#!/usr/bin/env python3
"""run_dehai.py -- validate the chain at DE-Hai (Hainich) against measurements.

WHY THIS SITE. The Zurich benchmark compares against MODELLED bands: the four
shortwave components shipped in the repository are 98.7% determined by (cloud
fraction, cos z), so they are a partitioning model's output, not instrument
readings. DE-Hai measures diffuse shortwave (SW_DIF) directly and reports
eddy-covariance GPP, so two links in the chain become testable against
observation rather than against another model:

    1  the direct/diffuse partition, against measured SW_DIF
    2  simulated GPP, against GPP_NT_VUT_REF and GPP_DT_VUT_REF

SIF itself remains unvalidated here -- DE-Hai carries no fluorescence
spectrometer -- but it is then driven by a radiation partition and a
photosynthesis rate that have each been checked against measurement, which is a
different claim from the aggregate consistency the 40-tower benchmark gives.

DATA. ICOS Carbon Portal, ICOSETC_DE-Hai_FLUXNET_FLUXMET_HH_*.csv, or the
FLUXNET2015 FULLSET product. Half-hourly; this aggregates to hourly.

    python run_dehai.py --csv ICOSETC_DE-Hai_FLUXNET_FLUXMET_HH_2000-2025_v1.3.csv

SITE. Hainich, Germany. Deciduous broadleaf (Fagus sylvatica), 51.0792 N,
10.4522 E, 430 m, UTC+1. Peak LAI 5-6, Vcmax25 for European beech 45-70.
"""
from __future__ import annotations
import argparse, os, sys, warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
import numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "..", "Forcing_Prep"))
from chain import photo_sif
from partition import solar_geometry, erbs_diffuse_fraction, SOLAR_CONST

LAT, LON, ELEV, GMT = 51.0792, 10.4522, 430.0, 1.0
MISSING = -9999.0

def read_icos(path):
    raw = np.genfromtxt(path, delimiter=",", names=True, dtype=None, encoding="utf-8")
    c = raw.dtype.names
    def col(*names):
        for n in names:
            if n in c:
                v = np.asarray(raw[n], float)
                return np.where(v <= MISSING + 1, np.nan, v)
        return None
    ts = np.asarray(raw["TIMESTAMP_START"], str)
    return dict(
        yr=np.array([int(s[0:4]) for s in ts]), mo=np.array([int(s[4:6]) for s in ts]),
        dy=np.array([int(s[6:8]) for s in ts]),
        hr=np.array([int(s[8:10]) + int(s[10:12])/60 for s in ts], float),
        SW_IN=col("SW_IN_F","SW_IN"), SW_DIF=col("SW_DIF"),
        PPFD_IN=col("PPFD_IN"), PPFD_DIF=col("PPFD_DIF"),
        TA=col("TA_F","TA"), VPD=col("VPD_F","VPD"), PA=col("PA_F","PA"),
        GPP_NT=col("GPP_NT_VUT_REF"), GPP_DT=col("GPP_DT_VUT_REF"),
        LAI=col("LAI"))

def calibrate_offset(SW, yr, mo, dy, hr):
    best = None
    for off in np.arange(-2.0, 2.01, 0.5):
        h, _ = solar_geometry(yr, mo, dy, hr+off, LAT, LON, GMT)
        cz = np.maximum(np.sin(h), 0.0)
        bad = int((np.nan_to_num(SW) > 20).__and__(cz <= 0).sum()
                  + ((cz > 0.3) & (np.nan_to_num(SW) <= 0)).sum())
        if best is None or bad < best[1]: best = (off, bad)
    return best

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--csv", required=True)
    ap.add_argument("--vmax", type=float, default=55.0)
    ap.add_argument("--laimax", type=float, default=5.5)
    ap.add_argument("--precoll", type=float, default=0.6)
    ap.add_argument("--kopt", type=float, default=0.5)
    a = ap.parse_args()
    D = read_icos(a.csv)
    n = len(D["yr"]); print(f"{n} records, {D['yr'].min()}-{D['yr'].max()}")
    for k in ("SW_DIF","GPP_NT","PPFD_IN"):
        print(f"  {k:8} {'present' if D[k] is not None else 'ABSENT'}"
              + (f", {np.isfinite(D[k]).mean()*100:.1f}% finite" if D[k] is not None else ""))

    off, bad = calibrate_offset(D["SW_IN"], D["yr"], D["mo"], D["dy"], D["hr"])
    print(f"\n  timestamp offset {off:+.1f} h ({bad} physically impossible hours)")
    h, r = solar_geometry(D["yr"], D["mo"], D["dy"], D["hr"]+off, LAT, LON, GMT)
    cz = np.maximum(np.sin(h), 0.0); I0 = SOLAR_CONST*r*cz
    SW = D["SW_IN"]; kt = np.where(I0 > 1, SW/np.maximum(I0, 1e-9), 0.0)

    # ---- TEST 1: direct/diffuse partition against MEASURED SW_DIF
    if D["SW_DIF"] is not None:
        obs_fd = D["SW_DIF"]/np.maximum(SW, 1e-9)
        m = np.isfinite(obs_fd) & (SW > 20) & (obs_fd >= 0) & (obs_fd <= 1)
        fd_erbs = erbs_diffuse_fraction(kt)
        rmse = lambda p: float(np.sqrt(np.nanmean((p[m]-obs_fd[m])**2)))
        print(f"\n  TEST 1  diffuse fraction vs MEASURED SW_DIF, n={m.sum()}")
        print(f"    Erbs (1982)          RMSE {rmse(fd_erbs):.4f}   bias "
              f"{np.nanmean((fd_erbs-obs_fd)[m]):+.4f}")
        from tcpart import TCPartition, predictors, _design_fd, fit_logit
        ktp, mm, lm, czp, pers = predictors(SW, h, r, np.exp(-ELEV/8434.0))
        idx = np.arange(n)                      # split by position, not by year:
        cut = idx[m][int(0.6*m.sum())] if m.sum() > 10 else n   # a year-based split
        tr = m & (idx <= cut)                   # is empty for single-year files
        b = fit_logit(_design_fd(ktp[tr], lm[tr], czp[tr], pers[tr]), obs_fd[tr], SW[tr])
        fd_fit = 1/(1+np.exp(-np.clip(_design_fd(ktp,lm,czp,pers)@b,-40,40)))
        te = m & ~tr
        print(f"    calibrated, held out RMSE "
              f"{float(np.sqrt(np.nanmean((fd_fit[te]-obs_fd[te])**2))):.4f}   "
              f"bias {np.nanmean((fd_fit-obs_fd)[te]):+.4f}   n={te.sum()}")
        fd_use = fd_fit
    else:
        print("\n  TEST 1 skipped: no SW_DIF column")
        fd_use = erbs_diffuse_fraction(kt)

    # ---- LAI: measured if present, else beech phenology
    if D["LAI"] is not None and np.isfinite(D["LAI"]).mean() > 0.3:
        LAI = np.where(np.isfinite(D["LAI"]), D["LAI"], 0.5); src = "measured"
    else:
        doy = np.array([sum([31,28,31,30,31,30,31,31,30,31,30,31][:mo-1])+dd
                        for mo, dd in zip(D["mo"], D["dy"])], float)
        g = np.clip((doy-110)/30, 0, 1)*np.clip((300-doy)/30, 0, 1)
        LAI = 0.3 + (a.laimax-0.3)*g; src = "beech phenology, DOY 110-300"
    print(f"\n  LAI: {src}, range {np.nanmin(LAI):.2f}-{np.nanmax(LAI):.2f}")

    # ---- TEST 2: simulated GPP against tower GPP
    PAR = 0.46*SW
    PARD = PAR*fd_use; PARB = PAR - PARD
    Fsun = (1-np.exp(-a.kopt*LAI))/np.maximum(a.kopt*LAI, 1e-9)
    Ta = D["TA"]; Ds = D["VPD"]*100.0 if D["VPD"] is not None else 1000.0
    Pre = D["PA"]*10.0 if D["PA"] is not None else 1013.0
    ok = np.isfinite(Ta) & np.isfinite(SW)
    Tc = np.clip(np.nan_to_num(Ta, nan=10.0), 0.1, 40)
    Dsc = np.clip(np.nan_to_num(Ds, nan=1000.0), 10, 6000)
    Prc = np.nan_to_num(Pre, nan=1013.0)
    lit = ok & (PAR > 5)
    A_s,_,F_s = photo_sif(np.where(lit, PARB/np.maximum(Fsun,1e-6)+PARD, 0), 400., Tc, Dsc, Prc, Vmax=a.vmax)
    A_h,_,F_h = photo_sif(np.where(lit, PARD, 0), 400., Tc, Dsc, Prc, Vmax=a.vmax)
    A_s,A_h,F_s,F_h = [np.where(lit, x, 0.0) for x in (A_s,A_h,F_s,F_h)]
    om = 0.87*(1.0 - min(max(a.precoll,0),1)); Kv = 0.5*np.sqrt(max(1-om,1e-6)); Ksv = a.kopt+Kv
    As = (1-np.exp(-a.kopt*LAI))/a.kopt; As_e = (1-np.exp(-Ksv*LAI))/Ksv
    Ah = LAI-As; Ah_e = (1-np.exp(-Kv*LAI))/Kv - As_e
    fs = As_e/np.maximum(As,1e-9); fh = Ah_e/np.maximum(Ah,1e-9)
    GPP_mod = A_s*As + A_h*Ah
    SIF_toc = F_s*As*fs + F_h*Ah*fh
    synth = os.path.basename(a.csv).lower().find("synth") >= 0
    if synth:
        print("\n  NOTE: filename contains 'synth'. TEST 2 compares against generated")
        print("  GPP and is a harness check only; the numbers carry no meaning.")
    for nm in ("GPP_NT","GPP_DT"):
        if D[nm] is None: continue
        g = D[nm]; mm2 = lit & np.isfinite(g) & (g > -5)
        if mm2.sum() < 100: continue
        bias = float(np.nanmean((GPP_mod-g)[mm2]))
        rm = float(np.sqrt(np.nanmean((GPP_mod-g)[mm2]**2)))
        rr = float(np.corrcoef(GPP_mod[mm2], g[mm2])[0,1])
        sl = np.linalg.lstsq(np.c_[g[mm2], np.ones(mm2.sum())], GPP_mod[mm2], rcond=None)[0][0]
        print(f"\n  TEST 2  simulated GPP vs tower {nm}, n={mm2.sum()}")
        print(f"    r2 {rr**2:.4f}   slope {sl:.3f}   bias {bias:+.3f}   RMSE {rm:.3f}  umol m-2 s-1")
        print(f"    mean modelled {GPP_mod[mm2].mean():.2f}   mean tower {g[mm2].mean():.2f}")

    # ---- daily SIF-GPP, in the 40-tower benchmark units
    day = D["yr"]*10000 + D["mo"]*100 + D["dy"]
    _, inv = np.unique(day, return_inverse=True); cnt = np.bincount(inv)
    dt_h = 0.5 if np.median(np.diff(D["hr"][:48])) < 0.75 else 1.0
    G = np.bincount(inv, weights=np.nan_to_num(GPP_mod))*3600*dt_h*12.011e-6
    S = np.bincount(inv, weights=np.nan_to_num(SIF_toc))/np.maximum(cnt,1)
    full = cnt >= (23/dt_h); G, S = G[full], S[full]
    o = (G > 0.2) & (S > 1e-4)
    if o.sum() > 100:
        sl, ic = np.linalg.lstsq(np.c_[S[o], np.ones(o.sum())], G[o], rcond=None)[0]
        rr = np.corrcoef(S[o], G[o])[0,1]
        print(f"\n  daily SIF-GPP, n={o.sum()} days")
        print(f"    slope {sl:.2f}   r2 {rr**2:.4f}   "
              f"{'IN RANGE' if 11.91<=sl<=68.59 else 'outside'} (11.91-68.59)")

if __name__ == "__main__":
    main()
