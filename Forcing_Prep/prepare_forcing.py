#!/usr/bin/env python3
"""prepare_forcing.py -- FLUXNET/ICOS hourly CSV -> a T&C-ready forcing .mat.

Emits every variable prova_Rural_Zurich.m loads:
    Date Pr Ta Ws ea esat Tdew Pre N SAB1 SAB2 SAD1 SAD2 PARB PARD
plus Lat Lon DeltaGMT Zbas.

The radiation partition is the part that needs care. Total shortwave is all a
tower gives; T&C wants four bands and two PAR streams. See tcpart.py for why
three fractions are needed rather than one.

    python prepare_forcing.py site.csv out.mat --lat 51.08 --lon 10.45 \
        --elev 430 --gmt 1 [--calib coeffs_zurich.json] [--fit-here]

--fit-here refits the coefficients on this site, for towers that measure
diffuse (SW_DIF). Without it the shipped Zurich coefficients are applied, which
is still better than a fixed 0.45 but is an extrapolation across climate.
"""
from __future__ import annotations
import argparse, json, sys
import numpy as np
from scipy.io import savemat
from partition import solar_geometry
from tcpart import TCPartition, predictors, _design_fd, _design_fv, _logistic

FLUX = dict(SW_IN="SW_IN_F", TA="TA_F", WS="WS_F", PA="PA_F", P="P_F",
            VPD="VPD_F", RH="RH", SW_DIF="SW_DIF")

def esat_hPa(T):                       # Magnus, T in C, returns Pa
    return 610.94*np.exp(17.625*T/(T+243.04))

def dewpoint(ea_Pa):
    lg = np.log(np.maximum(ea_Pa, 1e-3)/610.94)
    return 243.04*lg/(17.625-lg)

def calibrate_hour_offset(Rsw, yr, mo, dy, hr, lat, lon, gmt):
    """The offset minimising physically impossible hours. Silent 19% error if wrong."""
    best, tbl = None, []
    for off in np.arange(-2.0, 2.01, 0.5):
        h, _ = solar_geometry(yr, mo, dy, hr+off, lat, lon, gmt)
        cz = np.maximum(np.sin(h), 0.0)
        bad = int(((Rsw > 20) & (cz <= 0)).sum() + ((cz > 0.3) & (Rsw <= 0)).sum())
        tbl.append((off, bad))
        if best is None or bad < best[1]: best = (off, bad)
    return best[0], tbl

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("csv"); ap.add_argument("out")
    ap.add_argument("--lat", type=float, required=True)
    ap.add_argument("--lon", type=float, required=True)
    ap.add_argument("--elev", type=float, required=True)
    ap.add_argument("--gmt", type=float, required=True)
    ap.add_argument("--calib", default="coeffs_zurich.json")
    ap.add_argument("--fit-here", action="store_true")
    a = ap.parse_args()

    raw = np.genfromtxt(a.csv, delimiter=",", names=True, dtype=None, encoding="utf-8")
    cols = raw.dtype.names
    def col(key, req=True):
        n = FLUX.get(key, key)
        for c in (n, key):
            if c in cols: return np.asarray(raw[c], float)
        if req: sys.exit(f"missing column {n}")
        return None
    ts = np.asarray(raw["TIMESTAMP_START" if "TIMESTAMP_START" in cols else "TIMESTAMP"], str)
    yr = np.array([int(s[0:4]) for s in ts]); mo = np.array([int(s[4:6]) for s in ts])
    dy = np.array([int(s[6:8]) for s in ts]); hr = np.array([int(s[8:10]) + int(s[10:12])/60
                                                            for s in ts], float)
    for k in ("SW_IN","TA","WS","PA","P"):
        pass
    Rsw = np.maximum(col("SW_IN"), 0.0); Ta = col("TA"); Ws = np.maximum(col("WS"), 0.01)
    Pre = col("PA")*10.0                                   # kPa -> hPa
    Pr = np.maximum(col("P"), 0.0)
    es = esat_hPa(Ta)
    vpd = col("VPD", req=False)
    ea = es - vpd*100.0 if vpd is not None else es*col("RH")/100.0
    ea = np.clip(ea, 1.0, es)
    Tdew = dewpoint(ea)

    off, tbl = calibrate_hour_offset(Rsw, yr, mo, dy, hr, a.lat, a.lon, a.gmt)
    print(f"hour offset calibrated to {off:+.1f} h "
          f"({dict(tbl)[off]} inconsistent hours of {len(Rsw)})")
    h, r = solar_geometry(yr, mo, dy, hr+off, a.lat, a.lon, a.gmt)
    prr = np.exp(-a.elev/8434.0)

    swdif = col("SW_DIF", req=False)
    if a.fit_here and swdif is not None:
        print("refitting the diffuse fraction on this site's measured SW_DIF")
        M = TCPartition()
        fd = np.clip(swdif/np.maximum(Rsw, 1e-9), 0, 1)
        kt, m, lm, cz, pers = predictors(Rsw, h, r, prr)
        ok = Rsw > 20
        from tcpart import fit_logit
        M.b_fd = fit_logit(_design_fd(kt[ok], lm[ok], cz[ok], pers[ok]), fd[ok], Rsw[ok])
        C = json.load(open(a.calib))
        M.b_vdir = np.array(C["fvis_dir"]); M.b_vdif = np.array(C["fvis_dif"])
        M.b_pdir = np.array(C["par_dir"]);  M.b_pdif = np.array(C["par_dif"])
    else:
        C = json.load(open(a.calib))
        M = TCPartition()
        M.b_fd = np.array(C["f_diff"]); M.b_vdir = np.array(C["fvis_dir"])
        M.b_vdif = np.array(C["fvis_dif"]); M.b_pdir = np.array(C["par_dir"])
        M.b_pdif = np.array(C["par_dif"])
        print(f"applying calibration from {a.calib} (site: {C.get('site','?')})")
    B = M.predict(Rsw, h, r, prr)

    kt, _, _, cz, _ = predictors(Rsw, h, r, prr)
    N = np.clip(1.0 - kt/0.75, 0.0, 1.0)                  # crude, only used for LW
    Date = np.array([_datenum(y, m_, d_, hh) for y, m_, d_, hh in zip(yr, mo, dy, hr)])

    out = dict(Date=Date.reshape(-1,1), Pr=Pr.reshape(-1,1), Ta=Ta.reshape(-1,1),
               Ws=Ws.reshape(-1,1), ea=ea.reshape(-1,1), esat=es.reshape(-1,1),
               Tdew=Tdew.reshape(-1,1), Pre=Pre.reshape(1,-1), N=N.reshape(-1,1),
               Lat=a.lat, Lon=a.lon, DeltaGMT=a.gmt, Zbas=a.elev)
    for k in ("SAB1","SAB2","SAD1","SAD2","PARB","PARD"):
        out[k] = B[k].reshape(1,-1)
    savemat(a.out, out)
    s = sum(B[k] for k in ("SAB1","SAB2","SAD1","SAD2"))
    print(f"wrote {a.out}: {len(Rsw)} steps, {len(out)} variables")
    print(f"closure max |SAB1+SAB2+SAD1+SAD2 - Rsw| = {np.abs(s-Rsw).max():.2e}")

def _datenum(y, m, d, h):
    import datetime as dt
    t = dt.datetime(int(y), int(m), int(d)) + dt.timedelta(hours=float(h))
    return t.toordinal() + 366 + (t.hour*3600+t.minute*60)/86400.0

if __name__ == "__main__":
    main()
