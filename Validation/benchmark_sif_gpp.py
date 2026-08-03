#!/usr/bin/env python3
"""benchmark_sif_gpp.py -- the chain against the published 40-tower benchmark.

Reference: Zhang, Joiner, Alemohammad, Zhou, Gentine (2018), Biogeosciences 15,
5779-5800, "A global spatially contiguous solar-induced fluorescence (CSIF)
dataset using neural networks". Evaluating CSIF against GPP at 40 FLUXNET
tier-1 towers they report a regression slope spanning

    11.91 to 68.59   g C m-2 day-1 per mW m-2 nm-1 sr-1

with per-site r2 from 0.01 to 0.93 (median 0.64) and mean RMSE 1.67 g C m-2 d-1.

This aggregates the modelled hourly chain to those units and compares. Nothing
in the chain is fitted to the relationship.
"""
from __future__ import annotations
import argparse, os, sys
import numpy as np, scipy.io as sio
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chain import photo_sif, escape

LO, HI, R2_MED, RMSE_REF = 11.91, 68.59, 0.64, 1.67

def daily(forcing, LAI=None, Vmax=55.0, Kopt=0.5, omega=0.87, Ca=400.0, p_recoll=0.6):
    d = sio.loadmat(forcing); f = lambda k: np.asarray(d[k]).ravel().astype(float)
    PARB, PARD = f('PARB'), f('PARD')
    Ta, ea, es, Pre = f('Ta'), f('ea'), f('esat'), f('Pre')
    day = np.floor(f('D')).astype(np.int64)
    PAR = PARB + PARD; Ds = np.maximum(es-ea, 0)
    if LAI is None:                       # phenology from the site parameters
        import datetime as _dt
        from phenology import grass_LAI
        doy = np.array([(_dt.datetime.fromordinal(int(x)-366)).timetuple().tm_yday
                        for x in f('D')])
        LAI, _, _, _ = grass_LAI(Ta, day, doy)
    LAI = np.asarray(LAI, float)*np.ones_like(PAR)
    Fsun = (1-np.exp(-Kopt*LAI))/np.maximum(Kopt*LAI, 1e-9)
    PAR_sun = np.where(PAR > 0, PARB/np.maximum(Fsun,1e-6)+PARD, 0.0)
    valid = (Ta > -20) & (Ta < 45); lit = (PAR > 5) & valid
    Tc = np.clip(Ta, 0.1, 40)
    A_s,_,F_s = photo_sif(np.where(valid, PAR_sun, 0), Ca, Tc, Ds, Pre, Vmax=Vmax)
    A_h,_,F_h = photo_sif(np.where(valid, PARD, 0),    Ca, Tc, Ds, Pre, Vmax=Vmax)
    A_s,A_h,F_s,F_h = [np.where(lit, x, 0.0) for x in (A_s,A_h,F_s,F_h)]
    om_eff = omega*(1.0 - min(max(p_recoll,0.0),1.0))
    Kv = 0.5*np.sqrt(max(1-om_eff,1e-6)); Ksv = Kopt + Kv
    As = (1-np.exp(-Kopt*LAI))/Kopt
    As_e = (1-np.exp(-Ksv*LAI))/Ksv
    Ah = LAI - As
    Ah_e = (1-np.exp(-Kv*LAI))/Kv - As_e
    fs = As_e/np.maximum(As,1e-9); fh = Ah_e/np.maximum(Ah,1e-9)
    _, inv = np.unique(day, return_inverse=True)
    cnt = np.bincount(inv)
    GPP = np.bincount(inv, weights=A_s*As + A_h*Ah)*3600*12.011e-6   # g C m-2 d-1
    SIF = np.bincount(inv, weights=F_s*As*fs + F_h*Ah*fh)/np.maximum(cnt,1)
    full = cnt >= 23
    GPP, SIF = GPP[full], SIF[full]
    ok = (GPP > 0.2) & (SIF > 1e-4)
    sl, ic = np.linalg.lstsq(np.c_[SIF[ok], np.ones(ok.sum())], GPP[ok], rcond=None)[0]
    r = np.corrcoef(SIF[ok], GPP[ok])[0,1]
    rmse = float(np.sqrt(np.mean((GPP[ok]-(sl*SIF[ok]+ic))**2)))
    return dict(slope=sl, r2=r**2, rmse=rmse, n=int(ok.sum()),
                mGPP=GPP[ok].mean(), mSIF=SIF[ok].mean())

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--forcing", default="Inputs/Data_Run_Zurich_Fluntern.mat")
    ap.add_argument("--lai", type=float, default=None,
                    help="constant LAI; omit to use the site phenology")
    ap.add_argument("--vmax", type=float, default=55.0)
    ap.add_argument("--precoll", type=float, default=0.6)
    a = ap.parse_args()
    print("Benchmark: Zhang et al. 2018 BG, 40 FLUXNET tier-1 towers")
    print(f"  slope {LO}-{HI} g C m-2 d-1 per mW m-2 nm-1 sr-1;"
          f" r2 median {R2_MED}; RMSE {RMSE_REF}\n")
    R = daily(a.forcing, LAI=a.lai, Vmax=a.vmax, p_recoll=a.precoll)
    print(f"  modelled, n={R['n']} days, LAI={a.lai or 'site phenology'}, "
          f"Vmax={a.vmax}, p_recoll={a.precoll}")
    print(f"    slope {R['slope']:8.2f}   {'IN RANGE' if LO<=R['slope']<=HI else 'OUTSIDE'}")
    print(f"    r2    {R['r2']:8.4f}   {'above' if R['r2']>R2_MED else 'below'} the median")
    print(f"    RMSE  {R['rmse']:8.2f}   reference {RMSE_REF}")
    print(f"\n  {'Vmax':>5} {'slope':>8} {'r2':>7} {'RMSE':>7}  verdict")
    for V in (20,30,40,55,65,80,120):
        r = daily(a.forcing, LAI=a.lai, Vmax=V, p_recoll=a.precoll)
        print(f"  {V:5} {r['slope']:8.2f} {r['r2']:7.3f} {r['rmse']:7.2f}  "
              f"{'IN RANGE' if LO<=r['slope']<=HI else 'outside'}")
