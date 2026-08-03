#!/usr/bin/env python3
"""sif_gpp.py -- the emergent SIF-GPP relationship of the T&C chain.

Runs photosynthesis_biochemical.m's Farquhar + Lee et al. 2015 fluorescence
block, the sunlit/shaded scaling of Canopy_Resistence_An_Evolution.m line 164,
and SIF_Escape.m, on real forcing. Nothing in the chain is fitted to the
SIF-GPP relationship, so the relationship is a prediction.

    python sif_gpp.py [--forcing PATH] [--lai 4] [--vmax 55]
"""
from __future__ import annotations
import argparse, sys, os
import numpy as np, scipy.io as sio, datetime as dt
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from chain import photo_sif, escape

def load(path):
    d = sio.loadmat(path); f = lambda k: np.asarray(d[k]).ravel().astype(float)
    D = f('D'); n = len(D); mo = np.empty(n, int)
    for i, dv in enumerate(D):
        mo[i] = (dt.datetime.fromordinal(int(dv)-366)+dt.timedelta(days=float(dv)%1)).month
    return dict(PARB=f('PARB'), PARD=f('PARD'), Ta=f('Ta'), ea=f('ea'),
                esat=f('esat'), Pre=f('Pre'), mo=mo)

def run(F, LAI=4.0, Vmax=55.0, Kopt=0.5, omega=0.87, Ca=400.0, theta_v=0.0):
    PAR = F['PARB']+F['PARD']; Ds = np.maximum(F['esat']-F['ea'], 0)
    m = (PAR > 5) & (F['Ta'] > 0) & (F['Ta'] < 40)
    Fsun = (1-np.exp(-Kopt*LAI))/(Kopt*LAI)
    PAR_sun = F['PARB'][m]/max(Fsun,1e-6) + F['PARD'][m]
    PAR_shd = F['PARD'][m]
    A_s,_,Fl_s = photo_sif(PAR_sun, Ca, F['Ta'][m], Ds[m], F['Pre'][m], Vmax=Vmax)
    A_h,_,Fl_h = photo_sif(PAR_shd, Ca, F['Ta'][m], Ds[m], F['Pre'][m], Vmax=Vmax)
    fs, fh, As, Ah = escape(LAI, Kopt, theta_v, omega)
    GPP = A_s*As + A_h*Ah
    SIF = Fl_s*As*fs + Fl_h*Ah*fh
    ok = (GPP > 0.5) & (SIF > 1e-5)
    sl, ic = np.linalg.lstsq(np.c_[SIF[ok], np.ones(ok.sum())], GPP[ok], rcond=None)[0]
    r = np.corrcoef(SIF[ok], GPP[ok])[0,1]
    return dict(slope=sl, r2=r**2, n=int(ok.sum()), fesc_sun=fs, fesc_shd=fh,
                fesc_bulk=(As*fs+Ah*fh)/(As+Ah), SIF=SIF, GPP=GPP, ok=ok,
                mo=F['mo'][m], mean_SIF=SIF[ok].mean(), mean_GPP=GPP[ok].mean())

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--forcing", default="/home/claude/TeC/Inputs/Data_Run_Zurich_Fluntern.mat")
    ap.add_argument("--lai", type=float, default=4.0)
    ap.add_argument("--vmax", type=float, default=55.0)
    a = ap.parse_args()
    F = load(a.forcing)
    R = run(F, LAI=a.lai, Vmax=a.vmax)
    print(f"n = {R['n']}  LAI = {a.lai}  Vmax = {a.vmax}")
    print(f"  escape: sunlit {R['fesc_sun']:.4f}  shaded {R['fesc_shd']:.4f}  bulk {R['fesc_bulk']:.4f}")
    print(f"  GPP = {R['slope']:.2f} x SIF     r2 = {R['r2']:.4f}")
    print(f"  mean SIF {R['mean_SIF']:.3f} W m-2 sr-1 um-1   mean GPP {R['mean_GPP']:.2f} umol m-2 s-1")
    print("\n  SENSITIVITY")
    print(f"  {'LAI':>5} {'slope':>7} {'r2':>7} {'bulk fesc':>10}")
    for L in (1,2,3,4,6,8):
        r = run(F, LAI=L, Vmax=a.vmax)
        print(f"  {L:5} {r['slope']:7.2f} {r['r2']:7.3f} {r['fesc_bulk']:10.3f}")
    print(f"  {'Vmax':>5} {'slope':>7} {'r2':>7}")
    for V in (30,40,55,80,120):
        r = run(F, LAI=a.lai, Vmax=V)
        print(f"  {V:5} {r['slope']:7.2f} {r['r2']:7.3f}")
