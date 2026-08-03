#!/usr/bin/env python3
"""chain.py -- T&C's photosynthesis + SIF chain, transcribed, driven by real forcing.

Reproduces photosynthesis_biochemical.m (Farquhar + the Lee et al. 2015 GCB
fluorescence block) and Canopy_Resistence_An_Evolution.m line 164, then applies
SIF_Escape.m. Driven by the Zurich forcing shipped with the repo.

The test is the emergent SIF-GPP relationship. It is not fitted anywhere in the
chain: fiF comes from the NPQ parameterisation, GPP from the Farquhar solution.
If the coupled model is right, the slope should land in the range reported from
towers and OCO-2 at 757-760 nm.
"""
from __future__ import annotations
import numpy as np

def photo_sif(IPAR_W, Ca, Ts, Ds, Pre, Vmax=55.0, CT=3, a1=6.0, go=0.01,
              rjv=1.9, Oa=210.0, Do=1000.0):
    """Farquhar with the Lee et al. fluorescence block. IPAR_W in W/m2."""
    IPAR = IPAR_W*4.57                                   # umol photons /s/m2
    Tk = Ts + 273.15
    def arrh(c, dH): return np.exp(c - dH/(0.008314*Tk))
    Kc = arrh(38.05, 79.43); Ko = arrh(20.30, 36.38)
    GAM_s = arrh(19.02, 37.83)
    GAM = 0.5*np.exp(-3.3801 + 5220.0/(Tk*8.314))*Oa*Kc/Ko
    kT = np.exp(26.35 - 65.33/(0.008314*Tk))/(1 + np.exp((0.71*Tk - 220.0)/(0.008314*Tk)))
    Vm = Vmax*np.exp(26.35 - 65.33/(0.008314*Tk))/(1 + np.exp((0.65*Tk - 200.0)/(0.008314*Tk)))
    Jmax = Vmax*rjv; Jm = Jmax*kT
    FI = 0.081                                            # intrinsic quantum efficiency
    Q = FI*IPAR                                           # umolCO2 /s/m2
    th = 0.9
    J = (Q + Jm - np.sqrt(np.maximum((Q+Jm)**2 - 4*th*Q*Jm, 0)))/(2*th)
    Rdark = 0.015*Vm
    Cc = 0.7*Ca
    for _ in range(30):
        Cc = np.maximum(Cc, 1e-3)
        JC = Vm*(Cc - GAM)/(Cc + Kc*(1 + Oa/Ko))
        JE = J*(Cc - GAM)/(4*(Cc + 2*GAM))
        A = np.minimum(JC, JE)
        An = A - Rdark
        gsCO2 = go + a1*An*Pre/((Cc - GAM)*(1 + Ds/Do))
        gsCO2 = np.maximum(gsCO2, go)
        Cc_new = Ca - An*Pre/np.maximum(gsCO2, 1e-6)
        Cc = 0.5*Cc + 0.5*np.clip(Cc_new, 1e-3, Ca)
    # ---- Lee et al. 2015 fluorescence, verbatim from lines 270-291
    Jfe = A*(Cc + 2*GAM)/np.maximum(Cc - GAM, 1e-6) if CT == 3 else A
    fiP0 = FI*4
    fiP = fiP0*Jfe/np.maximum(Q, 1e-9)
    dls = np.clip(1 - fiP/fiP0, 0, 1)
    kf = 0.05
    kd = np.maximum(0.03*Ts + 0.0773, 0.087)
    kn = (6.2473*dls - 0.5944)*dls
    fiF = kf/(kf + kd + kn)*(1 - fiP)
    SIF = IPAR*fiF
    k = 0.0375*Vmax + 8.25
    F755 = np.where(IPAR > 0, SIF/k, 0.0)
    return np.maximum(A, 0.0), np.maximum(An, 0.0), np.maximum(F755, 0.0)

def escape(LAI, Kopt=0.5, theta_v=0.0, omega=0.87, G=0.5, p_recoll=0.6):
    """Escape fractions. omega_eff = omega*(1-p_recoll): a scattered photon is
    redirected isotropically, so the two-stream albedo over-credits directional
    escape. p_recoll = 1 recovers pure absorption."""
    omega_eff = omega*(1.0 - min(max(p_recoll, 0.0), 1.0))
    Kv = G/max(np.cos(theta_v), 1e-3)*np.sqrt(max(1-omega_eff, 1e-6))
    Ks = max(Kopt, 1e-6); Ksv = Ks + Kv
    As_e = (1-np.exp(-Ksv*LAI))/Ksv; As_t = (1-np.exp(-Ks*LAI))/Ks
    Ah_e = (1-np.exp(-Kv*LAI))/Kv - As_e; Ah_t = LAI - As_t
    return As_e/max(As_t,1e-9), Ah_e/max(Ah_t,1e-9), As_t, Ah_t
