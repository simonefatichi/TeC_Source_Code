#!/usr/bin/env python3
"""phenology.py -- prescribed grass LAI from the site's own T&C parameters.

MOD_PARAM_ZURICH_SMA.m gives, for the low-vegetation layer:
    aSE_L      = 2      grass species
    Tlo_L      = 0.0    mean temperature for leaf onset  [C]
    Tls_L      = NaN    no temperature-driven leaf shed
    dmg_L      = 20     days of maximum growth
    LAI_min_L  = 0.1    minimum LAI
    Sl_L       = 0.035  specific leaf area [m2/gC]

Holding LAI at 4 year-round contradicts LAI_min_L = 0.1 by a factor of forty in
dormancy, and assigns a full canopy to grass that is not there. In the SIF-GPP
regression that inflates GPP at the low-SIF end and flattens the slope, which is
the wrong direction for comparison against tower data.

This drives LAI from a running-mean temperature threshold at Tlo, with a
logistic build over dmg days and autumn senescence, bounded by LAI_min and a
prescribed peak. It is a prescribed phenology, not T&C's prognostic
VEGETATION_DYNAMIC, and it is used here only so the benchmark is run on a
canopy consistent with the site's own parameters.
"""
from __future__ import annotations
import numpy as np

def running_mean_daily(Ta_hourly, day_idx, window=7):
    u, inv = np.unique(day_idx, return_inverse=True)
    cnt = np.bincount(inv)
    Td = np.bincount(inv, weights=Ta_hourly)/np.maximum(cnt, 1)
    k = np.ones(window)/window
    Tsm = np.convolve(np.r_[np.repeat(Td[0], window), Td], k, mode="same")[window:]
    return u, Td, Tsm[:len(Td)]

def grass_LAI(Ta_hourly, day_idx, doy_hourly, Tlo=0.0, dmg=20, LAI_min=0.1,
              LAI_max=3.5, sen_doy=270, sen_len=45):
    """Daily LAI, mapped back onto the hourly index."""
    u, Td, Tsm = running_mean_daily(Ta_hourly, day_idx)
    _, inv = np.unique(day_idx, return_inverse=True)
    doy_d = np.bincount(inv, weights=doy_hourly)/np.maximum(np.bincount(inv), 1)
    LAI_d = np.full(len(u), LAI_min)
    grow = 0.0
    for i in range(len(u)):
        active = (Tsm[i] > Tlo) and (doy_d[i] < sen_doy)
        if active:
            grow = min(grow + 1.0/dmg, 1.0)
        else:
            grow = max(grow - 1.0/sen_len, 0.0)
        LAI_d[i] = LAI_min + (LAI_max - LAI_min)*grow
    return LAI_d[inv], LAI_d, u, doy_d
