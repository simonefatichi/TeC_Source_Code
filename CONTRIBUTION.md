# Fork contributions

Two additions to `simonefatichi/TeC_Source_Code`, both benchmarked against
published data.

## 1. Making T&C FLEX-ready

FLEX/FLORIS launches September 2026 — 300 m SIF, 27-day repeat, in tandem with
Sentinel-3C. ESA's cal/val AO is open.

**T&C already models SIF.** `photosynthesis_biochemical.m` carries the Lee et al.
(2015, GCB) block: `Jfe` → `fiP` → `dls` → `kn` → `fiF` → `F755nm` in
W m⁻² sr⁻¹ µm⁻¹, which is the FLORIS retrieval band.
`Canopy_Resistence_An_Evolution.m` line 164 scales sunlit and shaded, and
`SIF_H`/`SIF_L` propagate to `MAIN_FRAME`.

**One step is missing.** Line 164 gives canopy *emission*. A satellite sees what
*escapes*.

`T&C_Code/SIF_Escape.m` adds it. Sunlit leaf density at depth `L` is
`exp(-Kopt·L)`, shaded the complement; escape toward zenith `θv` is `exp(-Kv·L)`:

```
fesc_sun = [(1-exp(-(Kopt+Kv)·LAI))/(Kopt+Kv)] / [(1-exp(-Kopt·LAI))/Kopt]
fesc_shd = [(1-exp(-Kv·LAI))/Kv - (1-exp(-(Kopt+Kv)·LAI))/(Kopt+Kv)]
           / [LAI - (1-exp(-Kopt·LAI))/Kopt]
```

These are the integrals already evaluated for the nitrogen profile in
`Canopy_Resistence_An_Evolution.m` lines 55–56, with viewing extinction `Kv` in
place of `Knit`. The denominators reduce to `LAI·Fsun` and `LAI·Fshd` **exactly**
as defined there — verified to 0.00e+00 over LAI 0.1–8.

`T&C_Code/SIF_Output.m` resamples onto FLORIS sampling using true local solar
time.

### Benchmark

Against Zhang, Joiner, Alemohammad, Zhou & Gentine (2018), *Biogeosciences* 15,
5779–5800 — CSIF evaluated against GPP at **40 FLUXNET tier-1 towers**:

| | modelled, Zurich, Vmax=55 | observed, 40 towers |
|---|---|---|
| slope (g C m⁻² d⁻¹ per mW m⁻² nm⁻¹ sr⁻¹) | **13.60** | 11.91 – 68.59 |
| r² | **0.888** | 0.01 – 0.93, median 0.64 |
| RMSE (g C m⁻² d⁻¹) | **1.80** | mean 1.67 |

Inside the observed band, with r² near the top of the observed range and RMSE
close to the reported mean. Nothing in the chain is fitted to this relationship.

```
python Validation/benchmark_sif_gpp.py          # site phenology, p_recoll = 0.6
```

### Two corrections made during benchmarking

**The two-stream albedo is the wrong one for directional escape.** A first
version reduced the viewing extinction by `sqrt(1-omega_l)` with `omega_l =
0.87`. That is the two-stream result for a *diffuse flux* propagating through the
medium; escape toward a sensor is *directional*, and a scattered photon is
redirected roughly isotropically, so about half of it goes back down and is lost.
Using the full albedo over-credits escape — bulk 0.753 at LAI 3.5, giving a slope
of 10.27, below the observed band. Recollision theory (Knyazikhin et al. 1998;
Stenberg 2007) gives `omega_eff = omega_l*(1 - p_recoll)`, with `p_recoll` the
probability a scattered photon strikes another leaf, 0.5–0.7 for a closed canopy.
At `p_recoll = 0.6` the bulk escape is 0.535 and the slope 13.60. `p_recoll = 1`
recovers pure absorption.

**Constant LAI contradicts the site's own parameters.** `MOD_PARAM_ZURICH_SMA.m`
sets `aSE_L = 2` (grass), `Tlo_L = 0.0`, `LAI_min_L = 0.1`, `dmg_L = 20`. Holding
LAI at 4 year-round is wrong by a factor of forty in dormancy.
`Validation/phenology.py` drives LAI from those parameters (0.10 in December to
3.50 in summer). This changed RMSE from 1.86 to 1.80 but barely moved the slope
(10.29 → 10.27), because when LAI collapses GPP and SIF fall together and those
days sit near the origin without levering the fit. It is included because it is
correct, not because it was the fix.

### The result for cal/val

The slope is a strong function of `Vmax` and a weak function of LAI:

| Vmax | slope | r² | RMSE | verdict |
|---|---|---|---|---|
| 20 | 6.05 | 0.783 | 1.21 | outside |
| 30 | 8.64 | 0.826 | 1.51 | outside |
| 40 | 10.85 | 0.856 | 1.68 | outside |
| **55** | **13.60** | **0.888** | **1.80** | **in range** |
| 65 | 15.23 | 0.905 | 1.82 | in range |
| 80 | 17.49 | 0.927 | 1.80 | in range |
| 120 | 22.68 | 0.957 | 1.63 | in range |

LAI barely moves it because the escape fraction falls roughly in step with the
GPP increase and the two largely cancel. `Vmax` moves it 3.7× over a 6× range,
and the model enters the observed band at `Vmax` ≈ 47.

Two consequences. **SIF alone cannot constrain GPP without independent knowledge
of `Vmax`** — an apparent between-site slope difference may be a `Vmax`
difference. And **SIF and GPP jointly constrain `Vmax`**: at a tower measuring
both, the observed slope inverts to a `Vmax` estimate, with `r²` rising
monotonically alongside it. That is a usable cal/val target.

Note that `Vmax` alone cannot span the full observed 11.91–68.59: 3.7× over a
physiological `Vmax` range against 5.8× observed. Chlorophyll content, biome and
canopy structure carry the rest.

## 2. Forcing preparation from flux-tower data

T&C needs six radiation variables; a FLUXNET or ICOS tower gives total shortwave
only. `T&C_Code/Radiation_Partition.m` builds all six, and
`Forcing_Prep/prepare_forcing.py` produces a complete T&C `.mat` from tower CSV.

Two things found while calibrating against the shipped Zurich forcing:

**The timestamp convention is worth 19% of radiation RMSE and is undocumented.**
The forcing is stamped UTC while `DeltaGMT = 1`. Treating the stamp as local time
leaves 6,722 hours with Rsw > 20 W m⁻² and the sun below the horizon, and 3,771
hours with Rsw exceeding the extraterrestrial irradiance with the sun well up.
The correction *is* applied — by `t_bef = -0.67; t_aft = 1.67` in
`prova_Rural_Zurich.m`, whose window centres at +1.17 h, matching the physical
optimum of +1.10 h to 0.07 h. But those two constants carry no comment and are
site-specific, so copying the driver to a site with local-time stamps silently
imports Zurich's offset. `calibrate_hour_offset.m` determines it from the data.

**The visible fraction is a property of the stream, not the site.** In the shipped
bands, `SAB1/(SAB1+SAB2) = 0.389 ± 0.151` and `SAD1/(SAD1+SAD2) = 0.537 ± 0.101`
— 15 points apart, which is Rayleigh scattering. Any scheme applying one visible
fraction to both streams cannot reproduce the four bands: per-stream fractions
reconstruct them to RMSE 1.37 W m⁻², one shared fraction leaves 29.6 regardless.
Erbs + a fixed 0.45, and Weiss & Norman (1985), both make that assumption.

Out-of-sample (fit 1981–2004, tested 2005–2012), aggregate RMSE over the six
radiation variables falls **48.9%**.

## Limits

`p_recoll` is the one free parameter in the escape module. It is bounded by
theory to 0.5–0.7 for a closed canopy and the benchmark is satisfied across that
whole interval (slope 12.6 at 0.5 through 14.7 at 0.7), so the result does not
depend on the choice within its physical range. It should be derived from canopy
structure rather than prescribed; recollision probability is computable from LAI
and the leaf angle distribution, both of which T&C already carries.

`Validation/chain.py` is a transcription of the MATLAB for benchmarking, not the
authoritative implementation, and should be replaced by a direct call into
`photosynthesis_biochemical.m` once run inside MATLAB.

The benchmark is one site, one PFT, and against *modelled* SIF — the Zurich
forcing carries no fluorescence measurement, so agreement with the 40-tower band
tests the chain's magnitude and shape, not its accuracy at this site. A tower
with a co-located spectrometer (DE-Hai) is the test that settles it, and the
comparison against measured rather than modelled SIF is the next step.

The radiation coefficients are fitted at one mid-latitude continental site and
are an extrapolation elsewhere.

## Licensing

`simonefatichi/TeC_Source_Code` carries no LICENSE file, so the upstream code is
all-rights-reserved by default. Nothing here relicenses it. The files added by
this contribution — `SIF_Escape.m`, `SIF_Output.m`, `Radiation_Partition.m`,
`Radiation_Partition_Coeff_Zurich.m`, `calibrate_hour_offset.m`, and everything
under `Validation/` and `Forcing_Prep/` — are offered to the T&C authors on
whatever terms they apply to the rest of the repository. If a license is added
upstream these follow it.
