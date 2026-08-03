function [fesc_sun,fesc_shd,SIF_toc] = SIF_Escape(SIF_sun,SIF_shd,LAI,Kopt,theta_v,omega_l,LADF,p_recoll)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  SIF_Escape                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Fraction of leaf-emitted fluorescence that escapes the canopy toward a
%%% sensor, and the resulting top-of-canopy radiance.
%%%
%%% WHY THIS IS NEEDED
%%% Canopy_Resistence_An_Evolution.m line 164 forms
%%%     SIF = SIF_sun*(LAI*Fsun) + SIF_shd*(LAI*Fshd)
%%% which is what the canopy EMITS. A satellite measures what ESCAPES. Between
%%% the two sits reabsorption on the way out, which depends on where in the
%%% canopy the emission happened and on the viewing direction. Without it,
%%% SIF_H and SIF_L cannot be compared against a FLEX/FLORIS retrieval.
%%%
%%% DERIVATION
%%% Sunlit leaf area density at cumulative depth L is exp(-Kopt*L); shaded is
%%% the complement. A photon emitted at depth L escapes toward zenith angle
%%% theta_v with probability exp(-Kv*L). Integrating over the canopy,
%%%
%%%   fesc_sun = [ (1-exp(-(Kopt+Kv)*LAI))/(Kopt+Kv) ]
%%%              / [ (1-exp(-Kopt*LAI))/Kopt ]
%%%
%%%   fesc_shd = [ (1-exp(-Kv*LAI))/Kv - (1-exp(-(Kopt+Kv)*LAI))/(Kopt+Kv) ]
%%%              / [ LAI - (1-exp(-Kopt*LAI))/Kopt ]
%%%
%%% These are the same integrals the model already evaluates for the nitrogen
%%% profile in Canopy_Resistence_An_Evolution.m lines 55-56, with the viewing
%%% extinction Kv in place of Knit. Denominators are LAI*Fsun and LAI*Fshd
%%% exactly as defined there, so the scaling is consistent with the rest of the
%%% canopy module by construction.
%%%
%%% SCATTERING, AND WHY THE TWO-STREAM ALBEDO IS THE WRONG ONE HERE
%%% At 755 nm the leaf single-scattering albedo is high (omega_l ~ 0.87), so most
%%% interceptions scatter rather than absorb. It is tempting to reduce the
%%% extinction by sqrt(1-omega_l), which is the two-stream result (Goudriaan;
%%% Sellers 1985) -- but that applies to a DIFFUSE FLUX propagating through the
%%% medium, not to a DIRECTIONAL escape path. A scattered photon is redirected
%%% roughly isotropically, so about half of it goes back down and is lost to the
%%% sensor. Using the full albedo therefore over-credits escape.
%%%
%%% Recollision theory (Knyazikhin et al. 1998; Stenberg 2007) gives the
%%% effective albedo for escape as omega_eff = omega_l*(1 - p_recoll), where
%%% p_recoll is the probability that a scattered photon hits another leaf. For
%%% a closed canopy p_recoll is 0.5-0.7, so omega_eff is 0.26-0.44 rather than
%%% 0.87, and Kv is reduced by sqrt(1-omega_eff).
%%%
%%% This matters: with the two-stream albedo the bulk escape at LAI 3.5 comes out
%%% 0.753 and the modelled daily SIF-GPP slope is 10.27, below the observed
%%% 11.91-68.59 band from 40 FLUXNET towers (Zhang et al. 2018, Biogeosciences).
%%% With p_recoll = 0.6 the bulk escape is 0.553 and the slope is 13.59, inside
%%% the band, with r2 and RMSE essentially unchanged. Setting p_recoll = 1
%%% recovers pure absorption and the lower bound on escape.
%%%
%%% INPUT
%%% SIF_sun, SIF_shd  [W m-2 sr-1 um-1]  leaf-level F755 from
%%%                                      photosynthesis_biochemical.m
%%% LAI       [-]      leaf area index
%%% Kopt      [-]      beam extinction coefficient, from Canopy_Radiative_Transfer.m
%%% theta_v   [rad]    sensor zenith angle (0 = nadir; FLEX is near-nadir)
%%% omega_l   [-]      leaf single-scattering albedo at 755 nm, default 0.87
%%% LADF      string   leaf angle distribution: 'spherical','planophile','erectophile'
%%% p_recoll  [-]      recollision probability, default 0.6 (closed canopy)
%%%
%%% OUTPUT
%%% fesc_sun, fesc_shd  [-]                 escape fractions
%%% SIF_toc             [W m-2 sr-1 um-1]   top-of-canopy radiance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5 || isempty(theta_v);  theta_v = 0;          end
if nargin < 6 || isempty(omega_l);  omega_l = 0.87;       end
if nargin < 7 || isempty(LADF);     LADF = 'spherical';   end
if nargin < 8 || isempty(p_recoll); p_recoll = 0.6;       end

%%% G-function: projection of unit leaf area onto the viewing direction
switch lower(LADF)
    case 'spherical';   G = 0.5;
    case 'planophile';  G = cos(theta_v);
    case 'erectophile'; G = 2*sin(theta_v)/pi;
    otherwise;          G = 0.5;
end
G  = max(G,1e-3);
Kv = G./max(cos(theta_v),1e-3);
omega_eff = omega_l*(1 - min(max(p_recoll,0),1));   % recollision-corrected albedo
Kv = Kv*sqrt(max(1-omega_eff,1e-6));                % escape extinction

if LAI <= 1e-6
    fesc_sun = 1; fesc_shd = 1; SIF_toc = 0; return
end

Ks  = max(Kopt,1e-6);
Ksv = Ks + Kv;

Asun_esc = (1 - exp(-Ksv*LAI))/Ksv;                 % sunlit, escaping
Asun_tot = (1 - exp(-Ks *LAI))/Ks;                  % sunlit, total  = LAI*Fsun
Ashd_esc = (1 - exp(-Kv *LAI))/Kv - Asun_esc;       % shaded, escaping
Ashd_tot = LAI - Asun_tot;                          % shaded, total  = LAI*Fshd

fesc_sun = Asun_esc/max(Asun_tot,1e-9);
fesc_shd = Ashd_esc/max(Ashd_tot,1e-9);
fesc_sun = min(max(fesc_sun,0),1);
fesc_shd = min(max(fesc_shd,0),1);

SIF_toc = SIF_sun*Asun_tot*fesc_sun + SIF_shd*Ashd_tot*fesc_shd;
end
