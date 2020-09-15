%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Canopy_Radiative_Transfer  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Ref. Sellers, 1986 ; Bonan 1996; Oleson et al., 2004
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Iup,Idn,Kopt,om_vis_vg]=Canopy_Radiative_Transfer(PFT_opt,soil_alb,h_S,LAI,SAI,LAId,dw_SNO)
%%%%%%%%%%%%%%%%%%%%%%%%
%%% Two stream approximation 
%%%%%%% Dai 2004 Bonan 1996 Sellers 1985 Dickinson 1983 from Goudriaan 1977 
%%%% INPUT 
chiL=PFT_opt.chiL; %% Leaf Angle Distribution Parameter Ross, 1975
alf_lf_vis=PFT_opt.alf_lf_vis;
alf_lf_nir=PFT_opt.alf_lf_nir;
alf_st_vis =PFT_opt.alf_st_vis;
alf_st_nir =PFT_opt.alf_st_nir ;
alf_ld_vis =PFT_opt.alf_ld_vis;
alf_ld_nir =PFT_opt.alf_ld_nir ;
tau_lf_vis =PFT_opt.tau_lf_vis;
tau_lf_nir =PFT_opt.tau_lf_nir;
tau_st_vis =PFT_opt.tau_st_vis;
tau_st_nir =PFT_opt.tau_st_nir;
tau_ld_vis =PFT_opt.tau_ld_vis;
tau_ld_nir =PFT_opt.tau_ld_nir;
%%%-> ALBEDO_PROPERTIES 
soil_alb_dir_vis = soil_alb.dir_vis; 
soil_alb_dir_nir = soil_alb.dir_nir; 
soil_alb_dif_vis = soil_alb.dif_vis;
soil_alb_dif_nir = soil_alb.dif_nir;
%%% h_S 
%%% LAI 
%%% SAI 
fsnow= dw_SNO; %% [] Intercepted Snow Coefficient for H - Vegetation  
%%%% OUTPUT 
%Iup.dir_vis
%Idn.dir_vis
%Iup.dif_vis
%Idn.dif_vis
%Iup.dir_nir
%Idn.dir_nir
%Iup.dif_nir
%Idn.dif_nir
%Kopt 
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%[h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day] = SetSunVariables(Datam,DeltaGMT,Lon,Lat)
%%%%%%%%%%%%%%%%%%%
%%%  // Solar altitude or angle of radiation -- computed previously
%%%  // 'mu' is actually 'sinAlpha', the cosine of the zenith angle
%%%  // or, alternatively, the sine of the solar altitude
mu = sin(h_S); %mu is the cosine of the zenith angle of the incident beam
%%%%%%%%%%%%%%%%
%%% Estimate weighted leaf reflectance and transmittance
if (chiL > -0.4) && (chiL < 0.6)
    phi1 = 0.5 - 0.633*chiL - 0.33*chiL*chiL;
    phi2 = 0.877*(1-2*phi1);
    %%% Estimate the relative projected area of leaves and stems in acos(mu)
    G_mu = phi1 + phi2*mu;
    %%% Estimate the optical depth of direct beam per unit leaf and stem area
    Kopt = G_mu/mu;
    %%% The average inverse diffuse optical depth per unit leaf/stem area
    %%% Correction Dai et al., 2004
    if isequal(phi2,0) || isequal(phi1,0)
        if isequal(phi2,0)
            muav = 1/(2*phi1);
        else
            muav = 1/0.877;
        end
    else
        muav = (1/phi2)*(1 - (phi1/phi2)*log((phi1+phi2)/phi1));
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Get weighted canopy reflectance
omleaf = LAI/(LAI+SAI+LAId);
omstem = SAI/(LAI+SAI+LAId);
omlede = LAId/(LAI+SAI+LAId);
alf_vis = alf_lf_vis*omleaf + alf_st_vis*omstem + alf_ld_vis*omlede;
alf_nir = alf_lf_nir*omleaf + alf_st_nir*omstem + alf_ld_nir*omlede;
%%%%%%    
%%%Get weighted canopy transmittance
tau_vis = tau_lf_vis*omleaf + tau_st_vis*omstem +tau_ld_vis*omlede;
tau_nir = tau_lf_nir*omleaf + tau_st_nir*omstem +tau_ld_nir*omlede;
%%% Estimate CANOPY weighted scattering coeficient
om_vis_vg = alf_vis + tau_vis;
om_nir_vg = alf_nir + tau_nir;
%%% Estimate the single scattering albedo
tmp = (1 - (mu*phi1)/(mu*phi2+G_mu)*log((mu*phi1+mu*phi2+G_mu)/(mu*phi1)));
tmp = tmp*G_mu/(mu*phi2+G_mu);
alfs_mu_vis = 0.5*om_vis_vg*tmp;
alfs_mu_nir = 0.5*om_nir_vg*tmp;
%%% Estimate upscatter for DIFFUSE radiation
om_vis_bet_vg = 0.5*(om_vis_vg + (alf_vis-tau_vis)*((1+chiL)/2)*((1+chiL)/2));   
om_nir_bet_vg = 0.5*(om_nir_vg + (alf_nir-tau_nir)*((1+chiL)/2)*((1+chiL)/2));  
%%% Estimate upscatter for DIRECT BEAM radiation
om_vis_bet0_vg = alfs_mu_vis*(1 + muav*Kopt)/(muav*Kopt);
om_nir_bet0_vg = alfs_mu_nir*(1 + muav*Kopt)/(muav*Kopt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
om_vis_sn = 0.8; 
om_vis_bet_sn = 0.8*0.5; 
om_vis_bet0_sn = 0.8*0.5; 
om_nir_sn = 0.4; 
om_nir_bet_sn = 0.4*0.5; 
om_nir_bet0_sn = 0.4*0.5; 
%%%%%%%%%%%%%%%%%%%%%%%%
%%% NOTE: The above coefficients will have to be further
%%% used to account for snow presence on the canopy
%%% E.g. if snow occupies 'fsnow' fraction of the canopy:
%%% etc. See Bonan [1996],  section 2.1  or Oleson et al., [2004] section 3.1 for details
om_vis      = om_vis_vg*(1-fsnow)      + om_vis_sn*fsnow;
om_vis_bet  = om_vis_bet_vg*(1-fsnow)  + om_vis_bet_sn*fsnow;
om_vis_bet0 = om_vis_bet0_vg*(1-fsnow) + om_vis_bet0_sn*fsnow;
%%%%%%%%%%%%%%%
om_nir      = om_nir_vg*(1-fsnow)      + om_nir_sn*fsnow;
om_nir_bet  = om_nir_bet_vg*(1-fsnow)  + om_nir_bet_sn*fsnow;
om_nir_bet0 = om_nir_bet0_vg*(1-fsnow) + om_nir_bet0_sn*fsnow;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  %% Compute canopy UPWARD & DOWNWARD diffuse fluxes for incident
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  // -DIRECT BEAM- radiation flux (== VISIBLE range ==)
[H,sgm,s1]=ComputeRadCoeffs(om_vis,om_vis_bet,om_vis_bet0,soil_alb_dir_vis,muav,Kopt,LAI,SAI,LAId);
Iup.dir_vis = H(1)/sgm + H(2) + H(3);
Idn.dir_vis = (H(4)/sgm)*exp(-Kopt*(LAI+SAI+LAId)) + H(5)*s1 + H(6)/s1;
%%%% // Compute canopy UPWARD & DOWNWARD diffuse fluxes for incident
%%%  // -DIFFUSE BEAM- radiation flux (== VISIBLE range ==)
[H,sgm,s1]=ComputeRadCoeffs(om_vis,om_vis_bet,om_vis_bet0,soil_alb_dif_vis,muav,Kopt,LAI,SAI,LAId);
Iup.dif_vis = H(7) + H(8);
Idn.dif_vis = H(9)*s1 + H(10)/s1;
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%   // -DIRECT BEAM- radiation flux (== INFRARED range ==)
[H,sgm,s1]=ComputeRadCoeffs(om_nir,om_nir_bet,om_nir_bet0,soil_alb_dir_nir,muav,Kopt,LAI,SAI,LAId);
Iup.dir_nir = H(1)/sgm + H(2) + H(3);
Idn.dir_nir = (H(4)/sgm)*exp(-Kopt*(LAI+SAI+LAId)) + H(5)*s1 + H(6)/s1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  // -DIFFUSE BEAM- radiation flux (== INFRARED range ==)
[H,sgm,s1]=ComputeRadCoeffs(om_nir,om_nir_bet,om_nir_bet0,soil_alb_dif_nir,muav,Kopt,LAI,SAI,LAId);
Iup.dif_nir = H(7) + H(8);
Idn.dif_nir = H(9)*s1 + H(10)/s1;
%%%%%%%%%%%%%%%%%%%%%
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[H,sgm,s1]=ComputeRadCoeffs(om,om_bet,om_bet0,soil_alb,muav,Kopt,LAI,SAI,LAId)
%%%%%%%%%%%%%%%%%%%%%%%%
%%% Computes coefficients required for estimation of upward/downward
%%% diffuse fluxes per unit incident direct beam and diffuse flux
%%%%%%%%%%%%%%
  b = 1 - om + om_bet;
  c = om_bet;
  d = om_bet0*muav*Kopt;
  f = muav*Kopt*(om - om_bet0);
  h = sqrt(b*b - c*c)/muav;
  sgm = muav*Kopt*muav*Kopt + c*c - b*b;
  u1 = b - c/soil_alb;
  u2 = b - c*soil_alb;
  u3 = f + c*soil_alb;
  s1 = exp(-h*(LAI+SAI+LAId));
  s2 = exp(-Kopt*(LAI+SAI+LAId));
  p1 = b + muav*h;
  p2 = b - muav*h;
  p3 = b + muav*Kopt;
  p4 = b - muav*Kopt;
  d1 = p1*(u1 - muav*h)/s1 - p2*(u1 + muav*h)*s1;
  d2 = (u2 + muav*h)/s1 - (u2 - muav*h)*s1;
  %H(0) = sgm; // <-- to make things easier...
  H(1) = -d*p4 - c*f;
  H(2) =  ((d-H(1)/sgm*p3)*(u1-muav*h)/s1 - p2*(d-c-H(1)/sgm*(u1+muav*Kopt))*s2)/d1;
  H(3) = -((d-H(1)/sgm*p3)*(u1+muav*h)*s1 - p1*(d-c-H(1)/sgm*(u1+muav*Kopt))*s2)/d1;
  H(4) = -f*p3 - c*d;
  H(5) = -(H(4)/sgm*(u2+muav*h)/s1 + (u3-H(4)/sgm*(u2-muav*Kopt))*s2)/d2;
  H(6) =  (H(4)/sgm*(u2-muav*h)*s1 + (u3-H(4)/sgm*(u2-muav*Kopt))*s2)/d2;
  H(7) =  c*(u1-muav*h)/(d1*s1);
  H(8) = -c*(u1+muav*h)*s1/d1;
  H(9) =  (u2+muav*h)/(d2*s1);
  H(10) = -s1*(u2-muav*h)/d2;
  %H(11) = s1; // <-- to make things easier...
end 
