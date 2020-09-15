%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Albedo_OtherSurface_Properties    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[sur_alb,e_sur]=Albedo_Ice_Properties(alb_ice)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT 
%%% Ccover - Cwat Curb Crock %%% 
%%% h_S Solar Height 
%%% OUTPUT 
%sur_alb.dir_vis 
%sur_alb.dif_vis 
%sur_alb.dir_nir 
%sur_alb.dif_nir 
%%%%%%%%%%%%%%%%%%%%%
%%ALBEDO  dry_vis dry_nir sat_vis sat_nir 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
alb = alb_ice;
e_sur = 0.97;
%%% Visible 
sur_alb.dir_vis = alb; 
sur_alb.dif_vis = sur_alb.dir_vis;
%%%% NIR 
sur_alb.dir_nir = alb; 
sur_alb.dif_nir = sur_alb.dir_nir;
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
end 