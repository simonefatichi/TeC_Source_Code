%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Albedo_OtherSurface_Properties    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[sur_alb,e_sur]=Albedo_OS_Properties(Ccover,h_S,Par)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT 
%%% Ccover - Cwat Curb Crock %%% 
%%% h_S Solar Height 
%%% OUTPUT 
%sur_alb.dir_vis 
%sur_alb.dif_vis 
%sur_alb.dir_nir 
%sur_alb.dif_nir 
alb=NaN; 
%%%%%%%%%%%%%%%%%%%%%
if nargin < 3 
    Par = NaN; 
end 
%%ALBEDO  dry_vis dry_nir sat_vis sat_nir 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(Ccover,'Curb')
    alb = Par.alb; 
    e_sur = Par.e_sur; 
end 
if strcmp(Ccover,'Crock')
    alb = 0.25; 
    e_sur = 0.95; 
end 
if strcmp(Ccover,'Cdeb')
    alb = Par.alb; 
    e_sur = Par.e_sur ; 
end 
%%% Visible 
sur_alb.dir_vis = alb; 
sur_alb.dif_vis = sur_alb.dir_vis;
%%%% NIR 
sur_alb.dir_nir = alb; 
sur_alb.dif_nir = sur_alb.dir_nir;
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%
mu = sin(h_S); %mu is the cosine of the zenith angle of the incident beam
if strcmp(Ccover,'Cwat')
%%%%% Parameters from  Bonan et al., 1996 
    sur_alb.dir_vis = 0.06/(mu^(1.7) +0.15); 
    sur_alb.dir_nir = sur_alb.dir_vis; 
    sur_alb.dif_vis = 0.06; 
    sur_alb.dif_nir = 0.06; 
    e_sur = 0.96; 
end 
end 