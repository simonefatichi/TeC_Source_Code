%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Albedo_Soil_Properties    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[soil_alb,e_gr]=Albedo_Soil_Properties(OS,Color_Class)
%%%%%%%%%%%%% Ref. Lawrence and Chase (2007); Oleson et al., 2010; 
%%% INPUT 
%%% OS [] Water content frist Layer  
%%% Color_Class [0-20] 
%%% OUTPUT 
%soil_alb.dir_vis 
%soil_alb.dif_vis 
%soil_alb.dir_nir 
%soil_alb.dif_nir 
%%%%%%%%%%%%%%%%%%%%%
%%ALBEDO  dry_vis dry_nir sat_vis sat_nir 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
e_gr= 0.96; % Ground emissivity
%%%%% Parameters from Oleson et al., 2010 (updated color class) 
if Color_Class > 0
    ALBEDO=[0.36 0.61 0.25 0.50
        0.34 0.57 0.23 0.46
        0.32 0.53 0.21 0.42
        0.31 0.51 0.20 0.40
        0.30 0.49 0.19 0.38
        0.29 0.48 0.18 0.36
        0.28 0.45 0.17 0.34
        0.27 0.43 0.16 0.32
        0.26 0.41 0.15 0.30
        0.25 0.39 0.14 0.28
        0.24 0.37 0.13 0.26
        0.23 0.35 0.12 0.24
        0.22 0.33 0.11 0.22
        0.20 0.31 0.10 0.20
        0.18 0.29 0.09 0.18
        0.16 0.27 0.08 0.16
        0.14 0.25 0.07 0.14
        0.12 0.23 0.06 0.12
        0.10 0.21 0.05 0.10
        0.08 0.16 0.04 0.08];
    %%%%%%%%%%%%%%%%%%
    soil_alb_dry_vis = ALBEDO(Color_Class,1);
    soil_alb_dry_nir = ALBEDO(Color_Class,2);
    soil_alb_sat_vis = ALBEDO(Color_Class,3);
    soil_alb_sat_nir = ALBEDO(Color_Class,4);
else %%% General Class for Unknow Color  
    soil_alb_dry_vis = 0.22;
    soil_alb_dry_nir = 0.45;
    soil_alb_sat_vis = 0.11;
    soil_alb_sat_nir = 0.225;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dAlbedo = 0.11 - 0.4*OS; dAlbedo(dAlbedo<0)=0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Visible 
soil_alb.dir_vis = soil_alb_sat_vis + dAlbedo; 
soil_alb.dir_vis(soil_alb.dir_vis>soil_alb_dry_vis)=soil_alb_dry_vis;  
soil_alb.dif_vis = soil_alb.dir_vis;
%%%% NIR 
soil_alb.dir_nir = soil_alb_sat_nir+ dAlbedo; 
soil_alb.dir_nir(soil_alb.dir_nir>soil_alb_dry_nir)=soil_alb_dry_nir;
soil_alb.dif_nir = soil_alb.dir_nir;
%%%%%%%%%%%%%%%%%%
end 


