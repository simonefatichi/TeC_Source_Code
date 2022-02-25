%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Management_Parameter     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Mpar]=Vegetation_Management_Parameter()
%%%%%%%%%%%%%%%%
%%%%%%%% Default Options there is no management 
%%%%%%%%%%% Vegetation 
%%% Grazing + Cut 
Mpar.jDay_cut=NaN; %%%   [0 0 0] Doy of grass cut or period of grazing 
Mpar.LAI_cut=0; %%%  % LAI of grass after cut  or grazed biomass  [gC/m2]
%%%% Fruit harvesting (crop - related) 
Mpar.jDay_harv=NaN;  %%  [0 0 0]  Doy of fruit harvesting 
Mpar.B_harv=0; %% Harvest Biomass  [gC/m2]  "-1" stays for completly harvested 
%%%% Forest Logging / selective logging / scarification / windthrown  
Mpar.Date_log = NaN; %%% Date of logging 
Mpar.fract_log = 0 ; %%% Fraction of aboveground biomass logged 
%%% Fire effects 
Mpar.Date_fire = NaN; 
Mpar.fire_eff = 0; 
Mpar.funb_nit= 0.15; %% Fraction of unburned Nitrogen 
%%%%%%%%%% Option for Girdling 
Mpar.Date_girdling = NaN; 
Mpar.fract_girdling = 0; 
%%%%%%%%%%%% Options for Crops 
Mpar.Date_sowing = NaN; 
Mpar.Date_harvesting = NaN; 
Mpar.Crop_B=[0 0]; %%% [gC m-2] 
Mpar.Crop_crown = 1.0; 
%%% Option Surviving belowground / resrprouting
Mpar.fract_resprout = 0.2; 
%%%%%%%% Fraction of material left in the field  
Mpar.fract_left=1; %% Fraction_left of leaves and dead leaves 
Mpar.fract_left_fr=0; %% Fraction left of fruits 
Mpar.fract_left_AB = 0; %% Fraction left of harvested/fired wood aboveground 
%Mpar.fract_left_AB = 0.8; %% Fraction left of harvested/fired wood aboveground (for moderate Fire)  
Mpar.fract_left_BG = 1; %% Fraction left of harvested/fired wood belowground 
end 
