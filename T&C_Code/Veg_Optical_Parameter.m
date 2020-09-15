%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Veg_Optical_Parameter        %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[PFT_opt]=Veg_Optical_Parameter(PFT_Class)
%%%%%%%%Ref. Dorman and Sellers (1989) Asner et al. (1998)
%%%% Possible PFT_Class 
i=PFT_Class; 
%NET Temperate -1  
%NET Boreal -2 
%NDT Boreal -3 
%BET Tropical -4 
%BET temperate -5 
%BDT tropical -6 
%BDT temperate -7 
%BDT boreal -8 
%BES temperate -9 
%BDS temperate -10
%BDS boreal -11
%C3 arctic grass -12 
%C3 grass -13
%C4 grass -14
%Crop1 -15
%Crop2 -16
%Oil Palm - 17 
%%%%%%%%%%%%%%%%
%%%%%%  Oleson et al., 2010;  update of vegetation optical parameters
if i>0
    OPTICAL_PAR_VEG =[ 0.01 0.07 0.35 0.16 0.39 0.05 0.10 0.001 0.001
        0.01 0.07 0.35 0.16 0.39 0.05 0.10 0.001 0.001
        0.01 0.07 0.35 0.16 0.39 0.05 0.10 0.001 0.001
        0.10 0.10 0.45 0.16 0.39 0.05 0.25 0.001 0.001
        0.10 0.10 0.45 0.16 0.39 0.05 0.25 0.001 0.001
        0.01 0.10 0.45 0.16 0.39 0.05 0.25 0.001 0.001
        0.25 0.10 0.45 0.16 0.39 0.05 0.25 0.001 0.001
        0.25 0.10 0.45 0.16 0.39 0.05 0.25 0.001 0.001
        0.01 0.07 0.35 0.16 0.39 0.05 0.10 0.001 0.001
        0.25 0.10 0.45 0.16 0.39 0.05 0.25 0.001 0.001
        0.25 0.10 0.45 0.16 0.39 0.05 0.25 0.001 0.001
        -0.30 0.11 0.35 0.31 0.53 0.05 0.34 0.120 0.250
        -0.30 0.11 0.35 0.31 0.53 0.05 0.34 0.120 0.250
        -0.30 0.11 0.35 0.31 0.53 0.05 0.34 0.120 0.250
        -0.30 0.11 0.35 0.31 0.53 0.05 0.34 0.120 0.250
        -0.30 0.11 0.35 0.31 0.53 0.05 0.34 0.120 0.250
        -0.39 0.09 0.45 0.16 0.39 0.05 0.25 0.001 0.001];
else
    OPTICAL_PAR_VEG =NaN*ones(1,9);
    i=1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
PFT_opt.chiL=OPTICAL_PAR_VEG(i,1); %% Leaf Angle Distribution Parameter Ross, 1975
PFT_opt.alf_lf_vis=OPTICAL_PAR_VEG(i,2);  %reflectances
PFT_opt.alf_lf_nir=OPTICAL_PAR_VEG(i,3);
PFT_opt.alf_st_vis=OPTICAL_PAR_VEG(i,4);
PFT_opt.alf_st_nir=OPTICAL_PAR_VEG(i,5);
PFT_opt.alf_ld_vis=OPTICAL_PAR_VEG(i,4);
PFT_opt.alf_ld_nir=OPTICAL_PAR_VEG(i,5);
PFT_opt.tau_lf_vis=OPTICAL_PAR_VEG(i,6);  %transmittances
PFT_opt.tau_lf_nir=OPTICAL_PAR_VEG(i,7);
PFT_opt.tau_st_vis=OPTICAL_PAR_VEG(i,8);
PFT_opt.tau_st_nir=OPTICAL_PAR_VEG(i,9);
PFT_opt.tau_ld_vis=OPTICAL_PAR_VEG(i,8);
PFT_opt.tau_ld_nir=OPTICAL_PAR_VEG(i,9);
return
