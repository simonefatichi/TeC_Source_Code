%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Roughness                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[zom,zoh,disp_h,zom_H,zom_L,zoh_H,zoh_L,d_H,d_L,zom_other]=Roughness_New(D,ydepth,ICE_D,Cdeb,Deb_Par,Urb_Par,hc_H,hc_L,LAI_H,Ccrown_L,Cwat,Curb,Crock,Cice)
%%%%%%%%%%%%%%%%%%%%
%%REFERENCES %%   [Strack et al.,  2004]  --  [Brutsaert (1975)] --
%%% Mahat et al 2013
%%%INPUTS
% D snow depth [m]
% ydepth [m] Water ponding 
% ICE_D [m] ICE thickness 
% hc_H hc_L  Canopy height -- First and second layer [m]
% Cwat  
% Curb 
% Crock  
% Cice 
%%% OUTPUTS
%%%zom roughness eddy diffusivities for momentum [m]
%%%zoh roughness  eddy diffusivities for heat  [m]
%%% ROUGHNESS 
%%%%%%%% [Su 2002] new version of zom and zoh 
%%% [Wieringa 1992]  
zom_soil = 0.003; %% [m] bare soil roughness momentum 
zom_snow= 0.001; %% 0.001 %% [m] snow roughness momentum 
zom_ice = 0.001*(Cice==1); %% [m] ice roughenss momentum 
zom_wat = 0.0002*(Cwat == 1); %% %%0.0003 water  roughness momentum 
zom_urb = 0.123*Urb_Par.BuildH*(Curb > 0); %%% 0.3-2.5 [m] urban landscape roughness momentum 
zom_rock = 0.0003*(Crock == 1); %% [m] rock roughness momentum 
zom_debris = Deb_Par.zom*(Cdeb==1); %% [m] debris   roughness momentum 
%zom = hc*(0.23 - LAI^0.25/10 - (y-1)/67);
zom_H= 0.123*hc_H; %% vegetation roughness momentum [m] Brutsaert (1975)
zom_L =0.123*hc_L; %% vegetation roughness momentum [m] Brutsaert (1975)
%%%%%%%%%%%%%%%%%%
%%%% Displacement height 
%d = hc*(0.05 + LAI^0.02/2 + (y-1)/20);
d_L = 0.67*hc_L; 
d_H = 0.67*hc_H;
d_urban=0.67*(Urb_Par.BuildH)*(Curb); 
%%%%%%%%%%%%%%%%%
OPT_ROUGH = 0;
if OPT_ROUGH == 2
    %%%%%%% Computation of d and zom following Raupach 1992 and 1994 and
    %%%%%%% Verhoef et al 1997
    %%%% Only for overstory vegetation
    cd1 = 20.6;     %% ! taken from VEA1997
    cr  = 0.35;  %%    ! Value for trees, taken from VEA1997
    psicorr = 0.2;  %%  ! taken from VEA1997
    k =0.4; %% Von Karman constant
    csunder = 0.015; %% Vegetation
    cssoil = 0.003; %% Bare soil
    cs = csunder*Ccrown_L + cssoil*(1.0-Ccrown_L);
    frontala = LAI_H/2.0;
    gamraup = (cs+ cr*frontala).^(-0.5); gamraup(gamraup<3.33)=3.33;
    if LAI_H >0.01 
        d_H = hc_H.*(1 - (1-exp(-sqrt(cd1*LAI_H)))./sqrt(cd1*LAI_H));
    else
        d_H = 0.67*hc_H;
    end 
    zom_H = (hc_H - d_H)./(exp(k*gamraup - psicorr));
    %%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%
zom_L = (zom_L.*(max(0,1-D./hc_L)) + zom_snow*min(1,D./hc_L)).*(hc_L>0); %% [m] roughness comparison snow-vegetation  [Strack et al.,  2004]  
zom_H = (zom_H.*(max(0,1-D./hc_H)) + zom_snow*min(1,D./hc_H)).*(hc_H>0); %% [m] roughness comparison snow-vegetation  [Strack et al.,  2004]  
%%%%% 
zom_other =  zom_soil*(1-Cwat)*(1-Crock)*(1-Curb)*(1-Cice) + zom_wat + zom_urb*Curb + zom_rock*(1-Cice) + zom_ice ; % roughness eddy diffusivities for momentum [m]
zom_other =  zom_other*(1-Cdeb)+ zom_debris; %% % roughness [m] comparison Debris - Other 
if D > zom_other
    zom_other =  zom_snow; %%% [m] roughness comparison Other-Snow
end
%%%%%%
if ydepth > 0
    zom_wat = 0.0002;
    if (ydepth > zom_other) && (ydepth > D) && (ydepth > ICE_D) 
        zom_other = zom_wat;
    end
    zom_L = (zom_L.*(max(0,1-ydepth./hc_L)) + zom_wat*min(1,ydepth./hc_L)).*(hc_L>0); %% [m] roughness 
    zom_H = (zom_H.*(max(0,1-ydepth./hc_H)) + zom_wat*min(1,ydepth./hc_H)).*(hc_H>0); %% [m] roughness 
end
if Cice > 0 
     zom_L = (zom_L.*(max(0,1-ICE_D./hc_L)) + zom_ice*min(1,ICE_D./hc_L)).*(hc_L>0); %% [m] roughness 
     zom_H = (zom_H.*(max(0,1-ICE_D./hc_H)) + zom_ice*min(1,ICE_D./hc_H)).*(hc_H>0); %% [m] roughness 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Heat Roughness [m] 
zoh_L = zom_L*0.1; 
zoh_H=  zom_H*0.1; 
zoh_other= 0.1*zom_other; %% roughness  eddy diffusivities for heat  [m]  ???? [Brutsaert (1975)]
%%%%% PATCH SCALE ROUGHNESS 
zom= max(max(max(zom_H),max(zom_L)),zom_other); 
zoh= max(max(max(zoh_H),max(zoh_L)),zoh_other);
%%%%%%%%%%
disp_h = max(max(max(d_H),max(d_L)),d_urban);
return
%%%%%%%%%%%%%%%%%