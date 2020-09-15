%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Infiltration              %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[f,fpot]=Infiltration_2(Osat,Ohy,L,alpVG,nVG,lVG,Pe,Ks_Zs,O33,Ks_mac,Omac,alpVGM,nVGM,lVGM,Phy,s_SVG,bVG,SPAR,O,Dz,WIS,cosalp,Pond)
%%REFERENCES %%  Richards Equation ponded surface
%%%INPUTS
%%% Osat [] Saturation moisture 0 kPa
%%% OF Water Content Soil  Active Layer for Infiltration
%%% L % Slope of logaritimc tension-moisture curve
%%% Pe % Tension at air antry (bubbling pressure) [kPa]
%%% Ks_Zs  % saturation conductivty [mm/h]
%%% WIS = [mm/h] Water Incoming to Soil Layer
%%% Dz [mm] Distance from surface to half-layer
%%% OUTPUTS
%f  Infiltration rate [mm/h]
%fpot Potential Infiltration rate [mm/h]
%%%%%%%%%%%%%
switch SPAR
    %    %%% PAE [ suction negative]
    case 1
        PAE = 0;%(1/alpVG);  %% Suction [mm]
    case 2
        gw=9810;
        PAE = -1000*1000*Pe/(gw); %%[mm]% Suction at air antry (bubbling pressure)
    case 3
        PAE = 0;%(1/alpVG);  %% Suction [mm]
    case 4
        PAE = 0;%(1/alpVG);  %% Suction [mm]
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Pond > 0
    P0= Pond; %%[mm]
else
    P0=0;
end
[K,P]=Conductivity_Suction(SPAR,Ks_Zs,Osat,Ohy,L,Pe,O33,alpVG,nVG,lVG,Ks_mac,Omac,alpVGM,nVGM,lVGM,Phy,s_SVG,bVG,O);
P=-P; %% [mm] Water Potential
Khalf= 0.5*(K + Ks_Zs);
fpot=Khalf.*(1*cosalp - (-PAE + (P - P0) )./Dz); %% [mm/h] %%
%fpot=Khalf.*(1*cosalp - (P - P0)./Dz); %% [mm/h] %%
f= max(0, min(fpot,WIS)); % %% Infiltrazione rate [mm/h]
%%%%%%%%%
end