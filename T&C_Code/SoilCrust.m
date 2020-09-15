%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Soil Crust                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[KsC,LC,OsatC,OhyC,PeC,E]=SoilCrust(dth,ke,Etm1,rsd,Osat,Ohy,L,Pe,Ks)
%%REFERENCES %% 
%%% Mualem et al., 1990 ; Assouline and Mualem  1997  
%%% Assouline, 2004; Mualem and Assouline 1989 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% sandy loam --> loam   loam-silt  
%%% rsd %%  density density dry soil [kg/m^3] 
%%% Pr_TH [mm/h] 
%%% %ke  kinetic energy [J /m^2 h]  
%dth= %% [h]
%%% Parameter 
rhoR = 2650; %%[kg/m^3] density  soild soil 
drhoF = 400;%% 327-451  [kg/m^3]
dcF = 10; % [mm] 12.3 - 3.5
z =   3500;  %% [mm^2/J] 3050 - 4020 
eta = 7000;  %%[mm^2/J] 6100 - 8400  
%%%%% E Kinetic Energy 
h = 0.25; %% [mm] depth 
E = Etm1 + ke*dth/(10^6) ; %%% kinetic energy [J /mm^2]  
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Soil Crust Formation 
drho = drhoF*(1-exp(-eta*E)); %%%% [kg/m^3]
dc=dcF*(1-exp(-z*E)); %% [mm] 
gam =-log(0.001)/dc; gam(gam<0)=0;%%% [1/mm] 
rhoC = rsd + drho*exp(-gam*h);  %% [kg/m^3]
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%
%%% Modification Soil-Hydraulic Properties 
%%% Osat [] Saturation moisture 0 kPa
%%% Ohy [] Hygroscopic Moisture Evaporation cessation
%%% L % Slope of logaritimc tension-moisture curve
%%% Pe % Tension at air antry (bubbling pressure) [kPa]
%%% Ks  % saturation conductivty [mm/h]
C=0.25/1000; %%% [m^3/kg]  0.25 - 1.37 
OsatC = Osat - (rhoC-rsd)/rhoR; 
OhyC = Ohy*(1 + (rhoC-rsd)/rsd); 
PeC = Pe*(1 + (rhoC-rsd)/rsd).^3.72;
LC= L - C*(rhoC - rsd); LC(LC<0)=0; 
KsC = Ks*(((OsatC-OhyC)/(Osat-Ohy)).^2.5).*(Pe./PeC).*(Pe./PeC).*(LC*(1+L)./(L*(1+LC))).*(LC*(1+L)./(L*(1+LC))); 
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Harmonic or geometric mean of Ks
%%% Mean of Osat Ohy Pe L
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% simple mean 
%KsC = (KsC + Ks)/2; 
%%% geom. mean 
%KsC =  exp((log(KsC) + log(Ks))/2);% [mm/h] Crust Saturation Cond.
%%% harmonic mean 
%KsC =  2/(1/KsC + 1/Ks);% [mm/h] Crust Saturation Cond.
LC = 0.5*(LC + L);
OsatC = 0.5*(Osat + OsatC);
OhyC = 0.5*(Ohy + OhyC);
PeC = 0.5*(Pe + PeC);
%%%%%%%%%%%%%%%%%%%%
end 


