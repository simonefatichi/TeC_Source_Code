%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Erosion                  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[er,ke]=Erosion(dt,Pr_liq,hc_H,hc_L,K_usle,Ccrown,Cbare,Csno,Cfol_H,Cfol_L,CLitter,Dr_H,Dr_L)
%%%% Computation of splash erosion
%%REFERENCES %%  USLE  (Wischmeier & Smith,  1978) 
%%%INPUTS
%dt time step [s]
dth = dt/3600; %% [h]
%Pr_liq ; %%% [mm/h]  Precipitation Liquid 
%K_usle %%% Soil Erodibility factor [kg*h/J*mm] 
%Ccrown Cfol_L Cfol_H 
%Cbare
%Csno 
%%%%%%%%%%%%
%%% OUTPUTS
%er= % [kg/h m^2] %% [Erosion] 
%ke= % kinetic energy [J /m^2 h]   
%%%%%%%%%%%%%%%%%%%%%%%%%%%% Throughfall and Dripping  [Brandt 1989, 1990] 
Pr_TH= Pr_liq*Cbare*(1-Csno) + (1-Csno).*sum((1-CLitter).*Pr_liq.*Ccrown.*(1-Cfol_H).*(1-Cfol_L)); %% Direct Throughfall [mm/h] 
Pr_DR_Hv= (1-Csno)*(1-CLitter).*(Dr_H.*(1-Cfol_L));  % Dripping from vegetation H [mm/h] 
Pr_DR_Lv= (1-Csno)*(1-CLitter).*(Dr_L);  % Dripping from vegetation L [mm/h] 
Pr_ero = Pr_TH + sum(Pr_DR_Hv) + sum(Pr_DR_Lv) ; % Erosive Precipitation [mm/h] 
%%%%% Specific Kinetic Energy 
%%ke= 0.119 + 0.0873*log10(Pr); %% Specific Kinetic Energy of Rainfall [MJ/ha mm] %%% Foster et al 1981 
ke_TH = 8.95 + 8.44*log10(Pr_TH); %%[J/m^2 mm] 
ke_LD_Hv = 15.8*sqrt(hc_H)-5.87; %% [J/m^2 mm] 
ke_LD_Lv = 15.8*sqrt(hc_L)-5.87; %% [J/m^2 mm] 
ke_TH(ke_TH<0)=0; 
ke_LD_Hv(ke_LD_Hv<0)=0;
ke_LD_Lv(ke_LD_Lv<0)=0;
ke = ke_TH*Pr_TH + sum(ke_LD_Hv.*Pr_DR_Hv) +  sum(ke_LD_Lv.*Pr_DR_Lv); %%% kinetic energy [J /m^2 h]   
%%%%%%%%%%%%%%%%[ Morgan et al 1998 ]
%%%%%%%%%%%%% COMPUTATION 
I_ero= Pr_ero; %% Effective Rainfall Intensity [mm/h] 
%R_usle= ke.*I_ero; %%   Rainfall Erosivity [J mm/ m^2 h^2] 
er = ke*I_ero*K_usle; %%% [kg/h m^2] %% [Erosion] 
return 
