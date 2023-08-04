%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute New Water Potential      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Psi_s]=Salinity_Plant(Salt,Tdp,Tl,Psi_s,Psi_sto_00,Psi_sto_50,Osm_reg_Max,eps_root_base)
%%%INPUTS
%%% Salt = Salt Concentration [mol Salt/ m-3] 
R = 8.314472;%% Gas constant universal [J /K mol]
iv = 2; %% van't Hoff coefficient for NaCl
Tdp=Tdp+273.15; %% [K]
Tl =Tl+273.15; %% [K]
%%%%%%%%% Soil of Water Salinity 
%eps_root_base=0.9; 
eps_root = eps_root_base + 1e-4*Salt; %% [-] Filtration efficiency for Mangrove (Perri et al 2018, WRR)  
%%% Osmotic Potential in the soil  
Osm_Pot =  -(eps_root).*Salt*R*iv*Tdp; %% [J m-3] 
Osm_Pot = Osm_Pot*1e-6; %% [MPa] Osmotic Potential 
%%% Osmotic Potential in the leaves 
%%% Theoretically concentration of salt in the xylem water 
Osm_Pot_L =  -(1-eps_root).*Salt*R*iv*Tl; %% [J m-3] 
Osm_Pot_L = Osm_Pot_L*1e-6; %% [MPa] Osmotic Potential Leaf 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Psi_L = Psi_s + Osm_Pot; 
%%%%%%%%%% Osmoregulation 
%Osm_reg_Max=1.0; %% [MPa] 
b_os_reg=1; 
Osm_Reg = Osm_reg_Max*(Psi_L./Psi_sto_50).^b_os_reg; 
Osm_Reg(Osm_Reg > Osm_reg_Max) = Osm_reg_Max; 
%%%%%%%%%%
%%%%%%%%%%%%
Psi_s= Psi_s + Osm_Pot + Osm_Reg; %% [MPa]
Psi_L=Psi_s; 
%%% Water stress 
% Rgsws=0.02;
% p2= log((1 -Rgsws)/Rgsws)/(Psi_sto_00 - Psi_sto_50);%% [1/MPa]
% q2=-p2*Psi_sto_50; %%[-]
% Rgsw = 1./(1+exp(p2.*Psi_L+q2)); %% [fraction]
% fO=1-Rgsw; 
% %%%% 
% add_fclo=0.2;
% Rgsw=1-(fO-add_fclo); 
% Psi_target = (log((1 -Rgsw)/Rgsw) - q2)/p2; 
%%%%%%%%%%%%
end 


%% Salt Units 
% NaCl = 1 mg/l   0.0172 mol/l   =17.2 mol/m3 
% 15*10-4 dS/m ;;  ECe (2-8) 
% 
% 1 dS/m - 1000 uS/cm  = 640 mg/l (common salt)  530 mg/l (NaCl) 
% 
% 10 mol/m3 = 1 dS/m  // 580 mg/l 
% 35 PSU – 35.16 g/kg  - 1.02 kg/l -  594  (600) mol/m3  - 53-59 dS/m      
% - 55'000 us/cm seawater  --  35 ppt 

