%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Soil-Plant-Atmosphere-Continoum           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[CcF,An,rs,Rdark,F755nm,GAM,gsCO2]= photosynthesis_biochemical(Cc,IPAR,Csl,ra,rb,Ts,Pre,Ds,...
     Psi_L,Psi_sto_50,Psi_sto_00,...
    CT,Vmax,DS,Ha,FI,Oa,Do,a1,go,gmes,rjv)
%%%%%%%%%%%%%%%%%%%%%%%
%%% Leaf unit m^2 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% References 
%%%% References Collatz 1991- 1992 - Cox 2001 - Ivanov 2008  
%%% Dai et al 2004 Daly et al 2004 - Sellers 1996a,b  - Bonan et al 2011 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT
%%% Ci [umolCO2/molAIR] Leaf Interior CO2 concentration substomatal cavities 
%%% Cc [umolCO2/molAIR] Leaf Interior CO2 concentration mesophyll cells 
%%% IPAR [W/m^2] Intercepted Photosyntecially active Radiation
%%% Ca [ppm]-[umolCO2/molAIR] Atmospheric CO2 concentration
%%% ra =[s/m]
%%% rb =[s/m]
%%% Ts = Leaf temperature [°C]
%%% Ta  air temperature [°C]
%%% Pre = Atmospheric Pressure [mbar]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TP  --> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants   
%%% Vmax [umolCO2/ s m^2 ] Maximum Rubisco Capacity at 25°C 
%%% Ha =  [kJ/mol] Activation Energy of Vmax Plant Dependent 
%%% DS = [kJ / mol K]  entropy factor of Vmax Plant Dependent 
%%% FI Intrinsec quantum Efficiency [umolCO2/umolPhotons]
%%% Oa  Intercellular Partial Pressure Oxygen [umolO2/mol] -
%%% gsCO2  = [ umolCO2 / s m^2 ] Stomatal conductance to CO2 
%%% gmes =  Mesophyl Conductance [ molCO2 / s m^2 ] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ta=Ts; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUT
%%% CiF [umolCO2/molAIR] Leaf Interior  CO2 concentration
%%% An [gC/ s m^2]  Net Primary Productivity  -
%%% Cs %  [umolCO2/molAIR] Surface Leaf Concentration
%%% Rd % Dark Respiration [umolCO2/molAIR]
%%%%%%%%%%%%%%%%%%%%
Pre0 = 101325; %% [Pa] 
Tf = 273.15; %% [K] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%  Derivation of PAR Conversion to 
%lambda =0.55*(10^-6); %%[m] 0.4-0.7 [um]
%h = 6.626*10^-34; %% Planck Constant [J*s]
%Av = 6.02214*10^23; %% Avogadro Number
%c = 3*10^8; %% [m/s] Light Speed in Empty Space
%v = c/lambda; %% Frequency at Length lambda [1/s]
%%%% Energy of 1 Photon E = h*v
%E = h*v*Av*(10^-6); %%  [J/umolPhoton]  ---> 0.217
%b=1/E; ---> [umolPhoton/J]  4.59  [molPhoton/MJ]
%%%%%% Conversion Factors 
Pre = Pre*100; %% [Pa]
IPAR = IPAR*4.57 ; %%% [umolPhotons/s m^2] %%% [Dey 2004]
%%% Sellers et al 1996 
ra = ra*(0.0224*(Ta+273.15)*Pre0/(Tf*Pre))*10^-6; %% [m^2 s/umolH20]  %%% -> 
rb = rb*(0.0224*(Ta+273.15)*Pre0/(Tf*Pre))*10^-6; %% [m^2 s/umolH20]  %%% ->  
Cc = Cc*10^-6*Pre; %% [Pa] - Partial Pressure [Pa*molCO2/molAIR]
Oa = Oa*10^-6*Pre; %% [Pa]
Csl = Csl*10^-6*Pre; %% [Pa] -- Leaf surface CO2 concentration 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rmes = 1/(1e+6*gmes) ; %% [ s m^2 /umolCO2 ] Mesophyl Conductance 
go = go*10^6; %%%  [umolCO2 / s m^2] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% TEMPERATURE DEPENDENCE 
%%%% Maximum Rubisco Capacity Vm 
%%%%%%%%%% Temperature stress %%% Bonan et al., 2011; Kattge and Knorr 2007
Ts_k = Ts + 273.15; %%[K]
Tref = 25 + 273.15; %% [K] Reference Temperature
R =   0.008314; %%  [kJ·/ K mol] Gas Constant
%%%
ANS_TEMP=1; 
if  ANS_TEMP == 1 %% Kattge and Knorr 2007
    Hd =  200 ;% [kJ/mol]  Deactivation Energy -- constant
    %Ha = 72; [kJ/mol] Activation Energy - Plant Dependent
    %DS = 0.649; [kJ / mol K]  entropy factor - Plant Dependent
    %%%%%%% Mix of Temperature Function and High Temperature Inibition  %%%%
    kT =exp(Ha*(Ts_k-Tref)./(Tref*R*Ts_k)).*(1+exp((Tref*DS - Hd)/(Tref*R)))./(1+exp((Ts_k*DS-Hd)./(Ts_k*R)));
    Vm=Vmax*kT; %%% %% [umolCO2/ s m^2 ]
    %%%% Maximum Electron Transport Rate Jm
    %rjv= 1.97; %%[]
    Hd =  200;%  [kJ/mol]  Deactivation Energy -- constant
    Ha = 50;% [kJ/mol] Activation Energy
    DS = 0.646;% [kJ / mol K]  entropy factor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kT =exp(Ha*(Ts_k-Tref)./(Tref*R*Ts_k)).*(1+exp((Tref*DS - Hd)/(Tref*R)))./(1+exp((Ts_k*DS-Hd)./(Ts_k*R)));
    Jmax = Vmax*rjv;  %% [umol electrons/ s m^2 ]
    Jm= Jmax*kT; %%%% [umol electrons/ s m^2 ]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else %%% Bernacchi et al., 2001; 2003  Bonan et al., 2011;
    Hd = 149; %%  [kJ/mol]  Deactivation Energy -- constant
    Ha = 65.33; %%  [kJ/mol] Activation Energy - Plant Dependent
    DS = 0.485; %% [kJ / mol K]  entropy factor - Plant Dependent
    %%%%%%% Mix of Temperature Function and High Temperature Inibition  %%%%
    kT =exp(Ha*(Ts_k-Tref)./(Tref*R*Ts_k)).*(1+exp((Tref*DS - Hd)/(Tref*R)))./(1+exp((Ts_k*DS-Hd)./(Ts_k*R)));
    Vm=Vmax*kT; %%% %% [umolCO2/ s m^2 ]
    %%%% Maximum Electron Transport Rate Jm
    %rjv= 1.97; %%[]
    Hd = 152; %% [kJ/mol]  Deactivation Energy -- constant
    Ha = 43.5; %% [kJ/mol] Activation Energy
    DS = 0.495; %% [kJ / mol K]  entropy factor
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    kT =exp(Ha*(Ts_k-Tref)./(Tref*R*Ts_k)).*(1+exp((Tref*DS - Hd)/(Tref*R)))./(1+exp((Ts_k*DS-Hd)./(Ts_k*R)));
    Jmax = Vmax*rjv;  %% [umol electrons/ s m^2 ]
    Jm= Jmax*kT; %%%% [umol electrons/ s m^2 ]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%% Triose Phosphate Utilization
Ha = 53.1; %% [kJ/mol] Activation Energy
DS = 0.490; %%  [kJ / mol K]  entropy factor
Hd = 150.65; %% [kJ/mol]  Deactivation Energy
%%%%%%
TPU25=0.1182*Vmax; %% [umolCO2/ s m^2 ]
kT =exp(Ha*(Ts_k-Tref)./(Tref*R*Ts_k)).*(1+exp((Tref*DS - Hd)/(Tref*R)))./(1+exp((Ts_k*DS-Hd)./(Ts_k*R)));
TPU=TPU25*kT; %%  [umolCO2/ s m^2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if CT==4 
    %%%%% Temperature Dai et al 2004  - Sellers 1996b  - Bonan et al 2011 
    s1=0.3; %% [1/K] 
    s3=0.2; % [1/K] %% 0.3 Cox2001
    Tup = 40; %[°C]
    Tlow = 15; %%[°C]
    %%% Tup = 30-40 °C %%  Tlow = 5-15 °C (-inf)for C3  Sellers 1996b 
    f1T= 1/(1 +exp(s1*(Ts - Tup))); %%% Temperaure Function 1 for Maximum Rubisco Capacity
    f2T= 1/(1 +exp(s3*(Tlow - Ts)));%%% Temperaure Function 2 for Maximum Rubisco Capacity
    fT = 2.0^(0.1*(Ts-25)); 
    Vm = Vmax*fT*f1T*f2T; %% [umolCO2/ s m^2 ]
    %%%% 
    ke25 = 20000*Vmax; 
    ke = ke25*fT; 
end 

%%%% CO2 COMPENSATION POINT 
ANSG = 2; 
switch ANSG 
    case 0
        %%% Dai et al 2004- Cox 2001 Approach
        fT = 0.57^(0.1*(Ts-25));
        GAM = Oa/(2*2600*fT)*(CT==3); %% [Pa]
        %%%%%%%%%%%%%%%%%%%%%%%%%
    case 1
        %%%% Leuning Approach
        G0= 34.6;%%% 28  [umolCO2/molAIR]
        G1 = 0.0451; % 0.0509
        G2= 0.000347; % 0.001
        T0 = 293.2; %%[K]
        %%%%
        GAM = G0*(1 +G1*(Ts+273.15 -T0) + G2*(Ts+273.15-T0)^2); %%% [umolCO2/molAIR] -- CO2 Compensation point - Leuning 1995
        GAM = GAM*10^-6*Pre; %% [Pa] - Partial Pressure [Pa*molCO2/molAIR]
    case 2
        %%% Bonan et al 2011
        Ha = 37.83;   %% [kJ/mol] Activation Energy 
        kT= exp(Ha*(Ts_k-Tref)./(Tref*R*Ts_k));
        GAM25 = 42.75; %% %%[umol / mol]
        GAM25 = GAM25*10^-6*Pre; %%[Pa];
        GAM = GAM25*kT; %% [Pa] Michaelis-Menten Constant for C0_2
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


if CT == 3
    %%%%%%%%%%%%%%%%% Leuning 1995 - Dai et al 2004
    %%%%% Michaelis-Menten Constants for CO2 and O2
    %fT = 2.1^(0.1*(Ts-25));
    %Kc = 30*fT; %% [Pa] Michaelis-Menten Constant for C0_2
    %fT = 1.2^(0.1*(Ts-25));
    %Ko = 30000*fT; %% [Pa] Michaelis-Menten Constant for O_2
    %%%% Bonan et al 2011
    %%%%% Michaelis-Menten Constants for CO2 and O2
    Ha =79.43;  %% [kJ/mol] Activation Energy
    Kc25= 404.9; %%[umol / mol]
    Kc25= Kc25*10^-6*Pre; %%[Pa];
    kT= exp(Ha*(Ts_k-Tref)./(Tref*R*Ts_k));
    Kc = Kc25*kT; %% [Pa] Michaelis-Menten Constant for C0_2
    %%%
    Ha = 36.38;  %% [kJ/mol] Activation Energy
    Ko25 = 278.4 ; %%[mmol / mol]
    Ko25 = Ko25*10^-3*Pre; %%[Pa];
    kT= exp(Ha*(Ts_k-Tref)./(Tref*R*Ts_k));
    Ko = Ko25*kT; %% [Pa] Michaelis-Menten Constant for O_2
end

%%%%% Dark Respiration 
if CT == 3
    %%% Bonan et al 2011
    Ha = 46.39;
    DS = 0.490; 
    Hd = 150.65; 
    %%%%%%
    Rdark25 = 0.015*Vmax; 
    kT =exp(Ha*(Ts_k-Tref)./(Tref*R*Ts_k)).*(1+exp((Tref*DS - Hd)/(Tref*R)))./(1+exp((Ts_k*DS-Hd)./(Ts_k*R)));
    %%%
    Rdark = Rdark25*kT;
elseif CT ==4
    %%%%%%%%%% Dark Respiration Collatz 1991 - 1992 - Sellers 1996 Bonan et al 2011
    fT = 2.0^(0.1*(Ts-25));
    fT3= 1/(1 +exp(1.3*(Ts-55)));%%% Temperaure Function 3 for Respiration
    Rdark25 = 0.025*Vmax; 
    Rdark = Rdark25*fT*fT3; %%%  %%% [umolCO2/ s m^2 ] %% Leaf Maintenace Respiration / Dark Respiration
end


%%%%%%%%%%%%%%%%%%%%%%%% PHOTOSYNTHESIS FACTORS 
%%% Light - Electron Transport --> Farquhar and Wong 1984 - Leuning 1995 -
%%% Daly et al 2004 - Dai et al 2004  von Caemmerer 2000 
Q = FI*IPAR; %%  [umolCO2/ s m^2 ] %% Light Absorbed by Photosystem II in CO2 units 
d1= 0.7; d2= -(Q + Jm/4); d3= Q*Jm/4 ;  %% d1 = 0.95 Leuning 1995; d1 = 0.7 Bonan et al., 2011  
%%% ELectron Transport Rate 
J=min((-d2+sqrt(d2^2-4*d1*d3))/(2*d1),(-d2-sqrt(d2^2-4*d1*d3))/(2*d1)); % 
%J= min(roots([d1,d2,d3])); %% Smoothed Minimum between FI*IPAR and Jm  in [umolCO2/ s m^2 ]
if CT == 3
    %%% Rubisco Limited
    JC = Vm*(Cc -GAM)/(Cc + Kc*(1+Oa/Ko)); %%% Gross Assimilation Rate Limited by Rubisco% [umolCO2/ s m^2 ]
    %%% Light Limited
    JL = (J)*(Cc -GAM)/(Cc + 2*GAM); %%% Gross Assimilation Rate Limited by Light % [umolCO2/ s m^2 ]
    %%% Capacity of the leaf to export or utilize the products of photosynthesis
    JE = 3*TPU; %% Gross Assimilation Rate Limited by Export % [umolCO2/ s m^2 ]
elseif CT==4
    %%% Rubisco Limited
    JC = Vm; %%% Gross Assimilation Rate Limited by Rubisco% [umolCO2/ s m^2 ]
    %%% Light Limited
    JL = Q; %%% Gross Assimilation Rate Limited by Light % [umolCO2/ s m^2 ]
    %%% PEP Carboxylase Limited 
    JE = ke*Cc/Pre; 
    %JE =  20000*Vm*(Cc/Pre)*(CT==4); %% Gross Assimilation Rate Limited by Transport % [umolCO2/ s m^2 ]
end
%%%%%%%%%%%%%%%%%%%%%%%%%% Cox 2001 
%%% First Polynomium
if CT==3
    b1= 0.98; b2= -(JC+JL); b3= JC*JL ;  %% b1 =0.98  Collatz 1991  Selleers1996b; Bonan et al 2011 ;  b1 = 0.83 Cox 1998
elseif CT==4
    b1= 0.80; b2= -(JC+JL); b3= JC*JL ;  %%  b1 = 0.8 Bonan et al 2011;  for C4  0.83 Cox 2001
end
%JP= min(roots([b1,b2,b3])); %% Smoothed Minimum between JC and JE  [umolCO2/ s m^2 ]
JP=min((-b2+sqrt(b2^2-4*b1*b3))/(2*b1),(-b2-sqrt(b2^2-4*b1*b3))/(2*b1));
%%% Second Polynomium
if CT == 3
    c1 = 0.95; c2 = -(JP +JE); c3= JP*JE;  %%% c1 =0.95  Collatz 1991 Sellers1996b; Bonan et al 2011;   c1 =0.90 Cox1998
elseif CT == 4
    c1 = 0.95; c2 = -(JP +JE); c3= JP*JE;  %%% c1=0.95 Bonan et al 2011; for C4  0.93 Cox 2001
end
%A= min(roots([c1,c2,c3])); %%% %% Gross Assimilation Rate Potential % [umolCO2/ s m^2 ]
A=min((-c2+sqrt(c2^2-4*c1*c3))/(2*c1),(-c2-sqrt(c2^2-4*c1*c3))/(2*c1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% New Water Stress Function %%%%%
Rgsws=0.02;
p2= log((1 -Rgsws)/Rgsws)/(Psi_sto_00 - Psi_sto_50);%% [1/MPa]
q2=-p2*Psi_sto_50; %%[-]
Rgsw = 1./(1+exp(p2.*Psi_L+q2)); %% [fraction]
fO = (1-Rgsw) ;%%%%%%
fO(fO>1)=1;  fO(fO<0)=0; 
%%%%%%% Soil Moisture Stress
%%% Application to Vm-Jm or to Final A
%if (O > Owp) &&  (O < Oss)
%    fO = (O - Owp)/(Oss - Owp);
%else
%    if O <= Owp
%        fO=0;
%    end
%    if O >= Oss
%        fO=1;
%    end
%end
%%%%%%%%%%%%%%%%%%%% Solar-induced chlorophyll fluorescence (SIF)
%%% Lee et al 2015 GCB 
%%% Je is the actual electron transport rate calculated from the CO2 exchange data
if CT == 3
    Jfe = A*(Cc + 2*GAM)/(Cc -GAM);
elseif CT==4
    Jfe= A; %% [umolCO2/ s m^2 ]
end
%%%%%%%%%%%% 
fiP0= FI*4; %%% [umol Electrons/ umolPhotons]
fiP = fiP0*Jfe/Q;  %% [0.4 max - stress decrease ]  
%%%%%%%%%%%%%%%%%%%%
dls=1-fiP/fiP0; %% degree of light saturation 
%%%%%%%%%%
kf =0.05; 
kd = max(0.03*Ts+0.0773,0.087); 
kn = (6.2473*dls-0.5944)*dls; 
%kn = 5.01*((1+10)*dls.^1.93)./(10+dls.^1.93); 
%%%% 
fiF = kf/(kf+kd+kn)*(1-fiP);  % [umol Electrons/ umolPhotons]
SIF = IPAR*fiF; % %%% [umol electrons/s m^2]
%%% k theoretically a function of Vmax and Chlorophyll content 
k=0.0375*Vmax +8.25 ; %% [umol m-2 s-1 / W m-2 sr-1 um-1]  
F755nm =SIF/k; %% [W m-2 sr-1 um-1] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%% Cox 2001 
A = A*fO; %% Gross Assimilation Rate [umolCO2/ s m^2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
An = A - Rdark; % %% Net Assimilation Rate % [umolCO2/ s m^2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Ball-Woodrow-Berry --> Model  Dewar (2002) -- Correction Tuzet et al., 2003  
%gsCO2 = go + a1*An*Pre/((Cs-GAM)*(1+Ds/Do)); %%%  [umolCO2 / s m^2] -- Stomatal Conductance
gsCO2 = go + a1*An*Pre/((Cc-GAM)*(1+Ds/Do)); %%%  [umolCO2 / s m^2] -- Stomatal Conductance
gsCO2(gsCO2<go)=go; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rsCO2=1/gsCO2; %%% [ s m^2 / umolCO2 ] Stomatal resistence or Canopy 
%%%%%%%%%
%CiF = Cs - An*Pre/gsCO2; %%%%% [Pa] 
%CiF(CiF<0) = 0; 
CcF = Csl - An*Pre*(rsCO2 + rmes + 1.37*rb +ra); %%%%% [Pa] 
CcF(CcF<0) = 0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
rsH20 = (rsCO2/1.64)*(10^6); %%% [ s m^2 / molH20 ] Stomatal resistence or canopy 
An = (Csl - CcF)/(Pre*(rsCO2 + rmes + 1.37*rb + ra)); %%% Net Assimilation Rate  % [umolCO2/ s m^2 ]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CcF = CcF/(Pre*10^-6); %% [umolCO2 /molAIR ]
rs = rsH20*(Tf*Pre)/(0.0224*(Ts+273.15)*Pre0); %% [s/m]  Stomatal resistence or Canopy 
%An = An*12*(10^-6) ; %% [gC/ s m^2] Net Assimilation Rate  -
return
