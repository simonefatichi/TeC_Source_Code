%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Heat_fluxes                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%REFERENCES %%   --- Noilhan and Mafhouf 1996 ---
% --- Bertoldi et al., 2004 --- Deardorff (1978) --- Bras (1990) 
%--- Ivanov et al., 2008  Oleson et al., 2004 
function[H,QE,T_H,EIn_H]=Heat_fluxes_VegSnow(dt,...
    Ta,Ts,ea,Pre,Csno,...
    dw_H,dw_SNO,Ccrown,FsunH,FshdH,...
    LAI_H,SAI_H,In_H,...
    ra,rs_sunH,rs_shdH,rb_H,Vavail_plant_H)
%%%INPUTS
%dt time step [s]
dth = dt/3600; %% [h] 
%Ta = %% air temperature [°C] -- 
%Ts =  ; %% surface temperature [°C] --- 
%ea = %% vapor pressure [Pa]  --- 
%Pre = %pressure [hPa - mbar] --- 
Pre = Pre*100; %pressure [Pa] --- 
% Csno  Snow Cover [0-1]
% Crock  Rock Cover [0-1]
% Curb  Urban Cover [0-1]
% Cwat  Water Cover [0-1]
% Cbare Bare soil Cover [0-1] 
%dw_H=   Wet vegetation Fraction FirstLayer 
%dw_L=  Wet vegetation Fraction Second Layer 
%dw_SNO = % Snow Cover Fraction First Layer --- [Similitude with rain interception ] 
%In_H =  Intercption in First Layer [mm]
%In_L  = Interceptio in Second Layer [mm] 
%In_urb = Interceptio in Urban Landscape [mm] 
%In_rock = Interceptio in Rocks [mm] 
%SWE = Snow Water Equivalent 
%In_SWE  = Intercepted Snow [mm] 
%EvL_Zs, [%] Evaporation Layer fraction [1...m]
%RfH_Zs, [%] Root Fraction for First Layer of Vegetation [1...m]
%RfL_Zs  [%] Root Fraction for Second Layer of Vegetation [1...m] 
%Pr_liq -- Liquid Precipitation  [mm/h] 
%Pr_sno -- Snow Precipitation [mm/h] 
%ra =  %%% [s/m]  Aerodynamic resistence  Heat flux 
%r_soil =  soil surface resistence [s/m] 
%rs_H  = [s/m] Stomatal resistence First Layer  
%rs_L  = [s/m] Stomatal resistence Second Layer 
%rap_H  = [s/m] Undercanopy resitence First Layer 
%rap_L =  [s/m] Undercanopy resistence Second Layer 
%%%
%r_soil = ra/r_soil - ra;  %%[s/m] because input r_soil is the beta factor actually 
%%%%%
%%% OUTPUTS 
%H  Sensible Heat  [W/m^2]
%QE Latent Heat [W/m^2]
%T_H= Transpiration First Layer Canopy   %% [mm/h]
%T_L= Transpiration Second Layer Canopy %% [mm/h]
%EIn_H=  Evaporation from Interception First Layer Canopy   %% [mm/h] 
%EIn_L= Evaporation from Interception Second Layer Canopy %% [mm/h]
%EG =  Evaporation from Bare Soil %% [mm/h]
%ESN =   % Sublimation/Evaporation from Snow %% [mm/h]
%ESN_In % Sublimation/Evaporation from Intercepted Snow %% [mm/h]
%EWAT = % Sublimation/Evaporation from water %% % [mm/h]
%EIn_urb = % Evaporation from Interception in urban Landscape [mm/h]
%EIn_rock  % Evaporation from Interception in rocks [mm/h]
%Qv = %% [W/m^2] Advected Heat Flux in Precipiation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% HIP. ONLY ONE TEMPERATURE %%%%%%%%%%%%%%%%%%%%%%%%%
TvHsun=Ts;  TvHshd=Ts;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TacH = Ta;  TacL = Ta; 
%%%%%%%%%%%%%%%%%
qTa=0.622*ea/(Pre-0.378*ea);% Specifc humidity Air  [] 

%%% PARAMETERS 
ro = (Pre/(287.04*(Ta+273.15)))*(1-(ea/Pre)*(1-0.622)); %% dry air density [kg/m^3]
row = 1000; % water density [kg/m^3]
cwv = 1858 + 0.382*(Ta) + 4.22e-04*(Ta)^2 - 1.996e-07*(Ta)^3 ; % [J/Kg K] specific heat water vapor  
cp=1005 + ((Ta +23.15)^2)/3364; %% specific heat air  [J/kg K] 
L= 1000*(2501.3 - 2.361*(Ta)); %%% Latent heat vaporization/condensaition [J/kg]
cair = cp*(1-qTa) + qTa.*cwv; cp=cair; 

%%%%%%%%%%%%%%
%dw_H=  min(1,(In_H./In_max_H).^(2/3)).*(In_max_H>0); %% Wet vegetation Fraction First Layer  --- [Deardorff (1978)] 
%dw_SNO =  min(1,(In_SWE/In_max_SWE)^(2/3))*(In_max_SWE>0); %% Snow Cover Fraction First Layer --- [Similitude with rain wet fraction ] 
%dw_SNO =  min(1,In_SWE/In_max_SWE)*(In_max_SWE>0); %% Snow Cover Fraction First Layer --- [Lee and Mahrt 2004]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% FlUXES COMPUTATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SENSIBLE HEAT 
%%%% DT = ro*cp*(Ts-Ta); %%% Sensible Heat Potential [J/m^3]
%%% Sensible Heat from H_Vegetation and Intercepted In_H
H_Hsun = (Csno)*ro*cp*(TvHsun-Ta)./(ra + rb_H./(2*(LAI_H+SAI_H).*FsunH*(1-dw_SNO))); %% [W/m^2]
H_Hshd = (Csno)*ro*cp*(TvHshd-Ta)./(ra + rb_H./(2*(LAI_H+SAI_H).*FshdH*(1-dw_SNO))); %% [W/m^2]
H_Hsun(isnan(H_Hsun))=0; 
H_Hshd(isnan(H_Hshd))=0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CELL -QUANTITY
H = sum(Ccrown.*H_Hsun)+ sum(Ccrown.*H_Hshd) ;  %[W/m^2] -- Sensible Heat
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% LATENT HEAT 
qTa=0.622*ea/(Pre-0.378*ea);% Specifc humidity Air  [] 
%%% esat_s=611*exp(17.27*Ts/(237.3+Ts)); %%[Pa] vapor pressure saturation per Ts 
%%% qTs_sat=(0.622*esat_s)/(Pre-0.378*esat_s); %% Specific humidity at esat(Ts)  []
%%% Dq = ro*(qTs_sat-qTa); %% Potential Evaporation from Wet Surface
%%% [kg/m^3] 
esat_TvHsun=611*exp(17.27*TvHsun/(237.3+TvHsun)); %%[Pa] vapor pressure saturation per Ts 
qTvHsun_sat=(0.622*esat_TvHsun)/(Pre-0.378*esat_TvHsun); %% Specific humidity at esat(Ts)  []
esat_TvHshd=611*exp(17.27*TvHshd/(237.3+TvHshd)); %%[Pa] vapor pressure saturation per Ts 
qTvHshd_sat=(0.622*esat_TvHshd)/(Pre-0.378*esat_TvHshd); %% Specific humidity at esat(Ts)  []
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Dealing with dew %%% 
if qTa > qTvHsun_sat %% Dew formation 
    dw_H = dw_H*0 + 1;
end 
%%%%%%%%%%%%%%%%%%%%%%%%% POSSIBLES EVAPORATION 
%%% Evaporation from Interception 
%%Epot;  %%% [kg/m^2.s] 
%%%%%  Transpiration
Tpot_Hsun = (Csno)*ro*(qTvHsun_sat-qTa)./(ra + rb_H./(LAI_H.*FsunH*(1-dw_SNO).*(1-dw_H))+ ...
    rs_sunH./(LAI_H.*FsunH*(1-dw_SNO).*(1-dw_H)) ); %%% [kg/m^2.s] 
Tpot_Hshd = (Csno)*ro*(qTvHshd_sat-qTa)./(ra + rb_H./(LAI_H.*FshdH*(1-dw_SNO).*(1-dw_H))+ ...
    rs_shdH./(LAI_H.*FshdH*(1-dw_SNO).*(1-dw_H)) ); %%% [kg/m^2.s] 
Tpot_H= Tpot_Hsun + Tpot_Hshd; 
Tpot_H(isnan(Tpot_H))=0; 
%%% Evaporation from Interception 
Epot_In_H = (Csno)*( ro*(qTvHsun_sat-qTa)./(ra + rb_H./(LAI_H.*FsunH.*dw_H) ) + ... 
     ro*(qTvHshd_sat-qTa)./(ra + rb_H./(LAI_H.*FshdH.*dw_H)) );  %%% [kg/m^2.s] 
Epot_In_H(isnan(Epot_In_H))=0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Including Limitation on the intercepted evaporation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CELL -QUANTITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_H = Ccrown.*Tpot_H;%% Transpiration First Layer Canopy   %% [kg/m^2.s] 
T_H= min(T_H,Vavail_plant_H); %%% % [kg/m^2.s] 
%%%
EIn_H= min(Ccrown.*Epot_In_H,(In_H/dth)*(row/3600/1000)); % Evaporation from Interception First Layer Canopy   %% [kg/m^2.s] 
%%%%%%%%%%%%%%%%%%%%%%%%%% 
QE = (sum(T_H) + sum(EIn_H))*L ; %[W/m^2] Latent Heat Flux 
%%%%%%%%%%%%%%%%%%%%%%%%%%
T_H= T_H*(1000*3600/row); %% Transpiration First Layer Canopy   %% [mm/h]
EIn_H= EIn_H*(1000*3600/row); % Evaporation from Interception First Layer Canopy   %% [mm/h] 
%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%
end 