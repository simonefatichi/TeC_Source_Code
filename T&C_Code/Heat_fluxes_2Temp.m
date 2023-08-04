%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Heat_fluxes                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%REFERENCES %%   --- Noilhan and Mafhouf 1996 ---
% --- Bertoldi et al., 2004 --- Deardorff (1978) --- Bras (1990) 
%--- Ivanov et al., 2008  Oleson et al., 2004 
function[H_H,H_L,QE_H,QE_L,Qv_H,Qv_L,T_H,T_L,EIn_H,EIn_L,EG,ESN,ESN_In,EWAT,ELitter,EICE,EIn_urb,EIn_rock]=Heat_fluxes_2Temp(dt,...
    Ta,Ts_H,Ts_L,ea,Pre,Csno,Crock,Curb,Cwat,Cbare,Cice,Cicew,Csnow,CLitter,Cdeb,...
    dw_L,dw_H,dw_SNO,Ccrown,FsunH,FshdH,...
    FsunL,FshdL,LAI_H,LAI_L,SAI_H,SAI_L,...
    In_H,In_L,In_urb,In_rock,SWE,In_SWE,...
    Pr_liq,Pr_sno,ra,rs_sunH,rs_sunL,rs_shdH,rs_shdL,rb_H,rb_L,rap_H,rap_L,r_litter,...
    r_soil,b_soil,alp_soil,Vavail,Vavail_plant_H,Vavail_plant_L,WAT_avail,ICE_avail,In_Litter,alp_litter) 
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
Pr_liq=Pr_liq/(1000*3600); %% [m/s]
Pr_sno=Pr_sno/(1000*3600); %% [m/s]
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
Twat=Ts_L; Trock=Ts_L;  Turb=Ts_L; Tg=Ts_L;
Tice=Ts_H;Tsno=Ts_H;
TvHsun=Ts_H;  TvHshd=Ts_H;
TvLsun=Ts_H; TvLshd=Ts_H;
%%%%%%%%%%%%%%%%%%%%%%%%%%%
qTa=0.622*ea/(Pre-0.378*ea);% Specifc humidity Air  [] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%TacH = Ta;  TacL = Ta; 
%%%%%%%%%%%%%%%%%
%%% PARAMETERS 
ci = 2093; %%% [J/Kg K] specific heat ice 
cw =  4186;  % [J/Kg K] specific heat water 
cwv = 1858 + 0.382*(Ta) + 4.22e-04*(Ta)^2 - 1.996e-07*(Ta)^3 ; % [J/Kg K] specific heat water vapor  
ro = (Pre/(287.04*(Ta+273.15)))*(1-(ea/Pre)*(1-0.622)); %% wet air density [kg/m^3]
row = 1000; % water density [kg/m^3]
cp=1005 + ((Ta +23.15)^2)/3364; %% specific heat dry air  [J/kg K]
cair = cp*(1-qTa) + qTa.*cwv; cp=cair; 

Lf= 1000*333.7; %% [J/Kg] %% Latent heat melt/freezing 
L= 1000*(2501.3 - 2.361*(Ta)); %%% Latent heat vaporization/condensaition [J/kg]
Ls= Lf +L; %%% Latent Heat of sublimation  [J/kg]
%%%%%%%%%%%%%%
%dw_H=  min(1,(In_H./In_max_H).^(2/3)).*(In_max_H>0); %% Wet vegetation Fraction First Layer  --- [Deardorff (1978)] 
%dw_L=  min(1,(In_L./In_max_L).^(2/3)).*(In_max_L>0); %% Wet vegetation Fraction Second Layer 
%dw_SNO =  min(1,(In_SWE/In_max_SWE)^(2/3))*(In_max_SWE>0); %% Snow Cover Fraction First Layer --- [Similitude with rain wet fraction ] 
%dw_SNO =  min(1,In_SWE/In_max_SWE)*(In_max_SWE>0); %% Snow Cover Fraction First Layer --- [Lee and Mahrt 2004]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% FlUXES COMPUTATION 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SENSIBLE HEAT 
%%%% DT = ro*cp*(Ts-Ta); %%% Sensible Heat Potential [J/m^3]
%%% Sensible Heat from H_Vegetation and Intercepted In_H
H_Hsun = (1-Cice)*(1-Csno)*ro*cp*(TvHsun-Ta)./(ra + rb_H./(2*(LAI_H+SAI_H).*FsunH*(1-dw_SNO))); %% [W/m^2]
H_Hshd = (1-Cice)*(1-Csno)*ro*cp*(TvHshd-Ta)./(ra + rb_H./(2*(LAI_H+SAI_H).*FshdH*(1-dw_SNO))); %% [W/m^2]
H_Hsun(isnan(H_Hsun))=0; 
H_Hshd(isnan(H_Hshd))=0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Sensible Heat from L_Vegetation and Intercepted In_L
H_Lsun = (1-Cice)*(1-Csno)*ro*cp*(TvLsun-Ta)./(ra +  rap_H + rb_L./(2*(LAI_L+SAI_L).*FsunL)); %% [W/m^2]
H_Lshd = (1-Cice)*(1-Csno)*ro*cp*(TvLshd-Ta)./(ra +  rap_H + rb_L./(2*(LAI_L+SAI_L).*FshdL)); %% [W/m^2]
H_Lsun(isnan(H_Lsun))=0; 
H_Lshd(isnan(H_Lshd))=0; 
%%% Sensible Heat from Soil
H_bare = (1-Cice)*(1-Csno)*ro*cp*(Tg-Ta)/(ra); %% [W/m^2]
H_ground = (1-Cice)*(1-Csno)*ro*cp*(Tg-Ta)./(ra + rap_H + rap_L); %% [W/m^2]
%%% Sensible Heat from Snow or Water or Urban or Rock 
H_sno_Undveg =  Csno*ro*cp*(Tsno-Ta)./(ra + rap_H); 
H_sno_free =  Csno*ro*cp*(Tsno-Ta)/(ra); 
H_wat =  (1-Cicew)*(1-Csnow)*ro*cp*(Twat-Ta)/(ra); 
H_ice =  Cice*(1-Csno)*ro*cp*(Tice-Ta)/(ra); 
H_urb =  (1-Cice*(1-Cdeb))*(1-Csno)*ro*cp*(Turb-Ta)/(ra); 
H_rock = (1-Cice)*(1-Csno)*ro*cp*(Trock-Ta)/(ra);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CELL -QUANTITY
% H = sum(Ccrown.*H_Hsun)+ sum(Ccrown.*H_Hshd) + ...
%     sum(Ccrown.*H_Lsun)+ sum(Ccrown.*H_Lshd) +...
%     sum(Ccrown.*H_ground)+ H_bare*Cbare + ...
%     H_wat*Cwat + H_urb*Curb + H_rock*Crock +  H_urb*Cdeb +...
%     H_ice*Cice*(1 - (1-Cicew)*Cwat)*(1-Cdeb) + ... 
%     (sum(Ccrown.*(LAI_H+SAI_H))*dw_SNO + (1 - sum(Ccrown)- (1-Csnow)*Cwat))*H_sno_free+...
%     sum(Ccrown.*H_sno_Undveg);  %[W/m^2] -- Sensible Heat
%%%%%%%%%%%%%%
H_H=sum(Ccrown.*H_Hsun)+ sum(Ccrown.*H_Hshd) + ...
    sum(Ccrown.*H_Lsun)+ sum(Ccrown.*H_Lshd) +...
    H_urb*Cdeb +...
    H_ice*Cice*(1 - (1-Cicew)*Cwat)*(1-Cdeb) + ... 
    (sum(Ccrown.*(LAI_H+SAI_H))*dw_SNO + (1 - sum(Ccrown)- (1-Csnow)*Cwat))*H_sno_free+...
    sum(Ccrown.*H_sno_Undveg);  %[W/m^2]
H_L= sum(Ccrown.*H_ground)+ H_bare*Cbare + ...
    H_wat*Cwat + H_urb*Curb + H_rock*Crock ;  %[W/m^2]
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% LATENT HEAT 
%%% esat_s=611*exp(17.27*Ts/(237.3+Ts)); %%[Pa] vapor pressure saturation per Ts 
%%% qTs_sat=(0.622*esat_s)/(Pre-0.378*esat_s); %% Specific humidity at esat(Ts)  []
%%% Dq = ro*(qTs_sat-qTa); %% Potential Evaporation from Wet Surface [kg/m^3] 
esat_Tg=611*exp(17.27*Tg/(237.3+Tg)); %%[Pa] vapor pressure saturation per Ts 
qTg_sat=(0.622*esat_Tg)/(Pre-0.378*esat_Tg); %% Specific humidity at esat(Ts)  []
esat_Twat=611*exp(17.27*Twat/(237.3+Twat)); %%[Pa] vapor pressure saturation per Ts 
qTwat_sat=(0.622*esat_Twat)/(Pre-0.378*esat_Twat); %% Specific humidity at esat(Ts)  []
esat_Tice=611*exp(17.27*Tice/(237.3+Tice)); %%[Pa] vapor pressure saturation per Ts 
qTice_sat=(0.622*esat_Tice)/(Pre-0.378*esat_Tice); %% Specific humidity at esat(Ts)  []
esat_Trock=611*exp(17.27*Trock/(237.3+Trock)); %%[Pa] vapor pressure saturation per Ts 
qTrock_sat=(0.622*esat_Trock)/(Pre-0.378*esat_Trock); %% Specific humidity at esat(Ts)  []
esat_Turb=611*exp(17.27*Turb/(237.3+Turb)); %%[Pa] vapor pressure saturation per Ts 
qTurb_sat=(0.622*esat_Turb)/(Pre-0.378*esat_Turb); %% Specific humidity at esat(Ts)  []
esat_Tsno=611*exp(17.27*Tsno/(237.3+Tsno)); %%[Pa] vapor pressure saturation per Ts 
qTsno_sat=(0.622*esat_Tsno)/(Pre-0.378*esat_Tsno); %% Specific humidity at esat(Ts)  []
esat_TvHsun=611*exp(17.27*TvHsun/(237.3+TvHsun)); %%[Pa] vapor pressure saturation per Ts 
qTvHsun_sat=(0.622*esat_TvHsun)/(Pre-0.378*esat_TvHsun); %% Specific humidity at esat(Ts)  []
esat_TvHshd=611*exp(17.27*TvHshd/(237.3+TvHshd)); %%[Pa] vapor pressure saturation per Ts 
qTvHshd_sat=(0.622*esat_TvHshd)/(Pre-0.378*esat_TvHshd); %% Specific humidity at esat(Ts)  []
esat_TvLsun=611*exp(17.27*TvLsun/(237.3+TvLsun)); %%[Pa] vapor pressure saturation per Ts 
qTvLsun_sat=(0.622*esat_TvLsun)/(Pre-0.378*esat_TvLsun); %% Specific humidity at esat(Ts)  []
esat_TvLshd=611*exp(17.27*TvLshd/(237.3+TvLshd)); %%[Pa] vapor pressure saturation per Ts 
qTvLshd_sat=(0.622*esat_TvLshd)/(Pre-0.378*esat_TvLshd); %% Specific humidity at esat(Ts)  []
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% Dealing with dew %%% 
if qTa > qTvHshd_sat %% Dew formation 
    dw_H = dw_H*0 + 1;
    dw_L = dw_L*0 + 1;
    r_soil = 0;
    alp_soil = 1; 
    b_soil=1; 
    r_litter = 0*r_litter; 
    alp_litter = 1; 
end 
%%%%%%%%%%%%%%%%%%%%
if qTa > qTg_sat %% Dew formation 
    r_soil = 0;
    alp_soil = 1; 
    b_soil=1; 
    r_litter = 0*r_litter; 
    alp_litter = 1; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%% POSSIBLES EVAPORATION 
%%% Evaporation from Interception 
%%Epot;  %%% [kg/m^2.s] 
%%%%%  Transpiration
Tpot_Hsun = (1-Cice)*(1-Csno)*ro*(qTvHsun_sat-qTa)./(ra + rb_H./(LAI_H.*FsunH*(1-dw_SNO).*(1-dw_H))+ ...
    rs_sunH./(LAI_H.*FsunH*(1-dw_SNO).*(1-dw_H)) ); %%% [kg/m^2.s] 
Tpot_Hshd = (1-Cice)*(1-Csno)*ro*(qTvHshd_sat-qTa)./(ra + rb_H./(LAI_H.*FshdH*(1-dw_SNO).*(1-dw_H))+ ...
    rs_shdH./(LAI_H.*FshdH*(1-dw_SNO).*(1-dw_H)) ); %%% [kg/m^2.s] 
Tpot_H= Tpot_Hsun + Tpot_Hshd; 
Tpot_H(isnan(Tpot_H))=0; 
Tpot_Lsun = (1-Cice)*(1-Csno)*ro*(qTvLsun_sat-qTa)./(ra + rb_L./(LAI_L.*FsunL.*(1-dw_L))+...
    rs_sunL./(LAI_L.*FsunL.*(1-dw_L)) + rap_H  ); %%% [kg/m^2.s] 
Tpot_Lshd = (1-Cice)*(1-Csno)*ro*(qTvLshd_sat-qTa)./(ra + rb_L./(LAI_L.*FshdL.*(1-dw_L))+...
    rs_shdL./(LAI_L.*FshdL.*(1-dw_L)) + rap_H ); %%% [kg/m^2.s] 
Tpot_L= Tpot_Lsun + Tpot_Lshd; 
Tpot_L(isnan(Tpot_L))=0; 
%%%%%  Evaporation Bare Soil 
EGpot_ground =  (1-Cice)*(1-Csno)*ro*b_soil*(alp_soil*qTg_sat-qTa)./(ra+r_soil+rap_H+rap_L+r_litter);  %%% [kg/m^2.s] 
EGpot_bare =    (1-Cice)*(1-Csno)*ro*b_soil*(alp_soil*qTg_sat-qTa)/(ra+r_soil);  %%% [kg/m^2.s] 
%%%%%% Evaporation from Litter 
Epot_Litter =  (1-Cice)*(1-Csno)*ro*(alp_litter*qTg_sat-qTa)./(ra+rap_H+rap_L+r_litter);  %%% [kg/m^2.s] 
%%%%%  Evaporation from Water 
Epot_wat = (1-Cicew)*(1-Csnow)*ro*(qTwat_sat-qTa)./ra; % [kg/m^2.s]
%%%%%  Evaporation from Ice
Epot_ice = (1-Csno)*ro*(qTice_sat-qTa)./ra; % [kg/m^2.s]
%%%%%  Evaporation from rock
Epot_rock = (1-Cice)*(1-Csno)*ro*(qTrock_sat-qTa)./ra; % [kg/m^2.s]
%%%%%  Evaporation from urban 
Epot_urb = (1-Cice*(1-Cdeb))*(1-Csno)*ro*(qTurb_sat-qTa)./ra; % [kg/m^2.s]
%%% Evaporation from Snow 
Epot_snow_free = Csno*ro*(qTsno_sat-qTa)/(ra); % [kg/m^2.s] 
Epot_snow_Undveg = Csno*ro*(qTsno_sat-qTa)./(ra + rap_H); % [kg/m^2.s] 
%%% Evaporation from Interception 
Epot_In_H = (1-Cice)*(1-Csno)*(ro*(qTvHsun_sat-qTa)./(ra + rb_H./(LAI_H.*FsunH.*dw_H) ) + ... 
     ro*(qTvHshd_sat-qTa)./(ra + rb_H./(LAI_H.*FshdH.*dw_H) ) );  %%% [kg/m^2.s] 
Epot_In_L = (1-Cice)*(1-Csno)*(ro*(qTvLsun_sat-qTa)./(ra + rb_L./(LAI_L.*FsunL.*dw_L) + rap_H ) + ... 
     ro*(qTvLshd_sat-qTa)./(ra + rb_L./(LAI_L.*FshdL.*dw_L) + rap_H ) );  %%% [kg/m^2.s] 
Epot_In_H(isnan(Epot_In_H))=0; 
Epot_In_L(isnan(Epot_In_L))=0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Including Limitation on the intercepted evaporation 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% CELL -QUANTITY
%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_H = Ccrown.*Tpot_H;%% Transpiration First Layer Canopy   %% [kg/m^2.s] 
T_L=  Ccrown.*Tpot_L; %% Transpiration Second Layer Canopy %% [kg/m^2.s] 
EG =  Cbare*EGpot_bare + sum(Ccrown.*EGpot_ground); % Evaporation from Bare Soil %% [kg/m^2.s] 
EG = min(EG,Vavail); %% % [kg/m^2.s] 
ELitter = min(sum(CLitter.*Ccrown.*Epot_Litter),In_Litter*(row/3600/1000)); % [kg/m^2.s] 
T_H= min(T_H,Vavail_plant_H); %%% % [kg/m^2.s] 
T_L= min(T_L,Vavail_plant_L); %%% % [kg/m^2.s]                   
%%%
EIn_H= min(Ccrown.*Epot_In_H,(In_H/dth)*(row/3600/1000)); % Evaporation from Interception First Layer Canopy   %% [kg/m^2.s] 
EIn_L= min(Ccrown.*Epot_In_L,(In_L/dth)*(row/3600/1000)); % Evaporation from Interception Second Layer Canopy %% [kg/m^2.s] 
EWAT = min(Cwat*Epot_wat,(WAT_avail/dth)*(row/3600/1000)); % Evaporation from water %% % [kg/m^2.s] 
EICE = min((1-Cdeb)*Cice*(1 - (1-Cicew)*Cwat)*Epot_ice,(ICE_avail/dth)*(row/3600/1000)); % Evaporation from water %% % [kg/m^2.s] 
EIn_urb = min((Curb+Cdeb*(1-Curb))*Epot_urb,(In_urb/dth)*(row/3600/1000));% Evaporation from Interception in urban Landscape  [kg/m^2.s] 
EIn_rock = min(Crock*Epot_rock,(In_rock/dth)*(row/3600/1000)); % Evaporation from Interception in rocks [kg/m^2.s] 
ESN_In = min((dw_SNO*sum(Ccrown.*(LAI_H+SAI_H)))*Epot_snow_free,(In_SWE/dth)*(row/3600/1000)); % Evaporation from Snow Intercepted %% [kg/m^2.s] 
ESN = min((1-(1-Csnow)*Cwat-sum(Ccrown))*Epot_snow_free + ...
    sum(Ccrown.*Epot_snow_Undveg),(SWE/dth)*(row/3600/1000)); % Evaporation from Snow %% [kg/m^2.s] 
%%%%%%%%%%%%%%%%%%%%%%%%%% 
% QE = (sum(T_H) + sum(T_L) + sum(EIn_H) + sum(EIn_L) + EG + EWAT + ELitter + EIn_urb +  EIn_rock)*L +...
%     (ESN + ESN_In + EICE )*Ls ; %[W/m^2] Latent Heat Flux
if Curb>0
    QE_H=(sum(T_H) + sum(T_L) + sum(EIn_H) + sum(EIn_L))*L +...
        (ESN + ESN_In + EICE )*Ls ; %[W/m^2] Latent Heat Flux
    QE_L=(EG + EWAT + ELitter + EIn_urb +  EIn_rock)*L ; %[W/m^2] Latent Heat Flux
else
    QE_H=(sum(T_H) + sum(T_L) + sum(EIn_H) + sum(EIn_L) + EIn_urb)*L +...
        (ESN + ESN_In + EICE )*Ls ; %[W/m^2] Latent Heat Flux
    QE_L=(EG + EWAT + ELitter +  EIn_rock)*L ; %[W/m^2] Latent Heat Flux
end
%%%%%%%%%%%%%%%%%%%%%%%%%%
T_H= T_H*(1000*3600/row); %% Transpiration First Layer Canopy   %% [mm/h]
T_L= T_L*(1000*3600/row); %% Transpiration Second Layer Canopy %% [mm/h]
EIn_H= EIn_H*(1000*3600/row); % Evaporation from Interception First Layer Canopy   %% [mm/h] 
EIn_L= EIn_L*(1000*3600/row); % Evaporation from Interception Second Layer Canopy %% [mm/h]
EG = EG*(1000*3600/row);  % Evaporation from Bare Soil %% [mm/h]
ELitter = ELitter*(1000*3600/row);  % Evaporation from Litter %% [mm/h]
ESN = ESN*(1000*3600/row);  % Evaporation from Snow %% [mm/h]
ESN_In = ESN_In*(1000*3600/row);  % Evaporation from  Intercepted Snow %% [mm/h]
EWAT = EWAT*(1000*3600/row); % Evaporation from water %% % [mm/h]
EICE = EICE*(1000*3600/row); % Evaporation from water %% % [mm/h]
EIn_urb = EIn_urb*(1000*3600/row);% Evaporation from Interception in urban Landscape [mm/h]
EIn_rock = EIn_rock*(1000*3600/row);% Evaporation from Interception in rocks [mm/h]
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Qv Advected Heat Flux in Precipiation 
Qv_H = cw*Pr_liq*row*(max(Ta,0)-Ts_H) + ci*Pr_sno*row*(min(Ta,0)-Ts_H)  ; %% [W/m^2] Advected Heat Flux in Precipiation 
Qv_L = cw*Pr_liq*row*(max(Ta,0)-Ts_L) + ci*Pr_sno*row*(min(Ta,0)-Ts_L)  ; %% [W/m^2] Advected Heat Flux in Precipiation 
%%%%%%%%%%
end 