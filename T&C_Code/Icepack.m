%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute Ice  Balance       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[TsF,ICE,ICE_D,IP_wc,WR_IP,dQ,Qfm,Im]=Icepack(dt,Ts,Tstm1,ICEtm1,IP_wctm1,...
    Pr_liq,EICE,Rn,H,QE,G,Qv,Cwat,Ccrown,Cfol_H,Csno,Cicew,Ice_wc_sp,WR_SP)
%%%INPUTS
%dt time step [s]
dth = dt/3600; %% [h]
%%% OUTPUTS
%%% PARAMETERS
ci = 2093; %%% [J/Kg K] specific heat ice
row = 1000; % water density [kg/m^3]
roi = 916.2; %% ice density [kg/m^3]
Lf= 1000*333.7; %% [J/Kg] %% Latent heat melt/freezing
%L= 1000*(2501.3 - 2.361*(Ta)); %%% Latent heat vaporization/condensaition [J/kg]
%Ls= Lf +L; %%% Latent Heat of sublimation  [J/kg]
%%% WR_SP Water Released from Snowpack [mm]
%%%%%%%%%%%%%
%%%%%%%%% Liquid Precipitation on Ice
Pr_liq= Pr_liq*dth; %% [mm]
%%%%%%%%%%%%%%%%%%%%%%%%%%
if Csno == 0
    %%%%%%%%%%% ICEPACK ENERGETIC AND MASS BALANCE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ICE = ICEtm1 - EICE*dth;  %%% [mm] Icepack After Sublimation and new
    %IP_wc_L = IP_wctm1*(Tstm1>=0);  %%% Liquid in the icepack [mm]
    %IP_wc_F = IP_wctm1*(Tstm1 < 0); %%% Frozen in the icepack [mm]
    if ICE <= 0
        ICE = 0;
        Qfm = 0;
        dQ = Rn + Qv + Qfm - H - QE -G ;   % [W/m^2]  Variation Heat in the Snowpack
        TsF = 0;
        Im = 0;
    else
        %%%%%%%%%%
        ICEeb= min(2000,ICE); %%%[mm]  Maximum Ice for computing energy budget (no more than 2m)
        %%%% ENERGY BALANCE
        %Qfm = (Ts < 0)*(row*(IP_wc_L/1000)*Lf/dt) - (Ts >= 0)*(row*(IP_wc_F/1000)*Lf/dt); %%% [W/m^2]  %% Freezing/Melting Water in the icepack (not working)
        Qfm=0; %%%
        dQ = Rn + Qv + Qfm - H - QE -G;   % [W/m^2]  Variation Heat in the Icepack
        dTs = (1000*dQ*dt)/(ci*(ICEeb)*row); %% Variation temperature in the Icepack [°C]
        TsF = Tstm1 + dTs ; %%  Temperature Icepack [°C]
        Qcc =ci*(ICEeb)*TsF*row/1000; %%  Icepack Cold Content [J °C/m^2 K]
        %%%%% MASS BALANCE -- TEMPERATURE UPDATE
        if Qcc <= 0;
            Im =0 ; % Ice melt [mm]
        else
            TsF = 0;
            Im = Qcc/Lf; %  Ice melt [Kg/m^2]
            Im=  1000*(Im/row); %% Ice melt [mm]
            Im = min(Im,ICE);
            ICE =ICE - Im;
        end
    end
else
    %ICE = ICEtm1;
    %Im=0;
    %dQ=NaN; Qfm =NaN; TsF=0;
    %%%%%%%%%%
    ICE=ICEtm1; 
    ICEeb= min(2000,ICE); %%%[mm]  Maximum Ice for computing energy budget (no more than 2m)
    Qfm=0; %%%
    dQ = -G;   % [W/m^2]  Variation Heat in the Icepack
    dTs = (1000*dQ*dt)/(ci*(ICEeb)*row); %% Variation temperature in the Icepack [°C]
    TsF = Tstm1 + dTs ; %%  Temperature Icepack [°C]
    Qcc =ci*(ICEeb)*TsF*row/1000; %%  Icepack Cold Content [J °C/m^2 K]
    %%%%% MASS BALANCE -- TEMPERATURE UPDATE
    if Qcc <= 0;
        Im =0 ; % Ice melt [mm]
    else
        TsF = 0;
        Im = Qcc/Lf; %  Ice melt [Kg/m^2]
        Im=  1000*(Im/row); %% Ice melt [mm]
        Im = min(Im,ICE);
        ICE =ICE - Im;
    end
    
end
%%%
%%% ICE-DEPTH AND DENSITY
ICE_D= 0.001*ICE*(row/roi); %% Icepack [m]
%%%%%%%% LIQUID WATER SNOWPACK
IP_wc_max = Ice_wc_sp*ICE ; %% [mm] Ice water retention factor
KT= ICE/1000; %% [h] %% Vertical Glacier Flow Time Constant  V = 1.0 [m/h] 
IP_outflow = min(IP_wctm1,(IP_wctm1/KT));  %%% [mm/h] 
IP_wc =  IP_wctm1 + Im + WR_SP - IP_outflow*dth + Pr_liq*(1-Csno)*(1 - (1-Cicew)*Cwat - sum(Ccrown.*Cfol_H));%% Ice Pack Water Content [mm]
WR_IP = (IP_wc - IP_wc_max)*(IP_wc >= IP_wc_max); %% Water Released from Icepack filling [mm]
IP_wc =  IP_wc - WR_IP; %%[mm] %% Snow Pack Water Content [mm]
WR_IP = IP_outflow*dth + WR_IP; %% Total Water Released from Icepack [mm]
%%%%%%%%%%%%%%%%%%%%%%%%%
end




