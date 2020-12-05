%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute Snow Balance       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[TsF,SWE,D,ros,In_SWE,SP_wc,WR_SP,U_SWE,NIn_SWE,dQ,Qfm,t_sls,Sm,dQ2]=Snowpacks(dt,...
    Ta,Ts,Tstm1,Tdew,Ws,t_slstm1,SWEtm1,Dtm1,rostm1,SP_wctm1,In_SWEtm1,In_max_SWE,dw_SNO,...
    Pr_liq,Pr_sno,ESN,ESN_In,Rn,H,QE,G,Qv,Csnow,Ccrown,Cwat,Cfol_H,fpr,Pr_sno_day,Th_Pr_sno,ros_max1,ros_max2)
%%%INPUTS
%dt time step [s]
dth = dt/3600; %% [h]
%Ta = %% air temperature [°C] --
%Ts =  ; %% surface temperature [°C] --- before the computation
%Tstm1 =  ; %% surface temperature [°C] previous step
%Tdew [°C] dew point temperature
%Ws; %% [m/s] wind speed  ---
%t_slstm1 % time since last snowfall [s]
%SWEtm1 [mm]
%Dtm1 [m] snow depth
%rostm1 [kg/m^3] density of snow
%SP_wctm1 -- Snow Pack Water Content [mm]
% Atm1 %% Snow Albedo
%In_SWEtm1
%In_max_SWE,
%Pr_liq, %% [mm/h] Rain Fallen
%Pr_sno, %% [mm/h] Snow Fallen
%ESN =   % Sublimation/Evaporation from Snow %% [mm/h]
%ESN_In % Sublimation/Evaporation from Intercepted Snow %% [mm/h]
%Rn
%H
%QE
%G
%Qv,
%Csno
%fpr = frequency of precipitation in dt [0-1]
%%% OUTPUTS
%TsF =  ; %% surface temperature [°C] --- After the computation
%In_SWE =
%U_SWE =
%SWE %% [mm]
%SP_wc -- Snow Pack Water Content [mm]
%WR_SP-- Water Released from Snow Pack [mm]
%dQ = %%% [W/m^2]  Variation Heat in the Snowpack
%Qf =  %%% [W/m^2]  %% Release Heat of Freezing in the snowpack
% A %% Snow Albedo
% D %% [m] snow depth
% ros [kg/m^3] density of snow
%t_sls % time since last snowfall [s]
%%%%%%%%%%%%%
%%REFERENCES (Bras, 1990) (Williams and Tarboton 1999)
%%% PARAMETERS
ci = 2093; %%% [J/Kg K] specific heat ice
row = 1000; % water density [kg/m^3]
Lf= 1000*333.7; %% [J/Kg] %% Latent heat melt/freezing
%L= 1000*(2501.3 - 2.361*(Ta)); %%% Latent heat vaporization/condensaition [J/kg]
%Ls= Lf +L; %%% Latent Heat of sublimation  [J/kg]
%%%%%%%%%%%%%
%Th_Pr_sno = 3.0; %%% [mm/h] Threshold Intensity of snow to consider a New SnowFall
t_sls = (t_slstm1 + dth)*(Pr_sno_day <= Th_Pr_sno) ; %% Time since last snowfall [h]
%%%%%%%%%%%%%%
Pr_sno= Pr_sno*dth; %% [mm]
Pr_liq= Pr_liq*dth; %% [mm]
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[U_SWE,In_SWE,NIn_SWE]=Snow_Interception(dth,Pr_sno,Tdew,Ws,t_sls,ESN_In,In_max_SWE,In_SWEtm1);
%%% U_SWE unpack swe
%%% In_SWE interepted after evaporation
%%%%%%%%%%% SNOWPACKENERGY AND MASS BALANCE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SWE = SWEtm1 + (Pr_sno*(1-(1-Csnow)*Cwat)- NIn_SWE + U_SWE) - ESN*dth;  %%% [mm] Snow Water Equivalent After Snow and Evaporation
%%%%
f_sp_wc_involved = 5./SWE ; f_sp_wc_involved(f_sp_wc_involved>1)=1;
SP_wc_L = f_sp_wc_involved*SP_wctm1*(Tstm1 >= -0.01); %%% Snowpack Liquid in the snowpack [mm]
SP_wc_F = f_sp_wc_involved*SP_wctm1*(Tstm1 < -0.01); %%% Snowpack Frozen in the snowpack [mm]
if SWE <= 0
    SWE = 0;
    Qfm = 0;
    dQ = Rn + Qv + Qfm - H - QE -G;   % [W/m^2]  Variation Heat in the Snowpack
    TsF = 0;
    Sm1 = 0; Sm2 =0 ;
    dQ2 = dQ; 
else
    %%%%%%%%%%
    SWEeb= max(0.1,SWE+In_SWE); %%%[mm] Minumum snow-pack layer to compute energy budget temperature variations  
    SWEeb= min(2000,SWEeb); %%%[mm]  Maximum Snow for computing energy budget (no more than 2m)
    %SWEeb=SWE+In_SWE;
    %%%% ENERGY BALANCE
    Qfm = (Ts < -0.01)*(row*(SP_wc_L/1000)*Lf/dt) - (Ts >= -0.01)*(row*(SP_wc_F/1000)*Lf/dt); %%% [W/m^2]  %% Freezing/Melting Water in the snowpack
    dQ = Rn + Qv + Qfm - H - QE -G;   % [W/m^2]  Variation Heat in the Snowpack
    dTs = (1000*dQ*dt)/(ci*(SWEeb)*row); %% Variation temperature in the snowpack [°C]
    TsF = Tstm1 + dTs ; %%  Temperature Snowpack [°C]
    Qcc =ci*(SWEeb)*TsF*row/1000; %%  Snowpack Cold Content [J °C/m^2 K]
    %%%%% MASS BALANCE -- TEMPERATURE UPDATE
    if Qcc <= 0;
        Sm1 =0 ; Sm2=0; % Snow melt [mm]
        dQ2 =0 ; 
    else
        TsF = 0;
        Sm = Qcc/Lf; %  Snow melt [Kg/m^2]
        Sm=  1000*(Sm/row); %% Snow melt [mm]
        r = SWE/(SWE+In_SWE); %% Ratio allocation snowmelt between SWE and In_SWE
        Sm1 = min(r*Sm,SWE);
        Sm2 = min((1-r)*Sm,In_SWE);
        SWE =SWE -Sm1;
        In_SWE =  In_SWE - Sm2;
        %%%%%%%%%%%%%%%%%%%%%%%%%
        dQ2= Qcc/dt -  Lf*(Sm1+Sm2)*row/(1000*dt) ; %% [W/m^2] Energy Lost -- theroetically to increase TsF or G
        %Qcc =0; % Snowpack Cold Content [J °C/m^2 K]
    end
end
%%%
Sm = Sm1 + Sm2 ; %%% Snow melt intercepted and pack layer
%%% SNOW DEPTH AND DENSITY
[D,ros]=Snow_Depth(dt,Sm1,row,Ts,Ta,Pr_sno,SWE,SWEtm1,Dtm1,rostm1,ros_max1,ros_max2); %%% Snow depth and Density
%%%%%%%% LIQUID WATER SNOWPACK
SP_wc_max = SWE*(0.03*(ros >= 200) + (ros < 200)*(0.03 +(0.1-0.03)*(200-ros)/200)); %% [] snow water retention factor [Belair 2003]
%SP_wc_max=0;
SP_wc =  SP_wctm1 + Sm + Pr_liq*(1 - (1-Csnow)*Cwat - sum(Ccrown.*Cfol_H)*(1-dw_SNO)); %% Snow Pack Water Content [mm]
WR_SP = (SP_wc - SP_wc_max)*(SP_wc >= SP_wc_max); %% Water Released from Snowpack [mm]
SP_wc =  SP_wc - WR_SP; %%[mm] %% Snow Pack Water Content [mm]
%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%% CANOPY SNOW INTERCEPTION
%%REFERENCES [Hedstrom and Pomeroy 1998; Pomeroy et al., 1998 --- Gelfan et al 2004 ---
%Pomeroy et al 2002 --- Liston and Elder 2006 ]
%%%% Parameter --- Hedstrom and Pomeroy 1998;
function[U_SWE,In_SWE,NIn_SWE]=Snow_Interception(dth,Pr_sno,Tdew,Ws,t_sls,ESN_In,In_max_SWE,In_SWEtm1)
Cp=1;  %%
%a=1.157e-006;%% [1/s]
%ah=a*3600;%% [1/h]
ah=4.1e-3;%% [1/h]
if In_max_SWE <= 0
    U_SWE = 0;
    In_SWE =0;
    NIn_SWE =0;
else
    NIn_SWE = 0.7*(In_max_SWE - In_SWEtm1)*(1- exp(-Cp*Pr_sno/In_max_SWE)); %%% New interepted snowpack  [mm]
    In_SWE = In_SWEtm1 + NIn_SWE -ESN_In*dth; %%% First updated Interception  [mm]
    %%%%%
    if (Tdew > 0 && Ws > 0.5 ) %%%  Gelfan et al 2004
        U_SWE = In_SWE; %% Unload Intercepted Snowpack [mm]
        In_SWE = 0 ; %% Intercepted snowpack [mm]
    else
        %U_SWE = (1-exp(-a*dth))*exp(-a*t_sls)*In_SWE; %%% Unload Intercepted Snowpack [mm]
        %U_SWE = -((1-dth*ah/2)/(1+ah*dth/2)-1)*In_SWE;
        U_SWE = (1-exp(-ah*dth))*In_SWE;
        In_SWE = In_SWE - U_SWE; % Intercepted snowpack [mm]
        if In_SWE <= 0.01  %% [mm]
            U_SWE = U_SWE + In_SWE ;
            In_SWE = 0 ;
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SNOW DEPTH AND DENSITY FUNCTION
% References (Bras, 1990)  [Essery et al., 1999]
%[Hedstrom and Pomeray 1998]  [Douville et al 1995]  [Belair et al 2003]
function[D,ros]=Snow_Depth(dt,Sm,row,Ts,Ta,Pr_sno,SWE,SWEtm1,Dtm1,rostm1,ros_max1,ros_max2)
%%% INPUTS
% dt = [s]
% Sm snow melting
%row % water density [kg/m^3]
%Ta  % air temperature [°C]
% Pr_sno [mm]
% SWEtm1 [mm]
% Dtm1 [m]
% rostm1 [kg/m^3] density of snow
%%% OUTPUTS
% D [m] snow depth
% ros [kg/m^3] snow density
ANS = 2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch ANS
    case 1 %%%%%%%%% Compaction due to snowfall and New Snow mixing
        %ros_max = 1000; %% [kg/m^3]
        ros_max = 450; %% [kg/m^3] 1000
        if (SWEtm1 > 0) || (Pr_sno>0)
            %%PARAMETERS
            %tau_f = 0.24; %%
            %tau_f = 0.0072; %%
            %tau_1 = 86400; %% [s]
            %ros_n = 67.92 + 51.25*exp(Ta/2.59); %%  [Hedstrom and Pomeray 1998]
            ros_n = 1000*(0.05 + (((9/5)*Ta +32)> 0).*(((9/5)*Ta+32)./100).^2); %% [kg/m^3] new snow density [from Bras 1990]
            if rostm1 == 0
                Dn = 0.001*SWE*(row/ros_n); %% New Snowpack [m]
            else
                Dn = 0.001*Pr_sno*(row/ros_n); %% New Snowpack [m]
            end
            if rostm1 > 0
                %%% Snow Compaction ??Temperature  %%%% Douville et al 1995
                %if Ts == 0
                %    ros =  ros_max -(ros_max-rostm1)*exp(-tau_f*dt/tau_1);
                %else
                ros =rostm1 ;
                %end
                %%% Snow Compaction ??Temperature
                D = 0.001*row*SWEtm1/ros; %% First Update Snow Depth [m]
                %%% Anderson and Crawford, (1964)
                del_D = 0.0254*( ((Pr_sno/25.4)*(1000*Dtm1/25.4))/(SWEtm1/25.4))*((1000*Dtm1/254)^0.35) ; %%% [m]  Snow compaction new snow
                %%% Approximation D = Dtm1 and SWEtm1 = SWE
                D_sm = 0.001*row*Sm/ros; %% [m] melting snow
                D = D - del_D + Dn - D_sm; %% [m] new snow depth
                ros = row*(SWE/(D*1000)); %% New density of snow [kg/m^3]
            else
                D= Dn ;% [m] new snow depth
                ros = ros_n; %% [kg/m^3]
            end
            %%%%%%%
        else
            D = 0;
            ros = 0;
        end
        %%% Final Check
        D(D<0) = 0;  ros(ros < 0) = 0;
        D(ros > ros_max) =  0.001*row*SWE/ros_max;
        ros(ros > ros_max)= ros_max;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    case 2 %%% Belair 2003
        %%%%%%%%%%%%%%%
        %ros_max1=520; %600; %%% [kg/m^3]
        %ros_max2=320; %450; %%% [kg/m^3]
        if Sm >= 0.1 %%[mm]
            ros_max =( 1000/row)*(ros_max1- 20.47*(1-exp(-Dtm1/0.0673))/Dtm1); % [kg/m^3]
        else
            ros_max =( 1000/row)*(ros_max2- 20.47*(1-exp(-Dtm1/0.0673))/Dtm1);% [kg/m^3]
        end
        if (SWEtm1 > 0) || (Pr_sno>0)
            %%%PARAMETERS
            tau_f = 0.24; %%
            tau_1 = 86400; %% [s]
            %ros_n = 67.92 + 51.25*exp(Ta/2.59); %%  [Hedstrom and Pomeray 1998]
            ros_n = 1000*(0.05 + (((9/5)*Ta +32)> 0).*(((9/5)*Ta+32)./100).^2); %% [kg/m^3] new snow density [from Bras 1990]
            %ros_n = (1000/row)*(109+ (6*Ta) +26*sqrt(ua)); %% Brun et al 1989
            if rostm1 == 0
                Dn = 0.001*SWE*(row/ros_n); %% New Snowpack [m]
            else
                Dn = 0.001*Pr_sno*(row/ros_n); %% New Snowpack [m]
            end
            if rostm1 > 0
                %%% Snow Compaction
                if rostm1 < ros_max
                    ros1 =  ros_max -(ros_max-rostm1)*exp(-tau_f*dt/tau_1);
                else
                    ros1 = rostm1 ;
                end
                %%% Snow mixing with new snow
                ros = (ros_n*Pr_sno + SWEtm1*ros1)/(Pr_sno+SWEtm1); %% [kg/m^3]
                D= 0.001*SWE*(row/ros); %% Snowpack [m]
                %%%%%%
            else
                D= Dn ;% [m] new snow depth
                ros = ros_n; %% [kg/m^3]
            end
            %%%%%%%
        else
            D = 0;
            ros = 0;
        end
        %%% Final Check
        D(D<0) = 0;  ros(ros < 0) = 0;
end
end


