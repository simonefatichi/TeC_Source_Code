%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Partition Precipitation    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Pr_sno,Pr_liq]=Precipitation_partition(Pr,Ta,Tmin,Tmax,ea,Pre)
%%%INPUTS
%%% Precipitation [mm/h]
%%% Air Temperature  [°C]
%%% OUTPUTS
%Pr_sno  solid precipitation [mm/h]
%Pr_liq  liquid precipitation [mm/h]
%%%%%%%%%%%%%
OPT_Pr_Part=2;
%%%%%
switch OPT_Pr_Part
    case 1
        %%REFERENCES %%  Wigmosta et al., 1994
        %Tmin= -1.1; Tmax = 3.3; %% General
        %Tmin = -1.1 ; Tmax = 2.5 ; %% RCW optimized
        if (Ta > Tmin) && (Ta < Tmax)
            Pr_sno = Pr*(Tmax - Ta)/(Tmax- Tmin);
            Pr_liq = Pr- Pr_sno;
        end
        if Ta <= Tmin
            Pr_sno = Pr;
            Pr_liq= 0;
        end
        if Ta >= Tmax
            Pr_sno = 0;
            Pr_liq= Pr;
        end
    case 2
        %%%%%%
        % Ding, et , 2014.   J. Hydrol  http://dx.doi.org/10.1016/j.jhydrol.2014.03.038.
        %%%%%%%%%%%%%
        esat=611*exp(17.27*Ta./(237.3+Ta)); %% [Pa] Vapor pressure saturation
        U=ea./esat; %% Relative Humidity
        Laten= 1000*(2501.3 - 2.361*(Ta)); %%%% Latent heat vaporization/condensaition [J/kg]
        cp=1005 + ((Ta +23.15).^2)/3364; %% specific heat air  [J/kg K]
        gam=cp.*100.*Pre./(0.622*Laten); %% [Pa/C] psycrometric constant
        del=(4098*esat)./((237.3+Ta).^2); %% Pa/C
        Twb = Ta - ( esat - ea )./( gam + del);    % [C] %% Wet bulb temperature
        %%%%%
        %%%%%%
        %Wetbulb temperature following Sadeghi et.al 2013 DOI: 10.1175/JTECH-D-12-00191.1
        %         if Ta > 0
        %             a = 611;  b = 17.368; c = 238.88; gam = 6.42*10^-4;
        %         else
        %             a = 611; b = 17.966; c = 247.15; gam = 5.68*10^-4;
        %         end
        %         xi = (-3*10^-7)*Ta.^3 - (10^-5)*Ta.^2 +(2*10^-5).*Ta + 4.44*10^-2;%%empirical coefficient
        %         phi = xi + gam*Pre/10;%%empirical coefficient
        %         lam = 0.0014*exp(0.027*Ta);%%empirical coefficient
        %         psi = a - gam*Pre/10*Ta-ea/1000;%%empirical coefficient
        %         Twb = (-phi+sqrt(phi^2-4*lam*psi))/(2*lam); %wetbulb temperature
        %         esat = a*exp(b*(Ta)/(Ta+c));% saturation vapor pressure Pa
        %         U=ea./esat; %% Relative Humidity
        %%%%%%%%%%%%%%%%
        
        g= 9.81; %% [m/s^2] gravity acceleration
        P_Ref= 1013.25; %% [Pa] reference pressure
        Rd =287.05; %% [J/kg K] dry air gas constant
        %%%%% Reference Elevation
        Zref= -((Ta+15)/2+273.15)*(Rd/g)*log(Pre./P_Ref); %% [m]
        Zref =Zref/1000; % elevation, [m] => [km]
        
        %%% Thresholds for discrimination precipitation
        dT = 0.215 - 0.099*U + 1.018*U^2;
        dS = 2.374 - 1.634*U;
        T0 = -5.87 - 0.1042*Zref + 0.0885*Zref^2 + 16.06*U - 9.614*U^2;
        %%%%%%%%%%
        Tmin = T0;
        Tmax = T0;
        %%%%%%%%%%%
        if dT/dS > log(2)
            Tmin = T0 - dS*log(exp(dT./dS) - 2*exp(-dT./dS));
            Tmax = 2*T0 - Tmin;
        end
        %%% Calculation of solid fraction of precipitation
        solid_fraction = 1./(1 + exp((Twb - T0)./dS));    % sleet
        %%% Discrimination of rain, sleet, and snow
        if (Twb > Tmin) && (Twb < Tmax)
            % sleet
            Pr_sno = Pr.*solid_fraction;
            Pr_liq = Pr.*(1-solid_fraction);
        end
        if Twb <= Tmin
            % snow
            Pr_sno = Pr;
            Pr_liq= 0;
        end
        if Twb >= Tmax
            % rain
            Pr_sno = 0;
            Pr_liq= Pr;
        end
end
return
