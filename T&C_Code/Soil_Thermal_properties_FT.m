%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute Soil Thermal Properties     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[lanS,cv_Soil,CTt,Oice,O,Cw]=Soil_Thermal_properties_FT(Tdp,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
    Phy1,s_SVG,bVG,Osat,Ohy,Otot,OPT_FR_SOIL)
%%REFERENCES
%%% Noilhan and Mahfouf 1996 --- ISBA 2004
%%% Farouki (1981) --- Oleson et al., 2004 - Ivanov et al., 2008 --
%%% Boone et al., 2000 -- Johansen 1975 Oleson et al., 2004
%%%INPUTS
%rsd, %%  density density dry soil [kg/m^3]
%lan_dry, %%   Thermal conductivity dry soil [W/m K]
%lan_s, %% Thermal conductivity soil solid [W/m K]
%cv_s %% Volumetric heat capacity soil solid [J/m^3 K]
%Oice
%Oice_max
%Osat
%O
%Owp; %% Residual/Wilting Point Water Content
%%% OUTPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
n=length(Otot);
lanS = zeros(1,n);  %  [W/m K ] Thermal conductivity Soil
cv_Soil=zeros(1,n);  %  Volumetric heat capacity Soil  [J/m^3 K]
rsoil=zeros(1,n);  %   Soil Density [kg/m^3]
cs_Soil=zeros(1,n);  %  %%% [J/kg K]  %% Specific Heat Soil
Oice = zeros(1,n);
Pre = Pre*100; %% [Pa]
%%%%
Tdp=reshape(Tdp,1,n);
%CTt,  %  [K m^2/J] Total Thermal Capacity Soil
%lanS,  % [W/m K ] Thermal conductivity Soil
%cv_Soil %  Volumetric heat capacity Soil  [J/m^3 K]
%%%%%%%%%% THERMAL PROPERTIES SOIL
row=1000; %%[kg/m^3]
roi = 916.2; %% ice density [kg/m^3]
Lf= 1000*333.7; %% [J/Kg] %% Latent heat melt/freezing
g=9.81; %% m/s2
lan_wat = 0.58; %%% [W/m K ] Thermal conductivity water
lan_ice = 2.29; %%%  [W/m K ] Thermal conductivity ice
lan_air = 0.023;  %%%  [W/m K ] Thermal conductivity air
cv_w =  4186000;  % [J/m^3 K] Volumetric heat capacity water
cv_i = 916.2*2117;  % [J/m^3 K] Volumetric heat capacity ice

esat_T=611*exp(17.27*Tdp./(237.3+Tdp)); %%[Pa] vapor pressure saturation per Ts
ro = (Pre./(287.04*(Tdp+273.15))).*(1-(esat_T/Pre)*(1-0.622)); %%  air density [kg/m^3]
cp=1005 + ((Tdp +23.15).^2)/3364; %% specific heat dry air  [J/kg K]
cwv = 1858 + 0.382*(Tdp) + 4.22e-04*(Tdp).^2 - 1.996e-07*(Tdp).^3 ; % [J/Kg K] specific heat water vapor  
qTdp_sat=(0.622*esat_T)./(Pre-0.378*esat_T); %% Specific humidity at esat(Ts)  [-] g wat vap/ gAIR
cv_air = ro.*(cp + qTdp_sat.*cwv); % [J/m^3 K] Volumetric heat capacity air



%%%%% O == All water in the soil
O=Otot;
O(O<Ohy)=Ohy(O<Ohy); O(O>Osat)=Osat(O>Osat);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = Tdp + 273.15; %% [K]

switch SPAR
    case 1
        %%% Van-Genuchten, 1980 Corrected
        Se = 1; mVG= 1-1./nVG;
        h = (1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
        %%%%%%%%%%%%%% VG - Dall’Amico, 2010)
        Tcrit = 273.15 + g*(0.001*h)*(273.15)/Lf;   %%[K]
        hfrz =1000*(0.001*h + Lf/(g*273.15).*(T-Tcrit).*(T<Tcrit)); %%[mm]
        Ol = Ohy + (Osat-Ohy)./((1+abs(alpVG.*(hfrz)).^nVG).^mVG);%%  [-] Unfrozen water
        Oice_max = O-Ol; Oice_max(Oice_max<0)=0;
        %%%%%%%%%%%%
        Se = (O-Ohy)./(Osat-Ohy);
        Phy = Phy1*101.9368; %% [mm]
        P1 = (1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
        P2 = -Phy*exp(Se.*-abs(bVG)); % [mm]
        h = (Se<s_SVG).*P2 +(Se>=s_SVG).*P1; % [mm]
        %%%%%%%%%%%%%% VG - Dall’Amico, 2010)
        Tcrit = 273.15 + g*(0.001*h)*(273.15)/Lf;   %%[K]
        hfrz =1000*(0.001*h + Lf/(g*273.15).*(T-Tcrit).*(T<Tcrit)); %%[mm]
        Ol = Ohy + (Osat-Ohy)./((1+abs(alpVG.*(hfrz)).^nVG).^mVG);%%  [-] Unfrozen water
        Ol(Ol>Osat)=Osat(Ol>Osat); Ol(Ol<Ohy)=Ohy(Ol<Ohy);
        Ol = min(Ol,O);
        Oice = O-Ol; Oice=min(Oice,Oice_max);
        Se_l = (Ol-Ohy)./(Osat-Ohy);
        Cw = (-alpVG.*mVG.*(Osat-Ohy)./(1-mVG)).*(Se_l.^(1./mVG)).*((1-(Se_l).^(1./mVG)).^mVG); %% [1/mm]
        Cw=Cw.*(T<273.15).*(Oice_max>0);
    case 2
        Tf=273.15;
        B=1./L;
        A= exp(log(33)+B.*log(O33)); % Coefficient of moisture tension
        Psi = zeros(1,n); Cw = zeros(1,n); Ol = zeros(1,n);
        %%%%%%%%%%%%%%%%%%%%
        Pb =-1e+06*Pe./(9810); %%[mm]
        %%% Saxton and Rawls 2006 +  corrected with Dall’Amico, (2010) for energy conservation
        h=Pb;
        Tcrit = 273.15 + g*(0.001*h)*(273.15)/Lf;   %%[K]
        hfrz =1000*(0.001*h + Lf/(g*273.15).*(T-Tcrit).*(T<Tcrit)); %%[mm]
        hfrz=-hfrz*9810/1e+06; %%% [kPa] Tension
        for i=1:n
            if hfrz(i) < 33
                Ol(i) = O33(i)+(33-hfrz(i))*(Osat(i)-O33(i))/(33-Pe(i));
                if hfrz(i) < Pe(i)
                    Ol(i)=Osat(i);
                end
            else
                Ol(i) =(hfrz(i)/A(i))^(-1/B(i));
            end
        end
        Oice_max = O-Ol; Oice_max(Oice_max<0)=0;
        %%%%%%%%%%%%
        for i=1:n
            if O(i) < O33(i)
                Psi(i) = A(i)*(O(i)^-B(i)); %% [kPa]
            else
                Psi(i) =33-((O(i)-O33(i))*(33-Pe(i))/(Osat(i)-O33(i)));%% [kPa]
            end
        end
        h =-1e+06*Psi/(9810); %%[mm]% 
        %%%%%%%%%
        Tcrit = 273.15 + g*(0.001*h)*(273.15)/Lf;   %%[K]
        hfrz =1000*(0.001*h + Lf/(g*273.15).*(T-Tcrit).*(T<Tcrit)); %%[mm]
        hfrz=-hfrz*9810/1e+06; %%% [kPa] Tension
        for i=1:n
            if hfrz(i) < 33
                Ol(i) = O33(i)+(33-hfrz(i))*(Osat(i)-O33(i))/(33-Pe(i));
                if hfrz(i) < Pe(i)
                    Ol(i)=Osat(i);
                end
            else
                Ol(i) =(hfrz(i)/A(i))^(-1/B(i));
            end
        end 
        Ol(Ol>Osat)=Osat(Ol>Osat); Ol(Ol<Ohy)=Ohy(Ol<Ohy);
        Ol = min(Ol,O);
        Oice = O-Ol; Oice=min(Oice,Oice_max);
        %%%%
        for i=1:n
            if Ol(i) < O33(i)
                Psi(i) = A(i)*(Ol(i)^-B(i)); %% [kPa]
            else
                Psi(i) =33-((Ol(i)-O33(i))*(33-Pe(i))/(Osat(i)-O33(i)));%% [kPa]
            end
            if Psi(i) < 33
                Cw(i) = -(Osat(i)-O33(i))/(33-Pe(i)); %% [1/kPa]
            else
                Cw(i) =-L(i)*(A(i)^L(i))*Psi(i)^(-L(i)-1); %% [1/kPa]
            end
        end
        Cw =-9810*Cw/(1e+06); %%[1/mm]
        Cw=Cw.*(T<Tf).*(Oice_max>0);
    case {3,4}
        %%% Van-Genuchten, 1980
        Se = 1; mVG= 1-1./nVG;
        h = (1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
        %%%%%%%%%%%%%% VG - Dall’Amico, 2010)
        Tcrit = 273.15 + g*(0.001*h)*(273.15)/Lf;   %%[K]
        hfrz =1000*(0.001*h + Lf/(g*273.15).*(T-Tcrit).*(T<Tcrit)); %%[mm]
        Ol = Ohy + (Osat-Ohy)./((1+abs(alpVG.*(hfrz)).^nVG).^mVG);%%  [-] Unfrozen water
        Oice_max = O-Ol; Oice_max(Oice_max<0)=0;
        %%%%%%%%%%%%%%%%%%%
        Se = (O-Ohy)./(Osat-Ohy);
        h = (1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
        Tcrit = 273.15 + g*(0.001*h)*(273.15)/Lf;   %%[K]
        hfrz =1000*(0.001*h + Lf/(g*273.15).*(T-Tcrit).*(T<Tcrit)); %%[mm]
        Ol = Ohy + (Osat-Ohy)./((1+abs(alpVG.*(hfrz)).^nVG).^mVG);%%  [-] Unfrozen water
        Ol(Ol>Osat)=Osat(Ol>Osat); Ol(Ol<Ohy)=Ohy(Ol<Ohy);
        Ol = min(Ol,O);
        Oice = O-Ol; Oice=min(Oice,Oice_max);
        Se_l = (Ol-Ohy)./(Osat-Ohy);
        Cw = (-alpVG.*mVG.*(Osat-Ohy)./(1-mVG)).*(Se_l.^(1./mVG)).*((1-(Se_l).^(1./mVG)).^mVG); %% [1/mm]
        Cw=Cw.*(T<273.15).*(Oice_max>0);
    case 5
        Tf=273.15;
        Pb =-1e+06*Pe./(9810); %%[mm]
        %%% Campbell 1974 corrected with Dall’Amico, (2010) for energy conservation
        h=Pb;
        Tcrit = 273.15 + g*(0.001*h)*(273.15)/Lf;   %%[K]
        hfrz =1000*(0.001*h + Lf/(g*273.15).*(T-Tcrit).*(T<Tcrit)); %%[mm]
        Ol = Osat.*(Pb./hfrz).^L; %%  %  [-] Unfrozen water
        Oice_max = O-Ol; Oice_max(Oice_max<0)=0;
        %%%% using Clapeyron + Clapp and Hornberger is non-energy conservative
        % Ol = Osat.*(Lf/g*(Pb*0.001).*(T-Tf)./Tf).^(-L); %% Unfrozen water (Max)
        %%%%%%%%%%%%
        h=Pb.*(O./Osat).^(-1./L); %%
        Tcrit = 273.15 + g*(0.001*h)*(273.15)/Lf;   %%[K]
        hfrz =1000*(0.001*h + Lf/(g*273.15).*(T-Tcrit).*(T<Tcrit)); %%[mm]
        Ol = Osat.*(Pb./hfrz).^L; %%  %  [-] Unfrozen water
        Ol(Ol>Osat)=Osat(Ol>Osat); Ol(Ol<Ohy)=Ohy(Ol<Ohy);
        Ol = min(Ol,O);
        Oice = O-Ol; Oice=min(Oice,Oice_max);
        %%%%%%%%% Brooks-Corey
        %Se = (Ol-Ohy)/(Osat-Ohy);
        %h=Pb.*(Se).^(-1./L); %% %% [mm]
        %Cw=((Se).^(3+(2./L)))./((-h./(L.*(Osat-Ohy))).*Se.^(2+1./L));  %%% [1/mm]
        %%%%
        h=Pb.*(Ol./Osat).^(-1./L); %%
        Cw=-Osat.*L.*((Pb./h).^L)./h; %%% [1/mm]
        Cw=Cw.*(T<Tf).*(Oice_max>0);
end
%% Oice_max Oice
%%%%% O == Only liquid water in the soil
O = Otot - Oice;
Oliq = 0*O;
%%%% Each Soil layer
for i=1:n
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    if Tdp(i) > 0
        lan_sat = (lan_wat^Osat(i))*(lan_s(i)^(1-Osat(i))); %% Saturated Conductivity [W/m K]
        Ke = log(O(i)/Osat(i))+ 1; %% Kersten number
        Ke= Ke*(Ke>=0);
    else
        %Oice ; %% max ice for saturated soil
        Oliq(i) = Osat(i) - Oice_max(i); %%% Liquid water content at saturation
        lan_sat = (lan_wat^Oliq(i))*(lan_s(i)^(1-Osat(i)))*(lan_ice^(Oice_max(i))); %% Saturated Conductivity [W/m K]
        Ke = (Oliq(i)+Oice(i))/Osat(i);
    end
    %%%%
    if (O(i)+Oice(i))/Osat(i) > 1*10^-7
        lanS(i) = Ke*lan_sat + (1-Ke)*lan_dry(i) ; % [W/m K ] Thermal conductivity Soil
    else
        lanS(i) = lan_dry(i);%   [W/m K ] Thermal conductivity Soil
    end
    %%%%%%%%%%%%%%%%%%%%%%
    %cv_Soil(i) = cv_s(i)*(1-Osat(i)) + O(i)*cv_w + Oice(i)*cv_i  + (Osat(i)-O(i)-Oice(i))*cv_air(i);  %  Volumetric heat capacity Soil  [J/m^3 K]
    %%% Volumetric Heat Capacity of Air in the soil is irrelevant
    cv_Soil(i) = cv_s(i)*(1-Osat(i)) + O(i)*cv_w + Oice(i)*cv_i  ;  %  Volumetric heat capacity Soil  [J/m^3 K]
    rsoil(i)= rsd(i) + (O(i)-Ohy(i))*row + (Oice(i)-Ohy(i))*roi; %% Soil Density [kg/m^3]
    cs_Soil(i) = cv_Soil(i)/rsoil(i) ; %%% [J/kg K]  %% Specific Heat Soil
end
%%%%%%%%%%%%%%%%%%%%%%%%%
tau= 86400; %% [s] time constant
CTt=2*(sqrt(pi/(lanS(1)*cs_Soil(1)*rsoil(1)*tau))); %%  [K m^2/J] Total Thermal Capacity Soil
%%%%%%%%%%%%%%%%
%%%% No freezing component
if OPT_FR_SOIL==0
    O=Otot;
    Oice=Otot*0;
end
end