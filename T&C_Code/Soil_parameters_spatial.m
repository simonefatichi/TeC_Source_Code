%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute soil parameters    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Osat,L,Pe,Ks,O33,rsd,lan_dry,lan_s,cv_s,K_usle]=Soil_parameters_spatial(Psan,Pcla,Porg)
%%%INPUTS
%Psan = [%]
%Pcla = [%] 
%Porg= [%]  
%Prock = [%]
DF = 1; %% density factor 
Rw=0.0; %%% Weight fraction of gravel [g Gravel /g bulk soil] 
%%% OUTPUTS 
%%% Osat [] Saturation moisture 0 kPa 
%%% L % Slope of logaritimc tension-moisture curve 
%%% Pe % Tension at air antry (bubbling pressure) [kPa]
%%% Ks  % saturation conductivty [mm/h]
%%% O33 %% 33 kPa Moisture 
%%% rsd %%  density density dry soil [kg/m^3] 
%%% lan_dry=  Thermal conductivity dry soil [W/m K] 
%%% lan_s Thermal conductivity soil solid [W/m K] 
%%% cv_s =  Volumetric heat capacity soil solid [J/m^3 K]
%%% K_usle %%% Soil Erodibility factor [ton*h/MJ*mm]  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Check 
Psil=1-Psan-Pcla-Porg; 
if Psil < 0
    disp('SOIL PERCENTAGE INPUTS INCONSISTENT')
    return 
end
Porg=Porg*100; 
%%REFERENCES %% [Saxton and Rawls 2006] 
%%%%%%%%%%%%%%%%% PARAMETERIZATION 
O1500t = -0.024*Psan + 0.487*Pcla + 0.006*Porg + ...
    0.005*(Psan.*Porg) -0.013*(Pcla.*Porg) + 0.068*(Psan.*Pcla)+0.031; %% 
O1500= O1500t + 0.14*O1500t -0.02; %% 1500 kPa Moisture 
O33t= -0.251*Psan + 0.195*Pcla + 0.011*Porg + ...
    0.006*(Psan.*Porg) -0.027*(Pcla.*Porg) + 0.452*(Psan.*Pcla)+0.299;
O33= O33t +(1.283*(O33t).^2 -0.374*(O33t)-0.015); %%33 kPa Moisture 
Os_33t= 0.278*Psan + 0.034*Pcla +0.022*Porg -0.018*(Psan.*Porg)...
    -0.027*(Pcla.*Porg) -0.584*(Psan.*Pcla) + 0.078; 
Os_33=Os_33t +(0.636*Os_33t -0.107);%%% SAT-33 kPa Moisture  
%%%%%%%%%%%%%%%%%%%%%%%% PARAMETERS 
B=(log(1500)-log(33))./(log(O33)-log(O1500)); %% Coefficient of moisture tension 
A= exp(log(33)+B.*log(O33)); % Coefficient of moisture tension 
L=1./B; %% Slope of logaritimc tension-moisture curve 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
Osat= O33 + Os_33 - 0.097*Psan +0.043; %[] Saturation moisture 0 kPa
rsd=(1-Osat).*2650;  %% normal density dry soil  [kg/m^3] 
%%%%%%%%%%% DENSITY EFFECTS 
rsd_df=rsd.*DF; % adj density density dry soil [kg/m^3] 
Osat_df = 1 - (rsd_df/2650); %%% %[] Saturation moisture 0 kPa
O33_df = O33 -0.2*(Osat- Osat_df); %%33 kPa Moisture 
Os_33_df =  Osat_df - O33_df; % SAT-33 kPa Moisture 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Osat=Osat_df;  O33=O33_df; Os_33=Os_33_df;  rsd=rsd_df; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ks = 1930*(Osat-O33).^(3-L); %% saturation conductivty [mm/h]
Pet=  -21.67*Psan -27.93*Pcla -81.97*Os_33 + ...
    71.12*(Psan.*Os_33) + 8.29*(Pcla.*Os_33) + 14.05*(Psan.*Pcla)+ 27.16;
Pe= Pet + (0.02*Pet.^2 -0.113.*Pet -0.70); %% % Tension at air antry (bubbling pressure) [kPa]
Pe(Pe<0.5)=0.5; %%[kPa] 
%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Gravel Effects
alpha= rsd./2650; % matric soil density / gravel density
Rv = (alpha.*Rw)./(1-Rw.*(1-alpha)); % Volume fraction of gravel [g/cm^3]
rhoB=(rsd/1000).*(1-Rv) + (Rv*2.65); %Bulk soil density (matric plus gravel), [g/cm^3]
Osatb = Osat.*(1-Rv);
O33b =  O33.*(1-Rv);
Kb = Ks.*(1-Rw)./(1-Rw.*(1-3*alpha/2)); % Saturated conductivity bulk soil [mm/h]
%%%%%% After gravel effects 
rsd = rhoB*1000; %%[kg/m^3]
Osat = Osatb;%%[-]
O33=O33b;%%[-]
Ks=Kb; %%[mm/h]
%%%%%%%
%Salinity effects neglected 
%%%%%%%%%%%%%%%
%Chk = (Osat >= O33) & (O33 >= O1500) ;
%if not (Chk) 
%    disp ('SOIL PARAMETERS INCONSISTENT') 
%end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% THERMAL CHARCTERISTIC OF SOIL
%%REFERENCES %% [Farouki 1981] [De Vries 1963] [Oleson et al., 2004]
%%[Lawrence and Slater, 2008]  [Lawrence et al., 2018]
%%% OLD 
% rsd_s= 2700*(1-Osat); %%% % normal density dry soil for thermal propert. comput. [kg/m^3] 
% lan_dry= (0.135*rsd_s + 64.7)./(2700 - 0.947*rsd_s); %% Thermal conductivity dry soil [W/m K] 
% lan_s= (8.8*Psan+ 2.92*Pcla)./(Psan+Pcla); %% Thermal conductivity soil solid [W/m K] 
% cv_s = 1e+6*(2.128*Psan +2.385*Pcla)./(Psan+Pcla); %% Volumetric heat capacity soil solid [J/m^3 K]
%%%%%%%%%%%%%%%%%%%%
Osat_min= 0.489-0.1267.*Psan; %[] Saturated Water Content Mineral Soil 
%rsd_min=(1-Osat_min)*2700;  %% normal density dry soil  [kg/m^3] 
%fom = 0; % (rsd_min*Porg*0.58*0.01)/130; %%% Organic Carbon matter fraction [-] 
fom = (Osat-Osat_min)./(0.9-Osat_min); %%
fom(fom<0)=0; fom(fom>1)=1; 
lan_s_om = 0.25; %% Thermal conductivity  Organic matter [W/m K]
lan_dry_om = 0.05; %% Thermal conductivity  Organic matter dry [W/m K]
cv_s_om = 2.5*1e+6; %%% %% Volumetric heat capacity Organic matter [J/m^3 K]
rsd_s= 2700.*(1-Osat); %%% % normal density dry soil for thermal propert. comput. [kg/m^3] 
lan_dry_min = (0.135*rsd_s + 64.7)./(2700 - 0.947*rsd_s); %% Thermal conductivity dry soil [W/m K] 
lan_s_min= (8.8*Psan+ 2.92*Pcla)./(Psan+Pcla); %% Thermal conductivity soil solid [W/m K] 
lan_dry = (1-fom).*lan_dry_min + fom.*lan_dry_om; 
lan_s = (1-fom).*lan_s_min+ fom.*lan_s_om; 
cv_s_min = 1e+6*(2.128*Psan +2.385*Pcla)./(Psan+Pcla); %% Volumetric heat capacity soil solid [J/m^3 K]
cv_s = (1-fom).*cv_s_min + fom.*cv_s_om; 
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% K USLE Parameter 
%%REFERENCES %% [Williams 1995]  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Porg_c = Porg/1.72; %% Organic Carbon 
fsand= (0.2 + 0.3*exp(-25.6*Psan.*(1-Psil))); 
fcli = (Psil./(Pcla + Psil)).^0.3; 
forg = (1- (0.25*Porg_c)./(Porg_c + exp(3.72-2.95*Porg_c))); 
fhisand = (1- (0.7*(1-Psan))./((1-Psan) + exp(-5.51-22.9*(1-Psan))));
K_usle = fsand.*fcli.*forg.*fhisand ; %%% [ton*h/MJ*mm]  K_Usle 
K_usle = K_usle/1000; %%% [kg*h/J*mm] erosivity factor 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return 





