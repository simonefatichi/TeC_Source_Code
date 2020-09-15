%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Thermal Mod %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[dT]=SOIL_THERMAL_COND_FT(t,T,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
    Phy1,s_SVG,bVG,Osat,Ohy,Oicetm1,Otm1,dz,nz,Dz,G0,Gn,Tup,Tdown,AE,OPZ,OPT_FR_SOIL)
%%%%%%%%%%%%%%%%%%%%%%%%%%
T=reshape(T,1,length(T)); 
% T = % [K]
dT=zeros(nz,1); %[K s]
%%%%%%%%%%%%%%%%
g=9.81; %% m/s2
roi = 916.2; %% ice density [kg/m^3]
Lf= 1000*333.7; %% [J/Kg] %% Latent heat melt/freezing
%%%%%%%%%%%
[lanS,cv_Soil,~,~,~,Cw]=Soil_Thermal_properties_FT(T,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
    Phy1,s_SVG,bVG,Osat,Ohy,(Oicetm1+Otm1),OPT_FR_SOIL);
% %%% cv %  Volumetric heat capacity Soil  [J/m^3 K]
%capp = cv_Soil - roi*Lf*dOice_dT;%%% [J/m^3 K]
capp = cv_Soil + roi*Lf^2*(1000*Cw)./(g*(T+273.15));%% [J/m^3 K]
%capp = cv_Soil ;% %%% [J/m^3 K]
%capp=capp; 
capp(isnan(capp))=cv_Soil(isnan(capp));
%%%%%%
%%%%%
lanS_half(1:nz-1)= 0.5*(lanS(1:nz-1) + lanS(2:nz));
%%% Method of lines
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%OPZ 1 G0 Gn are inputs 
%%% Tup Tdown 
if OPZ == 2
    G0=-lanS_half(1).*((T(1)-Tup)./(0.001*dz(1)*0.5));  % [W/m^2]
    Gn=-lanS_half(nz-1).*((Tdown-T(nz))./(0.001*dz(nz)*0.5));  % [W/m^2]
end
%%% Tup Gn 
if OPZ == 3
    G0=-lanS_half(1).*((T(1)-Tup)./(0.001*dz(1)*0.5));  % [W/m^2]
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=-lanS_half.*((T(2:nz)-T(1:nz-1))./(0.001*Dz)); %% [W/m^2] %%% Flux positive downward from layer i+1 (above) to i (below)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% AE additional energy input to the layer 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% SOIL Thermal BALANCE
dT(1) = (1./(capp(1)*(0.001*dz(1))))*(G0 - q(1) + AE(1)) ;
for i =2:nz-1
    dT(i)=(1./(capp(i)*(0.001*dz(i))))*(q(i-1) - q(i) + AE(i)) ;
end
dT(nz) = (1./(capp(nz)*(0.001*dz(nz))))*(q(nz-1) - Gn + AE(nz)) ;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return