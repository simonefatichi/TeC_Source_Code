%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Thermal Mod %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [dT]=SOIL_THERMAL_COND_MOD(t,T,lan,cv,dz,nz,Dz,G0,Gn,Tup,Tdown,opt)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Method of lines
%%% Soil Water
T=reshape(T,1,length(T)); 
% T = % [K]
dT=zeros(nz,1); %[K s]
%%%%%%%%%%%
capp = cv; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opt == 1 
        G0=-lan(1).*((T(1)-Tup)./(0.001*Dz(1)));  %  [W °C/m^2 K]
else
        G0=G0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
q=-lan(2:nz).*((T(2:nz)-T(1:nz-1))./(0.001*Dz(2:nz))); %% [W °C/m^2 K] %%% Flux positive downward from layer i (above) to i+1 (below)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% SOIL Thermal BALANCE
dT(1) = (1./(capp(1)*(0.001*dz(1))))*(G0 - q(1)) ;
for i =2:nz-1
    dT(i)=(1./(capp(i)*(0.001*dz(i))))*(q(i-1) - q(i)) ;
end
dT(nz) = (1./(capp(nz)*(0.001*dz(nz))))*(q(nz-1) - Gn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return 