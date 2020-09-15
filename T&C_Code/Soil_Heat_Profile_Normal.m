%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute Soil Heat Flux     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[G0,T,Gn]=Soil_Heat_Profile_Normal(Tup,dt,Ttm1,ms,Zs,lanS,cv,Tdown,G0,Gn,OPZ)
%%%INPUTS
%%%%%% Boundary Mixed Condition
%Tdown = NaN; %% [K]
%G0 = NaN; %% [W/m^2]
%%%%
% OPZ case 1
% %%%% G0 Gn
% %%%%%%%%%%%%%
% OPZ case 2
% %%% Tup Tdown
% %%%%%%%%%%%%%%
% OPZ  case 3
% %%% Tup Gn
% OPZ case 4 
%%%% G0 Tdown 
% %%%%%%%%%%%%%%
dz= diff(Zs); %%%% [mm]  Thickness of the Layers
%%%%
if (OPZ == 1) || (OPZ == 3)
    Tdown= Ttm1(ms) - Gn*(0.001*dz(ms)*0.5)/lanS(ms); %%[K]
end
if (OPZ == 1) || (OPZ == 4)
    Tup = Ttm1(1) + G0*(0.001*dz(1)*0.5)/lanS(1);%% [K]
end
%%%%%%%%%%
Ttm1= [ Tup Ttm1 Tdown];
%%%%%%%
nit=10; %%% Number of Internal time step
%G0p = zeros(1,nit); Gnp = zeros(1,nit);
%%%
for kk=1:nit
 
    [Tout] = Heat_CN_Normal_Function(Ttm1,dt/nit,lanS,cv,Zs,Tup,Tdown,Gn,G0,OPZ); 
    
    Ttm1 = Tout;
    
%      %%%%     %%%%%%%%%%%%%%%
%      if  (OPZ == 2) || (OPZ == 3)
%          G0p(kk)=-lanS(1).*((Ttm1(2)-Tup)./(0.001*dz(1)*0.5));  % [W/m^2]
%      end
%      if  (OPZ == 2) || (OPZ == 4)
%          Gnp(kk)=-lanS(ms)*((Tdown-Ttm1(ms+1))./(0.001*dz(ms)*0.5));
%      end
     
end
T=Tout(2:end-1);% Temperature Layer %% [°C ]
T=T';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
if isnan(sum(T))
    %%%%%%%%%%%%%%%
    disp('NaN values in Soil Temperature')
    return
end
%%%%%%%%%%%%%%%%
% %%%%%%%%%%%%%%%
if  (OPZ == 2) || (OPZ == 3)
    G0=-lanS(1).*((T(1)-Tup)./(0.001*dz(1)*0.5));  % [W/m^2]
    %G0=mean(G0p);
end
if  (OPZ == 2) || (OPZ == 4)
    Gn=-lanS(ms)*((Tdown-T(ms))./(0.001*dz(end)*0.5));
    %Gn=mean(Gnp);
end
%Tdown= T(nz) - Gn*(0.001*dz(nz)*0.5)/lanS(nz); %%[°C][K]
return