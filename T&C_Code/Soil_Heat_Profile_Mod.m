%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute Soil Heat Flux     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[G0,T]=Soil_Heat_Profile_Mod(Ts,dt,Ttm1,dz,Dz,nz,lanS,cv_soil,Csno,Cice,G0)
%%%INPUTS
%dt = [s] time step 
%lanS =  %%  [W/m K ] Thermal conductivity Soil
%cv_soil  =  %%% [J/m^3 K] %  Volumetric heat capacity Soil
%dz  % [mm]  Thickness of the Layers
%nz % number of layers 
%Dz %
%%%%%% Boundary Mixed Condition
Tdown = NaN; %% [K]
%G0 = NaN; %% [W/m^2]
%%%%
if isnan(G0)
    opt=1; 
else
    opt=2; 
end 
%%%
if Csno == 1 || Cice == 1 
    Tup= 0;
else
    Tup =Ts; %% [°C]
end
%%%%
Gn = 0; %% [W/m^2]
%%%%%%%%%%
T_SPAN = [0 dt]; %% [s]
OPT_ST=odeset('AbsTol',0.05,'MaxStep',dt);%'NonNegative',ones(1,nz));
%%%%%%%%%%%%%%%%%%%
%%%%% Differential Equation Solving the Soil Moisture
%[Tout,Sout]=ode45(@SOIL_THERMAL_COND,T_SPAN,Stm1,OPT_ST,...
%    lanS',cv_soil',dz',nz,Dz',G0,Gn,Tup,Tdown,opt);
[tout,Tout]=ode45(@SOIL_THERMAL_COND,T_SPAN,Ttm1,OPT_ST,...
    lanS',cv_soil',dz',nz,Dz',G0,Gn,Tup,Tdown,opt);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T = Tout(end,:);   % Soil Layer Heat Content %% [J °C/m^2 K]
%%%%%%%%%%%%%%%%%%%%%%
if isnan(sum(T))
    %%%%%%%%%%%%%%%
    disp('NaN values in soil temperature')
    return
end
%T=S./(0.001*dz)./cv_soil; %% [°C]
%G=-lanS.*((T(2:nz)-T(1:nz-1))./(0.001*Dz));  % [W/m^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
G0=-lanS(1).*((T(1)-Tup)./(0.001*Dz(1)));  % [W/m^2]
%Tdown = T(nz) - Gn*(0.001*dz(nz)*0.5)/Par.lanS; %%[K]
 end
