%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute Soil Heat Flux     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[G,Tdp]=Soil_Heat(dt,Ts,Tstm1,Tdptm1,CTt)
%%%INPUTS
%dt time step [s]
%Ts =  ; %% surface temperature [°C]
%Tstm1 =  ; %% surface temperature [°C] previous step
%Tdptm1 = ; %% depth Temperature [°C] previous step
% CTt  [K m^2/J] Total Thermal Capacity Soil or Water 
%%% OUTPUTS
%Tdp
%G 
%%REFERENCES
%%% Noilhan and Mahfouf 1996 --- ISBA 2004
%%% Farouki (1981) --- Oleson et al., 2004 - Ivanov et al., 2008 --
%%% Boone et al., 2000 -- Johansen 1975
%%%%%%%%%% THERMAL PROPERTIES SOIL
tau= 86400; %% [s] time constant
%%%%%%%%%%%%%% COMPUTATION
dTs= Ts-Tstm1; %% Temperature Variation [°C]
%%% Noilhan & Planton (1989) Force Restore
%G = (1/CTt)*(2*pi*(Ts-Tdptm1)/tau + dTs/dt); %%  % Soil Heat Flux [W/m^2];
%Tdp = Tdptm1 + (dt/tau)*(Ts-Tdptm1); %% Depth Temperature [°C]
Tdp = (1/(1+dt/tau))*(Tdptm1 + (dt/tau)*Ts); %% Depth Temperature [°C]
G = (1/CTt)*(2*pi*(Ts-Tdp)/tau + dTs/dt); %%  % Soil Heat Flux [W/m^2];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
