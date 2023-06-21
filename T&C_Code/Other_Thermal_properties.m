%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute Other Surface Thermal Properties     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[lan,cv,CT]=Other_Thermal_properties(Cwat,Curb,Crock,ms)
%%REFERENCES
%%%INPUTS
%%% OUTPUTS
%Cv,  %  [K m^2/J] Total Thermal Capacity
%lan,  % [W/m K ] Thermal conductivity
%rho %     Density [kg/m^3]
%cs  %%% [J/kg K]  %% Specific Heat Soil
tau= 86400; %% [s] time constant
if Cwat == 1
    lan = 0.58; %%% [W/m K ] Thermal conductivity water
    rho = 1000; % [kg/m^3]
    cs =  4186 ; % [J/kg K] 
    cv =  rho*cs;  % [J/m^3 K]  Volumetric heat capcity water
    CT=2*(sqrt(pi/(lan*cs*rho*tau))); %%  [K m^2/J] Total Thermal Capacity
end
%%%
if Curb == 1
    lan = 0.92; %%% [W/m K ] Thermal conductivity
    rho = 2200; 
    cs =  700;
    cv =  rho*cs;  % [J/m^3 K] Volumetric heat capcity
    CT=2*(sqrt(pi/(lan*cs*rho*tau))); %%  [K m^2/J] Total Thermal Capacity
end
if Crock == 1
    lan = 3 ; % [W/m K ] Thermal conductivity
    rho = 2650; 
    cs =  920;  
    cv =  rho*cs;  % [J/m^3 K] Volumetric heat capcity
    CT=2*(sqrt(pi/(lan*cs*rho*tau))); %%  [K m^2/J] Total Thermal Capacity
end
lan=lan*ones(1,ms);
cv=cv*ones(1,ms); 
%% [W/m K ] Thermal conductivity
%Quartz mineral	3
%Rock, solid	2 - 7
%Rock, porous volcanic (Tuff)	0.5 - 2.5
%Asphalt	0.75
%Cement, portland	0.29
%Concrete, dense	1.0 - 1.8
%Concrete, stone	1.7
%Gravel	0.7
%Ice (0°C)	2.18
%Wood across the grain, white pine	0.12
%Wood across the grain, balsa	0.055
%Wood across the grain, yellow pine, timber	0.147
%% [kJ/kg K] Specific Heat Capacity 
%Asphalt	0.92
%Brick, common	0.9
%Cement dry	 1.55
%Concrete, stone		0.75
%Concrete, light	0.96
%Ice (0°C)		2.09
%Limestone	0.908
%Quartz mineral  0.8
%Rock salt	0.92
%Wood, oak		2
%Wood, white pine	2.5
%% [kg/m^3] Density 
%Asphalt 1700-2250 
%Concrete 2400 
%Rocks  2300-2900 
end