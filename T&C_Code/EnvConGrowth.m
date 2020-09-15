%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  ENVIRONMENTAL CONTROL GROWTH %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% References [Fatichi et al 2014]
function[GF]= EnvConGrowth(Ta,Bfac,gcoef)
%%%% INPUT
%Ta [°C]  Hourly temperature
%Bfac [0-1] Stess factor 
%gcoef = 0.1304; % [gC/m2 h]
gcoef = gcoef/24; %% % [gC/m2 h]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% ENVIRONMENTAL LIMITATION ON GROWTH %%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Temp. Function
Ts_k = Ta + 273.15; %%[K]
Tref = 20 + 273.15; %% [K] Reference Temperature
Ha = 76; %% [kJ/mol] Activation Energy - Plant Dependent
R =   0.008314; %%  [kJ·/ K mol] Gas Constant
Hd = 285; %% [kJ/mol]  Deactivation Energy -- constant
DS = 0.933; %% [kJ / mol K]  entropy factor - Plant Dependent
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kT =exp(Ha*(Ts_k-Tref)./(Tref*R*Ts_k)).*(1+exp((Tref*DS - Hd)/(Tref*R)))./(1+exp((Ts_k*DS-Hd)./(Ts_k*R)));
kT(Ta<5)=0;
%%%%%%%
%%% Water Limitation Function
%%% Bfac 
%%%%%%
GF=gcoef.*Bfac.*kT;
GF = sum(GF); %% [gC/m2 day]
%%%%%%%%%%%%%%%%
end