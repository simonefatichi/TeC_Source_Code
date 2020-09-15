%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction LongwaveFluxesOtherSurface    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Labs]=LongwaveFluxesOS(Ts,Latm,SvF,e_sur)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Bonan 1996; Ivanov 2008; 
%%%%%%% INPUT 
% Latm Incoming LongWave Radiation [W/m^2] 
% Ts surface radiative temperature [°C]
% SvF sky view factor 
% e_sur surface  emissivity 
%%%%%%% OUTPUT
%%Labs  Long Wave absorbed by the surface [W/m^2] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a_sur=e_sur; %% surface absorptivity
Ts_k = Ts +273.15; %% surface temperature  [K]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%missivity 
%%%%%%%%%%%%%%%%%%%%%%%%%%
sigmaSB = 5.6704e-8; % Stefan-Boltzmann constant [W/m^2 K4]  %%
%%%%%%%%%%%%%%%%%%%%%%%%
Lsur = e_sur*sigmaSB*(Ts_k)^4; %% Long Wave from the surface  [W/m^2] 
Labs = SvF*(a_sur*Latm - Lsur); %%% Long Wave absorbed by the surface [W/m^2] 
end 

