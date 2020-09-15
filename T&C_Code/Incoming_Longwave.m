%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  %%%%%%%%%%%%
%   Subfunction Incoming_Longwave            %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Latm]=Incoming_Longwave(Ta,ea,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Juszak and Pellicciotti 2013; Flerchinger et al 2009 
if N<=1 %% N is cloudiness
    %%%%%%%%%%%%%
    % Ta air temperature [°C]
    % ea vapor pressure [Pa]
    % N Cloudiness [0-1]
    Ta_k = Ta +273.15; %% air temperature  [K]
    %ea=ea; %% vapor pressure [Pa]
    sigmaSB = 5.6704e-8; % Stefan-Boltzmann constant [W/m^2 K4]  %%
    %%%%%%%% Compute emissivity clear sky
    %e_cs = 0.70 + 5.95*(1e-5).*(ea./100).*(exp(1500./Ta_k)); %% [Idso 1981]
    %e_cs = 1.24*((ea/100)./Ta_k).^(1/7); %% [Brutsaert 1975, Bertoldi et al., 2004]
    %e_cs = 0.23+ 0.484*(ea./Ta_k).^(1/8);%% [Konzelmann et al 1994]
    w=4.65*ea./Ta_k; %%[kg/m2] precipitable water [Prata 1996]
    %w= exp(0.07*Tdew-0.075)*10;%%%% [kg/m2] Iqbal 1983
    %e_cs = 1 -(1+w).*exp(-sqrt(1.2+3.*w)); %%% [Prata 1996]
    e_cs = (59.38+113.7*(Ta_k/273.16).^6 +96.96*sqrt(w./25))./(sigmaSB.*Ta_k.^4); %% Dilley and O'Brien 1998
    %%% Compute attenuation cloud cover
    %K= (1+0.29*N);%% [Murshunova 1966]
    %K = (1+0.17*N.^2); %%%%  [TVA, 1972]
    %K = (1-N.^4) + (0.952*N.^4)./e_cs; %%% [Konzelmann et al 1994]
    K= (1-0.84.*N) + 0.84*N./e_cs; %%% [Unsworth and Monteith 1975]
    %%%%
    %f8 = -0.6732 + 6.24*1e-3*Ta_k - 9.14*1e-6*Ta_k.^2 ;
    %e8z = 0.24 + (2.98*1e-6*(ea./1000).^2)*exp(3000./Ta_k);
    %t8 = 1- e8z*(1.4-0.4*e8z);
    %K= 1 + (t8*N*f8)./e_cs; %%% [Kimball et al., 1982] 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    Latm = K.*e_cs.*sigmaSB.*(Ta_k).^4; %% LONG_WAVE RADIATION Incoming  [W/m^2]
else %% N is the Longwawe radiation incoming
    Latm=N;
end
end