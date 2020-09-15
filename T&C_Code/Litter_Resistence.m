%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Litter Resistence           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%REFERENCES %%
function[r_litter,alp_litter]=Litter_Resistence(Ws,Ta,Pre,zatm,disp_h,zom,Sllit,BLit,In_Litter)
%%%INPUTS
%%% Blit = [kg/m2] Litter Biomass
%%% In_Litter =[mm] Interception in litter
%%% Air temperature [°C]
%%% Atm. Pressure [hPa]
%%% OUTPUTS
%r_litter=  [s/m] Litter surface resistance
Ts_k=Ta+273.15; %% Soil Temperature [K]
Da = (2.11*1e-5)*(((Ts_k)/273.15)^1.94)*(Pre*100/101325);%% [m^2/s] vapor molecular diffusivity
u = Ws; %% [m/s] wind speed  ---
k= 0.4; %% Von Karman Constant
d=disp_h; %% Zero plane displacement [m]
z=zatm; %% Measurement Height [m]
%%% Hypothesis Logaritmic distribution of wind speed
us =  k*u/log((z-d)/zom); %%% Friction Velocity  [m/s]
%%%%%%
%Sllit = 2 ; %%% [m2 Litter / kg DM]
%%%%% 
if BLit<=0
    r_litter=0;
    alp_litter=1;
else
    OPT=3;
    switch OPT
        case 1
            %%%% [Sakaguchi and Zeng 2009]
            Llitter= Sllit*BLit; % [Leaf Area Index of Litter]
            r_litter = (1-exp(-Llitter))/(0.004*us); %% [s/m]
            %%%%%%%%%%%%
        case 2
            %%% Sauer et al 1995
            r_litter = 1/(0.0035+0.0011*us); %% [s/m]
        case 3
            %%% Litter Biomass   [Puthena Cordery 1996] [Sato et al 2004]
            %Olit = 0.064*(U.^(-0.51)-1).^(-0.42); %% [-] Relative humidity Litter  Bristow et al 1986
            Thick = 1.6*BLit; %%% [cm] with Blit in [kg/m2]
            MaxStoCap = 0.8*BLit; %%% coeff. [0.1-1];  [mm] with BLit [kg/m2]
            MinStoCap = 0.1*BLit; %%% coeff. [0.05-0.35];  [mm] with BLit [kg/m2]
            minOlit = 0.4*MinStoCap./MaxStoCap;
            Llitter = Sllit*BLit; %%% [mm]  fraction [m2 Litter / m2 ground ]
            Olit= 0.4*(In_Litter/MaxStoCap);  Olit(Olit<minOlit)=minOlit;
            %%%% [Park et al 1998]
            %F0=0.1377*(Thick/100)^(1.662*Olit^0.5841);
            F0=0.2731*(Thick/100)^(1.609*Olit^0.3952);
            r_litter = F0/Da; %% [s/m] 
            r_litter = Llitter*r_litter; % [s/m] 
    end
    alp_litter=1;
end
return