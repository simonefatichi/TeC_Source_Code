%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Soil Resistence            %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%REFERENCES %%    Kondo et al 1990 %% Ivanov et al 2008  Cox et al 1999 %%
%%Laio et al., 2001; Oleson et al., 2008; Lawrence et al., 2010
%%% Haghighi et al 2013  %% Shahraeeni et al 2012
function[r_soil,b_soil,alp_soil]=Soil_Resistence(Ts,Pre,Ws,ea,q_runon,O,Ks,Osat,Ohy,L,Pe,O33,alpVG,nVG,lVG,Ks_mac,Omac,alpVGM,nVGM,lVGM,Phy,s_SVG,bVG,SPAR)
%%%INPUTS
%Ts = %% %% [°C] Soil temperature
%O = [] %% water content
%Osat [] Saturation Moisture
%%% Ohy [] Hygroscopic Moisture Evaporation cessation
Ts_k=Ts+273.15; %% Soil Temperature [K]
%%% OUTPUTS
%r_soil=  [s/m] soil resistence
%b_soil [0-1] beta factor
%alp_soil = [0-1] relative humidity
%%% PARAMETERS
row=1000; %%% [kg/m^3] water density
g= 9.81; %% [m/s^2] gravity acceleration
%Rd =287.05; %% [J/kg K] dry air gas constant
Rd = 461.5; %%  [J/kg K] water vapor gas constant
Da = (2.11*1e-5)*(((Ts_k)/273.15)^1.94)*(Pre*100/101325);%% [m^2/s] vapor molecular diffusivity
esat=611*exp(17.27*Ts./(237.3+Ts)); %% [Pa]
%%%%%%%%%
[Ko,Po]=Conductivity_Suction(SPAR,Ks,Osat,Ohy,L,Pe,O33,alpVG,nVG,lVG,Ks_mac,Omac,alpVGM,nVGM,lVGM,Phy,s_SVG,bVG,O) ;
Po(Po<0)=0; %% [mm]
alp_soil = exp(-Po*g/(1000*Rd*Ts_k)); %%% [-]
%%%%%%%%%%%%%%%%% From Haghighi et al 2013
Psz= (11.12*nVG^3.286)*1e-6; %% [m] 40-200 um  Size of the pores --  Particle Size/3
dm=2.26*1e-3/sqrt(Ws); %%% Shahraeeni et al 2012 [m] Boundary Layer Thickness
%%%%%%%%%%%%%%
%%%% From Lehmann et al 2018 
mVG= 1-1./nVG;
hc = -1/alpVG*((nVG-1)/nVG)^((1-2*nVG)/nVG); %%% [mm] 
%%%%%%%
ANSW=4;
switch ANSW
    case 1
        %%%%%%%%%%%%%%%%%%%%% Oleson 2008
        r_soil= exp(8.206 - 4.255*(O-Ohy)/(Osat-Ohy));
        %%%%%%%%% Constrain Resistence
        if (O <= Ohy)
            r_soil = Inf;
        end
        b_soil=1;
    case 2
        %%%%%%%%%%%%%%%%%%%%%%%
        %%Lee and Pielke [1992]
        b_soil = 0.25*(1-cos(O/O33*pi)).^2;
        if O >= O33
            b_soil=1;%[-]
        end
        if O<= Ohy
            b_soil=0;
        end
        r_soil = 0;
        %%%%%%%%%%%%%%%%%%%%%%
    case 3
        %%% Cathy Model Camporese et al [2010] Paniconi personal
        %%% communication
        r_soil=0;
        alpha=1;
        E = ro*(alpha*qsat-qa)/(ra+r_soil);  %%% [kg/m^2.s]
        Khalf= 0.5*(Kmin + Ko);
        Mev = -Khalf*(1- (Po-Pmin)/Dz);
        Ev = max(0,min(E,Mev)); %% [mm/h]
        r_soil = -ra + ro*(1*qsat-qa)/(Ev*row/1000/3600); %% [s/m]
    case 4
        %%%%%%%%%%%%%%%%%%%% Haghighi et al 2013
        gammap = (alp_soil*esat-ea)/(row*Rd*Ts_k); %% [-]
        if gammap < 0
            r_soil=0;
        else
            rsv = gammap/(4*Ko/(1000*3600)); %%% Internal soil viscous resistance [s/m]
            f_O= (2/pi)*(sqrt(pi/(4*O))-1)/sqrt(4*O); %% [-]
            rvbl = (dm +Psz*f_O)/Da ; %%% [s/m] viscous boundary layer resistance
            r_soil = rvbl + rsv;
        end
        %%%%%%%%% Constrain Resistence
        if (O <= Ohy)
            r_soil = Inf;
        end
        b_soil=1;
    case 5  %% %Only for Van-Genucthen option 
        if SPAR ~= 1 
            disp('ERROR: OPTION FOR SOIL RESISTANCE REQUIRES van-Genuchten soil hydraulic parameters')
        end 
        %%%%%%%%%%%%%%%%%%%% Lehmann et al 2018 ;; Haghighi et al 2013
        gammap = (alp_soil*esat-ea)/(row*Rd*Ts_k); %% [-]
        E0= Da/dm*gammap*1000/3600; %% [mm/h]
        Se = 1/((1+abs(alpVG*hc)^nVG)^mVG);
        Se(Se>1)=1;
        Khc= Ks*((Se)^(0.5))*(1-(1-(Se)^(1/mVG))^mVG)^2; %%% [mm/h]
        dhdz = 1 + E0/(4*Khc); %%%[-]
        %%%%%
        if gammap < 0
            r_soil=0;
        else
            rsv = gammap/(4*Ko*dhdz/(1000*3600)); %%% Internal soil viscous resistance [s/m]
            f_O= (2/pi)*(sqrt(pi/(4*O))-1)/sqrt(4*O); %% [-]
            rvbl = (dm +Psz*f_O)/Da ; %%% [s/m] viscous boundary layer resistance
            r_soil = rvbl + rsv;
        end
        %%%%%%%%% Constrain Resistence
        if (O <= Ohy)
            r_soil = Inf;
        end
        b_soil=1;
end
if q_runon > 0
    r_soil = 0;
    alp_soil = 1;
    b_soil=1;
end
return