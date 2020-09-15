%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  AVALANCHE_MODULE           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[SND,SWE,ros,S_exit]= AVALANCHES(DTM,cellsize,Area,Asur,Slo_top,SND,SWE,ros)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Reference: Bernhardt and Schulz (2010)
%%%%%% OUTPUT
%%% SND [m] n-coordinate 
%%% SWE [mm] v-coordinate 
%%% ros [kg/m^3]  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% INPUT
%%% DTM  [m a.s.l.]
%%% cellsize [m]
%%% Asur [-] % Effective Area / Projected Area
%%% Slo_top [-] fraction 
%%% SND [m] n-coordinate 
%%% SWE [mm] v-coordinate 
%%% ros [kg/m^3]  
%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameter
a=0.1012; %% [1/[°]]
C=99.05; %%[m]
ros_ava = 400; %% [kg/m^3] %% Density of avalanche deposition   Sovilla et al., [2006]
row=1000; %%[kg/m^3] 
%%% COMPUTE
S_exit=0;
SWM=SWE ; 
r=0; 
while sum(sum(SWM))>0.1 
    r=r+1; 
    %disp(sum(sum(SWM)));
    HSNO = DTM+SND.*cos(atan(Slo_top)); %%% Total Elevation Surface [m]
    [R_snow,SS] = dem_flow(HSNO,cellsize,cellsize);
    Slo_snow=Slope_Aspect_indexes(HSNO,cellsize,'mste'); %% [-]
    Slo_snow(SND==0)=0;
    %Slo_snow=(180/pi)*atan(SS); %% Slope in [°]
    Slo_snowG=(180/pi)*atan(Slo_snow); %% Slope in [°]
    Slo_snowG(Slo_snowG<0)=0;
    %%% Maximum Water Holding Capacity
    Snd_max=C*exp(-a*Slo_snowG); %%[m]
    Snd_max(Snd_max<0.05)=0.05; %%[m]
    %%%
    SNR=SND-Snd_max; SNR(SNR<0)=0; %%[m] Snow to Route
    fSNR = SNR./SND; %% [-] Fraction of snow to route
    %%%% SWE in the v-coordinate
    SWR=SWE.*fSNR; %% [mm] Snow water equivalent to route
    rhoR = max(ros,ros_ava); %%[kg/m^3] Density of moved material
    %%%% Correction on the maximum movement to not invert the slope 
    MDEP=(SWR./Asur).*(row./rhoR); %% Moved Depth %%[mm] 
    MDEP = min(MDEP,1000*cellsize*Slo_snow/2); %% Moved Depth corrected %%[mm] 
    MDEP(MDEP<0)=0;
    SWR = (MDEP.*Asur)./(row./rhoR); %% Corrected Snow water equivalent to route
    %%%%%%%
    T_flow_sno = flow_matrix(HSNO,R_snow,cellsize,cellsize); %% Flow Matrix %%%%
    [SWM]=Flow_Routing_Step2(HSNO,T_flow_sno,SWR); %%[mm] Moved Snow water equivalent
    [SWrosM]=Flow_Routing_Step2(HSNO,T_flow_sno,rhoR.*SWR); %%[kg/m^3 * mm] Moved Snow water equivalent *density 
    SWEn = SWM + (SWE - SWR) ; %% [mm] New SWE
    S_exit= S_exit + (nansum(nansum(SWE))- nansum(nansum(SWEn))); %% [mm] SWE exit domain 
    rosn = (SWrosM + (SWE-SWR).*ros)./SWEn; %% %%[kg/m^3] New Density
    %%% SND in the projected area n-coordinate
    SNDn= 0.001*(SWEn./Asur).*(row./rosn); %% Snowpack [m]
    SNDn(SNDn<0) = 0;  rosn(rosn < 0) = 0;
    SNDn(SWEn==0) = 0;  rosn(SWEn==0) = 0;
    %%%
    SND=SNDn;
    SWE=SWEn;
    ros=rosn;
end
S_exit = S_exit*(cellsize^2)/Area; %%[mm]
end 