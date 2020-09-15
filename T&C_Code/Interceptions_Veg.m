%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute Interceptions       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[In_H,In_L,In_Litter,Dr_H,Dr_L,WIS]=Interceptions_Veg(dt,...
    Ccrown,Cfol_H,Cfol_L,CLitter,Cbare,Csno,Cice,Crock,Curb,Cwat,dw_SNO,In_Ltm1,In_Htm1,In_Littertm1,...
    In_max_H,In_max_L,BLit,...
    Pr_liq,EIn_H,EIn_L,Elitter,q_runon,WR_SP,WR_IP,gc,Kc)
%%%INPUTS
%dt time step [s]
dth = dt/3600; %% [h] 
%Cveg_H 
%Cveg_L
%Csno,
%dw_SNO, 
%In_Ltm1 
%In_Htm1 ,
%In_max_H,
%In_max_L
%BLit [kg DM / m2 PFT]
%%%%
q_runon = q_runon*dth; %% [mm] 
Pr_liq= Pr_liq*dth; %% [mm]  
%EIn_H=  Evaporation from Interception First Layer Canopy   %% [mm/h] 
%EIn_L= Evaporation from Interception Second Layer Canopy %% [mm/h]
%%% OUTPUTS 
%In_H,
%In_L,
%Dr_H
%Dr_L
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% CANOPY RAIN INTERCEPTION
%%% %%REFERENCES [Deardorff (1978)---- Rutter et al.l, 1975 -- Mahfouf and Jacquemin 1989
%gc=3.7; %%% [1/mm]
%Kc=0.001*60*dth; %%%% [mm] -- Mahfouf and Jacquemin 1989
%%%% First Layer 
Pr_H = Ccrown.*Pr_liq; %%% [mm]  
In_H = In_Htm1 + (1-dw_SNO)*Pr_H.*Cfol_H - EIn_H*dth;  %In_H(In_H<0)=0; %% %% First updated Interception  [mm] 
SE_H=  (In_H -In_max_H).*(In_H>In_max_H); %% Storage Excess [mm]
In_H=  In_H -SE_H;  %%[mm]
Dr_H = Kc*exp(gc*(In_H - In_max_H)).*(In_H>0);   %%% Drainage  first layer 
In_H = In_H - Dr_H;  
Dr_H = Dr_H + In_H.*(In_H<0); In_H(In_H<0)=0; %% %% First updated Interception  [mm] 
Dr_H = Dr_H + SE_H; %% Dripping and Saturation Excess  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%% Second Layer 
Pr_L = (1-Cfol_H).*Pr_H + Dr_H; %%% [mm] 
In_L = In_Ltm1 + (1-Csno)*(1-Cice)*Pr_L.*Cfol_L - EIn_L*dth;  %In_L(In_L<0)=0; %% %% First updated Interception  [mm] 
SE_L=  (In_L -In_max_L).*(In_L>In_max_L); %% Storage Excess [mm]
In_L=  In_L -SE_L;  %%[mm]
Dr_L = Kc*exp(gc*(In_L - In_max_L)).*(In_L>0);   %%% Drainage  first layer 
In_L = In_L - Dr_L;
Dr_L = Dr_L + In_L.*(In_L<0); In_L(In_L<0)=0; %% %% First updated Interception  [mm]
Dr_L = Dr_L + SE_L; %% Dripping and Saturation Excess
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% WIS = [mm] Water Incoming to First Soil Layer
WIS_veg = (1-Csno)*(1-Cice)*sum((Pr_liq.*Ccrown.*(1-Cfol_H) + Dr_H).*(1-Cfol_L) + Dr_L) + ...
    + sum(Ccrown)*(WR_SP*(1-Cice)+WR_IP)*(1-Crock-Curb-Cwat)+ ...
    sum(Dr_H+Dr_L)*max(Csno,Cice) + sum(Ccrown)*q_runon*(1-Crock-Curb-Cwat); %%% [mm];
WIS_bs =  Pr_liq*Cbare*(1-Csno)*(1-Cice) + (1-sum(Ccrown))*(WR_SP*(1-Cice)+WR_IP)*(1-Crock-Curb-Cwat)+ ...
     (1-sum(Ccrown))*q_runon*(1-Crock-Curb-Cwat); %%% [mm];
%%% TOTAL 
WIS = WIS_veg + WIS_bs;  %%% [mm];
%WIS= (1-Csno)*(1-Cice)*sum((Pr_liq.*Ccrown.*(1-Cfol_H) + Dr_H).*(1-Cfol_L) + Dr_L) + ...
%    Pr_liq*Cbare*(1-Csno)*(1-Cice) + (WR_SP*(1-Cice)+WR_IP)*(1-Crock-Curb-Cwat)+ ...
%    sum(Dr_H+Dr_L)*max(Csno,Cice) + q_runon*(1-Crock-Curb-Cwat); %%% [mm];
%%%%%%%%%
%%% Litter Layer
if sum(BLit)>0
    %%%  [Puthena Cordery 1996] [Sato et al 2004]
    MaxStoCap = 0.8*BLit; %%% coeff. [0.1-1];  [mm] with BLit [kg/m2]
    MinStoCap = 0.1*BLit; %%% coeff. [0.05-0.35];  [mm] with BLit [kg/m2]
    In_max_Litter = sum(Ccrown.*(MaxStoCap - MinStoCap)); %% [mm]
    CLitter = sum(Ccrown.*CLitter); %% [-] %%% CLitter fraction cover by litter
    In_Litter = In_Littertm1 + CLitter*WIS_veg - Elitter*dth ; % First updated Interception  [mm]
    SE_Litter = (In_Litter -In_max_Litter)*(In_Litter>In_max_Litter);
    In_Litter =  In_Litter - SE_Litter;
    %%%% 
    WIS = SE_Litter +(1-CLitter)*WIS_veg + WIS_bs; %%% [mm];
else
    %WIS=WIS;
    In_Litter = 0;
end       
%%%%%%%%%%%%%%%%%%%%%
return 


