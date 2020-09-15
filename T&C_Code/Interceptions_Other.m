%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute Interceptions       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[In_urb,In_rock,SE_rock,SE_urb]=Interceptions_Other(dt,...
    Csno,Crock,Curb,Cice,...
    In_urbtm1,In_rocktm1,In_max_urb,In_max_rock,...
    Pr_liq,WR_SP,WR_IP,q_runon,EIn_urb,EIn_rock)
%%%INPUTS
%dt time step [s]
dth = dt/3600; %% [h] 
%In_urbtm1,
%In_rocktm1,
%In_max_H,
%In_max_L
%In_max_urb
%In_max_rock
%%%%
Pr_liq= Pr_liq*dth; %% [mm]  
%EIn_H=  Evaporation from Interception First Layer Canopy   %% [mm/h] 
%EIn_L= Evaporation from Interception Second Layer Canopy %% [mm/h]
%EIn_urb = % Evaporation from Interception in urban Landscape [mm/h]
%EIn_rock  % Evaporation from Interception in rocks [mm/h]
%%% OUTPUTS 
%SE_rock  Storage excess  rock 
%SE_urb Storage excess urban 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%% ROCK INTERCEPTION 
In_rock = In_rocktm1 + (1-Csno)*(1-Cice)*(Crock)*Pr_liq - EIn_rock*dth + ...
    (WR_SP*(1-Cice)+WR_IP)*Crock + q_runon*dth*Crock ; %In_rock(In_rock<0)=0; %% %% First updated Interception  [mm] 
SE_rock = (In_rock -In_max_rock)*(In_rock>In_max_rock); 
In_rock=  In_rock -(In_rock -In_max_rock)*(In_rock>In_max_rock);
%%%%% URBAN LANDSCAPE INTERCEPTION
In_urb = In_urbtm1 + (1-Csno)*(1-Cice)*(Curb)*Pr_liq - EIn_urb*dth + ...
    (WR_SP*(1-Cice)+WR_IP)*Curb + q_runon*dth*Curb;  %In_urb(In_urb<0)=0; %% %% First updated Interception  [mm]
SE_urb = (In_urb -In_max_urb)*(In_urb>In_max_urb);
In_urb=  In_urb -(In_urb -In_max_urb)*(In_urb>In_max_urb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return 


