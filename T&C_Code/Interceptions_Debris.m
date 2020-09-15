%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute Interceptions       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[In_deb,SE_deb]=Interceptions_Debris(dt,Csno,Cdeb,...
    In_debtm1,In_max_deb,...
    Pr_liq,WR_SP,EIn_deb)
%%%INPUT
%dt time step [s]
dth = dt/3600; %% [h] 
%In_debtm1,  [mm] 
%In_max_deb  [mm] 
%%%%
Pr_liq= Pr_liq*dth; %% [mm]  
%EIn_deb  % Evaporation from Interception in rocks [mm/h]
%%% OUTPUTS 
%SE_deb  Storage excess  rock  [mm] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
%%%%%% ROCK INTERCEPTION 
In_deb = In_debtm1 + (1-Csno)*(Cdeb)*Pr_liq - EIn_deb*dth + ...
    (WR_SP)*Cdeb  ; % %% %% First updated Interception  [mm] 
SE_deb = (In_deb -In_max_deb)*(In_deb>In_max_deb); 
In_deb=  In_deb -(In_deb -In_max_deb)*(In_deb>In_max_deb);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return 