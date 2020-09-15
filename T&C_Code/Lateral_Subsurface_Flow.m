%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Lateral Subsurface Flow     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Qi_out]=Lateral_Subsurface_Flow(Koh,dz,Slo_pot,aT,cosalp,sinalp,SN,I1)
%%REFERENCES %%  
%%%INPUTS
%LT = 1; %%[mm]  cell length 
%A = 1; %%[mm^2] Area cell 
%aT = A/LT; 
%%% OUTPUTS
%%%Qi_out Lateral subsurface flow [mm/h] 
%%%%%%%%%%%%%%%
To = Koh.*dz;  %%% Trasmissivity  unsaturated/saturated [mm^2/h]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%
%Qi_out = To.*Slo_pot/aT; %%%% [mm/h] 
%Qi_out = (To/aT).*(cosalp.*Slo_pot+sinalp); %%%% [mm/h] 
if SN==1
    %Qi_out = I1.*(To/aT).*cosalp.*(Slo_pot); %%%% [mm/h]
    Qi_out = I1.*(To/aT).*(sinalp); %%%% [mm/h]
else
    %Qi_out = (To/aT).*cosalp.*(Slo_pot); %%%% [mm/h]
    Qi_out = (To/aT).*(sinalp); %%%% [mm/h]
end
%%%%%%%%%%%%%%%%%%%%%%%%%
return 