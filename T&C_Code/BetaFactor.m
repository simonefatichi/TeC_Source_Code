%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Plant Water Status           %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% References 
function[Bfac,Bfac_all]= BetaFactor(Psi_l,PsiL00,PsiL50,PsiG50,PsiG99)
%%%% INPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PsiL50 = ;%%[MPa]  Water Potential at 50% loss conductivity
%PsiLs =  ; %% [MPa]  Water Potential at PLCs% loss conductivity
PLCs = 0.02;  %%[-]  specific loss of conductivity
%%% Leaf Water Status 
p= log((1 -PLCs)/PLCs)/(PsiL00 - PsiL50);%% [1/MPa]
q=-p*PsiL50; %%[-]
%%%%
PLC = 1./(1+exp(p.*Psi_l+q)); %% [fraction]
PLC(PLC>1)=1; PLC(PLC<0)=0;
Bfac=1-PLC.*(PLC>0.05); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Bfac = max(0,min(1,(Psi_l -PsiL99)/(PsiL50- PsiL99)));
%%%%
%%%% For Env. Control Growth and Allocation 
%%% Water Limitation Function
%PsiG50 =-0.45;%%[MPa]
%PsiG99 = -1.2; %% [MPa]
PLCs = 0.99;  %%[-]  specific loss of conductivity
p= log((1 -PLCs)/PLCs)/(PsiG99 - PsiG50);%% [1/MPa]
q=-p*PsiG50; %%[-]
PLC = 1./(1+exp(p.*Psi_l+q));  %% [fraction] Reduction Growth
PLC(PLC>1)=1; PLC(PLC<0)=0;
Bfac_all = 1-PLC; 
return
%%%%