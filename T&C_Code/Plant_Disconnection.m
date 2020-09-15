%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction - Plant Disconnection       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%REFERENCES
function[PLD]=Plant_Disconnection(Psi_x,Psi_l,Axyl,PsiL50,PsiL00,PsiX50)
%%%%%%%%%%%%%%%%%
%%% INPUTS
%%% OUTPUTS
%%%%%%% Plant Dead ---
PLD=0;
if Axyl > 0
    if Psi_x < PsiX50
        PLD = 1;
    end
else
    PLCs = 0.02;  %%[-]  specific loss of conductivity
    p= log((1 -PLCs)/PLCs)/(PsiL00 - PsiL50);%% [1/MPa]
    q=-p*PsiL50; %%[-]
    PLC_crit= 0.90; %% [fraction] 
    Psi_l_crit = 1/p*(log(1./PLC_crit -1)-q); %% [MPa]
    if Psi_l < Psi_l_crit
        PLD=1;
    end
end
return 