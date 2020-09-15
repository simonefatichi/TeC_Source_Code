%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction BIOGEOCHEMISTRY_DYNAMIC   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[BfixN]= Biogeo_Bio_fixation(AAET,LAI,B,REXMY,Tdp)
%%% AAET Annual Average ET [mm/yr]
OPT_BFIX = 3; 
switch OPT_BFIX
    case 1
        BfixN=0;
    case 2
        %%%%%%%%%%  Biological Nitrogen Fixation
        %%%% Thomas et al 2013
        if LAI < 1
            Nfix_sym = max(0,(0.0018*AAET - 0.0289)); %% [gN/m2 yr]
        else
            Nfix_sym = 0;
        end
        Nfix_nonsym = 0.0006*AAET+0.0117; %% [gN/m2 yr]
        BfixN = Nfix_sym + Nfix_nonsym;   % [gN/m2 yr]
        %%% Dickinson et al 2002
        %fT = exp(0.08*(Ts-25));
        %Ssb = 0.001296; %% % [gN/m^2 d]  microbiological fixation
        %BfixN = Ssb*Ccrown*fT*exp(-0.5*B(24));
        %%%% Zahele and Friend 2010
        if B(31)+B(32) > 2
            BfixN = 0;
        end
    case 3
        if REXMY(3) > 0 
        %%% Biological N-fixation
        a= -3.62; b=0.27; c=25.14;
        s=-30;
        Costfix = s*(exp(a + b*Tdp*(1 - 0.5*Tdp/c))-2); %% [gC/gN]
        BfixN = sum(REXMY)/Costfix; %% [gN/m2 d]
        else
            BfixN = 0 ; 
        end 
end
%%% BfixN is into ammonium NH4+ form (Fisher et al 2010 GB Cycl)
return