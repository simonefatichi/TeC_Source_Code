%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction ROOT EXUDATION        %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[RexmyT]= Root_Exudation(NPP,Broot,Bam,Bem,Creserve,rMc,rNc,rPc,rKc,ParEx,NupIm,NavlI,EM,Tdp)
%%%%%OUTPUT
%% RexmyT [gC/m2 day] Root exudation and export to mycorrhiza
%%%% INPUT
%%% NPP [gc/m2 day]
%%% Broot [gC/m2]
%%% rMc [-] 0.65-1.85 Nutrient limitation status
dexmy = ParEx.dexmy; %% [-] exudation/export fraction
OPT_EXU = ParEx.OPT_EXU;
bfix = ParEx.bfix;
dtd=1;
Bam(Bam<1e-5)=0; Bem(Bem<1e-5)=0; 
%%% EM [0/1] AM or EM mycorrhiza
%%% Tdp Soil temperature
%%%%-> All computations are for /m2 PFT not ground
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kc1 = ParEx.kc1; kc2 = ParEx.kc2; kc3 = ParEx.kc3;
kn1 = ParEx.kn1; kn2 = ParEx.kn2; kn3 = ParEx.kn3;
kp1 = ParEx.kp1; kp2 = ParEx.kp2; kp3 = ParEx.kp3;
kk1 = ParEx.kk1; kk2 = ParEx.kk2; kk3 = ParEx.kk3;
if ( OPT_EXU == 1 ) ||  sum(NupIm)==0
    %%% Root Exudation and transfer to Mychorriza
    if NPP>0
        Cexu = (-5.71*rMc + 6.71); Cexu(Cexu<1)=1;
        Rexmy =Cexu*dexmy*NPP; %%% [gC / m^2 d]
    else
        Rexmy=0;
    end
    Rexmy = max(0,min(0.01*Creserve/dtd,Rexmy));
    %%%% Fraction of  direct exudation; (1-ex)= mychorrizal export
    fract_ex = 0.2; %%% from Farrar et al 2008
    RexmyT(1) = Rexmy*fract_ex;
    RexmyT(2) = Rexmy*(1-fract_ex);
    RexmyT(3) = 0;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %CostactR = 7.485./(rMc-0.60) + 10./Broot; %% [gC/gN]
    %CostactEM = (0.6237./(rMc - 0.60 ) + 300./Broot); %% [gC/gN]
    %CostactAM = (1.382./(rMc - 0.60 ) + 50./Broot); %% [gC/gN]
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if OPT_EXU == 2
        if NPP>0
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Biological N-fixation
            a= -3.62; b=0.27; c=25.14;
            s=-30;
            Costfix = s*(exp(a + b*Tdp*(1 - 0.5*Tdp/c))-2); %% [gC/gN]
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %Nsoil = NavlI(1)*365.25; %% [gN/m2 yr] N-available last year
            %Psoil = NavlI(2)*365.25; %% [gP/m2 yr] P-available last year
            %Ksoil = NavlI(3)*365.25; %% [gK/m2 yr] K-available last year
            Nsoil = min(NavlI(1)*365.25,NupIm(1)*365.25); %% [gN/m2 yr] N-available last year
            Psoil = min(NavlI(2)*365.25,NupIm(2)*365.25); %% [gP/m2 yr] P-available last year
            Ksoil = min(NavlI(3)*365.25,NupIm(3)*365.25); %% [gK/m2 yr] K-available last year
            %Nsoil = NupIm(1)*365.25; %% [gN/m2 yr] N-uptake last year
            %Psoil = NupIm(2)*365.25; %% [gP/m2 yr] P-uptake  last year
            %Ksoil = NupIm(3)*365.25; %% [gK/m2 yr] K-uptake  last year
            %%%%
            CostactR_n =  kn1/Nsoil + kc1./Broot; %% [gC/gN]
            CostactEM_n = kn2/Nsoil + kc2./Broot; %% [gC/gN]
            CostactAM_n = kn3/Nsoil + kc3./Broot; %% [gC/gN]
            CostactR_p =  kp1/Psoil + kc1./Broot; %% [gC/gP]
            CostactEM_p = kp2/Psoil + kc2./Broot; %% [gC/gP]
            CostactAM_p = kp3/Psoil + kc3./Broot; %% [gC/gP]
            CostactR_k =  kk1/Ksoil + kc1./Broot; %% [gC/gK]
            CostactEM_k = kk2/Ksoil + kc2./Broot; %% [gC/gK]
            CostactAM_k = kk3/Ksoil + kc3./Broot; %% [gC/gK]
            %%%%
            if bfix == 1
                CostacqN = (1/Costfix  + 1/CostactR_n + EM/CostactEM_n + (1-EM)/CostactAM_n)^-1; %% [gC/gN]
            else
                CostacqN = (1/CostactR_n + EM/CostactEM_n + (1-EM)/CostactAM_n)^-1; %% [gC/gN]
            end
            CostacqP = (1/CostactR_p + EM/CostactEM_p + (1-EM)/CostactAM_p)^-1; %% [gC/gP]
            CostacqK = (1/CostactR_k + EM/CostactEM_k + (1-EM)/CostactAM_k)^-1; %% [gC/gK]
            %%%%
            [Rexmy,p] = max([CostacqN*NupIm(1),CostacqP*NupIm(2),CostacqK*NupIm(3)]); %%% [gC/ m2 day]
            Rexmy = max(0,min(0.01*Creserve/dtd,Rexmy));
            switch p
                case 1  %% N
                    if bfix == 1
                        Nui_bnf = Rexmy./Costfix ;
                    else
                        Nui_bnf = 0;
                    end
                    Nui_ex = Rexmy./CostactR_n ;
                    Nui_em = EM.*Rexmy./CostactEM_n ;
                    Nui_am = (1-EM).*Rexmy./CostactAM_n ;
                case 2 %% P
                    Nui_bnf = 0;
                    Nui_ex = Rexmy./CostactR_p ;
                    Nui_em = EM.*Rexmy./CostactEM_p ;
                    Nui_am = (1-EM).*Rexmy./CostactAM_p ;
                case 3 %% K
                    Nui_bnf = 0;
                    Nui_ex = Rexmy./CostactR_k ;
                    Nui_em = EM.*Rexmy./CostactEM_k ;
                    Nui_am = (1-EM).*Rexmy./CostactAM_k ;
            end
            Nui_tot = Nui_bnf + Nui_ex +  Nui_em +  Nui_am;
            RexmyT(1) = Rexmy*Nui_ex/Nui_tot ;   %% [gC / m^2 d]
            RexmyT(2) = Rexmy*(Nui_em+Nui_am)/Nui_tot; %% [gC / m^2 d]
            RexmyT(3) = Rexmy*Nui_bnf/Nui_tot; %% [gC / m^2 d]
        else
            RexmyT(1) = 0;
            RexmyT(2) = 0;
            RexmyT(3) = 0;
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if OPT_EXU == 3
        if NPP>0
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Biological N-fixation
            a= -3.62; b=0.27; c=25.14;
            s=-30;
            Costfix = s*(exp(a + b*Tdp*(1 - 0.5*Tdp/c))-2); %% [gC/gN]
            %%%%%%%%%%%%%%%%%%%%%%%%%%
            Nsoil = min(NavlI(1),NupIm(1)); %% [gN/m2] N-available last year
            Psoil = min(NavlI(2),NupIm(2)); %% [gP/m2] P-available last year
            Ksoil = min(NavlI(3),NupIm(3)); %% [gK/m2] K-available last year
            %%%%
            CostactR_n =  kn1/(Nsoil*Broot^1.5); %% [gC/gN]
            CostactEM_n = kn2/(Nsoil*Bem^0.5); %% [gC/gN]
            CostactAM_n = kn3/(Nsoil*Bam^0.5); %% [gC/gN]
            CostactR_p =  kp1/(Psoil*Broot^1.5); %% [gC/gP]
            CostactEM_p = kp2/(Psoil*Bem^0.5); %% [gC/gP]
            CostactAM_p = kp3/(Psoil*Bam^0.5); %% [gC/gP]
            CostactR_k =  kk1/(Ksoil*Broot^1.5); %% [gC/gK]
            CostactEM_k = kk2/(Ksoil*Bem^0.5); %% [gC/gK]
            CostactAM_k = kk3/(Ksoil*Bam^0.5); %% [gC/gK]
            %%%%
            if bfix == 1
                CostacqN = (1/Costfix  + 1/CostactR_n + EM/CostactEM_n + (1-EM)/CostactAM_n)^-1; %% [gC/gN]
            else
                CostacqN = (1/CostactR_n + EM/CostactEM_n + (1-EM)/CostactAM_n)^-1; %% [gC/gN]
            end
            CostacqP = (1/CostactR_p + EM/CostactEM_p + (1-EM)/CostactAM_p)^-1; %% [gC/gP]
            CostacqK = (1/CostactR_k + EM/CostactEM_k + (1-EM)/CostactAM_k)^-1; %% [gC/gK]
            %%%%
            [Rexmy,p] = max([CostacqN*NupIm(1),CostacqP*NupIm(2),CostacqK*NupIm(3)]); %%% [gC/ m2 day]
            Rexmy = max(0,min(0.01*Creserve/dtd,Rexmy));
            switch p
                case 1  %% N
                    if bfix == 1
                        Nui_bnf = Rexmy./Costfix ;
                    else
                        Nui_bnf = 0;
                    end
                    Nui_ex = Rexmy./CostactR_n ;
                    Nui_em = EM.*Rexmy./CostactEM_n ;
                    Nui_am = (1-EM).*Rexmy./CostactAM_n ;
                case 2 %% P
                    Nui_bnf = 0;
                    Nui_ex = Rexmy./CostactR_p ;
                    Nui_em = EM.*Rexmy./CostactEM_p ;
                    Nui_am = (1-EM).*Rexmy./CostactAM_p ;
                case 3 %% K
                    Nui_bnf = 0;
                    Nui_ex = Rexmy./CostactR_k ;
                    Nui_em = EM.*Rexmy./CostactEM_k ;
                    Nui_am = (1-EM).*Rexmy./CostactAM_k ;
            end
            Nui_tot = Nui_bnf + Nui_ex +  Nui_em +  Nui_am;
            RexmyT(1) = Rexmy*Nui_ex/Nui_tot ;   %% [gC / m^2 d]
            RexmyT(2) = Rexmy*(Nui_em+Nui_am)/Nui_tot; %% [gC / m^2 d]
            RexmyT(3) = Rexmy*Nui_bnf/Nui_tot; %% [gC / m^2 d]
        else
            RexmyT(1) = 0;
            RexmyT(2) = 0;
            RexmyT(3) = 0;
        end
    end
    if Broot == 0
        RexmyT(1:3) = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end