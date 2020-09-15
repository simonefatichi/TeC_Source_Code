%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Canopy_Resistence&NPP      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[rs_sun,rs_shd,Ci_sun,Ci_shd,An,Rdark,Lpho,SIF,DCi]=Canopy_Resistence_An_Evolution(PAR_sun,PAR_shd,LAI,...
    Kopt,Knit,Fsun,Fshd,Citm1_sun,Citm1_shd,...
    Ca,ra,rb,Ts,Ta,Pre,Ds,...
    Psi_L,Psi_sto_50,Psi_sto_99,...
    CT,Vmax,DS,Ha,FI,Oa,Do,a1,go,e_rel,e_relN,gmes,rjv,mSl,Sl,Opt_CR)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT
% Knit Plant type Nitrogen extinction coefficient
% LAI Leaf Area Index
%%% Ci [umolCO2/mol] Leaf Interior  CO2 concentration
%%% IPAR [W/m^2] Photosyntecially active Radiation intercepted
%%% Ca [ppm]-[umolCO2/mol] Atmospheric CO2 concentration
%%% ra =[s/m]
%%% Ts = surface temperature [°C]
%%% Pre = Atmospheric Pressure [mbar]
%%% Ds = Vapor Pressure Deficit [Pa]
%%% O = [] Water Content
%%% Owp - Wilting point Water content []
%%% Oss - Stomatal closure begin Water content []
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% TP  --> 'CT' == 3  'CT' ==  4  %% Photosyntesis Typology for Plants
%%% Vmax Maximum Rubisco Capacity [umolCO2/ s m^2 ]
%%% FI Intrinsec quantum Efficiency [umolCO2/umolPhotons]
%%% Oa Intercellular Partial Pressure Oxygen [umolO2/mol] -
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Leuning (1995) Model of Photosyntesis Parameters
%%% Do % [Pa] empirical Parameter 1
%%% a1 % []  empirical Parameter 2
%%% go  [mol / s m^2] minimum Stomatal Conductance
%%% e_rel Age relative efficiency of photosyhtesis
%%%% OUTPUT
%%% CiF [umolCO2/mol] Leaf Interior  CO2 concentration
%%% rc  [s/m]  Stomatal resistence
%%% An Net Assimiltation Rate  % [umolCO2/ s m^2 ]
%%% Rdark %  %% [umolCO2/ s m^2 ] %% Surface Leaf Concentration
Citm1_sun(Citm1_sun<200)=200;
Citm1_shd(Citm1_shd<200)=200; 
%%%%%%%%%%%%%%%%%%%%
ANSW_SCA=1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%% SCALING FROM LEAF TO CANOPY
Vmax=Vmax*e_rel*e_relN; %%% Relative efficiency for age
%IPAR = PAR_sun + PAR_shd ; %%%
switch ANSW_SCA
    case 1
        %%% To be recomputed for Vmax only for LAI and with Kopt to avoid issue with SAI LAIdead 
        FsunV = (1.0 - exp(-Kopt*(LAI)))/(Kopt*(LAI));
        FsunV(FsunV<0.01)=0; FsunV(FsunV>1)=1;
        FshdV = 1- FsunV;
        %%% Two big leaves with Kn
        Can_sun = (1-exp(-(Kopt+Knit)*LAI))/(Kopt + Knit);
        Can_shd = (1-exp(-Knit*LAI))/Knit -(1-exp(-(Kopt+Knit)*LAI))/(Kopt + Knit);
        %%% Two big leaves with Kn
        Vmax_sun= Vmax*Can_sun/(LAI*FsunV);
        Vmax_shd= Vmax*Can_shd/(LAI*FshdV);
        if FsunV == 0
            Vmax_sun = 0; 
        end 
    case 2
        if Fsun == 0 
            %%%% Scaling Sla and not Nitrogen  %%% Thornton and Zimmermann
            Sla_sun = NaN;
            Sla_shd = (Sl*LAI +0.5*mSl*LAI^2 )/(LAI);
        else
            FsunI = (1.0 - exp(-Kopt*(LAI)))./(Kopt*(LAI));
            FshdI = 1- FsunI;
            %%%%% Scaling Sla and not Nitrogen  %%% Thornton and Zimmermann
            Sla_sun = (mSl + Sl*Kopt - (mSl*(Kopt*LAI+1) + Sl*Kopt)*exp(-Kopt*LAI) )/(Kopt^2*LAI*FsunI);
            Sla_shd = (Sl*LAI +0.5*mSl*LAI^2 - Sla_sun*LAI*FsunI )/(LAI*FshdI);
        end
        Vmax_sun= Vmax*Sl./Sla_sun;
        Vmax_shd= Vmax*Sl./Sla_shd;
    case 3
        %%% Simple Big Leaf
        %Can_tot = Cans_sun + Can_shd;
        %Can_tot = (1-exp(-Knit*LAI))/Knit;
        %%%%% Simple Big Leaf
        % Vmax_C= Vmax*Can_tot;
end
%%%%%%%%%%%%%%%%%%
go_sun =go; %% minimum canopy conductance
rb_sun = rb; %% Canopy Boundary layer resistance
go_shd =go; %% minimum canopy conductance
rb_shd = rb ; %% Canopy Boundary layer resistance
%%%
gmes_sun =gmes; %%
gmes_shd =gmes; %%
%%%
PAR_sun = PAR_sun/(LAI*Fsun);
PAR_shd = PAR_shd/(LAI*Fshd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SUNLIT FRACTION
if Fsun > 0
    %Opt_CR = optimset('TolX',3);
    %[Ci_sun]=fzero(@CO2_Concentration,Citm1_sun,Opt_CR,PAR_sun,Ca,ra,rb_sun,Ts,Ta,Pre,Ds,...
    %    O,Owp,Oss,CT,Vmax_sun,Tup,Tlow,DS,Ha,FI,Oa,Do,a1,go_sun);
    %%%
    [Ci_sun]=fzero(@CO2_Concentration,Citm1_sun,Opt_CR,PAR_sun,Ca,ra,rb_sun,Ts,Pre,Ds,...
        Psi_L,Psi_sto_50,Psi_sto_99,CT,Vmax_sun,DS,Ha,FI,Oa,Do,a1,go_sun,gmes_sun,rjv);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[CiF_sun,An_sun,rc_sun,Rdark_sun]= PHOTOSYNTESIS(Ci_sun,PAR_sun,Ca,ra,rb_sun,Ts,Ta,Pre,Ds,...
    %    O,Owp,Oss,...
    %    CT,Vmax_sun,Tup,Tlow,DS,Ha,FI,Oa,Do,a1,go_sun);
    %%%%
    [CiF_sun,An_sun,rc_sun,Rdark_sun,SIF_sun]= photosynthesis_biochemical(Ci_sun,PAR_sun,Ca,ra,rb_sun,Ts,Pre,Ds,...
        Psi_L,Psi_sto_50,Psi_sto_99,...
        CT,Vmax_sun,DS,Ha,FI,Oa,Do,a1,go_sun,gmes_sun,rjv);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    Ci_sun= 0; CiF_sun =0; An_sun=0; Rdark_sun=0; rc_sun =Inf; SIF_sun = 0; 
end
%%%% SHADOWED FRACTION
if Fshd > 0
    %Opt_CR = optimset('TolX',3);
    %[Ci_shd]=fzero(@CO2_Concentration,Citm1_shd,Opt_CR,PAR_shd,Ca,ra,rb_shd,Ts,Ta,Pre,Ds,...
    %    O,Owp,Oss,CT,Vmax_shd,Tup,Tlow,DS,Ha,FI,Oa,Do,a1,go_shd);
    %%%%
    [Ci_shd]=fzero(@CO2_Concentration,Citm1_shd,Opt_CR,PAR_shd,Ca,ra,rb_shd,Ts,Pre,Ds,...
        Psi_L,Psi_sto_50,Psi_sto_99,CT,Vmax_shd,DS,Ha,FI,Oa,Do,a1,go_shd,gmes_shd,rjv);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %[CiF_shd,An_shd,rc_shd,Rdark_shd]= PHOTOSYNTESIS(Ci_shd,PAR_shd,Ca,ra,rb_shd,Ts,Ta,Pre,Ds,...
    %    O,Owp,Oss,...
    %    CT,Vmax_shd,Tup,Tlow,DS,Ha,FI,Oa,Do,a1,go_shd);
    %%%
    [CiF_shd,An_shd,rc_shd,Rdark_shd,SIF_shd]= photosynthesis_biochemical(Ci_shd,PAR_shd,Ca,ra,rb_shd,Ts,Pre,Ds,...
        Psi_L,Psi_sto_50,Psi_sto_99,...
        CT,Vmax_shd,DS,Ha,FI,Oa,Do,a1,go_shd,gmes_shd,rjv);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    Ci_shd=0; CiF_shd =0; An_shd=0; Rdark_shd=0; rc_shd =Inf; SIF_shd = 0; 
end
%%%%
DCi_sun = Ci_sun -CiF_sun;
DCi_shd = Ci_shd -CiF_shd;
DCi = (DCi_sun + DCi_shd)/2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
An = An_sun*(LAI*Fsun) + An_shd*(LAI*Fshd);
Rdark = Rdark_sun*(LAI*Fsun) + Rdark_shd*(LAI*Fshd);
SIF = SIF_sun*(LAI*Fsun) + SIF_shd*(LAI*Fshd);
%%%;
rs_sun = rc_sun;%% [s/m] stomatal resistence
rs_shd = rc_shd; %%% [s/m] stomatal resistence
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%lanp = 0.506 ; %% [J / umol CO2 ]  %% Drewry et al 2012  Nikolov et al 1995
lanp = 0.469 ; %% [J / umol CO2 ]  %% Blanken et al 1997
Lpho = (An+Rdark)*lanp; %% [W/m2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return