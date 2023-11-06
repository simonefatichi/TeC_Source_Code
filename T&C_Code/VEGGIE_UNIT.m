%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  VEGGIE_UNIT                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[LAI,B,NPP,ANPP,Rg,RA,Rms,Rmr,Rmc,PHE_S,dflo,AgeL,e_rel,e_relN,LAIdead,NBLeaf,Sr,Slf,Sfr,Swm,Sll,...
    Rexmy,Rootl,AgeDL,Bfac_day,Bfac_week,NPPm,Tsmm,NupIm,PAR_Im,NBL_Im,RB,FNC,Nreserve,Preserve,Kreserve,rNc,rPc,rKc,ManI]= VEGGIE_UNIT(Btm1,PHE_Stm1,dflotm1,AgeLtm1,AgeDLtm1,...
    Ta,Tdp,PAR,Psi_x,Psi_l,An,Rdark,NPPtm1,jDay,Datam,NPPI,TdpI,Bfac_weekI,NupI,NavlI,PAR_I,NBL_I,NBLeaftm1,....
    L_day,Lmax_day,Veg_Param_Dyn,cc,...
    Nreservetm1,Preservetm1,Kreservetm1,Nuptake,Puptake,Kuptake,FNCtm1,Se_bio,Tdp_bio,...
    ParEx,EM,Bam,Bem,Mpar,TBio,OPT_EnvLimitGrowth,OPT_VCA,OPT_VD,OPT_SoilBiogeochemistry)

%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% INPUT dth  = 1 [h]
%%%% LAItm1 LAI step before
%%%% Btm1 [gC / m^2] Carbon Pool step before
%%% dflotm1 [day] from Leaf onset step before
%%% AgeLtm1 [day] Average Age of Leaf step before
%%% PHE_Stm1 [#] Phenology State step before
%%%%%% Ta [°C] Air Temperature
%%%%%% Ts [°C] Soil Temperature
%%%%%% OV [] Water Conten Root Layer
%%%%%% An  [umolCO2/ s m^2 ] Net Assimilation
%%%%%% Rdark  [umolCO2/ s m^2 ] Dark Respiration
%%%%% Psi_x,
%%%%%% Psi_l
%%%%%% PsiL50,
%%%%% PsiL00,
%%%%% PsiX50
%%%% Osat
%%% i index of computation SVAT
%%% j index of computation VEGETATION
%%% NPPpr Net Primary Productivity
% Sl =  [m^2 gC] specific leaf area of  biomass
%gR = % [0.22 - 0.28] growth respiration  [] -- [Rg/(GPP-Rm)]
%r=  %% [0.066 -0.011]respiration rate at 10° [gC/gN d ]
%Ns = %%[ 50 -330] Sapwood Carbon Nitrogen  [gC/gN] Sapwood
%Nr =  %% [30- 60] [gC/gN] Fine root  Carbon Nitrogen
%drn =  %% turnover rate root  [1/d]
%dsn =  % [ normal transfer rate sapwood [1/d]
%dc_C [1/ d °C] -- [Factor of increasing mortality Cold]
%Tcold =  %% [°C] [ Cold stress for Leaf Fall  ]
%dd_max = %%  [1/d]  0.0250 -- 0.005-0.025 death maximum for drought
%aSE = %% Plant Type -- 1 Seasonal Plant --  0 Evergreen  -- 2 Grassspecies-- 3 Crops
%Trr =  %% Translocation rate [gC/m^2 d]
%Tlo = %% Temperature for Leaf onset
%Tls = %% Temperature for Leaf Shed
%%dmg = %%% 60 Tree 30 Grasses Length Period of Max Growth [d]
%%% Bfac_lo = %% Phenology Leaf Onset Water Stress
%%% Bfac_ls = %% Phenology Leaf Shed Water Stress
%%% LAI_min =
%%% age_cr= %% [day] Critical Leaf Age
%%% LtR, = Leaf to root ratio [0.75 -1.5]
%%% eps_ac = allocation to reserve parameter [0-1]
%%%%L_day, %%% Day light lenght [h]
%%% LDay_cr %% Threshold for senescence day light [h]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%
Sl = Veg_Param_Dyn.Sl(cc);
mSl = Veg_Param_Dyn.mSl(cc);
Stoich = Veg_Param_Dyn.Stoich(cc) ;
r = Veg_Param_Dyn.r(cc);
gR = Veg_Param_Dyn.gR(cc);
LtR = Veg_Param_Dyn.LtR(cc);
eps_ac = Veg_Param_Dyn.eps_ac(cc);
aSE = Veg_Param_Dyn.aSE(cc);
Trr = Veg_Param_Dyn.Trr(cc);
dd_max = Veg_Param_Dyn.dd_max(cc);
dc_C = Veg_Param_Dyn.dc_C(cc);
Tcold = Veg_Param_Dyn.Tcold(cc);
drn = Veg_Param_Dyn.drn(cc);
dsn = Veg_Param_Dyn.dsn(cc);
age_cr = Veg_Param_Dyn.age_cr(cc); 
Bfac_lo = Veg_Param_Dyn.Bfac_lo(cc); 
Bfac_ls = Veg_Param_Dyn.Bfac_ls(cc); 
Tlo = Veg_Param_Dyn.Tlo(cc); 
Tls = Veg_Param_Dyn.Tls(cc); 
mjDay = Veg_Param_Dyn.mjDay(cc); 
LDay_min = Veg_Param_Dyn.LDay_min(cc);
dmg = Veg_Param_Dyn.dmg(cc);
Mf = Veg_Param_Dyn.Mf(cc);
Wm =Veg_Param_Dyn.Wm(cc);
LAI_min = Veg_Param_Dyn.LAI_min(cc);
LDay_cr = Veg_Param_Dyn.LDay_cr(cc);
PsiG50 = Veg_Param_Dyn.PsiG50(cc);
PsiG99 = Veg_Param_Dyn.PsiG99(cc);
gcoef = Veg_Param_Dyn.gcoef(cc);
fab = Veg_Param_Dyn.fab(cc);
fbe = Veg_Param_Dyn.fbe(cc);
Klf = Veg_Param_Dyn.Klf(cc);
ff_r = Veg_Param_Dyn.ff_r(cc);
PAR_th = Veg_Param_Dyn.PAR_th(cc);
PsiL50=Veg_Param_Dyn.PsiL50(cc);
PsiL00=Veg_Param_Dyn.PsiL00(cc);
soCrop=Veg_Param_Dyn.soCrop(cc); 
Sl_emecrop=Veg_Param_Dyn.Sl_emecrop(cc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
jDay_cut = Mpar.jDay_cut;
LAI_cut = Mpar.LAI_cut;
jDay_harv = Mpar.jDay_harv;
B_harv = Mpar.B_harv;
%%%%  Introducing Girdling effects 
if  sum(abs(datenum(Datam(1),Datam(2),Datam(3),Datam(4),0,0)-Mpar.Date_girdling)<=0.49)>=1
    GirdOpt = Mpar.fract_girdling;
else
    GirdOpt =0;
end
%%%%%%%%%%%%%%%%% OUTPUT
%%% dflo [day] from Leaf onset
%%% AgeL [day] Average Age of Leaf
%%% PHE_S [#] Phenology State
%%%%
Btm1=squeeze(Btm1);
%%%%
dtd = 1; %% day
%%%
n2=30 ; %% 30 days
n3=7 ; %% 7 days
n4=365; %% 365 days
n5 =45; %% 15 days
n6=10; %% 10 days
%%%%%%%%%
%%% Refering to previous day -- Ta,Tdp,Psi_x,Psi_l,An,Rdark NPPpr
%%%%% Average 1 day
Tam = mean(Ta); Tam(Tam<-46)=-46;
Tsm = mean(Tdp); Tsm(Tsm<-46)=-46;
%%%%
[Bfac,Bfac_alloc]= BetaFactor(Psi_l,PsiL00,PsiL50,PsiG50,PsiG99);
Bfac_day = mean(Bfac);
Bfac_alloc = mean(Bfac_alloc);
%%%%
An = mean(An);
Rdark=mean(Rdark);
%%%%%%
%%%% Approx. Average 7 Days
NPPm = NPPI*(n3-1)/n3 + NPPtm1/n3;
Bfac_week = Bfac_weekI*(n3-1)/n3 + Bfac_day/n3;
%%%% Approx. Average 30 days
Tsmm = TdpI*(n2-1)/n2 + Tsm/n2;
%%%%%%%%%%
NupIm(1) = NupI(1)*(n4-1)/n4  + Nuptake/n4;
NupIm(2) = NupI(2)*(n4-1)/n4  + Puptake/n4;
NupIm(3) = NupI(3)*(n4-1)/n4  + Kuptake/n4;
%%%%%%%%%  Average PAR  30 and 45 days
PAR_Im(1)= PAR_I(1)*(n2-1)/n2 + mean(PAR)/n2;
PAR_Im(2)= PAR_I(2)*(n5-1)/n5 + mean(PAR)/n5;
PAR_Im(3)= PAR_I(3)*(n6-1)/n6 + (PAR_Im(1)-PAR_Im(2))/n6;
%%%%%  Average new leaf biomass 30 days
NBL_Im = NBL_I*(n2-1)/n2 + NBLeaftm1/n2;
%%%%%%%
if OPT_EnvLimitGrowth == 1
    [GF]= EnvConGrowth(Ta,Bfac_alloc,gcoef); %%% [gC /m2 day]
else
    GF=Inf; %% [gC /m2 day]
end
%%%%%%%%%%%%%%%%%%%%%
[FNC,e_relN,rNc,rPc,rKc,rMc,rNcR,Navailtm1,Pavailtm1,Kavailtm1]= Nutrients_Available(Btm1,FNCtm1,Stoich,Nreservetm1,Preservetm1,Kreservetm1,OPT_SoilBiogeochemistry);
%%%%%%%%%%%%%%%%%%%%%%
%%% Root Exudation and transfer to Mychorriza
[Rexmy]= Root_Exudation(NPPtm1,Btm1(3),Bam,Bem,Btm1(4),rMc,rNc,rPc,rKc,ParEx,NupIm,NavlI,EM,Tdp_bio); %%% % [gC / m^2 d]
%%%%%%%%%%%%%%%%%%%%
T_SPAN = [0 dtd];  %%%%%%%
%%%%%%%%%%%%%%%%%%%%%
sol=ode45(@VEGETATION_DYNAMIC,T_SPAN,Btm1,OPT_VD,...
    Tam,Tsm,An,Rdark,Bfac_day,Bfac_alloc,FNC,Se_bio,Tdp_bio,dtd,GF,...
    Sl,mSl,Sl_emecrop,Stoich,r,rNcR,gR,aSE,Trr,dd_max,dc_C,Tcold,drn,dsn,age_cr,PHE_Stm1,AgeLtm1,AgeDLtm1,LtR,eps_ac,...
    Mf,Wm,fab,fbe,Klf,ff_r,dmg,sum(Rexmy),NBLeaftm1,dflotm1,Navailtm1,Pavailtm1,Kavailtm1,soCrop,TBio,GirdOpt,OPT_EnvLimitGrowth,OPT_VCA);
%%%%%%%%%%%%%%%%%%
B = deval(sol,dtd);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[LAI,NPP,Rg,RA,Rms,Rmr,Rmc,ANPP,LAIdead,Sr,Slf,Sfr,Swm,Sll,NLeaf,NLeafdead,NBLeaf,Nreserve,Preserve,Kreserve]= VEG_DYN_RES(B,dtd,Btm1,Tam,Tsm,An,Rdark,...
    Sl,mSl,Sl_emecrop,Stoich,r,gR,aSE,AgeLtm1,AgeDLtm1,age_cr,dc_C,Tcold,Bfac_day,GF,dd_max,PHE_Stm1,dsn,drn,fab,fbe,Wm,Mf,Klf,NBLeaftm1,dmg,...
    Nreservetm1,Preservetm1,Kreservetm1,Nuptake,Puptake,Kuptake,rNcR,rNc,rPc,rKc,OPT_EnvLimitGrowth);
B(8)=0; %% Re-initialize B(8)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[PHE_S,dflo,AgeL,AgeDL]= PHENOLOGY_STATE(NLeaf,AgeLtm1,dtd,...
    LAIdead,NLeafdead,AgeDLtm1,...
    PHE_Stm1,LAI,aSE,age_cr,jDay,Tsmm,Bfac_day,Bfac_week,NPPm,PAR_Im(3),L_day,Bfac_lo,Bfac_ls,Tlo,Tls,mjDay,LDay_min,LDay_cr,dflotm1,dmg,PAR_th,LAI_min,jDay_cut,LAI_cut);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[e_rel]= RELATIVE_PC(AgeL,dflo,NBL_Im,B(1),age_cr,aSE,L_day,Lmax_day,jDay);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
RB(1:7)=0; %% Removed Live Leaves/ Sapwood/ Fine Roots /Carbohydrate Reserve /Fruit and Flower /Heartwood - Dead Sapwood /Standing Dead Leaves
%%%%%%%%%% Grass Cut and Grazing
if (aSE==2)  &&  not(isempty(intersect(jDay,jDay_cut)))
    [B,LAI,LAIdead,RB([1,7])]= Grass_Cut_Grazing(B,LAI,LAIdead,dtd,aSE,Sl,mSl,jDay,jDay_cut,LAI_cut) ;
end
%%%%%%%%%% Fruit Harvest
if not(isempty(intersect(jDay,jDay_harv)))
    [B,RB(5)]= Fruit_Harvest(B,dtd,jDay,jDay_harv,B_harv);
end
%%%%%%%%%% Forest Logging 
[B,RB,LAI,LAIdead,ManI]= Forest_Logging_Fire(B,RB,dtd,Sl,mSl,aSE,LAI,LAIdead,Datam,Mpar);
%%%%% Crop Management 
if (aSE==5)
    [B,RB,LAI,LAIdead,ManI,PHE_S,dflo,AgeL,AgeDL]= Crop_Planting_Harvest(B,RB,dtd,LAI,LAIdead,PHE_S,dflo,AgeL,AgeDL,Datam,Mpar);
end
%%%%%
%%%%%%%%% LAI Below Minimum
LAI(LAI<LAI_min)=0;
%%%  Root Properties
rroot =  0.5*1e-3 ;
[Rootl]=Root_properties(B(3),rroot);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%