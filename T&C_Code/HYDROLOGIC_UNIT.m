%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% HYDRO BUDGET FOR A UNIT  -DECOUPLED            %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[V,Vice,O,Oice,ZWT,OF,OS,OH,OL,Psi_s_H,Psi_s_L,Rd,Qi_out,WTR,...
    Rh,Lk,f,WIS,Ts,Pr_sno,Pr_liq,Csno,Cice,NDVI,rb_H,rb_L,rs_sunH,...
    rs_sunL,rs_shdH,rs_shdL,r_litter,...
    An_L,An_H,Rdark_L,Rdark_H,Ci_sunH,Ci_sunL,Ci_shdH,Ci_shdL,...
    rap_H,rap_L,r_soil,b_soil,alp_soil,ra,Rn,...
    H,QE,Qv,Lpho,T_H,T_L,EIn_H,EIn_L,EG,ESN,ESN_In,ELitter,EWAT,EICE,EIn_urb,EIn_rock,dw_SNO,...
    G,Gfin,Tdp,Tdpsnow,Tdeb,Tdamp,Tice,Tdp_H,Tdp_L,SWE,SND,ros,In_SWE,SP_wc,WR_SP,U_SWE,NIn_SWE,dQ,Qfm,t_sls,DQ,DT,...
    WAT,ICE,ICE_D,IP_wc,WR_IP,NIce,Cicew,Csnow,FROCK,Imelt,Smelt,...
    In_H,In_L,In_Litter,In_urb,In_rock,Dr_H,Dr_L,SE_rock,SE_urb,Lk_wat,Lk_rock,er,...
    gsr_H,Psi_x_H,Psi_l_H,Jsx_H,Jxl_H,Kleaf_H,Kx_H,Vx_H,Vl_H,...
    gsr_L,Psi_x_L,Psi_l_L,Jsx_L,Jxl_L,Kleaf_L,Kx_L,Vx_L,Vl_L,...
    fapar_H,fapar_L,SIF_H,SIF_L,...
    snow_alb,tau_sno,e_sno,Ws_under,dQVEG,HV,QEV,TsV,Ts_under,EK,POT]=HYDROLOGIC_UNIT(Vtm1,Oicetm1,aR,Zs,...
    EvL_Zs,Inf_Zs,Zinf,RfH_Zs,RfL_Zs,dz,Dz,ms,Kbot,Pr,Ta,Ds,Ws,zatm,Tstm1,Tstm1_under,IrD,dt,dth,ea,N,Pre,Tstm0,...
    LAI_H,SAI_H,LAI_L,SAI_L,LAIdead_H,LAIdead_L,Rrootl_H,Rrootl_L,BLit,Sllit,Kct,...
    Datam,DeltaGMT,Lon,Lat,t_bef,t_aft,...
    Ccrown,Cbare,Crock,Curb,Cwat,...
    Soil_Param,Interc_Param,SnowIce_Param,VegH_Param,VegL_Param,...
    Zs_deb,Deb_Par,...
    ZR95_H,ZR95_L,...
    SAB1,SAB2,SAD1,SAD2,PARB,PARD,SvF,SNDtm1,snow_albtm1,Color_Class,OM_H,OM_L,...
    PFT_opt_H,PFT_opt_L,hc_H,hc_L,d_leaf_H,d_leaf_L,...
    Ca,Oa,Citm1_sunH,Citm1_shdH,Citm1_sunL,Citm1_shdL,...
    e_rel_H,e_relN_H,e_rel_L,e_relN_L,...
    e_snotm1,In_Htm1,In_Ltm1,In_Littertm1,In_urbtm1,In_rocktm1,SWEtm1,In_SWEtm1,....
    Tdebtm1,Tdptm1,Tdpsnowtm1,Tdamptm1,Ticetm1,...
    WATtm1,ICEtm1,IP_wctm1,ICE_Dtm1,Cicewtm1,...
    Vxtm1_H,Vltm1_H,Vxtm1_L,Vltm1_L,Psi_xtm1_H,Psi_ltm1_H,Psi_xtm1_L,Psi_ltm1_L,...
    FROCKtm1,Krock,Ws_undertm1,...
    Tdew,t_slstm1,rostm1,SP_wctm1,fpr,Pr_sno_day,...
    Urb_Par,In_max_urb,In_max_rock,K_usle,tau_snotm1,Ta_day,Slo_top,Slo_pot,Asur,Ared,aTop,EKtm1,q_runon,Qi_in,...
    pow_dis,a_dis,Salt,...
    SPAR,SN,OPT_min_SPD,OPT_VegSnow,OPT_SoilTemp,OPT_PlantHydr,Opt_CR,Opt_ST,Opt_ST2,OPT_SM,OPT_STh,OPT_FR_SOIL,OPT_PH)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUTS
%%% INPUTS  /
%%% Conversion to Asur --- Computation in the n-coordinate
aTop= Asur*aTop/Ared; %% Area/Width for Asur accounting for rock content areal reduction (Ared)
Pr=Pr/Asur; %% [mm/h] su Asur
q_runon = q_runon/Asur; %% [mm/h] su Asur
Qi_in = Qi_in/Asur ; %% [mm/h] su Asur
In_Htm1=In_Htm1/Asur;%% [mm] su Asur
In_Ltm1=In_Ltm1/Asur;% [mm] su Asur
In_urbtm1=In_urbtm1/Asur;% [mm] su Asur
In_rocktm1=In_rocktm1/Asur;% [mm] su Asur
In_Littertm1  = In_Littertm1/Asur; %%[mm] su Asur
SWEtm1=SWEtm1/Asur;% [mm] su Asur
In_SWEtm1=In_SWEtm1/Asur;% [mm] su Asur
SP_wctm1=SP_wctm1/Asur;% [mm] su Asur
In_max_urb=In_max_urb/Asur;% [mm] su Asur
In_max_rock=In_max_rock/Asur;% [mm] su Asur
WATtm1=WATtm1/Asur;% [mm] su Asur
ICEtm1=ICEtm1/Asur;% [mm] su Asur
IP_wctm1=IP_wctm1/Asur;% [mm] su Asur
FROCKtm1=FROCKtm1/Asur;% [mm] su Asur
%%%% Topographic slopes sine - cosine
cosalp=cos(atan(Slo_top));  %% [-]
sinalp=sin(atan(Slo_top));  %% [-]
%%%%%%%% Water Logging
ydepth = 0.001*q_runon*dth; %%[m]  su Asur
%%%%
Pr_sno_day = sum(Pr_sno_day(2:end)*dth); %%[mm]
%%%% 
Ccrown_L=sum((ZR95_L>0).*Ccrown); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Soil_Parameters
Ks_Zs=Soil_Param.Ks_Zs;
Osat=Soil_Param.Osat;
Ohy =Soil_Param.Ohy;
L=Soil_Param.L;
Pe=Soil_Param.Pe;
O33=Soil_Param.O33;
alpVG=Soil_Param.alpVG;
nVG=Soil_Param.nVG;
Ks_mac=Soil_Param.Ks_mac;
Omac=Soil_Param.Omac;
alpVGM=Soil_Param.alpVGM;
nVGM =Soil_Param.nVGM;
rsd =Soil_Param.rsd;
lan_dry =Soil_Param.lan_dry;
lan_s =Soil_Param.lan_s;
cv_s =Soil_Param.cv_s;
s_SVG=Soil_Param.s_SVG;
bVG=Soil_Param.bVG;
lVG=Soil_Param.lVG;
lVGM=Soil_Param.lVGM;
Phy = Soil_Param.Phy; 
%%%%%%%%%% Interception Parameters
Sp_SN_In=Interc_Param.Sp_SN_In;
Sp_LAI_L_In=Interc_Param.Sp_LAI_L_In;
Sp_LAI_H_In=Interc_Param.Sp_LAI_H_In;
gcI=Interc_Param.gcI;
KcI=Interc_Param.KcI;
%%%%%% Snow-Ice Parameters
TminS=SnowIce_Param.TminS;
TmaxS=SnowIce_Param.TmaxS;
WatFreez_Th=SnowIce_Param.WatFreez_Th;
dz_ice=SnowIce_Param.dz_ice;
Th_Pr_sno=SnowIce_Param.Th_Pr_sno;
ros_max1=SnowIce_Param.ros_max1;
ros_max2=SnowIce_Param.ros_max2;
Ice_wc_sp=SnowIce_Param.Ice_wc_sp;
ros_Ice_thr=SnowIce_Param.ros_Ice_thr;
Aice=SnowIce_Param.Aice;
%%%%% Vegetation High Parameters
KnitH=VegH_Param.KnitH;
mSl_H=VegH_Param.mSl_H;
Sl_H=VegH_Param.Sl_H;
CT_H =VegH_Param.CT_H;
Vmax_H=VegH_Param.Vmax_H;
FI_H=VegH_Param.FI_H;
a1_H=VegH_Param.a1_H;
go_H=VegH_Param.go_H;
DSE_H=VegH_Param.DSE_H;
Ha_H=VegH_Param.Ha_H;
Do_H=VegH_Param.Do_H; 
gmes_H=VegH_Param.gmes_H;
rjv_H=VegH_Param.rjv_H;
Psi_sto_50_H=VegH_Param.Psi_sto_50_H;
Psi_sto_00_H=VegH_Param.Psi_sto_00_H;
Axyl_H=VegH_Param.Axyl_H;
PsiL50_H=VegH_Param.PsiL50_H;
PsiL00_H=VegH_Param.PsiL00_H;
Kleaf_max_H=VegH_Param.Kleaf_max_H;
Cl_H=VegH_Param.Cl_H;
Kx_max_H=VegH_Param.Kx_max_H;
PsiX50_H=VegH_Param.PsiX50_H;
Cx_H=VegH_Param.Cx_H;
Osm_reg_Max_H = VegH_Param.Osm_reg_Max_H;
eps_root_base_H = VegH_Param.eps_root_base_H; 
%%%%% Vegetation Low Parameters
KnitL=VegL_Param.KnitL;
mSl_L=VegL_Param.mSl_L;
Sl_L =VegL_Param.Sl_L;
CT_L=VegL_Param.CT_L;
Vmax_L =VegL_Param.Vmax_L;
FI_L =VegL_Param.FI_L;
a1_L=VegL_Param.a1_L;
go_L=VegL_Param.go_L;
DSE_L=VegL_Param.DSE_L;
Ha_L=VegL_Param.Ha_L;
Do_L=VegL_Param.Do_L;
gmes_L=VegL_Param.gmes_L;
rjv_L=VegL_Param.rjv_L;
Psi_sto_50_L=VegL_Param.Psi_sto_50_L;
Psi_sto_00_L=VegL_Param.Psi_sto_00_L;
Axyl_L=VegL_Param.Axyl_L;
PsiL50_L=VegL_Param.PsiL50_L;
PsiL00_L=VegL_Param.PsiL00_L;
Kleaf_max_L=VegL_Param.Kleaf_max_L;
Cl_L=VegL_Param.Cl_L;
Kx_max_L=VegL_Param.Kx_max_L;
PsiX50_L=VegL_Param.PsiX50_L;
Cx_L=VegL_Param.Cx_L;
Osm_reg_Max_L = VegL_Param.Osm_reg_Max_L;
eps_root_base_L = VegL_Param.eps_root_base_L; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%  Value of the Volume  -- V
%%% Initial Condition --->  Vtm1
%%%%%%%%
[Otm1,ZWTtm1,OFtm1,OStm1,Psi_stm1_H,Psi_stm1_L,gsr_Htm1,gsr_Ltm1,Exwat_Htm1,Exwat_Ltm1]=Soil_Water_MultiLayer(Vtm1,Zs,...
    dz,ms,Ccrown,Osat,Ohy,nVG,alpVG,lVG,Ks_Zs,L,Pe,O33,Ks_mac,Omac,alpVGM,nVGM,lVGM,s_SVG,bVG,Phy,SPAR,EvL_Zs,Inf_Zs,RfH_Zs,RfL_Zs,...
    Rrootl_H,Rrootl_L,PsiL50_H,PsiL50_L,PsiX50_H,PsiX50_L,Tstm1,Tdptm1,Psi_sto_00_H,Psi_sto_50_H,Psi_sto_00_L,Psi_sto_50_L,...
    Salt,Osm_reg_Max_H,Osm_reg_Max_L,eps_root_base_H,eps_root_base_L);
%%%%%%%%%%%%%%%
%%%%%% Thermal properties soil
if Crock ==1 || Curb ==1  || Cwat ==1
    [lan_oth,cv_oth,CTt_oth]=Other_Thermal_properties(Cwat,Curb,Crock,ms);
    CTt=[];
else
    [~,~,CTt]=Soil_Thermal_properties_FT(Tdptm1,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
        Phy,s_SVG,bVG,Osat,Ohy,(Otm1+Oicetm1),OPT_FR_SOIL);
    lan_oth=[];cv_oth=[];CTt_oth=[]; 
end
%%%%%
%%%%% Incoming Longwave  %%%%%%%%%%%
[Latm,N]=Incoming_Longwave(Ta,ea,N); % Latm Incoming LongWave Radiation [W/m^2]
%%%%% Sun height
[h_S]=SetSunVariables(Datam,DeltaGMT,Lon,Lat,t_bef,t_aft);
%%%%%%%%%%%%%%%%%%%%%%%%% Account for the water logging
% >--- OS [ Superficial Soil ]
[r_soil,b_soil,alp_soil]=Soil_Resistence(Tdptm1(1),Pre,Ws_undertm1,ea,q_runon,OStm1,Ks_Zs(1),Osat(1),Ohy(1),L(1),Pe(1),O33(1),alpVG(1),nVG(1),lVG(1),...
    Ks_mac(1),Omac(1),alpVGM(1),nVGM(1),lVGM(1),Phy,s_SVG(1),bVG(1),SPAR); %%
%%% Theoretically Ks_Zs should be reduced when soil is frozen 
%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%
%%%%%% RADIATION INPUT MANAGEMENT
Rsw.dir_vis = SAB1;
Rsw.dir_nir = SAB2;
Rsw.dif_vis = SAD1;
Rsw.dif_nir = SAD2;
PAR.dir = PARB;
PAR.dif = PARD;
clear SAB1 SAB2 SAD1 SAD2 PARB PARD
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% SVAT MANAGER
[Ts,Pr_sno,Pr_liq,Ws_under,Csno,Cice,Cfol_H,Cfol_L,CLitter,NDVI,rb_H,rb_L,rs_sunH,rs_sunL,rs_shdH,rs_shdL,r_litter,...
    An_L,An_H,Rdark_L,Rdark_H,Ci_sunH,Ci_sunL,Ci_shdH,Ci_shdL,rap_H,rap_L,ra,Rn,...
    H,QE,Qv,Lpho,T_H,T_L,EIn_H,EIn_L,EG,ESN,ESN_In,EWAT,EICE,ELitter,EIn_urb,EIn_rock,dw_SNO,...
    G,Gfin,Tdp,Tdpsnow,Oint,Oice,Tdamp,Tdeb,Tice,SWE,SND,ros,In_SWE,SP_wc,WR_SP,U_SWE,NIn_SWE,dQ,Qfm,t_sls,DQ,DT,...
    In_H,In_L,In_Litter,In_urb,In_rock,Dr_H,Dr_L,SE_rock,SE_urb,WIS,...
    ICE,ICE_D,IP_wc,WR_IP,NIce,Cicew,Csnow,...
    snow_alb,tau_sno,e_sno,...
    fapar_H,fapar_L,SIF_H,SIF_L,...
    gsr_H,Psi_x_H,Psi_l_H,Jsx_H,Jxl_H,Kleaf_H,Kx_H,Vx_H,Vl_H,...
    gsr_L,Psi_x_L,Psi_l_L,Jsx_L,Jxl_L,Kleaf_L,Kx_L,Vx_L,Vl_L,...
    dQVEG,HV,QEV,TsV,Ts_under,WAT,q_run,EG_dis,J_Hdis,J_Ldis,Imelt,Smelt]=SVAT_UNIT(Tstm0,Pr,Ta,Ds,Ws,zatm,Tstm1,Tstm1_under,dt,ea,Latm,N,Pre,...
    OStm1,Psi_stm1_H,Psi_stm1_L,gsr_Htm1,gsr_Ltm1,Exwat_Htm1,Exwat_Ltm1,....
    LAI_H,SAI_H,LAI_L,SAI_L,LAIdead_H,LAIdead_L,BLit,Sllit,Kct,TminS,TmaxS,...
    Sp_SN_In,Sp_LAI_L_In,Sp_LAI_H_In,h_S,...
    Ccrown,Cbare,Crock,Curb,Cwat,ydepth,Ccrown_L,...
    Rsw,PAR,SvF,SNDtm1,snow_albtm1,Color_Class,OM_H,OM_L,...
    PFT_opt_H,PFT_opt_L,hc_H,hc_L,d_leaf_H,d_leaf_L,...
    r_soil,b_soil,alp_soil,gcI,KcI,...
    KnitH,KnitL,mSl_H,Sl_H,mSl_L,Sl_L,Ca,Oa,Citm1_sunH,Citm1_shdH,CT_H,Citm1_sunL,Citm1_shdL,CT_L,Vmax_H,Vmax_L,FI_H,FI_L,a1_H,go_H,a1_L,go_L,...
    DSE_H,Ha_H,Do_H,e_rel_H,e_relN_H,DSE_L,Ha_L,Do_L,e_rel_L,e_relN_L,gmes_H,rjv_H,gmes_L,rjv_L,...
    e_snotm1,In_Htm1,In_Ltm1,In_Littertm1,In_urbtm1,In_rocktm1,SWEtm1,In_SWEtm1,....
    Oicetm1,Otm1,Tdptm1,Tdpsnowtm1,Tdamptm1,Ticetm1,dz,Zs,ms,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,Phy,s_SVG,bVG,Osat,Ohy,CTt,...
    lan_oth,cv_oth,CTt_oth,Zs_deb,Deb_Par,Tdebtm1,...
    Tdew,t_slstm1,rostm1,SP_wctm1,fpr,WatFreez_Th,dz_ice,...
    Pr_sno_day,Th_Pr_sno,ros_max1,ros_max2,...
    Urb_Par,In_max_urb,In_max_rock,tau_snotm1,Ta_day,Ws_undertm1,...
    EvL_Zs,RfH_Zs,RfL_Zs,Vtm1,WATtm1,...
    ICEtm1,IP_wctm1,ICE_Dtm1,Ice_wc_sp,ros_Ice_thr,Aice,Cicewtm1,...
    Vxtm1_H,Vltm1_H,Vxtm1_L,Vltm1_L,Psi_xtm1_H,Psi_ltm1_H,Psi_xtm1_L,Psi_ltm1_L,...
    Psi_sto_50_H,Psi_sto_00_H,ZR95_H,...
    Psi_sto_50_L,Psi_sto_00_L,ZR95_L,...
    Axyl_H,PsiL50_H,PsiL00_H,Kleaf_max_H,Cl_H,Kx_max_H,PsiX50_H,Cx_H,...
    Axyl_L,PsiL50_L,PsiL00_L,Kleaf_max_L,Cl_L,Kx_max_L,PsiX50_L,Cx_L,...
    OPT_min_SPD,OPT_VegSnow,OPT_SoilTemp,OPT_PlantHydr,Opt_CR,Opt_ST,Opt_ST2,OPT_STh,OPT_FR_SOIL,OPT_PH);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% All of the variable over Asur (normal coordinates) %%%%%%%%%%%%%%%%%%%%%%%%
%%% Pr_sno Pr_liq An_L An_H Rdark_H Rdark_L Rn H QE Qv T_H T_L EIn_H EIn_L
%%% ELitter WIS In_Litter
%%% EG ESN ESN_In EWAT EICE EIn_urb EIn_rock G SWE SND In_SWE SP_wc WR_SP Qfm
%%% In_H In_L In_urb In_rock Dr_H Dr_L ICE ICE_D IP_wc WR_IP NIce EG_dis THdis T_Ldis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if OPT_FR_SOIL == 1
    fQ= (Oice)./(Oice+Oint+1e-5-Ohy); %%[-]
    Kfro = 10.^(-5.*fQ).*Ks_Zs; %% Hansson et al 2004
    Ks_Zs=Kfro ;
end
%%%%%%%%%%%%%
q_runon = q_run;  %% [mm/h]
ydepth = q_runon*dth; %% [mm]  su Asur
WATtm1 = WAT;   %% [mm] %% It changes only in case of freezing lakes or in case Cwat < 1 
%%%%%%
WIS=WIS+IrD*dth; %%% Drip Irrigation added;; % [mm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(Crock ==1 || Curb ==1  || Cwat ==1)
    %%%%%%%%% NORMAL CASE OF SOIL %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  isnan(a_dis(1,1)) ||  Pr_liq == 0
        %%% Erosion
        [er,ke]=Erosion(dt,Pr_liq,hc_H,hc_L,K_usle,Ccrown,Cbare,Csno,Cfol_H,Cfol_L,CLitter,Dr_H,Dr_L); %%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Soil Sealing and Crust
        if (ke == 0) || isnan(a_dis(1,1));
            EK = 0; KsC=Ks_Zs(1); LC=L(1); OsatC=Osat(1); OhyC=Ohy(1); PeC=Pe(1);
        else
            [KsC,LC,OsatC,OhyC,PeC,EK]=SoilCrust(dth,ke,EKtm1,rsd(1),Osat(1),Ohy(1),L(1),Pe(1),Ks_Zs(1));
        end
        %%%% WIS = [mm] Water Incoming to First Soil Layer
        WIS= WIS/dth; %%% [mm/h]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [f]=Infiltration_2(OsatC,OhyC,LC,alpVG(1),nVG(1),lVG(1),PeC,KsC,O33(1),Ks_mac(1),Omac(1),alpVGM(1),nVGM(1),lVGM(1),Phy,s_SVG(1),bVG(1),SPAR,OStm1,Zinf,WIS,cosalp,ydepth);  %% [mm/h]  %% <--- OS First Layer
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Rh = (WIS-f)*dth; %%% [mm] Horton Runoff
        WIS=WIS*dth; %% [mm]
    else
        %%%% CASE with 5 min rainfall
        %%%%%%%% Rainfall Disaggregation
        dti=300; %%[s]
        dtih= dti/dt; %%[h]
        [Pr_liq_dis]=Disaggregator_Pr(Pr_liq*dth,pow_dis,a_dis);
        WIScon = WIS-sum(Pr_liq_dis); %% [mm]
        if WIScon < 0
            Pr_liq_dis = Pr_liq_dis*WIS./(sum(Pr_liq_dis)); % [mm]
            WIScon =0 ;
        end
        Pr_liq_dis=Pr_liq_dis/dtih; %%% disaggregated precipitation [mm/h]
        %%%%
        er=zeros(1,length(Pr_liq_dis));
        ke=zeros(1,length(Pr_liq_dis));
        f=zeros(1,length(Pr_liq_dis));
        WISvar=zeros(1,length(Pr_liq_dis));
        for jk=1:length(Pr_liq_dis)
            %%%%%% Erosion
            [er(jk),ke(jk)]=Erosion(dti,Pr_liq_dis(jk),hc_H,hc_L,K_usle,Ccrown,Cbare,Csno,Cfol_H,Cfol_L,CLitter,Dr_H,Dr_L);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Soil Sealing and Crust
            [KsC,LC,OsatC,OhyC,PeC,EK]=SoilCrust(dtih,ke(jk),EKtm1,rsd(1),Osat(1),Ohy(1),L(1),Pe(1),Ks_Zs(1));
            EKtm1=EK;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%% WIS = Water Incoming to First Soil Layer
            WISvar(jk)= Pr_liq_dis(jk) + WIScon/dth ; %%% [mm/h]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%% Only for SPAR == 2 - Saxton and  Rawls 
            [f(jk)]=Infiltration_2(OsatC,OhyC,LC,0,0,0,PeC,KsC,O33(1),Ks_mac(1),Omac(1),alpVGM(1),nVGM(1),lVGM(1),Phy,s_SVG(1),bVG(1),2,OStm1,Zinf,WISvar(jk),cosalp,ydepth);  %% [mm/h]  %% <--- OS First Layer
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        Rh = (WISvar-f)*dtih; %%% [mm] Horton Runoff
        %%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%
        Rh=sum(Rh); % [mm]
        f=sum(f*dtih)/dth;  % [mm/h]
        %WIS=sum(WISvar*dtih); %[mm]
        er=sum(er*dtih)/dth; %%  [kg/h m^2]
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Soil Bottom Leakage
    [Lk]=Leakage_Bottom(Otm1,Ks_Zs,Osat,Ohy,L,alpVG,nVG,lVG,Ks_mac,Omac,alpVGM,nVGM,lVGM,Kbot,ms,SPAR); %%% [mm/h]
    Lk=Lk*Ared; %% due to less available area
    %%%%%%%%%%
else
    %%%% No SOIL
    Rh = 0;  f = 0;  Lk = 0;  er = 0;  EK = 0; %WIS = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%  WATER
if Cwat > 0
    WAT = WATtm1 - EWAT*dth + (WR_SP*(1-Cice)+WR_IP)*Cwat + Cwat*dth*(Pr_liq*(1-Csnow)*(1-Cicew) + Pr_sno*(1-Cicew)*(1-Csnow)) + q_runon*dth*Cwat  + sum(Qi_in)*dth*Cwat; %%[mm]
    [Lk_wat]=Leakage_Rock(Krock,WAT,dth); %% [mm]
    WAT=WAT-Lk_wat;
else
    Lk_wat = 0;
    WAT = 0;
end
%%%%% URBAN
if Curb > 0
    [Lk_urb]=Leakage_Rock(Krock,In_urb,dth); %% [mm]
    In_urb=In_urb-Lk_urb;
else
    Lk_urb = 0;
end
%%%%%%%%%% ROCK  %%%%%%%%%%%%%%%%%
[Lk_rock]=Leakage_Rock(Krock,In_rock,dth); %% [mm/h]
In_rock = In_rock - Lk_rock; %%[mm]
%%% Fracture Rocks Storage
FROCK = FROCKtm1 + Lk + Lk_wat + Lk_urb + Lk_rock + Crock*sum(Qi_in)*dth; %%[mm]
%%%% Total Horton Runoff
Rh = Rh + SE_rock + SE_urb; %%[mm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Simplify Richards Model
%%% Initial Condition --->  V- Values
if Crock ==1 || Curb ==1 || Cwat ==1
    if OPT_FR_SOIL == 1
        V=(Oint-Ohy).*dz;
        Vice =(Oice).*dz;
        Rd_cryo=0; 
    else
        V=Vtm1;
        Vice=V*0;
        Rd_cryo=0; 
    end
    
else
    %%%%%%%%%%%%%%%%%%%
    %%% V over Asur and Ared (for rock content in the soil)
    %%% Operations over Ared
    f=f/Ared; Qi_in=Qi_in/Ared;  Lk=Lk/Ared;
    EG_dis=EG_dis/Ared; J_Hdis=J_Hdis/Ared; J_Ldis=J_Ldis/Ared;
    %%%%%%%%%%%%%%%%%%%
    V0 =Vtm1;
    T_SPAN = [0 dth]; 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if OPT_FR_SOIL == 1
        %%% Potentially reduced water content because of freezing 
        V0 = (Oint-Ohy).*dz; %%% [mm]
    else
        Oice = 0*Oice; 
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Param.Osat = Osat; Param.Ohy = Ohy; Param.O33 = O33;
    %Param.dz = dz; Param.Ks_Zs =Ks_Zs; Param.Dz = Dz; Param.ms = ms;
    %Param.L = L; Param.Pe = Pe; Param.aR =aR; Param.aTop = aTop;
    %Param.alpVG=alpVG; Param.nVG = nVG; Param.Zs=Zs;
    %Param.cosalp=cosalp; Param.sinalp=sinalp; Param.SN=SN;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ISeep =  ones(1,ms) ; %% Otm1 >= Osat-0.005 ;  %%% Decision to Seep
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [Tout,Vout]=ode23s(@SOIL_MOISTURES_RICH_COMP,T_SPAN,V0,OPT_SM,...
        Lk,f,EG_dis,sum(J_Hdis,1),sum(J_Ldis,1),...
        Qi_in,Slo_pot,ISeep,int32(SPAR),Osat,Ohy,O33,dz,Ks_Zs,Dz,int32(ms),L,Pe,aR,aTop,alpVG,nVG,lVG,...
        Ks_mac,Omac,alpVGM,nVGM,lVGM,Phy,s_SVG,bVG,Zs,cosalp,sinalp,int32(SN));
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %    [Tout,Vout]=ode23s(@SOIL_MOISTURES_RICH,T_SPAN,V0,OPT_SM,...
    %    Lk,f,EG_dis,sum(J_Hdis,1),sum(J_Ldis,1),...
    %    Qi_in,Slo_pot,ISeep,SPAR,Param);
    %%%%%%%%%%%%%%%%%%%%%%%%
    V = Vout(end,:); % [mm/h] su Asur
    if OPT_FR_SOIL == 1
        [V,~,Rd_cryo] = Cryosuction_stabilizer(Oice,V,dz,Osat,Ohy);
        Vice = (Oice).*dz;
    else
        Vice = (Oice).*dz;
        Rd_cryo=0; 
    end
        %%%%%%%%%%%%%%%%%%%%%%
    if isnan(sum(V))
        %%%%%%%%%%%%%%%
        disp('NaN values in the Volumes')
        return
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[O,ZWT,OF,OS,Psi_s_H,Psi_s_L,gsr_H,gsr_L,Exwat_H,Exwat_L,Rd,WTR,POT,OH,OL]=Soil_Water_MultiLayer(V,Zs,...
    dz,ms,Ccrown,Osat,Ohy,nVG,alpVG,lVG,Ks_Zs,L,Pe,O33,Ks_mac,Omac,alpVGM,nVGM,lVGM,s_SVG,bVG,Phy,SPAR,EvL_Zs,Inf_Zs,RfH_Zs,RfL_Zs,...
    Rrootl_H,Rrootl_L,PsiL50_H,PsiL50_L,PsiX50_H,PsiX50_L,Ts,Tdp,Psi_sto_00_H,Psi_sto_50_H,Psi_sto_00_L,Psi_sto_50_L,...
    Salt,Osm_reg_Max_H,Osm_reg_Max_L,eps_root_base_H,eps_root_base_L);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% 
[Tdp_H,Tdp_L]=RootZone_Temp(Tdp,RfH_Zs,RfL_Zs);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Qi_out, [mm] Lateral Flow outgoing  [1....m]
%Qi_in [mm] Lateral Flow incoming    [1....m]
Qi_out=zeros(1,ms);
if not(Crock ==1 || Curb ==1 || Cwat ==1)
    for jk = 2:length(Tout)
        %%%%%%%%%%
        Oor = 0.5*(Vout(jk,:)+Vout(jk-1,:))./dz + Ohy ; %%[-]
        I1 = Oor >= Osat-1e-5; I2 = Oor <= Ohy;
        Oor(I1)=Osat(I1); Oor(I2)=Ohy(I2)+1e-9;
        %%%%%%%%
        Ser = (Oor-Ohy)./(Osat-Ohy);
        mVG= 1-1./nVG;
        if SPAR == 1
            Kor= aR.*Ks_Zs.*((Ser).^(0.5)).*(1-(1-(Ser).^(1./mVG)).^mVG).^2; %%% [mm/h]
        else
            if SPAR == 2
                Kor = aR.*Ks_Zs.*(Oor./Osat).^(3+(2./L)); %%% <<<---
            elseif SPAR == 3
                [Kor]=Conductivity_Suction(SPAR,Ks_Zs,Osat,Ohy,L,Pe,O33,alpVG,nVG,lVG,Ks_mac,Omac,alpVGM,nVGM,lVGM,Phy,s_SVG,bVG,Oor);
            end
        end 
        %%%%%%%%%%%%%%%
        [Qi_out_t]=Lateral_Subsurface_Flow(Kor,dz,Slo_pot,aTop,cosalp,sinalp,SN,ISeep);
        Qi_out = Qi_out + Qi_out_t*(Tout(jk)-Tout(jk-1)); %%% [mm]
        %%%%%%%%%%%%%%%
    end
    Qi_out=Qi_out/dth; %% [mm/h]
    %%%%%%%%%%%%%%%  Volume Correction for Rd and WTR
    V(1) = V(1) +  WTR(2) - Rd;
    V(2:end-1)= V(2:end-1)+ (WTR(3:end) - WTR(2:end-1));
    V(end)= V(end) - WTR(end);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%% Volume Compensation -- Negative Value - %%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if sum(V < 0) > 0 
        EG=EG/Ared; Jsx_H=Jsx_H/Ared; Jsx_L=Jsx_L/Ared;
        [V,Jsx_H,Jsx_L,EG,Lk]=Volume_Correction(V,EvL_Zs,RfH_Zs,RfL_Zs,EG,Jsx_H,Jsx_L,Lk);
        EG=EG*Ared; Jsx_H=Jsx_H*Ared; Jsx_L=Jsx_L*Ared;
        %disp('Negative Volumes')
    end
end
%%%%%
Rd=Rd_cryo+Rd; %%% Recombine Saturation Excess Runoff 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% Correction becuase Qi_out is computed with approximations
%%%%%%%%%%%% discrete Tout
if (Slo_top > 0) && (sum(Qi_out) > 0) &&  not(Crock ==1 || Curb ==1 || Cwat ==1)
    dvol_correction = f*dth + sum(Vtm1) - sum(V) + sum(Oicetm1.*dz) - sum(Vice) - (EG/Ared)*dth - Lk*dth ...
        - sum(Qi_out)*dth -Rd -sum(Jsx_L/Ared).*dth -sum(Jsx_H/Ared).*dth  + sum(Qi_in)*dth ; %%[mm]
    Qi_out = Qi_out + (dvol_correction/dth)*Qi_out/sum(Qi_out); %%[mm/h] 
end
%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% Re-Transformation from Ared (rock content)
%Not retransformed --    WTR=WTR*Ared;  Qi_in=Qi_in*Ared;
Rd = Rd*Ared;
Qi_out=Qi_out*Ared;
Lk = Lk*Ared;
f=f*Ared;
%%%% Re-projected to the vertical coordinate
Qi_out = Qi_out*Asur;
Rh = Rh*Asur;
Rd = Rd*Asur;
Lk = Lk*Asur;
Pr_sno = Pr_sno*Asur;
Pr_liq = Pr_liq*Asur;
%Rn = Rn*Asur;
%H = H*Asur;
%QE= QE*Asur;
%Qv = Qv*Asur;
%G= G*Asur;
%Qfm = Qfm*Asur;
T_H=T_H*Asur;
T_L =T_L*Asur;
Jsx_H=Jsx_H*Asur;
Jsx_L =Jsx_L*Asur;
Jxl_H=Jxl_H*Asur;
Jxl_L =Jxl_L*Asur;
EIn_H = EIn_H*Asur;
EIn_L = EIn_L*Asur;
EG = EG*Asur;
ELitter = ELitter*Asur;
ESN= ESN*Asur;
ESN_In = ESN_In*Asur;
EWAT = EWAT*Asur;
EIn_urb = EIn_urb*Asur;
EIn_rock = EIn_rock*Asur;
SWE = SWE*Asur;
In_SWE = In_SWE*Asur;
SP_wc = SP_wc*Asur;
U_SWE=U_SWE*Asur; 
NIn_SWE=NIn_SWE*Asur; 
In_H=In_H*Asur;
In_L=In_L*Asur;
In_Litter=In_Litter*Asur;
In_urb=In_urb*Asur;
In_rock=In_rock*Asur;
Dr_H=Dr_H*Asur;
Dr_L=Dr_L*Asur;
SE_rock=SE_rock*Asur;
SE_urb=SE_urb*Asur;
er=er*Asur;
WAT=WAT*Asur;%
EICE= EICE*Asur;
ICE=ICE*Asur;%
IP_wc = IP_wc*Asur;
FROCK=FROCK*Asur;%
Lk_wat=Lk_wat*Asur;
Lk_rock=Lk_rock*Asur;
Imelt=Imelt*Asur;
Smelt=Smelt*Asur;

%%% Mantained in n-coordinate and accouting for rock content
%V,WTR,Vice 
%%%%%%%%%% Mantained in the projected area n-coordinate
%f,WIS,An_L,An_H,Rdark_L,Rdark_H,SND,WR_SP,dQ,DQ,DT,...
%Rn,H,QE,Qv,Lpho,G,Gfin,Qfm,dQVEG,HV,QEV,
%ICE_D,WR_IP
%Vx_H,Vl_H,Vx_L,Vl_L,...
%SIF_H,SIF_L
end