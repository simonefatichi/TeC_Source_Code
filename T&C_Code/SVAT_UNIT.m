%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% SVAT BUDGET FOR A UNIT              %%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Ts,Pr_sno,Pr_liq,Ws_under,Csno,Cice,Cfol_H,Cfol_L,CLitter,NDVI,rb_H,rb_L,rs_sunH,rs_sunL,rs_shdH,rs_shdL,r_litter,...
    An_L,An_H,Rdark_L,Rdark_H,Ci_sunH,Ci_sunL,Ci_shdH,Ci_shdL,rap_H,rap_L,ra,Rn,...
    H,QE,Qv,Lpho,T_H,T_L,EIn_H,EIn_L,EG,ESN,ESN_In,EWAT,EICE,ELitter,EIn_urb,EIn_rock,dw_SNO,...
    G,Gfin,Tdp,Tdpsnow,O,Oice,Tdamp,Tdeb,Tice,SWE,SND,ros,In_SWE,SP_wc,WR_SP,U_SWE,NIn_SWE,dQ,Qfm,t_sls,DQ,DT,...
    In_H,In_L,In_Litter,In_urb,In_rock,Dr_H,Dr_L,SE_rock,SE_urb,WIS,...
    ICE,ICE_D,IP_wc,WR_IP,NIce,Cicew,Csnow,...
    snow_alb,tau_sno,e_sno,...
    fapar_H,fapar_L,SIF_H,SIF_L,...
    gsr_H,Psi_x_H,Psi_l_H,Jsx_H,Jxl_H,Kleaf_H,Kx_H,Vx_H,Vl_H,...
    gsr_L,Psi_x_L,Psi_l_L,Jsx_L,Jxl_L,Kleaf_L,Kx_L,Vx_L,Vl_L,...
    dQVEG,HV,QEV,TsV,Ts_under,WAT,q_run,EG_dis,J_Hdis,J_Ldis,Imelt,Smelt]=SVAT_UNIT(Tstm0,Pr,Ta,Ds,Ws,zatm,Tstm1,Tstm1_under,dt,ea,Latm,N,Pre,...
    OS,Psi_s_H,Psi_s_L,gsr_H,gsr_L,Exwat_H,Exwat_L,...
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
    min_SPD,OPT_VegSnow,OPT_SoilTemp,OPT_PlantHydr,Opt_CR,Opt_ST,Opt_ST2,OPT_STh,OPT_FR_SOIL,OPT_PH)
%%% INPUTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OUTPUTS
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% SVAT COMPUTATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Initizialization Ccrown variables
Psi_x_H=zeros(1,length(Ccrown)); Psi_x_L=zeros(1,length(Ccrown));
Psi_l_H=zeros(1,length(Ccrown)); Psi_l_L=zeros(1,length(Ccrown));
Jsx_H=zeros(1,length(Ccrown)); Jsx_L=zeros(1,length(Ccrown));
Jxl_H=zeros(1,length(Ccrown)); Jxl_L=zeros(1,length(Ccrown));
Kleaf_H=zeros(1,length(Ccrown)); Kleaf_L=zeros(1,length(Ccrown));
Kx_H=zeros(1,length(Ccrown)); Kx_L=zeros(1,length(Ccrown));
Vx_H=zeros(1,length(Ccrown));  Vx_L=zeros(1,length(Ccrown));
Vl_H=zeros(1,length(Ccrown)); Vl_L=zeros(1,length(Ccrown));
%%%%%%%%%%%%
rb_H=zeros(1,length(Ccrown)); rb_L=zeros(1,length(Ccrown));
rs_sunH=zeros(1,length(Ccrown)); rs_sunL=zeros(1,length(Ccrown));
rs_shdH=zeros(1,length(Ccrown)); rs_shdL=zeros(1,length(Ccrown));
An_H=zeros(1,length(Ccrown)); An_L=zeros(1,length(Ccrown));
Ci_sunH=zeros(1,length(Ccrown)); Ci_sunL=zeros(1,length(Ccrown));
Ci_shdH=zeros(1,length(Ccrown)); Ci_shdL=zeros(1,length(Ccrown));
Rdark_H=zeros(1,length(Ccrown)); Rdark_L=zeros(1,length(Ccrown));
Lpho_H=zeros(1,length(Ccrown)); Lpho_L=zeros(1,length(Ccrown));
SIF_H=zeros(1,length(Ccrown)); SIF_L=zeros(1,length(Ccrown));
rap_H=zeros(1,length(Ccrown)); rap_L=zeros(1,length(Ccrown));
r_litter=zeros(1,length(Ccrown)); alp_litter=zeros(1,length(Ccrown));
Ws_under=zeros(1,length(Ccrown));
%%%%%%%%%%%%%
%%% Partition Pr_sno Pr_liq
[Pr_sno,Pr_liq]=Precipitation_partition(Pr,Ta,TminS,TmaxS,ea,Pre);
%%% Throughfall
[Cfol_H,Cfol_L,CLitter]=Throughfall(LAI_H,(SAI_H+LAIdead_H),LAI_L,(SAI_L+LAIdead_L),Sllit*BLit,Kct);
%%% Warm/Cold Hydrology
if Pr_sno > 0 || SWEtm1 > 0  || In_SWEtm1 > 0
    Csno = 1;
else
    Csno = 0;
end
%%% Glacied or not Glacied
if ICEtm1 > 0 || rostm1 >= ros_Ice_thr
    Cice = 1;
else
    Cice = 0;
end
%%%% Debris covered glacier or not
if length(Zs_deb)>2
    Cdeb = 1;
    if Cice ~= 1
        %disp('Error Debris without Ice')
        %%% Removing the debris as the ice is completely melted
        Zs_deb=[];
        Cdeb=0;
    end
else
    Cdeb = 0;
end
%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Freezing of water surfaces and snow above frozen water
if  (Cwat > 0 &&  WATtm1 > 0) && ((Tstm1 <=  WatFreez_Th) || (ICEtm1 >0))
    Cicew = 1; Cice = 1;
    if  Pr_sno > 0 || SWEtm1 > 0  || In_SWEtm1 > 0
        Csnow = 1 ; Csno = 1 ;
    else
        Csnow = 0;
    end
else
    Cicew = 0;
    Csnow = 0;
end
%%%%%%%%%%%%%%%%% 
%%% To prevent numerical instability
if  abs(Tstm0-Ta)>25
    Tstm0 = Tstm1;
end
if Tstm0 > 85 || Tstm0 < -100%  maximum initial point for surface temperature 85°C (Aminzadeh et al. 2023)
    Tstm0 = Ta-0.1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[In_max_SWE,In_max_L,In_max_H]=Maximum_Interception(Ccrown,LAI_L,LAI_H,(SAI_H+LAIdead_H),(SAI_L+LAIdead_L),...
    Ta,Sp_SN_In,Sp_LAI_L_In,Sp_LAI_H_In);
%%%%%%%%%%%%%%% Fraction of foliage wet or snow covered
dw_H=  min(1,(In_Htm1./In_max_H).^(2/3)).*(In_max_H>0); %% Wet vegetation Fraction First Layer  --- [Deardorff (1978)]
dw_L=  min(1,(In_Ltm1./In_max_L).^(2/3)).*(In_max_L>0); %% Wet vegetation Fraction Second Layer
dw_SNO =  min(1,In_SWEtm1/In_max_SWE)*(In_max_SWE>0); %% Snow Cover Fraction First Layer --- [Lee and Mahrt 2004]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% ShortWave Evalution %%%%%%%%%%%%%
[RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,PAR_sun_H,PAR_shd_H,...
    RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,PAR_sun_L,PAR_shd_L,FsunH,FshdH,...
    FsunL,FshdL,Kopt_H,Kopt_L,fapar_H,fapar_L,NDVI,ALB,Rabsb_sno,Rabsb_bare,Rabsb_urb,Rabsb_wat,Rabsb_rock,Rabsb_ice,Rabsb_deb,...
    soil_alb,e_gr,e_sur]=ShortwaveFluxes(Ccrown,Cbare,Crock,Curb,Cwat,Csno,Cice,...
    Rsw,PAR,SvF,dw_SNO,hc_H,hc_L,SNDtm1,ydepth,ICE_Dtm1,Cdeb,Deb_Par,Urb_Par,h_S,snow_albtm1,Aice,OS,Color_Class,OM_H,OM_L,...
    LAI_H,SAI_H,LAIdead_H,LAI_L,SAI_L,LAIdead_L,PFT_opt_H,PFT_opt_L);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Roughness --
[zom,zoh,disp_h,zom_H,zom_L,zoh_H,zoh_L,disp_h_H,disp_h_L,zom_under]=Roughness_New(SNDtm1,ydepth,ICE_Dtm1,Cdeb,Deb_Par,Urb_Par,hc_H,hc_L,LAI_H,Ccrown_L,Cwat,Curb,Crock,Cice);
%%%%%%%%%%%%%%%%%%%%
%%% Control unit of measurement -- Pre
if Pre < 100 || Pre > 1100
    disp('Error Atm. Pressure Units')
    return
end
%%%%%%%%%%%%%
%%%% Neutrel undercanopy resistence
[rap_H,rap_L,rb_H,rb_L]=Undercanopy_Leaf_Resistence2(Ws,Ta,Ta,Ccrown,hc_H,hc_L,...
    (LAI_H+SAI_H+LAIdead_H),(LAI_L+SAI_L+LAIdead_L),d_leaf_H,d_leaf_L,...
    zatm,disp_h,zom,zom_under,SNDtm1,disp_h_H,zom_H,disp_h_L,zom_L);
%%%%
for i=1:length(Ccrown)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Litter resistance %%%%%%%
    [r_litter(i),alp_litter(i)]=Litter_Resistence(Ws,Ta,Pre,zatm,disp_h,zom,Sllit,BLit(i),In_Littertm1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%% Neutrel undercanopy resistence
    %[rap_H(i),rap_L(i),rb_H(i),rb_L(i)]=Undercanopy_Leaf_Resistence(Ws,Ta,Ta,hc_H(i),hc_L(i),...
    %    (LAI_H(i)+SAI_H(i)+LAIdead_H(i)),(LAI_L(i)+SAI_L(i)+LAIdead_L(i)),d_leaf_H(i),d_leaf_L(i),...
    %    zatm,disp_h,zom,zom_under,SNDtm1,disp_h_H(i),zom_H(i),disp_h_L(i),zom_L(i));
    %%%%%%% Stomatal resistance
    if ((LAI_H(i) > 0) && (OPT_VegSnow==0)) || ((LAI_H(i)>0) && (OPT_VegSnow==1) && (dw_SNO <= 0.5))
        %%%
        ran = (1/((0.4^2)*Ws))*(log((zatm-disp_h)/zom))^2; %%% Neutral aerodynamic resistance  %%[s/m]
        %%%
        [rs_sunH(i),rs_shdH(i),Ci_sunH(i),Ci_shdH(i),An_H(i),Rdark_H(i),Lpho_H(i),SIF_H(i)]=Canopy_Resistence_An_Evolution(PAR_sun_H(i),PAR_shd_H(i),LAI_H(i),...
            Kopt_H(i),KnitH(i),FsunH(i),FshdH(i),Citm1_sunH(i),Citm1_shdH(i),...
            Ca,ran,rb_H(i),Ta,Ta,Pre,Ds,...
            Psi_ltm1_H(i),Psi_sto_50_H(i),Psi_sto_00_H(i),...
            CT_H(i),Vmax_H(i),DSE_H(i),Ha_H(i),FI_H(i),Oa,Do_H(i),a1_H(i),go_H(i),e_rel_H(i),e_relN_H(i),gmes_H(i),rjv_H(i),mSl_H(i),Sl_H(i),Opt_CR);
    else
        rs_sunH(i) = Inf; rs_shdH(i) = Inf; An_H(i) = 0; Rdark_H(i)=0;  Ci_sunH(i)=0;  Ci_shdH(i)=0; Lpho_H(i)=0; SIF_H(i) =0 ;
    end
    %%%%%%%%%%%%%%%%
    if (LAI_L(i) > 0) && (Csno == 0) && (Cice == 0) && (ydepth-0.15 < hc_L(i))
        %%%%%%%%%%%%
        %ran = (1/((0.4^2)*Ws))*(log((zatm-disp_h)/zom))^2 + rap_H(i); %%% Neutral aerodynamic resistance  %%[s/m]
        if hc_H(i) > 0
            %ran = (1/((0.4^2)*Ws_undertm1))*(log((hc_L(i)+2-disp_h_L)/zom_L))^2 ; %%% Neutral aerodynamic resistance  %%[s/m]
            ran = (1/((0.4^2)*Ws))*(log((zatm-disp_h)/zom))^2 + rap_H(i); %%% Neutral aerodynamic resistance  %%[s/m]
        else
            ran = (1/((0.4^2)*Ws))*(log((zatm-disp_h)/zom))^2 ; %%% Neutral aerodynamic resistance  %%[s/m]
        end
        %%%
        [rs_sunL(i),rs_shdL(i),Ci_sunL(i),Ci_shdL(i),An_L(i),Rdark_L(i),Lpho_L(i),SIF_L(i)]=Canopy_Resistence_An_Evolution(PAR_sun_L(i),PAR_shd_L(i),LAI_L(i),...
            Kopt_L(i),KnitL(i),FsunL(i),FshdL(i),Citm1_sunL(i),Citm1_shdL(i),...
            Ca,ran,rb_L(i),Ta,Ta,Pre,Ds,...
            Psi_ltm1_L(i),Psi_sto_50_L(i),Psi_sto_00_L(i),...
            CT_L(i),Vmax_L(i),DSE_L(i),Ha_L(i),FI_L(i),Oa,Do_L(i),a1_L(i),go_L(i),e_rel_L(i),e_relN_L(i),gmes_L(i),rjv_L(i),mSl_L(i),Sl_L(i),Opt_CR);
    else
        rs_sunL(i) = Inf; rs_shdL(i) = Inf; An_L(i) = 0; Rdark_L(i)=0; Ci_sunL(i)=0; Ci_shdL(i)=0;  Lpho_L(i)=0; SIF_L(i)=0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%% Energy spent in photosynthesis
Lpho = (Lpho_H + Lpho_L)*Ccrown'; %% [W/m2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Water available for Transpiration and Evaporation in a given time step
row=1000;  dth=1;
Vavail=(Vtm1/dth)*(row/3600/1000); %%% Volume available in the soil for evaporation [kg/m^2.s]
Vavail_mul=ones(length(Ccrown),1)*Vavail.*((Ccrown'/sum(Ccrown))*ones(1,length(Vtm1))); % Volume available in the different soil depth for plant [kg/m^2.s]
Vavail= sum(Vavail);
%%%% %%% Vavail,Vavail_plant_H,Vavail_plant_L [kg/m^2.s]
%%%% Exwat_H,Exwat_L [mm/h]
Vavail_plant_H =  sum(min(Vavail_mul.*RfH_Zs,Exwat_H*(row/3600/1000)),2)  ;  %%  [kg/m^2.s]
Vavail_plant_L =  sum(min(Vavail_mul.*RfL_Zs,Exwat_L*(row/3600/1000)),2) ;  %%  [kg/m^2.s]
%Vavail_plant_H = ((Vltm1_H + Vxtm1_H).*Ccrown) ; % [mm]
%Vavail_plant_L = ((Vltm1_L + Vxtm1_L).*Ccrown) ; % [mm]
%Vavail_plant_H =  (Vavail_plant_H'/dth)*(row/3600/1000) + sum(Vavail_mul.*RfH_Zs,2);  %%  [kg/m^2.s]
%Vavail_plant_L =  (Vavail_plant_L'/dth)*(row/3600/1000) + sum(Vavail_mul.*RfL_Zs,2);  %%  [kg/m^2.s]
Vavail_plant_H= Vavail_plant_H';
Vavail_plant_L= Vavail_plant_L';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%% ENERGY SOLUTION
if Csno == 1
    %lan_sno = lan_air + (7.75*1e-5*rostm1 + 1.105*1e-6*rostm1.^2).*(lan_ice-lan_air); %%  [W/m K ] Thermal conductivity snow
    lan_sno =  0.023 + (7.75*1e-5*rostm1 + 1.105*1e-6*rostm1.^2).*(2.29-0.023); %%  [W/m K ] Thermal conductivity snow
    %%%%%%%
    if Cice == 1
        if Cdeb == 1
            if SNDtm1<=min_SPD
                % specific heat capacity ice = 2093 J / kg K
                Gsno = lan_sno*(Tstm1-Tdebtm1(1))/(SNDtm1+0.001); %%% [W m-2]
                Gsno_max =2093*(SWEtm1+In_SWEtm1)*(Tstm1-Tdebtm1(1))/dt; %%  Maximum Flux [W /m^2 ]
                Gsno= sign(Gsno)*min([abs(Gsno),abs(Gsno_max)]);
            else
                ms_sno=length(Tdpsnowtm1);
                Gsno = lan_sno*(Tdpsnowtm1(ms_sno)-Tdebtm1(1))/(SNDtm1+0.001); %%% [W m-2]
                Gsno_max =2093*(SWEtm1+In_SWEtm1)*(Tdpsnowtm1(ms_sno)-Tdebtm1(1))/dt; %%  Maximum Flux [W /m^2 ]
                Gsno= sign(Gsno)*min([abs(Gsno),abs(Gsno_max)]);
            end

            %%% Snow over debris covered ice
            %%% Debris Heat Flux
            ms_deb = length(Tdebtm1);
            lan_deb =  Deb_Par.lan*ones(1,ms_deb) ;%%% [W/m K ] Thermal conductivity debris
            cv_deb =  Deb_Par.cs*Deb_Par.rho*ones(1,ms_deb); % [J/m^3 K]  Volumetric heat capcity debris

            [G,Tdeb,Gn]=Soil_Heat_Profile_Normal(NaN,dt,Tdebtm1,ms_deb,Zs_deb,lan_deb,cv_deb,Ticetm1,Gsno,NaN,4);
            %%%%
            Gice_max =2093*(ICEtm1)*(Tdeb(ms_deb)-Ticetm1(1))/dt; %%  Maximum Flux [W /m^2 ]
            Gn= sign(Gn)*min([abs(Gn),abs(Gice_max)]);

        else
            if SNDtm1<=min_SPD
                Gsno = lan_sno*(Tstm1-Ticetm1(1))/(SNDtm1+0.001); %%% [W m-2]
                Gsno_max =2093*(SWEtm1+In_SWEtm1)*(Tstm1-Ticetm1(1))/dt; %%  Maximum Flux [W /m^2 ]
                Gice_max =2093*(ICEtm1)*(Tstm1-Ticetm1(1))/dt; %%  Maximum Flux [W /m^2 ]
                Gsno= sign(Gsno)*min([abs(Gsno),abs(Gsno_max),abs(Gice_max)]);
            else
                ms_sno=length(Tdpsnowtm1);
                Gsno = lan_sno*(Tdpsnowtm1(ms_sno)-Ticetm1(1))/(SNDtm1+0.001); %%% [W m-2]
                Gsno_max =2093*(SWEtm1+In_SWEtm1)*(Tdpsnowtm1(ms_sno)-Ticetm1(1))/dt; %%  Maximum Flux [W /m^2 ]
                Gice_max =2093*(ICEtm1)*(Tdpsnowtm1(ms_sno)-Ticetm1(1))/dt; %%  Maximum Flux [W /m^2 ]
                Gsno= sign(Gsno)*min([abs(Gsno),abs(Gsno_max),abs(Gice_max)]);
            end
            %%% Snow over ice
            G = Gsno ;
        end
        %Tdamp = Tdamptm1;
        %Tdamp = Tdp(1);
    else
        if SNDtm1<=min_SPD
            Gsno =  lan_sno*(Tstm1-Tdptm1(1))/(SNDtm1+0.001); %%% [W m-2]
            Gsno_max =2093*(SWEtm1+In_SWEtm1)*(Tstm1-Tdptm1(1))/dt; %%  Maximum Flux [W /m^2 ]
            Gsno= sign(Gsno)*min([abs(Gsno),abs(Gsno_max)]);
        else
            ms_sno=length(Tdpsnowtm1);
            Gsno =  lan_sno*(Tdpsnowtm1(ms_sno)-Tdptm1(1))/(SNDtm1+0.001); %%% [W m-2]
            Gsno_max =2093*(SWEtm1+In_SWEtm1)*(Tdpsnowtm1(ms_sno)-Tdptm1(1))/dt; %%  Maximum Flux [W /m^2 ]
            Gsno= sign(Gsno)*min([abs(Gsno),abs(Gsno_max)]);
        end

        %%%%% Snow over other surfaces
        if Crock ==1 || Curb ==1  || Cwat ==1
            if OPT_SoilTemp==1
                [G,Tdp]=Soil_Heat_Profile_Normal(NaN,dt,Tdptm1,ms,Zs,lan_oth,cv_oth,NaN,Gsno,0,1);
                O=Otm1 ; Oice=Oicetm1;
                Tdamp = Tdp(1);
            else
                [G,Tdamp]=Soil_Heat(dt,0,Tstm1,Tdamptm1,CTt_oth);
                O=Otm1 ; Oice=Oicetm1;
                Tdp=Tdamp*ones(1,ms);
            end

        else
            %%%%% Snow cover over soil
            if OPT_SoilTemp==1
                [G,Tdp,O,Oice]=Soil_Heat_Profile_New(NaN,dt,Tdptm1,ms,dz,Zs,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
                    Phy,s_SVG,bVG,Osat,Ohy,Oicetm1,Otm1,Gsno,OPT_FR_SOIL,OPT_STh);
                Tdamp = Tdp(1);
            else
                %%% Snow over soil
                %%% Soil Heat Flux in Snow
                [G,Tdamp]=Soil_Heat(dt,0,Tstm1,Tdamptm1,CTt);
                Tdp=Tdamp*ones(1,ms);
                [~,~,~,Oice,O]=Soil_Thermal_properties_FT(Tdp,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
                    Phy,s_SVG,bVG,Osat,Ohy,(Otm1+Oicetm1),OPT_FR_SOIL);
            end
        end
    end
    Gfin=G;
    %%%% underneath surface cannot cool the snow
    %if G>0
    %    G=0;
    %end
    %%%%%
    if (dw_SNO <= 0.5) && (OPT_VegSnow==1) && (sum(hc_H)>0)
        %%%%%%%%%%% Uncovered Vegetation Surface Temperature
        [TsV]=fzero(@Surface_Temperature_VegSnow,Tstm0,Opt_ST,Tstm1,dt,Ta,ea,Latm,SvF,Pre,...
            Csno,Ccrown,...
            hc_H,hc_L,SNDtm1,LAI_H,LAI_L,(SAI_H+LAIdead_H),(SAI_L+LAIdead_L),...
            RabsbSun_vegH,RabsbShd_vegH,FsunH,FshdH,...
            FsunL,FshdL,e_snotm1,e_gr,...
            dw_H,dw_SNO,In_Htm1,...
            rs_sunH,rs_shdH,d_leaf_H,d_leaf_L,...
            zatm,disp_h,zom,zoh,zom_under,disp_h_H,zom_H,disp_h_L,zom_L,Ws,Lpho,Vavail_plant_H);
    else
        TsV = 0;
    end
    %%%%%
    for rr=1:2
        [Ts,~,exitflag]=fzero(@Surface_Temperature_Snow,Tstm0,Opt_ST,dt,Ta,ea,Latm,SvF,Pre,...
            Csno,Crock,Curb,Cbare,Ccrown,Cwat,Cice,Cfol_H,...
            hc_H,hc_L,SNDtm1,ydepth,ICE_Dtm1,Cdeb,LAI_H,LAI_L,(SAI_H+LAIdead_H),(SAI_L+LAIdead_L),...
            RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,...
            RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,FsunH,FshdH,...
            FsunL,FshdL,Rabsb_sno,Rabsb_bare,Rabsb_urb,Rabsb_wat,Rabsb_rock,Rabsb_ice,Rabsb_deb,...
            e_snotm1,e_gr,e_sur,Cicew,Csnow,CLitter,...
            dw_L,dw_H,dw_SNO,In_max_SWE,...
            In_Htm1,In_Ltm1,In_urbtm1,In_rocktm1,SWEtm1,In_SWEtm1,...
            Pr_liq,Pr_sno,rs_sunH,rs_sunL,rs_shdH,rs_shdL,d_leaf_H,d_leaf_L,r_litter,r_soil,b_soil,alp_soil,...
            Tstm1,G,Tdpsnowtm1,lan_sno,...
            zatm,disp_h,zom,zoh,zom_under,disp_h_H,zom_H,disp_h_L,zom_L,Ws,In_Littertm1,alp_litter,Pr_sno_day,Th_Pr_sno,ros_max1,ros_max2,...
            Tdew,t_slstm1,SWEtm1,SNDtm1,rostm1,SP_wctm1,In_SWEtm1,fpr,Vavail,Vavail_plant_H,Vavail_plant_L,WATtm1,ICEtm1,OPT_VegSnow,min_SPD,TsV);
        if exitflag>0
            break
        else
            if rr==1
                Tstm0=Ta-0.1;
            end
        end
    end
    %%%% %%%%%%%%%%%%%%%%%%%
    T2_flag=0; Ts_under=Tstm1_under; %% Case without 2 temperatures
else
    %%%% Ice without the snow
    if Cice > 0
        if Cdeb == 1
            %%% debris covered glacier
            ms_deb = length(Tdebtm1);
            [Ts]=fzero(@Surface_Temperature_Debris,Tstm0,Opt_ST,dt,Ta,ea,Latm,SvF,Pre,...
                Csno,Crock,Curb,Cbare,Ccrown,Cwat,Cice,...
                hc_H,hc_L,SNDtm1,ydepth,ICE_Dtm1,Cdeb,Zs_deb,LAI_H,LAI_L,(SAI_H+LAIdead_H),(SAI_L+LAIdead_L),...
                RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,...
                RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,FsunH,FshdH,...
                FsunL,FshdL,Rabsb_sno,Rabsb_bare,Rabsb_urb,Rabsb_wat,Rabsb_rock,Rabsb_ice,Rabsb_deb,...
                e_snotm1,e_gr,e_sur,Cicew,Csnow,CLitter,...
                dw_L,dw_H,dw_SNO,...
                In_Htm1,In_Ltm1,In_urbtm1,In_rocktm1,SWEtm1,In_SWEtm1,...
                Pr_liq,Pr_sno,rs_sunH,rs_sunL,rs_shdH,rs_shdL,d_leaf_H,d_leaf_L,r_litter,r_soil,b_soil,alp_soil,...
                Tdebtm1,Ticetm1,Deb_Par,ms_deb,...
                zatm,disp_h,zom,zoh,zom_under,disp_h_H,zom_H,disp_h_L,zom_L,Ws,In_Littertm1,alp_litter,...
                Vavail,Vavail_plant_H,Vavail_plant_L,WATtm1,ICEtm1);
        else
            Gice =  2.29*(Tstm1-Tdptm1(1))/(ICE_Dtm1+0.001); %%% [W m-2]
            Gice_max =2093*(ICEtm1)*(Tstm1-Tdptm1(1))/dt; %%  Maximum Flux [W /m^2 ]
            Gice= sign(Gice)*min([abs(Gice),abs(Gice_max)]);
            %%%%% Ice over other surfaces
            if Crock ==1 || Curb ==1  || Cwat ==1
                if OPT_SoilTemp==1
                    [G,Tdp]=Soil_Heat_Profile_Normal(NaN,dt,Tdptm1,ms,Zs,lan_oth,cv_oth,NaN,Gice,0,1);
                    O=Otm1 ; Oice=Oicetm1;
                    Tdamp = Tdp(1);
                else
                    [G,Tdamp]=Soil_Heat(dt,0,Tstm1,Tdamptm1,CTt_oth);
                    Tdp=Tdamp*ones(1,ms);  O=Otm1 ; Oice=Oicetm1;
                end
            else
                %%% Ice over  soil
                if OPT_SoilTemp==1
                    [G,Tdp,O,Oice]=Soil_Heat_Profile_New(NaN,dt,Tdptm1,ms,dz,Zs,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
                        Phy,s_SVG,bVG,Osat,Ohy,Oicetm1,Otm1,Gice,OPT_FR_SOIL,OPT_STh);
                    Tdamp = Tdp(1);
                else
                    %%% Ice over soil
                    %%% Soil Heat Flux in Snow
                    [G,Tdamp]=Soil_Heat(dt,0,Tstm1,Tdamptm1,CTt);
                    Tdp=Tdamp*ones(1,ms);
                    [~,~,~,Oice,O]=Soil_Thermal_properties_FT(Tdp,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
                        Phy,s_SVG,bVG,Osat,Ohy,(Otm1+Oicetm1),OPT_FR_SOIL);
                end

            end
            Gfin=G;
            %%%% underneath surface cannot cool the ice
            if G>0
                G=0;
            end
            [Ts]=fzero(@Surface_Temperature_Ice,Tstm0,Opt_ST,dt,Ta,ea,Latm,SvF,Pre,...
                Csno,Crock,Curb,Cbare,Ccrown,Cfol_H,Cwat,Cice,...
                hc_H,hc_L,SNDtm1,ydepth,ICE_Dtm1,Cdeb,LAI_H,LAI_L,(SAI_H+LAIdead_H),(SAI_L+LAIdead_L),...
                RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,...
                RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,FsunH,FshdH,...
                FsunL,FshdL,Rabsb_sno,Rabsb_bare,Rabsb_urb,Rabsb_wat,Rabsb_rock,Rabsb_ice,Rabsb_deb,...
                e_snotm1,e_gr,e_sur,Cicew,Csnow,CLitter,...
                dw_L,dw_H,dw_SNO,...
                In_Htm1,In_Ltm1,In_urbtm1,In_rocktm1,SWEtm1,In_SWEtm1,...
                Pr_liq,Pr_sno,rs_sunH,rs_sunL,rs_shdH,rs_shdL,d_leaf_H,d_leaf_L,r_litter,r_soil,b_soil,alp_soil,...
                Tstm1,G,...
                zatm,disp_h,zom,zoh,zom_under,disp_h_H,zom_H,disp_h_L,zom_L,Ws,In_Littertm1,alp_litter,Vavail,Vavail_plant_H,Vavail_plant_L,WATtm1,...
                ICEtm1,IP_wctm1,Ice_wc_sp);
            %%%%%%%%%%%%%%%%%%%
        end
        %%%% %%%%%%%%%%%%%%%%%%%
        T2_flag=0; Ts_under=Tstm1_under; %% Case without 2 temperatures
    else
        if  (Cbare == 1) || (Cwat == 1) || (Crock==1) || (Curb == 1) || isnan(Tstm1_under)
            [Ts]=fzero(@Surface_Temperature,Tstm0,Opt_ST,dt,Ta,ea,Latm,SvF,Pre,...
                Csno,Crock,Curb,Cbare,Ccrown,Cwat,Cice,...
                hc_H,hc_L,SNDtm1,ydepth,ICE_Dtm1,Cdeb,LAI_H,LAI_L,(SAI_H+LAIdead_H),(SAI_L+LAIdead_L),...
                RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,...
                RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,FsunH,FshdH,...
                FsunL,FshdL,Rabsb_sno,Rabsb_bare,Rabsb_urb,Rabsb_wat,Rabsb_rock,Rabsb_ice,Rabsb_deb,...
                e_snotm1,e_gr,e_sur,Cicew,Csnow,CLitter,...
                dw_L,dw_H,dw_SNO,...
                In_Htm1,In_Ltm1,In_urbtm1,In_rocktm1,SWEtm1,In_SWEtm1,...
                Pr_liq,Pr_sno,rs_sunH,rs_sunL,rs_shdH,rs_shdL,d_leaf_H,d_leaf_L,r_litter,r_soil,b_soil,alp_soil,...
                ms,dz,Zs,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,Phy,s_SVG,bVG,Osat,Ohy,...
                Oicetm1,Otm1,...
                Tstm1,Tdamptm1,Tdptm1,CTt,CTt_oth,OPT_FR_SOIL,OPT_STh,...
                zatm,disp_h,zom,zoh,zom_under,disp_h_H,zom_H,disp_h_L,zom_L,Ws,In_Littertm1,alp_litter,...
                Lpho,Vavail,Vavail_plant_H,Vavail_plant_L,WATtm1,ICEtm1,OPT_SoilTemp);
            %%%% %%%%%%%%%%%%%%%%%%%
            T2_flag=0; Ts_under=Tstm1_under; %% Case without 2 temperatures
        else
            Tinit=[Ta+0.1,Tstm1_under];
            for rr=1:3
                %[Ts2] = lsqnonlin(@Surface_Temperature_2Temp,[Tstm1,Tstm1_under],[],[],Opt_ST2,dt,Ta,ea,Latm,SvF,Pre,...
                [Ts2,~,exitflag] = fsolve(@Surface_Temperature_2Temp,Tinit,Opt_ST2,dt,Ta,ea,Latm,SvF,Pre,...
                    Csno,Crock,Curb,Cbare,Ccrown,Cwat,Cice,...
                    hc_H,hc_L,SNDtm1,ydepth,ICE_Dtm1,Cdeb,LAI_H,LAI_L,(SAI_H+LAIdead_H),(SAI_L+LAIdead_L),...
                    RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,...
                    RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,FsunH,FshdH,...
                    FsunL,FshdL,Rabsb_sno,Rabsb_bare,Rabsb_urb,Rabsb_wat,Rabsb_rock,Rabsb_ice,Rabsb_deb,...
                    e_snotm1,e_gr,e_sur,Cicew,Csnow,CLitter,...
                    dw_L,dw_H,dw_SNO,...
                    In_Htm1,In_Ltm1,In_urbtm1,In_rocktm1,SWEtm1,In_SWEtm1,...
                    Pr_liq,Pr_sno,rs_sunH,rs_sunL,rs_shdH,rs_shdL,d_leaf_H,d_leaf_L,r_litter,r_soil,b_soil,alp_soil,...
                    ms,dz,Zs,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,Phy,s_SVG,bVG,Osat,Ohy, ...
                    Oicetm1,Otm1,...
                    Tstm1,Tstm1_under,Tdamptm1,Tdptm1,CTt,CTt_oth,OPT_FR_SOIL,OPT_STh,...
                    zatm,disp_h,zom,zoh,zom_under,disp_h_H,zom_H,disp_h_L,zom_L,Ws,In_Littertm1,alp_litter,...
                    Lpho,Vavail,Vavail_plant_H,Vavail_plant_L,WATtm1,ICEtm1,OPT_SoilTemp);
                %%%%%%%%
                if exitflag>0
                    break
                else
                    if rr==1
                        Tinit=[Ta-0.1,Tstm1_under];
                    end
                    if rr==2
                        Tinit=[Tstm1,Tstm1_under];
                    end
                end
            end
            %%%%
            T2_flag=1;
            Ts=Ts2(1); Ts_under=Ts2(2);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(isreal(Ts)) ||  isnan(Ts)
    disp('Error Computation Ts')
    return
end
if abs(Ts) >= 120
    disp('Numerical instability on Ts')
    return
end
if abs(Ts) <= eps
    %%%% To prevent numerical instabilities
    Ts=0.0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% POST COMPUTATION %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Aerodynamic resistence
[ra]=Aerodynamic_Resistence(Ta,Ts,Pre,zatm,disp_h,zom,zoh,Ws,ea);
%%%%%%%%%% Undercanopy and Leaf resitence
[rap_H,rap_L,rb_H,rb_L,Ws_under]=Undercanopy_Leaf_Resistence2(Ws,Ta,Ts,Ccrown,hc_H,hc_L,...
    (LAI_H+SAI_H+LAIdead_H),(LAI_L+SAI_L+LAIdead_L),d_leaf_H,d_leaf_L,...
    zatm,disp_h,zom,zom_under,SNDtm1,disp_h_H,zom_H,disp_h_L,zom_L);
%for i=1:length(Ccrown)
%    %%%%%%
%    if  (hc_L(i) == 0) && (hc_H(i) == 0)
%        rap_H(i) = 0; rap_L(i) = 0 ;
%        rb_H(i)=0; rb_L(i)=0;
%        ums=Ws; ums2=0;
%    else
%        [rap_H(i),rap_L(i),rb_H(i),rb_L(i),ums,ums2]=Undercanopy_Leaf_Resistence(Ws,Ta,Ts,hc_H(i),hc_L(i),...
%            (LAI_H(i)+SAI_H(i)+LAIdead_H(i)),(LAI_L(i)+SAI_L(i)+LAIdead_L(i)),d_leaf_H(i),d_leaf_L(i),...
%            zatm,disp_h,zom,zom_under,SNDtm1,disp_h_H(i),zom_H(i),disp_h_L(i),zom_L(i));
%    end
%    if (hc_L(i)>0) && (hc_H(i)>0)  %%% Two vegetations
%        Ws_under(i) = ums2;
%    else
%        Ws_under(i)= ums;
%    end
%end
Ws_under = max(Ws_under);
%%%%
%%% Net Radiation
if (dw_SNO <= 0.5) && (OPT_VegSnow==1) && (sum(hc_H)>0) &&  (Csno == 1)
    [Rn]=Net_Radiation_Manager_SnowVeg(Ts,TsV,Latm,SvF,...
        Csno,Crock,Curb,Cwat,Cbare,Cice,Ccrown,...
        hc_L,SNDtm1,ydepth,ICE_Dtm1,LAI_H,LAI_L,SAI_H,SAI_L,...
        RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,...
        RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,FsunH,FshdH,...
        FsunL,FshdL,Rabsb_sno,Rabsb_bare,Rabsb_urb,Rabsb_wat,Rabsb_rock,Rabsb_ice,...
        e_snotm1,e_gr,e_sur,Cicew,Csnow);
else
    if T2_flag==1
        [Rn,Rn_under]=Net_Radiation_Manager_2Temp(Ts,Ts_under,Latm,SvF,...
            Csno,Crock,Curb,Cwat,Cbare,Cice,Ccrown,...
            hc_L,SNDtm1,ydepth,ICE_Dtm1,Cdeb,LAI_H,LAI_L,(SAI_H+LAIdead_H),(SAI_L+LAIdead_L),...
            RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,...
            RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,FsunH,FshdH,...
            FsunL,FshdL,Rabsb_sno,Rabsb_bare,Rabsb_urb,Rabsb_wat,Rabsb_rock,Rabsb_ice,Rabsb_deb,...
            e_snotm1,e_gr,e_sur,Cicew,Csnow);
        Rn=Rn+Rn_under; 
    else
        %%% Net Radiation  Normal conditions
        [Rn]=Net_Radiation_Manager(Ts,Latm,SvF,...
            Csno,Crock,Curb,Cwat,Cbare,Cice,Ccrown,...
            hc_L,SNDtm1,ydepth,ICE_Dtm1,Cdeb,LAI_H,LAI_L,(SAI_H+LAIdead_H),(SAI_L+LAIdead_L),...
            RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,...
            RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,FsunH,FshdH,...
            FsunL,FshdL,Rabsb_sno,Rabsb_bare,Rabsb_urb,Rabsb_wat,Rabsb_rock,Rabsb_ice,Rabsb_deb,...
            e_snotm1,e_gr,e_sur,Cicew,Csnow);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if T2_flag==1
    [H,H_under,QE,QE_under,Qv,Qv_under,T_H,T_L,EIn_H,EIn_L,EG,ESN,ESN_In,EWAT,ELitter,EICE,EIn_urb,EIn_rock]=Heat_fluxes_2Temp(dt,...
        Ta,Ts,Ts_under,ea,Pre,Csno,Crock,Curb,Cwat,Cbare,Cice,Cicew,Csnow,CLitter,Cdeb,...
        dw_L,dw_H,dw_SNO,Ccrown,FsunH,FshdH,...
        FsunL,FshdL,LAI_H,LAI_L,(SAI_H+LAIdead_H),(SAI_L+LAIdead_L),...
        In_Htm1,In_Ltm1,In_urbtm1,In_rocktm1,SWEtm1,In_SWEtm1,...
        Pr_liq,Pr_sno,ra,rs_sunH,rs_sunL,rs_shdH,rs_shdL,rb_H,rb_L,rap_H,rap_L,r_litter,...
        r_soil,b_soil,alp_soil,Vavail,Vavail_plant_H,Vavail_plant_L,WATtm1,ICEtm1,In_Littertm1,alp_litter);
    H=H+H_under; 
    QE=QE+QE_under; 
    Qv=Qv+Qv_under; 
else
    [H,QE,Qv,T_H,T_L,EIn_H,EIn_L,EG,ESN,ESN_In,EWAT,ELitter,EICE,EIn_urb,EIn_rock]=Heat_fluxes(dt,...
        Ta,Ts,ea,Pre,Csno,Crock,Curb,Cwat,Cbare,Cice,Cicew,Csnow,CLitter,Cdeb,...
        dw_L,dw_H,dw_SNO,Ccrown,FsunH,FshdH,...
        FsunL,FshdL,LAI_H,LAI_L,(SAI_H+LAIdead_H),(SAI_L+LAIdead_L),...
        In_Htm1,In_Ltm1,In_urbtm1,In_rocktm1,SWEtm1,In_SWEtm1,...
        Pr_liq,Pr_sno,ra,rs_sunH,rs_sunL,rs_shdH,rs_shdL,rb_H,rb_L,rap_H,rap_L,r_litter,...
        r_soil,b_soil,alp_soil,Vavail,Vavail_plant_H,Vavail_plant_L,WATtm1,ICEtm1,In_Littertm1,alp_litter);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if not(Csno == 1 || Cice == 1)
    %%%%% Thermal regime other surfaces
    if Crock ==1 || Curb ==1  || Cwat ==1
        if OPT_SoilTemp==1
            [G,Tdamp]=Soil_Heat(dt,Ts,Tstm1,Tdamptm1,CTt_oth);
            [Gfin,Tdp]=Soil_Heat_Profile_Normal(Ts,dt,Tdptm1,ms,Zs,lan_oth,cv_oth,NaN,NaN,0,3);
            O=Otm1 ; Oice=Oicetm1;
        else
            [G,Tdamp]=Soil_Heat(dt,Ts,Tstm1,Tdamptm1,CTt_oth);
            Gfin=G; Tdp=Tdamp*ones(1,ms);  O=Otm1 ; Oice=Oicetm1;
        end
    else
        if T2_flag==1
            %%%%% Thermal regime soil
            if OPT_SoilTemp==1
                [G,Tdamp]=Soil_Heat(dt,Ts_under,Tstm1_under,Tdamptm1,CTt);
                [Gfin,Tdp,O,Oice]=Soil_Heat_Profile_New(Ts_under,dt,Tdptm1,ms,dz,Zs,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
                    Phy,s_SVG,bVG,Osat,Ohy,Oicetm1,Otm1,NaN,OPT_FR_SOIL,OPT_STh);
            else
                [G,Tdamp]=Soil_Heat(dt,Ts_under,Tstm1_under,Tdamptm1,CTt);
                Gfin=G; Tdp=Tdamp*ones(1,ms);
                [~,~,~,Oice,O]=Soil_Thermal_properties_FT(Tdp,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
                    Phy,s_SVG,bVG,Osat,Ohy,(Otm1+Oicetm1),OPT_FR_SOIL);
            end
        else
            %%%%% Thermal regime soil
            if OPT_SoilTemp==1
                [G,Tdamp]=Soil_Heat(dt,Ts,Tstm1,Tdamptm1,CTt);
                [Gfin,Tdp,O,Oice]=Soil_Heat_Profile_New(Ts,dt,Tdptm1,ms,dz,Zs,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
                    Phy,s_SVG,bVG,Osat,Ohy,Oicetm1,Otm1,NaN,OPT_FR_SOIL,OPT_STh);
            else
                [G,Tdamp]=Soil_Heat(dt,Ts,Tstm1,Tdamptm1,CTt);
                Gfin=G; Tdp=Tdamp*ones(1,ms);
                [~,~,~,Oice,O]=Soil_Thermal_properties_FT(Tdp,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
                    Phy,s_SVG,bVG,Osat,Ohy,(Otm1+Oicetm1),OPT_FR_SOIL);
            end
        end
    end
end
%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Csno == 1
    %%%%%
    if SNDtm1 > min_SPD

        %%%% Snow cover with temperature profile and two layers

        [TsF,Tdpsnow,SWE,SND,ros,In_SWE,SP_wc,WR_SP,U_SWE,NIn_SWE,dQ,Qfm,t_sls,Smelt,Gres,dTres]=Snowpacks_2layers(dt,...
            Ta,Ts,Tstm1,Tdpsnowtm1,Tdew,Ws,t_slstm1,SWEtm1,SNDtm1,rostm1,SP_wctm1,In_SWEtm1,In_max_SWE,dw_SNO,...
            Pr_liq,Pr_sno,ESN,ESN_In,Rn,H,QE,G,Qv,Csnow,Ccrown,Cwat,Cfol_H,fpr,Pr_sno_day,Th_Pr_sno,ros_max1,ros_max2,lan_sno,min_SPD); 

    else
        %%%% Snow cover not deep enough to compute the temperature profile.
        %%%%
        [TsF,SWE,SND,ros,In_SWE,SP_wc,WR_SP,U_SWE,NIn_SWE,dQ,Qfm,t_sls,Smelt,Gres]=Snowpacks(dt,...
            Ta,Ts,Tstm1,Tdew,Ws,t_slstm1,SWEtm1,SNDtm1,rostm1,SP_wctm1,In_SWEtm1,In_max_SWE,dw_SNO,...
            Pr_liq,Pr_sno,ESN,ESN_In,Rn,H,QE,G,Qv,Csnow,Ccrown,Cwat,Cfol_H,fpr,Pr_sno_day,Th_Pr_sno,ros_max1,ros_max2);
        Tdpsnow=0*Tdpsnowtm1;
    end

    %%%%%%
    %%% Gres - additional heat left after the complete snow melting - currently lost
    %%%%%
    if (dw_SNO <= 0.5) && (OPT_VegSnow==1) && (sum(hc_H)>0)
        [raV]=Aerodynamic_Resistence(Ta,TsV,Pre,zatm,disp_h,zom,zoh,Ws,ea);
        %%%%
        [RnV]=Net_Radiation_Manager_VegSnow(TsV,Tstm1,Latm,SvF,...
            Csno,Ccrown,...
            hc_L,SNDtm1,LAI_H,LAI_L,(SAI_H+LAIdead_H),(SAI_L+LAIdead_L),...
            RabsbSun_vegH,RabsbShd_vegH,FsunH,FshdH,FsunL,FshdL,...
            e_snotm1,e_gr);
        %%%%%%%%%%%%%%
        [HV,QEV,T_HV,EIn_HV]=Heat_fluxes_VegSnow(dt,...
            Ta,TsV,ea,Pre,Csno,...
            dw_H,dw_SNO,Ccrown,FsunH,FshdH,...
            LAI_H,(SAI_H+LAIdead_H),In_Htm1,...
            raV,rs_sunH,rs_shdH,rb_H,Vavail_plant_H);
        %%%%%%%%%%%%%%%%%
        %T_HdisV=  T_HV'*ones(1,length(Vavail)).*RfH_Zs; %%% % [mm/h]
        %T_HdisV =min((T_HV'*ones(1,length(Vavail)).*RfH_Zs),Vavail_mul); %%%
        %%%%%%%%%%
        T_H=T_HV; EIn_H=EIn_HV; %T_Hdis=T_HdisV;
        dQVEG = RnV-HV-QEV-Lpho;
    else
        dQVEG=0; HV=0; QEV=0;
    end
    %%% Ice below snow
    if Cice == 1
        %%%
        if Cdeb == 1

            %%% Presence of debris cover == to rock for hydrology
            [In_deb,SE_deb]=Interceptions_Debris(dt,Csno,Cdeb,...
                In_urbtm1,In_max_urb,Pr_liq,WR_SP,EIn_urb);


            %%%% Icepacks
            [Tice,ICE,ICE_D,IP_wc,WR_IP,dQI,QfmI,Imelt]=Icepack(dt,NaN,Ticetm1,ICEtm1,IP_wctm1,...
                0,0,0,0,0,-Gn,0,Cwat,Ccrown,Cfol_H,Csno,Cicew,Ice_wc_sp,SE_deb);
            SE_deb = 0;

        else
            %%%% Icepacks
            [Tice,ICE,ICE_D,IP_wc,WR_IP,dQI,QfmI,Imelt]=Icepack(dt,NaN,Ticetm1,ICEtm1,IP_wctm1,...
                0,0,0,0,0,-G,0,Cwat,Ccrown,Cfol_H,Csno,Cicew,Ice_wc_sp,WR_SP);
            Tdeb=0;
        end

        if Crock ==1 || Curb ==1  || Cwat ==1
            [~,Tdp]=Soil_Heat_Profile_Normal(Tice,dt,Tdptm1,ms,Zs,lan_oth,cv_oth,NaN,NaN,0,3);
            O=Otm1 ; Oice=Oicetm1;
        else
            [~,Tdp,O,Oice]=Soil_Heat_Profile_New(Tice,dt,Tdptm1,ms,dz,Zs,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
                Phy,s_SVG,bVG,Osat,Ohy,Oicetm1,Otm1,NaN,OPT_FR_SOIL,OPT_STh);
        end
        Tdamp = Tdamptm1;
    else
        %%%%%%%% --> Parameter Icepack also when there isn't
        Tice = 0;
        ICE = 0;
        ICE_D = 0;
        IP_wc = 0 ;
        WR_IP = IP_wctm1 + ICEtm1 ;  %% All the terms should be equal to zero
        Imelt = 0;
        Tdeb=0;
    end
else
    %%% Ice without the snow
    if Cice == 1
        if Cdeb == 1
            %%% Debris Heat Flux
            lan_deb =  Deb_Par.lan*ones(1,ms_deb) ;%%% [W/m K ] Thermal conductivity debris
            cv_deb =  Deb_Par.cs*Deb_Par.rho*ones(1,ms_deb); % [J/m^3 K]  Volumetric heat capcity debris
            [G,Tdeb,Gn]=Soil_Heat_Profile_Normal(Ts,dt,Tdebtm1,ms_deb,Zs_deb,lan_deb,cv_deb,Ticetm1,NaN,NaN,2);
            %%% Heat Flux into the debris G
            %%% Heat Flux from the debris  Gn

            Gice_max =2093*(ICEtm1)*(Tdeb(ms_deb)-Ticetm1(1))/dt + 333700*(ICEtm1)/dt*(Tdeb(ms_deb)>Ticetm1(1)) ; %%  Maximum Flux [W /m^2 ]
            Gn= sign(Gn)*min([abs(Gn),abs(Gice_max)]);

            TsF = Tdeb(1) + G*(0.001*Zs_deb(2)*0.5)/lan_deb(1);%% [C]

            %%% Presence of debris cover == to rock for hydrology
            [In_deb,SE_deb]=Interceptions_Debris(dt,Csno,Cdeb,...
                In_urbtm1,In_max_urb,Pr_liq,0,EIn_urb);

            %%%% Icepacks
            [Tice,ICE,ICE_D,IP_wc,WR_IP,dQ,Qfm,Imelt]=Icepack(dt,NaN,Ticetm1,ICEtm1,IP_wctm1,...
                0,0,0,0,0,-Gn,0,Cwat,Ccrown,Cfol_H,Csno,Cicew,Ice_wc_sp,SE_deb);
            SE_deb = 0;

            %%%%%%%%%  Flux from ice to below
            if Crock ==1 || Curb ==1  || Cwat ==1
                [~,Tdp]=Soil_Heat_Profile_Normal(Tice,dt,Tdptm1,ms,Zs,lan_oth,cv_oth,NaN,NaN,0,3);
                O=Otm1 ; Oice=Oicetm1;
            else
                [~,Tdp,O,Oice]=Soil_Heat_Profile_New(Tice,dt,Tdptm1,ms,dz,Zs,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
                    Phy,s_SVG,bVG,Osat,Ohy,Oicetm1,Otm1,NaN,OPT_FR_SOIL,OPT_STh);
            end
            Tdamp = Tdamptm1;
            Gfin=G;
        else
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [TsF,ICE,ICE_D,IP_wc,WR_IP,dQ,Qfm,Imelt]=Icepack(dt,Ts,Tstm1,ICEtm1,IP_wctm1,...
                Pr_liq,EICE,Rn,H,QE,G,Qv,Cwat,Ccrown,Cfol_H,Csno,Cicew,Ice_wc_sp,0);
            Tice = TsF;
            Tdeb = 0;
        end
        %%%%%
        %%%%
        if Csno == 0 %%% --> Variables Snowpack also when there isn't
            SWE=0;
            WR_SP= Pr_sno*(1-Cwat)+ dw_SNO*sum(Ccrown)*Pr_liq + SP_wctm1 + SWEtm1 + In_SWEtm1 ; %% All the terms should be equal to zero
            SND=0;
            ros=0;
            In_SWE=0;
            SP_wc=0;
            U_SWE=0;
            t_sls =0;
            NIn_SWE=0;
            dQVEG=0; HV=0; QEV = 0;
            TsV=0;
            Smelt = 0;
            Tdpsnow=Tdpsnowtm1*0; 
        end
        %%%%%
    else
        %%%%%%%%%%
        %%% --> Variables Snowpack Icepack Debris and snow free vegetation also when there aren't
        %%%  TsF,SWE,D,ros,In_SWE,SP_wc,WR_SP,U_SWE,dQ,Qf,t_sls
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        TsF=0;
        SWE=0;
        WR_SP= Pr_sno*(1-Cwat)+ dw_SNO*sum(Ccrown)*Pr_liq + SP_wctm1 + SWEtm1 + In_SWEtm1 ; %% All the terms should be equal to zero
        SND=0;
        ros=0;
        In_SWE=0;
        SP_wc=0;
        U_SWE=0;
        dQ=0;
        Qfm=0;
        t_sls =0;
        NIn_SWE=0;
        dQVEG=0;  HV=0; QEV = 0;
        TsV=0;
        ICE = 0;
        ICE_D = 0;
        IP_wc = 0 ;
        WR_IP = IP_wctm1 + ICEtm1 ; %% All the terms should be equal to zero
        Imelt = 0;
        Smelt = 0;
        Tdpsnow=Tdpsnowtm1*0; 
        Tice = 0;
        Tdeb = 0;
    end
end
%%%%%
%%%%%%%%%%%%%%%%%%
if Csno == 1 || Cicew == 1
    [ICE,ICE_D,SWE,SND,SP_wc,WAT,NIce]=Snow2Ice(dt,Ts,ICE,SWE,SND,SP_wc,WATtm1,Cicew,ros,ros_Ice_thr,Cwat,Cicewtm1,dz_ice,WatFreez_Th);
else
    WAT=WATtm1;
    NIce= 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Treating Water Logging su Asur
q_run = 1000*ydepth/dth; %%%[mm/h]
if ydepth > 0
    if ELitter >= q_run
        ELitter = ELitter-q_run; %%[mm/h]
        EWAT = EWAT + q_run;
        WAT = WAT + q_run; 
        q_run = 0;
    else
        q_run = q_run - ELitter; %%[mm/h]
        EWAT = EWAT + ELitter;
        WAT = WAT + ELitter; 
        ELitter = 0;
    end
    if EG >= q_run
        EG = EG-q_run; %%[mm/h]
        EWAT = EWAT + q_run;
        WAT = WAT + q_run; 
        q_run = 0;
    else
        q_run = q_run - EG; %%[mm/h]
        EWAT = EWAT + EG;
        WAT = WAT + EG; 
        EG = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ENERGY BALANCE CONTROLS
if Csno > 0  || Cice >0
    %%%% Closure Energy for Surface Temperature (Snowpack)
    DT = TsF -Ts;  %% ---> 0
    if (Cdeb == 1) && (Csno == 0)
        DQ = Rn-H-QE-G+Qv +(-dQ+Gn);
    else
        DQ = Rn-H-QE-G+Qv+Qfm-dQ; %%%% ---> 0
    end
else
    %%%%% Closure Energy for Surface Temperature (Snow Free)
    DQ = Rn-H-QE-G+Qv-Lpho; %% ---> 0
    DT = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%
if (OPT_PlantHydr == 1)
    Vavail_plant_H =  sum(Vavail_mul.*RfH_Zs,2)*dth/(row/3600/1000) ;  %%  [mm]
    Vavail_plant_L =  sum(Vavail_mul.*RfL_Zs,2)*dth/(row/3600/1000);  %%  [mm]
    PLD_H=zeros(1,length(Ccrown)); PLD_L=zeros(1,length(Ccrown));
else
    %Ll= 1000*(2501.3 - 2.361*(Ta)); %%% Latent heat vaporization/condensaition [J/kg]
    %Tinit = (sum(T_H)+sum(T_L))/(1000*3600/1000); %%  % [kg/m^2.s]
    %T_Hdis= min((T_H'*ones(1,length(Vavail)).*RfH_Zs),(Vtm1/dth)); %%% % [mm/h]
    %T_Ldis= min((T_L'*ones(1,length(Vavail)).*RfL_Zs),(Vtm1/dth)); %%% %[mm/h]
    %T_H=sum(T_Hdis,2);  %[mm/h]
    %T_L=sum(T_Ldis,2);  %[mm/h]
    %Tfin = (sum(T_H)+sum(T_L))/(1000*3600/1000); %
    %DQ= DQ + (Tfin-Tinit)*Ll;
end
rho2 = 55555555; %% [mmolH20 /m^3]; Water density
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%% PLANT HYDRAULIC %%%%%%%%%%%%
for i=1:length(Ccrown)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (ZR95_H(i) ==0) || (OPT_PlantHydr == 0)
        if  (ZR95_H(i)> 0)
            Psi_l_H(i)=Psi_s_H(i);
            Psi_x_H(i)=Psi_s_H(i);
            Jsx_H(i)=T_H(i);
        else
            Psi_l_H(i)=0;
            Psi_x_H(i)=0;
            Jsx_H(i)=0;
            gsr_H(i)=0;
        end
        Jxl_H(i)=0;
        Kleaf_H(i)=0;
        Kx_H(i)=0;
        Vx_H(i)= Vxtm1_H(i);
        Vl_H(i)= Vltm1_H(i);
    else
        %%%%%
        [PLD_H(i)]=Plant_Disconnection(Psi_xtm1_H(i),Psi_ltm1_H(i),Axyl_H(i),PsiL50_H(i),PsiL00_H(i),PsiX50_H(i));
        %%%
        X0 =[Vxtm1_H(i)*rho2/1000;Vltm1_H(i)*rho2/1000];
        if PLD_H(i)== 0
            [Tout,Xout]=ode23s(@PLANT_HYDRAULIC_DIFF_NEW,[0 dt],X0,OPT_PH,...
                Ccrown(i),Psi_s_H(i),T_H(i),hc_H(i),LAI_H(i),Axyl_H(i),PsiL50_H(i),PsiL00_H(i),Kleaf_max_H(i),Cl_H(i),Sl_H(i),Kx_max_H(i),PsiX50_H(i),Cx_H(i),gsr_H(i),...
                dt,Vavail_plant_H(i));
        else  %% Plant wilted -- No status change --
            Xout=X0';
        end
        %%%%%
        [Psi_x_H(i),Psi_l_H(i),Jsx_H(i),Jxl_H(i),Kleaf_H(i),Kx_H(i),Vx_H(i),Vl_H(i)]=Plant_Hydraulic(Xout,X0,dth,Ccrown(i),T_H(i),Psi_s_H(i),hc_H(i),LAI_H(i),...
            Axyl_H(i),PsiL50_H(i),PsiL00_H(i),Kleaf_max_H(i),Kx_max_H(i),PsiX50_H(i),Sl_H(i),mSl_H(i),Cx_H(i),Cl_H(i));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (ZR95_L(i)==0) || (OPT_PlantHydr == 0)
        if (ZR95_L(i)>0)
            Psi_x_L(i)=Psi_s_L(i);
            Psi_l_L(i)=Psi_s_L(i);
            Jsx_L(i)=T_L(i);
        else
            Psi_x_L(i)=0;
            Psi_l_L(i)=0;
            Jsx_L(i)=0;
            gsr_L(i)=0;
        end
        Jxl_L(i)=0;
        Kleaf_L(i)=0;
        Kx_L(i)=0;
        Vx_L(i)= Vxtm1_L(i);
        Vl_L(i)= Vltm1_L(i);
    else
        %%%%%
        [PLD_L(i)]=Plant_Disconnection(Psi_xtm1_L(i),Psi_ltm1_L(i),Axyl_L(i),PsiL50_L(i),PsiL00_L(i),PsiX50_L(i));
        %%%
        X0 =[Vxtm1_L(i)*rho2/1000;Vltm1_L(i)*rho2/1000];
        if PLD_L(i)== 0
            [Tout,Xout]=ode23s(@PLANT_HYDRAULIC_DIFF_NEW,[0 dt],X0,OPT_PH,...
                Ccrown(i),Psi_s_L(i),T_L(i),hc_L(i),LAI_L(i),Axyl_L(i),PsiL50_L(i),PsiL00_L(i),Kleaf_max_L(i),Cl_L(i),Sl_L(i),Kx_max_L(i),PsiX50_L(i),Cx_L(i),gsr_L(i),...
                dt,Vavail_plant_L(i));
        else   %% Plant wilted -- No status change --
            Xout=X0';
        end
        %%%%%
        [Psi_x_L(i),Psi_l_L(i),Jsx_L(i),Jxl_L(i),Kleaf_L(i),Kx_L(i),Vx_L(i),Vl_L(i)]=Plant_Hydraulic(Xout,X0,dth,Ccrown(i),T_L(i),Psi_s_L(i),hc_L(i),LAI_L(i),...
            Axyl_L(i),PsiL50_L(i),PsiL00_L(i),Kleaf_max_L(i),Kx_max_L(i),PsiX50_L(i),Sl_L(i),mSl_L(i),Cx_L(i),Cl_L(i));
    end
end
%%%%%%%% Distributed Sink %%%%%%%%%%%%%%%%%%%%%
EG_dis = EG*EvL_Zs; %% % [mm/h]
%T_Hdis=  T_H'*ones(1,length(Vavail)).*RfH_Zs; %%% % [mm/h]
%T_Ldis=  T_L'*ones(1,length(Vavail)).*RfL_Zs; %%% % [mm/h]
%EG_dis = min(EG*EvL_Zs,Vavail); %% %
%T_Hdis= min((T_H'*ones(1,length(Vavail)).*RfH_Zs),Vavail_mul); %%% % % [kg/m^2.s]
%T_Ldis= min((T_L'*ones(1,length(Vavail)).*RfL_Zs),Vavail_mul); %%% % % [kg/m^2.s]
%%%%%%%%%%
J_Hdis=  Jsx_H'*ones(1,length(Vavail_mul)).*RfH_Zs; %%% % [mm/h]
J_Ldis=  Jsx_L'*ones(1,length(Vavail_mul)).*RfL_Zs; %%% % [mm/h]
%%%%%% For numerical reasons only
for rj=1:length(Ccrown)
    if  sum(J_Hdis(rj,:))>0
        J_Hdis(rj,:) = Jsx_H(rj)./sum(J_Hdis(rj,:)).*J_Hdis(rj,:);% [mm/h]
    end
    if  sum(J_Ldis(rj,:))>0
        J_Ldis(rj,:) = Jsx_L(rj)./sum(J_Ldis(rj,:)).*J_Ldis(rj,:);% [mm/h]
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% INTERCEPTION and WATER TO THE SOIL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[In_H,In_L,In_Litter,Dr_H,Dr_L,WIS]=Interceptions_Veg(dt,...
    Ccrown,Cfol_H,Cfol_L,CLitter,Cbare,Csno,Cice,Crock,Curb,Cwat,dw_SNO,In_Ltm1,In_Htm1,In_Littertm1,...
    In_max_H,In_max_L,BLit,...
    Pr_liq,EIn_H,EIn_L,ELitter,q_run,WR_SP,WR_IP,gcI,KcI);
%%%%%%%%%%%%%%
if Cdeb == 1
    [~,In_rock,SE_rock]=Interceptions_Other(dt,...
        Csno,Crock,Curb,Cice,...
        In_urbtm1,In_rocktm1,In_max_urb,In_max_rock,...
        Pr_liq,WR_SP,WR_IP,q_run,EIn_urb,EIn_rock);
    In_urb = In_deb; SE_urb =SE_deb;
else
    [In_urb,In_rock,SE_rock,SE_urb]=Interceptions_Other(dt,...
        Csno,Crock,Curb,Cice,...
        In_urbtm1,In_rocktm1,In_max_urb,In_max_rock,...
        Pr_liq,WR_SP,WR_IP,q_run,EIn_urb,EIn_rock);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Updated Snow Albedo
if Csno > 0
    [snow_alb,tau_sno,e_sno]=Albedo_Snow_Properties(dt,SWE,h_S,Ts,Ta,SWEtm1,tau_snotm1,snow_albtm1,Th_Pr_sno,Pr_sno_day,Aice,Deb_Par,Cdeb,Cice,Ta_day,Pr_sno,Pr_liq,ros,N);
else
    snow_alb.dir_vis = ALB;
    snow_alb.dir_nir = ALB;
    snow_alb.dif_vis = ALB;
    snow_alb.dif_nir = ALB;
    tau_sno = 0; e_sno = e_gr;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return