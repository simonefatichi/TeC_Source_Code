%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% MAIN_FRAME OF HBM-VEG %%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INITIZIALIZATION VARIABLES
%%%%%--->>  INIZIALIZATION VARIABLES  <<--- %%%%%%%%%%%%%%%%%%%%%%%
%%% j time dt = 1 day  %%%%%%%%%%%%%%%%
%%% i time dt = 1h
%%% ms soil layer
%%% cc Crown Area number
%%% NN time step
dtd = 1; %% [day]
dth = dt/3600; %% [hour]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
V=zeros(NN,ms);
O=zeros(NN,ms);
Qi_out=zeros(NN,ms);
WTR=zeros(NN,ms);
POT=zeros(NN,ms);
Tdp=zeros(NN,ms);
%Sdp=zeros(NN,ms);
Vice=zeros(NN,ms);
Oice=zeros(NN,ms);
%%%%%%%%%%%%%%%%%%%%%%%%%%%
OF=zeros(NN,1); OS=zeros(NN,1);
ZWT=zeros(NN,1);
Pr_sno=zeros(NN,1); Pr_liq=zeros(NN,1);
er=zeros(NN,1);
Ts=zeros(NN,1);ra=zeros(NN,1);r_soil=zeros(NN,1);
b_soil=zeros(NN,1); Tice=zeros(NN,1);
Rn=zeros(NN,1);H=zeros(NN,1);
QE=zeros(NN,1);Qv=zeros(NN,1); Lpho=zeros(NN,1);
EG=zeros(NN,1); ESN=zeros(NN,1);ESN_In=zeros(NN,1);
EWAT=zeros(NN,1);EIn_urb=zeros(NN,1);
EIn_rock=zeros(NN,1); dw_SNO=zeros(NN,1);
G=zeros(NN,1); SWE=zeros(NN,1);
SND=zeros(NN,1);%snow_alb=zeros(NN,1);
ros=zeros(NN,1);In_SWE=zeros(NN,1);SP_wc=zeros(NN,1);
WR_SP=zeros(NN,1);U_SWE=zeros(NN,1);NIn_SWE=zeros(NN,1);
dQ=zeros(NN,1);Qfm=zeros(NN,1);t_sls=zeros(NN,1);
DQ=zeros(NN,1);DT=zeros(NN,1);
In_urb=zeros(NN,1);In_rock=zeros(NN,1);
SE_rock=zeros(NN,1);SE_urb=zeros(NN,1);
Csno=zeros(NN,1);WIS=zeros(NN,1);
Lk=zeros(NN,1);f=zeros(NN,1);
Rh=zeros(NN,1);Rd=zeros(NN,1);
NDVI=zeros(NN,1);alp_soil=zeros(NN,1);
tau_sno=zeros(NN,1);e_sno=zeros(NN,1);ALB=zeros(NN,1);
EK=zeros(NN,1);
dQVEG=zeros(NN,1);TsV=zeros(NN,1);
HV=zeros(NN,1);QEV=zeros(NN,1);
Cice=zeros(NN,1);
Lk_wat=zeros(NN,1); Lk_rock=zeros(NN,1);
EICE=zeros(NN,1);
WAT=zeros(NN,1);ICE=zeros(NN,1);ICE_D=zeros(NN,1);
WR_IP=zeros(NN,1);NIce=zeros(NN,1);Cicew=zeros(NN,1);
IP_wc=zeros(NN,1);
Csnow=zeros(NN,1);FROCK=zeros(NN,1);
Imelt=zeros(NN,1);Smelt=zeros(NN,1);
Tdamp=zeros(NN,1); Gfin=zeros(NN,1);
In_Litter=zeros(NN,1); ELitter=zeros(NN,1);
Ws_under=zeros(NN,1);
Ts_under=NaN*ones(NN,1);
Tdpsnow = zeros(NN,5);
if  not(exist('IrD','var'))
    IrD=zeros(NN,1);
end
if  not(exist('Salt','var'))
    Salt=zeros(NN,1); %%% Salt = Salt Concentration [mol Salt/ m-3] 
end
NetWatWet=zeros(NN,1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
r_litter=zeros(NN,cc);
%%%
In_H=zeros(NN,cc);In_L=zeros(NN,cc);
OH=zeros(NN,cc);OL=zeros(NN,cc);
Psi_s_H=zeros(NN,cc); Psi_s_L=zeros(NN,cc);
Tdp_H=zeros(NN,cc);Tdp_L=zeros(NN,cc);
rap_H=zeros(NN,cc);rap_L=zeros(NN,cc);
Dr_H=zeros(NN,cc);Dr_L=zeros(NN,cc);
T_H=zeros(NN,cc);T_L=zeros(NN,cc);EIn_H=zeros(NN,cc);EIn_L=zeros(NN,cc);
rb_H=zeros(NN,cc);rb_L=zeros(NN,cc);
An_H=zeros(NN,cc);Rdark_H=zeros(NN,cc);
An_L=zeros(NN,cc);Rdark_L=zeros(NN,cc);
Ci_sunH=zeros(NN,cc);Ci_sunL=zeros(NN,cc);
Ci_shdH=zeros(NN,cc);Ci_shdL=zeros(NN,cc);
rs_sunH=zeros(NN,cc);rs_sunL=zeros(NN,cc);
rs_shdH=zeros(NN,cc);rs_shdL=zeros(NN,cc);
%%%%%
Vx_H=zeros(NN,cc);Vl_H=zeros(NN,cc);
Vx_L=zeros(NN,cc);Vl_L=zeros(NN,cc);
Psi_x_H=zeros(NN,cc);Psi_l_H=zeros(NN,cc);
Psi_x_L=zeros(NN,cc);Psi_l_L=zeros(NN,cc);
gsr_H=zeros(NN,cc);
Jsx_H=zeros(NN,cc); Jxl_H=zeros(NN,cc);
Kleaf_H=zeros(NN,cc);Kx_H=zeros(NN,cc);
gsr_L=zeros(NN,cc);
Jsx_L=zeros(NN,cc);Jxl_L=zeros(NN,cc);
Kleaf_L=zeros(NN,cc);Kx_L=zeros(NN,cc);
fapar_H=zeros(NN,cc);fapar_L=zeros(NN,cc);
SIF_H=zeros(NN,cc);SIF_L=zeros(NN,cc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NCP = 8; %% Number of Carbon Pool
NNd = ceil(NN/24)+1;
LAI_L=zeros(NNd,cc); B_L=zeros(NNd,cc,NCP); NPP_L=zeros(NNd,cc);Rg_L=zeros(NNd,cc);
RA_L=zeros(NNd,cc);Rmc_L=zeros(NNd,cc);ANPP_L=zeros(NNd,cc);
Rms_L=zeros(NNd,cc);Rmr_L=zeros(NNd,cc); PHE_S_L=zeros(NNd,cc); dflo_L=zeros(NNd,cc);
AgeL_L=zeros(NNd,cc);e_rel_L=ones(NNd,cc);SAI_L=zeros(NNd,cc); hc_L=zeros(NNd,cc); e_relN_L=ones(NNd,cc); BA_L=zeros(NNd,cc);
%%%
LAI_H=zeros(NNd,cc); B_H=zeros(NNd,cc,NCP);NPP_H=zeros(NNd,cc);Rg_H=zeros(NNd,cc);
RA_H=zeros(NNd,cc);Rms_H=zeros(NNd,cc);Rmr_H=zeros(NNd,cc); ANPP_H=zeros(NNd,cc);
PHE_S_H=zeros(NNd,cc); dflo_H=zeros(NNd,cc);Rmc_H=zeros(NNd,cc);
AgeL_H=zeros(NNd,cc);e_rel_H=ones(NNd,cc);SAI_H=zeros(NNd,cc); hc_H=zeros(NNd,cc); e_relN_H=ones(NNd,cc); BA_H=zeros(NNd,cc);
%%%
Rrootl_H=zeros(NNd,cc); Rrootl_L=zeros(NNd,cc);
Bfac_dayH=ones(NNd,cc); Bfac_weekH=ones(NNd,cc); NPPI_H=zeros(NNd,cc); TdpI_H=zeros(NNd,cc);
Bfac_dayL=ones(NNd,cc); Bfac_weekL=ones(NNd,cc); NPPI_L=zeros(NNd,cc); TdpI_L=zeros(NNd,cc);
%%%
Sr_H=zeros(NNd,cc);  Slf_H=zeros(NNd,cc);
Sfr_H=zeros(NNd,cc); Swm_H=zeros(NNd,cc);  Sll_H=zeros(NNd,cc);
Sr_L=zeros(NNd,cc); Slf_L=zeros(NNd,cc);
Sfr_L=zeros(NNd,cc); Swm_L=zeros(NNd,cc); Sll_L=zeros(NNd,cc);
LAIdead_H=zeros(NNd,cc); LAIdead_L=zeros(NNd,cc);
Rexmy_H=zeros(NNd,cc,3); AgeDL_H=zeros(NNd,cc); Nreserve_H=zeros(NNd,cc);
Preserve_H=zeros(NNd,cc); Kreserve_H=zeros(NNd,cc);
rNc_H=zeros(NNd,cc);rPc_H=zeros(NNd,cc);rKc_H=zeros(NNd,cc);
Rexmy_L=zeros(NNd,cc,3); AgeDL_L=zeros(NNd,cc); Nreserve_L=zeros(NNd,cc);
Preserve_L=zeros(NNd,cc); Kreserve_L=zeros(NNd,cc);
rNc_L=zeros(NNd,cc);rPc_L=zeros(NNd,cc);rKc_L=zeros(NNd,cc);
NBLeaf_H =zeros(NNd,cc); PARI_H=zeros(NNd,cc,3) ; NBLI_H=zeros(NNd,cc);
NBLeaf_L =zeros(NNd,cc);PARI_L=zeros(NNd,cc,3) ; NBLI_L=zeros(NNd,cc);
%%%%%
NupI_H=zeros(NNd,cc,3);
NupI_L=zeros(NNd,cc,3);
NavlI=zeros(NNd,3);
%%%
RB_L=zeros(NNd,cc,7);
RB_H=zeros(NNd,cc,7);
NuLit_H =zeros(NNd,cc,3);
NuLit_L =zeros(NNd,cc,3);
%%%%%
BLit=zeros(NNd,cc);
%%%%%
Nuptake_H=zeros(NNd,cc);
Puptake_H=zeros(NNd,cc);
Kuptake_H=zeros(NNd,cc);
FNC_H=zeros(NNd,cc);
TexC_H=zeros(NNd,cc);TexN_H=zeros(NNd,cc);TexP_H=zeros(NNd,cc);TexK_H=zeros(NNd,cc);
TNIT_H=zeros(NNd,cc);TPHO_H=zeros(NNd,cc);TPOT_H=zeros(NNd,cc);
SupN_H=zeros(NNd,cc);SupP_H=zeros(NNd,cc);SupK_H=zeros(NNd,cc);
%%%%
Nuptake_L=zeros(NNd,cc);
Puptake_L=zeros(NNd,cc);
Kuptake_L=zeros(NNd,cc);
FNC_L=zeros(NNd,cc);
TexC_L=zeros(NNd,cc);TexN_L=zeros(NNd,cc);TexP_L=zeros(NNd,cc);TexK_L=zeros(NNd,cc);
TNIT_L=zeros(NNd,cc);TPHO_L=zeros(NNd,cc);TPOT_L=zeros(NNd,cc);
SupN_L=zeros(NNd,cc);SupP_L=zeros(NNd,cc);SupK_L=zeros(NNd,cc);
%%%%
ISOIL_H=zeros(NNd,cc,18);
ISOIL_L=zeros(NNd,cc,18);
ManIH = zeros(cc,1);
ManIL = zeros(cc,1);

AgrHarNut =  zeros(NNd,3);
%%%%%%%%%%%%%%%%%
jDay=zeros(NNd,1);L_day=zeros(NNd,1);
%%%%
Ccrown_t =ones(NNd,cc);
AgePl_H =zeros(NNd,cc); AgePl_L =zeros(NNd,cc);
Tden_H =zeros(NNd,cc); Tden_L =zeros(NNd,cc);
TBio_Ht =zeros(NNd,cc); TBio_Lt =zeros(NNd,cc);
ZR95_Ht =zeros(NNd,cc); ZR95_Lt =zeros(NNd,cc);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if  not(exist('OPT_BG','var'))
    OPT_SoilBiogeochemistry = 0;
else
    if  OPT_BG == 0
        OPT_SoilBiogeochemistry = 0;
    else
        OPT_SoilBiogeochemistry = 1;
    end 
end
%%%%%%%%%%%%%
if OPT_SoilBiogeochemistry == 1
    P=zeros(NNd,55);
    R_litter=zeros(NNd,1);
    R_microbe=zeros(NNd,1);
    R_litter_sur=zeros(NNd,1);
    R_ew=zeros(NNd,1);
    VOL=zeros(NNd,1);
    N2flx = zeros(NNd,1);
    Min_N = zeros(NNd,1);
    Min_P = zeros(NNd,1);
    R_bacteria= zeros(NNd,1);
    RmycAM = zeros(NNd,1);
    RmycEM = zeros(NNd,1);
    Prod_B = zeros(NNd,1);
    Prod_F = zeros(NNd,1);
    BfixN = zeros(NNd,1);
    LitFirEmi =  zeros(NNd,2);
    LEAK_NH4 = zeros(NNd,1);
    LEAK_NO3 = zeros(NNd,1);
    LEAK_P = zeros(NNd,1);
    LEAK_K = zeros(NNd,1);
    LEAK_DOC = zeros(NNd,1);
    LEAK_DON = zeros(NNd,1);
    LEAK_DOP = zeros(NNd,1);
    R_NH4 = zeros(NNd,1);
    R_NO3 = zeros(NNd,1);
    R_P = zeros(NNd,1);
    R_K = zeros(NNd,1);
    R_DOC = zeros(NNd,1);
    R_DON = zeros(NNd,1);
    R_DOP = zeros(NNd,1);
    RexmyI=zeros(NNd,3);
end

%%%%%%%%% CALL PARAMETERS AND INITIAL CONDITION
run(PARAM_IC);
Restating_parameters;
if length(Oice)==1
    Oice=zeros(NN,ms);
end
%%%%%%%%%%
Tdeb =zeros(NN,max(1,length(Zs_deb)-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Check_Land_Cover_Fractions(Crock,Curb,Cwat,Cbare,Ccrown);
CcrownFIX = Ccrown;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Lateral Contribution
q_runon=zeros(NN,1); %%[mm/h]
Qi_in=zeros(NN,ms); %%[mm/h]
%%%%%%%%%%%%%%%%%%%%%%%%%
tic ;
CK1=zeros(NN,1);  CK2=zeros(NN,1);
%profile on
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
j=1; %%% time dt = 1 day  %%%%%%%%%%%%%%%%
Tstm0= Ts(1);
%%% i time dt = 1h
%%% ms soil layer
%%% cc Crown Area present
%%%%%%%%%%%%%%%% NUMERICAL METHODS OPTIONS
%Opt_CR = optimset('TolX',3);
%Opt_ST = optimset('TolX',0.1);
%Opt_ST = optimset('TolX',0.1,'Display','iter');
Opt_CR = optimset('TolFun',1);%,'UseParallel','always');
Opt_ST = optimset('TolFun',0.1);%,'UseParallel','always');
Opt_ST2 = optimset('TolFun',0.1,'Display','off');
OPT_SM=  odeset('AbsTol',0.05,'MaxStep',dth);
OPT_VD=  odeset('AbsTol',0.05);
OPT_PH= odeset('AbsTol',0.01);
OPT_STh = odeset('AbsTol',0.02);
OPT_VegSnow = 1;
OPT_SoilTemp = 1;
OPT_FR_SOIL = 1; 
OPT_min_SPD = Inf; %% [m] minimum snow pack depth to have a multilayer snow 
%%%%
if  not(exist('OPT_VCA','var'))
    OPT_VCA = 0;
    OPT_ALLOME = 0;
end
if  not(exist('OPT_DROOT','var'))
    OPT_DROOT = 0;
end
if  not(exist('OPT_PlantHydr','var'))
    OPT_PlantHydr = 0;
end
if  not(exist('OPT_EnvLimitGrowth','var'))
    OPT_EnvLimitGrowth = 0;
end
if  not(exist('OPT_WET','var'))
    OPT_WET = 0;
else
    if  not(exist('Wlev','var'))
        Wlev=zeros(NN,1);
    end
    Wlevm1 = Rd(1)+Rh(1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=2:NN
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if  (mod(i,1000) == 0) || (i == 2)
        disp('Iter:'); disp(i);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    pdind = [max(1,i-24):i-1]; %% previous day indexes
    %%%
    if (Datam(i,4)==1)
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        j=j+1; [jDay(j)]= JULIAN_DAY(Datam(i,:));
        [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day(j)]= SetSunVariables(Datam(i,:),DeltaGMT,Lon,Lat,t_bef,t_aft);
        clear h_S delta_S zeta_S T_sunrise T_sunset
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%% Biogeochemistry
        [Se_bio,Se_fc,Psi_bio,Tdp_bio,VSUM,VTSUM]=Biogeo_environment(Tdp(pdind,:),O(pdind,:),V(pdind,:),Soil_Param,Phy,SPAR,Bio_Zs);% sum(V(i,:))
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        if OPT_SoilBiogeochemistry == 1
            IS= Ccrown*squeeze(ISOIL_L(j-1,:,:)) + Ccrown*squeeze(ISOIL_H(j-1,:,:));
            REXMY= Ccrown*squeeze(Rexmy_L(j-1,:,:)) + Ccrown*squeeze(Rexmy_H(j-1,:,:));
            FireA = 1*((sum(ManIH==-5) + sum(ManIL==-5)) > 0);
            [P(j,:),LEAK_NH4(j),LEAK_NO3(j),LEAK_P(j),LEAK_K(j),LEAK_DOC(j),LEAK_DON(j),LEAK_DOP(j),...
                R_NH4(j),R_NO3(j),R_P(j),R_K(j),R_DOC(j),R_DON(j),R_DOP(j),...
                Nuptake_H(j,:),Puptake_H(j,:),Kuptake_H(j,:),Nuptake_L(j,:),Puptake_L(j,:),Kuptake_L(j,:),RexmyI(j,:),...
                R_litter(j),R_microbe(j),R_litter_sur(j),R_ew(j),VOL(j),N2flx(j),Min_N(j),Min_P(j),...
                R_bacteria(j),RmycAM(j),RmycEM(j),Prod_B(j),Prod_F(j),BfixN(j),NavlI(j,:),LitFirEmi(j,:)]=BIOGEO_UNIT(P(j-1,:),IS,Zbio,sum(Bio_Zs.*rsd),PHs,Tdp_bio,mean(Ta(pdind)),Psi_bio,Se_bio,Se_fc,VSUM,VTSUM,...
                Ccrown,Bio_Zs,RfH_Zs,RfL_Zs,sum(Lk(pdind)),sum(Rd(pdind)),sum(Rh(pdind)),sum(Pr(pdind)),sum(T_H(pdind,:),1),sum(T_L(pdind,:),1),B_H(j-1,:,3),B_L(j-1,:,3),LAI_H(j-1,:),LAI_L(j-1,:),...
                SupN_H(j-1,:),SupP_H(j-1,:),SupK_H(j-1,:),SupN_L(j-1,:),SupP_L(j-1,:),SupK_L(j-1,:),...
                REXMY,RexmyI(j-1,:),ExEM,NavlI(j-1,:),Pcla,Psan,...
                B_IO,jDay(j),FireA,0);
            %%%%%
            %%%%%
            BLit(j,:)= 0.002*sum(P(j,1:5))*Ccrown; %% %%[kg DM / m2]
            Bam =  P(j,20); %%[gC/m2]
            Bem =  P(j,21); %%[gC/m2]
        else
            %%%%%
            BLit(j,:)= 0.0 ; %  0.100 ; %% %%[kg DM / m2]
            %%%
            Nuptake_H(j,:)= 0.0;
            Puptake_H(j,:)= 0.0;
            Kuptake_H(j,:)= 0.0; %% [gK/m^2 day]
            %%%%
            Nuptake_L(j,:)= 0.0;  %% [gN/m^2 day]
            Puptake_L(j,:)= 0.0;
            Kuptake_L(j,:)= 0.0;
            %%%%
            NavlI(j,:)=[0 0 0];
            %P(j,:)=0.0;
            %LEAK_NH4(j)= 0.0; LEAK_NO3(j)= 0.0; LEAK_P(j)= 0.0; LEAK_K(j) = 0.0; LEAK_DOC(j)= 0.0;
            %R_litter(j)= 0.0; R_microbe(j)= 0.0; R_ew(j)=0; R_litter_sur(j)= 0.0; VOL(j)= 0.0; N2flx(j)= 0.0; BfixN(j)=0; LitFirEmi(j,:)=0;
            %LEAK_DOP(j)= 0.0; LEAK_DON(j)= 0.0; Min_N(j)=0; Min_P(j) =0; R_bacteria(j)=0; RmycAM(j)=0; RmycEM(j)=0;
            %RexmyI(j,:)=0.0;
            Bam=NaN; Bem=NaN;
        end


        %%% projected area n-coordinate
        for cc=1:length(Ccrown)
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%
            if (ZR95_H(cc) > 0) || (ZRmax_H(cc) > 0)
                [LAI_H(j,cc),B_H(j,cc,:),NPP_H(j,cc),ANPP_H(j,cc),Rg_H(j,cc),RA_H(j,cc),Rms_H(j,cc),Rmr_H(j,cc),Rmc_H(j,cc),PHE_S_H(j,cc),...
                    dflo_H(j,cc),AgeL_H(j,cc),e_rel_H(j,cc),e_relN_H(j,cc),LAIdead_H(j,cc),NBLeaf_H(j,cc),Sr_H(j,cc),Slf_H(j,cc),Sfr_H(j,cc),Swm_H(j,cc),Sll_H(j,cc),Rexmy_H(j,cc,:),Rrootl_H(j,cc),...
                    AgeDL_H(j,cc),Bfac_dayH(j,cc),Bfac_weekH(j,cc),NPPI_H(j,cc),TdpI_H(j,cc),NupI_H(j,cc,:),PARI_H(j,cc,:),NBLI_H(j,cc),RB_H(j,cc,:),FNC_H(j,cc),Nreserve_H(j,cc),Preserve_H(j,cc),Kreserve_H(j,cc),...
                    rNc_H(j,cc),rPc_H(j,cc),rKc_H(j,cc),ManIH(cc)]= VEGGIE_UNIT(B_H(j-1,cc,:),PHE_S_H(j-1,cc),dflo_H(j-1,cc),AgeL_H(j-1,cc),AgeDL_H(j-1,cc),...
                    Ta(pdind),Tdp_H(pdind,cc),PARB(pdind)+PARD(pdind),Psi_x_H(pdind,cc),Psi_l_H(pdind,cc),An_H(pdind,cc),Rdark_H(pdind,cc),NPP_H(j-1,cc),jDay(j),Datam(i,:),...
                    NPPI_H(j-1,cc),TdpI_H(j-1,cc),Bfac_weekH(j-1,cc),NupI_H(j-1,cc,:),NavlI(j,:),PARI_H(j-1,cc,:),NBLI_H(j-1,cc),NBLeaf_H(j-1,cc),...
                    L_day(j),Lmax_day,VegH_Param_Dyn,cc,...
                    Nreserve_H(j-1,cc),Preserve_H(j-1,cc),Kreserve_H(j-1,cc),Nuptake_H(j,cc),Puptake_H(j,cc),Kuptake_H(j,cc),FNC_H(j-1,cc),Se_bio,Tdp_bio,...
                    ParEx_H(cc),ExEM,Bam,Bem,Mpar_H(cc),TBio_Ht(j-1,cc),OPT_EnvLimitGrowth,OPT_VCA,OPT_VD,OPT_SoilBiogeochemistry);
                %%%%%%%%%%%%%%%%%%
                [TexC_H(j,cc),TexN_H(j,cc),TexP_H(j,cc),TexK_H(j,cc),TNIT_H(j,cc),TPHO_H(j,cc),TPOT_H(j,cc),NuLit_H(j,cc,:),Nreserve_H(j,cc),Preserve_H(j,cc),Kreserve_H(j,cc),...
                    SupN_H(j,cc),SupP_H(j,cc),SupK_H(j,cc),ISOIL_H(j,cc,:)]= Plant_Exports(B_H(j,cc,:),B_H(j-1,cc,:),NuLit_H(j-1,cc,:),...
                    Slf_H(j,cc),Sfr_H(j,cc),Swm_H(j,cc),Sll_H(j,cc),Sr_H(j,cc),Rexmy_H(j,cc,:),Stoich_H(cc),Mpar_H(cc),fab_H(cc),fbe_H(cc),RB_H(j,cc,:),...
                    Nreserve_H(j,cc),Preserve_H(j,cc),Kreserve_H(j,cc),rNc_H(j,cc),rPc_H(j,cc),rKc_H(j,cc),ManIH(cc),OPT_SoilBiogeochemistry);
                %%%%%%%%%%%%%%%%
                
                %%%%%%%%%%%%%%%%%%%%%% Change in height SAI and Ccrown
                SAI_H(j,cc) =SAI_H(j-1,cc);
                if aSE_H(cc) == 2
                    [hc_H(j,cc)] = GrassHeight(LAI_H(j,cc),LAIdead_H(j,cc));
                elseif aSE_H(cc) == 5
                    [hc_H(j,cc),SAI_H(j,cc),B_H(j,:,:),Ccrown,Nreserve_H(j-1:j,:),Preserve_H(j-1:j,:),Kreserve_H(j-1:j,:),AgrHarNut(j,:)] = CropHeightType(LAI_H(j,cc),LAIdead_H(j,cc),cc,B_H(j,:,:),...
                        Ccrown,Nreserve_H(j,:),Preserve_H(j,:),Kreserve_H(j,:),ManIH,Mpar_H,VegH_Param_Dyn,OPT_SoilBiogeochemistry);
                    %%%%
                else
                    hc_H(j,cc)= hc_H(j-1,cc); %%%[m]
                    if OPT_VCA == 1
                        %%%%%%
                        [Ccrown_t(j,cc),hc_H(j,cc),SAI_H(j,cc),BA_H(j,cc),Tden_H(j,cc),AgePl_H(j,cc),TBio_Ht(j,cc)]=Vegetation_Structural_Attributes(dtd,...
                            Ccrown_t(j-1,cc),B_H(j,cc,:),fab_H(cc),Tden_H(j-1,cc),AgePl_H(j-1,cc),OPT_ALLOME);
                        %%%%
                        Ccrown(cc) = CcrownFIX(cc)*Ccrown_t(j,cc);
                        B_H(j,cc,:)= B_H(j,cc,:)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
                        Nreserve_H(j,cc)=Nreserve_H(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
                        Preserve_H(j,cc)=Preserve_H(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
                        Kreserve_H(j,cc)=Kreserve_H(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
                    end
                end
                if OPT_DROOT == 1 
                    %%%%%%%%%
                    [ZR95_H(cc),RfH_Zs(cc,:)]=Root_Depth_Dynamics(CASE_ROOT,B_H(j,cc,:),B_H(j-1,cc,:),Rrootl_H(j,cc),Zs,ZR95_H(cc),ZR50_H(cc),ZRmax_H(cc),...
                        Bfac_dayH(j,cc),Psan,Tdp_H(pdind,cc),O(pdind,:),Soil_Param,a_root_H(cc));
                end
                ZR95_Ht(j,cc)=ZR95_H(cc); 
                
            else
                LAI_H(j,cc)=0;B_H(j,cc,:)=0;NPP_H(j,cc)=0;ANPP_H(j,cc)=0;Rg_H(j,cc)=0;RA_H(j,cc)=0;Rms_H(j,cc)=0;
                Rmr_H(j,cc)=0;PHE_S_H(j,cc)=1;dflo_H(j,cc)=0;AgeL_H(j,cc)=0;e_rel_H(j,cc)=0; e_relN_H(j,cc)=0;
                SAI_H(j,cc)=0; hc_H(j,cc)=0;
                LAIdead_H(j,cc) =LAIdead_H(j-1,cc);
                Sr_H(j,cc)=0;Slf_H(j,cc)=0;Sfr_H(j,cc)=0;Swm_H(j,cc)=0; Rrootl_H(j,cc)=0;
                Rexmy_H(j,cc,:)=0;AgeDL_H(j,cc)=0;FNC_H(j,cc)=1; Nreserve_H(j,cc)=0;Preserve_H(j,cc)=0;Kreserve_H(j,cc)=0;
                rNc_H(j,cc)=1;rPc_H(j,cc)=1;rKc_H(j,cc)=1;
                Bfac_dayH(j,cc)=0; Bfac_weekH(j,cc)=0; NPPI_H(j,cc)=0; TdpI_H(j,cc)=0; NupI_H(j,cc,:)=0; RB_H(j,cc,:)=0;
                TexC_H(j,cc)=0;TexN_H(j,cc)=0;TexP_H(j,cc)=0;TexK_H(j,cc)=0;TNIT_H(j,cc)=0;TPHO_H(j,cc)=0;TPOT_H(j,cc)=0;
                ISOIL_H(j,cc,:)=0; NuLit_H(j,cc,:)=0;
            
            end
            %%%%%%%%%%%%%%%%%
            if (ZR95_L(cc) > 0) || (ZRmax_L(cc) > 0)
                [LAI_L(j,cc),B_L(j,cc,:),NPP_L(j,cc),ANPP_L(j,cc),Rg_L(j,cc),RA_L(j,cc),Rms_L(j,cc),Rmr_L(j,cc),Rmc_L(j,cc),PHE_S_L(j,cc),...
                    dflo_L(j,cc),AgeL_L(j,cc),e_rel_L(j,cc),e_relN_L(j,cc),LAIdead_L(j,cc),NBLeaf_L(j,cc),Sr_L(j,cc),Slf_L(j,cc),Sfr_L(j,cc),Swm_L(j,cc),Sll_L(j,cc),Rexmy_L(j,cc,:),Rrootl_L(j,cc),...
                    AgeDL_L(j,cc),Bfac_dayL(j,cc),Bfac_weekL(j,cc),NPPI_L(j,cc),TdpI_L(j,cc),NupI_L(j,cc,:),PARI_L(j,cc,:),NBLI_L(j,cc),RB_L(j,cc,:),FNC_L(j,cc),Nreserve_L(j,cc),Preserve_L(j,cc),Kreserve_L(j,cc),...
                    rNc_L(j,cc),rPc_L(j,cc),rKc_L(j,cc),ManIL(cc)]= VEGGIE_UNIT(B_L(j-1,cc,:),PHE_S_L(j-1,cc),dflo_L(j-1,cc),AgeL_L(j-1,cc),AgeDL_L(j-1,cc),...
                    Ta(pdind),Tdp_L(pdind,cc),PARB(pdind)+PARD(pdind),Psi_x_L(pdind,cc),Psi_l_L(pdind,cc),An_L(pdind,cc),Rdark_L(pdind,cc),NPP_L(j-1,cc),jDay(j),Datam(i,:),...
                    NPPI_L(j-1,cc),TdpI_L(j-1,cc),Bfac_weekL(j-1,cc),NupI_L(j-1,cc,:),NavlI(j,:),PARI_L(j-1,cc,:),NBLI_L(j-1,cc),NBLeaf_L(j-1,cc),...
                    L_day(j),Lmax_day,VegL_Param_Dyn,cc,...
                    Nreserve_L(j-1,cc),Preserve_L(j-1,cc),Kreserve_L(j-1,cc),Nuptake_L(j,cc),Puptake_L(j,cc),Kuptake_L(j,cc),FNC_L(j-1,cc),Se_bio,Tdp_bio,...
                    ParEx_L(cc),ExEM,Bam,Bem,Mpar_L(cc),TBio_Lt(j-1,cc),OPT_EnvLimitGrowth,OPT_VCA,OPT_VD,OPT_SoilBiogeochemistry);
                %%%%%%%%%%%%%%%%%%%
                [TexC_L(j,cc),TexN_L(j,cc),TexP_L(j,cc),TexK_L(j,cc),TNIT_L(j,cc),TPHO_L(j,cc),TPOT_L(j,cc),NuLit_L(j,cc,:),Nreserve_L(j,cc),Preserve_L(j,cc),Kreserve_L(j,cc),...
                    SupN_L(j,cc),SupP_L(j,cc),SupK_L(j,cc),ISOIL_L(j,cc,:)]= Plant_Exports(B_L(j,cc,:),B_L(j-1,cc,:),NuLit_L(j-1,cc,:),...
                    Slf_L(j,cc),Sfr_L(j,cc),Swm_L(j,cc),Sll_L(j,cc),Sr_L(j,cc),Rexmy_L(j,cc,:),Stoich_L(cc),Mpar_L(cc),fab_L(cc),fbe_L(cc),RB_L(j,cc,:),...
                    Nreserve_L(j,cc),Preserve_L(j,cc),Kreserve_L(j,cc),rNc_L(j,cc),rPc_L(j,cc),rKc_L(j,cc),ManIL(cc),OPT_SoilBiogeochemistry);
                
                %%%%%%%%%%%%%%%%%%%%%% Change in height SAI and Ccrown
                SAI_L(j,cc) =SAI_L(j-1,cc);
                if aSE_L(cc) == 2
                    [hc_L(j,cc)] = GrassHeight(LAI_L(j,cc),LAIdead_L(j,cc));
                elseif aSE_L(cc) == 5
                    [hc_L(j,cc),SAI_L(j,cc),B_L(j,:,:),Ccrown,Nreserve_L(j,:),Preserve_L(j,:),Kreserve_L(j,:),AgrHarNut(j,:)] = CropHeightType(LAI_L(j,cc),LAIdead_L(j,cc),cc,B_L(j,:,:),...
                        Ccrown,Nreserve_L(j,:),Preserve_L(j,:),Kreserve_L(j,:),ManIL,Mpar_L,VegL_Param_Dyn,OPT_SoilBiogeochemistry);
                else
                    hc_L(j,cc)= hc_L(j-1,cc); %%%[m]
                    if OPT_VCA == 1
                        %%%%%%
                        [Ccrown_t(j,cc),hc_L(j,cc),SAI_L(j,cc),BA_L(j,cc),Tden_L(j,cc),AgePl_L(j,cc),TBio_Lt(j,cc)]=Vegetation_Structural_Attributes(dtd,...
                            Ccrown_t(j-1,cc),B_L(j,cc,:),fab_L(cc),Tden_L(j-1,cc),AgePl_L(j-1,cc),OPT_ALLOME);
                        %%%%
                        Ccrown(cc) = CcrownFIX(cc)*Ccrown_t(j,cc);
                        B_L(j,cc,:)= B_L(j,cc,:)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
                        Nreserve_L(j,cc)=Nreserve_L(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
                        Preserve_L(j,cc)=Preserve_L(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
                        Kreserve_L(j,cc)=Kreserve_L(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
                    end
                end

                if OPT_DROOT == 1 
                    %%%%%%%%%
                    [ZR95_L(cc),RfL_Zs(cc,:)]=Root_Depth_Dynamics(CASE_ROOT,B_L(j,cc,:),B_L(j-1,cc,:),Rrootl_L(j,cc),Zs,ZR95_L(cc),ZR50_L(cc),ZRmax_L(cc),...
                        Bfac_dayL(j,cc),Psan,Tdp_L(pdind,cc),O(pdind,:),Soil_Param,a_root_L(cc));
                end
                ZR95_Lt(j,cc)=ZR95_L(cc); 

            else
                LAI_L(j,cc)=0;B_L(j,cc,:)=0;NPP_L(j,cc)=0;ANPP_L(j,cc)=0;Rg_L(j,cc)=0;RA_L(j,cc)=0;Rms_L(j,cc)=0;
                Rmr_L(j,cc)=0;PHE_S_L(j,cc)=1;dflo_L(j,cc)=0;AgeL_L(j,cc)=0;e_rel_L(j,cc)=0; e_relN_L(j,cc)=0;
                SAI_L(j,cc)=0; hc_L(j,cc)=0;
                LAIdead_L(j,cc) =LAIdead_L(j-1,cc);
                Sr_L(j,cc)=0;Slf_L(j,cc)=0;Sfr_L(j,cc)=0;Swm_L(j,cc)=0;  Rrootl_L(j,cc)=0;
                Rexmy_L(j,cc,:)=0;AgeDL_L(j,cc)=0;FNC_L(j,cc)=1;Nreserve_L(j,cc)=0;Preserve_L(j,cc)=0;Kreserve_L(j,cc)=0;
                rNc_L(j,cc)=1;rPc_L(j,cc)=1;rKc_L(j,cc)=1;
                Bfac_dayL(j,cc)=0; Bfac_weekL(j,cc)=0; NPPI_L(j,cc)=0; TdpI_L(j,cc)=0;  NupI_L(j,cc,:)=0; RB_L(j,cc,:)=0;
                TexC_L(j,cc)=0;TexN_L(j,cc)=0;TexP_L(j,cc)=0;TexK_L(j,cc)=0;TNIT_L(j,cc)=0;TPHO_L(j,cc)=0;TPOT_L(j,cc)=0;
                ISOIL_L(j,cc,:)=0; NuLit_L(j,cc,:)=0;
            end
            %%%%%%%%%%%%%%%%%
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if OPT_VCA == 0  %%% Solutions for Crops aSE==5 
            Ccrown_t(j,:) = Ccrown;
            Cbare = 1 - sum(Ccrown) - Cwat - Curb - Crock;
            Check_Land_Cover_Fractions(Crock,Curb,Cwat,Cbare,Ccrown);
        end
    end
    %%%%%%% %%%%%%%

    %%%%%%% %%%%%%%
    if OPT_VCA == 1
        if((OPT_ALLOME == 1)&& (ZR95_L(2))>0) %%% Ad hoc solution for oil palm
            Ccrown_t(j,2) = 1-Ccrown_t(j,1);  % understory of Oil Palm
            Ccrown(2) = Ccrown_t(j,2);
            %%% Theoretically to preserve mass // but this is considered
            %%% lost understory in the  oil palm 
            %B_L(j,cc,:)= B_L(j,cc,:)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
            %Nreserve_L(j,cc)=Nreserve_L(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
            %Preserve_L(j,cc)=Preserve_L(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
            %Kreserve_L(j,cc)=Kreserve_L(j,cc)*Ccrown_t(j-1,cc)/Ccrown_t(j,cc);
        end
        Cbare = 1 - sum(Ccrown);
        if zatm < (max(max(hc_H),max(hc_L))+2)
            zatm = zatm + 2;
        end
    end
    
    Qi_in(i,:)=Qi_out(i-1,:);
    q_runon(i)=0; %Rd(i-1)+Rh(i-1);
    if OPT_WET == 1 %%% Wetland option with water standing in the surface
        q_runon(i)=Wlev(i);
        if isnan(Wlev(i))
            q_runon(i)=Wlevm1;
        end
        NetWatWet(i) = q_runon(i) - Wlevm1 ;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [V(i,:),Vice(i,:),O(i,:),Oice(i,:),ZWT(i),OF(i),OS(i),OH(i,:),OL(i,:),Psi_s_H(i,:),Psi_s_L(i,:),Rd(i),Qi_out(i,:),WTR(i,:),...
        Rh(i),Lk(i),f(i),WIS(i),Ts(i),Pr_sno(i),Pr_liq(i),Csno(i),Cice(i),NDVI(i),rb_H(i,:),rb_L(i,:),rs_sunH(i,:),...
        rs_sunL(i,:),rs_shdH(i,:),rs_shdL(i,:),r_litter(i,:),...
        An_L(i,:),An_H(i,:),Rdark_L(i,:),Rdark_H(i,:),Ci_sunH(i,:),Ci_sunL(i,:),Ci_shdH(i,:),Ci_shdL(i,:),...
        rap_H(i,:),rap_L(i,:),r_soil(i),b_soil(i),alp_soil(i),ra(i),Rn(i),...
        H(i),QE(i),Qv(i),Lpho(i),T_H(i,:),T_L(i,:),EIn_H(i,:),EIn_L(i,:),EG(i),ESN(i),ESN_In(i),ELitter(i),EWAT(i),EICE(i),EIn_urb(i),EIn_rock(i),dw_SNO(i),...
        G(i),Gfin(i),Tdp(i,:),Tdpsnow(i,:),Tdeb(i,:),Tdamp(i),Tice(i),Tdp_H(i,:),Tdp_L(i,:),SWE(i),SND(i),ros(i),In_SWE(i),SP_wc(i),WR_SP(i),U_SWE(i),NIn_SWE(i),dQ(i),Qfm(i),t_sls(i),DQ(i),DT(i),...
        WAT(i),ICE(i),ICE_D(i),IP_wc(i),WR_IP(i),NIce(i),Cicew(i),Csnow(i),FROCK(i),Imelt(i),Smelt(i),...
        In_H(i,:),In_L(i,:),In_Litter(i),In_urb(i),In_rock(i),Dr_H(i,:),Dr_L(i,:),SE_rock(i),SE_urb(i),Lk_wat(i),Lk_rock(i),er(i),...
        gsr_H(i,:),Psi_x_H(i,:),Psi_l_H(i,:),Jsx_H(i,:),Jxl_H(i,:),Kleaf_H(i,:),Kx_H(i,:),Vx_H(i,:),Vl_H(i,:),...
        gsr_L(i,:),Psi_x_L(i,:),Psi_l_L(i,:),Jsx_L(i,:),Jxl_L(i,:),Kleaf_L(i,:),Kx_L(i,:),Vx_L(i,:),Vl_L(i,:),...
        fapar_H(i,:),fapar_L(i,:),SIF_H(i,:),SIF_L(i,:),...
        snow_alb,tau_sno(i),e_sno(i),Ws_under(i),dQVEG(i),HV(i),QEV(i),TsV(i),Ts_under(i),EK(i),POT(i,:)]=HYDROLOGIC_UNIT(V(i-1,:),Oice(i-1,:),aR,Zs,...
        EvL_Zs,Inf_Zs,Zinf,RfH_Zs,RfL_Zs,dz,Dz,ms,Kbot,Pr(i),Ta(i),Ds(i),Ws(i),zatm,Ts(i-1),Ts_under(i-1),IrD(i),dt,dth,ea(i),N(i),Pre(i),Tstm0,...
        LAI_H(j,:),SAI_H(j,:),LAI_L(j,:),SAI_L(j,:),LAIdead_H(j,:),LAIdead_L(j,:),Rrootl_H(j,:),Rrootl_L(j,:),BLit(j,:),Sllit,Kct,...
        Datam(i,:),DeltaGMT,Lon,Lat,t_bef,t_aft,...
        Ccrown,Cbare,Crock,Curb,Cwat,...
        Soil_Param,Interc_Param,SnowIce_Param,VegH_Param,VegL_Param,...
        Zs_deb,Deb_Par,...
        ZR95_H,ZR95_L,...
        SAB1(i),SAB2(i),SAD1(i),SAD2(i),PARB(i),PARD(i),SvF,...
        SND(i-1),snow_alb,Color_Class,OM_H,OM_L,...
        PFT_opt_H,PFT_opt_L,hc_H(j,:),hc_L(j,:),d_leaf_H,d_leaf_L,...
        Ca(i),Oa,Ci_sunH(i-1,:),Ci_shdH(i-1,:),Ci_sunL(i-1,:),Ci_shdL(i-1,:),...
        e_rel_H(j,:),e_relN_H(j,:),e_rel_L(j,:),e_relN_L(j,:),...
        e_sno(i-1),In_H(i-1,:),In_L(i-1,:),In_Litter(i-1),In_urb(i-1),In_rock(i-1),SWE(i-1),In_SWE(i-1),....
        Tdeb(i-1,:),Tdp(i-1,:),Tdpsnow(i-1,:),Tdamp(i-1),Tice(i-1),...
        WAT(i-1),ICE(i-1),IP_wc(i-1),ICE_D(i-1),Cicew(i-1),...
        Vx_H(i-1,:),Vl_H(i-1,:),Vx_L(i-1,:),Vl_L(i-1,:),Psi_x_H(i-1,:),Psi_l_H(i-1,:),Psi_x_L(i-1,:),Psi_l_L(i-1,:),...
        FROCK(i-1),Krock,Ws_under(i-1),...
        Tdew(i),t_sls(i-1),ros(i-1),SP_wc(i-1),fpr,Pr_sno(pdind),...
        Urb_Par,In_max_urb,In_max_rock,K_usle,tau_sno(i-1),Ta(pdind),Slo_top,Slo_pot,Asur,Ared,aTop,EK(i-1),q_runon(i),Qi_in(i,:),...
        pow_dis,a_dis,Salt(i),...
        SPAR,SN,OPT_min_SPD,OPT_VegSnow,OPT_SoilTemp,OPT_PlantHydr,Opt_CR,Opt_ST,Opt_ST2,OPT_SM,OPT_STh,OPT_FR_SOIL,OPT_PH);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    Tstm0 =2*Ts(i)-Ts(i-1);
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    STOT= SAB1(i)+SAB2(i)+SAD1(i)+SAD2(i);
    ALB(i)= SAB1(i)/STOT*snow_alb.dir_vis + SAD1(i)/STOT*snow_alb.dif_vis + ...
        SAB2(i)/STOT*snow_alb.dir_nir + SAD2(i)/STOT*snow_alb.dif_nir;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if OPT_WET == 1 
        Wlevm1 = Rd(i)+Rh(i); 
    end 
    %%% v-coordinate
    CK1(i) = f(i)*dth*Asur*Ared + sum(V(i-1,:) - V(i,:))*Asur*Ared + sum(Vice(i-1,:) - Vice(i,:))*Asur*Ared - EG(i)*dth - Lk(i)*dth ...
        - sum(Qi_out(i,:))*dth -Rd(i) -sum(Jsx_L(i,:)).*dth -sum(Jsx_H(i,:)).*dth  + sum(Qi_in(i,:))*dth  ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%close(bau)
Computational_Time =toc;
%profile off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('COMPUTATIONAL TIME [h] ')
disp(Computational_Time/3600)
disp(' COMPUTATIONAL TIME [ms/cycle] ')
disp(1000*Computational_Time/NN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PROF1 = profile('info');
%profile('status')
%profile viewer
if OPT_SoilBiogeochemistry == 1
    NEE = -(NPP_H+NPP_L)*Ccrown' + R_litter + R_microbe + R_ew; %% [gC/m2 day]
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% MASS BALANCE CHECK
dV= (V(1,:)-V(end,:))*Asur*Ared +(Vice(1,:)-Vice(end,:))*Asur*Ared;
Ck = sum((Pr_liq(2:end)+Pr_sno(2:end))*dth)+ sum(IrD*dth) + sum(dV)...
    -sum(EG*dth)  -sum(ELitter*dth) -sum(sum(T_L)*dth)  -sum(sum(T_H)*dth) -sum(sum(EIn_L)*dth)  -sum(sum(EIn_H)*dth) -sum(Rh)...
    -sum(sum(Qi_out*dth)) -sum(Rd) + sum(sum(Qi_in*dth)) + sum(q_runon*dth)  - sum(EWAT*dth) ...
    - sum(ESN*dth) -sum(ESN_In*dth) -sum(EIn_urb*dth)- sum(EIn_rock*dth) + ...
    (SWE(1)-SWE(end)) + (In_SWE(1) - In_SWE(end)) + (SP_wc(1) -SP_wc(end))  + ...
    sum(In_H(1,:) - In_H(end,:)) + sum((In_L(1,:) - In_L(end,:)))+ (In_Litter(1) - In_Litter(end)) + (In_urb(1) - In_urb(end)) + (In_rock(1) - In_rock(end)) + ...
    - sum(EICE*dth) + (IP_wc(1) -IP_wc(end)) + (ICE(1) -ICE(end)) + (WAT(1) -WAT(end)) + (FROCK(1) -FROCK(end)) + ...
    Asur*( (Vx_H(1,:)-Vx_H(end,:)) + (Vl_H(1,:)-Vl_H(end,:))  + (Vx_L(1,:)-Vx_L(end,:)) + (Vl_L(1,:)-Vl_L(end,:)) )*Ccrown' ; %%%[mm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp(Ck);
disp(mean(DQ));
% T_H, T_L  EG, EIn_urb, EIn_rock, [mm/h]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% ENERGY TRANSFORMATION CHECK
Laten= 1000*(2501.3 - 2.361*(Ta)); %%%
ETen = (QE)*1000*3600./(reshape(Laten,size(QE))*1000);  %% [mm/h]
ET =  sum(T_H+EIn_H,2) + sum(T_L+EIn_L,2) +  EG +  ELitter + ESN + ESN_In + EWAT +  EICE+ EIn_urb + EIn_rock ;  %% [mm/h]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%% CARBON BALANCE CHECK
for kj=1:cc
    if RA_H(end,kj) == 0
        dB_H= squeeze((B_H(1,kj,:)-B_H(end-1,kj,:)));
    else
        dB_H= squeeze((B_H(1,kj,:)-B_H(end,kj,:)));
    end
    if RA_L(end,kj) == 0
        dB_L= squeeze((B_L(1,kj,:)-B_L(end-1,kj,:)));
    else
        dB_L= squeeze((B_L(1,kj,:)-B_L(end,kj,:)));
    end
    %%%%%%%%
 %CkC_H(kj) =  sum(dB_H)+ sum(NPP_H(:,j)) - sum(Swm_H(:,j))-sum(Sfr_H(:,j))-sum(Sr_H(:,j))-sum(Slf_H(:,j)) - sum(Rexmy_H(:,j,:));%
    CkC_H(kj) =  sum(dB_H)+ sum(1.0368*An_H(:,kj)/24)- sum(Rmr_H(:,kj)) -sum(Rmc_H(:,kj)) -sum(Rms_H(:,kj)) -sum(Rg_H(:,kj))...
        -sum(TexC_H(:,kj));
    %- sum(Swm_H(:,kj))-sum(Sfr_H(:,kj))-sum(Sr_H(:,kj))-sum(Slf_H(:,kj)) - sum(Rexmy_H(:,kj,:));%
    % CkC_H(kj) =  sum(dB_H)+ sum(1.0368*An_H(:,kj)/24)- sum(Rmr_H(:,kj)) -sum(Rmc_H(:,kj)) -sum(Rms_H(:,kj)) -sum(Rg_H(:,kj))...
    %     -sum(RB_H(:,kj,:),[1,3])- sum(Swm_H(:,kj))-sum(Sfr_H(:,kj))-sum(Sr_H(:,kj))-sum(Slf_H(:,kj))-sum(Rexmy_H(:,kj,:),[1,3]);
    CkC_L(kj) =  sum(dB_L)+ sum(1.0368*An_L(:,kj)/24)- sum(Rmr_L(:,kj)) -sum(Rmc_L(:,kj)) -sum(Rms_L(:,kj)) -sum(Rg_L(:,kj))...
        -sum(TexC_L(:,kj));
    %- sum(Swm_L(:,kj))-sum(Sfr_L(:,kj))-sum(Sr_L(:,kj))-sum(Slf_L(:,kj)) - sum(Rexmy_L(:,kj,:));%
    %    CkC_L(kj) =  sum(dB_L)+ sum(1.0368*An_L(:,kj)/24)- sum(Rmr_L(:,kj)) -sum(Rmc_L(:,kj)) -sum(Rms_L(:,kj)) -sum(Rg_L(:,kj))...
    %    -sum(RB_L(:,kj,:),[1,3])- sum(Swm_L(:,kj))-sum(Sfr_L(:,kj))-sum(Sr_L(:,kj))-sum(Slf_L(:,kj))-sum(Rexmy_L(:,kj,:),[1,3]);
end
CkC_ALL = sum(CkC_H)+sum(CkC_L);

for kj=1:cc
    %%%%%%%% NUTRIENT BALANCE CHECK
    if RA_H(end,kj) == 0
        ed=length(Nreserve_H)-1;
    else
        ed=length(Nreserve_H);
    end
    dNres_H=Nreserve_H(1,kj) - Nreserve_H(ed,kj);
    CkN_H(kj)= dNres_H -  sum(TexN_H(:,kj))  + (TNIT_H(2,kj) - TNIT_H(ed,kj)) + sum(Nuptake_H(:,kj));
    dPres_H=Preserve_H(1,kj) - Preserve_H(ed,kj);
    CkP_H(kj)= dPres_H -  sum(TexP_H(:,kj))  + (TPHO_H(2,kj) - TPHO_H(ed,kj)) + sum(Puptake_H(:,kj));
    dKres_H=Kreserve_H(1,kj) - Kreserve_H(ed,kj);
    CkK_H(kj)= dKres_H -  sum(TexK_H(:,kj))  + (TPOT_H(2,kj) - TPOT_H(ed,kj)) + sum(Kuptake_H(:,kj));
    %%%%
    dNres_L=Nreserve_L(1,kj) - Nreserve_L(ed,kj);
    CkN_L(kj)= dNres_L -  sum(TexN_L(:,kj))  + (TNIT_L(2,kj) - TNIT_L(ed,kj)) + sum(Nuptake_L(:,kj));
    dPres_L=Preserve_L(1,kj) - Preserve_L(ed,kj);
    CkP_L(kj)= dPres_L -  sum(TexP_L(:,kj))  + (TPHO_L(2,kj) - TPHO_L(ed,kj)) + sum(Puptake_L(:,kj));
    dKres_L=Kreserve_L(1,kj) - Kreserve_L(ed,kj);
    CkK_L(kj)= dKres_L -  sum(TexK_L(:,kj))  + (TPOT_L(2,kj) - TPOT_L(ed,kj)) + sum(Kuptake_L(:,kj));
end

%%%% Difference between carbon and nutrient ISOIL and Tex (0 except for
%%%% management cases)
%%% To add LitFirEmi to the litter budget
IS = zeros(18, 1);
Tex = 0;
for jj = 1:length(Ccrown_t)
    IS = IS + squeeze(ISOIL_L(jj, :, :))'*Ccrown_t(jj, :)' + squeeze(ISOIL_H(jj, :, :))'*Ccrown_t(jj, :)';
    Tex = Tex + TexC_L(jj,:)*Ccrown_t(jj, :)' + TexC_H(jj,:)*Ccrown_t(jj, :)';
end

CkExC = Tex - sum(IS(1:9));

clear dB_H  dB_L dNres_H  dNres_L  dPres_H dKres_H  dPres_L dKres_L ed
clear Se_bio Psi_bio Tdp_bio VSUM VTSUM  IS Tex 
clear  Tstm0 snow_alb  kj m p  pdind STOT  PARAM_IC
clear Bam Bem FireA Se_fc REXMY

%Ckf=0*Ts;
%for i=2:NN
%dV= V(1,:)-V(i,:) +(Vice(1,:)-Vice(i,:));
%Ckf(i) = sum((Pr_liq(1:i)+Pr_sno(1:i))*dth) + sum(dV)...
%    -sum(EG(1:i)*dth)  -sum(ELitter(1:i)*dth) -sum(sum(T_L(1:i))*dth)  -sum(sum(T_H(1:i))*dth) -sum(sum(EIn_L(1:i))*dth)  -sum(sum(EIn_H(1:i))*dth) -sum(Rh(1:i))...
%    -sum(sum(Qi_out(1:i,:)*dth)) -sum(Rd(1:i)) + sum(sum(Qi_in(1:i,:)*dth)) + sum(q_runon(1:i)*dth) - sum(EWAT(1:i)*dth)...
%    - sum(ESN(1:i)*dth) -sum(ESN_In(1:i)*dth) -sum(EIn_urb(1:i)*dth)- sum(EIn_rock(1:i)*dth) + ...
%    (SWE(1)-SWE(i)) + (In_SWE(1) - In_SWE(i)) + (SP_wc(1) -SP_wc(i))  + ...
%    sum(In_H(1,:) - In_H(i,:)) + sum((In_L(1,:) - In_L(i,:))) + (In_urb(1) - In_urb(i)) + (In_Litter(1) - In_Litter(i))  + (In_rock(1) - In_rock(i))+...
%    - sum(EICE(1:i)*dth) + (IP_wc(1) -IP_wc(i)) + (ICE(1) -ICE(i)) + (WAT(1) -WAT(i)) + (FROCK(1) -FROCK(i)) + ...
%    Asur*( (Vx_H(1,:)-Vx_H(i,:)) + (Vl_H(1,:)-Vl_H(i,:))  + (Vx_L(1,:)-Vx_L(i,:)) + (Vl_L(1,:)-Vl_L(i,:)) )*Ccrown'; %%%[mm]
%end


%for i=2:NN
%Ck1(i)= (T_L(i)-Jsx_L(i))+(Vl_L(i)-Vl_L(i-1)) ;
%end


%for i=2:NN
%Cc(i)=Jsx_L(i)+Lk(i)+EG(i) +( sum(V(i,:))-sum(V(i-1,:))) -f(i);
%end

