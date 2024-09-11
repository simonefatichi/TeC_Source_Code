%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% MAIN_FRAME SPATIAL CHLORIS %%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[LAI_H,B_H,NPP_H,ANPP_H,Rg_H,RA_H,Rms_H,Rmr_H,Rmc_H,PHE_S_H,...
    dflo_H,AgeL_H,e_rel_H,e_relN_H,...
    LAI_L,B_L,NPP_L,ANPP_L,Rg_L,RA_L,Rms_L,Rmr_L,Rmc_L,PHE_S_L,...
    dflo_L,AgeL_L,e_rel_L,e_relN_L,...
    SAI_H,hc_H,SAI_L,hc_L,...
    LAIdead_H,NBLeaf_H,Sr_H,Slf_H,Sfr_H,Sll_H,Swm_H,Rexmy_H,NupI_H,NuLit_H,...
    LAIdead_L,NBLeaf_L,Sr_L,Slf_L,Sfr_L,Sll_L,Swm_L,Rexmy_L,NupI_L,NuLit_L,...
    Rrootl_H,AgeDL_H,Bfac_dayH,Bfac_weekH,NPPI_H,TdpI_H,PARI_H,NBLI_H,RB_H,FNC_H,Nreserve_H,Preserve_H,Kreserve_H,rNc_H,rPc_H,rKc_H,ManIH,...
    Rrootl_L,AgeDL_L,Bfac_dayL,Bfac_weekL,NPPI_L,TdpI_L,PARI_L,NBLI_L,RB_L,FNC_L,Nreserve_L,Preserve_L,Kreserve_L,rNc_L,rPc_L,rKc_L,ManIL,...
    TexC_H,TexN_H,TexP_H,TexK_H,TNIT_H,TPHO_H,TPOT_H,SupN_H,SupP_H,SupK_H,ISOIL_H,...
    TexC_L,TexN_L,TexP_L,TexK_L,TNIT_L,TPHO_L,TPOT_L,SupN_L,SupP_L,SupK_L,ISOIL_L,...
    BA_H,Tden_H,AgePl_H,BA_L,Tden_L,AgePl_L,Ccrown_t]=VEGETATION_MODULE_PAR(cc_max,Ccrown,ZR_H,ZR_L,B_Htm1,...
    PHE_S_Htm1,dflo_Htm1,AgeL_Htm1,AgeDL_Htm1,...
    Ta24,PAR24,Tdp_H24,Psi_x_H24,Psi_l_H24,An_H24,Rdark_H24,NPP_Htm1,jDay,Datam,...
    NPPI_Htm1,TdpI_Htm1,Bfac_weekHtm1,...
    Stoich_H,aSE_H,VegH_Param_Dyn,...
    Nreserve_Htm1,Preserve_Htm1,Kreserve_Htm1,Nuptake_H,Puptake_H,Kuptake_H,FNC_Htm1,Tden_Htm1,AgePl_Htm1,...
    fab_H,fbe_H,ParEx_H,Mpar_H,TBio_H,SAI_Htm1,hc_Htm1,...
    B_Ltm1,PHE_S_Ltm1,dflo_Ltm1,AgeL_Ltm1,AgeDL_Ltm1,...
    Tdp_L24,Psi_x_L24,Psi_l_L24,An_L24,Rdark_L24,NPP_Ltm1,...
    NPPI_Ltm1,TdpI_Ltm1,Bfac_weekLtm1,...
    NupI_Htm1,NupI_Ltm1,NuLit_Htm1,NuLit_Ltm1,NBLeaf_Htm1,NBLeaf_Ltm1,...
    PARI_Htm1,NBLI_Htm1,PARI_Ltm1,NBLI_Ltm1,...
    Stoich_L,aSE_L,VegL_Param_Dyn,...
    NavlI,Bam,Bem,Ccrown_t_tm1,...
    Nreserve_Ltm1,Preserve_Ltm1,Kreserve_Ltm1,Nuptake_L,Puptake_L,Kuptake_L,FNC_Ltm1,Tden_Ltm1,AgePl_Ltm1,...
    fab_L,fbe_L,ParEx_L,Mpar_L,TBio_L,SAI_Ltm1,hc_Ltm1,...
    ExEM,Lmax_day,L_day,Se_bio,Tdp_bio,OPT_EnvLimitGrowth,OPT_VD,OPT_VCA,OPT_ALLOME,OPT_SoilBiogeochemistry)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% %% To be adjusted
dtd=1; %% day
Ccrown_t = Ccrown_t_tm1;
%%% currently not working from crop rotations Ccrown fixed on the external parameter file

%%%%%% VEGETATION MODULE
hc_Htm1 =squeeze(hc_Htm1);
hc_Ltm1 =squeeze(hc_Ltm1);
AgeL_Htm1 =squeeze(AgeL_Htm1);
AgeL_Ltm1 =squeeze(AgeL_Ltm1);
B_Htm1 =squeeze(B_Htm1);
B_Ltm1 =squeeze(B_Ltm1);
PHE_S_Htm1 =squeeze(PHE_S_Htm1);
PHE_S_Ltm1 =squeeze(PHE_S_Ltm1);
SAI_Htm1 =squeeze(SAI_Htm1);
SAI_Ltm1 =squeeze(SAI_Ltm1);
dflo_Htm1 =squeeze(dflo_Htm1);
dflo_Ltm1 =squeeze(dflo_Ltm1);
AgeDL_Htm1 =squeeze(AgeDL_Htm1);
NPP_Htm1 =squeeze(NPP_Htm1);
NPPI_Htm1 =squeeze(NPPI_Htm1);
TdpI_Htm1 =squeeze(TdpI_Htm1);
Bfac_weekHtm1 =squeeze(Bfac_weekHtm1);
Nreserve_Htm1 =squeeze(Nreserve_Htm1);
Preserve_Htm1 =squeeze(Preserve_Htm1);
Kreserve_Htm1 =squeeze(Kreserve_Htm1);
Nuptake_H =squeeze(Nuptake_H);
Puptake_H =squeeze(Puptake_H);
Kuptake_H =squeeze(Kuptake_H);
FNC_Htm1 =squeeze(FNC_Htm1);
AgeDL_Ltm1 =squeeze(AgeDL_Ltm1);
NPP_Ltm1 =squeeze(NPP_Ltm1);
NPPI_Ltm1 =squeeze(NPPI_Ltm1);
TdpI_Ltm1 =squeeze(TdpI_Ltm1);
Bfac_weekLtm1 =squeeze(Bfac_weekLtm1);
Nreserve_Ltm1 =squeeze(Nreserve_Ltm1);
Preserve_Ltm1 =squeeze(Preserve_Ltm1);
Kreserve_Ltm1 =squeeze(Kreserve_Ltm1);
Nuptake_L =squeeze(Nuptake_L);
Puptake_L =squeeze(Puptake_L);
Kuptake_L =squeeze(Kuptake_L);
FNC_Ltm1 =squeeze(FNC_Ltm1);
NupI_Htm1 =squeeze(NupI_Htm1);
NupI_Ltm1 = squeeze(NupI_Ltm1);
NuLit_Htm1 = squeeze(NuLit_Htm1);
NuLit_Ltm1 = squeeze(NuLit_Ltm1);
PARI_Htm1 = squeeze(PARI_Htm1);
PARI_Ltm1 = squeeze(PARI_Ltm1);
Tden_Htm1 = squeeze(Tden_Htm1);
Tden_Ltm1 = squeeze(Tden_Ltm1);
AgePl_Htm1 = squeeze(AgePl_Htm1);
AgePl_Ltm1 = squeeze(AgePl_Ltm1);
%%%%%%%%%%%
Ta24 =squeeze(Ta24);
PAR24 =squeeze(PAR24);
Tdp_H24 =squeeze(Tdp_H24);
Tdp_L24 =squeeze(Tdp_L24);
Rdark_H24 =squeeze(Rdark_H24);
Rdark_L24 =squeeze(Rdark_L24);
An_H24 =squeeze(An_H24);
An_L24 =squeeze(An_L24);
Psi_x_H24 =squeeze(Psi_x_H24);
Psi_l_H24 =squeeze(Psi_l_H24);
Psi_x_L24 =squeeze(Psi_x_L24);
Psi_l_L24 =squeeze(Psi_l_L24);
%%%%%%%%%%
if cc_max == 1
    An_H24=An_H24';
    An_L24=An_L24';
    B_Htm1=B_Htm1';
    B_Ltm1=B_Ltm1';
    Psi_x_H24=Psi_x_H24';
    Psi_l_H24=Psi_l_H24';
    Psi_x_L24=Psi_x_L24';
    Psi_l_L24=Psi_l_L24';
    Tdp_H24=Tdp_H24';
    Tdp_L24=Tdp_L24';
    Rdark_H24=Rdark_H24';
    Rdark_L24=Rdark_L24';
    NupI_Htm1 = NupI_Htm1';
    NupI_Ltm1 = NupI_Ltm1';
    NuLit_Htm1 = NuLit_Htm1';
    NuLit_Ltm1 =  NuLit_Ltm1';
    PARI_Htm1 = PARI_Htm1';
    PARI_Ltm1= PARI_Ltm1';
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NCP = 8; %% Number of Carbon Pool
%cc=length(Ccrown); %%% Number of PFT within cell
%cc_max ; %%% Max number of PFT
%%%%
LAI_H=zeros(1,cc_max); B_H=zeros(cc_max,NCP);NPP_H=zeros(1,cc_max);Rg_H=zeros(1,cc_max);
RA_H=zeros(1,cc_max);Rms_H=zeros(1,cc_max);Rmr_H=zeros(1,cc_max); ANPP_H=zeros(1,cc_max);
PHE_S_H=zeros(1,cc_max); dflo_H=zeros(1,cc_max);Rmc_H=zeros(1,cc_max);
AgeL_H=zeros(1,cc_max);e_rel_H=zeros(1,cc_max);e_relN_H=zeros(1,cc_max);SAI_H=zeros(1,cc_max); hc_H=zeros(1,cc_max);
LAIdead_H=zeros(1,cc_max); Sr_H=zeros(1,cc_max);  Sll_H=zeros(1,cc_max); Slf_H=zeros(1,cc_max);Sfr_H=zeros(1,cc_max);Swm_H=zeros(1,cc_max);
Rexmy_H=zeros(cc_max,3);Rrootl_H=zeros(1,cc_max);
AgeDL_H=zeros(1,cc_max);Bfac_dayH=zeros(1,cc_max);Bfac_weekH=zeros(1,cc_max);NPPI_H=zeros(1,cc_max);TdpI_H=zeros(1,cc_max);
RB_H=zeros(cc_max,7);FNC_H=zeros(1,cc_max);Nreserve_H=zeros(1,cc_max);Preserve_H=zeros(1,cc_max);Kreserve_H=zeros(1,cc_max);
TexC_H=zeros(1,cc_max);TexN_H=zeros(1,cc_max);TexP_H=zeros(1,cc_max);TexK_H=zeros(1,cc_max);
TNIT_H=zeros(1,cc_max);TPHO_H=zeros(1,cc_max);TPOT_H=zeros(1,cc_max);
SupN_H=zeros(1,cc_max);SupP_H=zeros(1,cc_max);SupK_H=zeros(1,cc_max);ISOIL_H=zeros(cc_max,18);
rNc_H=zeros(1,cc_max);rPc_H=zeros(1,cc_max);rKc_H=zeros(1,cc_max);ManIH=zeros(1,cc_max);
NupI_H=zeros(cc_max,3); NuLit_H=zeros(cc_max,3); NBLeaf_H=zeros(1,cc_max);
PARI_H=zeros(cc_max,3);NBLI_H=zeros(1,cc_max);
BA_H=zeros(1,cc_max); Tden_H=zeros(1,cc_max); AgePl_H=zeros(1,cc_max); TBio_Ht=zeros(1,cc_max);
%%%%
LAI_L=zeros(1,cc_max); B_L=zeros(cc_max,NCP);NPP_L=zeros(1,cc_max);Rg_L=zeros(1,cc_max);
RA_L=zeros(1,cc_max);Rms_L=zeros(1,cc_max);Rmr_L=zeros(1,cc_max); ANPP_L=zeros(1,cc_max);
PHE_S_L=zeros(1,cc_max); dflo_L=zeros(1,cc_max);Rmc_L=zeros(1,cc_max);
AgeL_L=zeros(1,cc_max);e_rel_L=zeros(1,cc_max);e_relN_L=zeros(1,cc_max);SAI_L=zeros(1,cc_max); hc_L=zeros(1,cc_max);
LAIdead_L=zeros(1,cc_max); Sr_L=zeros(1,cc_max);Sll_L=zeros(1,cc_max);  Slf_L=zeros(1,cc_max);Sfr_L=zeros(1,cc_max);Swm_L=zeros(1,cc_max);
Rexmy_L=zeros(cc_max,3);Rrootl_L=zeros(1,cc_max);
AgeDL_L=zeros(1,cc_max);Bfac_dayL=zeros(1,cc_max);Bfac_weekL=zeros(1,cc_max);NPPI_L=zeros(1,cc_max);TdpI_L=zeros(1,cc_max);
RB_L=zeros(cc_max,7);FNC_L=zeros(1,cc_max);Nreserve_L=zeros(1,cc_max);Preserve_L=zeros(1,cc_max);Kreserve_L=zeros(1,cc_max);
TexC_L=zeros(1,cc_max);TexN_L=zeros(1,cc_max);TexP_L=zeros(1,cc_max);TexK_L=zeros(1,cc_max);
TNIT_L=zeros(1,cc_max);TPHO_L=zeros(1,cc_max);TPOT_L=zeros(1,cc_max);
SupN_L=zeros(1,cc_max);SupP_L=zeros(1,cc_max);SupK_L=zeros(1,cc_max);ISOIL_L=zeros(cc_max,18);
rNc_L=zeros(1,cc_max);rPc_L=zeros(1,cc_max);rKc_L=zeros(1,cc_max);ManIL=zeros(1,cc_max);
NupI_L=zeros(cc_max,3);  NuLit_L=zeros(cc_max,3); NBLeaf_L=zeros(1,cc_max);
PARI_L=zeros(cc_max,3);NBLI_L=zeros(1,cc_max);
BA_L=zeros(1,cc_max); Tden_L=zeros(1,cc_max); AgePl_L=zeros(1,cc_max); TBio_Lt=zeros(1,cc_max);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% ,
for cc=1:length(Ccrown)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ZR_H(cc) > 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [LAI_H(cc),B_H(cc,:),NPP_H(cc),ANPP_H(cc),Rg_H(cc),RA_H(cc),Rms_H(cc),Rmr_H(cc),Rmc_H(cc),PHE_S_H(cc),...
            dflo_H(cc),AgeL_H(cc),e_rel_H(cc),e_relN_H(cc),LAIdead_H(cc),NBLeaf_H(cc),Sr_H(cc),Slf_H(cc),Sfr_H(cc),Swm_H(cc),Sll_H(cc),Rexmy_H(cc,:),Rrootl_H(cc),...
            AgeDL_H(cc),Bfac_dayH(cc),Bfac_weekH(cc),NPPI_H(cc),TdpI_H(cc),NupI_H(cc,:),PARI_H(cc,:),NBLI_H(cc),RB_H(cc,:),FNC_H(cc),Nreserve_H(cc),Preserve_H(cc),Kreserve_H(cc),...
            rNc_H(cc),rPc_H(cc),rKc_H(cc),ManIH(cc)]= VEGGIE_UNIT(B_Htm1(cc,:),PHE_S_Htm1(cc),dflo_Htm1(cc),AgeL_Htm1(cc),AgeDL_Htm1(cc),...
            Ta24,Tdp_H24(cc,:),PAR24,Psi_x_H24(cc,:),Psi_l_H24(cc,:),An_H24(cc,:),Rdark_H24(cc,:),NPP_Htm1(cc),jDay,Datam,...
            NPPI_Htm1(cc),TdpI_Htm1(cc),Bfac_weekHtm1(cc),NupI_Htm1(cc,:),NavlI,PARI_Htm1(cc,:),NBLI_Htm1(cc),NBLeaf_Htm1(cc),...
            L_day,Lmax_day,VegH_Param_Dyn,cc,...
            Nreserve_Htm1(cc),Preserve_Htm1(cc),Kreserve_Htm1(cc),Nuptake_H(cc),Puptake_H(cc),Kuptake_H(cc),FNC_Htm1(cc),Se_bio,Tdp_bio,...
            ParEx_H(cc),ExEM,Bam,Bem,Mpar_H(cc),TBio_H(cc),OPT_EnvLimitGrowth,OPT_VCA,OPT_VD,OPT_SoilBiogeochemistry);

        %%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%
        [TexC_H(cc),TexN_H(cc),TexP_H(cc),TexK_H(cc),TNIT_H(cc),TPHO_H(cc),TPOT_H(cc),NuLit_H(cc,:),Nreserve_H(cc),Preserve_H(cc),Kreserve_H(cc),...
            SupN_H(cc),SupP_H(cc),SupK_H(cc),ISOIL_H(cc,:)]= Plant_Exports(B_H(cc,:),B_Htm1(cc,:),NuLit_Htm1(cc,:),...
            Slf_H(cc),Sfr_H(cc),Swm_H(cc),Sll_H(cc),Sr_H(cc),Rexmy_H(cc,:),Stoich_H(cc),Mpar_H(cc),fab_H(cc),fbe_H(cc),RB_H(cc,:),Nreserve_H(cc),Preserve_H(cc),Kreserve_H(cc),...
            rNc_H(cc),rPc_H(cc),rKc_H(cc),ManIH(cc),OPT_SoilBiogeochemistry);
        %%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%% Change in height SAI and Ccrown
        SAI_H(cc) =SAI_Htm1(cc);
        %%%%
        if aSE_H(cc) == 2
            [hc_H(cc)] = GrassHeight(LAI_H(cc),LAIdead_H(cc));

        elseif aSE_H(cc) == 5

            [hc_H(cc),SAI_H(cc),B_H(:,:),~,Nreserve_H(cc),Preserve_H(cc),Kreserve_H(cc),AgrHarNut(:)] = CropHeightType(LAI_H(cc),LAIdead_H(cc),cc,B_H(:,:),...
                Ccrown,Nreserve_H(cc),Preserve_H(cc),Kreserve_H(cc),ManIH,Mpar_H,VegH_Param_Dyn,OPT_SoilBiogeochemistry);

        else
            hc_H(cc)= hc_Htm1(cc); %%%[m]
            if OPT_VCA == 1
                %%%%%%
                [Ccrown_t(cc),hc_H(cc),SAI_H(cc),BA_H(cc),Tden_H(cc),AgePl_H(cc),TBio_Ht(cc)]=Vegetation_Structural_Attributes(dtd,...
                    Ccrown_t_tm1(cc),B_H(cc,:),fab_H(cc),Tden_Htm1(cc),AgePl_Htm1(cc),OPT_ALLOME);
                %%%%
                Ccrown(cc) = CcrownFIX(cc)*Ccrown_t(cc);
                TBio_H(cc) = TBio_Ht(cc);
                B_H(cc,:)= B_H(cc,:)*Ccrown_t_tm1(cc)/Ccrown_t(cc);
                Nreserve_H(cc)=Nreserve_H(cc)*Ccrown_t_tm1(cc)/Ccrown_t(cc);
                Preserve_H(cc)=Preserve_H(cc)*Ccrown_t_tm1(cc)/Ccrown_t(cc);
                Kreserve_H(cc)=Kreserve_H(cc)*Ccrown_t_tm1(cc)/Ccrown_t(cc);
            end
        end

    else
        LAI_H(cc)=0;B_H(cc,:)=0;NPP_H(cc)=0;ANPP_H(cc)=0;Rg_H(cc)=0;RA_H(cc)=0;Rms_H(cc)=0;
        Rmr_H(cc)=0;PHE_S_H(cc)=0;dflo_H(cc)=0;AgeL_H(cc)=0;e_rel_H(cc)=0; e_relN_H(cc)=0;
        SAI_H(cc)=0; hc_H(cc)=0;
        LAIdead_H(cc) =0; NupI_H(cc,:)=0; ManIH(cc)=0;
        Sr_H(cc)=0;Slf_H(cc)=0;Sll_H(cc)=0;Sfr_H(cc)=0;Swm_H(cc)=0; Rrootl_H(cc)=0;
        Rexmy_H(cc,:)=0;AgeDL_H(cc)=0;FNC_H(cc)=1; Nreserve_H(cc)=0;Preserve_H(cc)=0;Kreserve_H(cc)=0;
        rNc_H(cc)=0;rPc_H(cc)=0;rKc_H(cc)=0;
        Bfac_dayH(cc)=0; Bfac_weekH(cc)=0; NPPI_H(cc)=0; TdpI_H(cc)=0;
        TexC_H(cc)=0;TexN_H(cc)=0;TexP_H(cc)=0;TexK_H(cc)=0;TNIT_H(cc)=0;TPHO_H(cc)=0;TPOT_H(cc)=0;
        SupN_H(cc)=1;SupP_H(cc)=1;SupK_H(cc)=1;ISOIL_H(cc,:)=0;
        BA_H(cc) = 0;Tden_H(cc) = 0;  AgePl_H(cc) = 0;
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ZR_L(cc) > 0
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        [LAI_L(cc),B_L(cc,:),NPP_L(cc),ANPP_L(cc),Rg_L(cc),RA_L(cc),Rms_L(cc),Rmr_L(cc),Rmc_L(cc),PHE_S_L(cc),...
            dflo_L(cc),AgeL_L(cc),e_rel_L(cc),e_relN_L(cc),LAIdead_L(cc),NBLeaf_L(cc),Sr_L(cc),Slf_L(cc),Sfr_L(cc),Swm_L(cc),Sll_L(cc),Rexmy_L(cc,:),Rrootl_L(cc),...
            AgeDL_L(cc),Bfac_dayL(cc),Bfac_weekL(cc),NPPI_L(cc),TdpI_L(cc),NupI_L(cc,:),PARI_L(cc,:),NBLI_L(cc),RB_L(cc,:),FNC_L(cc),Nreserve_L(cc),Preserve_L(cc),Kreserve_L(cc),...
            rNc_L(cc),rPc_L(cc),rKc_L(cc),ManIL(cc)]= VEGGIE_UNIT(B_Ltm1(cc,:),PHE_S_Ltm1(cc),dflo_Ltm1(cc),AgeL_Ltm1(cc),AgeDL_Ltm1(cc),...
            Ta24,Tdp_L24(cc,:),PAR24,Psi_x_L24(cc,:),Psi_l_L24(cc,:),An_L24(cc,:),Rdark_L24(cc,:),NPP_Ltm1(cc),jDay,Datam,...
            NPPI_Ltm1(cc),TdpI_Ltm1(cc),Bfac_weekLtm1(cc),NupI_Ltm1(cc,:),NavlI,PARI_Ltm1(cc,:),NBLI_Ltm1(cc),NBLeaf_Ltm1(cc),...
            L_day,Lmax_day,VegL_Param_Dyn,cc,...
            Nreserve_Ltm1(cc),Preserve_Ltm1(cc),Kreserve_Ltm1(cc),Nuptake_L(cc),Puptake_L(cc),Kuptake_L(cc),FNC_Ltm1(cc),Se_bio,Tdp_bio,...
            ParEx_L(cc),ExEM,Bam,Bem,Mpar_L(cc),TBio_L(cc),OPT_EnvLimitGrowth,OPT_VCA,OPT_VD,OPT_SoilBiogeochemistry);
        %%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%
        [TexC_L(cc),TexN_L(cc),TexP_L(cc),TexK_L(cc),TNIT_L(cc),TPHO_L(cc),TPOT_L(cc),NuLit_L(cc,:),Nreserve_L(cc),Preserve_L(cc),Kreserve_L(cc),...
            SupN_L(cc),SupP_L(cc),SupK_L(cc),ISOIL_L(cc,:)]= Plant_Exports(B_L(cc,:),B_Ltm1(cc,:),NuLit_Ltm1(cc,:),...
            Slf_L(cc),Sfr_L(cc),Swm_L(cc),Sll_L(cc),Sr_L(cc),Rexmy_L(cc,:),Stoich_L(cc),Mpar_L(cc),fab_L(cc),fbe_L(cc),RB_L(cc,:),Nreserve_L(cc),Preserve_L(cc),Kreserve_L(cc),...
            rNc_L(cc),rPc_L(cc),rKc_L(cc),ManIL(cc),OPT_SoilBiogeochemistry);
        %%%%%%%%%%%%%%%%,
        %%%%%%%%%%%%%%%%%%%%%% Change in height SAI and Ccrown
        SAI_L(cc) =SAI_Ltm1(cc);
        %%%%
        if aSE_L(cc) == 2
            [hc_L(cc)] = GrassHeight(LAI_L(cc),LAIdead_L(cc));
        elseif aSE_L(cc) == 5

            [hc_L(cc),SAI_L(cc),B_L(:,:),~,Nreserve_L(cc),Preserve_L(cc),Kreserve_L(cc),AgrHarNut(:)] = CropHeightType(LAI_L(cc),LAIdead_L(cc),cc,B_L(:,:),...
                Ccrown,Nreserve_L(cc),Preserve_L(cc),Kreserve_L(cc),ManIL,Mpar_L,VegL_Param_Dyn,OPT_SoilBiogeochemistry);

        else
            hc_L(cc)= hc_Ltm1(cc); %%%[m]
            if OPT_VCA == 1
                %%%%%%
                [Ccrown_t(cc),hc_L(cc),SAI_L(cc),BA_L(cc),Tden_L(cc),AgePl_L(cc),TBio_Lt(cc)]=Vegetation_Structural_Attributes(dtd,...
                    Ccrown_t_tm1(cc),B_L(cc,:),fab_L(cc),Tden_Ltm1(cc),AgePl_Ltm1(cc),OPT_ALLOME);
                %%%%
                Ccrown(cc) = CcrownFIX(cc)*Ccrown_t(cc);
                TBio_L(cc) = TBio_Lt(cc);
                B_L(cc,:)= B_L(cc,:)*Ccrown_t_tm1(cc)/Ccrown_t(cc);
                Nreserve_L(cc)=Nreserve_L(cc)*Ccrown_t_tm1(cc)/Ccrown_t(cc);
                Preserve_L(cc)=Preserve_L(cc)*Ccrown_t_tm1(cc)/Ccrown_t(cc);
                Kreserve_L(cc)=Kreserve_L(cc)*Ccrown_t_tm1(cc)/Ccrown_t(cc);
            end
        end

    else
        LAI_L(cc)=0;B_L(cc,:)=0;NPP_L(cc)=0;ANPP_L(cc)=0;Rg_L(cc)=0;RA_L(cc)=0;Rms_L(cc)=0;
        Rmr_L(cc)=0;PHE_S_L(cc)=0;dflo_L(cc)=0;AgeL_L(cc)=0;e_rel_L(cc)=0;e_relN_L(cc)=0;
        SAI_L(cc)=0; hc_L(cc)=0;
        LAIdead_L(cc) =0;NupI_L(cc,:)=0; ManIL(cc)=0;
        Sr_L(cc)=0;Slf_L(cc)=0;Sll_L(cc)=0;Sfr_L(cc)=0;Swm_L(cc)=0; Rrootl_L(cc)=0;
        Rexmy_L(cc,:)=0;AgeDL_L(cc)=0;FNC_L(cc)=1; Nreserve_L(cc)=0;Preserve_L(cc)=0;Kreserve_L(cc)=0;
        rNc_L(cc)=0;rPc_L(cc)=0;rKc_L(cc)=0;
        Bfac_dayL(cc)=0; Bfac_weekL(cc)=0; NPPI_L(cc)=0; TdpI_L(cc)=0;
        TexC_L(cc)=0;TexN_L(cc)=0;TexP_L(cc)=0;TexK_L(cc)=0;TNIT_L(cc)=0;TPHO_L(cc)=0;TPOT_L(cc)=0;
        SupN_L(cc)=1;SupP_L(cc)=1;SupK_L(cc)=1;ISOIL_L(cc,:)=0;
        BA_L(cc) = 0;Tden_L(cc) = 0;  AgePl_L(cc) = 0;
    end
end
return