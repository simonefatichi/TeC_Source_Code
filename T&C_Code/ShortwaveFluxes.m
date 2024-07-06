%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction ShortwaveFluxes              %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,PAR_sun_H,PAR_shd_H,...
    RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,PAR_sun_L,PAR_shd_L,FsunH,FshdH,...
    FsunL,FshdL,Kopt_H,Kopt_L,fapar_H,fapar_L,NDVI,ALB,Rabsb_sno,Rabsb_bare,Rabsb_urb,Rabsb_wat,Rabsb_rock,Rabsb_ice,Rabsb_deb,...
    soil_alb,e_gr,e_sur]=ShortwaveFluxes(Ccrown,Cbare,Crock,Curb,Cwat,Csno,Cice,...
    Rsw,PAR,SvF,dw_SNO,hc_H,hc_L,SnoDep,ydepth,IceDep,Cdeb,Deb_Par,Urb_Par,h_S,snow_alb,Aice,OS,Color_Class,OM_H,OM_L,...
    LAI_H,SAI_H,LAId_H,LAI_L,SAI_L,LAId_L,PFT_opt_H,PFT_opt_L)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT
%%%%% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% -- %%% alb Rsw --->  .dir_vis .dir_nir .dif_vis .dif_nir
%%%%%%% ALBEDOs Computation %%%%%%%%
[soil_alb,e_gr]=Albedo_Soil_Properties(OS,Color_Class);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Absorbed Radiation Not Vegetated Patches % -Other Surfaces - %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Rock
if (Crock > 0)  && (h_S > 0)
    [sur_alb,e_sur.R]=Albedo_OS_Properties('Crock',h_S);
    [Rabsb_rock,Rrfl_vis_rock,Rrfl_nir_rock]=ShortwaveFluxesOS(Rsw,SvF,sur_alb); %%% [W/m^2]
else
    Rabsb_rock=0; e_sur.R=1; Rrfl_vis_rock= 0; Rrfl_nir_rock=0;
end
%%%% Urban
if (Curb > 0)  && (h_S > 0)
    [sur_alb,e_sur.U]=Albedo_OS_Properties('Curb',h_S,Urb_Par);
    [Rabsb_urb,Rrfl_vis_urb,Rrfl_nir_urb]=ShortwaveFluxesOS(Rsw,SvF,sur_alb); %%% [W/m^2]
else
    Rabsb_urb=0;e_sur.U=1; Rrfl_vis_urb=0;Rrfl_nir_urb=0;
end
%%%%% Water
if ((Cwat > 0) || (ydepth >0)) && (h_S > 0)
    [sur_alb,e_sur.W]=Albedo_OS_Properties('Cwat',h_S);
    [Rabsb_wat,Rrfl_vis_wat,Rrfl_nir_wat]=ShortwaveFluxesOS(Rsw,SvF,sur_alb); %%% [W/m^2]
else
    Rabsb_wat=0;  e_sur.W =1; Rrfl_vis_wat=0; Rrfl_nir_wat=0;
    sur_alb = NaN;
end
%%% Glacier
if (Cice > 0) && (h_S > 0)
    if Cdeb == 0
        [ice_alb,e_sur.I]=Albedo_Ice_Properties(Aice);
        [Rabsb_ice,Rrfl_vis_ice,Rrfl_nir_ice]=ShortwaveFluxesOS(Rsw,SvF,ice_alb); %%% [W/m^2]
        Rabsb_deb=0; e_sur.D =1; Rrfl_vis_deb=0; Rrfl_nir_deb=0;
    else
        [sur_alb,e_sur.D]=Albedo_OS_Properties('Cdeb',h_S,Deb_Par);
        [Rabsb_deb,Rrfl_vis_deb,Rrfl_nir_deb]=ShortwaveFluxesOS(Rsw,SvF,sur_alb); %%% [W/m^2]
        Rabsb_ice=0; e_sur.I =1; Rrfl_vis_ice=0; Rrfl_nir_ice=0;
    end
else
    Rabsb_ice=0; e_sur.I=1; Rrfl_vis_ice=0; Rrfl_nir_ice=0;
    ice_alb = NaN;
    Rabsb_deb=0; e_sur.D =1; Rrfl_vis_deb=0; Rrfl_nir_deb=0;
end
%%%%%%%%%%%%%%%%%%%%
%%% Bare Soil
if (Cbare > 0) && (h_S > 0)
    [Rabsb_bare,Rrfl_vis_bare,Rrfl_nir_bare]=ShortwaveFluxesOS(Rsw,SvF,soil_alb); %%% [W/m^2]
else
    Rabsb_bare=0;Rrfl_vis_bare=0;Rrfl_nir_bare=0;
end
%%% Snow Cover
if (Csno > 0) && (h_S > 0)
    if sum(snow_alb.dir_vis + snow_alb.dif_vis+snow_alb.dir_nir + snow_alb.dif_nir)==0
        snow_alb=soil_alb;
    end
    [Rabsb_sno,Rrfl_vis_sno,Rrfl_nir_sno]=ShortwaveFluxesOS(Rsw,SvF,snow_alb); %%% [W/m^2]
else
    Rabsb_sno=0; Rrfl_vis_sno=0; Rrfl_nir_sno=0;
end
%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc=length(Ccrown);
RabsbSun_vegH=zeros(1,cc);
RabsbShd_vegH=zeros(1,cc);
Rabsb_soiH=zeros(1,cc);
Rrfl_vis_vegH=zeros(1,cc);Rrfl_nir_vegH=zeros(1,cc);PAR_sun_H=zeros(1,cc);PAR_shd_H=zeros(1,cc);
RabsbSun_vegL=zeros(1,cc);RabsbShd_vegL=zeros(1,cc);
Rabsb_soiL=zeros(1,cc);
Rrfl_vis_vegL=zeros(1,cc);Rrfl_nir_vegL=zeros(1,cc);PAR_sun_L=zeros(1,cc);PAR_shd_L=zeros(1,cc);
FsunH=zeros(1,cc);FshdH=zeros(1,cc);
FsunL=zeros(1,cc);FshdL=zeros(1,cc);
Kopt_H=zeros(1,cc);Kopt_L=zeros(1,cc);
fapar_H=zeros(1,cc);fapar_L=zeros(1,cc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Within Canopy - Clumping Factor %%%
LAI_H = LAI_H.*OM_H;
SAI_H = SAI_H.*OM_H;
LAId_H = LAId_H.*OM_H;
LAI_L = LAI_L.*OM_L;
SAI_L = SAI_L.*OM_L ;
LAId_L = LAId_L.*OM_L;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:cc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% SHORTWAVE %%%%%%%%%%%%%%%%%
    %%% Normal
    [RabsbSun_vegH(i),RabsbShd_vegH(i),Rabsb_soiH(i),Rrfl_vis_vegH(i),Rrfl_nir_vegH(i),PAR_sun_H(i),PAR_shd_H(i),...
        RabsbSun_vegL(i),RabsbShd_vegL(i),Rabsb_soiL(i),Rrfl_vis_vegL(i),Rrfl_nir_vegL(i),PAR_sun_L(i),PAR_shd_L(i),FsunH(i),FshdH(i),...
        FsunL(i),FshdL(i),Kopt_H(i),Kopt_L(i)]=ShortwaveFluxesVEG(Rsw,PAR,SvF,...
        LAI_H(i),SAI_H(i),LAId_H(i),LAI_L(i),SAI_L(i),LAId_L(i),PFT_opt_H(i),PFT_opt_L(i),snow_alb,soil_alb,h_S,Csno,dw_SNO,hc_L(i),SnoDep,Rabsb_sno);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Water Logging
    if ydepth > 0 && Csno == 0
        dw_WatH= min(1,ydepth./hc_H(i));
        dw_WatL= min(1,ydepth./hc_L(i));
        %%% 
        [RabsbSun_vegH(i),RabsbShd_vegH(i),Rabsb_soiH(i),Rrfl_vis_vegH(i),Rrfl_nir_vegH(i),PAR_sun_H(i),PAR_shd_H(i),...
            RabsbSun_vegL(i),RabsbShd_vegL(i),Rabsb_soiL(i),Rrfl_vis_vegL(i),Rrfl_nir_vegL(i),PAR_sun_L(i),PAR_shd_L(i),FsunH(i),FshdH(i),...
            FsunL(i),FshdL(i),Kopt_H(i),Kopt_L(i)]=ShortwaveFluxes_SUBM_VEG(Rsw,PAR,SvF,...
            LAI_H(i),SAI_H(i),LAId_H(i),LAI_L(i),SAI_L(i),LAId_L(i),PFT_opt_H(i),PFT_opt_L(i),sur_alb,h_S,dw_WatH,dw_WatL,hc_H(i),hc_L(i),ydepth);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    %%% Glacier
    if Cice == 1 && Csno == 0
        if Cdeb == 0
            if (LAI_H(i)+SAI_H(i)+LAId_H(i))>0 
                dw_ICE= min(1,IceDep./hc_H(i)); 
            else
                dw_ICE=1; 
            end 
            [RabsbSun_vegH(i),RabsbShd_vegH(i),Rabsb_soiH(i),Rrfl_vis_vegH(i),Rrfl_nir_vegH(i),PAR_sun_H(i),PAR_shd_H(i),...
                RabsbSun_vegL(i),RabsbShd_vegL(i),Rabsb_soiL(i),Rrfl_vis_vegL(i),Rrfl_nir_vegL(i),PAR_sun_L(i),PAR_shd_L(i),FsunH(i),FshdH(i),...
                FsunL(i),FshdL(i),Kopt_H(i),Kopt_L(i)]=ShortwaveFluxesVEG(Rsw,PAR,SvF,...
                LAI_H(i),SAI_H(i),LAId_H(i),LAI_L(i),SAI_L(i),LAId_L(i),PFT_opt_H(i),PFT_opt_L(i),ice_alb,soil_alb,h_S,Cice,dw_ICE,hc_L(i),IceDep,Rabsb_ice);
        else
            %%% Glacier with debris 
            [RabsbSun_vegH(i),RabsbShd_vegH(i),Rabsb_soiH(i),Rrfl_vis_vegH(i),Rrfl_nir_vegH(i),PAR_sun_H(i),PAR_shd_H(i),...
                RabsbSun_vegL(i),RabsbShd_vegL(i),Rabsb_soiL(i),Rrfl_vis_vegL(i),Rrfl_nir_vegL(i),PAR_sun_L(i),PAR_shd_L(i),FsunH(i),FshdH(i),...
                FsunL(i),FshdL(i),Kopt_H(i),Kopt_L(i)]=ShortwaveFluxesVEG(Rsw,PAR,SvF,...
                LAI_H(i),SAI_H(i),LAId_H(i),LAI_L(i),SAI_L(i),LAId_L(i),PFT_opt_H(i),PFT_opt_L(i),sur_alb,soil_alb,h_S,Cdeb,1,hc_L(i),IceDep,Rabsb_deb);
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    if (h_S > 0)
        fapar_H(i) = (PAR_sun_H(i)+PAR_shd_H(i))./(PAR.dir+PAR.dif);
        fapar_L(i) = (PAR_sun_L(i)+PAR_shd_L(i))./(PAR.dir+PAR.dif);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% COMPUTATION NDVI %%%%%%%%%%%%%%
if  (h_S > 0) && ((Rsw.dir_vis + Rsw.dif_vis+Rsw.dir_nir + Rsw.dif_nir)>0)
    if sum(LAI_H+SAI_H+LAId_H) > 0
        Rup_vis = (1-Csno)*(1-Cice)*(Cwat*Rrfl_vis_wat + Crock*Rrfl_vis_rock + Curb*Rrfl_vis_urb + Cbare*Rrfl_vis_bare) + ...
            (1-Csno)*Cice*((1-Cdeb)*Rrfl_vis_ice + Cdeb*Rrfl_vis_deb) + (1-sum(Ccrown))*Csno*Rrfl_vis_sno;
        Rup_nir = (1-Csno)*(1-Cice)*(Cwat*Rrfl_nir_wat + Crock*Rrfl_nir_rock + Curb*Rrfl_nir_urb + Cbare*Rrfl_nir_bare) + ...
            (1-Csno)*Cice*((1-Cdeb)*Rrfl_nir_ice + Cdeb*Rrfl_nir_deb) + (1-sum(Ccrown))*Csno*Rrfl_nir_sno;
    else
        Rup_vis = (1-Csno)*(1-Cice)*(Cwat*Rrfl_vis_wat + Crock*Rrfl_vis_rock + Curb*Rrfl_vis_urb + Cbare*Rrfl_vis_bare) + ...
            (1-Csno)*Cice*((1-Cdeb)*Rrfl_vis_ice + Cdeb*Rrfl_vis_deb) + Csno*Rrfl_vis_sno;
        Rup_nir = (1-Csno)*(1-Cice)*(Cwat*Rrfl_nir_wat + Crock*Rrfl_nir_rock + Curb*Rrfl_nir_urb + Cbare*Rrfl_nir_bare) + ...
            (1-Csno)*Cice*((1-Cdeb)*Rrfl_nir_ice + Cdeb*Rrfl_nir_deb) + Csno*Rrfl_nir_sno;
    end
    %%%
    if sum(LAI_H+SAI_H+LAId_H) > 0
        Rup_vis_veg = sum(Ccrown.*Rrfl_vis_vegH);
        Rup_nir_veg = sum(Ccrown.*Rrfl_nir_vegH);
        if sum(LAI_L+SAI_L+LAId_L) > 0  && (Csno == 0)
            Rup_vis_veg = Rup_vis_veg + sum(Ccrown.*Rrfl_vis_vegL);
            Rup_nir_veg = Rup_nir_veg + sum(Ccrown.*Rrfl_nir_vegL);
        end
    else
        if sum(LAI_L+SAI_L+LAId_L) > 0  && (Csno == 0)
            Rup_vis_veg = sum(Ccrown.*Rrfl_vis_vegL);
            Rup_nir_veg = sum(Ccrown.*Rrfl_nir_vegL);
        else
            Rup_vis_veg = 0;
            Rup_nir_veg = 0;
        end
    end
    rvis= (Rup_vis+ Rup_vis_veg)/(Rsw.dir_vis + Rsw.dif_vis);
    rnir= (Rup_nir+ Rup_nir_veg)/(Rsw.dir_nir + Rsw.dif_nir);
    NDVI=(rnir-rvis)/(rnir+rvis);
    ALB = rvis*(Rsw.dir_vis + Rsw.dif_vis)/(Rsw.dir_vis + Rsw.dif_vis+Rsw.dir_nir + Rsw.dif_nir) + rnir*(Rsw.dir_nir + Rsw.dif_nir)/(Rsw.dir_vis + Rsw.dif_vis+Rsw.dir_nir + Rsw.dif_nir);
else
    NDVI=0;
    ALB=0;
end
%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%
function[RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,Rrfl_vis_vegH,Rrfl_nir_vegH,PAR_sun_H,PAR_shd_H,...
    RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,Rrfl_vis_vegL,Rrfl_nir_vegL,PAR_sun_L,PAR_shd_L,FsunH,FshdH,...
    FsunL,FshdL,Kopt_H,Kopt_L]=ShortwaveFluxesVEG(Rsw,PAR,SvF,...
    LAI_H,SAI_H,LAId_H,LAI_L,SAI_L,LAId_L,PFT_opt_H,PFT_opt_L,snow_alb,soil_alb,h_S,Csno,dw_SNO,hc_L,SnoDep,Rabsb_sno)
%%%%%%%% INPUT
%%%%% INPUT
%%% Variables:
% Rsw  shortwave [W/m^2]
% PAR PAR radition [W/m^2]
% ShF shadow factor
% SvF sky view factor
%%% LAI_H SAI_H LAI_L SAI_L
%%% PFT_opt_H PFT_opt_L
%%% h_S
%%% Csno
%%% dw_SNO
%%% hc_L
%%% SnoDep
%%% snow_alb
%%% soil_alb
%Rsw.dir_vis= Rsw.dir_vis;
%Rsw.dir_nir= Rsw.dir_nir;
%%%
%Rsw.dif_vis = Rsw.dif_vis;
%Rsw.dif_nir = Rsw.dif_nir;
%%%
%PAR.dir = PAR.dir;
%PAR.dif = PAR.dif;
%%%%%%%% OUTPUT
%%%    Rabsb_vegH,Rabsb_soiH,NDVI_H,PAR_sun_H,PAR_shd_H,FsunH,FshdH
%%%    Rabsb_vegL,Rabsb_soiL,NDVI_L,PAR_sun_L,PAR_shd_L,FsunL,FshdL
%%%%%%%%%%%%%%%%%%%%%%%
if h_S <= 0 %%%%% NIGHT
    RabsbSun_vegH=0; RabsbShd_vegH=0;Rabsb_soiH=0;Rrfl_vis_vegH=0;Rrfl_nir_vegH=0;PAR_sun_H=0;PAR_shd_H=0;FsunH=0; FshdH=1; Kopt_H=Inf;
    RabsbSun_vegL=0; RabsbShd_vegL=0;Rabsb_soiL=0;Rrfl_vis_vegL=0;Rrfl_nir_vegL=0;PAR_sun_L=0;PAR_shd_L=0;FsunL=0; FshdL=1; Kopt_L=Inf;
else
    %%% DAYLIGHT-TIME
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%  Iup Idn alb Rsw --->  .dir_vis .dir_nir .dif_vis .dif_nir
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (LAI_H+SAI_H+LAId_H)>0 && (LAI_L+SAI_L+LAId_L)>0  %%%%%%% BOTH VEGETATION ARE PRESENT
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (SnoDep < hc_L)  %%%  Normal Condition and Height of snow -- Less than Vegetation L
            if Csno > 0
                [Iup_L,Idn_L,Kopt_L,om_vis_vg_L]=Canopy_Radiative_Transfer(PFT_opt_L,snow_alb,h_S,LAI_L,SAI_L,LAId_L,Csno);
            else
                [Iup_L,Idn_L,Kopt_L,om_vis_vg_L]=Canopy_Radiative_Transfer(PFT_opt_L,soil_alb,h_S,LAI_L,SAI_L,LAId_L,0);
            end
            [Iup_H,Idn_H,Kopt_H,om_vis_vg_H]=Canopy_Radiative_Transfer(PFT_opt_H,Iup_L,h_S,LAI_H,SAI_H,LAId_H,dw_SNO);
            %%%%%%%%%%%%%%%%%%
            PARL.dir=PAR.dir*(exp(-Kopt_H*(LAI_H+SAI_H+LAId_H)));
            PARL.dif=PAR.dif*Idn_H.dif_vis + PAR.dir*Idn_H.dir_vis;
            %%%%%%%%%%%%%%%%%%%%%%%%
            RswL.dir_vis =Rsw.dir_vis*(exp(-Kopt_H*(LAI_H+SAI_H+LAId_H)));
            RswL.dir_nir =Rsw.dir_nir*(exp(-Kopt_H*(LAI_H+SAI_H+LAId_H)));
            RswL.dif_vis =Rsw.dif_vis*Idn_H.dif_vis + Rsw.dir_vis*Idn_H.dir_vis;
            RswL.dif_nir= Rsw.dif_nir*Idn_H.dif_nir + Rsw.dir_nir*Idn_H.dir_nir;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,PAR_sun_H,PAR_shd_H,FsunH,FshdH,Rrfl_vis_vegH,Rrfl_nir_vegH]=ShortwaveFluxesVEG_unit(Rsw,...
                PAR,LAI_H,SAI_H,LAId_H,Kopt_H,om_vis_vg_H,Iup_L,Iup_H,Idn_H,SvF);
            if Csno > 0
                [RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,PAR_sun_L,PAR_shd_L,FsunL,FshdL,Rrfl_vis_vegL,Rrfl_nir_vegL]=ShortwaveFluxesVEG_unit(RswL,...
                    PARL,LAI_L,SAI_L,LAId_L,Kopt_L,om_vis_vg_L,snow_alb,Iup_L,Idn_L,1);
            else
                [RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,PAR_sun_L,PAR_shd_L,FsunL,FshdL,Rrfl_vis_vegL,Rrfl_nir_vegL]=ShortwaveFluxesVEG_unit(RswL,...
                    PARL,LAI_L,SAI_L,LAId_L,Kopt_L,om_vis_vg_L,soil_alb,Iup_L,Idn_L,1);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Rabsb_soiH = 0;
        else %%%%   Vegetation L cover by snow
            [Iup_H,Idn_H,Kopt_H,om_vis_vg_H]=Canopy_Radiative_Transfer(PFT_opt_H,snow_alb,h_S,LAI_H,SAI_H,LAId_H,dw_SNO);
            %%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,PAR_sun_H,PAR_shd_H,FsunH,FshdH,Rrfl_vis_vegH,Rrfl_nir_vegH]=ShortwaveFluxesVEG_unit(Rsw,...
                PAR,LAI_H,SAI_H,LAId_H,Kopt_H,om_vis_vg_H,snow_alb,Iup_H,Idn_H,SvF);
            %%%%%%%%%%%%%%%
            RabsbSun_vegL=0;
            RabsbShd_vegL=0;
            Rabsb_soiL=0;
            Rrfl_vis_vegL=0;Rrfl_nir_vegL=0; Kopt_L = 0;
            PAR_sun_L=0;PAR_shd_L=0;FsunL=0;FshdL=1;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        if (LAI_H+SAI_H+LAId_H)>0  %%% ONLY HIGH VEGETATION
            if Csno > 0 %% snow cover
                [Iup_H,Idn_H,Kopt_H,om_vis_vg_H]=Canopy_Radiative_Transfer(PFT_opt_H,snow_alb,h_S,LAI_H,SAI_H,LAId_H,dw_SNO);
                %%%
                [RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,PAR_sun_H,PAR_shd_H,FsunH,FshdH,Rrfl_vis_vegH,Rrfl_nir_vegH]=ShortwaveFluxesVEG_unit(Rsw,...
                    PAR,LAI_H,SAI_H,LAId_H,Kopt_H,om_vis_vg_H,snow_alb,Iup_H,Idn_H,SvF);
            else %%% snow free
                [Iup_H,Idn_H,Kopt_H,om_vis_vg_H]=Canopy_Radiative_Transfer(PFT_opt_H,soil_alb,h_S,LAI_H,SAI_H,LAId_H,dw_SNO);
                %%%
                [RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,PAR_sun_H,PAR_shd_H,FsunH,FshdH,Rrfl_vis_vegH,Rrfl_nir_vegH]=ShortwaveFluxesVEG_unit(Rsw,...
                    PAR,LAI_H,SAI_H,LAId_H,Kopt_H,om_vis_vg_H,soil_alb,Iup_H,Idn_H,SvF);
            end
            RabsbSun_vegL=0;
            RabsbShd_vegL=0;
            Rabsb_soiL=0;
            Rrfl_vis_vegL=0;Rrfl_nir_vegL=0; Kopt_L =0 ;
            PAR_sun_L=0;PAR_shd_L=0;FsunL=0;FshdL=0;
            %%%
        else
            if (LAI_L+SAI_L+LAId_L)>0   %%% ONLY LOW VEGETATION
                if (SnoDep < hc_L)  %%%  Height of snow -- Less than Vegetation L
                    if Csno > 0
                        [Iup_L,Idn_L,Kopt_L,om_vis_vg_L]=Canopy_Radiative_Transfer(PFT_opt_L,snow_alb,h_S,LAI_L,SAI_L,LAId_L,Csno);
                        %%%
                        [RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,PAR_sun_L,PAR_shd_L,FsunL,FshdL,Rrfl_vis_vegL,Rrfl_nir_vegL]=ShortwaveFluxesVEG_unit(Rsw,...
                            PAR,LAI_L,SAI_L,LAId_L,Kopt_L,om_vis_vg_L,snow_alb,Iup_L,Idn_L,SvF);
                    else  %%%%% snow free
                        [Iup_L,Idn_L,Kopt_L,om_vis_vg_L]=Canopy_Radiative_Transfer(PFT_opt_L,soil_alb,h_S,LAI_L,SAI_L,LAId_L,0);
                        %%%
                        [RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,PAR_sun_L,PAR_shd_L,FsunL,FshdL,Rrfl_vis_vegL,Rrfl_nir_vegL]=ShortwaveFluxesVEG_unit(Rsw,...
                            PAR,LAI_L,SAI_L,LAId_L,Kopt_L,om_vis_vg_L,soil_alb,Iup_L,Idn_L,SvF);
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else %%%%   Vegetation L covered by snow
                    Rabsb_soiL=Rabsb_sno;
                    %%%
                    RabsbSun_vegL=0;
                    RabsbShd_vegL=0;
                    Rrfl_vis_vegL=0;Rrfl_nir_vegL=0;Kopt_L = 0;
                    PAR_sun_L=0;PAR_shd_L=0;FsunL=0;FshdL=1;
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                RabsbSun_vegH=0;
                RabsbShd_vegH=0;
                Rabsb_soiH=0;
                Rrfl_vis_vegH = 0; Rrfl_nir_vegH =0; Kopt_H=0;
                PAR_sun_H=0;PAR_shd_H=0;FsunH=0;FshdH=0;
            else
                %%%  No Vegetation %%%%%%%%%%%%%%%%%%%%%%%%%%
                RabsbSun_vegL=0;
                RabsbShd_vegL=0;
                Rabsb_soiL=0;
                Rrfl_vis_vegL=0;Rrfl_nir_vegL=0; Kopt_L =0 ;
                PAR_sun_L=0;PAR_shd_L=0;FsunL=0;FshdL=0;
                RabsbSun_vegH=0;
                RabsbShd_vegH=0;
                Rabsb_soiH=0;
                Rrfl_vis_vegH =0; Rrfl_nir_vegH =0; Kopt_H=0;
                PAR_sun_H=0;PAR_shd_H=0;FsunH=0;FshdH=0;
            end
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%
function[RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,Rrfl_vis_vegH,Rrfl_nir_vegH,PAR_sun_H,PAR_shd_H,...
    RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,Rrfl_vis_vegL,Rrfl_nir_vegL,PAR_sun_L,PAR_shd_L,FsunH,FshdH,...
    FsunL,FshdL,Kopt_H,Kopt_L]=ShortwaveFluxes_SUBM_VEG(Rsw,PAR,SvF,...
    LAI_H,SAI_H,LAId_H,LAI_L,SAI_L,LAId_L,PFT_opt_H,PFT_opt_L,sur_alb,h_S,dw_WatH,dw_WatL,hc_H,hc_L,ydepth)
%%%%%%%%  Light attenuation in water -- // very simplified
%%%%
Kdwat=2.0;   %% 0.4-3.2 m-1  Light attenuation coefficient in water 
%%%%% INPUT
%%% Variables:
% Rsw  shortwave [W/m^2]
% PAR PAR radition [W/m^2]
% ShF shadow factor
% SvF sky view factor
%%% LAI_H SAI_H LAI_L SAI_L
%%% PFT_opt_H PFT_opt_L
%%% h_S
%%% Cwat
%%% dw_wat
%%% hc_L
%%% ydepth
%%% sur_alb
%%% soil_alb
%Rsw.dir_vis= Rsw.dir_vis;
%Rsw.dir_nir= Rsw.dir_nir;
%%%
%Rsw.dif_vis = Rsw.dif_vis;
%Rsw.dif_nir = Rsw.dif_nir;
%%%
%PAR.dir = PAR.dir;
%PAR.dif = PAR.dif;
%%%%%%%% OUTPUT
%%%    Rabsb_vegH,Rabsb_soiH,NDVI_H,PAR_sun_H,PAR_shd_H,FsunH,FshdH
%%%    Rabsb_vegL,Rabsb_soiL,NDVI_L,PAR_sun_L,PAR_shd_L,FsunL,FshdL
%%%%%%%%%%%%%%%%%%%%%%%
if h_S <= 0 %%%%% NIGHT
    RabsbSun_vegH=0; RabsbShd_vegH=0;Rabsb_soiH=0;Rrfl_vis_vegH=0;Rrfl_nir_vegH=0;PAR_sun_H=0;PAR_shd_H=0;FsunH=0; FshdH=1; Kopt_H=Inf;
    RabsbSun_vegL=0; RabsbShd_vegL=0;Rabsb_soiL=0;Rrfl_vis_vegL=0;Rrfl_nir_vegL=0;PAR_sun_L=0;PAR_shd_L=0;FsunL=0; FshdL=1; Kopt_L=Inf;
else
    %%% DAYLIGHT-TIME
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%  Iup Idn alb Rsw --->  .dir_vis .dir_nir .dif_vis .dif_nir
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (LAI_H+SAI_H+LAId_H)>0 && (LAI_L+SAI_L+LAId_L)>0  %%%%%%% BOTH VEGETATION ARE PRESENT
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if (ydepth < hc_L)  %%%  Height of water -- Less than Vegetation L
            [Iup_L,Idn_L,Kopt_L,om_vis_vg_L]=Canopy_Radiative_Transfer(PFT_opt_L,sur_alb,h_S,LAI_L,SAI_L,LAId_L,dw_WatL);
            %%%%% 
            [Iup_H,Idn_H,Kopt_H,om_vis_vg_H]=Canopy_Radiative_Transfer(PFT_opt_H,Iup_L,h_S,LAI_H,SAI_H,LAId_H,dw_WatH);
            %%%%%%%%%%%%%%%%%%
            PARL.dir=PAR.dir*(exp(-Kopt_H*(LAI_H+SAI_H+LAId_H)));
            PARL.dif=PAR.dif*Idn_H.dif_vis + PAR.dir*Idn_H.dir_vis;
            %%%%%%%%%%%%%%%%%%%%%%%%
            RswL.dir_vis =Rsw.dir_vis*(exp(-Kopt_H*(LAI_H+SAI_H+LAId_H)));
            RswL.dir_nir =Rsw.dir_nir*(exp(-Kopt_H*(LAI_H+SAI_H+LAId_H)));
            RswL.dif_vis =Rsw.dif_vis*Idn_H.dif_vis + Rsw.dir_vis*Idn_H.dir_vis;
            RswL.dif_nir= Rsw.dif_nir*Idn_H.dif_nir + Rsw.dir_nir*Idn_H.dir_nir;
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,PAR_sun_H,PAR_shd_H,FsunH,FshdH,Rrfl_vis_vegH,Rrfl_nir_vegH]=ShortwaveFluxesVEG_unit(Rsw,...
                PAR,LAI_H,SAI_H,LAId_H,Kopt_H,om_vis_vg_H,Iup_L,Iup_H,Idn_H,SvF);
            %%%
            [RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,PAR_sun_L,PAR_shd_L,FsunL,FshdL,Rrfl_vis_vegL,Rrfl_nir_vegL]=ShortwaveFluxesVEG_unit(RswL,...
                PARL,LAI_L,SAI_L,LAId_L,Kopt_L,om_vis_vg_L,sur_alb,Iup_L,Idn_L,1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Rabsb_soiH = 0;
        else %%%%   Vegetation L covered by water 
            [Iup_H,Idn_H,Kopt_H,om_vis_vg_H]=Canopy_Radiative_Transfer(PFT_opt_H,sur_alb,h_S,LAI_H,SAI_H,LAId_H,dw_WatH);
            %%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            [RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,PAR_sun_H,PAR_shd_H,FsunH,FshdH,Rrfl_vis_vegH,Rrfl_nir_vegH]=ShortwaveFluxesVEG_unit(Rsw,...
                PAR,LAI_H,SAI_H,LAId_H,Kopt_H,om_vis_vg_H,sur_alb,Iup_H,Idn_H,SvF);
            %%%%%%%%%%%%%%%
            %%%% Low vegetation covered by water 
            [Iup_L,Idn_L,Kopt_L,om_vis_vg_L]=Canopy_Radiative_Transfer(PFT_opt_L,sur_alb,h_S,LAI_L,SAI_L,LAId_L,dw_WatL);
            %%%
            Atten = exp(-Kdwat*(ydepth - hc_L)); 
            %%%
            PARL.dir=PAR.dir*(exp(-Kopt_H*(LAI_H+SAI_H+LAId_H)));
            PARL.dif=PAR.dif*Idn_H.dif_vis + PAR.dir*Idn_H.dir_vis;
            RswL.dir_vis =Rsw.dir_vis*(exp(-Kopt_H*(LAI_H+SAI_H+LAId_H)));
            RswL.dir_nir =Rsw.dir_nir*(exp(-Kopt_H*(LAI_H+SAI_H+LAId_H)));
            RswL.dif_vis =Rsw.dif_vis*Idn_H.dif_vis + Rsw.dir_vis*Idn_H.dir_vis;
            RswL.dif_nir= Rsw.dif_nir*Idn_H.dif_nir + Rsw.dir_nir*Idn_H.dir_nir;
            %%%
            PARL.dir=PARL.dir*Atten; 
            PARL.dif=PARL.dif*Atten; 
            RswL.dir_vis=RswL.dir_vis*Atten; 
            RswL.dir_nir=RswL.dir_nir*Atten; 
            RswL.dif_vis=RswL.dif_vis*Atten; 
            RswL.dif_nir=RswL.dif_nir*Atten; 
            %%%%% 
            [RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,PAR_sun_L,PAR_shd_L,FsunL,FshdL,Rrfl_vis_vegL,Rrfl_nir_vegL]=ShortwaveFluxesVEG_unit(RswL,...
                PARL,LAI_L,SAI_L,LAId_L,Kopt_L,om_vis_vg_L,sur_alb,Iup_L,Idn_L,1);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Rabsb_soiH = 0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        if (LAI_H+SAI_H+LAId_H)>0  %%% ONLY HIGH VEGETATION
            %%%%
            if (ydepth < hc_H)
                [Iup_H,Idn_H,Kopt_H,om_vis_vg_H]=Canopy_Radiative_Transfer(PFT_opt_H,sur_alb,h_S,LAI_H,SAI_H,LAId_H,dw_WatH);
                %%%
                [RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,PAR_sun_H,PAR_shd_H,FsunH,FshdH,Rrfl_vis_vegH,Rrfl_nir_vegH]=ShortwaveFluxesVEG_unit(Rsw,...
                    PAR,LAI_H,SAI_H,LAId_H,Kopt_H,om_vis_vg_H,sur_alb,Iup_H,Idn_H,SvF);
                %%%%%%%%%
            else %%%   Vegetation H covered by water
                [Iup_H,Idn_H,Kopt_H,om_vis_vg_H]=Canopy_Radiative_Transfer(PFT_opt_H,sur_alb,h_S,LAI_H,SAI_H,LAId_H,dw_WatH);
                %%%%
                Atten = exp(-Kdwat*(ydepth - hc_H));
                PAR.dir=PAR.dir*Atten;
                PAR.dif=PAR.dif*Atten;
                Rsw.dir_vis=Rsw.dir_vis*Atten;
                Rsw.dir_nir=Rsw.dir_nir*Atten;
                Rsw.dif_vis=Rsw.dif_vis*Atten;
                Rsw.dif_nir=Rsw.dif_nir*Atten;
                %%%
                [RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,PAR_sun_H,PAR_shd_H,FsunH,FshdH,Rrfl_vis_vegH,Rrfl_nir_vegH]=ShortwaveFluxesVEG_unit(Rsw,...
                    PAR,LAI_H,SAI_H,LAId_H,Kopt_H,om_vis_vg_H,sur_alb,Iup_H,Idn_H,SvF);
            end
            RabsbSun_vegL=0;
            RabsbShd_vegL=0;
            Rabsb_soiL=0;
            Rrfl_vis_vegL=0;Rrfl_nir_vegL=0; Kopt_L =0 ;
            PAR_sun_L=0;PAR_shd_L=0;FsunL=0;FshdL=0;
            %%%
        else
            if (LAI_L+SAI_L+LAId_L)>0   %%% ONLY LOW VEGETATION
                if (ydepth < hc_L)  %%%  Height of water -- Less than Vegetation L
                        [Iup_L,Idn_L,Kopt_L,om_vis_vg_L]=Canopy_Radiative_Transfer(PFT_opt_L,sur_alb,h_S,LAI_L,SAI_L,LAId_L,dw_WatL);
                        %%%
                        [RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,PAR_sun_L,PAR_shd_L,FsunL,FshdL,Rrfl_vis_vegL,Rrfl_nir_vegL]=ShortwaveFluxesVEG_unit(Rsw,...
                            PAR,LAI_L,SAI_L,LAId_L,Kopt_L,om_vis_vg_L,sur_alb,Iup_L,Idn_L,SvF);
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                else %%%%   Vegetation L covered by water
                    [Iup_L,Idn_L,Kopt_L,om_vis_vg_L]=Canopy_Radiative_Transfer(PFT_opt_L,sur_alb,h_S,LAI_L,SAI_L,LAId_L,dw_WatL);
                    %%%%
                    Atten = exp(-Kdwat*(ydepth - hc_L));
                    PAR.dir=PAR.dir*Atten;
                    PAR.dif=PAR.dif*Atten;
                    Rsw.dir_vis=Rsw.dir_vis*Atten;
                    Rsw.dir_nir=Rsw.dir_nir*Atten;
                    Rsw.dif_vis=Rsw.dif_vis*Atten;
                    Rsw.dif_nir=Rsw.dif_nir*Atten;
                    %%%
                    [RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,PAR_sun_L,PAR_shd_L,FsunL,FshdL,Rrfl_vis_vegL,Rrfl_nir_vegL]=ShortwaveFluxesVEG_unit(Rsw,...
                        PAR,LAI_L,SAI_L,LAId_L,Kopt_L,om_vis_vg_L,sur_alb,Iup_L,Idn_L,SvF);
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                end
                RabsbSun_vegH=0;
                RabsbShd_vegH=0;
                Rabsb_soiH=0;
                Rrfl_vis_vegH = 0; Rrfl_nir_vegH =0; Kopt_H=0;
                PAR_sun_H=0;PAR_shd_H=0;FsunH=0;FshdH=0;
            else
                %%%  No Vegetation %%%%%%%%%%%%%%%%%%%%%%%%%%
                RabsbSun_vegL=0;
                RabsbShd_vegL=0;
                Rabsb_soiL=0;
                Rrfl_vis_vegL=0;Rrfl_nir_vegL=0; Kopt_L =0 ;
                PAR_sun_L=0;PAR_shd_L=0;FsunL=0;FshdL=0;
                RabsbSun_vegH=0;
                RabsbShd_vegH=0;
                Rabsb_soiH=0;
                Rrfl_vis_vegH =0; Rrfl_nir_vegH =0; Kopt_H=0;
                PAR_sun_H=0;PAR_shd_H=0;FsunH=0;FshdH=0;
            end
        end
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function[RabsbSun_veg,RabsbShd_veg,Rabsb_soi,PAR_sun,PAR_shd,Fsun,Fshd,Rrfl_vis_veg,Rrfl_nir_veg]=ShortwaveFluxesVEG_unit(Rsw,...
    PAR,LAI,SAI,LAId,Kopt,om_vis_vg,sur_alb,Iup,Idn,SvF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% INPUT
%%% Variables:   Rsw.dir_vis, Rsw.dir_nir,
%%%              Rsw.dif_vis, Rsw.dif_nir
%%%              PAR.dir  PAR.dif
%%% LAI SAI Kopt om_vis_vg
%sur_alb.dir_vis =
%sur_alb.dif_vis =
%sur_alb.dir_nir =
%sur_alb.dif_nir =
soil_alb = sur_alb;
%%% Iup Idn
Iup.dir_vis = SvF*Iup.dir_vis;
Iup.dif_vis = SvF*Iup.dif_vis;
Iup.dir_nir = SvF*Iup.dir_nir;
Iup.dif_nir = SvF*Iup.dif_nir;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% OUTPUT
%%% SHORTWAVE RADIATION BALANCE (at the element scale):
%%%   -- Absorbed by vegetation   =  Rabsb_veg
%%%   -- Absorbed by soil         =  Rabsb_soi
%%%   -- Reflected by vegetation  =  Rrefl_veg
%%%   -- Reflected by soil        =  Rrefl_soi
%%% NDVI
%%% PAR_sun
%%% PAR_shd
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Kopt < 0
    disp('Optical depth is less than zero!')
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Kmpt = Kopt*sqrt(1-om_vis_vg); %%% Sellers, 1985
%Kmpt = Kopt; %%% Dai et al., 2004
%%% TO CHECK FOR LOW VEGETATION
%%%%%% Shaded and Sunlit Fraction of Canopy Fsun Fshd *something on
Fsun = (1.0 - exp(-Kmpt*(LAI+SAI+LAId)))/(Kmpt*(LAI+SAI+LAId));
Fsun(Fsun<0.01)=0; Fsun(Fsun>1)=1;
Fshd = 1- Fsun;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%// Compute absorption coefficients for the DIRECT BEAM fluxes
%%%// == VISIBLE range ==
Idir_vis = 1 - Iup.dir_vis - (1-soil_alb.dir_vis)*Idn.dir_vis -...
    (1-soil_alb.dir_vis)*exp(-Kopt*(LAI+SAI+LAId));
%%// == INFRARED range ==
Idir_nir = 1 - Iup.dir_nir - (1-soil_alb.dir_nir)*Idn.dir_nir -...
    (1-soil_alb.dir_nir)*exp(-Kopt*(LAI+SAI+LAId));
%%% // Compute absorption coefficients for the DIFFUSE fluxes
%%%	// == VISIBLE range ==
Idif_vis = 1 - Iup.dif_vis - (1-soil_alb.dif_vis)*Idn.dif_vis;
%%%	// == INFRARED range ==
Idif_nir = 1 - Iup.dif_nir - (1-soil_alb.dif_nir)*Idn.dif_nir;
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The total solar radiation ABSORBED by SUNLIT leaves
%%% NOTE: It is assumed that the sunlit leaves absorb the direct beam
%%%       radiation, all leaves absorb diffuse radiation, and that
%%%       leaves absorb LAI/(LAI+SAI) of the absorbed radiation
%%% == PAR range only ==
if (Fsun > 0) && (LAI > 0)
    PAR_sun = (PAR.dir*Idir_vis + Fsun*PAR.dif*Idif_vis);
    PAR_sun = PAR_sun*(LAI/(LAI+SAI+LAId));
else
    PAR_sun = 0;
end
%%% The total solar radiation ABSORBED by SHADED leaves
if (Fshd > 0) && (LAI > 0)
    if  (Fsun == 0)
        tmp1 = PAR.dir*Idir_vis;
        Fshd = 1.0;
    else
        tmp1 = 0;
    end
    PAR_shd = (tmp1 + Fshd*PAR.dif*Idif_vis); %
    PAR_shd = PAR_shd*(LAI/(LAI+SAI+LAId));
else
    PAR_shd = 0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% ------------------------------------------------------
%%% The total solar radiation ABSORBED by the VEGETATION
%%% Variables:   Rsw.dir_vis, Rsw.dir_nir,
%%%              Rsw.dif_vis, Rsw.dif_nir
%%%              PAR.dir  PAR.dif
%%% == VISIBLE and INFRARED range ==
RabsSun_vis_veg = Rsw.dir_vis*Idir_vis + Fsun*Rsw.dif_vis*Idif_vis;
RabsShd_vis_veg = Fshd*Rsw.dif_vis*Idif_vis;
RabsSun_nir_veg = Rsw.dir_nir*Idir_nir + Fsun*Rsw.dif_nir*Idif_nir;
RabsShd_nir_veg = Fshd*Rsw.dif_nir*Idif_nir;
RabsbSun_veg = RabsSun_vis_veg + RabsSun_nir_veg ;
RabsbShd_veg = RabsShd_vis_veg + RabsShd_nir_veg ;
%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The total solar radiation ABSORBED by the SOIL under the canopy
%%% == VISIBLE and INFRARED range ==
Rabs_vis_soi = Rsw.dir_vis*exp(-Kopt*(LAI+SAI+LAId))*(1-soil_alb.dir_vis) +...
    (Rsw.dir_vis*Idn.dir_vis + Rsw.dif_vis*Idn.dif_vis)*(1-soil_alb.dif_vis);
Rabs_nir_soi = Rsw.dir_nir*exp(-Kopt*(LAI+SAI+LAId))*(1-soil_alb.dir_nir) +...
    (Rsw.dir_nir*Idn.dir_nir + Rsw.dif_nir*Idn.dif_nir)*(1-soil_alb.dif_nir);
Rabsb_soi = Rabs_vis_soi +Rabs_nir_soi;
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% The solar radiation REFLECTED back to the ATMOSPHERE by VEGETATION
%%% == VISIBLE and INFRARED range ==
Rrfl_vis_veg = Rsw.dir_vis*Iup.dir_vis + Rsw.dif_vis*Iup.dif_vis;
Rrfl_nir_veg = Rsw.dir_nir*Iup.dir_nir + Rsw.dif_nir*Iup.dif_nir;
%%% Compute the Normalized Difference Vegetation Index (NDVI)
%tmp1 = (Rrfl_vis_veg)/(Rsw.dir_vis + Rsw.dif_vis);
%tmp2 = (Rrfl_nir_veg)/(Rsw.dir_nir + Rsw.dif_nir);
%NDVI = (tmp2 - tmp1)/(tmp2 + tmp1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end