%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction BIOGEOCHEMISTRY_DYNAMIC   %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% References Manzoni and Porporato, 2007; 2009;  Dickinson et al 2002;
%%% Xu-Ri and Prentice 2008; Li et al., 1992; 2000
%%% Porporato et al., 2003; Zahele and Friend 2010;
%%% Orwin et al 2011 ; Wang et al 2012
function[dB,R_litter,R_microbe,R_litter_sur,R_ew,VOL,N2flx,Min_N,Min_P,R_bacteria,RmycAM,RmycEM,Prod_B,Prod_F]= BIOGEOCHEMISTRY_DYNAMIC3(t,B,ZBIOG,rsd,...
    ISOIL,Ts,Ta,Psi_s,PH,Se,Se_fc,...
    FertN,DepN,BfixN,FertP,DepP,FertK,DepK,...
    NH4_Uptake,NO3_Uptake,P_Uptake,K_Uptake,LEAK_DOC,LEAK_NH4,LEAK_NO3,LEAK_P,LEAK_K,LEAK_DON,LEAK_DOP,Tup_P,Tup_K,ExEM,Pcla,Psan,BiogeoPar,SC_par,IMAN,opt_cons_CUE)
%%% OUTPUT
%VOL ;
%N2flx ;
%R_microbe
%R_litter
%%%% INPUT
dtd=1; %% [day]
Psi_s = Psi_s*1000; %% [kPa] Soil Water Potential
%Lk = [mm/d]
%V=[mm]
%%%%%%%%%%%
%ZBIOG = 0.25; %%[m]
%rsd =1300; %% density dry soil [kg/m^3]
rsd=rsd*1000; %% [gsoil /m3]
%%%% Carbon Inputs %% [gC/gsoil day]
IS.C_met_sur_lit =  ISOIL(1)/(ZBIOG*rsd) + IMAN(1)/(ZBIOG*rsd);
IS.C_str_sur_lit_lig = ISOIL(2)/(ZBIOG*rsd) + IMAN(2)/(ZBIOG*rsd);
IS.C_str_sur_lit_nlig = ISOIL(3)/(ZBIOG*rsd)+ IMAN(3)/(ZBIOG*rsd);
%%%%
IS.C_wod_sur_lit_lig =  ISOIL(4)/(ZBIOG*rsd);
IS.C_wod_sur_lit_nlig =  ISOIL(5)/(ZBIOG*rsd);
%%%%
IS.C_met_ssr_lit =  ISOIL(6)/(ZBIOG*rsd);
IS.C_str_ssr_lit_lig = ISOIL(7)/(ZBIOG*rsd);
IS.C_str_ssr_lit_nlig = ISOIL(8)/(ZBIOG*rsd);
%%%
IS.C_myc = ISOIL(9)/(ZBIOG*rsd);
%%%% Nitrogen Inputs %  [gN/gsoil day]
IS.N_sur_lit = ISOIL(10)/(ZBIOG*rsd)+ IMAN(4)/(ZBIOG*rsd);
IS.N_wod_lit = ISOIL(11)/(ZBIOG*rsd);
IS.N_ssr_lit = ISOIL(12)/(ZBIOG*rsd);
%%%% Phosoporus Inputs   [gP/gsoil day]
IS.P_sur_lit = ISOIL(13)/(ZBIOG*rsd) + IMAN(5)/(ZBIOG*rsd);
IS.P_wod_lit = ISOIL(14)/(ZBIOG*rsd);
IS.P_ssr_lit = ISOIL(15)/(ZBIOG*rsd);
%%%% Potassium Inputs   [gK/gsoil day]
IS.K_sur_lit =  ISOIL(16)/(ZBIOG*rsd) + IMAN(6)/(ZBIOG*rsd) ;
IS.K_wod_lit =  ISOIL(17)/(ZBIOG*rsd);
IS.K_ssr_lit = ISOIL(18)/(ZBIOG*rsd);
%%%%%%%%%%%%%%%%%
%%%% Tectonic Uplift %% [gX/gsoil d]
Tup_P= Tup_P/(ZBIOG*rsd);
Tup_K= Tup_K/(ZBIOG*rsd);
%%%%%%%% Plant and Mychorrizal uptake  %% [gX/gsoil d]
NH4_Uptake= NH4_Uptake/(ZBIOG*rsd);
NO3_Uptake=NO3_Uptake/(ZBIOG*rsd);
P_Uptake=P_Uptake/(ZBIOG*rsd);
K_Uptake=K_Uptake/(ZBIOG*rsd);
%%%%%%%%%%%%%% Leakage of Nutrient Compounds
LEAK_NH4 = LEAK_NH4/(ZBIOG*rsd) ; %% [gN/gsoil d]
LEAK_NO3 = LEAK_NO3/(ZBIOG*rsd) ; %% [gN/gsoil d]
LEAK_P = LEAK_P/(ZBIOG*rsd);   %% [gP/gsoil d]
LEAK_K = LEAK_K/(ZBIOG*rsd);  %% [gK/gsoil d]
LEAK_DOC = LEAK_DOC/(ZBIOG*rsd);
LEAK_DOP = LEAK_DOP/(ZBIOG*rsd);
LEAK_DON = LEAK_DON/(ZBIOG*rsd);
%%% METABOLIC (fast)  Proteins  Carbohydrate Reserve (Starch Fructants)
%%% Chlorophyll
%%% STRUCTURAL (Cellulose, hemicelllulose, pectin Lignin Tannin Polyphenolic
%%% Lipids Cutine Suberine
%%% OUTPUT
%%% dB [gX /m^2 d]
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% PARAMETERS
%%%%%%% CARBON POOL %%%%%%%%%%%
%%% B1 Above-ground Litter Metabolic
%%% B2 Above-ground Litter Structural - Cellulose/Hemicellulose
%%% B3 Above-ground Litter Structura - Lignin
%%% B4 Above-ground Woody  - Cellulose/Hemicellulose
%%% B5 Above-ground Woody - Lignin
%%% B6 Below-ground Litter Metabolic
%%% B7 Below-ground Litter Structural - Cellulose/Hemicellulose
%%% B8 Below-ground Litter Structura - Lignin
%%% B9  SOM-POC- lignin
%%% B10 SOM-POC -Cellulose/Hemicellulose
%%% B11 SOM-MOC
%%% B12 DOC - for bacteria
%%% B13 DOC - for fungi
%%% B14 Enzyme for decomposition of POC-Bact
%%% B15 Enzyme for decomposition of POC-Fung
%%% B16 Enzyme for decomposition of MOC-Bact
%%% B17 Enzyme for decomposition of MOC-Fung
%%% B18 Bacteria pool
%%% B19 Fungi saprotrophic
%%% B20 AM-Mycorrhizal - C
%%% B21 EM-Mycorrhizal - C
%%% B22 Earthworms - C
%%%%%% NITROGEN POOL
%%% B23 Nitrogen Above-ground Litter
%%% B24 Nitrogen Above-ground Woody
%%% B25 Nitrogen Below-ground Litter
%%% B26 Nitrogen SOM
%%% B27 Nitrogen Bacteria
%%% B28 Nitrogen Fungi
%%% B29 AM Mycorrhizal - N
%%% B30 EM Mycorrhizal - N
%%% B31 Nitrogen Ione Ammonium NH4+
%%% B32 Nitrogen Nitrate NO3-
%%% B33 DON
%%% B34 Earthworms - N
%%%%%% PHOSPHORUS POOL
%%% B35 phosphorus Above-ground Litter
%%% B36 phosphorus Above-ground Woody
%%% B37 phosphorus Below-ground Litter
%%% B38 phosphorus SOM
%%% B39 phosphorus Bacteria
%%% B40 phosphorus Fungi
%%% B41 AM - Mycorrhizal - P
%%% B42 EM - Mycorrhizal - P
%%% B43 phosphorus Mineral
%%% B44 phosphorus primary
%%% B45 phosphorus secondary
%%% B46 phosphorus occluded
%%% B47 DOP
%%%%%% POTASSIUM POOL
%%% B48 Potassium Above-ground Litter
%%% B49 Potassium  Above-ground Woody
%%% B50 Potassium  Below-ground Litter
%%% B51 Potassium SOM
%%% B52 Potassium  Mineral  solution
%%% B53 Potassium  exchangeable
%%% B54 Potassium fixed or non-exchangeable
%%% B55  Potassium in the lattice of certain primary minerals
%%%%%%%%
%%%%%%%%%%%%
B= B/(ZBIOG*rsd); %% [gX /g soil]
%%%%%%%% subtract quantities uptaken for internal computations
B(31)= B(31)- NH4_Uptake;%  Nitrogen Ione Ammonium NH4+
B(32)=B(32)-NO3_Uptake;%%   Nitrogen Nitrate NO3-
B(43)=B(43)-P_Uptake;%  phosphorus Mineral
B(52)=B(52)-K_Uptake;%% Potassium  Mineral  solution
B(B<=0)=1e-15;
%%%% Decay Constant on a 40°C basis (Kirschbaum 2002)
k_met_sur = BiogeoPar.k_met_sur; %%[1/day]
k_str_sur = BiogeoPar.k_str_sur; %% %[1/day]
k_met_ssr = BiogeoPar.k_met_ssr;
k_str_ssr = BiogeoPar.k_str_ssr;
k_wod_sur = BiogeoPar.k_wod_sur; %% [1/92 - 1/460]
%%%%%%%%%%%%%%%%%%%%
%%%%%% Respiration Coefficient for Litter [-]
%rr_met_sur=0.55;
%rr_str_sur=0.45;
%rr_wod_sur=0.55;
%rr_met_ssr=0.55;
%rr_str_ssr=0.55;
%%%%%%% Critical Ratios
R_CN_bac = BiogeoPar.R_CN_bac ; %%[C/N]  8.6
R_CN_fun = BiogeoPar.R_CN_fun ;
R_CN_mycAM = BiogeoPar.R_CN_myc;
R_CN_mycEM = BiogeoPar.R_CN_myc;
R_CN_ew = BiogeoPar.R_CN_ew;
R_CP_bac= BiogeoPar.R_CP_bac; %%[C/P] 60
R_CP_fun= BiogeoPar.R_CP_fun;
R_CP_mycAM = BiogeoPar.R_CP_myc;
R_CP_mycEM = BiogeoPar.R_CP_myc;
%%%%%%%%%%%%%%%%%%
%ExEM = [0-1] Fraction of EM vs AM
%%% Parameters
gd =BiogeoPar.gd; % [0.2-0.8] fraction of dead B allocated to D
fd =BiogeoPar.fd; % [0.2-0.8]  fraction of decomposed P allocated to D
%Ecb =0.47; % carbon use efficiency [-]
%Ecf =0.47;
Vdb = BiogeoPar.Vdb; %maximum specific uptake rate of D for growth of B [gC-D/gC-Bday]
Vdf = BiogeoPar.Vdf;
Kdb =BiogeoPar.Kdb ; % [0.00014 - 0.00038] half-saturation constant of uptake of D for growth of B [ g C/g soil]
Kdf =BiogeoPar.Kdf ;
mrb = SC_par(1)*BiogeoPar.mrb; % [ 1.20*1e-3 - 30*1e-3] specific maintenance factor or rate [ g C /g C day]
mrf = SC_par(2)*BiogeoPar.mrf;
mram = SC_par(2)*BiogeoPar.mram;
mrem = SC_par(2)*BiogeoPar.mrem;
Vpc =  BiogeoPar.Vpc; % maximum specific decomposition rate for P by EP [gC-P/gC-EP day]
Vpl =  BiogeoPar.Vpl;
Kpl =BiogeoPar.Kpl; % half-saturation constant for decomposition of P [gC / gsoil]
Kpc =BiogeoPar.Kpc ;
Vm = BiogeoPar.Vm;  %maximum specific decomposition rate for M by [gC-M /gC-EM day]
Km =BiogeoPar.Km; %half-saturation constant for decomposition of [ gC /gsoil];
rem = BiogeoPar.rem ; % [0.0075 - 0.075] turnover rate of EM [g C /g C day]
rep =BiogeoPar.rep ; %  [0.0075 - 0.075] turnover rate of EP [g C /g C day]
pepb = SC_par(3)*BiogeoPar.pepb; % [0.0031 - 0.031] fraction of mR for production of EP
pemb =SC_par(3)*BiogeoPar.pemb;% [0.0031 - 0.031] fraction of mR for production of EM
pepf = SC_par(3)*BiogeoPar.pepf; %
pemf =SC_par(3)*BiogeoPar.pemf;%
pepem = SC_par(3)*BiogeoPar.pepem; %
pemem =SC_par(3)*BiogeoPar.pemem;%
%%%%%%%%%% Leaching Coefficients
lambda_C = BiogeoPar.lambda_C;
lambda_N = BiogeoPar.lambda_N;
lambda_P = BiogeoPar.lambda_P;
lambda_K = BiogeoPar.lambda_K;
f_org_lea = BiogeoPar.f_org_lea; %% Fraction of Organic leaching
%%%% Earthworms parameter
%kse = BiogeoPar.kse ; % half-saturation constant for consumption of litter
VmaxE=  SC_par(4)*BiogeoPar.VmaxE ;%%% [gC / gC-food d] Curry and Schmidt 2007
rmanEa = BiogeoPar.rmanEa ; % %% [gC / gC-EW d]
rmanEp = BiogeoPar.rmanEp ; % %% [gC / gC-EW d]
dew = BiogeoPar.dew ; % %% [gC / gC-EW d]
fabsE = BiogeoPar.fabsE ; %% [-]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%  Soil Moisture regulation %%%%%%%%%%%%%
%%%% Moyano et al 2013; Manzoni et al 2012
%%% Psi_s in [kPa]
if Psi_s == 0
    fSM_litter =0 ; fSM_microbe = 0; Aew=0;
else
    Psi_s(Psi_s>-10)=-10.1;
    Psi_opt = -3;  Psi_th = -15800;  alpha = 1.47;
    fSM_microbe  = 1 - ((log(Psi_s)-log(Psi_opt))./(log(Psi_th)-log(Psi_opt))).^alpha;
    Psi_opt = -10;  Psi_th = -28800;  alpha = 1.0;
    fSM_litter =  1 - ((log(Psi_s)-log(Psi_opt))./(log(Psi_th)-log(Psi_opt))).^alpha;
    %%%%
    fSM_microbe(fSM_microbe <0)=0;
    fSM_litter(fSM_litter<0)=0;
    %%%%
    Aew = 0.000008575*exp(11.67*Se./Se_fc)*(Se>0.2);
    Aew(Se>0.9)=0.002729*Se(Se>0.9).^(-56.03); Aew(Aew>1)=1;
end
%if Se <= Sfc
%    fSM = Se/Sfc;
%else
%    fSM =  Sfc/Se;  %% [-] Porporato et al 2003
%end
%%% Wang et al 2012b
%PHopt=6; Psen=2;
%fPH= exp(-((PH-PHopt)/Psen)^2) ;
%%% Orwin et al 2011  [min 2 4.5 optm 7.5 12 max]
if PH<=4.5
    fPH = 0.4*PH - 0.8;
else
    if PH >=7.5
        fPH = -0.2222*PH + 2.666;
    else
        fPH=1;
    end
end
fPH(fPH>1)=1; fPH(fPH<0)=0;
%%%%%%%%% Clay Dependence of half-saturation and fraction of dead B (Wieder et al. 2015, GMD, Six et al 2002 PLSO)
fClay =   0.6940+1.36*(Pcla);
%fClay = (2.582*exp(-2*sqrt(Pcla))).^(-1);
fClay2 = (1.2*exp(-0.8*(Pcla)));
%%%%%%%%%%%%%%%%%%%%%% Temperature Effects Kinetiks
R= 8.314*1e-3 ; %[kJ mol-1 K-1]
%V=Vref*exp(-Ea/R*(1/(Ts+273.15) - 1/(Tref+273.15)));
Tref=BiogeoPar.T1;
%%%%%%%
Vpc = Vpc*exp(-37/R*(1/(Ts+273.15) - 1/(Tref+273.15)))*fSM_microbe*fPH;
Vpl = Vpl*exp(-53/R*(1/(Ts+273.15) - 1/(Tref+273.15)))*fSM_microbe*fPH;
Vm = Vm*exp(-47/R*(1/(Ts+273.15) - 1/(Tref+273.15)))*fSM_microbe*fPH;
Vdb = Vdb*exp(-47/R*(1/(Ts+273.15) - 1/(Tref+273.15)))*fSM_microbe*fPH;
Vdf = Vdf*exp(-47/R*(1/(Ts+273.15) - 1/(Tref+273.15)))*fSM_microbe*fPH;
Kpl = Kpl*exp(-30/R*(1/(Ts+273.15) - 1/(Tref+273.15)))*fClay;
Kpc = Kpc*exp(-30/R*(1/(Ts+273.15) - 1/(Tref+273.15)))*fClay;
Km = Km*exp(-30/R*(1/(Ts+273.15) - 1/(Tref+273.15)))*fClay;
Kdb= Kdb*exp(-30/R*(1/(Ts+273.15) - 1/(Tref+273.15)));
Kdf= Kdf*exp(-30/R*(1/(Ts+273.15) - 1/(Tref+273.15)));
mrb=mrb*exp(-20/R*(1/(Ts+273.15) - 1/(Tref+273.15)));
mrf=mrf*exp(-20/R*(1/(Ts+273.15) - 1/(Tref+273.15)));
mram=mram*exp(-20/R*(1/(Ts+273.15) - 1/(Tref+273.15)));
mrem=mrem*exp(-20/R*(1/(Ts+273.15) - 1/(Tref+273.15)));
%%%%
gd = gd*fClay2; gd(gd>1)=1; 
%%%% Control of available space for accumulating MOC 
fd = ((1-fd)*((B(11)*1000)./(4.825*(100-Psan*100).^0.6287)).^3.322+fd); fd(fd>1)=1; 
%%%%%% Soil moisture accessibility of DOC
%%% Unresolved issue if to use or not
fsm= 1;% (Se./Se_fc); fsm(fsm>1)=1; %%[-]
%%%%
%rem
%pemb pemf
%rep
%pepb pepf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Scaling of enzyme production with Ninorg or Norg/Ninorg 
%ROI=(B(:,26)+B(:,27)+B(:,28)+B(:,29)+B(:,30)+B(:,33))./(B(:,31)+B(:,32)); %%[-] 
knp=1;  %knp =  1-exp(-ROI/500); 
%%%% Scaling enzyme production with harsh conditions ;; Fundamental!!
%%%%%
%Cpep = 5.1285e-004*(B(18))^(-0.78); Cpep(Cpep>10)=10;
%Cpep = 4.2433e-004*(B(18))^(-0.85); Cpep(Cpep>10)=10;
%Cpep = 2.736e-004*(B(18))^(-0.85); Cpep(Cpep>6)=6;
%Cpep = 5.3409e-004*(B(18))^(-0.75); Cpep(Cpep>6)=6;
Cpep = 0.0078*(B(18))^(-0.5); Cpep(Cpep>6)=6;
%Cpep = 3.5e-004*(B(18))^(-0.85); Cpep(Cpep>5)=5;
pepb = knp*Cpep*pepb; %
pemb = knp*Cpep*pemb;%
%Cpep=0.0070*(B(19))^(-0.6); Cpep(Cpep>10)=10;
%Cpep = 9.0331e-004*(B(19))^(-0.85); Cpep(Cpep>10)=10;
%Cpep = 6.419e-004*(B(19))^(-0.85); Cpep(Cpep>3.5)=3.5;
%Cpep = 0.0016*(B(19))^(-0.75); Cpep(Cpep>3.5)=3.5;
Cpep = 0.0134*(B(19))^(-0.5); Cpep(Cpep>3.0)=3.0;
%Cpep = 8e-004*(B(19))^(-0.85); Cpep(Cpep>3)=3;
pepf = knp*Cpep*pepf; %
pemf = knp*Cpep*pemf;%
%%%% enzyme production from EM fungi 
pepem = knp*pepem; %
pemem = knp*pemem; %
%%%%%
%%%%%%
%pepb = rep*B(14)./(mrb*B(18)); pepb(pepb<0.0031)=0.0031; pepb(pepb>0.031)=0.031;
%pemb = rem*B(16)./(mrb*B(18)); pemb(pemb<0.0031)=0.0031; pemb(pemb>0.031)=0.031;
%pepf = rep*B(15)./(mrf*B(19)); pepf(pepf<0.0031)=0.0031; pepf(pepf>0.031)=0.031;
%pemf = rem*B(17)./(mrf*B(19)); pemf(pemf<0.0031)=0.0031; pemf(pemf>0.031)=0.031;
%%%%%%%% Litter Decoomposition on a 40°C basis
Tref4=BiogeoPar.T4;
fT2 = exp(3.36*(Ts-Tref4)/(Ts+31.79)); %% Kirschbaum 1995 2000 2002
fT1 = exp(3.36*(Ta-Tref4)/(Ta+31.79)); %% Kirschbaum 1995 2000 2002
fT2(Ts <-31)=0; fT1(Ta <-31)=0;
%fT = 2^((Ts-30)/10); %%% [-] %% Foley et al 1995 ;   Krinner et al., 2005

%%%%%%%%%%%%% Concentration in the soil
%%%%%%% C/N Ratios
%%% Parton et all 1988 - Krinner et al 2005  Orwin et al 2011
r_met_str_CN = 5; %% Ratio between metabolic and structural [C/N]
N_met_sur = B(23)/(1+(B(2)+B(3))/(r_met_str_CN*(B(1)))); %%
N_str_sur= B(23)-N_met_sur;
f_str_sur = N_str_sur./B(23);
N_met_ssr = B(25)/(1+(B(7)+B(8))/(r_met_str_CN*(B(6)))); %%[gN/gP]
N_str_ssr = B(25)-N_met_ssr;
f_str_ssr = N_str_ssr./B(25);
%%%%
CN_met_sur_lit= (B(1))/N_met_sur; %% [gC/gN]
CN_str_sur_lit= (B(2)+B(3))/N_str_sur; %% [gC/gN]
CN_wod_lit= (B(4)+B(5))/B(24); %% [gC/gN]
CN_met_ssr_lit= (B(6))/N_met_ssr; %% [gC/gN]
CN_str_ssr_lit= (B(7)+B(8))/N_str_ssr; %% [gC/gN]
CN_som = (B(9)+B(10)+B(11))/B(26); %% [gC/gN]
CN_bac = B(18)/B(27); %% [gC/gN]
CN_fun = B(19)/B(28); %% [gC/gN]
%CN_ew = B(22)/B(34); %% [gC/gN]
CN_mycAM = B(20)/B(29); %% [gC/gN]
CN_mycEM = B(21)/B(30); %% [gC/gN]
%%%  C/P Ratios
CP_met_sur_lit= (B(1))/((1-f_str_sur)*B(35)); %% [gC/gP]
CP_str_sur_lit= (B(2)+B(3))/(f_str_sur*B(35)); %% [gC/gP]
CP_wod_lit= (B(4)+B(5))/B(36); %% [gC/gP]
CP_met_ssr_lit= (B(6))/((1-f_str_ssr)*B(37)); %% [gC/gP]
CP_str_ssr_lit= (B(7)+B(8))/(f_str_ssr*B(37)); %% [gC/gP]
CP_som = (B(9)+B(10)+B(11))/(B(38)); %% [gC/gP]
CP_bac = B(18)/B(39); %% [gC/gP]
CP_fun = B(19)/B(40); %% [gC/gP]
CP_mycAM = B(20)/B(41); %% [gC/gP]
CP_mycEM = B(21)/B(42); %% [gC/gP]
%%% C/K Ratios
CK_met_sur_lit= (B(1))/((1-f_str_sur)*B(48)); %% [gC/gK]
CK_str_sur_lit= (B(2)+B(3))/(f_str_sur*B(48)); %% [gC/gK]
CK_wod_lit= (B(4)+B(5))/B(49); %% [gC/gK]
CK_met_ssr_lit= (B(6))/((1-f_str_ssr)*B(50)); %% [gC/gK]
CK_str_ssr_lit= (B(7)+B(8))/(f_str_ssr*B(50)); %% [gC/gK]
CK_som = (B(9)+B(10)+B(11))/B(51); %% [gC/gK]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% CUE-Litter and CUE-Microbe/Fungi
%%% Sinsabaguh et al 2013 ; Manzoni et al 2008
%rLO=1./L_cn; rb=1/B_cn;
%E= 0.45*rLO.^0.76./rb;
CUEmax= BiogeoPar.CUEmax;
%E = CUEmax./(1 +0.015*L_cn);
%TER_cn = 2.33*CN.^0.78;
%TER_cp = 2.91*CP.^0.83;
Ae=BiogeoPar.Ae;
%E = min(CUEmax,min(Ae*B_cb./TER_cn,Ae*B_cp./TER_cp));
%%%%%%%%%%%%%%%%%% SC_par(4)
CUE_met_sur=CUEmax./(1 +0.015*CN_met_sur_lit);%1-0.55;
CUE_str_sur=CUEmax./(1 +0.015*CN_str_sur_lit);%1-0.45;
CUE_wod_sur=CUEmax./(1 +0.015*CN_wod_lit);%1-0.55;
CUE_met_ssr=CUEmax./(1 +0.015*CN_met_ssr_lit);%1-0.55;
CUE_str_ssr=CUEmax./(1 +0.015*CN_str_ssr_lit);%1-0.55;
%%%%%
TER_cn = 2.33*CN_som.^0.78;
TER_cp = 2.91*CP_som.^0.83;
%%% Microbe
Ecb = BiogeoPar.Ecb; % carbon use efficiency [-]
%Ecb = min(CUEmax,min(Ae*CN_bac./TER_cn,Ae*CP_bac./TER_cp));
%%% Fungi
Ecf =BiogeoPar.Ecf;
%Ecf =  min(CUEmax,min(Ae*CN_fun./TER_cn,Ae*CP_fun./TER_cp));
%%%%%%%%%%%%%%%%%%%%
%%% Temperature effects on CUE
Tref=BiogeoPar.T2;
mtem = BiogeoPar.mtem; %%[°C] [0 - 0.016] Allison et al 2010; Li et al 2014
rr_met_sur = 1 - min(CUEmax,(CUE_met_sur + mtem*(Ts-Tref)));
rr_str_sur= 1 - min(CUEmax,(CUE_str_sur + mtem*(Ts-Tref)));
rr_wod_sur=1 - min(CUEmax,(CUE_wod_sur + mtem*(Ts-Tref)));
rr_met_ssr=1 - min(CUEmax,(CUE_met_ssr + mtem*(Ts-Tref)));
rr_str_ssr=1 - min(CUEmax,(CUE_str_ssr + mtem*(Ts-Tref)));
Ecb = min(CUEmax,Ecb + mtem*(Ts-Tref));
Ecf = min(CUEmax,Ecf + mtem*(Ts-Tref));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%opt_cons_CUE=0;
if opt_cons_CUE == 1
    Ecf = BiogeoPar.Ecf;
    Ecb = BiogeoPar.Ecb;
    rr_met_sur = BiogeoPar.rr_met_sur;
    rr_str_sur=   BiogeoPar.rr_str_sur;
    rr_wod_sur=  BiogeoPar.rr_wod_sur;
    rr_met_ssr= BiogeoPar.rr_met_ssr;
    rr_str_ssr=  BiogeoPar.rr_str_ssr;
end
%%%%%%%% Decomposition rates  %%[1/day]
D_met_sur_lit = k_met_sur*fT1;
D_str_sur_lit = k_str_sur*fT1*exp(-5*B(3)/(B(2)+B(3)));
D_met_ssr_lit = k_met_ssr*fT2*fSM_litter;
D_str_ssr_lit = k_str_ssr*fT2*fSM_litter*exp(-5*B(8)/(B(7)+B(8)));
D_wod_sur_lit = k_wod_sur*fT1*exp(-5*B(5)/(B(4)+B(5)));

%%% Decomposition Fluxes  [gC/gsoil  day]
D1=D_met_sur_lit*B(1);
D2=D_str_sur_lit*B(2);
D3=D_str_sur_lit*B(3);
D4=D_wod_sur_lit*B(4);
D5=D_wod_sur_lit*B(5);
D6=D_met_ssr_lit*B(6);
D7=D_str_ssr_lit*B(7);
D8=D_str_ssr_lit*B(8);
%%% Input to SOC
Ip_som_poc_lig = (1-rr_str_sur)*D3 + (1-rr_wod_sur)*D5 +  (1-rr_str_ssr)*D8 ; %%% [gC/gsoil d]
Ip_som_poc_cel = (1-rr_met_sur)*D1 + (1-rr_str_sur)*D2 +  (1-rr_met_ssr)*D6  + (1-rr_str_ssr)*D7  + (1-rr_wod_sur)*D4 ; %% [gC/gsoil d]
%%% leaching-C
Id = lambda_C*Ip_som_poc_cel;  %% [gC/gsoil d]
Ip_som_poc_cel = (1-lambda_C)*Ip_som_poc_cel;
%%%%
R_litter=  (rr_str_sur)*D3 + (rr_wod_sur)*D5 +  (rr_str_ssr)*D8 + (rr_met_sur)*D1 + (rr_str_sur)*D2 +  (rr_met_ssr)*D6  + (rr_str_ssr)*D7  + (rr_wod_sur)*D4; % [gC/gsoil d]
R_litter_sur =   (rr_str_sur)*D3 + (rr_wod_sur)*D5  + (rr_wod_sur)*D4 + (rr_met_sur)*D1 + (rr_str_sur)*D2 ; % [gC/gsoil d]
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Soil Carbon Fluxes  %% [gC/gsoil d]
F1b = (1/Ecb)*(Vdb+mrb)*(fsm*B(12))*B(18)/(Kdb+(fsm*B(12))); %% Flux from DOC to Bacteria
F1f = (1/Ecf)*(Vdf+mrf)*(fsm*B(13))*B(19)/(Kdf+(fsm*B(13))); %% Flux from DOC to Fungi
F2c_b =  Vpc*B(14)*B(10)/(Kpc + B(10)); %% Flux from POC-Cellulose to DOC - Bacteria
F2c_f =  Vpc*B(15)*B(10)/(Kpc + B(10)); %% Flux from POC-Cellulose to DOC - Fungi
F2l =  Vpl*B(15)*B(9)/(Kpl + B(9));  %% Flux from POC-Lignin to DOC  (Fungi only)
F3_b =   Vm*B(16)*B(11)/(Km + B(11)); %%% Flux from MOC to DOC-Bacteria
F3_f =   Vm*B(17)*B(11)/(Km + B(11)); %%% Flux from MOC to DOC-Fungi
F4b = (1/Ecb - 1)*Vdb*B(18)*(fsm*B(12))/(Kdb+(fsm*B(12)));  %% Respiration growth Bacteria
F4f = (1/Ecf - 1)*Vdf*B(19)*(fsm*B(13))/(Kdf+(fsm*B(13)));  %% Respiration growth  Fungi
F5b = (1/Ecb - 1)*mrb*B(18)*(fsm*B(12))/(Kdb+(fsm*B(12)));  %% Respiration Maintenance Bacteria
F5f = (1/Ecf - 1)*mrf*B(19)*(fsm*B(13))/(Kdf+(fsm*B(13)));  %% Respiration Maintenance  Fungi
F8b = (1-pepb -pemb)*mrb*B(18);  %% Microbe Bacteria  mortality
F8f = (1-pepf -pemf)*mrf*B(19); %% Fungi  mortality
F8am = mram*B(20);  %% AM  mortality
F8em = (1- pepem - pemem)*mrem*B(21);  %% EM  mortality
F9ep_b= pepb*mrb*B(18);  %% Enzyme Production
F9em_b= pemb*mrb*B(18);
F9ep_f =pepf*mrf*B(19);
F9em_f= pemf*mrf*B(19);
F9ep_em =pepem*mrem*B(21);
F9em_em= pemem*mrem*B(21);
F10ep_b= rep*B(14); %% Enzyme Turnover
F10ep_f= rep*B(15); %% Enzyme Turnover
F10em_b= rem*B(16);
F10em_f= rem*B(17);
%%%%%%%%%%%%%%%%%%%%%
Resp_mycAM = mram*B(20); %%% Respiration maintenance AM-Mycorrhizal
Resp_mycEM = mrem*B(21); %%% Respiration maintenance EM-Mycorrhizal
R_microbe =  F4b+F4f+F5f+F5b + Resp_mycAM + Resp_mycEM; % [gC/gsoil d]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Rbf = (B(12)/(B(12)+B(13))); %% Ratio of DOC - Bacteria to B+F
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Earthworms
fpal = 1.6 - (0.04*CN_som); fpal(fpal<0)=0;
Tref=BiogeoPar.T3; Tref_k=Tref+273.15;
Ts_k=Ts+273.15;
fT3 =exp(50*(Ts_k-Tref_k)./(Tref_k*R*Ts_k)).*(1+exp((Tref_k*0.649 - 200)/(Tref_k*R)))./(1+exp((Ts_k*0.649-200)./(Ts_k*R))); %% Assimilation
fT4 = exp(-20/R*(1/(Ts+273.15) - 1/(Tref+273.15))); %% Respiration
fClay3 =  2.938*exp(-4.82*Pcla) + -2.509*exp(-11.82*Pcla); fClay3(Pcla<0.15)=1;
%Few  = fabsE*B(22)*B(10)/(kse + B(10))*VmaxE*fpal*fT3*Aew; %%% [gC/ gsoil d] POC consumed by earthworms
Few  = fabsE*B(10)*VmaxE*fpal*fT3*fClay3*fPH*Aew; %%% % [gC/ gsoil d] POC consumed by earthworms
RespEw_m =  Aew*rmanEa*B(22)*fT4 + (1-Aew)*rmanEp*B(22)*fT4 ; %%% Maintenance Respiration
RespEw_g = (1-CUEmax)*Few ; %% Growth Respiration
%CUEapp  = CUEmax*((Few + RespEw_m)./Few);
RespE = RespEw_g + RespEw_m ; %% [gC/gsoil d]  Total EW Respiration
M_ew = dew*B(22)*fT4; %% % [gC/gsoil d]  Death of Earthworms
FewN =  (Few - RespE - M_ew)/R_CN_ew; %% [gN/gsoil d]  N-requirement
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Nitrogen Fluxes  %%% Kirschbaum et al 2002
if CN_bac < R_CN_bac
    %%% Minerelazation
    Nmin_b = B(27) - B(18)/R_CN_bac;  %% [gN/gsoil]
    Nimm_b=0;
else
    %%% Immobilization
    Nimm_b = B(18)/R_CN_bac - B(27);  %% [gN/gsoil]
    %if Nimm_b > 0.90*B(31)
    %    Nimm_b = 0.90*B(31);
    %end
    Nmin_b=0;
end
if CN_fun < R_CN_fun
    %%% Minerelazation
    Nmin_f = B(28) - B(19)/R_CN_fun;  %% [gN/gsoil]
    Nimm_f=0;
else
    %%% Immobilization
    Nimm_f = B(19)/R_CN_fun - B(28);  %% [gN/gsoil]
    Nmin_f=0;
    %if Nimm_f > 0.90*B(31)
    %    Nimm_f = 0.90*B(31);
    %end
end
if ExEM < 1
    if CN_mycAM < R_CN_mycAM
        %%% Minerelazation
        Nmin_am = B(29) - B(20)/R_CN_mycAM;  %% [gN/gsoil]
        Nimm_am=0;
    else
        %%% Immobilization
        Nimm_am = B(20)/R_CN_mycAM - B(29);  %% [gN/gsoil]
        Nmin_am=0;
        %if Nimm_am > 0.90*B(31)
        %    Nimm_am = 0.90*B(31);
        %end
    end
else
    Nimm_am=0; Nmin_am=0;
end
if ExEM > 0
    if CN_mycEM < R_CN_mycEM
        %%% Minerelazation
        Nmin_em = B(30) - B(21)/R_CN_mycEM;  %% [gN/gsoil]
        Nimm_em=0;
    else
        %%% Immobilization
        Nimm_em = B(21)/R_CN_mycEM - B(30);  %% [gN/gsoil]
        Nmin_em=0;
        %if Nimm_em > 0.90*B(31)
        %    Nimm_em = 0.90*B(31);
        %end
    end
else
    Nimm_em=0; Nmin_em=0;
end
%%% Input to SOC - Nitrogen
Ip_som_nit = (D1)/CN_met_sur_lit + (D2+D3)/CN_str_sur_lit  + (D4+D5)/CN_wod_lit  +  (D6)/CN_met_ssr_lit +(D7+D8)/CN_str_ssr_lit;
Ip_lea_nit = (lambda_N)*Ip_som_nit;  %% leaching of mineral-N
Ip_som_nit =  (1-lambda_N)*Ip_som_nit;
%%%%%%
N_imm_fix_b = 0.011*(B(31))/dtd; %% [gN/gsoil day] %% Continous immobilization due to chemical stabilization
N_imm_fix_f = 0.011*(B(32))/dtd; %% [gN/gsoil day] %% Continous immobilization due to chemical stabilization
%%%%
NMIN = (Nmin_b + Nmin_f + Nmin_am + Nmin_em)/dtd; %% [gN/gsoil day]
NIMM = (Nimm_b + Nimm_f + Nimm_am + Nimm_em  + N_imm_fix_f*dtd + N_imm_fix_b*dtd)/dtd; %% [gN/gsoil day]
NMIN_IN = NMIN - NIMM; %% [gN/gsoil day]
%%%%%%%%%%%%%
if  B(31) + NMIN_IN < 0
    Norg_imm = -NMIN_IN - B(31);
    %NIMM_M = B(31) + NMIN;
    B(31)=0;
    %Nimm_b = Nimm_b*NIMM_M/NIMM;
    %Nimm_f = Nimm_f*NIMM_M/NIMM;
    %Nimm_am = Nimm_am*NIMM_M/NIMM;
    %Nimm_em = Nimm_em*NIMM_M/NIMM;
    %N_imm_fix_f = N_imm_fix_f*NIMM_M/NIMM;
    %N_imm_fix_b = N_imm_fix_b*NIMM_M/NIMM;
else
    Norg_imm = 0;
end
%NIMM = (Nimm_b + Nimm_f + Nimm_am + Nimm_em  + N_imm_fix_f*dtd + N_imm_fix_b*dtd)/dtd; %% [gN/gsoil day]
%NMIN_IN = NMIN - NIMM; %% [gN/gsoil day]
%%%%%%% Porporato et al., 2003
NH4_in = 0.90*NMIN_IN; %%% [gN/gsoil d] Mineralization - Immoblization NH4
NO3_in = 0.10*NMIN_IN; %%% [gN/gsoil d]
%%%% External source, Biological fixation, deposition, fertilization
%%%  added as NH4 only to simplifiy -->  Dickinson et al., 2002
Incoming_N = (FertN + DepN + BfixN)/(ZBIOG*rsd);  %% [gN/ gsoil d]
%%%% Volatilization of Ammonia NH3
kv = BiogeoPar.kv; %% [1/d]
VOL = kv*B(31);  %% [gN/gsoil d]  Dickinson et al., 2002
%%%%%%% Nitrification (Aerobic)
KNmax = BiogeoPar.KNmaxn; %% [1/d] % Dickinson et al., 2002
ftemp=power((70-Ts)/(70-38),12)*exp(12*(Ts-38)/(70-38)); %%% Xu-Ri and Prentice 2008
fmoist = Se*(1-Se)/(0.25);  %% Dickinson et al., 2002
NO3flx = KNmax*ftemp*fmoist*B(31); %% [gN/gsoil d]
%%%%%%% Denitrification (Anaerobic)
%Kc = BiogeoPar.Kc/(ZBIOG*rsd); %% [gC/ gsoil] %% Li et al., 1992
%Kn = BiogeoPar.Kn/(ZBIOG*rsd); %% [gN/ gsoil] % Li et al., 1992
ftemp=exp(308.56*(1/68.02-1/(Ts+46.02))); % Xu-Ri and Prentice 2008
%fmoist=Se^(1/L); %% % Dickinson et al., 2002
fmoist=Se^(1/0.5); %% % Dickinson et al., 2002
%N2flx = (B(13)/(B(13) + Kc))*ftemp*fmoist*(B(28)/(B(28) + Kn));% [gN/gsoil d]
KNmax= BiogeoPar.KNmaxd; %% [1/d] % Dickinson et al., 2002
N2flx = KNmax*ftemp*fmoist*B(32); %% [gN/gsoil d]
%%%%%
S_DON = 0.01*B(33); %% Stabilization DON
%%%%%%%%%%%%%%%%%%%%%%%%
%Phosphorus Fluxes
if CP_bac < R_CP_bac
    %%% Minerelazation
    Pmin_b = B(39) - B(18)/R_CP_bac;  %% [gP/gsoil]
    Pimm_b=0;
else
    %%% Immobilization
    Pimm_b = B(18)/R_CP_bac - B(39);  %% [gP/gsoil]
    Pmin_b=0;
    %if Pimm_b > B(43)
    %    Pimm_b = B(43);
    %end
end
if CP_fun < R_CP_fun
    %%% Minerelazation
    Pmin_f = B(40) - B(19)/R_CP_fun;  %% [gP/gsoil]
    Pimm_f=0;
else
    %%% Immobilization
    Pimm_f = B(19)/R_CP_fun - B(40);  %% [g/gsoil]
    Pmin_f =0;
    %if Pimm_f > B(43)
    %    Pimm_f = B(43);
    %end
end
if ExEM < 1
    if CP_mycAM < R_CP_mycAM
        %%% Minerelazation
        Pmin_am = B(41) - B(20)/R_CP_mycAM;  %% [gP/gsoil]
        Pimm_am=0;
    else
        %%% Immobilization
        Pimm_am = B(20)/R_CP_mycAM - B(41);  %% [gP/gsoil]
        Pmin_am =0;
        %if Pimm_am > B(43)
        %    Pimm_am = B(43);
        %end
    end
else
    Pimm_am=0;Pmin_am =0;
end
if ExEM > 0
    if CP_mycEM < R_CP_mycEM
        %%% Minerelazation
        Pmin_em = B(42) - B(21)/R_CP_mycEM;  %% [gP/gsoil]
        Pimm_em=0;
    else
        %%% Immobilization
        Pimm_em = B(21)/R_CP_mycEM - B(42);  %% [g/gsoil]
        Pmin_em =0;
        %if Pimm_em > B(43)
        %    Pimm_em = B(43);
        %end
    end
else
    Pimm_em=0;Pmin_em =0;
end
%%% Input to SOC - Phosporus
Ip_som_pho = (D1)/CP_met_sur_lit + (D2+D3)/CP_str_sur_lit  + (D4+D5)/CP_wod_lit  +  (D6)/CP_met_ssr_lit +(D7+D8)/CP_str_ssr_lit;
Ip_lea_pho = (lambda_P)*Ip_som_pho;
Ip_som_pho =  (1-lambda_P)*Ip_som_pho;
%%%%%%
K1p= BiogeoPar.K1p; %%% [1/day]
K2p= BiogeoPar.K2p; %% [1/day]
K3p= BiogeoPar.K3p; %% [1/day]
K4p= BiogeoPar.K4p; %% [1/day]
PMIN = (Pmin_b + Pmin_f + Pmin_em + Pmin_am)/dtd ; PIMM =  (Pimm_b + Pimm_f + Pimm_em + Pimm_am)/dtd;  %% [gP/gsoil d]
PMIN_IN = (PMIN - PIMM );  %% [gP/gsoil d]
%%%
if  B(43) + PMIN_IN < 0
    Porg_imm = -PMIN_IN - B(43);
    %PIMM_M = B(43) + PMIN;
    B(43)=0;
    %Pimm_b = Pimm_b*PIMM_M/PIMM;
    %Pimm_f = Pimm_f*PIMM_M/PIMM ;
    %Pimm_am = Pimm_am*PIMM_M/PIMM ;
    %Pimm_em = Pimm_em*PIMM_M/PIMM ;
else
    Porg_imm = 0;
end
%PIMM =  Pimm_b + Pimm_f + Pimm_em + Pimm_am;
%PMIN_IN = (PMIN - PIMM);  %% [gP/gsoil d]
%%%%
Psec_exch = K1p*B(43) - K2p*B(45);
Ppri_exch = K4p*B(44); %% %% [gP/gsoil day]
Pocc_exch = K3p*B(45);
%%%
S_DOP = 0.01*B(47); %% Stabilization DOP
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%% Input to SOC - Potassium
Ip_som_pot = (D1)/CK_met_sur_lit + (D2+D3)/CK_str_sur_lit  + (D4+D5)/CK_wod_lit  +  (D6)/CK_met_ssr_lit +(D7+D8)/CK_str_ssr_lit;
Ip_lea_pot = (lambda_K)*Ip_som_pot;
Ip_som_pot =  (1-lambda_K)*Ip_som_pot;
%%%%%%
Kdk=BiogeoPar.Kdk;%%%% [1/day]
Kak= BiogeoPar.Kak; %% [1/day]
%%%%
K1k= BiogeoPar.K1k ; % []; %%% [1/day]
K2k= BiogeoPar.K2k; % []; %% [1/day]
K3k= BiogeoPar.K3k; %% [1/day]
Kex_sol =  Kdk*B(53)- Kak*B(52); %% Absrorptio desporption
Kfix_sol = K2k*B(54) - K1k*B(52);  %% From fixed to solution
Kmin_rel = K3k*B(55) ; %% %% [gK/gsoil day] Conversion mineral
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Incoming_P = (FertP + DepP)/(ZBIOG*rsd);  %% [gP/ gsoil d]
Incoming_K = (FertK + DepK)/(ZBIOG*rsd);  %% [gK/ gsoil d]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% MASS BALANCE %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%% Carbon [gC/g soil d]
dB(1)= IS.C_met_sur_lit -  D1; %%% B1  Above-ground Litter Metabolic
dB(2) = IS.C_str_sur_lit_nlig  - D2 ;  %%% B2  Above-ground Litter Structural - Cellulose/Hemicellulose
dB(3)= IS.C_str_sur_lit_lig - D3 ; %%% B3  Above-ground Litter Structura - Lignin
dB(4)= IS.C_wod_sur_lit_nlig - D4 ;  %%% B4  Above-ground Woody  - Cellulose/Hemicellulose
dB(5)= IS.C_wod_sur_lit_lig - D5 ;  %%% B5  Above-ground Woody - Lignin
dB(6)= IS.C_met_ssr_lit - D6; %%% B6 Below-ground Litter Metabolic
dB(7)= IS.C_str_ssr_lit_nlig - D7 ;  %%% B7 Below-ground Litter Structural - Cellulose/Hemicellulose
dB(8)= IS.C_str_ssr_lit_lig - D8 ; %%% B8 Below-ground Litter Structura - Lignin
%%%%
dB(9)=  Ip_som_poc_lig  -F2l ; %%% B9 SOM-POC- lignin
dB(10)= Ip_som_poc_cel - F2c_b - F2c_f  + (1-gd)*(F8b+F8f+F8am+F8em+M_ew) - Few ;  %%% B10 SOM-POC -  Cellulose/Hemicellulose
dB(11)= (1-fd)*(F2l + F2c_b + F2c_f) - F3_b - F3_f;  %%% B11 SOM-MOC
dB(12)= Rbf*Id + fd*(F2c_b) + Rbf*gd*(F8b+F8f+F8am+F8em+M_ew) + F3_b + (F10ep_b + F10em_b) - F1b - Rbf*LEAK_DOC ;  %%% B12 DOC - for bacteria
dB(13)= (1-Rbf)*Id + fd*(F2l + F2c_f) + (1-Rbf)*gd*(F8b+F8f+F8am+F8em+M_ew) + F3_f + (F10ep_f + F10em_f) - F1f  - (1-Rbf)*LEAK_DOC ;  %%% B13 DOC - Fungi
dB(14)= F9ep_b + Rbf*F9ep_em - F10ep_b;  %%% B14 Enzyme for decomposition of POC-Bact
dB(15)= F9ep_f + (1-Rbf)*F9ep_em  - F10ep_f;  %%% B15 Enzyme for decomposition of POC-Fung
dB(16)= F9em_b + Rbf*F9em_em - F10em_b; %%% B16 Enzyme for decomposition of MOC-Bact
dB(17)= F9em_f + (1-Rbf)*F9em_em  - F10em_f;  %%% B17 Enzyme for decomposition of MOC-Fung
dB(18)= F1b - (F4b+F5b) - F8b - (F9ep_b + F9em_b) ;  %%% B18 Bacteria pool
dB(19)= F1f - (F4f+F5f) - F8f - (F9ep_f + F9em_f) ;  %%% B19 Fungi saprotrophic
dB(20)= (1-ExEM)*IS.C_myc - F8am - Resp_mycAM; %%% B20 AM-Mycorrhizal
dB(21)= (ExEM)*IS.C_myc - F8em - Resp_mycEM - (F9ep_em + F9em_em)  ; %%% B21 EM-Mycorrhizal
dB(22)= Few - RespE - M_ew; %% B22 Earthworm -C
%%%%%%
%%%%%%%%%%%%%% Nitrogen [gN/ gsoil  d]
%%%%%% NITROGEN POOL
dB(23)=  IS.N_sur_lit - (D1)/CN_met_sur_lit - (D2+D3)/CN_str_sur_lit ; %%% B23 Nitrogen Above-ground Litter
dB(24)=  IS.N_wod_lit - (D4+D5)/CN_wod_lit;%%% B24 Nitrogen Above-ground Woody
dB(25)=  IS.N_ssr_lit - (D6)/CN_met_ssr_lit -(D7+D8)/CN_str_ssr_lit;%%% B25 Nitrogen Below-ground Litter
dB(26)=  Ip_som_nit - (fd*(F2l + F2c_b + F2c_f) + (F3_b + F3_f))/CN_som  + (F8b/CN_bac + F8f/CN_fun + F8am/CN_mycAM + F8em/CN_mycEM) -FewN + S_DON - Norg_imm; %%% B26 Nitrogen SOM
dB(27)=  (1-lambda_N)*(fd*(F2c_b) + (F3_b))/CN_som  - F8b/CN_bac + (N_imm_fix_b*dtd+Nimm_b-Nmin_b)/dtd;  %%% B27 Nitrogen Bacteria
dB(28)=  (1-lambda_N)*(fd*(F2l + F2c_f) + (F3_f))/CN_som  - F8f/CN_fun + (N_imm_fix_f*dtd+Nimm_f-Nmin_f)/dtd; %%% B28 Nitrogen Fungi S
dB(29)=  (Nimm_am-Nmin_am)/dtd - F8am/CN_mycAM ; %% B29 Nitrogen AM Mycorrhizal
dB(30)=  (Nimm_em-Nmin_em)/dtd - F8em/CN_mycEM ; %% B30 Nitrogen EM Mycorrhizal
dB(31)=  Incoming_N + (1-f_org_lea)*Ip_lea_nit  +  NH4_in  - NO3flx  - NH4_Uptake - LEAK_NH4 - VOL + 0.9*Norg_imm;  %%% B31 Nitrogen Ione Ammonium NH4+
dB(32)=  NO3_in + NO3flx - NO3_Uptake - LEAK_NO3 - N2flx + 0.1*Norg_imm ;%%% B32 Nitrogen Nitrate NO3-
dB(33) = f_org_lea*Ip_lea_nit + (lambda_N)*(fd*(F2l + F2c_b + F2c_f) + (F3_b + F3_f))/CN_som  - LEAK_DON - S_DON;  %%% B33 DON
dB(34) = FewN; %%  B34 Earthworm -N
%%%%%%
%%%%%%%%%%%%%% Phosphorus [gP/gsoil  d]
%%%%%% PHOSPHORUS POOL
dB(35)=  IS.P_sur_lit - (D1)/CP_met_sur_lit - (D2+D3)/CP_str_sur_lit ; %%%B35 phosphorus Above-ground Litter
dB(36)=  IS.P_wod_lit- (D4+D5)/CP_wod_lit;%%%B36 phosphorus Above-ground Woody
dB(37)=  IS.P_ssr_lit - (D6)/CP_met_ssr_lit -(D7+D8)/CP_str_ssr_lit;%%% B37 phosphorus Below-ground Litter
dB(38)=  Ip_som_pho - (fd*(F2l + F2c_b + F2c_f) + (F3_b + F3_f))/CP_som  + (F8b/CP_bac + F8f/CP_fun + F8am/CP_mycAM + F8em/CP_mycEM) + S_DOP - Porg_imm ; %%% B38 phosphorus SOM
dB(39)=  (1-lambda_P)*(fd*(F2c_b) + (F3_b))/CP_som - F8b/CP_bac + (Pimm_b - Pmin_b)/dtd  ;  %%%B39 phosphorus Bacteria
dB(40)=  (1-lambda_P)*(fd*(F2l + F2c_f) + (F3_f))/CP_som  - F8f/CP_fun + (Pimm_f - Pmin_f)/dtd; %%% B40 phosphorus Fungi
dB(41)= (Pimm_am - Pmin_am)/dtd - F8am/CP_mycAM ;%% B41 phosphorus AM- Mycorrhizal
dB(42)= (Pimm_em - Pmin_em)/dtd - F8em/CP_mycEM ;%% B42 phosphorus EM- Mycorrhizal
dB(43)=  Incoming_P + (1-f_org_lea)*Ip_lea_pho + PMIN_IN - P_Uptake- Psec_exch + Ppri_exch  - LEAK_P + Porg_imm ;  %%% B43 phosphorus Mineral
dB(44)=  Tup_P-Ppri_exch;%%%  B44 phosphorus primary
dB(45)= Psec_exch - Pocc_exch  ;  %%% B45 phosphorus secondary
dB(46)=  Pocc_exch  ;%%%% B46  phosphorus occluded
dB(47) = f_org_lea*Ip_lea_pho + (lambda_P)*(fd*(F2l + F2c_b + F2c_f) + (F3_b + F3_f))/CP_som   -LEAK_DOP - S_DOP;  %%% B47 DOP
%%%%%%
%%%%%%%%%%%%%% Potassium [gK/gsoil d]
%%%%%%%%% POTASSIUM POOL
dB(48)=  IS.K_sur_lit - (D1)/CK_met_sur_lit - (D2+D3)/CK_str_sur_lit ; %% B48 Potassium Above-ground Litter
dB(49)=  IS.K_wod_lit- (D4+D5)/CK_wod_lit;%%%B49 Potassium  Above-ground Woody
dB(50)=  IS.K_ssr_lit - (D6)/CK_met_ssr_lit -(D7+D8)/CK_str_ssr_lit;%%% B50 Potassium  Below-ground Litter
dB(51)=  Ip_som_pot -  (fd*(F2l + F2c_b + F2c_f) + (F3_b + F3_f))/CK_som ; %%%B51 Potassium SOM
dB(52)=  Incoming_K + Ip_lea_pot + (fd*(F2l + F2c_b + F2c_f) + (F3_b + F3_f))/CK_som  - K_Uptake  + Kmin_rel + Kfix_sol + Kex_sol -LEAK_K;  %%% B52 phosphorus Mineral solution
dB(53)=  -Kex_sol; %%%% B53 the exchangeable K
dB(54)=  -Kfix_sol ;  %%%% B54  the fixed or non-exchangeble K
dB(55)= Tup_K -Kmin_rel ;  %%%% B55 K in the lattice of certain primary minerals (Kl).
%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
dB= dB*(ZBIOG*rsd); %% [gX / m2 day]
R_litter_sur = R_litter_sur*(ZBIOG*rsd); %% %% [gC/ m2 day]
R_litter = R_litter*(ZBIOG*rsd); %% %% [gC/ m2 day]
R_microbe = R_microbe*(ZBIOG*rsd); %% %% [gC/ m2 day]
R_bacteria = (F4b+F5b)*(ZBIOG*rsd); %% %% [gC/ m2 day]
RmycAM = Resp_mycAM*(ZBIOG*rsd) ; %% [gC/ m2 day]
RmycEM = Resp_mycEM*(ZBIOG*rsd) ; %% [gC/ m2 day]
Prod_B=F1b*(ZBIOG*rsd);  %% [gC/ m2 day]
Prod_F=F1f*(ZBIOG*rsd);  %% [gC/ m2 day]
R_ew = RespE*(ZBIOG*rsd); %% %% [gC/ m2 day]
Min_N = NMIN_IN*(ZBIOG*rsd); %% %% [gN/ m2 day]
Min_P = PMIN_IN*(ZBIOG*rsd); %% %% [gP/ m2 day]
VOL=VOL*(ZBIOG*rsd);
N2flx=N2flx*(ZBIOG*rsd);
return
