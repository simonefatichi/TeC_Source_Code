%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  PLANT EXPORTS   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[]= Biogeo_Spin_Up()

jj=0;
for i=1:24:length(Ta)
    jj=jj+1;
    pdind = [max(1,i-24):i-1];
    [Se_bio(jj),Se_fc,Psi_bio(jj),Tdp_bio(jj),VSUM(jj),VTSUM(jj)]=Biogeo_environment(Tdp(pdind,:),O(pdind,:),V(pdind,:),Osat,Ohy,nVG,alpVG,Ks_Zs,L,Pe,O33,Ofc,SPAR,Bio_Zs);
    Aew(jj) = 0.000008575*exp(11.67*Se_bio(jj)./Se_fc)*(Se_bio(jj)>0.2); 
end
Aew(Aew>1)=1;  Se = (log(nanmean(Aew)/0.000008575)/11.67)*Se_fc;

%%%%%
ZBIOG = 0.20; %%[m]
rsd =1214; %% density dry soil [kg/m^3]
PH=6;  %% Soil - Ph
Zs=0.8; %% [m]
Se_fc =0.5994; 

Lat=   61.8475; Lon=  61.8475;
Ts= 5.8070;
Ta=  4.239 ;
%%%%
Psi_s =-0.0253; %%[MPa]
Se=   0.4501;
Lk=    0.5379; %%[mm/day]
V=153.2; %% [mm]
T = 0.771; %%[mm/day]
VT= 13.4750 ;%% [mm]
Broot = 395; %% [gC-Root/m2]
RexmyI =[0 0 0]; %% [gC m2 /day]
Ccrown=1;


%%%% INPUT  %% [gC/m^2 day]
%%%% Carbon Inputs
IS(1)= 0.254643151899803; %IS.C_met_sur_lit =  [gC/m^2 day]
IS(2)=   0.010614314347318;%IS.C_str_sur_lit_lig =
IS(3)= 0.069697565585479;%IS.C_str_sur_lit_nlig =
%%%%
IS(4)=  0.149738390292539;%IS.C_wod_sur_lit_lig =
IS(5)= 0.449215170877616;%IS.C_wod_sur_lit_nlig =
%%%%
IS(6)=  0.377784324905609;%IS.C_met_ssr_lit =
IS(7)=  0.057826996473149;%IS.C_str_ssr_lit_lig =
IS(8)=  0.239552991552177;%IS.C_str_ssr_lit_nlig =
%%%
IS(9)= 0.025026470782282;  %IS.C_myc
%%%% Nitrogen Inputs
IS(10)= 0.005645913077114; %IS.N_sur_lit =  %  [gN/m^2 day]
IS(11)=  0.002285480693939;%IS.N_wod_lit =
IS(12)=   0.008763794965459;%IS.N_ssr_lit =   % [gN/m^2 day]
%%%% Phosoporus Inputs
IS(13)=   3.136618376174205e-004;%IS.P_sur_lit =  % [gP/m^2 day]
IS(14)=  1.269711496632632e-004;%IS.P_wod_lit = ;
IS(15)=   4.868774980810796e-004;%IS.P_ssr_lit = ;   % [gP/m^2 day]
%%%% Potassium Inputs
IS(16)= 0.001881971025705;%IS.K_sur_lit =  ; % [gK/m^2 day]
IS(17)= 7.618268979795798e-004;%IS.K_wod_lit = ;
IS(18)=  0.002921264988486;%IS.K_ssr_lit = ;   % [gK/m^2 day]


%%%% INITIALIZATION  [gX /m^2]


%%%%%%% CARBON POOL %%%%%%%%%%%
Btm1(1)=90; %%% B1 Above-ground Litter Metabolic
Btm1(2)=160; %%% B2 Above-ground Litter Structural - Cellulose/Hemicellulose
Btm1(3)=25; %%% B3 Above-ground Litter Structura - Lignin
Btm1(4)=0; %%% B4 Above-ground Woody  - Cellulose/Hemicellulose
Btm1(5)=0; %%% B5 Above-ground Woody - Lignin
Btm1(6)=20; %%% B6 Below-ground Litter Metabolic
Btm1(7)=188; %%% B7 Below-ground Litter Structural - Cellulose/Hemicellulose
Btm1(8)=51; %%% B8 Below-ground Litter Structura - Lignin
Btm1(9)= 80; %%% B9  SOM-POC- lignin
Btm1(10) = 400; %%% B10 SOM-POC -Cellulose/Hemicellulose
Btm1(11) = 4000; %%% B11 SOM-MOC 
Btm1(12) = 20; %%%B12 DOC - for bacteria 
Btm1(13) = 40; %B13 DOC - for fungi 
Btm1(14) = 0.23; %%% B14 Enzyme for decomposition of POC-Bact 
Btm1(15) = 0.23; %%% B15 Enzyme for decomposition of POC-Fung
Btm1(16) = 0.23; %% B16 Enzyme for decomposition of MOC-Bact
Btm1(17) = 0.23; %%% B17 Enzyme for decomposition of MOC-Fung
Btm1(18) = 30; %%% B18 Bacteria pool
Btm1(19) = 90; %%%  B19 Fungi saprotrophic 
Btm1(20) = 10; %%% B20 AM-Mycorrhizal - C 
Btm1(21) = 10; %%%  B21 EM-Mycorrhizal - C 
Btm1(22) = 4; %%%  B22 Earthworms - C 
%%%
%%%%%% NITROGEN POOL
Btm1(23) = 3.7; %%% B23 Nitrogen Above-ground Litter
Btm1(24) = 0;%%% B24 Nitrogen Above-ground Woody
Btm1(25) = 1.8; %%% B25 Nitrogen Below-ground Litter
Btm1(26) = 650; %%% B26 Nitrogen SOM
Btm1(27) = 5.7;%%% B27 Nitrogen Bacteria 
Btm1(28) = 6;%%% B28 Nitrogen Fungi 
Btm1(29) = 6;%%% B29 AM Mycorrhizal - N 
Btm1(30) = 6;%%% B30 EM Mycorrhizal - N  
Btm1(31) = 0.87;%%% B31 Nitrogen Ione Ammonium NH4+
Btm1(32) = 0.039;%%%B32 Nitrogen Nitrate NO3-
Btm1(33) = 0.1; %%% B33 DON
Btm1(34) = 0.4; %%% B34 Earthworms - N 
%%%%%% PHOSPHORUS POOL
Btm1(35) = 0.21; %%% B35 phosphorus Above-ground Litter
Btm1(36) = 0; %%% B36 phosphorus Above-ground Woody
Btm1(37) = 0.097; %%% B37 phosphorus Below-ground Litter
Btm1(38) = 67; %%% B38 phosphorus SOM
Btm1(39) = 0.66;%%%%%% B39 phosphorus Bacteria 
Btm1(40) = 0.75;%%%% B40 phosphorus Fungi 
Btm1(41) = 0.75;%%% B41 AM - Mycorrhizal - P
Btm1(42) = 0.75;%%%% B42 EM - Mycorrhizal - P
Btm1(43) = 0.032;%%% B43 phosphorus Mineral
Btm1(44) = 250;%%%%% B44 phosphorus primary
Btm1(45) = 0.7;%% B45 phosphorus secondary
Btm1(46) = 15;%%%%% B46 phosphorus occluded
Btm1(47) = 0.1; %%%% B47 DOP
%%%%%% POTASSIUM POOL
Btm1(48) = 1.26;%%%% B48 Potassium Above-ground Litter
Btm1(49) = 0; %%% B36 Potassium  Above-ground Woody
Btm1(50) = 0.58;%%%% B50 Potassium  Below-ground Litter
Btm1(51) = 21; %%%% B51 Potassium SOM
Btm1(52) = 0.26; %%%% B52 Potassium  Mineral  solution
Btm1(53) = 0.26; %%%% B53 Potassium  exchangeable
Btm1(54) = 2.6; %%% B54 Potassium fixed or non-exchangeable 10- 80
Btm1(55) = 400; %%%% B55  Potassium in the lattice of certain primary minerals  130-9500
%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% External Inputs
DepN= 0.001578082191781; %
DepP= 9.202049856316554e-006; %  gP/m2 yr
DepK= 1.106849315068493e-004;
FertN=0*ones(1,366);
FertP=0*ones(1,366);
FertK=0*ones(1,366);
Upl=0.0; %% Soil production [mm/yr]
HIST=0; 
[B_IO]=Biogeochemistry_IO(Zs*1000,Lat,Lon,Upl,HIST,FertN,FertP,FertK,DepN,DepP,DepK);
%[B_IO]=Biogeochemistry_IO(Zs*1000,Lat,Lon,Upl);
ExEM = 0; 
B_IO.SC_par= [ 1 1 1 1]; % (Maint Bact , Maint Fungi,  EP, EW)
jDay=1;
FertN=B_IO.FertN(jDay);
DepN=B_IO.DepN;
FertP=B_IO.FertP(jDay);
DepP=B_IO.DepP;
FertK=B_IO.FertK(jDay);
DepK=B_IO.DepK;
Tup_P=B_IO.Tup_P;
Tup_K=B_IO.Tup_K;
SC_par=B_IO.SC_par;

opt_cons_CUE=0;
[BiogeoPar]=Biogeochemistry_Parameter(opt_cons_CUE);

t=1:3650;
B=zeros(length(t),55);
R_litter=zeros(length(t),1);
R_litter_sur=zeros(length(t),1);
R_microbe=zeros(length(t),1);
R_ew=zeros(length(t),1);
VOL=zeros(length(t),1);
BfixN = zeros(length(t),1);
Min_N = zeros(length(t),1);
Min_P = zeros(length(t),1);
R_bacteria= zeros(length(t),1);
RmycAM = zeros(length(t),1);
RmycEM = zeros(length(t),1);
N2flx = zeros(length(t),1);
NH4_Uptake=zeros(length(t),1);
NO3_Uptake = zeros(length(t),1);
P_Uptake=zeros(length(t),1);
K_Uptake = zeros(length(t),1);
LEAK_NH4 = zeros(length(t),1);
LEAK_NO3 = zeros(length(t),1);
LEAK_P = zeros(length(t),1);
LEAK_K = zeros(length(t),1);
LEAK_DOC = zeros(length(t),1);
LEAK_DON = zeros(length(t),1);
LEAK_DOP = zeros(length(t),1);

for j=1:length(t);
    %%%%%%%%%
    
    [LEAK_NH4(j),LEAK_NO3(j),LEAK_P(j),LEAK_K(j),LEAK_DOC(j),LEAK_DON(j),LEAK_DOP(j)]= Biogeo_Leakage(Btm1,Lk,V,BiogeoPar);
    %%%%
    [NH4_Uptake(j),NO3_Uptake(j),P_Uptake(j),K_Uptake(j)]= Biogeo_uptake(Btm1,Broot,Ts,T,VT,Ccrown,ExEM,BiogeoPar);
    %%%%
    [BfixN(j)]= Biogeo_Bio_fixation(0,0,Btm1,RexmyI,Ts);
    %
    %%%
    [dB,R_litter(j),R_microbe(j),R_litter_sur(j),R_ew(j),VOL(j),N2flx(j),Min_N(j),Min_P(j),R_bacteria(j),RmycAM(j),RmycEM(j)]= BIOGEOCHEMISTRY_DYNAMIC3(t(j),Btm1,ZBIOG,rsd,IS,Ts,Ta,Psi_s,PH,Se,Se_fc,FertN,DepN,BfixN(j),FertP,DepP,FertK,DepK,...
        NH4_Uptake(j),NO3_Uptake(j),P_Uptake(j),K_Uptake(j),LEAK_DOC(j),LEAK_NH4(j),LEAK_NO3(j),LEAK_P(j),LEAK_K(j),LEAK_DON(j),LEAK_DOP(j),Tup_P,Tup_K,ExEM,BiogeoPar,SC_par,opt_cons_CUE);
    
    %%%%
    if isreal(sum(dB))==0 || isnan(sum(dB)) == 1
        break
    end
    %%%
    B(j,:)=Btm1+dB;
    Btm1=B(j,:);
    %%%
end



%%%%%%%% SOIL BIOGEOCHEMISTRY BALANCE CHECK
IS=IS*length(t);
C_exp = sum(IS(1:9));
N_exp = sum(IS(10:12)) + DepN*length(t) + FertN*length(t);
P_exp = sum(IS(13:15)) + DepP*length(t) + Tup_P*length(t) + FertP*length(t) ;
K_exp = sum(IS(16:18)) + DepK*length(t) + Tup_K*length(t) + FertK*length(t) ;
dP_soil = B(1,:) - B(end,:);
C_out = sum(LEAK_DOC)+sum(R_litter)+sum(R_microbe)+sum(R_ew);
N_out = sum(sum(NH4_Uptake +NO3_Uptake)) + sum(LEAK_NH4) + sum(LEAK_NO3) + sum(LEAK_DON) + sum(VOL) +sum(N2flx);
P_out = sum(sum(P_Uptake)) + sum(LEAK_DOP) + sum(LEAK_P);
K_out = sum(sum(K_Uptake)) +sum(LEAK_K);


CkC_s = sum(dP_soil([1:22])) + C_exp - C_out;
CkN_s= sum(dP_soil([23:34]))+ N_exp - N_out;
CkP_s= sum(dP_soil([35:47]))+ P_exp - P_out;
CkK_s= sum(dP_soil([48:55]))+ K_exp - K_out;

CN_is = sum(IS(1:9))./sum(IS(10:12)); 
CP_is = sum(IS(1:9))./sum(IS(13:15));
CK_is = sum(IS(1:9))./sum(IS(16:18)); 


figure(106)
subplot(3,1,1)
plot(B(:,1),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('gC/m^2')
plot(B(:,2),'g','LineWidth', 1.5);
plot(B(:,3),'r','LineWidth', 1.5);
legend('Ab. Met','Ab Str','Ab Str Lig')
subplot(3,1,2)
plot(B(:,6),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('gC/m^2')
plot(B(:,7),'g','LineWidth', 1.5);
plot(B(:,8),'r','LineWidth', 1.5);
legend('Be. Met','Be Str','Be Str Lig')
subplot(3,1,3)
plot(B(:,4),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('gC/m^2')
plot(B(:,5),'g','LineWidth', 1.5);
legend('Ab. Wood','Ab Wood Lign')

figure(107)
subplot(2,2,1)
plot(B(:,9),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,10),'k','LineWidth', 1.5);
plot(B(:,11),'y','LineWidth', 1.5);
legend('SOM POC Lign','SOM POC - Cell','SOM MOC')
subplot(2,2,2)
plot(B(:,12),'k','LineWidth', 1.5);
hold on; grid on;
plot(B(:,13),'r','LineWidth', 1.5);
legend('DOC-B','DOC-F')
subplot(2,2,3)
plot(B(:,18),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,19),'k','LineWidth', 1.5);
plot(B(:,20),'r','LineWidth', 1.5);
plot(B(:,21),'b','LineWidth', 1.5);
plot(B(:,22),'y','LineWidth', 1.5);
title('Carbon Pool')
xlabel('Days'); ylabel('gC/m^2')
legend('Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza','Earthworms')
subplot(2,2,4)
plot(B(:,14),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,15),'k','LineWidth', 1.5);
plot(B(:,16),'r','LineWidth', 1.5);
plot(B(:,17),'b','LineWidth', 1.5);
title('Carbon Pool')
xlabel('Days'); ylabel('gC/m^2')
legend('EM-POC-B','EM-POC-F','EM-MOC-B','EM-MOC-F')

figure(108)
subplot(2,2,1)
plot(B(:,23),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('gN/m^2')
plot(B(:,24),'g','LineWidth', 1.5);
plot(B(:,25),'r','LineWidth', 1.5);
legend('Ab. Lit','Ab Wod','Be Lit')
subplot(2,2,2)
plot(B(:,26),'g','LineWidth', 1.5);
hold on; grid on;
legend('SOM')
subplot(2,2,3)
plot(B(:,27),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,28),'k','LineWidth', 1.5);
plot(B(:,29),'r','LineWidth', 1.5);
plot(B(:,30),'b','LineWidth', 1.5);
title('Nitrogen Pool')
xlabel('Days'); ylabel('gN/m^2')
legend('Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza')
subplot(2,2,4)
plot(B(:,31),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,32),'k','LineWidth', 1.5);
plot(B(:,33),'b','LineWidth', 1.5);
title('Nitrogen Pool')
xlabel('Days'); ylabel('gN/m^2')
legend('NH4+ ','NO3-','DON')

figure(109)
subplot(3,2,1)
plot(B(:,35),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('gP/m^2')
plot(B(:,36),'g','LineWidth', 1.5);
plot(B(:,37),'r','LineWidth', 1.5);
legend('Ab. Lit','Ab Wod','Be Lit')
subplot(3,2,2)
plot(B(:,38),'g','LineWidth', 1.5);
hold on; grid on;
legend('SOM')
subplot(3,2,3)
plot(B(:,39),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,40),'k','LineWidth', 1.5);
plot(B(:,41),'r','LineWidth', 1.5);
plot(B(:,42),'b','LineWidth', 1.5);
title('Phosporus Pool')
xlabel('Days'); ylabel('gP/m^2')
legend('Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza')
subplot(3,2,4)
plot(B(:,43),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,47),'r','LineWidth', 1.5);
title('Phosporus Pool')
xlabel('Days'); ylabel('gP/m^2')
legend('Mineral','DOP')
subplot(3,2,5)
plot(B(:,44),'k','LineWidth', 1.5);
hold on; grid on;
title('Phosporus Pool')
xlabel('Days'); ylabel('gP/m^2')
legend('Primary Material')
subplot(3,2,6)
plot(B(:,46),'k','LineWidth', 1.5);
hold on; grid on;
plot(B(:,45),'g','LineWidth', 1.5);
xlabel('Days'); ylabel('gP/m^2')
legend('Occluded','Secondary')


figure(110)
subplot(3,2,1)
plot(B(:,48),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('gK/m^2')
plot(B(:,49),'g','LineWidth', 1.5);
plot(B(:,50),'r','LineWidth', 1.5);
legend('Ab. Lit','Ab Wod','Be Lit')
subplot(3,2,2)
plot(B(:,51),'g','LineWidth', 1.5);
hold on; grid on;
legend('SOM')
subplot(3,2,3)
plot(B(:,52),'g','LineWidth', 1.5);
hold on; grid on;
title('Potassium Pool')
xlabel('Days'); ylabel('gK/m^2')
legend('Mineral Solution ')
subplot(3,2,4)
plot(B(:,53),'g','LineWidth', 1.5);
hold on; grid on;
plot(B(:,54),'k','LineWidth', 1.5);
title('Potassium Pool')
xlabel('Days'); ylabel('gK/m^2')
legend('Excheangeable ','Non-Excheangeable')
subplot(3,2,5)
plot(B(:,55),'g','LineWidth', 1.5);
hold on; grid on;
legend('Primary Minerals')



figure(111)
plot((B(:,2)+B(:,3)+B(:,1))./B(:,23),'g','LineWidth', 1.5);
hold on; grid on;
plot((B(:,5)+B(:,4))./B(:,24),'m','LineWidth', 1.5);
plot((B(:,6)+B(:,7)+B(:,8))./B(:,25),'k','LineWidth', 1.5);
plot((B(:,9)+B(:,10)+B(:,11))./B(:,26),'b','LineWidth', 1.5);
plot(B(:,18)./B(:,27),'r','LineWidth', 1.5);
plot(B(:,19)./B(:,28),'y','LineWidth', 1.5);
plot(B(:,20)./B(:,29),'c','LineWidth', 1.5);
plot(B(:,21)./B(:,30),'Color',[0.168 0.50586 0.3372],'LineWidth', 1.5);
plot(B(:,22)./B(:,34),'Color',[0.06 0.7 0.6],'LineWidth', 1.5);
title('C:N Ratio')
xlabel('Days'); ylabel('C:N')
legend('AG Litter','AG Wood','BG Litter','SOM','Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza','Earthworms')
%%%%%%%%%%%%%%%%%%%%%%%%%


figure(112)
plot((B(:,2)+B(:,3)+B(:,1))./B(:,35),'g','LineWidth', 1.5);
hold on; grid on;
plot((B(:,5)+B(:,4))./B(:,36),'m','LineWidth', 1.5);
plot((B(:,6)+B(:,7)+B(:,8))./B(:,37),'k','LineWidth', 1.5);
plot((B(:,9)+B(:,10)+B(:,11))./B(:,38),'b','LineWidth', 1.5);
plot(B(:,18)./B(:,39),'r','LineWidth', 1.5);
plot(B(:,19)./B(:,40),'y','LineWidth', 1.5);
plot(B(:,20)./B(:,41),'c','LineWidth', 1.5);
plot(B(:,21)./B(:,42),'Color',[0.168 0.50586 0.3372],'LineWidth', 1.5);
title('C:P Ratio')
xlabel('Days'); ylabel('C:P')
legend('AG Litter','AG Wood','BG Litter','SOM','Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza')
%%%%%%%%%%%%%%%%%%%%%%%%%


figure(113)
plot((B(:,2)+B(:,3)+B(:,1))./B(:,48),'g','LineWidth', 1.5);
hold on; grid on;
plot((B(:,5)+B(:,4))./B(:,50),'m','LineWidth', 1.5);
plot((B(:,6)+B(:,7)+B(:,8))./B(:,37),'k','LineWidth', 1.5);
plot((B(:,9)+B(:,10)+B(:,11))./B(:,51),'b','LineWidth', 1.5);
title('C:K Ratio')
xlabel('Days'); ylabel('C:K')
legend('AG Litter','AG Wood','BG Litter','SOM')
%%%%%%%%%%%%%%%%%%%%%%%%%

figure(114)
subplot(2,2,1)
plot(R_litter,'r','LineWidth', 1.5);
hold on; grid on;
plot(R_litter-R_litter_sur,'m','LineWidth', 1.5);
plot(R_microbe,'b','LineWidth', 1.5);
title('Respiration Het.')
plot(R_ew,'g','LineWidth', 1.5);
xlabel('Days'); ylabel('[gC/m2 day]')
legend('Litter','Litter below','Microbe','Earthworms')
subplot(2,2,2)
plot(VOL,'r','LineWidth', 1.5);
hold on; grid on;
plot(N2flx,'b','LineWidth', 1.5);
title('N - Fluxes')
xlabel('Days'); ylabel('[gN/m2 day]')
legend('NH_4 Vol.','N_2')
subplot(2,2,3)
plot(Min_N,'k','LineWidth', 1.5);
hold on ;  grid on 
title('N - Fluxes')
xlabel('Days'); ylabel('[gN/m2 day]')
legend('Min-N')
subplot(2,2,4)
plot(Min_P,'k','LineWidth', 1.5);
hold on ;  grid on 
title('P - Fluxes')
xlabel('Days'); ylabel('[gP/m2 day]')
legend('Min-P')
%%%%%%%%%%%%%%%%%%%%%%%%%


figure(115)
subplot(3,1,1)
plot(NH4_Uptake,'r','LineWidth', 1.5);
hold on; grid on;
plot(NO3_Uptake,'b','LineWidth', 1.5);
title('N Uptake.')
xlabel('Days'); ylabel('[gN/m2 day]')
legend('NH_4','NO_3')
subplot(3,1,2)
plot(P_Uptake,'g','LineWidth', 1.5);
hold on; grid on;
title('P Uptake')
xlabel('Days'); ylabel('[gP/m2 day]')
subplot(3,1,3)
plot(K_Uptake,'m','LineWidth', 1.5);
hold on; grid on;
title('K Uptake')
xlabel('Days'); ylabel('[gK/m2 day]')

%%%%%%%%%%%%%%%%%%%%%%%%%

figure(116)
subplot(2,2,1)
plot(LEAK_DOC,'b','LineWidth', 1.5);
hold on; grid on;
title('DOC Leaching')
xlabel('Days'); ylabel('[gC/m2 day]')
subplot(2,2,2)
plot(LEAK_NH4,'r','LineWidth', 1.5);
hold on; grid on;
plot(LEAK_NO3,'b','LineWidth', 1.5);
hold on; grid on;
plot(LEAK_DON,'m','LineWidth', 1.5);
title('N Leaching')
xlabel('Days'); ylabel('[gN/m2 day]')
legend('NH_4','NO_3','DON')
subplot(2,2,3)
plot(LEAK_P,'g','LineWidth', 1.5);
hold on; grid on;
plot(LEAK_DOP,'m','LineWidth', 1.5);
title('P Leaching')
legend('P','PON')
xlabel('Days'); ylabel('[gP/m2 day]')
subplot(2,2,4)
plot(LEAK_K,'m','LineWidth', 1.5);
hold on; grid on;
title('K Leaching')
xlabel('Days'); ylabel('[gK/m2 day]')
%%%%%%%%%%%%%%%%%%%%%%%%%


Lkday=Lk;
figure(120)
subplot(2,2,1)
plot(LEAK_DOC./(Lkday)*1000,'b','LineWidth', 1.5);
hold on; grid on;
title('DOC Conc.')
xlabel('Days'); ylabel('[mg/l]')
subplot(2,2,2)
plot(LEAK_NH4./(Lkday)*1000,'r','LineWidth', 1.5);
hold on; grid on;
plot(LEAK_NO3./(Lkday)*1000,'b','LineWidth', 1.5);
hold on; grid on;
plot(LEAK_DON./(Lkday)*1000,'m','LineWidth', 1.5);
title('N Conc.')
xlabel('Days'); ylabel('[mg/l]')
legend('NH_4','NO_3','DON')
subplot(2,2,3)
plot(LEAK_P./(Lkday)*1000*1000,'g','LineWidth', 1.5);
hold on; grid on;
plot(LEAK_DOP./(Lkday)*1000*1000,'m','LineWidth', 1.5);
title('P Conc.')
legend('P','DOP')
xlabel('Days');  ylabel('[ug/ l]')
subplot(2,2,4)
plot(LEAK_K./(Lkday)*1000,'m','LineWidth', 1.5);
hold on; grid on;
title('K Conc.')
xlabel('Days'); ylabel('[mg/ l]')


figure(121)
subplot(3,2,1)
plot((B(:,18)+B(:,19)+B(:,20) + B(:,21))./(B(:,9)+B(:,10)+B(:,11)+B(:,12)+B(:,13)),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[-]')
legend('Microbial/Substrate')
subplot(3,2,2)
plot((B(:,22))./(B(:,18)+B(:,19)+B(:,20) + B(:,21)),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[-]')
legend('Earthworms/Microbial')
subplot(3,2,3)
plot(sum(B(:,6:21),2)./(ZBIOG*rsd),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('SOC [gC /kg soil]')
subplot(3,2,4)
plot(sum(B(:,6:21),2),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('SOC [gC /m^2]')
subplot(3,2,5)
plot(((B(:,19) + B(:,20) + B(:,21))./B(:,18)),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[-]')
legend('Fungi/Bacteria')
subplot(3,2,6)
plot( B(:,19)./(B(:,20)+B(:,21)),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[-]')
legend('Saprotrophic/Mycorrhiza')

Rlitter_sub= R_litter-R_litter_sur; 
rb = R_bacteria./R_microbe; 
rf = 1 -rb; 
TSR = Rlitter_sub + R_microbe + R_ew; 
figure(122)
subplot(1,1,1)
plot((R_bacteria + rb.*Rlitter_sub)./TSR,'r','LineWidth', 1.5);
hold on; grid on;
plot((rf.*R_microbe + rf.*Rlitter_sub)./TSR,'b','LineWidth', 1.5);
title('Respiration Het.')
plot(R_ew./TSR,'g','LineWidth', 1.5);
xlabel('Days'); ylabel('[%]')
legend('Bacteria','Fungi','Earthworms')

return
