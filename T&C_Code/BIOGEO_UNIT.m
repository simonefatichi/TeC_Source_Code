%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[P,LEAK_NH4,LEAK_NO3,LEAK_P,LEAK_K,LEAK_DOC,LEAK_DON,LEAK_DOP,R_NH4,R_NO3,R_P,R_K,R_DOC,R_DON,R_DOP,Nuptake_H,Puptake_H,Kuptake_H,Nuptake_L,Puptake_L,Kuptake_L,RexmyI,...
    R_litter,R_microbe,R_litter_sur,R_ew,VOL,N2flx,Min_N,Min_P,R_bacteria,RmycAM,RmycEM,Prod_B,Prod_F,BfixN,NavlI,LitFirEmi]= BIOGEO_UNIT(Ptm1,IS,ZBIOG,rsd,PH,Ts,Ta,Psi_s,Se,Se_fc,V,VT,Ccrown,Bio_Zs,RfH_Zs,RfL_Zs,...
    Lk,Rd,Rh,Pr,T_H,T_L,Broot_H,Broot_L,LAI_H,LAI_L,...
    SupN_H,SupP_H,SupK_H,SupN_L,SupP_L,SupK_L,Rexmy,RexmyI,ExEM,NavlI,Pcla,Psan,...
    B_IO,jDay,FireA,AAET)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% INPUT
ZBIOG=ZBIOG/1000; %% [m]
%%%%%%%
%%% Mean value of an exponential distribution MV = Zperc/(-ln(1-perc)); %%
%%%%%
LAI=(LAI_H+LAI_L)*Ccrown'; %% [m2/m2]
%Broot =(Broot_H+Broot_L)*Ccrown'; % [gC/m2]
%%%%%%%%%%
ns = 365; 
RexmyI(1) = RexmyI(1)*(ns-1)/ns  + Rexmy(1)/ns; %%% Root exudation [gC/m2 day] 
RexmyI(2) = RexmyI(2)*(ns-1)/ns  + Rexmy(2)/ns; %% C export AM/EM [gC/m2 day] 
RexmyI(3) = RexmyI(3)*(ns-1)/ns  + Rexmy(3)/ns; %% C for Bfix [gC/m2 day] 
%%%%%%%%%%
cfTL=(sum((ones(length(Ccrown),1)*(Bio_Zs>0)).*RfL_Zs,2))';
cfTH=(sum((ones(length(Ccrown),1)*(Bio_Zs>0)).*RfH_Zs,2))';
T_H= (cfTH.*T_H); %% [mm/day] Transp from ZBIOG depth 
T_L= (cfTL.*T_L); %% [mm/day] Transp from ZBIOG depth 
%%%%
cc=length(Ccrown); 
NH4_Uptake_H=zeros(1,cc);NH4_Uptake_L=zeros(1,cc);
NO3_Uptake_H=zeros(1,cc);NO3_Uptake_L=zeros(1,cc);
P_Uptake_H=zeros(1,cc);P_Uptake_L=zeros(1,cc);
K_Uptake_H=zeros(1,cc);K_Uptake_L=zeros(1,cc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FertN=B_IO.FertN(jDay);
DepN=B_IO.DepN;
FertP=B_IO.FertP(jDay);
DepP=B_IO.DepP; 
FertK=B_IO.FertK(jDay); 
DepK=B_IO.DepK; 
Tup_P=B_IO.Tup_P; 
Tup_K=B_IO.Tup_K; 
SC_par = B_IO.SC_par; 
%%%%%%
if Rh>0
    f_run=1-Rh/Pr; f_run(f_run<0)=0; %% 1- fraction of direct infiltration excess runoff
    DepN=DepN*f_run;
    DepP=DepP*f_run;
    DepK=DepK*f_run;
end
%%%%%%%%%%%%%%%%%%%%%%%
ManF = B_IO.ManF(jDay);  %%% [gC /m2 day] 
N_Man =B_IO.N_Man; %% Manure  [gC/gN]
P_Man =B_IO.P_Man; % Manure  [gC/gP]
K_Man =B_IO.K_Man;% Manure  [gC/gK]
Lig_fr_Man =B_IO.Lig_fr_Man; %% Lignin fraction in Manure  [g Lignin / g DM] 
%%%%%%%%%%%%%
frac_to_metabolic_Man = 0.85 - 0.018*(N_Man*2*Lig_fr_Man); 
frac_to_metabolic_Man(frac_to_metabolic_Man<0)=0;
IMAN(1)= frac_to_metabolic_Man*(ManF) ; %% met_sur_lit  [gC/m^2 d]
IMAN(2)= (1-frac_to_metabolic_Man)*(ManF)*Lig_fr_Man ;%%  str_sur_lit_lig [gC/m^2 d]
IMAN(3)=  (1-frac_to_metabolic_Man)*(ManF)*(1-Lig_fr_Man);  %% str_sur_lit_nlig [gC/m^2 d]
IMAN(4)=  ManF./N_Man; % [gN/m^2 day]  
IMAN(5)=  ManF./P_Man; % [gP/m^2 day]
IMAN(6)=  ManF./K_Man; % [gK/m^2 day]
%%%%%%%%%%%%%%%%%%%%%%%%
opt_cons_CUE=1;
[BiogeoPar]=Biogeochemistry_Parameter(opt_cons_CUE);
%%%%%%%%%
[LEAK_NH4,LEAK_NO3,LEAK_P,LEAK_K,LEAK_DOC,LEAK_DON,LEAK_DOP]= Biogeo_Leakage(Ptm1,Lk+Rd,V,BiogeoPar);
%%%%
for cc=1:length(Ccrown)
    if Broot_H(cc) > 0 
    [NH4_Uptake_H(cc),NO3_Uptake_H(cc),P_Uptake_H(cc),K_Uptake_H(cc)]= Biogeo_uptake(Ptm1,Broot_H(cc),Ts,T_H(cc),VT,Ccrown(cc),ExEM,BiogeoPar);
    end 
    if Broot_L(cc) > 0 
    [NH4_Uptake_L(cc),NO3_Uptake_L(cc),P_Uptake_L(cc),K_Uptake_L(cc)]= Biogeo_uptake(Ptm1,Broot_L(cc),Ts,T_L(cc),VT,Ccrown(cc),ExEM,BiogeoPar);
    end 
end
%%%%
NH4_Uptake=sum(NH4_Uptake_H.*(1-SupN_H) + NH4_Uptake_L.*(1-SupN_L)); 
NO3_Uptake=sum(NO3_Uptake_H.*(1-SupN_H) + NO3_Uptake_L.*(1-SupN_L)); 
P_Uptake=sum(P_Uptake_H.*(1-SupP_H) + P_Uptake_L.*(1-SupP_L)); 
K_Uptake=sum(K_Uptake_H.*(1-SupK_H) + K_Uptake_L.*(1-SupK_L)); 
%%%%
[BfixN]= Biogeo_Bio_fixation(AAET,LAI,Ptm1,RexmyI,Ts);
%%%%%
t=[];
[dP,R_litter,R_microbe,R_litter_sur,R_ew,VOL,N2flx,Min_N,Min_P,R_bacteria,RmycAM,RmycEM,Prod_B,Prod_F]= BIOGEOCHEMISTRY_DYNAMIC3(t,Ptm1,ZBIOG,rsd,IS,Ts,Ta,Psi_s,PH,Se,Se_fc,FertN,DepN,BfixN,FertP,DepP,FertK,DepK,...
    NH4_Uptake,NO3_Uptake,P_Uptake,K_Uptake,LEAK_DOC,LEAK_NH4,LEAK_NO3,LEAK_P,LEAK_K,LEAK_DON,LEAK_DOP,Tup_P,Tup_K,ExEM,Pcla,Psan,BiogeoPar,SC_par,IMAN,opt_cons_CUE);
%%%%%%
if isreal(sum(dP))==0 || isnan(sum(dP)) == 1 
    disp('NaN in Biogeochemistry Pools')
    return
end
%%%
P=Ptm1+dP;
if (sum(P>(10^15)))>1 
    disp('Issue in the Biogeochemistry Pools')
    return
end
%%%%%%%%%%%%%%%%%%%%%%%%%
rRd=Rd./(Lk+Rd); rRd(isnan(rRd))=0; 
R_NH4=rRd*LEAK_NH4; LEAK_NH4=LEAK_NH4-R_NH4;
R_NO3=rRd*LEAK_NO3;LEAK_NO3=LEAK_NO3-R_NO3;
R_P=rRd*LEAK_P;LEAK_P=LEAK_P-R_P;
R_K=rRd*LEAK_K;LEAK_K=LEAK_K-R_K;
R_DOC=rRd*LEAK_DOC;LEAK_DOC=LEAK_DOC-R_DOC;
R_DON=rRd*LEAK_DON;LEAK_DON=LEAK_DON-R_DON;
R_DOP=rRd*LEAK_DOP; LEAK_DOP=LEAK_DOP-R_DOP;
if Rh>0
    R_NH4=R_NH4+DepN*(1-f_run);
    R_P=R_P+DepP*(1-f_run);
    R_K=R_K+DepK*(1-f_run);
end
%%%%%  Passing External Uptakes in units of [./m2 PFT]
Nuptake_H=(NH4_Uptake_H+NO3_Uptake_H).*(1-SupN_H)./Ccrown;
Puptake_H=P_Uptake_H.*(1-SupP_H)./Ccrown;
Kuptake_H=K_Uptake_H.*(1-SupK_H)./Ccrown;
%%%%
Nuptake_L=(NH4_Uptake_L+NO3_Uptake_L).*(1-SupN_L)./Ccrown;
Puptake_L=P_Uptake_L.*(1-SupP_L)./Ccrown;
Kuptake_L=K_Uptake_L.*(1-SupK_L)./Ccrown;
%%%%
Nuptake_H(Ccrown == 0) = 0;
Nuptake_L(Ccrown == 0) = 0;
Puptake_H(Ccrown == 0) = 0;
Puptake_L(Ccrown == 0) = 0;
Kuptake_H(Ccrown == 0) = 0;
Kuptake_L(Ccrown == 0) = 0;
%%%%% Updating Mineral Nutrient in the soil mean of last 365 days 
n4=365; 
NavlI(1) = NavlI(1)*(n4-1)/n4  + (P(31)+P(32))/n4; 
NavlI(2) = NavlI(2)*(n4-1)/n4  + P(43)/n4; 
NavlI(3) = NavlI(3)*(n4-1)/n4  + P(52)/n4; 
%%%%%%%%%%%%%%%%
if FireA == 1 
    [P,LitFirEmi]=Litter_Fire(P,FireA); 
else
    LitFirEmi=[0 0]; 
end 
return