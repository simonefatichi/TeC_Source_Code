%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  VEGETATION_DYNAMIC_RESULTS %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% References [Cox 2001] Bonan et al., 2003 Krinner et al.,  2005 Ivanov
%%% et al., 2008  Sitch et al 2003 White et al 2000  Knorr 2000  Arora and
%%% Boer 2003
function[LAI,NPP,Rg,RA,Rms,Rmr,Rmc,ANPP,LAIdead,Sr,Slf,Sfr,Swm,Sll,NLeaf,NLeafdead,NBLeaf,Nreserve,Preserve,Kreserve]= VEG_DYN_RES(B,dtd,Btm1,Tam,Tsm,An,Rdark,...
    Sl,mSl,Sl_emecrop,St,r,gR,aSE,AgeL,AgeDL,age_cr,dc_C,Tcold,Bfac,GF,dd_max,PHE_S,dsn,drn,fab,fbe,Wm,Mf,Klf,NBLeaftm1,dmg,Nreservetm1,Preservetm1,Kreservetm1,...
    Nuptake,Puptake,Kuptake,rNcR,rNc,rPc,rKc,OPT_EnvLimitGrowth)
%%%% INPUT
%Sl specific leaf area of  biomass [m^2 / kg C]
%Tam [°C]  Air daily temperature
%Tsm [°C]  Soil Daily temperature
% Tam [°C]  Mean Daily Temperature
% Tsm [°C]  Soil Daily Temperature
% An  [umolCO2/ s m^2]  Net Assimilation Rate
% Rdark  % [umolCO2/ s m^2]  Leaf Maintenace Respiration / Dark Respiration
% gR growth respiration  [] -- [Rg/(GPP-Rm)]
% r respiration rate at 10° [gC/gN d ]
% Ns [gC/gN] Sapwood
% Nr  [gC/gN] Fine root
%%%%%%% CARBON POOL %%%%%%%%%%%
%%% B1 Leaves - Grass
%%% B2 Sapwood
%%% B3 Fine Root
%%% B4 Carbohydrate Reserve
%%% B5 Fruit and Flower
%%% B6 Heartwood - Dead Sapwood
%%% B7 Leaves - Grass -- Standing Dead
%%%%%%%%%%%%%% Stochiometry  - Concentration -
Nl = St.Nl;  %%%  [gC/gN]
Ns = St.Ns;  %%%  [gC/gN]
Nr = St.Nr ; %%%  [gC/gN]
Nc= Ns; %%%  [gC/gN] Carbohydrate Reserve Carbon Nitrogen
Nf = St.Nf; %%%  [gC/gN]
Nh = St.Nh;  %%%  [gC/gN]
%%%
Phol = St.Phol;  %%%  [gC/gP]
Phos = St.Phos;  %%%  [gC/gP]
Phor = St.Phor ; %%%  [gC/gP]
Phoc= Phos; %%%  [gC/gP] Carbohydrate Reserve Carbon Phosophorus
Phof = St.Phof; %%%  [gC/gP]
Phoh = St.Phoh;  %%%  [gC/gP]
%%%
Kpotl = St.Kpotl;  %%%  [gC/gK]
Kpots = St.Kpots;  %%%  [gC/gK]
Kpotr = St.Kpotr ; %%%  [gC/gK]
Kpotc= Kpots; %%%  [gC/gK] Carbohydrate Reserve Carbon Potassium
Kpotf = St.Kpotf; %%%  [gC/gK]
Kpoth = St.Kpoth;  %%%  [gC/gK]
ftransL = St.ftransL;
ftransR = St.ftransR;
%%%%%%%%%%%%%%
if mSl == 0
    LAItm1 = Sl*(Btm1(1)); %% Leaf area index for green biomass

    if aSE==5  %% for crops
         %%% theoretically should be using the AgeLtm1 approximated as -1 day 
        Sl_n=Sl+(1-(AgeL-1)./dmg)*(Sl_emecrop)*((AgeL-1)<dmg); 
        LAItm1 = Sl_n*Btm1(1); %% Leaf area index for green biomass
    end

    LAIdeadtm1 = Sl*(Btm1(7)); %% Lead area index standing dead biomass
else
    LAItm1 = Sl*((exp(mSl*Btm1(1))-1)/mSl);
    LAIdeadtm1 = Sl*((exp(mSl*Btm1(7))-1)/mSl);
end
%%%%%%%%%%
if mSl == 0
    LAI = Sl*(B(1)); %% Leaf area index for green biomass

    if aSE==5  %% for crops
        Sl_n=Sl+(1-AgeL./dmg)*(Sl_emecrop)*(AgeL<dmg);
        LAI = Sl_n*B(1); %% Leaf area index for green biomass
    end

    LAIdead = Sl*(B(7)); %% Lead area index standing dead biomass
else
    LAI = Sl*((exp(mSl*B(1))-1)/mSl);
    LAIdead = Sl*((exp(mSl*B(7))-1)/mSl);
end
%%%%%%%%%%%%%% Maintenance and Growth Respiration
GPP = 1.0368*(An + Rdark); %% [gC / m^2 d]  Gross Primary Productivty  --> A
%%%%Ref--  Sitch 2003 Ivanov 2008 Ruimy 1996 Ryan 1991
gTam = exp(308.56*(1/56.02 - 1/(Tam+46.02)));
gTsm = exp(308.56*(1/56.02 - 1/(Tsm+46.02)));
Rms = fab*r*0.5*(B(2)+Btm1(2))*gTam/Ns +  fbe*r*0.5*(B(2)+Btm1(2))*gTsm/Ns; %%% [gC / m^2 d]
Rmr = rNcR*r*0.5*(B(3)+Btm1(3))*gTsm/Nr; %%% [gC / m^2 d]
if (aSE == 2)
    Rmc = rNcR*r*0.5*(B(4)+Btm1(4))*gTsm/Nc; %%% [gC / m^2 d]
else
    Rmc = fab*rNcR*r*0.5*(B(4)+Btm1(4))*gTam/Nc + fbe*rNcR*r*0.5*(B(4)+Btm1(4))*gTsm/Nc; %%% [gC / m^2 d]
end
Rm = Rms + Rmr + Rmc + 1.0368*Rdark; %%% [gC / m^2 d] Maintenance Respiration
%Rg = max(0,gR*GPP); %%%   Growth Respiration  [gC / m^2 d] only for GPP>0
Rg = max(0,gR*(GPP-Rm)); %%%   Growth Respiration  [gC / m^2 d] only for GPP>0
RA = Rg + Rm; %% Autotrphic Respiration  [gC / m^2 d]
NPP = GPP-RA; %% Net Primary Productivity [gC / m^2 d]  NPP = An - Rg -Rmr -Rms -Rmc
%%%%%% Turnover output --> dl dr ds Sll Ss Sr
%%%%%%%%%%%%%% Leaf Shed whit Age
switch aSE
    case 0
        dla= AgeL/((age_cr)^2); %% [1/d] Mortality for normal leaf age
    case 1
        dla= min(0.99,(1/age_cr)*(AgeL/age_cr).^4); %% [1/d] Mortality for normal leaf age
    case 2
        dla= min(1/age_cr,AgeL/((age_cr)^2)); %% [1/d] Mortality for normal grass age
    case 3
        if NPP > 0
            dlaK=(NBLeaftm1)/B(1)*age_cr; %% [-]
        else
            dlaK=1;
        end
        %%%%
        dla= dlaK*AgeL/((age_cr)^2); %% [1/d] Mortality for normal leaf age
        if AgeL>age_cr
            dla=max(0.25/age_cr,dla);
        end 
    case 5
        %dla= min(1/age_cr,AgeL/((age_cr)^2)); %% [1/d] Mortality as in grass
        dla=(1/age_cr)*(0.5*tanh(10*(AgeL/age_cr)-7)+0.5); %% [1/d] Mortality for crop
end
%%%%%%
%%% Leaf Mortality to Cold Stress  Linear [Cox 2001]
dc= (dc_C*(Tcold - Tam))*(Tam<Tcold); %% [1/d]
%%% Leaf Mortality to Drought  [Arora and Boer 2005]
dd = dd_max*(1-Bfac)^3*(Tsm>0.0);  %% [1/d]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dl = dla +dc +dd; %%% % [1/d] Leaf Mortality Senescence Cold Drought
if (aSE == 1) && (PHE_S == 1)
    dl = 1;  %% [1/d]
end
%%%%%%%%%%
dr = drn; %% [1/d] Fine root Turnover
ds = dsn; %%% [1/d]  Sapwood Conversion to Heartwood
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if aSE == 2
    fT = 2.0^(0.1*(Tam-20));
    Klf= Klf*fT;
end
dld= min(0.99,(Klf)*(AgeDL*Klf).^4); %% [1/d]
%%%
%Rexmy % Root Exudation and transfer to Mychorriza %% % [gC / m^2 d] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sll= dl*0.5*(B(1)+Btm1(1)); % Sll Mortality factor for leaves
Ss= ds*0.5*(B(2)+Btm1(2)); % Ss Converted factor for sapwood
Sr= dr*0.5*(B(3)+Btm1(3)); % Sr Mortality factor for fine roots
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Sfr = Mf*0.5*(B(5) + Btm1(5)); %% Sfr Fruit maturation % [gC/m^2 d]
Swm = Wm*0.5*(B(6) + Btm1(6)); %% Wood mortality % [gC/m^2 d]
Slf = dld*0.5*(B(7) + Btm1(7));  % Litterfall factor
%%%% Environemtnal Constraints Growth  GF  [gC /m2 day]
%if OPT_EnvLimitGrowth == 1
%    if NPP > GF ; %% [gC /m^2 d]
%        f_red = (GF/NPP); %%% [-]
%        NPP=f_red*NPP;
%    end
%end
%%%%%%%%%%%%%%%%% STOICHIOMETRIC Constraints Growth
if NPP > 0
    %%% B8 Idling Respiration
    WResp=B(8); %% [gC /m^2 d]
    Rmc= Rmc + WResp;
    RA = RA + WResp ; %% Autotrphic Respiration  [gC / m^2 d]
    NPP= NPP - WResp; %% % [gC /m^2 d]
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if (aSE == 2)
    ANPP= (B(1) -Btm1(1) + (B(5) -Btm1(5)))/dtd + Sll +Sfr; %%% [gC/m^2 d]
else
    ANPP= (B(1) -Btm1(1) + fab*(B(2) - Btm1(2)) + fab*(B(4) - Btm1(4)) + (B(5) -Btm1(5)) )/dtd + Sll + fab*Ss + Sfr; %%% [gC/m^2 d]
end
%%%%% New Biomass leaf 
NBLeaf =  (B(1)- Btm1(1)+Sll*dtd)*((B(1)- Btm1(1)+Sll*dtd)>0); %% New Biomass Leaf [gC/m2]  
%%%%% New Leaves
if mSl == 0
    Leafsen = Sl*Sll*dtd;
    
else
    Leafsen = (Sl +mSl*(LAI+LAItm1)/4)*Sll*dtd;
end
NLeaf= (LAI-LAItm1+Leafsen)*((LAI-LAItm1+Leafsen)>0); %% New Leaf [] LAI
%%% New dead Leaves
if mSl == 0
    Leaffall = Sl*Slf*dtd;
else
    Leaffall =  (Sl +mSl*(LAIdeadtm1+LAIdead)/4)*Slf*dtd;
end
NLeafdead= (LAIdead-LAIdeadtm1+Leaffall)*((LAIdead-LAIdeadtm1+Leaffall)>0); %% New dead Leaf [] LAI
%%%%% OTHER POOLS
%%%  Nitrogen - Reserve Pool
%%% Phosphoros - Reserve Pool
%%%  Potassium/Micronutrient Reserve Pool
%%%%% Exports  Nitrogen Plant [gN/m^2]
%(1-ftransL)*Slf./Nl + Sfr/Nf; % [gN/m^2 day]
%(1-ftransR)*Sr./Nr  ;  % [gN/m^2 day]
%fab*Swm/Nh;
%fbe*Swm/Nh;
%%%%% Exports  Phosporus Plant [gP/m^2]
%(1-ftransL)*Slf./Phol + Sfr/Phof; % [gP/m^2 day]
%(1-ftransR)*Sr./Phor  ;  % [gP/m^2 day]
%fab*Swm/Phoh;
%fbe*Swm/Phoh;
%%%%% Exports Potassium Plant [gK/m^2]
%(1-ftransL)*Slf./Kpotl + Sfr/Kpotf; % [gK/m^2 day]
%(1-ftransR)*Sr./Kpotr  ;  % [gK/m^2 day]
%fab*Swm/Kpoth;
%fbe*Swm/Kpoth;
%%%% Exports from plants (standing leaves not included)
TexN = rNc*(1-ftransL*(rNc<=1))*Sll./Nl +  rNc*Sfr/Nf +  rNc*(1-ftransR*(rNc<=1))*Sr./Nr  + Swm/Nh ; % [gN/m^2 day]
TexP = rPc*(1-ftransL*(rPc<=1))*Sll./Phol + rPc*Sfr/Phof + rPc*(1-ftransR*(rPc<=1))*Sr./Phor  + Swm/Phoh ;
TexK = rKc*(1-ftransL*(rKc<=1))*Sll./Kpotl + rKc*Sfr/Kpotf + rKc*(1-ftransR*(rKc<=1))*Sr./Kpotr  + Swm/Kpoth ;
%%%%%%%%
TNITtm1 =  Btm1(1)/Nl + Btm1(2)/Ns + Btm1(3)/Nr + Btm1(4)/Nc + Btm1(5)/Nf  + Btm1(6)./Nh  ;  %% Total Nitrogen Plant [gN/m^2 PFT]
TPHOtm1 =  Btm1(1)/Phol + Btm1(2)/Phos + Btm1(3)/Phor + Btm1(4)/Phoc + Btm1(5)/Phof  + Btm1(6)./Phoh  ;  %% Total Phosporus Plant [gP/m^2 PFT]
TPOTtm1 =  Btm1(1)/Kpotl + Btm1(2)/Kpots + Btm1(3)/Kpotr + Btm1(4)/Kpotc + Btm1(5)/Kpotf  + Btm1(6)./Kpoth  ;  %% Total Potassium Plant [gK/m^2 PFT]
%%%%%
TNIT =  B(1)/Nl + B(2)/Ns + B(3)/Nr + B(4)/Nc + B(5)/Nf  + B(6)./Nh  ;  %% Total Nitrogen Plant [gN/m^2]
TPHO =  B(1)/Phol + B(2)/Phos + B(3)/Phor + B(4)/Phoc + B(5)/Phof  + B(6)./Phoh  ;  %% Total Phosporus Plant [gP/m^2]
TPOT =  B(1)/Kpotl + B(2)/Kpots + B(3)/Kpotr + B(4)/Kpotc + B(5)/Kpotf  + B(6)./Kpoth  ;  %% Total Potassium Plant [gK/m^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Fictious Variable when equal to compute nutrient stress
Nreserve = Nreservetm1 + Nuptake*dtd + (TNITtm1 - TNIT) - TexN*dtd; %%% [gN/m^2]
Preserve = Preservetm1 + Puptake*dtd + (TPHOtm1 - TPHO) - TexP*dtd; %%% [gP/m^2]
Kreserve = Kreservetm1 + Kuptake*dtd + (TPOTtm1 - TPOT) - TexK*dtd; %%% [gK/m^2]
%%% 
return
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
