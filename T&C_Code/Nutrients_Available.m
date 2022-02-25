%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction NUTRIENT AVAILABILITY  %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[FNC,e_relN,rNc,rPc,rKc,rMc,rNcR,Navail,Pavail,Kavail]= Nutrients_Available(B,FNCtm1,St,Nreserve,Preserve,Kreserve,OPT_SoilBiogeochemistry)  
%%%%%%%%%
%%%%-> All computations are for /m2 PFT not ground 
%%%% INPUT
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
%%%%%
FiS = St.FiS; %%[-] ; %%% Factor to increment nutrient reserve buffer 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Total Nutrient Content 
%TNIT =  B(1)/Nl + B(2)/Ns + B(3)/Nr + B(4)/Nc + B(5)/Nf  + B(6)./Nh  ;  %% Total Nitrogen Plant [gN/m^2]
%TPHO =  B(1)/Phol + B(2)/Phos + B(3)/Phor + B(4)/Phoc + B(5)/Phof  + B(6)./Phoh  ;  %% Total Phosporus Plant [gP/m^2]
%TPOT =  B(1)/Kpotl + B(2)/Kpots + B(3)/Kpotr + B(4)/Kpotc + B(5)/Kpotf  + B(6)./Kpoth  ;  %% Total Potassium Plant [gK/m^2]
%%%%%%%%% Non-Structural Nutrient Cotent (Leaves, Fine Root, Reserve, Fruits )
TNITns =  B(1)/Nl +  B(3)/Nr + B(4)/Nc + B(5)/Nf   ;  %% Total Nitrogen Plant [gN/m^2]
TPHOns =  B(1)/Phol +  B(3)/Phor + B(4)/Phoc + B(5)/Phof  ;  %% Total Phosporus Plant [gP/m^2]
TPOTns =  B(1)/Kpotl +  B(3)/Kpotr + B(4)/Kpotc  + B(5)/Kpotf ;  %% Total Potassium Plant [gK/m^2]
%%%% 
Nsto = FiS*B(4)*(1/Nl-1/Nc); %% Flexible storage Nitrogen Plant [gN/m^2]
Psto = FiS*B(4)*(1/Phol-1/Phoc); %%
Ksto = FiS*B(4)*(1/Kpotl-1/Kpotc); %%
%%%%%%%%% Relative Concentration of Nutrients 
if Nreserve < 0
    rNc =  1 + Nreserve./TNITns ; 
else
    if Nreserve > Nsto
        rNc = 1 + (Nreserve - Nsto)./TNITns ; 
    else %% Nreserve < Nsto
        rNc = 1;
    end
end
if Preserve < 0
    rPc =  1 + Preserve./TPHOns ; 
else
    if Preserve > Psto
        rPc = 1 + (Preserve - Psto)./TPHOns ; 
    else %% +
        rPc = 1;
    end
end
if Kreserve < 0
    rKc =  1 + Kreserve./TPOTns ; 
else
    if Kreserve > Ksto
        rKc = 1 + (Kreserve - Ksto)./TPOTns ; 
    else %% 
        rKc = 1;
    end
end
if sum(B) == 0
    rNc=1; rPc=1; rKc=1; 
end 
if OPT_SoilBiogeochemistry == 1 
    if rNc < 0.64 || rPc <0.64 || rKc <0.64 %% || rNc > 1.85 || rPc > 1.85 || rKc > 1.85
        disp('Error Nutrient Concentrations outside allowed bounds!!')
        return
    end
end 
rMc = min([rNc,rPc,rKc]); %% Overall nutrient limitation used by exudation; 
%%%%%%%%%%%%%%%%%%%%%%%%% Minimum Nutriennt concentration is 0.65 of Non-struct nutrients 
Navail= Nreserve + 0.35*TNITns ;  %% Available Nitrogen Plant [gN/m^2]
Pavail= Preserve + 0.35*TPHOns ;  %% Available Phosporus Plant [gP/m2] 
Kavail= Kreserve + 0.35*TPOTns;   %% Available Potassium Plant [gK/m2] 
Navail(Navail<0)=0; 
Pavail(Pavail<0)=0; 
Kavail(Kavail<0)=0; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
n4=90; %% 90 days 
FNC = FNCtm1*(n4-1)/n4 + (1-(rNc<0.66))/n4;  %% [-] Nitrogen Stress Factor
m= 0.5; % smoothing factor 
e_relN = m*(rNc-1)+1; %% Photosynthesis N efficiency 
rNcR = min(1.85,m*(rNc-1)+1); %%% Factor to increase respiration due to change in N-content 
if e_relN > 100 && OPT_SoilBiogeochemistry == 1
    disp('Error unrealistic Photosynthesis N efficiency!!') 
    return 
end 
%%%%%%%%%%%%%%%%%%%%%%%%
if OPT_SoilBiogeochemistry == 0 
    e_relN=1;  FNC = 1; rNcR = 1; 
    rNc=1; rPc=1; rKc=1; rMc = 1; 
end 
end 