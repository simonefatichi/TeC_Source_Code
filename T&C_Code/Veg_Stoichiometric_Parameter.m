%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Stoichiometric_Parameter     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Stoich]=Veg_Stoichiometric_Parameter(Nl)
%%%%%%
%Nl = [15 -42 ]; %[gC/gN ] Leaf Carbon-Nitrogen ratio 
%Ns  = [50 50];%%%[ 50-330] Sapwood Carbon Nitrogen  [gC/gN] Sapwood
%N = [58 58]; %% [30- 60] [gC/gN] Fine root  Carbon Nitrogen
%Nc = Ns 
Stoich.Nl = Nl;  %% Leaf Carbon Nitrogen 
Stoich.Ns= Nl/0.145;%%%[ 50 -330] Sapwood Carbon Nitrogen  [gC/gN] Sapwood
Stoich.Nr= Nl/0.860; %%% [30- 60] [gC/gN] Fine root  Carbon Nitrogen
Stoich.Nf = Nl;  %%%  [gC/gN]  Fruit/Reproduction  
Stoich.Nh = Nl/0.145;  %%%  [gC/gN] Heartwood/Dead sapwood 
%%%% 
Pl = Nl*14;
Stoich.Phol = Pl; %%%  [gC/gP]
Stoich.Phos = Pl/0.145;
Stoich.Phor = Pl/0.860; 
Stoich.Phof = Pl; 
Stoich.Phoh= Pl/0.145; 
%%%%%
Kl = Nl*2;
Stoich.Kpotl = Kl;  %%  [gC/gK]
Stoich.Kpots  = Kl/0.145;
Stoich.Kpotr = Kl/0.20;
Stoich.Kpotf = Kl;
Stoich.Kpoth = Kl/0.145;
%%%%%%%%%%%%%%%% Max. Translocation Rates of nutrients 
Stoich.ftransR = 0.2; 
Stoich.ftransL = 0.6; 
%%%%%% 
Stoich.FiS=1; %%[-] ; %%% Factor to increment or decrement nutrient reserve buffer 
%%%%%%
%%%%% Lignin Fraction 
%% Poorter 1994; Chave et al 2009 ; Roumet et al 2015 ; Fortunel et al 2009
Stoich.Lig_fr_l= 0.15; %% Lignin fraction in leaves  [g Lignin / g DM] 
Stoich.Lig_fr_fr=0.10; %%% Lignin fraction in fruit [g Lignin / g DM] 
Stoich.Lig_fr_h= 0.25; %% Lignin fraction in wood [g Lignin / g DM] 
Stoich.Lig_fr_r= 0.10; %% Lignin fraction in fine roots [g Lignin / g DM] 
return
