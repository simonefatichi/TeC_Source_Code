%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Litter and Fire                         %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[B,LitFirEmi]=Litter_Fire(B,FireA)
fire_eff= FireA;  
funb_nit = 0.15; % Mpar.funb_nit; %% Fraction of unburned Nitrogen in Litter 
%%%%%%%%%%%%%%
%%% B1 Above-ground Litter Metabolic
%%% B2 Above-ground Litter Structural - Cellulose/Hemicellulose
%%% B3 Above-ground Litter Structura - Lignin
%%% B4 Above-ground Woody  - Cellulose/Hemicellulose
%%% B5 Above-ground Woody - Lignin
%%% B23 Nitrogen Above-ground Litter
%%% B24 Nitrogen Above-ground Woody
%%% B35 phosphorus Above-ground Litter
%%% B36 phosphorus Above-ground Woody
%%% B48 Potassium Above-ground Litter
%%% B49 Potassium  Above-ground Woody
%%%%%%%%%%
%%% B31 Nitrogen Ione Ammonium NH4+
%%% B43 phosphorus Mineral
%%% B52 Potassium  Mineral  solution
%%%%  Mycorrhizal
%%% B20 AM-Mycorrhizal - C
%%% B21 EM-Mycorrhizal - C
%%% B29 AM Mycorrhizal - N
%%% B30 EM Mycorrhizal - N
%%% B41 AM - Mycorrhizal - P
%%% B42 EM - Mycorrhizal - P
%%%%%%%
LitFirEmi(1)=fire_eff*(B(1)+B(2)+B(3)+B(4)+B(5)); %% Burned Litter C 
B(1)=B(1)*(1-fire_eff); 
B(2)=B(2)*(1-fire_eff); 
B(3)=B(3)*(1-fire_eff); 
B(4)=B(5)*(1-fire_eff); 
B(5)=B(5)*(1-fire_eff); 
%%%%%
B(31)=  B(31) + (funb_nit)*(B(23)+B(24))*fire_eff; 
LitFirEmi(2)= (1-funb_nit)*(B(23)+B(24))*fire_eff; %%% Burned Litter Nitrogen 
B(23)=B(23)*(1-fire_eff); 
B(24)=B(24)*(1-fire_eff); 
%%%%%%%%
B(43)= B(43)+ (B(35) + B(36))*fire_eff; %%% P-burned 
B(35)=B(35)*(1-fire_eff); 
B(36)=B(36)*(1-fire_eff); 
%%%%%%%%
B(52)= B(52)+ (B(48) + B(49))*fire_eff; %%% K-burned 
B(48)=B(48)*(1-fire_eff); 
B(49)=B(49)*(1-fire_eff); 
%%%%%%% Mycorrhizal - effect 
%%%%%%%
B(B<0)=0;
end 