%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Throughfall                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Cfol_H,Cfol_L,CLitter]=Throughfall(LAI_H,SAI_H,LAI_L,SAI_L,Llitter,Kc)
%%% References [Van Dijk and Bruinzeel (2001)] ; [Ramirez and Senarath
%%% 2000] [Mahfouf and Jacquemin (1989)] 
%%%INPUTS
%%% LAI [] Leaf Area Index 
%%% SAI [] Stem Area Index 
%%% Llitter [] Litter Area Index 
%%% OUTPUTS
%Cfol_H [0-1] Foliage Cover H Vegetation 
%Cfol_L [0-1] Foliage Cover L Vegetation 
%Clitter [0-1] Litter Cover   
%%%%%%%%
%Kc = 0.75; 
%%%%%%%%%%%%%%%%%%% Vegetation Cover First and Second Layer %%%%%%%%%%%%
Cfol_H= 1 - exp(-Kc*(LAI_H+SAI_H)); %First Layer 
Cfol_L= 1 - exp(-Kc*(LAI_L+SAI_L)); %Second Layer 
CLitter = 1 - exp(-Kc*(Llitter)); % Litter Layer 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return 