%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Evaporation_layers   %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[EvL_Zs]=Evaporation_layers(Zs,Zdes)
%%%% INPUT
%Zs, [mm] Depth Layers [1....m]
%Zdes [mm] Depth of desorption 
%%% OUTPUT
%EvL_Zs, [%] Evaporation Layer fraction [1...m]
n=length(Zs)-1;
EvL_Zs= zeros(1,n);
if Zdes < Zs(1) 
    disp('ERROR FIRST LAYER TOO DEPTH') 
    return 
end 
if Zdes > Zs(n+1) 
    disp('ERROR LAST LAYER TOO SHALLOW') 
    return 
end 
NiL= sum(Zs < Zdes); %% Number interested Layer 
dz= diff(Zs); %%%% [mm]  Thickness of the Layers
dz(NiL)= Zdes - Zs(NiL); %%% Thickness last layer 
dz=dz(1:NiL); %% Interested Layer 
EvL_Zs(1:NiL) = dz/Zdes; 
end 
