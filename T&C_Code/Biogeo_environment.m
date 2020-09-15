%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Temperature Humidity Biogeochemistry Zone
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Se_bio,Se_fc,Psi_bio,Tdp_bio,Vtot,VT]=Biogeo_environment(Tdp,O,V,Soil_Param,Phy,SPAR,Bio_Zs)
O=mean(O,1); Tdp=mean(Tdp,1);  V=mean(V,1); 
%%%%%%%
Ks_Zs=Soil_Param.Ks_Zs; 
Osat=Soil_Param.Osat;
Ohy =Soil_Param.Ohy; 
L=Soil_Param.L;
Pe=Soil_Param.Pe; 
O33=Soil_Param.O33; 
alpVG=Soil_Param.alpVG; 
nVG=Soil_Param.nVG; 
lVG=Soil_Param.lVG; 
Ofc=Soil_Param.Ofc; 
Ks_mac=Soil_Param.Ks_mac; 
Omac=Soil_Param.Omac; 
alpVGM=Soil_Param.alpVGM;
nVGM =Soil_Param.nVGM;
lVGM =Soil_Param.lVGM;
s_SVG =Soil_Param.s_SVG;
bVG =Soil_Param.bVG;
%%%%%%%%
%%%% Biogeochemistry Layer Water Content
Obio = sum(Bio_Zs.*O);
Ohy = sum(Bio_Zs.*Ohy); 
Osat = sum(Bio_Zs.*Osat); 
Ofc = sum(Bio_Zs.*Ofc); 
Obio(Obio<=Ohy)= Ohy + 1e-5; 
Obio(Obio>=Osat)= Osat - 1e-5; 
Se_bio = (Obio-Ohy)./(Osat-Ohy);%%
Se_fc = (Ofc-Ohy)./(Osat-Ohy);%%
[Ko,Psi_bio]=Conductivity_Suction(SPAR,sum(Bio_Zs.*Ks_Zs),Osat,Ohy,...
    sum(Bio_Zs.*L),sum(Bio_Zs.*Pe),sum(Bio_Zs.*O33),sum(Bio_Zs.*alpVG),sum(Bio_Zs.*nVG),sum(Bio_Zs.*lVG),...
     sum(Bio_Zs.*Ks_mac),sum(Bio_Zs.*Omac), sum(Bio_Zs.*alpVGM),sum(Bio_Zs.*nVGM),sum(Bio_Zs.*lVGM),Phy,sum(Bio_Zs.*s_SVG),sum(Bio_Zs.*bVG),Obio);
Tdp_bio = sum(Bio_Zs.*Tdp); %% [°C]
Psi_bio=-(Psi_bio/1000)*1000*9.81/1e+6; %%[MPa]
Vtot = sum(V); %% [mm] Volume of entire depth 
VT = sum(Bio_Zs.*V); %% [mm] Volume in the Biogeochemistry layer
end
