%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Water Contents into Soil Dunne Runoff  Lateral Subsurface flow 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[O,ZWT,OF,OS,Psi_s_H,Psi_s_L,gsr_H,gsr_L,Exwat_H,Exwat_L,Rd,WTR,POT,OH,OL]=Soil_Water_MultiLayer(V,Zs,...
    dz,n,Ccrown,Osat,Ohy,nVG,alpVG,lVG,Ks_Zs,L,Pe,O33,Ks_mac,Omac,alpVGM,nVGM,lVGM,s_SVG,bVG,Phy1,SPAR,EvL_Zs,Inf_Zs,RfH_Zs,RfL_Zs,...
    Rrootl_H,Rrootl_L,PsiL50_H,PsiL50_L,PsiX50_H,PsiX50_L,Ts,Tdp,Psi_sto_00_H,Psi_sto_50_H,Psi_sto_00_L,Psi_sto_50_L,...
    Salt,Osm_reg_Max_H,Osm_reg_Max_L,eps_root_base_H,eps_root_base_L)
%%% INPUT
%%%V -- Volume [mm] in the Layer 
%%% Ks --- Saturation Conductivity Vertical [mm/h]
%%% Kh, -- Saturation Conductivity Horizontal [mm/h] 
%Damp_Zs, --- Depth Soil Heat exchange fraction [1...m]
%Zdes,...   -- Depth for Evaporation Process [mm]  
%Osat, Water content Saturation []  
%Ohy, Water Content Hygrosopic [] 
%Owp, Water Content wilting point []
%EvL_Zs, [%] Evaporation Layer fraction [1...m]
%RfH_Zs, [%] Root Fraction for First Layer of Vegetation [1...m]x[1....cc]
%RfL_Zs  [%] Root Fraction for Second Layer of Vegetation [1...m]x[1....cc]
%dz= diff(Zs); %%%% [mm]  Thickness of the Layers
%n = length(V); 
%Zs, [mm] Depth Layers [1....m+1] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUT  
%O, [] Water Content [1...m] Layer  
%OF, [] Water Content for infiltration 
%OS, [] Water Content for evaporation 
%OZdep, [] Water Content for soil heat exchange 
%OH, [] Water Content for First Layer of Vegetation
%OL,[] Water Content for Second Layer of Vegetation 
%Rd, [mm] Dunne Runoff 
%WTR [mm] Water Table Rise 
%POT [mm] Soil Water Potential 
%Psi_s_H,  [MPa] Soil Water Potential  for First Layer of Vegetation
%Psi_s_L, [MPa] Soil Water Potential  for Second Layer of Vegetation
%Ko_H, [mm/h]
%Ko_L, [mm/h]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%% Definition of -- OF OS OZdep OH  OL
%%% Zs(i)  dz(i)  RfL_Z(i)  RfH_Z(i)  WEG(i)
O=ones(1,n); %%% Water Content [] n Layer ---
%%%%%%% Definition Layer Depth %%%% and Water Content
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%
WTR= zeros(1,n); 
for i=n:-1:1
    if i == n
        O(i) = (V(i)/dz(i))+ Ohy(i); %%%%  Water Content [] All layers
        WTR(i) = (O(i) - Osat(i))*dz(i)*(O(i) > Osat(i)); %%% [mm] Water Table Rise 
    else
        O(i) = (V(i)+ WTR(i+1))/dz(i) + Ohy(i) ; %%%%  Water Content [] All layers
        WTR(i) = (O(i) - Osat(i))*dz(i)*(O(i) > Osat(i)); %%% [mm] Excess From the Reservoir - Below
    end
    if O(i) < Ohy(i)
        O(i) = Ohy(i);
    end
    if O(i) > Osat(i)
        O(i)=Osat(i);
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%
i=n;
while (O(i) > Osat(i)-1e-5) %% Computation precision 
    i=i-1; 
    if i == 0
        break
    end 
end
ZWT=Zs(i+1); 
PHead = (Zs(1:n)+dz/2)-ZWT;
PHead(PHead<0)=0; 
%%%  Compute First Potential 
switch SPAR
    case 1
        %%% Hydraulic Head
        Se = (O-Ohy)./(Osat-Ohy);
        mVG= 1-1./nVG;
        %%%%%%%%%%%%%%%
        Phy = Phy1*101.9368; %% [mm] 
        P1 = (1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
        P2 = -Phy*exp(Se.*-abs(bVG)); % [mm]
        P = (Se<s_SVG).*P2 +(Se>=s_SVG).*P1; % [mm]
        P(Se==0)=P2(Se==0); 
        POT = PHead + P; %%% [mm]
    case 2
        %%% Hydraulic Head
        Ptem=zeros(1,n);
        for jk=1:n
            [~,Ptem(jk)]=Conductivity_Suction(2,Ks_Zs(jk),Osat(jk),Ohy(jk),L(jk),Pe(jk),O33(jk),alpVG(jk),nVG(jk),lVG(jk),...
                Ks_mac(jk),Omac(jk),alpVGM(jk),nVGM(jk),lVGM(jk),Phy1,s_SVG(jk),bVG(jk),O(jk));
        end
        POT = PHead - Ptem ; %%% [mm]
    case 3
        Ptem=zeros(1,n);
        for jk=1:n
            [~,Ptem(jk)]=Conductivity_Suction(3,Ks_Zs(jk),Osat(jk),Ohy(jk),L(jk),Pe(jk),O33(jk),alpVG(jk),nVG(jk),lVG(jk),...
                Ks_mac(jk),Omac(jk),alpVGM(jk),nVGM(jk),lVGM(jk),Phy1,s_SVG(jk),bVG(jk),O(jk));
        end
        POT = PHead - Ptem ; %%% [mm]
    case 4
        %%% Hydraulic Head
        Se = (O-Ohy)./(Osat-Ohy);
        mVG= 1-1./nVG;
        POT = PHead + (1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
end
%%%%%%%%%%%%%%%%%%%%%%%%
if i==n
    ZWT =Zs(i+1);
else
    if i == 0
        ZWT = 0;
    else
        if POT(i+1)>0 %% Usual case
            %%% Find Water Table Depth
            CZ = Zs(1:n)+dz/2;
            ZWT = CZ(i+1)+ (CZ(i)-CZ(i+1))*(-POT(i+1))/(POT(i)-POT(i+1)); %% [mm]
        else %%% Otherwise Soil Heterog. in Pe creates lower POT below WT  
            ZWT =Zs(i+1);
        end
        %%% Recompute Potential
        PHead = (Zs(1:n)+dz/2)-ZWT;
        PHead(PHead<0)=0;
        %POT = PHead + (1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
    end
end
%%%%%%% Recompute Potential
switch SPAR
    case 1
        %%% Hydraulic Head
        Se = (O-Ohy)./(Osat-Ohy);
        mVG= 1-1./nVG;
        %%%%%%%%%%%%%%%
        Phy = Phy1*101.9368; %% [mm]
        P1 = (1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
        P2 = -Phy*exp(Se.*-abs(bVG)); % [mm]
        P = (Se<s_SVG).*P2 +(Se>=s_SVG).*P1; % [mm]
        P(Se==0)=P2(Se==0); 
        POT = PHead + P; %%% [mm]
    case 2
        %%% Hydraulic Head
        Ptem=zeros(1,n);
        for jk=1:n
            [~,Ptem(jk)]=Conductivity_Suction(2,Ks_Zs(jk),Osat(jk),Ohy(jk),L(jk),Pe(jk),O33(jk),alpVG(jk),nVG(jk),lVG(jk),...
                Ks_mac(jk),Omac(jk),alpVGM(jk),nVGM(jk),lVGM(jk),Phy1,s_SVG(jk),bVG(jk),O(jk));
        end
        POT = PHead - Ptem ; %%% [mm]
    case 3
        Ptem=zeros(1,n);
        for jk=1:n
            [~,Ptem(jk)]=Conductivity_Suction(3,Ks_Zs(jk),Osat(jk),Ohy(jk),L(jk),Pe(jk),O33(jk),alpVG(jk),nVG(jk),lVG(jk),...
                Ks_mac(jk),Omac(jk),alpVGM(jk),nVGM(jk),lVGM(jk),Phy1,s_SVG(jk),bVG(jk),O(jk));
        end
        POT = PHead - Ptem ; %%% [mm]
    case 4
        %%% Hydraulic Head
        Se = (O-Ohy)./(Osat-Ohy);
        mVG= 1-1./nVG;
        POT = PHead + (1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
end
%%%%%%%%%%%%%%%%%%%%%%%%
Rd = WTR(1); %%% [mm]  Dunne Runoff
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Evaporation Layer Water Content 
OS = sum(EvL_Zs.*O);%% Evaporation Bare Soil WC []
%%%%%%%%%%%%%%%%%% First Layer 
OF = sum(Inf_Zs.*O); %%% Infiltration Water Content []
%%%%%%%%%%%%%%%%% VEGETATION 
[cc,n]=size(RfH_Zs); 
OH=zeros(1,cc); OL=zeros(1,cc);
Psi_s_H=zeros(1,cc); Psi_s_L=zeros(1,cc);
for i=1:cc
    OH(i) = sum(RfH_Zs(i,:).*O);
    OL(i) = sum(RfL_Zs(i,:).*O);
    [~,Psi_s_H(i)]=Conductivity_Suction(SPAR,sum(RfH_Zs(i,:).*Ks_Zs),sum(RfH_Zs(i,:).*Osat),sum(RfH_Zs(i,:).*Ohy),...
        sum(RfH_Zs(i,:).*L),sum(RfH_Zs(i,:).*Pe),sum(RfH_Zs(i,:).*O33),sum(RfH_Zs(i,:).*alpVG),sum(RfH_Zs(i,:).*nVG),sum(RfH_Zs(i,:).*lVG),...
        sum(RfH_Zs(i,:).*Ks_mac),sum(RfH_Zs(i,:).*Omac), sum(RfH_Zs(i,:).*alpVGM),sum(RfH_Zs(i,:).*nVGM),sum(RfH_Zs(i,:).*lVGM),...
        Phy1,sum(RfH_Zs(i,:).*s_SVG),sum(RfH_Zs(i,:).*bVG),OH(i));
    [~,Psi_s_L(i)]=Conductivity_Suction(SPAR,sum(RfL_Zs(i,:).*Ks_Zs),sum(RfL_Zs(i,:).*Osat),sum(RfL_Zs(i,:).*Ohy),...
        sum(RfL_Zs(i,:).*L),sum(RfL_Zs(i,:).*Pe),sum(RfL_Zs(i,:).*O33),sum(RfL_Zs(i,:).*alpVG),sum(RfL_Zs(i,:).*nVG),sum(RfL_Zs(i,:).*lVG),...
         sum(RfL_Zs(i,:).*Ks_mac),sum(RfL_Zs(i,:).*Omac), sum(RfL_Zs(i,:).*alpVGM),sum(RfL_Zs(i,:).*nVGM),sum(RfL_Zs(i,:).*lVGM),...
          Phy1,sum(RfL_Zs(i,:).*s_SVG),sum(RfL_Zs(i,:).*bVG),OL(i));
    %%%%%%%%%%%
end
Psi_s_H=-(Psi_s_H/1000)*1000*9.81/1e+6; %%[MPa]
Psi_s_L=-(Psi_s_L/1000)*1000*9.81/1e+6; %%[MPa]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Salinity effect
if Salt>0
    [Tdp_H,Tdp_L]=RootZone_Temp(Tdp,RfH_Zs,RfL_Zs);
    for i=1:length(Ccrown)
        [Psi_s_H(i)]=Salinity_Plant(Salt,Tdp_H(i),Ts,Psi_s_H(i),Psi_sto_00_H(i),Psi_sto_50_H(i),Osm_reg_Max_H(i),eps_root_base_H(i));
        [Psi_s_L(i)]=Salinity_Plant(Salt,Tdp_L(i),Ts,Psi_s_L(i),Psi_sto_00_L(i),Psi_sto_50_L(i),Osm_reg_Max_L(i),eps_root_base_L(i));
    end
end  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rho2 = 55555555; %% [mmolH20 /m^3]; Water density
%rcyl= 2.0*1e-3 ;%%%  [m] radius cylinder of soil to which root has access to
rroot= 0.5*1e-3 ;%% [0.5-6 *10^-3] [m] root radius
Psi_s= zeros(n,1); 
gsr_L=zeros(cc,n);
gsr_H=zeros(cc,n);
%%%%% layer by layer analysis 
for jk=1:n
    [Ko,Psi_s(jk)]=Conductivity_Suction(SPAR,Ks_Zs(jk),Osat(jk),Ohy(jk),L(jk),Pe(jk),O33(jk),alpVG(jk),nVG(jk),lVG(jk),...
        Ks_mac(jk),Omac(jk),alpVGM(jk),nVGM(jk),lVGM(jk),Phy1,s_SVG(jk),bVG(jk),O(jk)); 
    for i=1:cc
        [gsr_L(i,jk)]= root_soil_Conductance(Ko,RfL_Zs(i,jk).*Rrootl_L(i).*Ccrown(i),rroot,dz(jk)); % % [mmol H20 / m^2 ground s MPa]
        [gsr_H(i,jk)]= root_soil_Conductance(Ko,RfH_Zs(i,jk).*Rrootl_H(i).*Ccrown(i),rroot,dz(jk)); % [mmol H20 / m^2 ground s MPa]
    end
end
Psi_s=-(Psi_s/1000)*1000*9.81/1e+6; %%[MPa]
%%%%%%%%%%
%Psi_minH = nanmin(PsiX50_H,PsiL50_H); 
%Psi_minL = nanmin(PsiX50_L,PsiL50_L); 
Psi_minH = min([PsiX50_H;PsiL50_H],[],'omitnan'); 
Psi_minL = min([PsiX50_L;PsiL50_L],[],'omitnan');  
Exwat_L = gsr_L/rho2*1000*3600.*(-Psi_minL'*ones(1,n) + ones(cc,1)*(Psi_s')); %%  %% [mm m2 / m2 ground h ] %% Max extractable water 
Exwat_H = gsr_H/rho2*1000*3600.*(-Psi_minH'*ones(1,n) + ones(cc,1)*(Psi_s')); %%  %% [mm m2 / m2 ground h ] %% Max extractable water 
Exwat_L(Exwat_L<0)=0;
Exwat_H(Exwat_H<0)=0;
gsr_L = sum(gsr_L,2); % [mmol H20 / m^2 ground s MPa]
gsr_H = sum(gsr_H,2); % [mmol H20 / m^2 ground s MPa]
return
