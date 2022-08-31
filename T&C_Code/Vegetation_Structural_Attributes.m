%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Vegetation Structural Evolution    %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Ccrown_t,Hc,SAI,BA,Tden,AgePl,TBio]=Vegetation_Structural_Attributes(dt,Ccrown_tm1,B,fab,Tden_tm1,AgePl_tm1,OPT_ALLOME)
%%%%%%%%%%%%%%%%
%%% INPUT
%Ccrown_tm1 %%[-]
TBio = (B(1) + B(2) + B(3) + B(4) + B(6)); %%% [gC/m2 VEG] Total Biomass
TBio = 2*TBio/100 ;  %[ton DM / ha VEG]
Babove = fab*(B(2) + B(4) + B(6));  %%% [gC/m2 VEG] Stem Aboveground Biomass
Babove2 = (B(1) + fab*(B(2) + B(4) + B(6)));  %%% [gC/m2 VEG] Leaf+ Stem Aboveground Biomass
%AgePl [yr] Age_Plantation
%Tden= 160; %% tree density [palm/ha]
%%% OUTPUT
%Ccrown_t
%%%%%%%%%%%%%%%%%%%%%%%%%%%
AgePl = AgePl_tm1 + dt; %%% [days]
%%%%%%%%%%%
switch OPT_ALLOME
    case 1
        %%%% OPTION OPAL
        %%% Per unit of area 
        Babove2 = (B(1) + fab*(B(2) + B(4) + B(6)))*Ccrown_tm1;  %%% [gC/m2 area] Leaf+ Stem Aboveground Biomass
        %%%% OIL PALM
        Tden = Tden_tm1; %%% [number /ha]
        %%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%
        %%% Option 1
        %CD = (1/0.2733)*log((AgePl/365.25)/0.7344) ;
        %Hc = 1.1444*exp(0.2455*CD);
        %%% Option 2
        %BMp = 2*Babove2/100/Tden; %%% [Mg DM /tree]
        %Hc=(BMp./0.1839).^(1/0.766);  % PALM HEIGHT, h [m] vs Palm BIOMASS (total), BM [MgDM/palm] - Khasanah et al. 2015
        %%% Option 3
        BMp = 1000*2*Babove2/100/Tden; %%% [kg DM /tree]
        Hc = (BMp+7.0872)/71.797; % PALM HEIGHT, h [m] Asari et al 2013
        %%%%
        CD=1/0.5736.*log(Hc/0.0161);   % [m] h trunk
        %%%%
        CA = pi.*CD.^2/4; %% [m2 per tree]
        CA = CA.*Tden/10000; %%[ ha/ha]
        CA(CA>1)=1;
        Ccrown_t = CA;
        %%%%%%%%%%%%%%%%%%%
        D = 0.2;%%% [m] diamter tree
        fv = 0.7; %%%  fraction of vertical stem
        SAI = ((1-fv)*D*Hc + fv*pi*D^2/4)*(Tden/10000); %%[m2 /m2]
        SAI =SAI/Ccrown_t; %%[m2 SAI /m2 PFT]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Stem area
        BA = (pi*D^2/4)*(Tden);  %[m2 /ha]
    case 2
        %%%% Generic Forest  - Currently no change in Crown-Area
        Tot_Biomass = max(0,TBio) ;  %25:1:500; %275; %% [MgDM / ha] %%%  [ton DM / ha ]
        TB = Tot_Biomass*1000; %% [kg DM / ha ]
        %%% M = [kg DM / tree ]
        %%% N= [tree /ha]
        %%% TB = M*N; [kgDM / ha]
        %%%% Self thinning law  Yoda et al 1963 Wolf et al 2011 GCB
        %M=alpha*N.^beta;
        alpha = 10^(6.39); %% Ang. 6.22 Gym 6.62
        beta = -1.38; %% Ang. -1.32 -1.45  %%% Theory: -4/3 -3/2
        %%%%%
        %M=alpha*N.^beta;
        %N= (1/alpha)^(beta^-1).*M.^(beta^-1);
        %TB=alpha*N.^(beta+1);
        Tden = (1/alpha)^((beta+1)^-1).*TB.^((beta+1)^-1); %% [n°ind / ha ]
        %%%%%
        %%%%%%%%%% Allometric Scaling of tree
        %%%%%%%
        w_den = 450*1000; %%[gDM /m^3]
        R_C = 0.5; %% [gC/gDM] Ratio Carbon-Dry Matter in trunk
        fv = 0.7; %%%  fraction of vertical stem
        %%%%%%%%% Starting from Biomass --->
        %%%%
        Vtree = Babove./w_den./R_C./(Tden/10000) ; %%% [m^3]
        %%%%
        Mtree = Vtree*w_den*1000; %% [kg-DM]
        %%Mtree [kg-DM tree]  Aboveground mass tree
        %D = (1/8659)*Mtree^0.4134; %%% [m]  Cannell 1982
        %%% geometrical consideation
        %D = (6./(40*pi).*Vtree).^(1/(2+0.5)); %%% [m] Diameter
        D = (6./(54.64*pi).*Vtree).^(1/(2+0.7034)); %%% [m] Diameter
        %%%%%%%%%%%
        %Hc = 40*D.^0.5; %% [m] %%% Sitch et al., 2003  Sato eto al., 2007 West et al., 2009
        Hc = 54.64*D.^0.7034; %% [m] %%% Schepaschenko
        %%%%
        Acrown = 228.6*D^1.53; %%% BAAD Database
        %Acrown = 150*D.^1.6; %%[m^2 PFT]  %%% Sitch et al., 2003  Sato eto al., 2007 West et al., 2009
        CA = (Tden/10000).*Acrown; %% [m^2 PFT / m^2 ground]
        CA(CA<0.01)=0.01;
        CA(CA>1)=1;
        %Ccrown_t = CA;
        Ccrown_t = 1;
        %%%%%
        SAI = ((1-fv)*D*Hc + fv*pi*D^2/4)*(Tden/10000);  %[m2 /m2]
        SAI =SAI/Ccrown_t; %%[m2 SAI /m2 PFT]
        %%% Stem area / Basal Area
        BA = (pi*D^2/4)*(Tden);  %[m2 /ha]
    case 3
        %%% Specific forest with known allometric relation
    case 4
        %%% Eucalyptos Regnans
        AgeE =AgePl/365.25; %% yr
        if AgeE>1.5
            Hc= 1./(0.012807792+(6.1243119*log(AgeE)./AgeE.^2));  %% [m]
            %TBio = 2./(0.0030239638+(2.4875134*log(AgeE)./AgeE.^2)); %% Aboveground Biomass - DM [ton DM/ha]
            %BaboveD= 2*Babove/100; %% [ton DM/ha]
            %x=(2-BaboveD*0.0030239638)./(2.4875134*BaboveD);
            %Hc= 1./(0.012807792+(6.1243119*x));  %% [m]
        else
            Hc=0.3+0.6/1.5*AgeE;
            %TBio=  4.4/1.5*AgeE;
        end
        BA= NaN;
        Tden = NaN;
        SAI=0.2;
        Ccrown_t=1;
end
%%%%%
end
