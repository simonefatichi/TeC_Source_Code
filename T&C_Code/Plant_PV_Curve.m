%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction - Plant Psi-V Curve         %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [Psi] = Plant_PV_Curve(V,C,TYP,LAI,Axyl,H,Sl)

%%% if C length == 1 
%% C --> Cx or Cl, 
%Cl  = 12'00;  %%%  %  [mmolH20 / m^2 leaf MPa]
%Cx = 150; %%% [kg / m^3 sapwood MPa]
%%% if C length == 2 
%% C  [Px12 Px50] or [Eleaf LDMC] 

rho= 55.5; %% [mmolH20 /cm^3] Water density
%%%  Sapwood 
fwat = 0.30; % ; %% Free Water in Xylem  %%% Free Water is less 0.25 --
n_sap = 0.60; %% 0.55-0.65 [-] Overall Sapwood-porosity  Density drw-wood
%Axyl = 20  ; %% % [cm^2 stem /m^2 PFT]
%H= 10; %% height [m]
%%% Leaf 
%LAI = 1 ; %% [Leaf Area]
%Sl = 0.020 %% [m^2 leaf / gC ]
LDMC = 0.40 ; %% [ gDM / g leaf Fresh] %%% Vile et al., 2005   [0.22-0.50] Scoffoni et al., 2011  Niinemets 2001
LMA =2./Sl; %% [gDM / m^2 Leaf ]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(TYP,'X')
    TYPN=1;
end
if strcmp(TYP,'L')
    TYPN=2;
    if length(C) == 2 
        LDMC = C(2);
    end 
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch TYPN
    case  1  %%% XYLEM 
        Vx=V;%% [mmolH20/ m2 PFT]
        Vwat = rho*fwat*Axyl*H*100; %% whole plant [mmolH20  / m^2 PFT ]
        RWC = (Vx)./(Vwat); %%%[-]
        RWC(RWC>1)=1; RWC(RWC==0)=1e-5;
        if length(C) == 1
            Cx=C; %%% Capacitance %[kg H20 /m^3 sapwood MPa]
            a = -1000/(Cx*(1 +(1-n_sap)./fwat));
            b = 0;
        else
            Px12= C(1);
            Px50= C(2);
            %%% Barnard et al., 2011 PCE  doi:10.1111/j.1365-3040.2010.02269.x
            %Px50=-2.5 ; % -1.0; %%%
            %Px12=-0.24*2.5; % -0.14;  %%%
            %%%%
            b= (Px12 -0.24*Px50)/(0.12*(Px50-Px12)); %%[-]
            a= Px50*(2+b);   %%[MPa]

        end
        
        %%%%
        Psi = a*(1-RWC)./(1+b*(1-RWC)); %%% [MPa]
        Psi(Psi<-10) = -10;
        %%%%%%%%
        %VxF_Vwat = 1 +(1-n_sap)./fwat; %% [-]
        %%% Recompute Capacitance
        %Cap = (0.001*( -a./((b*Psi-a).^2)))./VxF_Vwat*1e+6; %%[kg H20 /m^3 sapwood MPa]
        
        
    case 2 %%% LEAF 
        Vl=V;% %% [mmolH20/ m2 PFT]
        %%%%% Maximum Leaf Volume
        Vleaf =  LAI*LMA*(1- LDMC)/LDMC; %%% [gH20 /m^2 PFT]
        RWC = Vl./(Vleaf*rho) ;  %%% [-]
        RWC(RWC>1)=1; RWC(RWC==0)=1e-5;
        
        if length(C) == 1
            %%% Cl [500 - 3000]%  Leaf capacitance [mmolH20 / m^2 leaf MPa]
            Cl=C*LAI; %%%    Leaf capacitance Plant scale [mmolH20 / m^2 PFT MPa]
            Eleaf = rho*Vleaf./Cl; %%% [MPa] 
            %Psi =  -(rho*Vleaf)*(1 - RWC)./C;
        else
            Eleaf= C(1);
            %C = rho*Vleaf./Eleaf; %% %   [mmolH20 leaf / m^2 leaf Mpa]
        end
        Psi =  -Eleaf*(1 - RWC);
        Psi(Psi<-10) = -10;    
        
end

Psi(Psi>0) = 0;

end