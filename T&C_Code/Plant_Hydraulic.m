%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction - Plant Hydraulic           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%REFERENCES
function[Psi_x,Psi_l,Jsx,Jxl,Kleaf,Kx,Vx,Vl]=Plant_Hydraulic(Xout,X0,dth,Ccrown,T,Psi_s,hc,LAI,Axyl,PsiL50,PsiL00,Kleaf_max,Kx_max,PsiX50,Sl,mSl,Cx,Cl)

%%%%%%%%%%%%%%%%%
%Psi_x =Xout(end,1); 
%Psi_l = Xout(end,2);
Vx = Xout(end,1); Vl =Xout(end,2); %% [mmolH20/ m2 PFT]
Vxtm1=X0(1);  Vltm1=X0(2);  %%% [mmolH20/ m2 PFT]

% Plant Pressure Volume curve
[Psi_x] = Plant_PV_Curve(Vx,Cx,'X',LAI,Axyl,hc,Sl); 
[Psi_l] = Plant_PV_Curve(Vl,Cl,'L',LAI,Axyl,hc,Sl); 


Vxtm1(Vxtm1<0) = 0; Vltm1(Vltm1<0)=0;
Vx(Vx<0) = 0; Vl(Vl<0)=0;

%%% INPUTS
% T Transpiration [mm/h]
% hc Canopy height [m]
% LAI  Leaf Area [m2 LAI / m2 PFT ]
% dt [s]
% Psi_s [MPa] Soil water potential felt by plant
% Vx %%% [mmolH20/ m2 PFT]  Volume water xylem
% Vl %%% [mmolH20/ m2 PFT] Volume water leaves
% Psi_ltm1 % [MPa]  Water potential in the leaves
% Psi_xtm1 % [MPa] Water potential in the xylem
%gsr  % soil to root conductance [mmol H20 / m^2 PFT s MPa]
%%% OUTPUTS
% Vx %%% [mm / m2 PFT]  Volume water xylem
% Vl %%% [mm/ m2 PFT] Volume water leaves
%Psi_x,  % [MPa] Water potential in the xylem
%Psi_l, % [MPa]  Water potential in the leaves
%Jsx,%% [mmolH20 /m^2 PFT s ] Soil to Xylem flux
%Jxl, %% [mmolH20 /m^2 PFT s ] Xylem to leaves flux
%Kleaf, %% [mmolH20/ MPa s m^2 PFT ] Leaf conductance
%Kx, %% [mmolH20 /m^2 PFT s MPa] Xylem Conductance

%%% Parameters
row = 1000; % water density [kg/m^3]
rho2 = 55555555; %% [mmolH20 /m^3]; Water density
%%%%%
T = T./Ccrown; %% Transpiration Canopy  %% [mm m2 / m2 PFT h]
Vx = Vx/rho2*1000; %%  %% [mm m2 / m2 PFT ]
Vl =  Vl/rho2*1000; %%  %% [mm m2 / m2 PFT ]
Vxtm1 =  Vxtm1/rho2*1000; %%  %% [mm m2 / m2 PFT ]
Vltm1 =  Vltm1/rho2*1000; %%  %% [mm m2/ m2 PFT ]
%%%%%%%%%%%%%%%%%%%%%%%%
[PLD]=Plant_Disconnection(Psi_x,Psi_l,Axyl,PsiL50,PsiL00,PsiX50); 
%%%%%%%%%%%%%%%%
if LAI == 0 || PLD == 1
    if Axyl > 0
        Psi_l= Psi_s; 
        Psi_x= Psi_s;
    else
        Psi_l = Psi_s;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%Axyl = 15.0 ; %% [cm^2 stem /m^2 PFT]
Axyl =  Axyl/10000; %% [m2 stem/m2 PFT]
%%% Leaf
%PsiL50 =  -0.7 ;%%[MPa]  Water Potential at 50% loss conductivity
%PsiLs =  -1.6; %% [MPa]  Water Potential at PLCs% loss conductivity
PLCs = 0.02;  %%[-]  specific loss of conductivity
%Kleaf_max = 10 ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
%Cl  = 1;  %%%  %  [mmolH20 / m^2 leaf MPa]
%Cl = Cl*LAI ; %% [mmolH20 / m^2 PFT MPa]
%%%%%%%%%%
%%% gsr ; %% [mmol H20 / m^2 PFT s MPa]
%%% Xylem
%Kx_max = 80000;  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
%PsiX50 = -3.5; %%[MPa]  Water Potential at 50% loss conductivity
PsiXe = PsiX50/1.42; %% [MPa] (Typically at 12% of Xylem, Meinzer et al., 2010)    Meinzer et al., 2009
%Cx = 150; %%% [kg / m^3 sapwood MPa]
%Cx = Cx/row*rho2*hc*Axyl; %%  % [mmolH20 / m^2 PFT MPa]

%%% Leaf conductance
p= log((1 -PLCs)/PLCs)/(PsiL00 - PsiL50);%% [1/MPa]
q=-p*PsiL50; %%[-]
PLC = 1./(1+exp(p.*Psi_l+q)); %% [fraction]
PLC(PLC>1)=1; PLC(PLC<0)=0;
Kleaf = (1-PLC).*Kleaf_max ;%%% %%%% [mmolH20/MPa s m^2 leaf ]
Kleaf = Kleaf*LAI; %% [mmolH20/ MPa s m^2 PFT ]
%%%% Stem - Conductance
p=2/(PsiXe - PsiX50); %% [1/MPa]
q=-p*PsiX50; %%[-]
PLC = 1./(1+exp((p*Psi_x)+q)); %% [fraction]
PLC(PLC>1)=1; PLC(PLC<0)=0;
Kx = (1-PLC).*Kx_max; %%% %%%% [mmolH20/ m MPa s ]  Conductivity specific
Kx = Kx./hc*Axyl ; %% [mmolH20 /m^2 PFT s MPa] Conductance
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if PLD == 0
    if Axyl > 0
        %%%%% Fluxes
        Jxl = T + (Vl-Vltm1)/dth ; %%  [mm m2 /m^2 PFT h]
        Jsx = Jxl +(Vx-Vxtm1)/dth ; %%  [mm m2 /m^2 PFT h]
    else
        Jsl = T + (Vl-Vltm1)/dth ; %%  [mm m2 /m^2 PFT h]
        Jsx = Jsl;
        Jxl = Jsl;
    end
else
    Jxl = T; 
    Jsx = T;
end
Jxl=Jxl*Ccrown; %%% [mm/h]
Jsx=Jsx*Ccrown; %%% [mm/h]
return