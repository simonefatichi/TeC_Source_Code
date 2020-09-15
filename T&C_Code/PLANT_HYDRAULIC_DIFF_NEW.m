%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction - Plant Hydraulic           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%REFERENCES
function [dX]=PLANT_HYDRAULIC_DIFF_NEW(t,X,Ccrown,Psi_s,T,hc,LAI,Axyl,PsiL50,PsiL00,Kleaf_max,Cl,Sl,Kx_max,PsiX50,Cx,gsr,dt,V_avail)
%%% INPUTS
dX=zeros(2,1);
Vx=X(1);  Vl=X(2);
%%%%%%%%%%%%%%
%%%%%%%%%%%%%
% T Transpiration [mm/h]
% hc Canopy height [m]
% LAI  Leaf Area [m2 LAI / m2 PFT ]
% dt [s]
% Ccrown [m2 PFT / m2 Ground ]
% Psi_s [MPa] Soil water potential felt by plant
% Vxtm1 % [mmolH20 /m^2 PFT ]  Volume water xylem
% Vltm1 % [mmolH20 /m^2 PFT ] Volume water leaves
% Psi_ltm1 % [MPa]  Water potential in the leaves
% Psi_xtm1 % [MPa] Water potential in the xylem
% gsr  % soil to root conductance [mmol H20 / m^2 PFT s MPa]
%%% OUTPUTS
% Vx % [mmolH20 /m^2 PFT ]  Volume water xylem
% Vl % [mmolH20 /m^2 PFT ] Volume water leaves
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
V_avail = V_avail./Ccrown/(1000/row);% % [kg/m^2 PFT ]
V_avail = V_avail/row*rho2; %% %%  [mmolH20/m^2 PFT] 
%%%%%
T = T./Ccrown/(1000*3600/row);%% Transpiration Canopy  %% [kg/m^2 PFT s]
T = T/row*rho2; %%  [mmolH20/m^2 PFT s]
%%%
%Axyl = 15.0 ; %% [cm^2 stem /m^2 PFT]
%Axyl =  Axyl/10000; %% [m2 stem/m2 PFT]
%%% Leaf
%PsiL50 =  -0.7 ;%%[MPa]  Water Potential at 50% loss conductivity
%PsiLs =  -1.6; %% [MPa]  Water Potential at PLCs% loss conductivity
PLCs = 0.02;  %%[-]  specific loss of conductivity
%Kleaf_max = 10 ; %%  %%%  [mmolH20 m^2 leaf s /MPa]
%Cl  = 1;  %%%  %  [mmolH20 / m^2 leaf MPa]
%Cl = Cl*LAI ; %% [mmolH20 / m^2 PFT MPa]
%%%%%%%%%%
%%% gsr ; %% [mmol H20 / m^2 PFT s MPa]
%%%%%%% 
%%% Xylem
%Kx_max = 80000;  %%5550-555550 [mmolH20 /m s MPa]  Xylem Conductivity specific for water;
%PsiX50 = -3.5; %%[MPa]  Water Potential at 50% loss conductivity
PsiXe = PsiX50/1.42; %% [MPa] (Typically at 12% of Xylem, Meinzer et al., 2010)    Meinzer et al., 2009
%Cx = 150; %%% [kg / m^3 sapwood MPa]
%Cx = Cx/row*rho2*hc*Axyl; %%  % [mmolH20 / m^2 PFT MPa]
%%%%%%% 

% Plant Pressure Volume curve
[Psi_x] = Plant_PV_Curve(Vx,Cx,'X',LAI,Axyl,hc,Sl); 
[Psi_l] = Plant_PV_Curve(Vl,Cl,'L',LAI,Axyl,hc,Sl); 

%%% Leaf conductance
p= log((1 -PLCs)/PLCs)/(PsiL00 - PsiL50);%% [1/MPa]
q=-p*PsiL50; %%[-]
PLC = 1./(1+exp(p.*Psi_l+q)); %% [fraction]
PLC(PLC>1)=1; PLC(PLC<0)=0; 
Kleaf = (1-PLC).*Kleaf_max ;%%% %%%% [mmolH20/MPa s m^2 leaf ]
Kleaf = Kleaf*LAI; %% [mmolH20/ MPa s m^2 PFT ]
%Cl = Cl*(PLC<0.90);
%%%% Stem - Conductance
p=2/(PsiXe - PsiX50); %% [1/MPa]
q=-p*PsiX50; %%[-]
PLC = 1./(1+exp((p*Psi_x)+q)); %% [fraction]
PLC(PLC>1)=1; PLC(PLC<0)=0;
Kx = (1-PLC).*Kx_max; %%% %%%% [mmolH20/ m MPa s ]  Conductivity specific
Kx = Kx./hc*Axyl ; %% [mmolH20 /m^2 PFT s MPa] Conductance
%%%%%%%%
%Cx = Cx*(1-PLC); 
%%%%%%
%%%%%%% 
if Axyl > 0
    %%%%% Fluxes
    KM = 2*(Kx*gsr)/(Kx+gsr);
    Jsx = min(V_avail/dt,(KM)*(Psi_s - Psi_x)); %% [mmolH20 /m^2 PFT s ]
    Jsx(Jsx<0) = 0;
    if LAI == 0 
        %Jxl = 0;
        %dX(1) = (Jsx-Jxl)/Cx; %  [MPa/s]
        dX(1) = (Jsx - Jxl);  
        dX(2) = (Jxl - T); %% [mmolH20/ m2 PFT s];
    else
        KM = 2*(Kx*Kleaf)/(Kx+Kleaf);
        Jxl =KM*(Psi_x - Psi_l); %% [mmolH20/m^2 PFT  s]
        %Jsx = ((2*Kx*gsr)/(2*Kx + gsr))*(Psi_s - Psi_xtm1); %% [mmolH20 /m^2 PFT s ]
        %Jxl = ((2*Kx*Kleaf)/(2*Kx+Kleaf))*(Psi_xtm1 - Psi_ltm1); %% [mmolH20 /m^2 PFT  s]
        %%%%
        %dX(1) = (Jsx-Jxl)/Cx; % [MPa / s ]
        %dX(2) = (Jxl - T)/Cl;%  [MPa / s ]
        %Psi_x = Psi_xtm1 + dt*(Jsx - Jxl)/Cx;
        %Psi_l = Psi_ltm1 + dt*(Jxl - T)/Cl;
        %%%%%%%%%%%%
        %Vx = Vxtm1 + dt*(Jsx - Jxl); %% 
        %Vl = Vltm1 + dt*(Jxl - T); %% 
        dX(1) = (Jsx - Jxl); % %  [mmolH20/m^2 PFT s]
        dX(2) = (Jxl - T); %% [mmolH20/ m2 PFT s];
    end
else   
    if LAI == 0 
        Jsl=0;
        dX(2) = (Jsl - T); %%% [mmolH20/ m2 PFT s];
    else
        %Jsl = ((2*Kleaf*gsr)/(2*Kleaf + gsr))*(Psi_s - Psi_ltm1); %% [mmolH20 /m^2 PFT s ]
        %Psi_l = Psi_ltm1 + dt*(Jsl - T)/Cl;
        %Vl = Vltm1 + dt*(Jsl - T)/rho2; %% [m3/ m2 PFT];
        KM=2*(Kleaf*gsr)/(Kleaf + gsr); 
        Jsl = min(V_avail/dt,(KM)*(Psi_s - Psi_l)); %% [mmolH20 /m^2 PFT s ]
        %dX(2) = (Jsl - T)/Cl; %  [MPa/s]
        dX(2) = (Jsl - T); %% [mmolH20/ m2 PFT s];
        %%%
        %Jsx = Jsl; Jxl = Jsl;
        %Psi_x = Psi_l;
        %Vx = 0;
    end
end
%%%%%%%%%%%%%%%
%%%%%%--> Death of Plant 
if Psi_x < PsiX50
    dX=zeros(2,1);
end 
return