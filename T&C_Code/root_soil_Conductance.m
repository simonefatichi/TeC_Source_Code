%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Root to Soil Conductivity                 %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Ksr]= root_soil_Conductance(Ks,Rl,rcyl,rroot,Zr)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Rl = 5300 ; %%% root length index [m root / m^2 ground]
%rcyl = 4.25*1e-3 ;%%%  [m] radius cylinder of soil to which root has access to
%rroot = 1*1e-3 ;%% [0.5-6 *10^-3] [m] root radius
Zr = Zr/1000 ; %% [m] Rooting depth
row= 1000; %% [kg/m^3] water density
g= 9.81; %% [m/s^2] gravity acceleration
rho= 55555; %% [mmolH20 / kg] Water density
%%%
CF = 10^6/(row*g); %% [m/MPa]
%%%%%%%%%%%%
%%%% Limit Hydraulic Conductivtiy to a minimum Value 
%Ks(Ks< 1e-06) = 1e-06; %% [mm/h] 
%%%%
%%% PsiS = [MPa] Soil Water Potential
% Ks [mm/h]
Ks=Ks/3600000; %% [m/s]
%Ks=Ks*row*rho; %%[mmol H20 /m^2 ground s]
%%%%%%%%%%%% Soil to Root conductance  --> Ksr
OPT=3; 
switch OPT
    case 1
        %%%%%%%%%%%% ; Daly et al 2004 ; Manzoni et al 2014; Katul et al 2003; 
        %%% Assumption grs = Ks/Lrs ;
        %RAI = Rl*rroot ; %%   Root Area Index  [m^2 root / m^2 ground]
        %Lrs = sqrt(2*Zr/Rl); %% [m]
        %gsr = Ks*sqrt(RAI./(2*rroot*Zr));  %%% [1/s] %%
        gsr = Ks*sqrt(Rl./(2*Zr));  %%% [1/s] %%
        gsr = gsr*CF; %% [m / s MPa]
        Ksr = gsr*row*rho;%% [mmol H20 / m^2 ground s MPa]
    case 2
        %%%% Janott et al., 2011
        Kr = 5*1e-08./CF ; % Radial conductivity or root [1/s]
        Lrs = rcyl; % radial thickenss of the rhizosphere [m]
        gsr =sqrt(Ks*Kr/Lrs); %% [1/s]
        Ksr = gsr*CF*row*rho ; % [mmol H20 / m^2 ground s MPa]
    case 3
        %%%%%%%%%%%% Holtta et al 2009  Sperry et al., 1998
        %gsr = Ks*Rl*(2*pi)/(log(rcyl/rroot)); %%[1/s]
        %%%%%%%%%%%%%%% %%%%%%%%%%
        %%%%%%%%%%  (Deckyman et al 2008; Newman 1969 -- ANAFORE model
        gsr =  Ks*Rl*(2*pi)/(log(rcyl/rroot)); %%[1/s]
        %%%%%%%%%%%%%%%%%%%%%%%%%
        Ksr = gsr*CF*row*rho ; % [mmol H20 / m^2 ground s MPa]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
end
%%% RADIAL --> Flow --> Krad*dP*row = [kg/m2 root s] or Krad*dP*row*RAI [kg/m2 ground s]
%Krad = 5*10^-8; %% % [1/s] %%  radial conductivity root Amenu and Kumar 2008
%Krad = 15*1e-8 ;%% [m /Pa s] radial conductivity root  -- North, G.B., Peterson, C.A. (2005) in Vascular transport in plants, Water
%%% Krad = 10^-9 - 7*10^-7 ;%% [m /MPa s] radial conductivity root Mendel et a 2002 WRR,   Huang and Nobel, 1994
%Krad = 0.3 - 20 *10^-8 ;%% [m /MPa s] radial conductivity root  Steudle and Peterson 1998
%Krad = 5 - 13 *10^-8 ;%% [m /MPa s] radial conductivity root  Bramley et al 2009 
%Krad = 2*10^-11 - 10^-9; %% % [1/s] %%  radial conductivity Schneider et al 2010
%Krad = 2*10^-9; %% % [1/s] %%  Javaux et al 2010
%Krad = 2*10^-7 -- 2*10^-5 [m /Mpa s] %% Draye et  al 2010 
%Krad= 0.5--2*10^-7  [m /Mpa s] %% Doussan et  al 2006
%Krad= 10^-9--10^-7  [m /Mpa s] %% Doussan et  al 1998 

%%%  AXIAL --> Flow = Kax/dL*dP*row ;; [kg / s]   
% Kax/dL*dP*row/(rroot*dL) ;; [kg/m^2 root /s] 
% Kax/dL*dP*row/(rroot*dL)*RAI ;; [kg/m^2 ground /s] 
%%% Kax = 0.2 ; % mm2/s   Root Axial  % Amenu and Kumar 2008
%%% Kax = 5*10^-11 - 4.2*10^-10 ;%% [m4 /MPa s] axial conductivity root Mendel et a 2002 WRR,   
%%% Kax = 2-6*10^-9 ;%% [m3 /MPa s] Bramley et al 2009 
%%% Kax = 2*10^-12 - 5*10^-9  ;%% [m4 /MPa s] Pierret et al 2006  
%%% Kax = 1*10^-12 - 1*10^-9  ;%% [m3 / s] Schneider et al 2010
%%%% Kax =5^10^-13-5*10^-12; %% % [m3 /s] %%  Javaux et al 2010
%Kax= 5*10^-11 -- 5*10^-9 [m4 /Mpa s] %% Draye et  al 2010 
%Kax= 5*10^-11 -- 1*10^-8 [m4 /Mpa s] %% Doussan et  al 2006 


  
 