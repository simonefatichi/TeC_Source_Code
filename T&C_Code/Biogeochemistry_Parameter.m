%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Biogeochemistry_Parameter     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[BiogeoPar]=Biogeochemistry_Parameter(opt_cons_CUE)
%%%%%%%%%%%%%%
%%%% Decay Constant on a 40°C basis (Kirschbaum 2002)
BiogeoPar.k_met_sur = 1/12.5; %%[1/day]
BiogeoPar.k_str_sur = 1/46; %% %[1/day]
BiogeoPar.k_met_ssr =  1/10;
BiogeoPar.k_str_ssr = 1/37;
BiogeoPar.k_wod_sur = 1/150; %% [1/92 - 1/460]
%%%%%%%%%%%%%%%%%%%%
%%%%%% Respiration Coefficient for Litter [-]  Parton et al 1987 
%rr_met_sur=0.55;
%rr_str_sur=0.45;
%rr_wod_sur=0.55;
%rr_met_ssr=0.55;
%rr_str_ssr=0.55;
%%%%%%% Critical Ratios
BiogeoPar.R_CN_bac = 5.2; % 7 ;%  %%[C/N]  8.6
BiogeoPar.R_CN_fun = 6.5; % 12 ;% 
BiogeoPar.R_CN_myc = 18;
BiogeoPar.R_CN_ew = 10; 
BiogeoPar.R_CP_bac = 16; % 35; %%[C/P] 60
BiogeoPar.R_CP_fun = 40;% 80; %% 120 
BiogeoPar.R_CP_myc = 120; % 120;
%%%%%%%%%%%%%%%%%%
%%% Parameters  %%% Wang et al 2012  Waring-Averill Ecol. Let 2013 
BiogeoPar.gd =0.2;% 0.3;% 0.3;% [0.2-0.8] fraction of dead B allocated to D
BiogeoPar.fd =0.4;% 0.5; % [0.2-0.8]  fraction of decomposed P allocated to D
%Ecm =0.47; % carbon use efficiency [-]
%Ecf =0.47;
BiogeoPar.Vdb = 0.04; %0.120; %0.010;  % [0.0024-0.060] maximum specific uptake rate of DOC for growth of B [gC-D/gC-B day]
BiogeoPar.Vdf = 0.02; %0.025;% 0.012;  % 
BiogeoPar.Kdb = 0.00026; % [0.00014 - 0.00038] half-saturation constant of uptake of DOC for growth of B [ g C/g soil]
BiogeoPar.Kdf = 0.00026;
BiogeoPar.mrb =  5.0*1e-3; % 15.0*1e-3; % [ 1.20*1e-3 - 30*1e-3] specific maintenance factor or rate of B [ g C /g C day]
BiogeoPar.mrf =  2.0*1e-3; % 5.0*1e-3;
BiogeoPar.mram = 1.2*1e-3; % 1.5*1e-3; specific maintenance factor of AM 
BiogeoPar.mrem = 1.2*1e-3; % 1.5*1e-3; specific maintenance factor of EM 
BiogeoPar.Vpc =200;% 95; % 120; %% [4.8 - 792] maximum specific decomposition rate for POC-C by EP [gC-P/gC-EP day]
BiogeoPar.Vpl =23; % 8; % 23;% 
BiogeoPar.Vm = 100; % 48;% 24;  %[1.2 -528] maximum specific decomposition rate for MOC by [gC-M /gC-EM day]
BiogeoPar.Kpl =0.05; % [0.023 - 0.094] half-saturation constant for decomposition of POC-L [gC / gsoil]
BiogeoPar.Kpc =0.05; 
BiogeoPar.Km =0.25; % [0.057 - 0.95] half-saturation constant for decomposition of MOC [ gC /gsoil];
BiogeoPar.rem = 0.018; % 0.025; % [0.0075 - 0.075] turnover rate of EM [g C/g C day]
BiogeoPar.rep = 0.018; % 0.025; %  [0.0075 - 0.075] turnover rate of EP [g C/g C day]
BiogeoPar.pepb = 0.012; % [0.0031 - 0.031] fraction of mR for production of EP
BiogeoPar.pemb = 0.005;% [0.0031 - 0.031] fraction of mR for production of EM
BiogeoPar.pepf = 0.006; % [0.0031 - 0.031] fraction of mR for production of EP
BiogeoPar.pemf = 0.010;% [0.0031 - 0.031] fraction of mR for production of EM
BiogeoPar.pepem = 0.006;
BiogeoPar.pemem = 0.010; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
BiogeoPar.Ecf = 0.27;% 0.47; %%% Sinsbaugh et al 2016 
BiogeoPar.Ecb = 0.27;% 0.47;
%%% Sinsabaguh et al 2013 ; Manzoni et al 2008
BiogeoPar.CUEmax= 0.6;
BiogeoPar.Ae=1; %% (not used)

%%% Temperature effects on CUE
%Tref=20;
BiogeoPar.T1= 12;  %%%%%% Temperature Effects Kinetiks
BiogeoPar.T2= 20;  %%% Temperature effects on CUE
BiogeoPar.T3= 15; % Temperature Earthworms
BiogeoPar.T4= 40; %%% Ref. Temp. Litter Decoomposition
%%%%%%%%
%%%%%
BiogeoPar.mtem = -0.008; %%[°C] [0 - 0.016] Allison et al 2010; Li et al 2014
if opt_cons_CUE == 1
    BiogeoPar.rr_met_sur = 0.55;
    BiogeoPar.rr_str_sur=  0.45;
    BiogeoPar.rr_wod_sur=  0.45;
    BiogeoPar.rr_met_ssr= 0.55;
    BiogeoPar.rr_str_ssr= 0.55;
end
%%%%%%%%%%%%%% Leaching Coefficients 
BiogeoPar.lambda_C=0.0015;  %%% 0 up to 2500 mm/yr 0.5 at 4000 mm/yr linear between 
BiogeoPar.lambda_N=0.0015;
BiogeoPar.lambda_P=0.0001;
BiogeoPar.lambda_K=0.90;
BiogeoPar.f_org_lea=1; %% Fraction of Organic leaching 
%%%%%%%%%%%%%%%%%% Earthworm parameters 
%BiogeoPar.kse =  5e-004 ; % [gC/ gsoil] half-saturation constant for consumption of food %
%BiogeoPar.VmaxE = 0.038 ; %[gC / gC-EW d] Ref earthworm growth 15°C
BiogeoPar.VmaxE = 0.0050 ; %%  [gC / gC-food d]
BiogeoPar.rmanEa= 0.0300 ; %% [gC / gC-EW d] Maint. respiration 15°C 
BiogeoPar.rmanEp= 0.0146 ; %% [gC / gC-EW d] Maint. respiration 15°C 
BiogeoPar.dew = 0.010; % %% [gC / gC-EW d] 
BiogeoPar.fabsE =0.15; % [-] fraction of absorbed food 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Volatilization of Ammonia NH3
BiogeoPar.kv = 8.6*1e-5; %% [1/d]
%%%%%%% Nitrification (Aerobic)
BiogeoPar.KNmaxn = 0.086; %% [1/d] % Dickinson et al., 2002
%%%%%%% Denitrification (Anaerobic)
BiogeoPar.Kc = 17.5; %% [gC/ gsoil] %% Li et al., 1992 (not used)
BiogeoPar.Kn = 82.5; %% [gN/ gsoil] % Li et al., 1992 (not used)
BiogeoPar.KNmaxd= 0.217; %% [1/d] % Dickinson et al., 2002
%%%%%% Parton et al 1997 
%%%%%%
BiogeoPar.K1p= 1/600; %%% [1/day]
BiogeoPar.K2p= 1/13500; %% [1/day]
BiogeoPar.K3p= 1/(30*1e+6); %% [1/day]
BiogeoPar.K4p= 1/(4.38*1e+6);% 1/300000; %% [1/day]
%%%%%
%%%%%% Myself 
BiogeoPar.Kdk=1/2;%%%% [1/day] desorption exchangeable
BiogeoPar.Kak= 1/2; %% [1/day] absorption exchangeable
%%%%
BiogeoPar.K1k= 1/600 ; % []; %%% [1/day] Fixation of  mineral K
BiogeoPar.K2k= 1/13500; % []; %% [1/day]  Release of fixed K 
BiogeoPar.K3k= 1/350000; %% [1/day] Release primary K 
%%%%%%%%%%%%%%%%%%
%%%%%%%%%%% Solubility coefficients for leakage and passive uptake 
BiogeoPar.aNH4 =0.05;  %%0.1 CLM4
BiogeoPar.aNO3 = 1;  
BiogeoPar.aP = 0.005;
BiogeoPar.aK = 1; 
BiogeoPar.aDOC = 0.04; 
BiogeoPar.aDON= 1;
BiogeoPar.aDOP= 1;
%%%%%%%%%%%%%%% Uptake coefficients 
%%% Root EM AM  
%%%--< These are dependent on the exudation/mycorrhiza allocations costs 
% Cex_root = 5; 
% kc1 = 10; kc2 = 300; kc3 = 50 ;
% kn1 = 600; kn2 = 50; kn3 = 100;
%kp1 = kn1/8 ; kp2 = kn2/8;  kp3 = kn3/8;
%kk1 = kn1/2.5; kk2 = kn2/2.5; kk3 = kn3/2.5;
%%%%%%%
% BiogeoPar.vem  = BiogeoPar.mrm/kc2; %  5e-006 ;  %%[gN /gC^2 day]
% BiogeoPar.Kvem = (kn2/kc2/365); % 45.6e-005 ; %%[gN /gC ] 
% BiogeoPar.vam =  BiogeoPar.mrm/kc3;% 3e-005 ;  %%[gN /gC^2 day]
% BiogeoPar.Kvam = (kn3/kc3/365);% 0.0055 ; %%[gN /gC ] 
% BiogeoPar.vr = (1/kc1)*Cex_root/365;% 0.0014; %%[gN /gC day]
% BiogeoPar.Kvr = (kn1/kc1/365);% 0.1644 ; %%[gN /gC ] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kn1 = 7000; kn2 = 3.5; kn3 = 2.0;
Zr = 0.25; % [m]
rroot =  0.5*1e-3 ; % 3.3*1e-4 ;%% [0.5-6 *10^-3] [m] root radius
rmyc_am =  5*1e-6 ; % [m] hyphae radius
rmyc_em =  5*1e-6 ; % [m] hyphae radius
rho_r = 122*1000; %% [gC / m^3] Root density  Jackson et al., 1997
rho_m = 200*1000; %% [gC / m^3]
BiogeoPar.g_r = (sqrt(2)*(rroot.^2)*(Zr.^1.5)*(rho_r.^1.5));
BiogeoPar.g_em = (sqrt(2)*(rmyc_em.^2)*(Zr.^1.5)*(rho_m.^1.5));
BiogeoPar.g_am = (sqrt(2)*(rmyc_am.^2)*(Zr.^1.5)*(rho_m.^1.5));
BiogeoPar.vr =   2.3e-7; %(Cex_root/365)*BiogeoPar.g_r/kn1 ;%% [m^2/day]
BiogeoPar.vem =  2.3e-7; %BiogeoPar.g_em*BiogeoPar.mrm/kn2 ;%% [m^2/day]
BiogeoPar.vam =  2.3e-7; %BiogeoPar.g_am*BiogeoPar.mrm/kn3 ; %% [m^2/day]

%%%%%% Uptake coefficients 
%BiogeoPar.vNH4 = 1.5*1e-4  ;  %%[gN /gC day]  2.33*1e-3  5.14*1e-6 
%BiogeoPar.vNO3 =  1.5*1e-4 ;  
%BiogeoPar.vP =  1.5*1e-4 ;
%BiogeoPar.vK =  1.5*1e-4 ;
%%%%%%%%%% half-constants  for uptakes 
%BiogeoPar.KnNH4 = 0.5; %% [0.83-1.0]  [gN/m2] 
%BiogeoPar.KnNO3 = 0.2; 
%BiogeoPar.KnP=0.05;
%BiogeoPar.KnK=0.1; 
return 


