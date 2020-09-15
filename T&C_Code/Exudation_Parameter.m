%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Exudation_Parameter     %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[ParEx]=Exudation_Parameter(bfix)
%%%%%%%%%%%%%%%%
%%%%%%%% Default Options
ParEx.dexmy=0.04; %%%   %% [-]  Root exudation fraction of NPP 
ParEx.OPT_EXU=3; %% Method to compute Root exudation (1) constant fraction NPP (2) Cost evaluation of Nutrient uptake - FUN2.0 (3) Cost evaluation based on Uptake Physics  
ParEx.bfix=bfix; %%% Biological fixers (1/0)
%%%%%
if ParEx.OPT_EXU == 2
    ParEx.kc1 = 10; ParEx.kc2 = 300; ParEx.kc3 = 50 ;
    ParEx.kn1 = 600; ParEx.kn2 = 50; ParEx.kn3 = 100;
    ParEx.kp1 = ParEx.kn1/8 ; ParEx.kp2 = ParEx.kn2/8;  ParEx.kp3 = ParEx.kn3/8;
    ParEx.kk1 = ParEx.kn1/2.5; ParEx.kk2 = ParEx.kn2/2.5; ParEx.kk3 = ParEx.kn3/2.5;
end
if ParEx.OPT_EXU == 3
    ParEx.kc1 = NaN; ParEx.kc2 = NaN; ParEx.kc3 =NaN ;
    ParEx.kn1 = 7000; ParEx.kn2 = 3.50; ParEx.kn3 = 2.00;
    ParEx.kp1 = ParEx.kn1/15 ; ParEx.kp2 = ParEx.kn2/15;  ParEx.kp3 = ParEx.kn3/15;
    ParEx.kk1 = ParEx.kn1/2.5; ParEx.kk2 = ParEx.kn2/2.5; ParEx.kk3 = ParEx.kn3/2.5;
end
