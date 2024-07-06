%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% OUTPUT WRITING   %%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Mass_balance_Variable
if t==2
    Qlat_in_tgtm1 = 0;
    q_runon_tgtm1 = 0;
    Q_channel_tgtm1 = 0;
else
    V_tgtm1 = V_tg;
    Vice_tgtm1 = Vice_tg; 
    SWE_tgtm1 = SWE_tg;
    In_tgtm1 =In_tg;
    ICE_tgtm1 =  ICE_tg; %%
    WAT_tgtm1 =  WAT_tg; %%
    FROCK_tgtm1 = FROCK_tg; %%
    %%%%---
    Qlat_in_tgtm1 = Qlat_in_tg;
    q_runon_tgtm1 = q_runon_tg;
    Q_channel_tgtm1 = Q_channel_tg;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SPATIAL AVERAGE over the watershed
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Pr_tg = mean(Pr_S(Kinde)); %%% [mm]
Ta_tg = mean(Ta_S(Kinde)); %%% [°C]
Ws_tg = mean(Ws_S(Kinde)); %%% [m/s]
Ds_tg = mean(Ds_S(Kinde)); %%% [°C]
ea_tg = mean(ea_S(Kinde)); %%% [Pa]
N_tg = mean(N_S(Kinde)); %%% [-]
Pre_tg = mean(Pre_S(Kinde)); %%% [mbar]
Tdew_tg = mean(Tdew_S(Kinde)); %%% [°C]
Rsw_space = SAB1_S+ SAB2_S + SAD1_S+ SAD2_S;
Rsw_tg = mean(Rsw_space(Kinde)); %%% [W/m^2]
PAR_space = PARB_S + PARD_S; %%% %% [W/m^2]
PAR_tg= mean(PAR_space(Kinde));%% [W/m^2]
Ca_tg = mean(Ca_S(Kinde)); %% [ppm]
%%%%%
Ts_tg  =  mean(Ts(Kinde));%% [°C]
Tdamp_tg  = mean(Tdamp(Kinde));%% [°C]
Csno_tg = mean(Csno(Kinde)); %%[]
Cice_tg = mean(Cice(Kinde)); %%[]
Csnow_tg = mean(Csnow(Kinde)); %%[]
Cicew_tg = mean(Cicew(Kinde)); %%[]
Pr_sno_tg  =  mean(Pr_sno(Kinde));%%[mm]
Pr_liq_tg = mean(Pr_liq(Kinde)); %%[mm]
%%%%
Rn_tg=mean(Rn(Kinde));%% [W/m^2]
H_tg=mean(H(Kinde));%% [W/m^2]
G_tg=mean(G(Kinde));%% [W/m^2]
Gfin_tg=mean(Gfin(Kinde)); % [W/m^2]
QE_tg=mean(QE(Kinde));%% [W/m^2]
Qv_tg = mean(Qv(Kinde)); %% [W/m^2]
Qfm_tg =  mean(Qfm(Kinde));% [W/m^2]
%%%%
SWE_tg = mean(SWE(Kinde)); %% [mm]
SND_tg = mean(SND(Kinde)); %% [m]
WR_SP_tg = mean(WR_SP(Kinde)); %%[mm]
U_SWE_tg = mean(U_SWE(Kinde)); %%[mm]
NIn_SWE_tg  =  mean(NIn_SWE(Kinde));%% [mm]
dw_SNO_tg=mean(dw_SNO(Kinde)); %%[]
ros_tg = mean(ros(Kinde)); %% [kg/m^3]
In_SWE_tg =  mean(In_SWE(Kinde)); %%[mm]
SP_wc_tg =  mean(SP_wc(Kinde));%%[mm]
%%%%
ICE_tg = mean(ICE(Kinde)); %% [mm]
ICE_D_tg = mean(ICE_D(Kinde)); %% [m]
WR_IP_tg = mean(WR_IP(Kinde)); %%[mm]
IP_wc_tg =  mean(IP_wc(Kinde));%%[mm]
NIce_tg  =  mean(NIce(Kinde));%% [mm]
%%%
Imelt_tg = mean(Imelt(Kinde)); %% [mm]
Smelt_tg = mean(Smelt(Kinde)); %% [mm]
Tice_tg =  mean(Tice(Kinde)); %% [C]
%%%%%%%%%
T_H_space = sum(T_H,2);
T_L_space = sum(T_L,2);
T_H_tg = mean(T_H_space(Kinde)); %%%[mm/h]
T_L_tg = mean(T_L_space(Kinde)); %%%[mm/h]
T_tg= T_L_tg +T_H_tg  ; %%% [mm/h]
EIn_H_space = sum(EIn_H,2);
EIn_L_space = sum(EIn_L,2);
EIn_L_tg  = mean(EIn_L_space(Kinde)); %%% [mm/h]
EIn_H_tg =  mean(EIn_H_space(Kinde)); %%% [mm/h]
EIn_tg =EIn_H_tg+ EIn_L_tg; %%% [mm/h]
EG_tg =mean(EG(Kinde)); %%% [mm/h]
ESN_tg =mean(ESN(Kinde)+ESN_In(Kinde)); %%% [mm/h]
EWAT_tg = mean(EWAT(Kinde)); %%[mm/h]
EICE_tg = mean(EICE(Kinde)); %%[mm/h]
Dr_H_space = sum(Dr_H,2);
Dr_L_space = sum(Dr_L,2);
Dr_L_tg  = mean(Dr_L_space(Kinde)); %%% [mm]
Dr_H_tg =  mean(Dr_H_space(Kinde)); %%% [mm]
%%%%%%%%%%%%%%%%
EIn_urb_tg = mean(EIn_urb(Kinde)); %%[mm/h]
EIn_rock_tg = mean(EIn_rock(Kinde)); %%[mm/h]
SE_rock_tg = mean(SE_rock(Kinde)); %%[]
SE_urb_tg = mean(SE_urb(Kinde)); %%[]
%%%
In_H_space = sum(In_H,2);
In_L_space = sum(In_L,2);
In_urb_tg =  mean(In_urb(Kinde)); %%[mm]
In_rock_tg =  mean(In_rock(Kinde)); %%[mm]
In_tg = mean(In_H_space(Kinde) +  In_L_space(Kinde) +  SP_wc(Kinde) + In_SWE(Kinde) + In_urb(Kinde) + In_rock(Kinde)+ IP_wc(Kinde) ); %%% [mm]
Inveg_tg = mean(In_H_space(Kinde) +  In_L_space(Kinde));
%%%
WAT_tg= mean(WAT(Kinde)) ;%%%[mm]
FROCK_tg = mean(FROCK(Kinde)); %%%[mm]
%%%
f_tg = mean(f(Kinde)*dth); %%[mm]
WIS_tg = mean(WIS(Kinde)); %%[mm]
Rd_tg =mean(Rd(Kinde)); %% [mm]
Rh_tg =mean(Rh(Kinde)); %%[mm]
Lk_tg =mean(Lk(Kinde)*dth); %% [mm]
Lk_wat_tg =mean(Lk_wat(Kinde)*dth); %% [mm]
Lk_rock_tg =mean(Lk_rock(Kinde)*dth); %% [mm]
%%%
OF_tg = mean(OF(Kinde)); %%[]
OS_tg = mean(OS(Kinde)); %%[]
ZWT_tg= mean(ZWT(Kinde)); %% [mm]
Fract_sat_tg = sum(ZWT(Kinde)==0)/num_cell;%% [-]
%%%%%%%%
Qi_in_space = sum(Qi_in*dth,2);
Qi_out_space = sum(Qi_out*dth,2);
Qlat_in_tg=  mean(Qi_in_space(Kinde)); %% [mm]
Qlat_out_tg =  mean(Qi_out_space(Kinde)); %% [mm]
q_runon_tg = mean(q_runon(Kinde)*dth); %% [mm]
Q_channel_tg =  mean(Q_channel(Kinde)); %% [mm]
%%%%%%%%%
V_space =  sum(V,2); %%%
Vice_space =  sum(Vice,2); %%%
O_space =  V_space./Zs_OUT + Ohy_OUT ; %%[-]
V_tg = mean(Asur(Kinde).*V_space(Kinde)); %% [mm]
Vice_tg = mean(Asur(Kinde).*Vice_space(Kinde)); %% [mm]
O_tg = mean(O_space(Kinde)); %% [-]
%%%%%%
ra_tg  =  mean(ra(Kinde));%% [s/m]
r_soil_tg  =  mean(r_soil(Kinde));%% [s/m]
alp_soil_tg = mean(alp_soil(Kinde)); %%[-]
%%%%
Tdp_space =  mean(Tdp,2); %%%
Tdp_tg = mean(Tdp_space(Kinde)); %% % [°C]
Tdpsnow_space =  mean(Tdpsnow,2); %%%
Tdpsnow_tg = mean(Tdpsnow_space(Kinde)); %% % [°C]
%%%%
er_tg =  mean(er(Kinde)); %%  % [kg/s m^2]
%%%
TsVEG_tg  =  mean(TsVEG(Kinde));%% [°C]
Ts_under_tg  =  mean(Ts_under(Kinde));%% [°C]
DQ_S_tg = mean(DQ_S(Kinde)); %%[]
DT_S_tg = mean(DT_S(Kinde)); %%[]
dQ_S_tg = mean(dQ_S(Kinde)); %%[]
CK1_tg=  mean(CK1(Kinde));%% [mm]
%CK2_tg =  mean(CK2(Kinde));%% [mm]
%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% SIMPLE MASS CHECK CONTROL
if sum(SNn)>=1
    CKt = (V_tgtm1 - V_tg) + (Vice_tgtm1 - Vice_tg) + Pr_tg - EG_tg*dth - T_tg*dth - EIn_tg*dth ...
        - ESN_tg*dth - EIn_urb_tg*dth - EWAT_tg*dth - EIn_rock_tg*dth - EICE_tg*dth ...
        + (SWE_tgtm1 -SWE_tg) + (In_tgtm1 -In_tg) ...
        +  (ICE_tgtm1 -ICE_tg) +  (WAT_tgtm1 -WAT_tg)  +  (FROCK_tgtm1 -FROCK_tg) ...
        + Qlat_in_tgtm1 + q_runon_tgtm1 + Q_channel_tgtm1  ...
        - Qlat_in_tg - Q_exit - Qsub_exit - q_runon_tg - Q_channel_tg -Swe_exit;
else
    CKt = (V_tgtm1 - V_tg) + (Vice_tgtm1 - Vice_tg) + Pr_tg - EG_tg*dth - T_tg*dth - EIn_tg*dth ...
        - ESN_tg*dth - EIn_urb_tg*dth - EWAT_tg*dth - EIn_rock_tg*dth - EICE_tg*dth ...
        + (SWE_tgtm1 -SWE_tg) + (In_tgtm1 -In_tg) ...
        +  (ICE_tgtm1 -ICE_tg) +  (WAT_tgtm1 -WAT_tg)  +  (FROCK_tgtm1 -FROCK_tg) ...
        + Qlat_in_tgtm1 + q_runon_tgtm1   ...
        - Qlat_in_tg - Q_exit - Qsub_exit - q_runon_tg  -Swe_exit;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Ckk1=  Pr_tg  + q_runon_tgtm1 - EIn_tg*dth - WIS_tg + (In_tgtm1 -In_tg); %%Not correct
%Ckk2 = WIS_tg - f_tg - Rh_tg;
%Ckk3 = Rh_tg - Q_exit - q_runon_tg;
%Ckk4 = f_tg + (V_tgtm1 - V_tg) - EG_tg*dth - T_tg*dth - Lk_tg ...
%  - Qlat_in_tg -Rd_tg + Qlat_in_tgtm1 - Qsub_exit ;
%Ckk5 = Pr_sno_tg*dth + Pr_liq_tg*dth*Csno_tg ...
%    - ESN_tg*dth + (SWE_tgtm1 -SWE_tg) ...
%    + (In_tgtm1 -In_tg) -sum(Dr_H_tg)*Csno_tg ; %% Not correct
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t==2
    tit{1}=strcat('OUTPUT_',TITLE_SAVE,'_AVG.dat');
    fid(1)=fopen(tit{1},'a');
end
fprintf(fid(1),'%g\t',Pr_tg); %1
fprintf(fid(1),'%g\t',Ta_tg); %2
fprintf(fid(1),'%g\t',Ws_tg); %3
fprintf(fid(1),'%g\t',Ds_tg); %4
fprintf(fid(1),'%g\t',ea_tg); %5
fprintf(fid(1),'%g\t',N_tg);  %6
fprintf(fid(1),'%g\t',Pre_tg); %7
fprintf(fid(1),'%g\t',Tdew_tg);%8
fprintf(fid(1),'%g\t',Rsw_tg); %9
fprintf(fid(1),'%g\t',PAR_tg);%10
fprintf(fid(1),'%g\t',Ca_tg);%%11
%%%%%
fprintf(fid(1),'%g\t',Ts_tg);%12
fprintf(fid(1),'%g\t',Tdamp_tg);%13
fprintf(fid(1),'%g\t',Csno_tg);%14
fprintf(fid(1),'%g\t',Cice_tg);%15
fprintf(fid(1),'%g\t',Csnow_tg);%16
fprintf(fid(1),'%g\t',Cicew_tg);%17
fprintf(fid(1),'%g\t',Pr_sno_tg);%18
fprintf(fid(1),'%g\t',Pr_liq_tg);%19
%%%%
fprintf(fid(1),'%g\t',Rn_tg);%20
fprintf(fid(1),'%g\t',H_tg);%21
fprintf(fid(1),'%g\t',G_tg);%22
fprintf(fid(1),'%g\t',Gfin_tg);%23
fprintf(fid(1),'%g\t',QE_tg);%24
fprintf(fid(1),'%g\t',Qv_tg);%25
fprintf(fid(1),'%g\t',Qfm_tg);%26
%%%%
fprintf(fid(1),'%g\t',SWE_tg);%27
fprintf(fid(1),'%g\t',SND_tg);%28
fprintf(fid(1),'%g\t',WR_SP_tg);%29
fprintf(fid(1),'%g\t',dw_SNO_tg);%30
fprintf(fid(1),'%g\t',ros_tg);%31
fprintf(fid(1),'%g\t',In_SWE_tg);%32
fprintf(fid(1),'%g\t',SP_wc_tg);%33
%%%%%%%%%
fprintf(fid(1),'%g\t',ICE_tg);%34
fprintf(fid(1),'%g\t',ICE_D_tg);%35
fprintf(fid(1),'%g\t',WR_IP_tg);%36
fprintf(fid(1),'%g\t',NIce_tg);%37
fprintf(fid(1),'%g\t',IP_wc_tg);%38
%%%%%%%%%
fprintf(fid(1),'%g\t',T_H_tg);%39
fprintf(fid(1),'%g\t',T_L_tg );%%40
fprintf(fid(1),'%g\t',EIn_H_tg);%%41
fprintf(fid(1),'%g\t',EIn_L_tg);%%42
fprintf(fid(1),'%g\t',EG_tg);%%43
fprintf(fid(1),'%g\t',ESN_tg);%%44
fprintf(fid(1),'%g\t',EWAT_tg);%%45
fprintf(fid(1),'%g\t',EICE_tg);%%46
fprintf(fid(1),'%g\t',Dr_L_tg);%47
fprintf(fid(1),'%g\t',Dr_H_tg);%%48
%%%%%%%%
fprintf(fid(1),'%g\t',EIn_urb_tg);%%49
fprintf(fid(1),'%g\t',EIn_rock_tg);%%50
fprintf(fid(1),'%g\t',In_tg);%%51
fprintf(fid(1),'%g\t',Inveg_tg);%%52
%%%
fprintf(fid(1),'%g\t',WAT_tg);%%53
fprintf(fid(1),'%g\t',FROCK_tg);%%54
%%%
fprintf(fid(1),'%g\t',f_tg);%%55
fprintf(fid(1),'%g\t',WIS_tg);%%56
fprintf(fid(1),'%g\t',Rd_tg);%%57
fprintf(fid(1),'%g\t',Rh_tg);%%58
fprintf(fid(1),'%g\t',Lk_tg);%%59
fprintf(fid(1),'%g\t',Lk_wat_tg);%%60
fprintf(fid(1),'%g\t',Lk_rock_tg);%%61
fprintf(fid(1),'%g\t',OF_tg);%%62
fprintf(fid(1),'%g\t',OS_tg);%%63
fprintf(fid(1),'%g\t',ZWT_tg);%%64
fprintf(fid(1),'%g\t',Fract_sat_tg);%%65
%%%
fprintf(fid(1),'%g\t',Qlat_in_tg);%%66
fprintf(fid(1),'%g\t',Qlat_out_tg);%%67
fprintf(fid(1),'%g\t',q_runon_tg);%%68
fprintf(fid(1),'%g\t',Q_channel_tg);%%69
fprintf(fid(1),'%g\t',V_tg);%%70
fprintf(fid(1),'%g\t',O_tg);%%71
fprintf(fid(1),'%g\t',Q_exit);%%72
fprintf(fid(1),'%g\t',Qsub_exit);%%73
fprintf(fid(1),'%g\t',Swe_exit);%%74
%%%%
fprintf(fid(1),'%g\t',r_soil_tg);%75
fprintf(fid(1),'%g\t',alp_soil_tg);%76
fprintf(fid(1),'%g\t',ra_tg );%77
fprintf(fid(1),'%g\t',Tdp_tg );%78
%%%%%%%%%%%%%%%%
fprintf(fid(1),'%g\t',er_tg);%%79
fprintf(fid(1),'%g\t',TsVEG_tg);%%80
fprintf(fid(1),'%g\t',Ts_under_tg);%%81
fprintf(fid(1),'%g\t',Tdpsnow_tg);%%82
fprintf(fid(1),'%g\t',DQ_S_tg);%%83
fprintf(fid(1),'%g\t',DT_S_tg);%%84
fprintf(fid(1),'%g\t',dQ_S_tg);%%85
%%%%%%%%%
fprintf(fid(1),'%g\t',Imelt_tg);%%86
fprintf(fid(1),'%g\t',Smelt_tg);%%87
fprintf(fid(1),'%g\t',Tice_tg);%%88
fprintf(fid(1),'%g\t',Vice_tg);%%89
%%%
fprintf(fid(1),'%g\t',CK1_tg);%%90
fprintf(fid(1),'%g\t',t);%%91
fprintf(fid(1),'%g\t\n',CKt);%%92
%%%%
if t==N_time_step
    fclose(fid(1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SPATIAL STANDARD DEVIATION
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
std_Pr_tg = std(Pr_S(Kinde)); %%% [mm]
std_Ta_tg = std(Ta_S(Kinde)); %%% [°C]
std_Ws_tg = std(Ws_S(Kinde)); %%% [m/s]
std_Ds_tg = std(Ds_S(Kinde)); %%% [°C]
std_ea_tg = std(ea_S(Kinde)); %%% [Pa]
std_N_tg = std(N_S(Kinde)); %%% [-]
std_Pre_tg = std(Pre_S(Kinde)); %%% [mbar]
std_Tdew_tg = std(Tdew_S(Kinde)); %%% [°C]
std_Rsw_tg = std(Rsw_space(Kinde)); %%% [W/m^2]
std_PAR_tg= std(PAR_space(Kinde));%% [W/m^2]
std_Ca_tg= std(Ca_S(Kinde));%% [W/m^2]
std_Ts_tg  =  std(Ts(Kinde));%% [°C]
std_Tdamp_tg  = std(Tdamp(Kinde));%% [°C]
std_Pr_sno_tg  =  std(Pr_sno(Kinde));%% [mm]
std_Pr_liq_tg = std(Pr_liq(Kinde)); %%[mm]
%%%%
%%%%
std_Rn_tg=std(Rn(Kinde));%% [W/m^2]
std_H_tg=std(H(Kinde));%% [W/m^2]
std_G_tg=std(G(Kinde));%% [W/m^2]
std_Gfin_tg=std(Gfin(Kinde));%% [W/m^2]
std_QE_tg=std(QE(Kinde));%% [W/m^2]
std_Qv_tg = std(Qv(Kinde)); %% [W/m^2]
std_Qfm_tg =  std(Qfm(Kinde));% [W/m^2]
%%%%
std_SWE_tg = std(SWE(Kinde)); %% [mm]
std_SND_tg = std(SND(Kinde)); %% [m]
std_WR_SP_tg = std(WR_SP(Kinde)); %%[]
std_U_SWE_tg = std(U_SWE(Kinde)); %%[]
std_NIn_SWE_tg  =  std(NIn_SWE(Kinde));%% [m/s]
std_dw_SNO_tg=std(dw_SNO(Kinde)); %%[]
std_ros_tg = std(ros(Kinde)); %% [kg/m^3]
std_In_SWE_tg =  std(In_SWE(Kinde)); %%[mm]
std_SP_wc_tg =  std(SP_wc(Kinde));%%[mm]
%%%%
std_ICE_tg = std(ICE(Kinde)); %% [mm]
std_ICE_D_tg = std(ICE_D(Kinde)); %% [m]
std_WR_IP_tg = std(WR_IP(Kinde)); %%[]
std_NIce_tg  =  std(NIce(Kinde));%%
std_IP_wc_tg =  std(IP_wc(Kinde));%%[mm]
%%%
std_T_H_tg = std(T_H_space(Kinde)); %%%[mm/h]
std_T_L_tg = std(T_L_space(Kinde)); %%%[mm/h]
std_EIn_L_tg  = std(EIn_L_space(Kinde)); %%% [mm/h]
std_EIn_H_tg =  std(EIn_H_space(Kinde)); %%% [mm/h]
std_EG_tg =std(EG(Kinde)); %%% [mm/h]
std_ESN_tg =std(ESN(Kinde)+ESN_In(Kinde)); %%% [mm/h]
std_EWAT_tg = std(EWAT(Kinde)); %%[mm/h]
std_EICE_tg = std(EICE(Kinde)); %%[mm/h]
%%%%%%%%%%%%%%%%
std_EIn_urb_tg = std(EIn_urb(Kinde)); %%[mm/h]
std_EIn_rock_tg = std(EIn_rock(Kinde)); %%[mm/h]
std_In_tg = std(In_H_space(Kinde) +  In_L_space(Kinde) +  SP_wc(Kinde) + In_SWE(Kinde) +  In_urb(Kinde) + In_rock(Kinde)+ IP_wc(Kinde) ); %%% [mm]
std_Inveg_tg = std(In_H_space(Kinde) +  In_L_space(Kinde));
%%%
std_WAT_tg = std(WAT(Kinde)); %% [mm]
std_FROCK_tg = std(FROCK(Kinde)); %% [mm]
%%%
std_f_tg = std(f(Kinde)*dth); %%[mm]
std_WIS_tg = std(WIS(Kinde)); %%[mm]
std_Rd_tg =std(Rd(Kinde)); %% [mm]
std_Rh_tg =std(Rh(Kinde)); %%[mm]
std_Lk_wat_tg =std(Lk_wat(Kinde)*dth); %% [mm]
std_Lk_rock_tg =std(Lk_rock(Kinde)*dth); %% [mm]
std_Lk_tg =std(Lk(Kinde)*dth); %% [mm]
std_OF_tg = std(OF(Kinde)); %%[]
std_OS_tg = std(OS(Kinde)); %%[]
std_ZWT_tg= std(ZWT(Kinde)); %% [mm]
%%%%%%%%
std_Qlat_in_tg=  std(Qi_in_space(Kinde)); %% [mm]
std_Qlat_out_tg =  std(Qi_out_space(Kinde)); %% [mm]
std_q_runon_tg = std(q_runon(Kinde)*dth); %% [mm]
std_Q_channel_tg = std(Q_channel(Kinde)); %% [mm]
%%%%%%%%%
std_O_tg =std(O_space(Kinde)); %%
std_V_tg =std(V_space(Kinde)); %%
std_Tdp_tg = std(Tdp_space(Kinde));
std_TsVEG_tg  = std(TsVEG(Kinde));%% [°C]
%%%%%%%%
std_er_tg =  std(er(Kinde)); %%  % [kg/s m^2]
std_r_soil_tg  =  std(r_soil(Kinde));%% [°C]
std_alp_soil_tg = std(alp_soil(Kinde)); %%[]
std_ra_tg  =  std(ra(Kinde));%% [m/s]
%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t==2
    tit4{1}=strcat('OUTPUT_',TITLE_SAVE,'_STD.dat');
    fid4(1)=fopen(tit4{1},'a');
end
%%%%%%%%%
fprintf(fid4(1),'%g\t',std_Pr_tg); %1
fprintf(fid4(1),'%g\t',std_Ta_tg); %2
fprintf(fid4(1),'%g\t',std_Ws_tg); %3
fprintf(fid4(1),'%g\t',std_Ds_tg); %4
fprintf(fid4(1),'%g\t',std_ea_tg); %5
fprintf(fid4(1),'%g\t',std_N_tg);  %6
fprintf(fid4(1),'%g\t',std_Pre_tg); %7
fprintf(fid4(1),'%g\t',std_Tdew_tg);%8
fprintf(fid4(1),'%g\t',std_Rsw_tg); %9
fprintf(fid4(1),'%g\t',std_PAR_tg);%10
fprintf(fid4(1),'%g\t',std_Ca_tg);%11
%%%%%
fprintf(fid4(1),'%g\t',std_Ts_tg);%12
fprintf(fid4(1),'%g\t',std_Tdamp_tg);%13
fprintf(fid4(1),'%g\t',std_Pr_sno_tg);%14
fprintf(fid4(1),'%g\t',std_Pr_liq_tg);%15
%%%%%%%%%
fprintf(fid4(1),'%g\t',std_Rn_tg);%16
fprintf(fid4(1),'%g\t',std_H_tg);%17
fprintf(fid4(1),'%g\t',std_G_tg);%18
fprintf(fid4(1),'%g\t',std_Gfin_tg);%19
fprintf(fid4(1),'%g\t',std_QE_tg);%20
fprintf(fid4(1),'%g\t',std_Qv_tg);%21
fprintf(fid4(1),'%g\t',std_Qfm_tg);%22
%%%%
fprintf(fid4(1),'%g\t',std_SWE_tg);%23
fprintf(fid4(1),'%g\t',std_SND_tg);%24
fprintf(fid4(1),'%g\t',std_WR_SP_tg);%25
fprintf(fid4(1),'%g\t',std_dw_SNO_tg);%26
fprintf(fid4(1),'%g\t',std_ros_tg);%27
fprintf(fid4(1),'%g\t',std_In_SWE_tg);%28
fprintf(fid4(1),'%g\t',std_SP_wc_tg);%29
%%%%%%%%%
fprintf(fid4(1),'%g\t',std_ICE_tg);%30
fprintf(fid4(1),'%g\t',std_ICE_D_tg);%31
fprintf(fid4(1),'%g\t',std_WR_IP_tg);%32
fprintf(fid4(1),'%g\t',std_NIce_tg);%33
fprintf(fid4(1),'%g\t',std_IP_wc_tg);%34
%%%%%%%%%
fprintf(fid4(1),'%g\t',std_T_H_tg);%35
fprintf(fid4(1),'%g\t',std_T_L_tg );%%36
fprintf(fid4(1),'%g\t',std_EIn_H_tg);%%37
fprintf(fid4(1),'%g\t',std_EIn_L_tg);%%38
fprintf(fid4(1),'%g\t',std_EG_tg);%%39
fprintf(fid4(1),'%g\t',std_ESN_tg);%%40
fprintf(fid4(1),'%g\t',std_EWAT_tg);%%41
fprintf(fid4(1),'%g\t',std_EICE_tg);%%42
fprintf(fid4(1),'%g\t',std_EIn_urb_tg);%%43
fprintf(fid4(1),'%g\t',std_EIn_rock_tg);%%44
fprintf(fid4(1),'%g\t',std_In_tg);%%45
fprintf(fid4(1),'%g\t',std_Inveg_tg);%%46
%%%%
fprintf(fid4(1),'%g\t',std_WAT_tg);%%47
fprintf(fid4(1),'%g\t',std_FROCK_tg);%%48
%%%%
fprintf(fid4(1),'%g\t',std_f_tg);%%49
fprintf(fid4(1),'%g\t',std_WIS_tg);%%50
fprintf(fid4(1),'%g\t',std_Rd_tg);%%51
fprintf(fid4(1),'%g\t',std_Rh_tg);%%52
fprintf(fid4(1),'%g\t',std_Lk_tg);%%53
fprintf(fid4(1),'%g\t',std_Lk_wat_tg);%%54
fprintf(fid4(1),'%g\t',std_Lk_rock_tg);%%55
fprintf(fid4(1),'%g\t',std_OF_tg);%%56
fprintf(fid4(1),'%g\t',std_OS_tg);%%57
fprintf(fid4(1),'%g\t',std_ZWT_tg);%%58
%%%
fprintf(fid4(1),'%g\t',std_Qlat_in_tg);%%59
fprintf(fid4(1),'%g\t',std_Qlat_out_tg);%%60
fprintf(fid4(1),'%g\t',std_q_runon_tg);%%61
fprintf(fid4(1),'%g\t',std_Q_channel_tg);%%62
fprintf(fid4(1),'%g\t',std_V_tg);%%63
fprintf(fid4(1),'%g\t',std_O_tg);%%64
fprintf(fid4(1),'%g\t',std_Tdp_tg);%%65
fprintf(fid4(1),'%g\t',std_TsVEG_tg);%%66
%%%
fprintf(fid4(1),'%g\t',std_er_tg);%%67
fprintf(fid4(1),'%g\t',std_r_soil_tg);%68
fprintf(fid4(1),'%g\t',std_alp_soil_tg);%69
fprintf(fid4(1),'%g\t\n',std_ra_tg );%70
%%%
if t==N_time_step
    fclose(fid4(1));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


OPT_VAR_SM_ON=0;


if OPT_VAR_SM_ON == 1 
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%% SPATIAL VARIABILITY OF SOIL MOISTURE
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% V, f in n-coordinate
    %%% O_tg O_tg_prec dth O_space O_space_prec  Qi_in_prec  std_O_tg
    f_ver_tg = mean(Asur(Kinde).*f(Kinde)*dth); %%[mm]
    %%%
    T_H_tg_n = mean(T_H_space(Kinde)./Asur(Kinde)); %%%[mm/h]
    T_L_tg_n = mean(T_L_space(Kinde)./Asur(Kinde)); %%%[mm/h]
    EG_tg_n =mean(EG(Kinde)./Asur(Kinde)); %%% [mm/h]
    Rd_tg_n =mean(Rd(Kinde)./Asur(Kinde)); %% [mm]
    Lk_tg_n =mean(Lk(Kinde)./Asur(Kinde)*dth); %% [mm]
    Qlat_out_tg_n =  mean(Qi_out_space(Kinde)./Asur(Kinde)); %% [mm]
    %%%%%%%%%%%
    V_space_prec = sum(V_prec,2);
    O_space_prec =  V_space_prec./Zs_OUT + Ohy_OUT ; %%[-]
    %V_tg_prec = mean(Asur(Kinde).*V_space_prec(Kinde)); %% [mm]
    O_tg_prec = mean(O_space_prec(Kinde)); %% [-]
    Qi_in_space_prec = sum(Qi_in_prec*dth,2);
    Qlat_in_tg_prec =  mean(Qi_in_space_prec(Kinde)); %% [mm]
    %%%
    Qlat_in_tg_prec_n=  mean(Qi_in_space_prec(Kinde)./Asur(Kinde)); %% [mm]
    %%%%%
    Qi_in_prec = Qi_in ; %%%
    V_prec = V ;
    %V_space =  sum(V,2); %%%
    %O_space =  V_space./Zs_OUT + Ohy_OUT ; %%[-]
    %V_tg = mean(Asur(Kinde).*V_space(Kinde)); %% [mm]
    %O_tg = mean(O_space(Kinde)); %% [-]
    %%%%%%%%%% Computation in n-coordinate
    Te1 = (O_space-O_space_prec)/dth - (O_tg - O_tg_prec)/dth; %%[1/h]
    Te2 = 2*(O_space-O_tg); %% [-]
    dOspacep2dt = mean(Te2(Kinde).*Te1(Kinde));
    %dOspacedt2 =  mean( (1./Zs_OUT(Kinde)).*(Te2(Kinde).*(Asur.*f(Kinde)-f_ver_tg/dth) - Te2(Kinde).*(EG(Kinde)-EG_tg) - Te2(Kinde).*(T_H_space(Kinde)+T_L_space(Kinde)-T_H_tg-T_L_tg) + ...
    %     -Te2(Kinde).*(Lk(Kinde)-Lk_tg/dth) + Te2(Kinde).*(Qi_in_space_prec(Kinde)/dth-Qlat_in_tg_prec/dth) - Te2(Kinde).*(Qi_out_space(Kinde)/dth-Qlat_out_tg/dth) +...
    %    -Te2(Kinde).*(Rd(Kinde)/dth-Rd_tg/dth)) );
    %%%%%%%%%
    %Cov_f  = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(f(Kinde)-f_tg/dth) ));
    %Cov_Eg =  mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(EG(Kinde)-EG_tg) ));
    %Cov_T =  mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(T_H_space(Kinde)+T_L_space(Kinde)-T_H_tg-T_L_tg) ));
    %Cov_Lk = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(Lk(Kinde)-Lk_tg/dth) ));
    %Cov_Qin = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(Qi_in_space_prec(Kinde)/dth-Qlat_in_tg_prec/dth) ));
    %Cov_Qout = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(Qi_out_space(Kinde)/dth-Qlat_out_tg/dth) ));
    %Cov_Rd = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(Rd(Kinde)/dth-Rd_tg/dth) ));
    %%%%%%%
    Cov_f  = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(f(Kinde)-f_tg/dth) ));
    Cov_Eg =  mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(EG(Kinde)-EG_tg_n) ));
    Cov_T =  mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(T_H_space(Kinde)./Asur+T_L_space(Kinde)./Asur-T_H_tg_n-T_L_tg_n) ));
    Cov_Lk = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(Lk(Kinde)./Asur-Lk_tg_n/dth) ));
    Cov_Qin = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(Qi_in_space_prec(Kinde)/dth./Asur-Qlat_in_tg_prec_n/dth) ));
    Cov_Qout = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(Qi_out_space(Kinde)/dth./Asur-Qlat_out_tg_n/dth) ));
    Cov_Rd = mean( (1./Zs_OUT(Kinde)).*( Te2(Kinde).*(Rd(Kinde)/dth./Asur-Rd_tg_n/dth) ));
    %%%%%%
    dOspacep2dt2 = Cov_f - Cov_Eg - Cov_T - Cov_Lk + Cov_Qin - Cov_Qout - Cov_Rd;
    A_var = Cov_f - Cov_Eg  - Cov_Lk + Cov_Qin - Cov_Qout - Cov_Rd;
    B_var = Cov_T;
    %%%%%%%%%%%%%%%% %% Budget in v - Coordinate
    dOspacedt =  (V_tg-V_tgtm1) ;
    dOspacedt2 = f_ver_tg  - EG_tg*dth - T_tg*dth - Lk_tg +  Qlat_in_tgtm1 - Qlat_out_tg - Rd_tg ;
    Ckk4 = (V_tgtm1 - V_tg) + f_ver_tg  - EG_tg*dth - T_tg*dth - Lk_tg +  Qlat_in_tgtm1 + (- Qlat_in_tg  +  Qlat_in_tg) - Qlat_out_tg - Rd_tg ;
    %Ckk44 = f_tg + (V_tgtm1 - V_tg) - EG_tg*dth - T_tg*dth - Lk_tg - Qlat_in_tg -Rd_tg + Qlat_in_tgtm1 - Qsub_exit ;
    A_mu =  mean(1./Zs_OUT(Kinde)).*(f_ver_tg  - EG_tg*dth - Lk_tg +  Qlat_in_tgtm1 - Qlat_out_tg - Rd_tg);
    B_mu =  mean(1./Zs_OUT(Kinde)).*(T_tg*dth);
    %%%%%%%%%%%%%
    OCv_tg = std_O_tg/O_tg;
    dCv_1 =  (A_var/(2*O_tg*std_O_tg));
    dCv_2 =  A_mu*OCv_tg/O_tg;
    dCv_3 =  B_mu*OCv_tg/O_tg;
    dCv_4 =  (B_var/(2*O_tg*std_O_tg));
    dCvdt = dCv_1 - dCv_2 + dCv_3 - dCv_4;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    deltaTSO = (O_space - O_tg)./O_tg ; %%% [-] Mean Relative difference
    deltaTSO2 = (deltaTSO).^2 ; %%% [-] Square of Mean Relative difference
    %%%%Ohy_OUT Osat_OUT
    SEspace = (O_space-Ohy_OUT)./(Osat_OUT-Ohy_OUT);
    SE_tg = mean( SEspace(Kinde) ); %% Effective saturation [-]
    std_SE_tg =std(SEspace(Kinde)); %%
    deltaTSS = (SEspace - SE_tg)./SE_tg ; %%% [-] Mean Relative difference
    deltaTSS2 = (deltaTSS).^2 ; %%% [-] Square of Mean Relative difference
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [WL_Mdep1]=Evaporation_layers(Zs,Mdep1); %%% Weighiting layers
    OMdep1 = (O*WL_Mdep1');%% Weighting layers 1
    OMdep1_tg = mean(OMdep1(Kinde));
    std_OMdep1_tg =std(OMdep1(Kinde)); %%
    [WL_Mdep2]=Evaporation_layers(Zs,Mdep2); %%% Weighiting layers
    OMdep2 = (O*WL_Mdep2');%% Weighting layers 2
    OMdep2_tg = mean(OMdep2(Kinde));
    std_OMdep2_tg =std(OMdep2(Kinde)); %%
    [WL_Mdep3]=Evaporation_layers(Zs,Mdep3); %%% Weighiting layers
    OMdep3 = (O*WL_Mdep3');%% Weighting layers 3
    OMdep3_tg = mean(OMdep3(Kinde));
    std_OMdep3_tg =std(OMdep3(Kinde)); %%
    %%%%%%%%%%%%%%%%
    
    if t==2
        tit8{1}=strcat('OUTPUT_',TITLE_SAVE,'_VAR_SM.dat');
        fid8(1)=fopen(tit8{1},'a');
    end
    fprintf(fid8(1),'%g\t',dOspacep2dt); %1
    fprintf(fid8(1),'%g\t',Cov_f); %2
    fprintf(fid8(1),'%g\t',Cov_Eg); %3
    fprintf(fid8(1),'%g\t',Cov_T); %4
    fprintf(fid8(1),'%g\t',Cov_Lk); %5
    fprintf(fid8(1),'%g\t',Cov_Qin); %6
    fprintf(fid8(1),'%g\t',Cov_Qout); %7
    fprintf(fid8(1),'%g\t',Cov_Rd); %8
    fprintf(fid8(1),'%g\t',dOspacep2dt2); %9
    fprintf(fid8(1),'%g\t',A_var); %10
    fprintf(fid8(1),'%g\t',B_var); %11
    %%%%%%
    fprintf(fid8(1),'%g\t',dOspacedt); %12
    fprintf(fid8(1),'%g\t',dOspacedt2); %13
    fprintf(fid8(1),'%g\t',Ckk4); %14
    fprintf(fid8(1),'%g\t',A_mu); %15
    fprintf(fid8(1),'%g\t',B_mu); %16
    %%%%%
    fprintf(fid8(1),'%g\t',dCv_1); %17
    fprintf(fid8(1),'%g\t',dCv_2); %18
    fprintf(fid8(1),'%g\t',dCv_3); %19
    fprintf(fid8(1),'%g\t',dCv_4); %20
    fprintf(fid8(1),'%g\t',dCvdt); %21
    fprintf(fid8(1),'%g\t',SE_tg); %22
    fprintf(fid8(1),'%g\t',std_SE_tg); %23
    fprintf(fid8(1),'%g\t',OMdep1_tg); %24
    fprintf(fid8(1),'%g\t',std_OMdep1_tg); %25
    fprintf(fid8(1),'%g\t',OMdep2_tg); %26
    fprintf(fid8(1),'%g\t',std_OMdep2_tg); %27
    fprintf(fid8(1),'%g\t',OMdep3_tg); %28
    fprintf(fid8(1),'%g\t',std_OMdep3_tg); %29
    fprintf(fid8(1),'%g\t\n',t);%%30
    
    %%%%
    if t==N_time_step
        fclose(fid8(1));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SPATIAL AVERAGE over the Vegetation type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for ijki=1:cc_max
    for ievc=1:length(EVcode)
        riev= find(ksv==EVcode(ievc));
        %%%%%%
        OH_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*OH(riev,ijki)); %%[]
        std_OH_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*OH(riev,ijki)); %%[]
        OL_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*OL(riev,ijki)); %%[]
        std_OL_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*OL(riev,ijki)); %%[]
        An_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*An_H(riev,ijki)); %%[]
        std_An_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*An_H(riev,ijki)); %%[]
        An_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*An_L(riev,ijki)); %%[]
        std_An_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*An_L(riev,ijki)); %%[]
        Tdp_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*Tdp_H(riev,ijki)); %%[]
        std_Tdp_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*Tdp_H(riev,ijki)); %%[]
        Tdp_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*Tdp_L(riev,ijki)); %%[]
        std_Tdp_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*Tdp_L(riev,ijki)); %%[]
        Rdark_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*Rdark_H(riev,ijki)); %%[]
        std_Rdark_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*Rdark_H(riev,ijki)); %%[]
        Rdark_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*Rdark_L(riev,ijki)); %%[]
        std_Rdark_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*Rdark_L(riev,ijki)); %%[]
        LAI_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*LAI_H(riev,ijki)); %%[]
        std_LAI_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*LAI_H(riev,ijki)); %%[]
        LAI_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*LAI_L(riev,ijki)); %%[]
        std_LAI_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*LAI_L(riev,ijki)); %%[]
        NPP_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*NPP_H(riev,ijki)); %%[]
        std_NPP_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*NPP_H(riev,ijki)); %%[]
        NPP_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*NPP_L(riev,ijki)); %%[]
        std_NPP_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*NPP_L(riev,ijki)); %%[]
        ANPP_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*ANPP_H(riev,ijki)); %%[]
        std_ANPP_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*ANPP_H(riev,ijki)); %%[]
        ANPP_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*ANPP_L(riev,ijki)); %%[]
        std_ANPP_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*ANPP_L(riev,ijki)); %%[]
        RA_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*RA_H(riev,ijki)); %%[]
        std_RA_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*RA_H(riev,ijki)); %%[]
        RA_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*RA_L(riev,ijki)); %%[]
        std_RA_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*RA_L(riev,ijki)); %%[]
        Rg_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*Rg_H(riev,ijki)); %%[]
        std_Rg_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*Rg_H(riev,ijki)); %%[]
        Rg_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*Rg_L(riev,ijki)); %%[]
        std_Rg_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*Rg_L(riev,ijki)); %%[]
        LAIdead_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*LAIdead_H(riev,ijki)); %%[]
        std_LAIdead_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*LAIdead_H(riev,ijki)); %%[]
        LAIdead_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*LAIdead_L(riev,ijki)); %%[]
        std_LAIdead_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*LAIdead_L(riev,ijki)); %%[]
        hc_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*hc_H(riev,ijki)); %%[]
        std_hc_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*hc_H(riev,ijki)); %%[]
        hc_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*hc_L(riev,ijki)); %%[]
        std_hc_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*hc_L(riev,ijki)); %%[]
        AgeL_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*AgeL_H(riev,ijki)); %%[]
        std_AgeL_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*AgeL_H(riev,ijki)); %%[]
        AgeL_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*AgeL_L(riev,ijki)); %%[]
        std_AgeL_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*AgeL_L(riev,ijki)); %%[]
        SAI_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*SAI_H(riev,ijki)); %%[]
        std_SAI_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*SAI_H(riev,ijki)); %%[]
        SAI_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*SAI_L(riev,ijki)); %%[]
        std_SAI_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*SAI_L(riev,ijki)); %%[]
        PHE_S_H_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*PHE_S_H(riev,ijki)); %%[]
        std_PHE_S_H_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*PHE_S_H(riev,ijki)); %%[]
        PHE_S_L_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*PHE_S_L(riev,ijki)); %%[]
        std_PHE_S_L_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*PHE_S_L(riev,ijki)); %%[]
        %Llitter_tg(ievc,ijki) =  mean(Ccrown_OUT(ievc,ijki)*Llitter(riev,ijki)); %%[]
        %std_Llitter_tg(ievc,ijki) = std(Ccrown_OUT(ievc,ijki)*Llitter(riev,ijki)); %%[]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %,B_H(ij,:,:),,,
        %,B_L(ij,:,:),
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t==2
    for ijki=1:cc_max
        for ievc=1:length(EVcode)
            if Ccrown_OUT(ievc,ijki)>0
                tit2{ievc,ijki}=strcat('OUTPUT_',TITLE_SAVE,'_AVG_PFT_',num2str(ijki),'_code_',num2str(EVcode(ievc)),'.dat');
                fid2(ievc,ijki)=fopen(tit2{ievc,ijki},'a');
            end
        end
    end
end
for ijki=1:cc_max
    for ievc=1:length(EVcode)
        if Ccrown_OUT(ievc,ijki)>0
            fprintf(fid2(ievc,ijki),'%g\t',OH_tg(ievc,ijki)); %1
            fprintf(fid2(ievc,ijki),'%g\t',OL_tg(ievc,ijki)); %2
            fprintf(fid2(ievc,ijki),'%g\t',An_H_tg(ievc,ijki));%3
            fprintf(fid2(ievc,ijki),'%g\t',An_L_tg(ievc,ijki));%4
            fprintf(fid2(ievc,ijki),'%g\t',Rdark_H_tg(ievc,ijki));%5
            fprintf(fid2(ievc,ijki),'%g\t',Rdark_L_tg(ievc,ijki));%6
            fprintf(fid2(ievc,ijki),'%g\t',Tdp_H_tg(ievc,ijki));%7
            fprintf(fid2(ievc,ijki),'%g\t',Tdp_L_tg(ievc,ijki));%8
            %%%
            fprintf(fid2(ievc,ijki),'%g\t',LAI_H_tg(ievc,ijki));%9
            fprintf(fid2(ievc,ijki),'%g\t',LAI_L_tg(ievc,ijki));%10
            fprintf(fid2(ievc,ijki),'%g\t',NPP_H_tg(ievc,ijki));%11
            fprintf(fid2(ievc,ijki),'%g\t',NPP_L_tg(ievc,ijki));%12
            fprintf(fid2(ievc,ijki),'%g\t',ANPP_H_tg(ievc,ijki));%13
            fprintf(fid2(ievc,ijki),'%g\t',ANPP_L_tg(ievc,ijki));%14
            fprintf(fid2(ievc,ijki),'%g\t',RA_H_tg(ievc,ijki));%15
            fprintf(fid2(ievc,ijki),'%g\t',RA_L_tg(ievc,ijki));%16
            fprintf(fid2(ievc,ijki),'%g\t',Rg_H_tg(ievc,ijki));%17
            fprintf(fid2(ievc,ijki),'%g\t',Rg_L_tg(ievc,ijki));%18
            fprintf(fid2(ievc,ijki),'%g\t',PHE_S_H_tg(ievc,ijki));%19
            fprintf(fid2(ievc,ijki),'%g\t',PHE_S_L_tg(ievc,ijki));%20
            fprintf(fid2(ievc,ijki),'%g\t',AgeL_H_tg(ievc,ijki));%21
            fprintf(fid2(ievc,ijki),'%g\t',AgeL_L_tg(ievc,ijki));%22
            fprintf(fid2(ievc,ijki),'%g\t',SAI_H_tg(ievc,ijki));%23
            fprintf(fid2(ievc,ijki),'%g\t',SAI_L_tg(ievc,ijki));%24
            fprintf(fid2(ievc,ijki),'%g\t',LAIdead_H_tg(ievc,ijki));%25
            fprintf(fid2(ievc,ijki),'%g\t',LAIdead_L_tg(ievc,ijki));%26
            fprintf(fid2(ievc,ijki),'%g\t',hc_H_tg(ievc,ijki));%27
            fprintf(fid2(ievc,ijki),'%g\t',hc_L_tg(ievc,ijki));%28
            %fprintf(fid2(ievc,ijki),'%g\t\n',Llitter_tg(ievc,ijki));%29
        end
    end
end

%%%
if t==N_time_step
    for ijki=1:cc_max
        for ievc=1:length(EVcode)
            if Ccrown_OUT(ievc,ijki)>0
                fclose(fid2(ievc,ijki));
            end
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% SPATIAL STD over the vegetation type
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t==2
    for ijki=1:cc_max
        for ievc=1:length(EVcode)
            if Ccrown_OUT(ievc,ijki)>0
                tit3{ievc,ijki}=strcat('OUTPUT_',TITLE_SAVE,'_STD_PFT_',num2str(ijki),'_code_',num2str(EVcode(ievc)),'.dat');
                fid3(ievc,ijki)=fopen(tit3{ievc,ijki},'a');
            end
        end
    end
end
%%%%%%%%%%
for ijki=1:cc_max
    for ievc=1:length(EVcode)
        if Ccrown_OUT(ievc,ijki)>0
            fprintf(fid3(ievc,ijki),'%g\t',std_OH_tg(ievc,ijki)); %1
            fprintf(fid3(ievc,ijki),'%g\t',std_OL_tg(ievc,ijki)); %2
            fprintf(fid3(ievc,ijki),'%g\t',std_An_H_tg(ievc,ijki));%3
            fprintf(fid3(ievc,ijki),'%g\t',std_An_L_tg(ievc,ijki));%4
            fprintf(fid3(ievc,ijki),'%g\t',std_Rdark_H_tg(ievc,ijki)); %5
            fprintf(fid3(ievc,ijki),'%g\t',std_Rdark_L_tg(ievc,ijki));%6
            fprintf(fid3(ievc,ijki),'%g\t',std_Tdp_H_tg(ievc,ijki)); %7
            fprintf(fid3(ievc,ijki),'%g\t',std_Tdp_L_tg(ievc,ijki));%8
            %%%
            fprintf(fid3(ievc,ijki),'%g\t',std_LAI_H_tg(ievc,ijki));%9
            fprintf(fid3(ievc,ijki),'%g\t',std_LAI_L_tg(ievc,ijki));%10
            fprintf(fid3(ievc,ijki),'%g\t',std_NPP_H_tg(ievc,ijki));%11
            fprintf(fid3(ievc,ijki),'%g\t',std_NPP_L_tg(ievc,ijki));%12
            fprintf(fid3(ievc,ijki),'%g\t',std_ANPP_H_tg(ievc,ijki));%13
            fprintf(fid3(ievc,ijki),'%g\t',std_ANPP_L_tg(ievc,ijki));%14
            fprintf(fid3(ievc,ijki),'%g\t',std_RA_H_tg(ievc,ijki));%15
            fprintf(fid3(ievc,ijki),'%g\t',std_RA_L_tg(ievc,ijki));%16
            fprintf(fid3(ievc,ijki),'%g\t',std_Rg_H_tg(ievc,ijki));%17
            fprintf(fid3(ievc,ijki),'%g\t',std_Rg_L_tg(ievc,ijki));%18
            fprintf(fid3(ievc,ijki),'%g\t',std_PHE_S_H_tg(ievc,ijki)); %19
            fprintf(fid3(ievc,ijki),'%g\t',std_PHE_S_L_tg(ievc,ijki)); %20
            fprintf(fid3(ievc,ijki),'%g\t',std_AgeL_H_tg(ievc,ijki));%21
            fprintf(fid3(ievc,ijki),'%g\t',std_AgeL_L_tg(ievc,ijki)); %22
            fprintf(fid3(ievc,ijki),'%g\t',std_SAI_H_tg(ievc,ijki)); %23
            fprintf(fid3(ievc,ijki),'%g\t',std_SAI_L_tg(ievc,ijki)); %24
            fprintf(fid3(ievc,ijki),'%g\t',std_LAIdead_H_tg(ievc,ijki)); %25
            fprintf(fid3(ievc,ijki),'%g\t',std_LAIdead_L_tg(ievc,ijki)); %26
            fprintf(fid3(ievc,ijki),'%g\t',std_hc_H_tg(ievc,ijki)); %27
            fprintf(fid3(ievc,ijki),'%g\t',std_hc_L_tg(ievc,ijki)); %28
            %fprintf(fid3(ievc,ijki),'%g\t\n',std_Llitter_tg(ievc,ijki)); %29
        end
    end
end
%%%
if t==N_time_step
    for ijki=1:cc_max
        for ievc=1:length(EVcode)
            if Ccrown_OUT(ievc,ijki)>0
                fclose(fid3(ievc,ijki));
            end
        end
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% TEMPORAL AVERAGE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if t==2
    toutp=2;
else
    toutp=toutp+1;
end
%%%
if toutp==2
    %%%%%
    Pr_spatial=Pr_S;
    Ta_spatial=Ta_S;
    Ws_spatial=Ws_S;
    Ds_spatial=Ds_S;
    ea_spatial=ea_S;
    N_spatial=N_S;
    Pre_spatial=Pre_S;
    Tdew_spatial=Tdew_S;
    Rsw_spatial=Rsw_space;
    PAR_spatial=PAR_space;
    Ca_spatial = Ca_S;
    %%%
    Ts_spatial=Ts;
    Tdamp_spatial=Tdamp;
    Csno_spatial=Csno;
    Cice_spatial=Cice;
    Csnow_spatial=Csnow;
    Cicew_spatial=Cicew;
    Pr_sno_spatial=Pr_sno;
    Pr_liq_spatial=Pr_liq;
    %%%%
    Rn_spatial=Rn;
    H_spatial=H;
    G_spatial=G;
    Gfin_spatial=Gfin;
    QE_spatial=QE;
    Qv_spatial=Qv;
    Qfm_spatial=Qfm;
    %%%%
    SWE_spatial = SWE;
    SND_spatial = SND ;
    WR_SP_spatial = WR_SP;
    U_SWE_spatial = U_SWE;
    dw_SNO_spatial = dw_SNO;
    ros_spatial = ros;
    In_SWE_spatial=  In_SWE;
    SP_wc_spatial =  SP_wc ;
    SWE_avalanched_spatial = SWE_avalanched; 
    %%%%
    ICE_spatial = ICE;
    ICE_D_spatial = ICE_D ;
    WR_IP_spatial = WR_IP;
    NICe_spatial = NIce;
    IP_wc_spatial =  IP_wc ;
    %%%%
    Imelt_spatial = Imelt; 
    Smelt_spatial = Smelt; 
    Tice_spatial = Tice;     
    T_H_spatial=T_H_space;
    T_L_spatial=T_L_space;
    EIn_H_spatial=EIn_H_space;
    EIn_L_spatial=EIn_L_space;
    EG_spatial = EG;
    ESN_spatial = ESN + ESN_In;
    EWAT_spatial = EWAT ;
    EICE_spatial = EICE ;
    Dr_H_spatial=Dr_H_space;
    Dr_L_spatial=Dr_L_space;
    EIn_rock_spatial = EIn_rock;
    EIn_urb_spatial = EIn_urb;
    SE_rock_spatial =  SE_rock;
    SE_urb_spatial =  SE_urb;
    %%%%
    In_spatial =  In_H_space + In_L_space + SP_wc + In_SWE  + In_urb + In_rock + IP_wc;
    Inveg_spatial = In_H_space + In_L_space;
    In_rock_spatial = In_rock;
    WAT_spatial = WAT;
    FROCK_spatial = FROCK ;
    %%%%%
    f_spatial = f;
    WIS_spatial = WIS;
    Rd_spatial =Rd;
    Rh_spatial =Rh;
    Lk_spatial =Lk;
    Lk_rock_spatial =Lk_rock;
    Lk_wat_spatial =Lk_wat;
    OF_spatial =OF;
    OS_spatial = OS;
    ZWT_spatial=ZWT;
    %%%%%%%%
    Qlat_in_spatial=  Qi_in_space;
    Qlat_out_spatial =  Qi_out_space;
    q_runon_spatial = q_runon;
    Q_channel_spatial = Q_channel;
    %%%%%%%%%
    O_spatial=O_space;
    V_spatial=V_space;
    Vice_spatial=Vice_space;
    Tdp_spatial = Tdp_space;
    Tdpsnow_spatial=Tdpsnow_space;
    SAT_spatial = (ZWT==0);
    TsVEG_spatial = TsVEG;
    Ts_under_spatial = Ts_under;
    %%%%%%%%%%%%%%%%
    er_spatial =  er ;
    DQ_S_spatial =DQ_S;
    DT_S_spatial = DT_S;
    dQ_S_spatial =dQ_S;
    CK1_spatial =  CK1;
    %CK2_spatial = CK2;
    %%%%
    OH_spatial = OH;
    OL_spatial = OL;
    An_H_spatial = An_H;
    An_L_spatial=  An_L ;
    Rdark_H_spatial =Rdark_H;
    Rdark_L_spatial = Rdark_L;
    Tdp_H_spatial =Tdp_H;
    Tdp_L_spatial = Tdp_L;
    Ci_sunH_spatial =Ci_sunH;
    Ci_sunL_spatial = Ci_sunL;
    Ci_shdH_spatial =Ci_shdH;
    Ci_shdL_spatial = Ci_shdL;
    LAI_H_spatial = LAI_H;
    LAI_L_spatial = LAI_L;
    GPP_H_spatial = NPP_H+RA_H;
    GPP_L_spatial = NPP_L +RA_L;
    NPP_H_spatial = NPP_H;
    NPP_L_spatial = NPP_L;
    ANPP_H_spatial = ANPP_H;
    ANPP_L_spatial = ANPP_L;
    B_H_spatial = B_H;
    B_L_spatial = B_L;
    RA_H_spatial =RA_H;
    RA_L_spatial = RA_L;
    Rg_H_spatial =Rg_H;
    Rg_L_spatial = Rg_L;
    Rms_H_spatial =Rms_H;
    Rms_L_spatial = Rms_L;
    Rmr_H_spatial =Rmr_H;
    Rmr_L_spatial = Rmr_L;
    Rmc_H_spatial =Rmc_H;
    Rmc_L_spatial = Rmc_L;
    SAI_H_spatial = SAI_H;
    SAI_L_spatial = SAI_L;
    hc_H_spatial = hc_H;
    hc_L_spatial = hc_L;
    LAIdead_H_spatial = LAIdead_H;
    LAIdead_L_spatial = LAIdead_L;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
else
    %%%%%
    Pr_spatial = Pr_spatial + Pr_S;
    Ta_spatial = ((toutp-2)*Ta_spatial + Ta_S)/(toutp-1) ;
    Ws_spatial = ((toutp-2)*Ws_spatial + Ws_S)/(toutp-1) ;
    Ds_spatial = ((toutp-2)*Ds_spatial + Ds_S)/(toutp-1) ;
    ea_spatial = ((toutp-2)*ea_spatial + ea_S)/(toutp-1) ;
    N_spatial = ((toutp-2)*N_spatial + N_S)/(toutp-1) ;
    Pre_spatial = ((toutp-2)*Pre_spatial + Pre_S)/(toutp-1) ;
    Tdew_spatial = ((toutp-2)*Tdew_spatial + Tdew_S)/(toutp-1) ;
    Rsw_spatial = ((toutp-2)*Rsw_spatial + Rsw_space)/(toutp-1) ;
    PAR_spatial = ((toutp-2)*PAR_spatial + PAR_space)/(toutp-1) ;
    Ca_spatial = ((toutp-2)*Ca_spatial + Ca_S)/(toutp-1) ;
    %%%
    Ts_spatial = ((toutp-2)*Ts_spatial + Ts)/(toutp-1) ;
    Tdamp_spatial = ((toutp-2)*Tdamp_spatial + Tdamp)/(toutp-1) ;
    Csno_spatial = Csno_spatial + Csno ;
    Csnow_spatial = Csnow_spatial + Csnow ;
    Cice_spatial = Cice_spatial + Cice ;
    Cicew_spatial = Cicew_spatial + Cicew ;
    Pr_sno_spatial = Pr_sno_spatial + Pr_sno;
    Pr_liq_spatial = Pr_liq_spatial + Pr_liq;
    %%%
    Rn_spatial = ((toutp-2)*Rn_spatial + Rn)/(toutp-1) ;
    H_spatial = ((toutp-2)*H_spatial + H)/(toutp-1) ;
    G_spatial = ((toutp-2)*G_spatial + G)/(toutp-1) ;
    Gfin_spatial = ((toutp-2)*Gfin_spatial + Gfin)/(toutp-1) ;
    QE_spatial = ((toutp-2)*QE_spatial + QE)/(toutp-1) ;
    Qv_spatial = ((toutp-2)*Qv_spatial + Qv)/(toutp-1) ;
    Qfm_spatial = ((toutp-2)*Qfm_spatial + Qfm)/(toutp-1) ;
    %%%%
    SWE_spatial = ((toutp-2)*SWE_spatial + SWE)/(toutp-1) ;
    SND_spatial = ((toutp-2)*SND_spatial + SND)/(toutp-1) ;
    WR_SP_spatial = ((toutp-2)*WR_SP_spatial + WR_SP)/(toutp-1) ;
    U_SWE_spatial = ((toutp-2)*U_SWE_spatial + U_SWE)/(toutp-1) ;
    dw_SNO_spatial = ((toutp-2)*dw_SNO_spatial + dw_SNO)/(toutp-1) ;
    ros_spatial = ((toutp-2)*ros_spatial + ros)/(toutp-1) ;
    In_SWE_spatial= ((toutp-2)*In_SWE_spatial + In_SWE)/(toutp-1) ;
    SP_wc_spatial = ((toutp-2)*SP_wc_spatial + SP_wc)/(toutp-1) ;
    SWE_avalanched_spatial = ((toutp-2)*SWE_avalanched_spatial + SWE_avalanched)/(toutp-1) ;   
    %%%%
    ICE_spatial = ((toutp-2)*ICE_spatial + ICE)/(toutp-1) ;
    ICE_D_spatial= ((toutp-2)*ICE_D_spatial + ICE_D)/(toutp-1) ;
    WR_IP_spatial = ((toutp-2)*WR_IP_spatial + WR_IP)/(toutp-1) ;
    NICe_spatial = ((toutp-2)*NICe_spatial + NIce)/(toutp-1) ;
    IP_wc_spatial = ((toutp-2)*IP_wc_spatial + IP_wc)/(toutp-1) ;
    Imelt_spatial = ((toutp-2)*Imelt_spatial + Imelt)/(toutp-1) ;
    Smelt_spatial = ((toutp-2)*Smelt_spatial + Smelt)/(toutp-1) ;
    Tice_spatial = ((toutp-2)*Tice_spatial + Tice)/(toutp-1) ;  
    %%%%
    T_H_spatial = ((toutp-2)*T_H_spatial + T_H_space)/(toutp-1) ;
    T_L_spatial = ((toutp-2)*T_L_spatial + T_L_space)/(toutp-1) ;
    EIn_H_spatial = ((toutp-2)*EIn_H_spatial + EIn_H_space)/(toutp-1) ;
    EIn_L_spatial = ((toutp-2)*EIn_L_spatial + EIn_L_space)/(toutp-1) ;
    EG_spatial = ((toutp-2)*EG_spatial + EG)/(toutp-1) ;
    ESN_spatial = ((toutp-2)*ESN_spatial + ESN + ESN_In)/(toutp-1) ;
    EWAT_spatial = ((toutp-2)*EWAT_spatial + EWAT)/(toutp-1) ;
    EICE_spatial = ((toutp-2)*EICE_spatial + EICE)/(toutp-1) ;
    Dr_H_spatial = ((toutp-2)*Dr_H_spatial + Dr_H_space)/(toutp-1) ;
    Dr_L_spatial = ((toutp-2)*Dr_L_spatial + Dr_L_space)/(toutp-1) ;
    EIn_rock_spatial = ((toutp-2)*EIn_rock_spatial + EIn_rock)/(toutp-1) ;
    EIn_urb_spatial = ((toutp-2)*EIn_urb_spatial + EIn_urb)/(toutp-1) ;
    SE_rock_spatial = ((toutp-2)*SE_rock_spatial + SE_rock)/(toutp-1);
    SE_urb_spatial = ((toutp-2)*SE_urb_spatial + SE_urb)/(toutp-1);
    %%%%%
    In_spatial = ((toutp-2)*In_spatial + In_H_space + In_L_space + SP_wc + In_SWE + In_urb + In_rock + IP_wc)/(toutp-1) ;
    Inveg_spatial = ((toutp-2)*Inveg_spatial + In_H_space + In_L_space)/(toutp-1) ;
    In_rock_spatial =  ((toutp-2)*In_rock_spatial + In_rock)/(toutp-1) ;
    FROCK_spatial = ((toutp-2)*FROCK_spatial + FROCK)/(toutp-1) ;
    WAT_spatial = ((toutp-2)*WAT_spatial + WAT)/(toutp-1) ;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    f_spatial = ((toutp-2)*f_spatial +  f)/(toutp-1) ;
    WIS_spatial = ((toutp-2)*WIS_spatial +  WIS)/(toutp-1) ;
    Rd_spatial =((toutp-2)*Rd_spatial +  Rd)/(toutp-1) ;
    Rh_spatial =((toutp-2)*Rh_spatial +  Rh)/(toutp-1) ;
    Lk_spatial =((toutp-2)*Lk_spatial +  Lk)/(toutp-1) ;
    Lk_rock_spatial =((toutp-2)*Lk_rock_spatial +  Lk_rock)/(toutp-1) ;
    Lk_wat_spatial =((toutp-2)*Lk_wat_spatial +  Lk_wat)/(toutp-1) ;
    OF_spatial =((toutp-2)*OF_spatial +  OF)/(toutp-1) ;
    OS_spatial = ((toutp-2)*OS_spatial +  OS)/(toutp-1) ;
    ZWT_spatial=((toutp-2)*ZWT_spatial +  ZWT)/(toutp-1) ;
    %%%%%%%%
    Qlat_in_spatial=  ((toutp-2)*Qlat_in_spatial +  Qi_in_space)/(toutp-1) ;
    Qlat_out_spatial =  ((toutp-2)*Qlat_out_spatial +  Qi_out_space)/(toutp-1) ;
    q_runon_spatial = ((toutp-2)*q_runon_spatial +  q_runon)/(toutp-1) ;
    Q_channel_spatial = ((toutp-2)*Q_channel_spatial +  Q_channel)/(toutp-1) ; 
    %%%%%%%%%
    V_spatial= ((toutp-2)*V_spatial + V_space)/(toutp-1) ;
    Vice_spatial=((toutp-2)*Vice_spatial + Vice_space)/(toutp-1) ;
    O_spatial = ((toutp-2)*O_spatial + O_space)/(toutp-1) ;
    SAT_spatial = SAT_spatial + (ZWT==0) ;
    Tdp_spatial = ((toutp-2)*Tdp_spatial + Tdp_space)/(toutp-1) ;
    TsVEG_spatial = ((toutp-2)*TsVEG_spatial + TsVEG)/(toutp-1) ;
    Tdpsnow_spatial = ((toutp-2)*Tdpsnow_spatial + Tdpsnow_space)/(toutp-1) ;
    Ts_under_spatial = ((toutp-2)*Ts_under_spatial + Ts_under)/(toutp-1) ;
    %%%%%%%%%
    er_spatial = ((toutp-2)*er_spatial + er)/(toutp-1) ;
    DQ_S_spatial = DQ_S_spatial + DQ_S;
    DT_S_spatial = DT_S_spatial + DT_S;
    dQ_S_spatial = dQ_S_spatial + dQ_S;
    CK1_spatial = CK1_spatial + CK1;
    %CK2_spatial = CK2_spatial + CK2;
    %%%%%
    OH_spatial = ((toutp-2)*OH_spatial + OH)/(toutp-1) ;
    OL_spatial = ((toutp-2)*OL_spatial + OL)/(toutp-1) ;
    An_H_spatial = ((toutp-2)*An_H_spatial + An_H)/(toutp-1) ;
    An_L_spatial = ((toutp-2)*An_L_spatial + An_L)/(toutp-1) ;
    Rdark_H_spatial = ((toutp-2)*Rdark_H_spatial + Rdark_H)/(toutp-1) ;
    Rdark_L_spatial = ((toutp-2)*Rdark_L_spatial + Rdark_L)/(toutp-1) ;
    Ci_sunH_spatial = ((toutp-2)*Ci_sunH_spatial + Ci_sunH)/(toutp-1) ;
    Ci_sunL_spatial = ((toutp-2)*Ci_sunL_spatial + Ci_sunL)/(toutp-1) ;
    Ci_shdH_spatial = ((toutp-2)*Ci_shdH_spatial + Ci_shdH)/(toutp-1) ;
    Ci_shdL_spatial = ((toutp-2)*Ci_shdL_spatial + Ci_shdL)/(toutp-1) ;
    LAI_H_spatial = ((toutp-2)*LAI_H_spatial + LAI_H)/(toutp-1) ;
    LAI_L_spatial = ((toutp-2)*LAI_L_spatial + LAI_L)/(toutp-1) ;
    GPP_H_spatial = ((toutp-2)*GPP_H_spatial + NPP_H+RA_H)/(toutp-1) ;
    GPP_L_spatial = ((toutp-2)*GPP_L_spatial + NPP_L+RA_L)/(toutp-1) ;
    NPP_H_spatial = ((toutp-2)*NPP_H_spatial + NPP_H)/(toutp-1) ;
    NPP_L_spatial = ((toutp-2)*NPP_L_spatial + NPP_L)/(toutp-1) ;
    ANPP_H_spatial = ((toutp-2)*ANPP_H_spatial + ANPP_H)/(toutp-1) ;
    ANPP_L_spatial = ((toutp-2)*ANPP_L_spatial + ANPP_L)/(toutp-1) ;
    B_H_spatial = ((toutp-2)*B_H_spatial + B_H)/(toutp-1) ;
    B_L_spatial = ((toutp-2)*B_L_spatial + B_L)/(toutp-1) ;
    RA_H_spatial = ((toutp-2)*RA_H_spatial + RA_H)/(toutp-1) ;
    RA_L_spatial = ((toutp-2)*RA_L_spatial + RA_L)/(toutp-1) ;
    Rg_H_spatial = ((toutp-2)*Rg_H_spatial + Rg_H)/(toutp-1) ;
    Rg_L_spatial = ((toutp-2)*Rg_L_spatial + Rg_L)/(toutp-1) ;
    Rms_H_spatial = ((toutp-2)*Rms_H_spatial + Rms_H)/(toutp-1) ;
    Rms_L_spatial = ((toutp-2)*Rms_L_spatial + Rms_L)/(toutp-1) ;
    Rmr_H_spatial = ((toutp-2)*Rmr_H_spatial + Rmr_H)/(toutp-1) ;
    Rmr_L_spatial = ((toutp-2)*Rmr_L_spatial + Rmr_L)/(toutp-1) ;
    Rmc_H_spatial = ((toutp-2)*Rmc_H_spatial + Rmc_H)/(toutp-1) ;
    Rmc_L_spatial = ((toutp-2)*Rmc_L_spatial + Rmc_L)/(toutp-1) ;
    SAI_H_spatial = ((toutp-2)*SAI_H_spatial +  SAI_H)/(toutp-1) ;
    SAI_L_spatial = ((toutp-2)*SAI_L_spatial +  SAI_L)/(toutp-1) ;
    hc_H_spatial = ((toutp-2)*hc_H_spatial +  hc_H)/(toutp-1) ;
    hc_L_spatial = ((toutp-2)*hc_L_spatial +  hc_L)/(toutp-1) ;
    LAIdead_H_spatial = ((toutp-2)*LAIdead_H_spatial +  LAIdead_H)/(toutp-1) ;
    LAIdead_L_spatial = ((toutp-2)*LAIdead_L_spatial +  LAIdead_L)/(toutp-1) ;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if  length(intersect(t,tstore))==1 ||  t==N_time_step
    Title_save = strcat('OUTPUT_',TITLE_SAVE,'_SPATIAL_',num2str(t));
    save(Title_save,'DTM','cellsize','xllcorner','yllcorner','x_cell','y_cell','m_cell','n_cell','Pr_spatial',...
        'Ta_spatial','Ws_spatial','Ds_spatial','ea_spatial','N_spatial','Pre_spatial','Tdew_spatial',...
        'Rsw_spatial','PAR_spatial','Ca_spatial','Ts_spatial','Tdamp_spatial','Csno_spatial','Pr_sno_spatial',...
        'Csnow_spatial','Cice_spatial','Cicew_spatial',...
        'Pr_liq_spatial','Rn_spatial','H_spatial','G_spatial','Gfin_spatial','QE_spatial','Qv_spatial','Qfm_spatial',...
        'SWE_spatial','SND_spatial','WR_SP_spatial','ros_spatial','dw_SNO_spatial','U_SWE_spatial','In_SWE_spatial','SP_wc_spatial',...
        'ICE_spatial','ICE_D_spatial','WR_IP_spatial','NICe_spatial','IP_wc_spatial','SWE_avalanched_spatial',...
        'T_H_spatial','T_L_spatial','EIn_H_spatial','EIn_L_spatial','EG_spatial','ESN_spatial','EWAT_spatial','EICE_spatial','EIn_rock_spatial',...
        'EIn_urb_spatial','SE_rock_spatial','SE_urb_spatial',...
        'Dr_H_spatial','Dr_L_spatial','In_spatial','Inveg_spatial','In_rock_spatial','FROCK_spatial','WAT_spatial',...
        'f_spatial','WIS_spatial','Rd_spatial','Rh_spatial',...
        'Lk_spatial','Lk_rock_spatial','Lk_wat_spatial','OF_spatial','OS_spatial','ZWT_spatial','Qlat_in_spatial',...
        'Qlat_out_spatial','q_runon_spatial','Q_channel_spatial','O_spatial','V_spatial','SAT_spatial','Tdp_spatial','TsVEG_spatial',...
        'Ts_under_spatial','Tdpsnow_spatial',...
        'er_spatial','DQ_S_spatial','DT_S_spatial','dQ_S_spatial',...
        'Tice_spatial','Imelt_spatial','Smelt_spatial','Vice_spatial',...
        'CK1_spatial',...
        'OH_spatial','OL_spatial','An_H_spatial','An_L_spatial','Rdark_H_spatial',...
        'Rdark_L_spatial','Ci_sunH_spatial','Ci_sunL_spatial','Ci_shdH_spatial','Ci_shdL_spatial','LAI_H_spatial','LAI_L_spatial','GPP_H_spatial',...
        'GPP_L_spatial','NPP_H_spatial','NPP_L_spatial','ANPP_H_spatial','ANPP_L_spatial','B_H_spatial',...
        'B_L_spatial','RA_H_spatial','RA_L_spatial','Rg_H_spatial','Rg_L_spatial','Rms_H_spatial',...
        'Rms_L_spatial','Rmr_H_spatial','Rmr_L_spatial','Rmc_H_spatial','Rmc_L_spatial','SAI_H_spatial',...
        'SAI_L_spatial','hc_H_spatial','hc_L_spatial','LAIdead_H_spatial','LAIdead_L_spatial');
    %%%%%%%%%%
    toutp = 1;
    %%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% if t==2
%     tou_ch_add=2;
% else
%     tou_ch_add=tou_ch_add+1;
% end
% %%%
% if tou_ch_add==2
%     Q_ch_added_spatial = Q_ch_added;
% else
%     Q_ch_added_spatial = Q_ch_added_spatial + Q_ch_added;
% end
% if  mod(t,24) == 0
%     %%%%%%%%%%%
%     Title_save = strcat('OUTPUT_',TITLE_SAVE,'_Q_CH_ADD_',num2str(t));
%     save(Title_save,'Q_ch_added_spatial');
%     %%%%%%%%%%
%     tou_ch_add = 1;
%     %%%%%%%%%%
% end


%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% TRACKED PIXELS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ipo=1:npoint
    ij = sub2ind(size(DTM),Yout(ipo),Xout(ipo));
    if t==2
        tit5{1,ipo}=strcat('OUTPUT_',TITLE_SAVE,'_PIXEL_',num2str(ij),'.dat');
        fid5(1,ipo)=fopen(tit5{1,ipo},'a');
    end
    fprintf(fid5(1,ipo),'%g\t',Pr_S(ij));%1
    fprintf(fid5(1,ipo),'%g\t',Ta_S(ij));%2
    fprintf(fid5(1,ipo),'%g\t',Ws_S(ij));%3
    fprintf(fid5(1,ipo),'%g\t',Ds_S(ij));%4
    fprintf(fid5(1,ipo),'%g\t',ea_S(ij));%5
    fprintf(fid5(1,ipo),'%g\t',N_S(ij));%6
    fprintf(fid5(1,ipo),'%g\t',Pre_S(ij));%7
    fprintf(fid5(1,ipo),'%g\t',Tdew_S(ij));%8
    fprintf(fid5(1,ipo),'%g\t',Rsw_space(ij));%9
    fprintf(fid5(1,ipo),'%g\t',PAR_space(ij));%10
    fprintf(fid5(1,ipo),'%g\t',Ca_S(ij));%11
    %%%%
    fprintf(fid5(1,ipo),'%g\t',Ts(ij));%12
    fprintf(fid5(1,ipo),'%g\t',Tdamp(ij));%13
    fprintf(fid5(1,ipo),'%g\t',Csno(ij));%14
    fprintf(fid5(1,ipo),'%g\t',Cice(ij));%15
    fprintf(fid5(1,ipo),'%g\t',Csnow(ij));%16
    fprintf(fid5(1,ipo),'%g\t',Cicew(ij));%17
    fprintf(fid5(1,ipo),'%g\t',Pr_sno(ij));%18
    fprintf(fid5(1,ipo),'%g\t',Pr_liq(ij));%19
    %%%%
    fprintf(fid5(1,ipo),'%g\t',Rn(ij));%20
    fprintf(fid5(1,ipo),'%g\t',H(ij));%21
    fprintf(fid5(1,ipo),'%g\t',G(ij));%22
    fprintf(fid5(1,ipo),'%g\t',Gfin(ij));%23
    fprintf(fid5(1,ipo),'%g\t',QE(ij));%24
    fprintf(fid5(1,ipo),'%g\t',Qv(ij));%25
    fprintf(fid5(1,ipo),'%g\t',Qfm(ij));%26
    %%%%
    fprintf(fid5(1,ipo),'%g\t',SWE(ij));%27
    fprintf(fid5(1,ipo),'%g\t',SND(ij));%28
    fprintf(fid5(1,ipo),'%g\t',WR_SP(ij));%29
    fprintf(fid5(1,ipo),'%g\t',U_SWE(ij));%30
    fprintf(fid5(1,ipo),'%g\t',NIn_SWE(ij));%31
    fprintf(fid5(1,ipo),'%g\t',dw_SNO(ij));%32
    fprintf(fid5(1,ipo),'%g\t',ros(ij));%33
    fprintf(fid5(1,ipo),'%g\t',In_SWE(ij));%34
    fprintf(fid5(1,ipo),'%g\t',SP_wc(ij));%35
    fprintf(fid5(1,ipo),'%g\t',ICE(ij));%36
    fprintf(fid5(1,ipo),'%g\t',ICE_D(ij));%37
    fprintf(fid5(1,ipo),'%g\t',WR_IP(ij));%38
    fprintf(fid5(1,ipo),'%g\t',IP_wc(ij));%39
    fprintf(fid5(1,ipo),'%g\t',NIce(ij));%40
    %%%
    fprintf(fid5(1,ipo),'%g\t',EG(ij));%41
    fprintf(fid5(1,ipo),'%g\t',ESN(ij));%42
    fprintf(fid5(1,ipo),'%g\t',ESN_In(ij));%43
    fprintf(fid5(1,ipo),'%g\t',EWAT(ij));%44
    fprintf(fid5(1,ipo),'%g\t',EICE(ij));%45
    fprintf(fid5(1,ipo),'%g\t',EIn_urb(ij));%46
    fprintf(fid5(1,ipo),'%g\t',EIn_rock(ij));%47
    fprintf(fid5(1,ipo),'%g\t',SE_rock(ij));%48
    fprintf(fid5(1,ipo),'%g\t',SE_urb(ij)); %49
    fprintf(fid5(1,ipo),'%g\t',In_urb(ij));%50
    fprintf(fid5(1,ipo),'%g\t',In_rock(ij)); %51
    fprintf(fid5(1,ipo),'%g\t',WAT(ij)); %52
    fprintf(fid5(1,ipo),'%g\t',FROCK(ij)); %53
    %%%
    fprintf(fid5(1,ipo),'%g\t',f(ij));%54
    fprintf(fid5(1,ipo),'%g\t',WIS(ij));%55
    fprintf(fid5(1,ipo),'%g\t',Rd(ij));%56
    fprintf(fid5(1,ipo),'%g\t',Rh(ij));%57
    fprintf(fid5(1,ipo),'%g\t',Lk(ij));%58
    fprintf(fid5(1,ipo),'%g\t',Lk_wat(ij));%59
    fprintf(fid5(1,ipo),'%g\t',Lk_rock(ij));%60
    %%%
    fprintf(fid5(1,ipo),'%g\t',OF(ij));%61
    fprintf(fid5(1,ipo),'%g\t',OS(ij));%62
    fprintf(fid5(1,ipo),'%g\t',ZWT(ij));%63
    fprintf(fid5(1,ipo),'%g\t',q_runon(ij));%64
    fprintf(fid5(1,ipo),'%g\t',Q_channel(ij));%65
    %%%
    fprintf(fid5(1,ipo),'%g\t',r_soil(ij));%66
    fprintf(fid5(1,ipo),'%g\t',alp_soil(ij));%67
    fprintf(fid5(1,ipo),'%g\t',ra(ij));%68
    %%%%
    fprintf(fid5(1,ipo),'%g\t',er(ij));%69
    %%%
    fprintf(fid5(1,ipo),'%g\t',DQ_S(ij));%70
    fprintf(fid5(1,ipo),'%g\t',DT_S(ij));%%71
    fprintf(fid5(1,ipo),'%g\t',dQ_S(ij));%72
    %%%
    fprintf(fid5(1,ipo),'%g\t',QpointH(ipo));%73
    fprintf(fid5(1,ipo),'%g\t',QpointC(ipo));%%74
    fprintf(fid5(1,ipo),'%g\t',UpointH(ipo));%75
    fprintf(fid5(1,ipo),'%g\t',UpointC(ipo));%76
    %%%
    fprintf(fid5(1,ipo),'%g\t\n',CK1(ij));%77
    %fprintf(fid5(1,ipo),'%g\t\n',CK2(ij));%78
    if t==N_time_step
        fclose(fid5(1,ipo));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%
    if t==2
        tit6{ipo}=strcat('OUTPUT_',TITLE_SAVE,'_PIXEL_SOIL_',num2str(ij),'.dat');
        fid6(ipo)=fopen(tit6{ipo},'a');
    end
    for kjk = 1:ms_max
        fprintf(fid6(ipo),'%g\t',V(ij,kjk));
        fprintf(fid6(ipo),'%g\t',O(ij,kjk));
        fprintf(fid6(ipo),'%g\t',Tdp(ij,kjk));
        fprintf(fid6(ipo),'%g\t',Qi_in(ij,kjk));
        fprintf(fid6(ipo),'%g\t',Qi_out(ij,kjk));
    end
    fprintf(fid6(ipo),'%g\t\n',0*CK1(ij));
    if t==N_time_step
        fclose(fid6(ipo));
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%
    if t==2
        for ijki=1:cc_max
            tit7{1+ijki,ipo}=strcat('OUTPUT_',TITLE_SAVE,'_PIXEL_',num2str(ij),'_PFT_',num2str(ijki),'.dat');
            fid7(1+ijki,ipo)=fopen(tit7{1+ijki,ipo},'a');
        end
    end
    for ijki = 1:cc_max
        fprintf(fid7(1+ijki,ipo),'%g\t',T_H(ij,ijki)); %1
        fprintf(fid7(1+ijki,ipo),'%g\t',T_L(ij,ijki)); %2
        fprintf(fid7(1+ijki,ipo),'%g\t',EIn_H(ij,ijki)); %3
        fprintf(fid7(1+ijki,ipo),'%g\t',EIn_L(ij,ijki)); %4
        fprintf(fid7(1+ijki,ipo),'%g\t',Dr_H(ij,ijki));  %5
        fprintf(fid7(1+ijki,ipo),'%g\t',Dr_L(ij,ijki));  %6
        fprintf(fid7(1+ijki,ipo),'%g\t',In_H(ij,ijki));  %7
        fprintf(fid7(1+ijki,ipo),'%g\t',In_L(ij,ijki));  %8
        fprintf(fid7(1+ijki,ipo),'%g\t',Tdp_H(ij,ijki));   %9
        fprintf(fid7(1+ijki,ipo),'%g\t',Tdp_L(ij,ijki));  %10
        fprintf(fid7(1+ijki,ipo),'%g\t',OH(ij,ijki));   %11
        fprintf(fid7(1+ijki,ipo),'%g\t',OL(ij,ijki));  %12
        fprintf(fid7(1+ijki,ipo),'%g\t',rb_H(ij,ijki)); %13
        fprintf(fid7(1+ijki,ipo),'%g\t',rb_L(ij,ijki)); %14
        %%%%
        fprintf(fid7(1+ijki,ipo),'%g\t',rs_sunH(ij,ijki)); %15
        fprintf(fid7(1+ijki,ipo),'%g\t',rs_sunL(ij,ijki)); %16
        fprintf(fid7(1+ijki,ipo),'%g\t',rs_shdH(ij,ijki)); %17
        fprintf(fid7(1+ijki,ipo),'%g\t',rs_shdL(ij,ijki)); %18
        fprintf(fid7(1+ijki,ipo),'%g\t',rap_H(ij,ijki)); %19
        fprintf(fid7(1+ijki,ipo),'%g\t',rap_L(ij,ijki)); %20
        %%%
        fprintf(fid7(1+ijki,ipo),'%g\t',An_H(ij,ijki)); %21
        fprintf(fid7(1+ijki,ipo),'%g\t',An_L(ij,ijki)); %22
        fprintf(fid7(1+ijki,ipo),'%g\t',Rdark_H(ij,ijki)); %23
        fprintf(fid7(1+ijki,ipo),'%g\t',Rdark_L(ij,ijki)); %24
        %%%
        fprintf(fid7(1+ijki,ipo),'%g\t',Ci_sunH(ij,ijki)); %25
        fprintf(fid7(1+ijki,ipo),'%g\t',Ci_sunL(ij,ijki)); %26
        fprintf(fid7(1+ijki,ipo),'%g\t',Ci_shdH(ij,ijki)); %27
        fprintf(fid7(1+ijki,ipo),'%g\t',Ci_shdL(ij,ijki)); %28
        %%%
        fprintf(fid7(1+ijki,ipo),'%g\t',LAI_H(ij,ijki)); %29
        fprintf(fid7(1+ijki,ipo),'%g\t',LAI_L(ij,ijki)); %30
        fprintf(fid7(1+ijki,ipo),'%g\t',NPP_H(ij,ijki)); %31
        fprintf(fid7(1+ijki,ipo),'%g\t',NPP_L(ij,ijki)); %32
        fprintf(fid7(1+ijki,ipo),'%g\t',ANPP_H(ij,ijki)); %33
        fprintf(fid7(1+ijki,ipo),'%g\t',ANPP_L(ij,ijki)); %34
        fprintf(fid7(1+ijki,ipo),'%g\t',RA_H(ij,ijki)); %35
        fprintf(fid7(1+ijki,ipo),'%g\t',RA_L(ij,ijki)); %36
        fprintf(fid7(1+ijki,ipo),'%g\t',Rg_H(ij,ijki)); %37
        fprintf(fid7(1+ijki,ipo),'%g\t',Rg_L(ij,ijki)); %38
        fprintf(fid7(1+ijki,ipo),'%g\t',Rms_H(ij,ijki)); %39
        fprintf(fid7(1+ijki,ipo),'%g\t',Rms_L(ij,ijki)); %40
        fprintf(fid7(1+ijki,ipo),'%g\t',Rmr_H(ij,ijki)); %41
        fprintf(fid7(1+ijki,ipo),'%g\t',Rmr_L(ij,ijki)); %42
        fprintf(fid7(1+ijki,ipo),'%g\t',Rmc_H(ij,ijki)); %43
        fprintf(fid7(1+ijki,ipo),'%g\t',Rmc_L(ij,ijki)); %44
        fprintf(fid7(1+ijki,ipo),'%g\t',PHE_S_H(ij,ijki)); %45
        fprintf(fid7(1+ijki,ipo),'%g\t',PHE_S_L(ij,ijki)); %46
        fprintf(fid7(1+ijki,ipo),'%g\t',dflo_H(ij,ijki)); %47
        fprintf(fid7(1+ijki,ipo),'%g\t',dflo_L(ij,ijki)); %48
        fprintf(fid7(1+ijki,ipo),'%g\t',AgeL_H(ij,ijki)); %49
        fprintf(fid7(1+ijki,ipo),'%g\t',AgeL_L(ij,ijki)); %50
        fprintf(fid7(1+ijki,ipo),'%g\t',SAI_H(ij,ijki));  %51
        fprintf(fid7(1+ijki,ipo),'%g\t',SAI_L(ij,ijki));  %52
        fprintf(fid7(1+ijki,ipo),'%g\t',LAIdead_H(ij,ijki));  %53
        fprintf(fid7(1+ijki,ipo),'%g\t',LAIdead_L(ij,ijki));  %54
        fprintf(fid7(1+ijki,ipo),'%g\t',hc_H(ij,ijki)); %55
        fprintf(fid7(1+ijki,ipo),'%g\t',hc_L(ij,ijki)); %56
        fprintf(fid7(1+ijki,ipo),'%g\t',e_rel_L(ij,ijki)); %57
        fprintf(fid7(1+ijki,ipo),'%g\t',e_rel_H(ij,ijki)); %58
        fprintf(fid7(1+ijki,ipo),'%g\t',B_H(ij,ijki,1)); %59
        fprintf(fid7(1+ijki,ipo),'%g\t',B_H(ij,ijki,2)); %60
        fprintf(fid7(1+ijki,ipo),'%g\t',B_H(ij,ijki,3)); %61
        fprintf(fid7(1+ijki,ipo),'%g\t',B_H(ij,ijki,4)); %62
        fprintf(fid7(1+ijki,ipo),'%g\t',B_H(ij,ijki,5)); %63
        fprintf(fid7(1+ijki,ipo),'%g\t',B_H(ij,ijki,6)); %64
        fprintf(fid7(1+ijki,ipo),'%g\t',B_H(ij,ijki,7)); %65
        fprintf(fid7(1+ijki,ipo),'%g\t',B_L(ij,ijki,1)); %66
        fprintf(fid7(1+ijki,ipo),'%g\t',B_L(ij,ijki,2)); %67
        fprintf(fid7(1+ijki,ipo),'%g\t',B_L(ij,ijki,3)); %68
        fprintf(fid7(1+ijki,ipo),'%g\t',B_L(ij,ijki,4)); %69
        fprintf(fid7(1+ijki,ipo),'%g\t',B_L(ij,ijki,5)); %70
        fprintf(fid7(1+ijki,ipo),'%g\t',B_L(ij,ijki,6)); %71
        fprintf(fid7(1+ijki,ipo),'%g\t',B_L(ij,ijki,7)); %72
        fprintf(fid7(1+ijki,ipo),'%g\t\n',0*CK1(ij)); %73
    end
    if t==N_time_step
        for ijki=1:cc_max
            fclose(fid7(1+ijki,ipo));
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% TRACKING RESERVOIRS 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if not(isempty(RES_ID_List))
    for ire=1:length(RES_ID_List)
        if t==2
            tit8{1,ire}=strcat('OUTPUT_',TITLE_SAVE,'_RESERVOIRS_',num2str(RES_ID_List(ire)),'.dat');
            fid8(1,ire)=fopen(tit8{1,ire},'a');
        end
        fprintf(fid8(1,ire),'%g\t',Q_out_Res(ire));%1 %% [mm/h] 
        fprintf(fid8(1,ire),'%g\t',H_Res(ire));%2 %%[m] 
        fprintf(fid8(1,ire),'%g\t',VOL_Res(ire));%3 %% [m3] 
        fprintf(fid8(1,ire),'%g\t\n',0*Q_out_Res(ire)); %
        if t==N_time_step
            fclose(fid8(1,ire));
        end
    end
end