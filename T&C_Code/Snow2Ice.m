%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute Snow-Ice Conversion  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[ICE,ICE_D,SWE,SND,SP_wc,WAT,NIce]=Snow2Ice(dt,Ts,ICEtm1,SWEtm1,SNDtm1,SP_wctm1,WATtm1,Cicew,ros,ros_Ice_thr,Cwat,Cicewtm1,dz_ice,WatFreez_Th)
%%%INPUTS
%%% OUTPUTS
row = 1000; % water density [kg/m^3]
roi = 916.2; %% ice density [kg/m^3]
dth=dt/3600; %% [h]
%dz_ice [mm / h]  %% Conversion rate 
%ros_Ice_thr [kg/m^3]
%%%%
%%% Transition at 65m in 200 yr (Cuffrey - Peterson 2010) %%%
%%% Conversion rate: 0.0371 [mm/h]
%%%% Too abrupt transition
if ros > ros_Ice_thr
    %NIce = SWEtm1;
    NIce = min(SWEtm1,0.037*dth); %%[mm]
else
    NIce = 0;
end
ICE = ICEtm1 + NIce;  %%% [mm] Icepack After Sublimation and new
ICE_D= 0.001*ICE*(row/roi); %% New Icepack Depth [m]
SWE = SWEtm1 - NIce; %%[mm] Snow depth after compaction
%%%%
%%%%% Treating SND SP_wc ros ;
if ros>0
    del_D = 0.001*row*NIce/ros; %%% height reduction becuase of frozen snow
    SND = SNDtm1 - del_D; %%[m] new-snow
else
    SND = SNDtm1;
end
SP_wc=SP_wctm1;
%%%
%if Cicew == 1 &&  Cicewtm1 == 0
%NIce = NIce + Cwat*dz_ice*(roi/row); %%% [mm]
if Cicew == 1 
    %%% Yang et al., 2002 
    if (Ts < 0) && (SNDtm1 < 0.1)
        NIce = NIce + Cwat*dz_ice*dth; %% [mm]
    elseif (SNDtm1 > 0.1)
        NIce = NIce + Cwat*(dz_ice/5)*dth; %% [mm]
    end 
    NIce = min(NIce,WATtm1);
    ICE = ICE + NIce;  %%% [mm]
    ICE_D= 0.001*ICE*(row/roi); %% New Icepack Depth [m]
    WAT = WATtm1 - NIce ;
else
    WAT = WATtm1;
end

end
