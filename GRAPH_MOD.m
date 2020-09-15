%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% GRAPHIC RESULTS HYDROLOGICAL AND VEGETATION MODEL
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
T_H = sum(T_H,2);
T_L = sum(T_L,2);
EIn_H = sum(EIn_H,2);
EIn_L = sum(EIn_L,2);
Dr_H = sum(Dr_H,2);
Dr_L = sum(Dr_L,2);
Rexmy_H = squeeze(sum(Rexmy_H,3));
Rexmy_L = squeeze(sum(Rexmy_L,3));
%OH = OH*(ones(NN,1)*(Ccrown/sum(Ccrown)));
%OL = OL*(ones(NN,1)*(Ccrown/sum(Ccrown)));
%rs_H = rs_H*((Ccrown/sum(Ccrown))');
%rs_L = rs_L*((Ccrown/sum(Ccrown))');
rb_H = rb_H*((Ccrown/sum(Ccrown))');
rb_L = rb_L*((Ccrown/sum(Ccrown))');
rap_H = rap_H*((Ccrown/sum(Ccrown))');
rap_L = rap_L*((Ccrown/sum(Ccrown))');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Area = 1;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUAL RESULT
prompt={'Precipitation [mm] ',...
    'Ground Evaporation [mm] ',...
    'Transpiration High Vegetation [mm] ',...
    'Transpiration Low Vegetation [mm] ',...
    'Sublimation/Evaporation Snow [mm] ',...
    'Saturation Excess (Dunne) [mm]',...
    'Infiltration Excess (Horton) [mm]',...
    'Evaporation from Interception [mm] ',...
    'Evaporation from Urban Soil [mm] ',...
    'Evaporation from Rocks [mm] ',...
    'Evaporation from Litter [mm] ',...
    'Recharge [mm] ',...
    'Net Lateral Subsurface [mm] ',...
    'Soil Water Content Variation [mm]',...
    'Snowpack Water Content Variation [mm]',...
    'TOTAL BALANCE [mm]'};
defo={num2str(sum((sum((Pr_liq+Pr_sno)*dth).*Area))/sum(Area)),...
    num2str(sum((sum(EG*dth).*Area))/sum(Area)),...
    num2str(sum((sum(T_H*dth).*Area))/sum(Area)),...
    num2str(sum((sum(T_L*dth).*Area))/sum(Area)),...
    num2str(sum((sum((ESN+ESN_In)*dth).*Area))/sum(Area)),...
    num2str(sum((sum(Rd).*Area))/sum(Area)),...
    num2str(sum((sum(Rh).*Area))/sum(Area)),...
    num2str(sum((sum((EIn_H+EIn_L)*dth).*Area))/sum(Area)),...
    num2str(sum((sum((EIn_urb*dth)).*Area))/sum(Area)),...
    num2str(sum((sum((EIn_rock*dth)).*Area))/sum(Area)),...
    num2str(sum((sum((ELitter*dth)).*Area))/sum(Area)),...
    num2str(sum((sum(Lk*dth).*Area))/sum(Area)),...
    num2str(sum(((sum(sum(Qi_in*dth))-sum(sum(Qi_out*dth))).*Area))/sum(Area)),...
    num2str(sum((sum(dV).*Area))/sum(Area)),...
    num2str(sum((SWE(end,:).*Area))/sum(Area) - sum((SWE(1,:).*Area))/sum(Area)),...
    num2str(sum((sum(Ck).*Area))/sum(Area))};
datiii1=inputdlg(prompt,'MAIN VARIABLES HydrologicalBalanceModel',1,defo,'on');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:length(Ccrown)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUAL RESULT II
    prompt={'GPP Gross Primary Productivity [gC/m^2] ',...
        'NPP Net Primary Productivity [gC/m^2] ',...
        'ANPP Net Primary Productivity [gC/m^2] ',...
        'RA Autrophic Respiration [gC/m^2] ',...
        'Rg Growth Respiration [gC/m^2] ',...
        'Rml Leaf Maintenance Respiration [gC/m^2] ',...
        'Rexmy Exudation and Root Export [gC/m^2] ',...
        'Rms Sapwood Maintenance Respiration [gC/m^2] ',...
        'Rms Carbohydrate Maintenance Respiration [gC/m^2] ',...
        'Rmr Root Maintenance Respiration  [gC/m^2]'};
    defo={num2str(sum((sum((NPP_H(:,i)+RA_H(:,i))*dtd).*Area))/sum(Area)),...
        num2str(sum((sum(NPP_H(:,i)*dtd).*Area))/sum(Area)),...
        num2str(sum((sum(ANPP_H(:,i)*dtd).*Area))/sum(Area)),...
        num2str(sum((sum(RA_H(:,i)*dtd).*Area))/sum(Area)),...
        num2str(sum((sum(Rg_H(:,i)*dtd).*Area))/sum(Area)),...
        num2str(sum((sum((RA_H(:,i)-Rg_H(:,i)-Rms_H(:,i)-Rmr_H(:,i)-Rmc_H(:,i))*dtd).*Area))/sum(Area)),...
        num2str(sum((sum(Rexmy_H(:,i)*dtd).*Area))/sum(Area)),...
        num2str(sum((sum(Rms_H(:,i)*dtd).*Area))/sum(Area)),...
        num2str(sum((sum(Rmc_H(:,i)*dtd).*Area))/sum(Area)),...
        num2str(sum((sum(Rmr_H(:,i)*dtd).*Area))/sum(Area))};
    tit=strcat('CARBON BALANCE HIGH-VEGETATION - CROWN ',num2str(i));
    datiii2=inputdlg(prompt,tit,1,defo,'on');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% VISUAL RESULT III
    prompt={'GPP Gross Primary Productivity [gC/m^2] ',...
        'NPP Net Primary Productivity [gC/m^2] ',...
        'ANPP Net Primary Productivity [gC/m^2] ',...
        'RA Autrophic Respiration [gC/m^2] ',...
        'Rg Growth Respiration [gC/m^2] ',...
        'Rml Leaf Maintenance Respiration [gC/m^2] ',...
        'Rexmy Exudation and Root Export [gC/m^2] ',...
        'Rms Sapwood Maintenance Respiration [gC/m^2] ',...
        'Rms Carbohydrate Maintenance Respiration [gC/m^2] ',...
        'Rmr Root Maintenance Respiration  [gC/m^2]'};
    defo={num2str(sum((sum((NPP_L(:,i)+RA_L(:,i))*dtd).*Area))/sum(Area)),...
        num2str(sum((sum(NPP_L(:,i)*dtd).*Area))/sum(Area)),...
        num2str(sum((sum(ANPP_L(:,i)*dtd).*Area))/sum(Area)),...
        num2str(sum((sum(RA_L(:,i)*dtd).*Area))/sum(Area)),...
        num2str(sum((sum(Rg_L(:,i)*dtd).*Area))/sum(Area)),...
        num2str(sum((sum((RA_L(:,i)-Rg_L(:,i)-Rms_L(:,i)-Rmr_L(:,i)-Rmc_L(:,i))*dtd).*Area))/sum(Area)),...
        num2str(sum((sum(Rexmy_L(:,i)*dtd).*Area))/sum(Area)),...
        num2str(sum((sum(Rms_L(:,i)*dtd).*Area))/sum(Area)),...
        num2str(sum((sum(Rmc_L(:,i)*dtd).*Area))/sum(Area)),...
        num2str(sum((sum(Rmr_L(:,i)*dtd).*Area))/sum(Area))};
    tit=strcat('CARBON BALANCE LOW-VEGETATION CROWN',num2str(i));
    datiii3=inputdlg(prompt,tit,1,defo,'on');
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
NN=1:NN; NNd=1:NNd;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch_met = 0;
switch_basic = 0;
switch_snow = 0;
switch_ice = 0;
switch_soil = 0;
switch_veg = 1;
switch_summary = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if switch_met == 1
    figure(1)
    subplot(3,1,1);
    set(gca,'FontSize',11);
    plot(NN,Ta(1:length(NN)),'r','LineWidth', 1.5);
    hold on; grid on;
    title('Air  Temperature ')
    ylabel('[°C]')
    subplot(3,1,2);
    plot(NN,ea(1:length(NN)),'m','LineWidth', 1.5);
    hold on; grid on;
    plot(NN,esat(1:length(NN)),':g','LineWidth', 1.5);
    title('Vapor Pressure ')
    ylabel('[Pa]')
    legend('ea','ea_{sat}')
    %datetick('x',3)
    subplot(3,1,3);
    hold on; grid on;
    plot(NN,SAB1(1:length(NN))+SAB2(1:length(NN))+...
        SAD1(1:length(NN))+SAD2(1:length(NN)),'k','LineWidth', 1.5);
    title('Incoming Short-Wave Radiation ')
    xlabel('Date'); ylabel('[W/m^2]')
    %datetick('x',3)
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
if switch_basic == 1
    figure(2)
    subplot(1,1,1);
    set(gca,'FontSize',11);
    plot(NN,OS,'b','LineWidth', 1.5);
    hold on; grid on;
    for i = 1:length(Ccrown)
        plot(NN,OH(:,i),'g','LineWidth', 1.5);
        plot(NN,OL(:,i),'--k','LineWidth', 1.5);
    end
    title(' \theta  Soil Moisture')
    ylabel('[ ]')
    legend('Evaporation Layer','H-Vegetation Layer','L-Vegetation Layer')
    figure(3)
    subplot(1,1,1);
    %plot(NN,Epot,':y','LineWidth', 1.5);
    hold on; grid on;
    plot(NN,EG,'b','LineWidth', 1.5);
    plot(NN,sum(T_H,2),'r','LineWidth', 1.5);
    plot(NN,sum(T_L,2),'m','LineWidth', 1.5);
    title('Evaporation and Transpiration')
    xlabel('Hour'); ylabel('[mm/h]')
    %legend('E-Virtual','Ground','H-Vegetation','L-Vegetation')
    legend('Ground','H-Vegetation','L-Vegetation')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(4)
    subplot(2,1,1);
    set(gca,'FontSize',11);
    plot(NN,Ts,'r','LineWidth', 1.5);
    hold on; grid on;
    plot(NN,Tdp,'b','LineWidth', 1.5);
    plot(NN,Ta(1:length(NN)),':g','LineWidth', 1.5);
    title('Temperature')
    ylabel('[°C]')
    legend('Surface Radiative','Soil at Damp. depth','Air Temperature')
    subplot(2,1,2);
    hold on; grid on;
    plot(NN,Rn,':k','LineWidth', 1.5);
    hold on; grid on;
    plot(NN,QE,'g','LineWidth', 1.5);
    plot(NN,H,'b','LineWidth', 1.5);
    plot(NN,G,'m','LineWidth', 1.5);
    plot(NN,Qv,'y','LineWidth', 1.5);
    legend('Observed','Simulated')
    title('HEAT FLUXES ')
    xlabel('Hour'); ylabel('[W/m^2]')
    legend('Rn','LE','H','G','Qv')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ra(ra>250)=NaN; r_soil(r_soil>4000)=NaN;
    rs_sunH(rs_sunH>4000)=NaN; rs_sunL(rs_sunL>4000)=NaN;
    rs_shdH(rs_shdH>4000)=NaN; rs_shdL(rs_shdL>4000)=NaN;
    rb_H(rb_H>300)=NaN; rb_L(rb_L>300)=NaN;
    rap_H(rap_H>700)=NaN; rap_L(rap_L>700)=NaN;
    figure(5)
    subplot(4,1,1);
    set(gca,'FontSize',11);
    plot(NN,ra,'r','LineWidth', 1.5);
    hold on; grid on;
    title('Aerodynamic Resistance  ra')
    subplot(4,1,2);
    set(gca,'FontSize',11);
    plot(NN,rb_H,'m','LineWidth', 1.5);
    hold on; grid on;
    plot(NN,rb_L,'k','LineWidth', 1.5);
    title('Boundary Leaf Resistance ')
    legend(' H-Vegetation rb',' L-Vegetation rb');
    ylabel('[s/m]')
    subplot(4,1,3);
    set(gca,'FontSize',11);
    hold on; grid on;
    plot(NN,rap_H,'r','LineWidth', 1.5);
    hold on; grid on;
    plot(NN,rap_L,'b','LineWidth', 1.5);
    title('Undercanopy Resistance ')
    legend('Soil rap',' L-Vegetation rap');
    xlabel('Hour'); ylabel('[s/m]')
    axis([ 0 length(NN) 0 600])
    subplot(4,1,4);
    set(gca,'FontSize',11);
    hold on; grid on;
    plot(NN,r_soil,'r','LineWidth', 1.5);
    hold on; grid on;
    title('Soil Resistance ')
    legend('Soil Resistance ');
    xlabel('Hour'); ylabel('[s/m]')
    axis([ 0 length(NN) 0 4000])
    figure(6)
    for cc=1:length(Ccrown)
        subplot(2,1,1)
        plot(NN,rs_sunH(:,cc),'b','LineWidth', 1.5);
        hold on; grid on;
        plot(NN,rs_shdH(:,cc),'--g','LineWidth', 1.5);
        ylabel('[s/m]'); xlabel('Hour');
        title('r_{s,sun}  - H Vegetation ');
        axis([ 0 length(NN) 0 1000])
        legend('H-Veg rs_{shd}','H-Veg rs_{sun}');
        subplot(2,1,2)
        plot(NN,rs_sunL(:,cc),'b','LineWidth', 1.5);
        hold on; grid on;
        plot(NN,rs_shdL(:,cc),'--g','LineWidth', 1.5);
        ylabel('[s/m]'); xlabel('Hour');
        title('r_{s,sun} - L Vegetation ');
        legend('L-Veg rs_{shd}','L-Veg rs_{sun}');
        axis([ 0 length(NN) 0 1000])
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(7)
    subplot(2,1,1);
    set(gca,'FontSize',11);
    plot(NN,In_H,'r','LineWidth', 1.5);
    hold on; grid on;
    plot(NN,In_L,'b','LineWidth', 1.5);
    legend('Interception H-Veg','Interception L-Veg')
    ylabel('mm')
    subplot(2,1,2);
    plot(NN,EIn_H,'r','LineWidth', 1.5);
    hold on; grid on;
    plot(NN,EIn_L,'b','LineWidth', 1.5);
    legend(' H-Veg',' L-Veg')
    title('Evaporation from Interception');
    xlabel('Hour'); ylabel('[mm/h]')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    nnd2=length(NNd)-2;
    T_Hd= sum(reshape(T_H(1:nnd2*24),24,nnd2));
    T_Ld= sum(reshape(T_L(1:nnd2*24),24,nnd2));
    EGd= sum(reshape(EG(1:nnd2*24),24,nnd2));
    figure(8)
    subplot(1,1,1);
    set(gca,'FontSize',11);
    plot(NNd(1:end-2),T_Hd,'r','LineWidth', 1.5);
    hold on; grid on;
    plot(NNd(1:end-2),T_Ld,'m','LineWidth', 1.5);
    title('Moisture Outfluxes')
    xlabel('days')
    plot(NNd(1:end-2),EGd,':g','LineWidth', 1.5);
    ylabel('[mm/day]')
    legend('H-Vegetation','L-Vegetation','Ground Evaporation')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
if switch_snow == 1
    figure(9)
    set(gca,'FontSize',11);
    title('Cover Fractions');
    for i = 1:length(Ccrown)
        plot([NN(1) NN(end)],[Ccrown(i) Ccrown(i)],':b','LineWidth', 1.5);
        hold on; grid on;
    end
    plot([NN(1) NN(end)],[Curb Curb],'g','LineWidth', 1.5);
    plot(NN,Csno,'c','LineWidth', 1.5);
    plot([NN(1) NN(end)],[Curb Curb],'g','LineWidth', 1.5);
    plot([NN(1) NN(end)],[Crock Crock],'y','LineWidth', 1.5);
    plot([NN(1) NN(end)],[Cwat Cwat],':b','LineWidth', 1.5);
    legend('Crown Area','','Fraction Bare Soil','Snow','Urban','Rocks','Water')
    axis([ 1 NN(end) -0.1 1.1 ])
    ylabel('[]')
    %%%%%%%%%%%%%%%%%%%%%%
end
if switch_basic == 1
    figure(10)
    set(gca,'FontSize',11);
    subplot(2,1,1);
    plot(NN,Pr_liq,'r','LineWidth', 1.5);
    hold on; grid on;
    plot(NN,Pr_sno,'b','LineWidth', 1.5);
    legend(' Liquid Pr',' Snow Pr')
    title('Precipitation');
    ylabel('[mm/h]')
    subplot(2,1,2);
    plot(NN,Dr_H,'g','LineWidth', 1.5);
    hold on; grid on;
    plot(NN,Dr_L,'m','LineWidth', 1.5);
    legend(' Drainage H-Veg',' Drainage L-Veg')
    title('Drainage');
    xlabel('Hour'); ylabel('[mm]')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(11)
    set(gca,'FontSize',11);
    subplot(2,1,1);
    plot(NN,Rh,'r','LineWidth', 1.5);
    hold on; grid on;
    plot(NN,f,'b','LineWidth', 1.5);
    plot(NN,Rd,'g','LineWidth', 1.5);
    plot(NN,sum(Qi_out-Qi_in,2),'y','LineWidth', 1.5);
    legend(' Horton Runoff',' Infiltration','Dunne Runoff','Lateral Flow')
    title('Water Fluxes');
    ylabel('[mm]')
    subplot(2,1,2);
    plot(NN,Lk,'r','LineWidth', 1.5);
    grid on;
    legend('Deep Acquifer Recharge')
    title('Groundwater');
    xlabel('Hour'); ylabel('[mm/h]')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
if switch_soil == 1
    figure(12)
    set(gca,'FontSize',11);
    for i=1:ms
        subplot(2,1,1);
        LW=1;%+m*0.2;
        plot(NN,V(:,i),'b','LineWidth', LW);
        hold on; grid on;
        title('Volumes');
        ylabel('[mm]')
        subplot(2,1,2);
        plot(NN,-Zs'*ones(1,length(NN)),'g','LineWidth', LW);
        hold on; grid on;
        title('Depths');
        xlabel('Hour'); ylabel('[mm]')
    end
    plot(NN,-ZWT,'c','LineWidth', 2);
    plot(NN,-Zdes,'k','LineWidth', 2);
    figure(13)
    set(gca,'FontSize',11);
    subplot(2,1,1)
    plot(NN,-ZWT,'b','LineWidth', 1.5);
    hold on; grid on;
    ylabel('[mm]'); xlabel('Hour');
    title('Water Table Depth');
    subplot(2,1,2)
    imagesc(O'); ylabel('Layer')
    xlabel('Time Hours');% impixelinfo ;
    colormap('hot');
    colorbar;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if switch_snow == 1
    figure(14)
    set(gca,'FontSize',11);
    subplot(3,1,1);
    plot(NN,SND,'r','LineWidth', 1.5);
    hold on; grid on;
    title('Snow Depth');
    ylabel('[m]')
    subplot(3,1,2);
    plot(NN,SWE,'g','LineWidth', 1.5);
    hold on; grid on;
    title('SWE');
    ylabel('[mm]')
    subplot(3,1,3);
    plot(NN,ros,'g','LineWidth', 1.5);
    hold on; grid on;
    title('Snow Density');
    xlabel('Hour'); ylabel('[kg/m^3]')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    figure(15)
    subplot(2,1,1);
    plot(NN,SP_wc,'g','LineWidth', 1.5);
    hold on; grid on;
    title('Water inside Snowpack');
    ylabel('[mm]')
    subplot(2,1,2);
    plot(NN,Qfm,'b','LineWidth', 1.5);
    hold on; grid on;
    title('Freezing/Melting Water in the Snowpack Heat');
    xlabel('Hour'); ylabel('[W/m^2]')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    figure(16)
    set(gca,'FontSize',11);
    subplot(3,1,1);
    plot(NN,WR_SP,'r','LineWidth', 1.5);
    hold on; grid on;
    title('Water from Snowpack');
    ylabel('[mm]')
    subplot(3,1,2);
    plot(NN,ESN,'g','LineWidth', 1.5);
    hold on; grid on;
    plot(NN,ESN_In,'b','LineWidth', 1.5);
    title('Snow Evaporation');
    legend(' Snowpack Evap.','Intercepted Snow Evap.')
    ylabel('[mm/h]')
    subplot(3,1,3);
    plot(NN,U_SWE,'k','LineWidth', 1.5);
    hold on; grid on;
    plot(NN,In_SWE,'r','LineWidth', 1.5);
    title('Intercepted Snow');
    xlabel('Hour'); ylabel('[mm]')
    legend(' Unload Snow','Intercepted Snow')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if switch_ice == 1
    figure(17)
    set(gca,'FontSize',11);
    subplot(2,1,1);
    plot(NN,ICE_D,'r','LineWidth', 1.5);
    hold on; grid on;
    title('Ice Depth');
    ylabel('[m]')
    subplot(2,1,2);
    plot(NN,ICE,'g','LineWidth', 1.5);
    hold on; grid on;
    title('ICE water');
    ylabel('[mm]')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
    figure(18)
    subplot(3,1,1);
    plot(NN,IP_wc,'g','LineWidth', 1.5);
    hold on; grid on;
    title('Water inside Icepack');
    ylabel('[mm]')
    set(gca,'FontSize',11);
    subplot(3,1,2);
    plot(NN,WR_IP,'r','LineWidth', 1.5);
    hold on; grid on;
    title('Water released from Icepack');
    ylabel('[mm]')
    subplot(3,1,3);
    plot(NN,EICE,'g','LineWidth', 1.5);
    hold on; grid on;
    title('Ice Evaporation');
    ylabel('[mm/h]')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if switch_veg == 1
    for cc=1:length(Ccrown)
        figure(100*cc + 1)
        set(gca,'FontSize',11);
        subplot(2,1,1)
        plot(NN,An_H(:,cc),'b','LineWidth', 1.5);
        hold on; grid on;
        %plot(NN,An_H(:,cc)+Rdark_H(:,cc),'g','LineWidth', 1.5);
        plot(NN,Rdark_H(:,cc),'r','LineWidth', 1.5)
        ylabel('[\mu mol CO_2 / m^2 s ]'); xlabel('Hour');
        title('Assimilation Rate H - Vegetation');
        %legend('Net Assimiliation Rate','Gross Assimilation Rate','Foliage Respiration')
        legend('Net Assimiliation Rate','Foliage Respiration')
        subplot(2,1,2)
        plot(NN,An_L(:,cc),'b','LineWidth', 1.5);
        hold on; grid on;
        %plot(NN,An_L(:,cc)+Rdark_L(:,cc),'g','LineWidth', 1.5);
        plot(NN,Rdark_L(:,cc),'r','LineWidth', 1.5);
        ylabel('[\mu mol CO_2 / m^2 s ]'); xlabel('Hour');
        title('Assimilation Rate L - Vegetation');
        %legend('Net Assimiliation Rate','Gross Assimilation Rate','Foliage Respiration')
        legend('Net Assimiliation Rate','Foliage Respiration')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(100*cc + 3)
        set(gca,'FontSize',11);
        subplot(3,1,1)
        plot(NNd,LAI_H(:,cc),'b','LineWidth', 1.5);
        hold on; grid on;
        plot(NNd,LAI_L(:,cc),'r','LineWidth', 1.5);
        plot(NNd,LAIdead_H(:,cc),'--b','LineWidth', 1.5);
        hold on; grid on;
        plot(NNd,LAIdead_L(:,cc),'--r','LineWidth', 1.5);
        ylabel('[LAI]'); xlabel('Day');
        title('LAI - Leaf Area Index');
        legend('H_Vegetation','L_Vegetation')
        subplot(3,1,2)
        plot(NNd,NPP_H(:,cc),'b','LineWidth', 1.5);
        hold on; grid on;
        plot(NNd,NPP_L(:,cc),'r','LineWidth', 1.5);
        plot(NNd,NPP_H(:,cc)+RA_H(:,cc),':b','LineWidth', 1.5);
        plot(NNd,NPP_L(:,cc)+RA_L(:,cc),':r','LineWidth', 1.5);
        ylabel('[gC / m^2 d ]'); xlabel('Day');
        title('GPP/NPP - Gross/Net Primary Productivity');
        legend('H_Vegetation','L_Vegetation')
        subplot(3,1,3)
        plot(NNd,AgeL_H(:,cc),'b','LineWidth', 1.5);
        hold on; grid on;
        plot(NNd,AgeL_L(:,cc),'r','LineWidth', 1.5);
        ylabel('[Days]'); xlabel('Day');
        title('Average Leaf Age');
        legend('H-Vegetation','L-Vegetation')
        figure(100*cc + 4)
        set(gca,'FontSize',11);
        subplot(2,1,1)
        plot(NNd,PHE_S_H(:,cc),'b','LineWidth', 1.5);
        hold on; grid on;
        ylabel('[#]');
        title('PHENOLOGY STATE');
        for k=2:length(PHE_S_H(:,cc))
            if (PHE_S_H(k,cc) == 2 && PHE_S_H(k-1,cc) == 1) || (PHE_S_H(k,cc) == 1 && PHE_S_H(k-1,cc) == 4)
                text(NNd(k),PHE_S_H(k,cc),[num2str(jDay(k))],'hor','left','vert','bottom','FontSize',8.5);
            end
        end
        legend('H-Vegetation')
        subplot(2,1,2)
        plot(NNd,PHE_S_L(:,cc),'r','LineWidth', 1.5);
        hold on; grid on;
        ylabel('[#]'); xlabel('Day');
        title('PHENOLOGY STATE');
        for k=2:length(PHE_S_L(:,cc))
            if (PHE_S_L(k,cc) == 2 && PHE_S_L(k-1,cc) == 1) || (PHE_S_L(k,cc) == 1 && PHE_S_L(k-1,cc) == 4)
                text(NNd(k),PHE_S_L(k,cc),[num2str(jDay(k))],'hor','left','vert','bottom','FontSize',8.5);
            end
        end
        legend('L-Vegetation')
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        figure(100*cc + 5)
        set(gca,'FontSize',11);
        subplot(2,1,1)
        plot(NNd,B_H(:,cc,1),'g','LineWidth', 1.5);
        hold on; grid on;
        plot(NNd,B_H(:,cc,2),'b','LineWidth', 1.5);
        plot(NNd,B_H(:,cc,3),'k','LineWidth', 1.5);
        plot(NNd,B_H(:,cc,4),'y','LineWidth', 1.5);
        title('Carbon Pool H_{VEG}')
        xlabel('Days'); ylabel('gC/m^2')
        legend('Foliage','Sapwood','Fine Roots','Carbohydrate Reserve')
        subplot(2,1,2)
        plot(NNd,B_L(:,cc,1),'g','LineWidth', 1.5);
        hold on; grid on;
        plot(NNd,B_L(:,cc,2),'b','LineWidth', 1.5);
        plot(NNd,B_L(:,cc,3),'k','LineWidth', 1.5);
        plot(NNd,B_L(:,cc,4),'y','LineWidth', 1.5);
        title('Carbon Pool L_{VEG}')
        xlabel('Days'); ylabel('gC/m^2')
        legend('Foliage','Sapwood','Fine Roots','Carbohydrate Reserve')
        set(gca,'FontSize',11);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%%%%%%%%%%%%%%%%%%%%
        figure(100*cc + 6)
        set(gca,'FontSize',11);
        subplot(2,1,1)
        plot(NNd,TdpI_H(:,cc),'b','LineWidth', 1.5);
        hold on; grid on;
        plot(NNd,TdpI_L(:,cc),'r','LineWidth', 1.5);
        ylabel('°C')
        title('Soil Temperature 30 days Average')
        subplot(2,1,2)
        hold on; grid on;
        plot(NNd,Bfac_dayH(:,cc),'g','LineWidth', 1.5);
        plot(NNd,Bfac_dayL(:,cc),':r','LineWidth', 1.5);
        xlabel('Hours'); ylabel('\beta factor []')
        legend('Root Zone H_{VEG} ','Root Zone L_{VEG} ')
        title('\beta factor daily []')
    end
    %%%%
    figure(19)
    set(gca,'FontSize',11);
    subplot(1,1,1)
    plot(NNd,(LAI_H+LAI_L)*Ccrown','k','LineWidth',1.5);
    hold on; grid on;
    ylabel('[LAI]'); xlabel('Day');
    title('LAI - Leaf Area Index');
    legend('Total weighted')
    for ij=1:length(NNd)
        if mod(ij,365)==0
            plot([ij ij],[min((LAI_H+LAI_L)*Ccrown')-0.1,max((LAI_H+LAI_L)*Ccrown')+0.1],'--y')
        end
    end
    %figure(16)
    %subplot(1,1,1)
    %plot(NN,NDVI,'g','LineWidth', 1.5);
    %hold on; grid on;
    %xlabel('Hours'); ylabel('NDVI')
    %title('NDVI')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
if switch_summary == 1
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    clear Pr_yr
    Yrs=year(Date);
    r=0;
    for i=min(Yrs):max(Yrs)
        if length(find(Yrs==i))>350*24
            r=r+1;
            Pr_yr(r)=sum(Pr(find(Yrs==i)));
            ET_yr(r)=sum(sum(T_H(find(Yrs==i),:),2)+sum(T_L(find(Yrs==i),:),2) + EG(find(Yrs==i)) +  ...
                sum(EIn_H(find(Yrs==i),:),2)+sum(EIn_L(find(Yrs==i),:),2)+EIn_urb(find(Yrs==i))+EIn_rock(find(Yrs==i)) + ...
                ESN(find(Yrs==i)) + ESN_In(find(Yrs==i)) );
            Lk_yr(r)=sum(Lk(find(Yrs==i)));
            T_yr(r)=sum(sum(T_H(find(Yrs==i),:),2)+sum(T_L(find(Yrs==i),:),2));
            VPD_yr(r)=mean(Ds(find(Yrs==i)));
        end
    end
    Yrs=year(Date(1):1:Date(end)+1);
    r=0;
    GPP_H=(NPP_H+RA_H);
    GPP_L=(NPP_L+RA_L);
    Yrs=Yrs(1:length(GPP_H));
    for i=min(Yrs):max(Yrs)
        if length(find(Yrs==i))>350
            r=r+1;
            GPP_yr(r)=  sum((GPP_H(find(Yrs==i),:)+GPP_L(find(Yrs==i),:))*Ccrown');
            NPP_yr(r)= sum((NPP_H(find(Yrs==i),:)+NPP_L(find(Yrs==i),:))*Ccrown');
            ANPP_yr(r)=  sum((ANPP_H(find(Yrs==i),:)+ANPP_L(find(Yrs==i),:))*Ccrown');
            Yrs_yr(r)= i;
        end
    end
    figure(20)
    subplot(1,2,1)
    plot(Pr_yr,ET_yr,'xk','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Pr [mm]'); ylabel('ET [mm]')
    title('Annual Value')
    subplot(1,2,2)
    plot(Pr_yr,Lk_yr,'xk','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Pr [mm]'); ylabel('Recharge [mm]')
    title('Annual Value')
    figure(18)
    subplot(1,3,1)
    plot(Pr_yr,GPP_yr,'xm','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Pr [mm]'); ylabel('GPP [gC/m^2]')
    title('Annual Value')
    subplot(1,3,2)
    plot(Pr_yr,NPP_yr,'xm','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Pr [mm]'); ylabel('NPP [gC/m^2]')
    title('Annual Value')
    subplot(1,3,3)
    plot(Pr_yr,ANPP_yr,'xm','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Pr [mm]'); ylabel('ANPP [gC/m^2]')
    title('Annual Value')
    figure(21)
    subplot(2,1,1)
    plot(Yrs_yr,GPP_yr,'-k','LineWidth', 1.5);
    hold on; grid on;
    plot(Yrs_yr,NPP_yr,'-r','LineWidth', 1.5);
    plot(Yrs_yr,ANPP_yr,'-g','LineWidth', 1.5);
    xlabel('Year'); ylabel('[gC/m^2]')
    title('Vegetation Productivities')
    legend('GPP','NPP','ANPP')
    subplot(2,1,2)
    plot(Yrs_yr,Pr_yr,'-k','LineWidth', 1.5);
    hold on; grid on;
    plot(Yrs_yr,ET_yr,'-g','LineWidth', 1.5);
    plot(Yrs_yr,T_yr,'-r','LineWidth', 1.5);
    plot(Yrs_yr,Lk_yr,'-m','LineWidth', 1.5);
    xlabel('Year'); ylabel('[mm]')
    title('Hydrologic Components')
    legend('Pr','ET','Transp.','Rech.')
    figure(22)
    %%%% Huang et al 2015 GCB
    subplot(3,2,1)
    plot(Yrs_yr,GPP_yr./ET_yr,'--k','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Year'); ylabel('[gC/m^2 mm]')
    title('EWUE = GPP/ET')
    subplot(3,2,2)
    plot(Yrs_yr,GPP_yr./T_yr,'--k','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Year'); ylabel('[gC/m^2 mm]')
    title('WUE_T = GPP/T')
    subplot(3,2,3)
    plot(Yrs_yr,0.001*GPP_yr.*VPD_yr./ET_yr,'--k','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Year'); ylabel('[gC kPa/m^2 mm]')
    title('IWUE = GPP*VPD/ET')
    subplot(3,2,4)
    plot(Yrs_yr,GPP_yr./T_yr,'--k','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Year'); ylabel('[gC/m^2 mm]')
    title('WUE_{leaf} = Ag/T')
    subplot(3,2,5)
    plot(Yrs_yr,ANPP_yr./Pr_yr,'--k','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Year'); ylabel('[gC/m^2 mm]')
    title('RAIN USE EFFICIENCY (Huxman et al., 2004)')
    subplot(3,2,6)
    plot(Yrs_yr,ET_yr./(ET_yr+Lk_yr),'--k','LineWidth', 1.5);
    hold on; grid on;
    xlabel('Year'); ylabel('[-]')
    title('HORTON INDEX')
    
    
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% EG T_H T_L ESN ESN_In Rd Rh EIn_H+EIn_L+EIn_urb+EIn_rock Qsub Qi
    Nyr_plot = length(Date)/8766;
    clear Tm ESNm EGm Einm Qim Rhm Rdm Prm Lkm TT
    PRECIP = (Pr_liq+Pr_sno)*dth;
    TRASP = (T_H + T_L)*dth;
    EINT = (ELitter+EIn_H+EIn_L+EIn_urb+EIn_rock)*dth;
    QPER = sum(Qi_out-Qi_in,2) ;
    ESNOW = (EICE+ESN + ESN_In)*dth;
    %ET=ET*dth;
    for j=1:12
        %%%%%%%%%%%%%%%%%%%%%%%%%
        Tm(j)=sum(TRASP(Datam(:,2)==j));
        ESNm(j)=sum(ESNOW(Datam(:,2)==j));
        EGm(j)=sum(EG(Datam(:,2)==j));
        Einm(j)=sum(EINT(Datam(:,2)==j));
        Qim(j)=sum(QPER(Datam(:,2)==j));
        Rhm(j)=sum(Rh(Datam(:,2)==j));
        Rdm(j)=sum(Rd(Datam(:,2)==j));
        Prm(j)=sum(PRECIP(Datam(:,2)==j));
        Lkm(j)=sum(Lk(Datam(:,2)==j));
        %ETm(j)=sum(ET(Datam(:,2)==j));
        %Qsubm(j)= Prm(j)-Tm(j)-EGm(j) - ESNm(j)- Einm(j) -Qim(j) -Rhm(j) - Rdm(j); %%%
        TT(j)=Tm(j)+ESNm(j)+EGm(j)+Einm(j)+Qim(j)+Lkm(j)+Rhm(j)+Rdm(j);
        
        %%%%%%%%%%%%%%%%%%%%%%
    end
    figure(1001)
    set(gca,'FontSize',9);
    mmm=[0,1:12,13];
    fill(mmm,[0 Rhm+Rdm+Qim+ESNm+Einm+Lkm+Tm+EGm 0]/Nyr_plot,'y')
    hold on ;
    fill(mmm,[0 Rhm+Rdm+Qim+ESNm+Einm+Lkm+Tm 0]/Nyr_plot,'r')
    fill(mmm,[0 Rhm+Rdm+ESNm+Einm+Lkm+Qim 0]/Nyr_plot,'g')
    fill(mmm,[0 Rhm+Rdm+ESNm+Einm 0]/Nyr_plot,'c')
    fill(mmm,[0 Rhm+Rdm 0]/Nyr_plot,'b')
    plot(1:12,Prm/Nyr_plot,'o--b','LineWidth', 1.5);
    %plot(1:12,ETm/Nyr_plot,'o--r','LineWidth', 1.5);
    grid on;
    xlim([1 12])
    %axis([1 12 0 max(Prm/Nyr_plot)+12])
    legend('Soil Evap.','Transp.','Rec. + Lat. Flow','Interc. Evap.','Runoff','Precip.')
    ylabel('[mm]'); xlabel('Month')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1002)
    set(gca,'FontSize',9);
    mmm=[0,1:12,13];
    fill(mmm,[0 Rhm+Rdm+Qim+ESNm+Einm+Lkm+Tm+EGm 0]/Nyr_plot,'y')
    hold on ;
    fill(mmm,[0 Rhm+Rdm+Lkm+Qim 0]/Nyr_plot,'r')
    fill(mmm,[0 Rhm+Rdm+Qim 0]/Nyr_plot,'b')
    fill(mmm,[0 Rhm+Rdm 0]/Nyr_plot,'g')
    fill(mmm,[0 Rhm 0]/Nyr_plot,'c')
    plot(1:12,Prm/Nyr_plot,'o--b','LineWidth', 1.5);
    grid on;
    xlim([1 12])
    %axis([1 12 0 max(Prm/Nyr_plot)+12])
    legend('Evapotransp.','Rec.','Lat. flow','Sat. Excess','Infilt. Excess','Precip.')
    ylabel('[mm]'); xlabel('Month')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1003)
    set(gca,'FontSize',9);
    mmm=[0,1:12,13];
    fill(mmm,[0 (Rhm+Rdm+Qim+ESNm+Einm+Lkm+Tm+EGm)./TT 0],'y')
    hold on ;
    fill(mmm,[0 (Rhm+Rdm+Qim+ESNm+Einm+Lkm+Tm)./TT 0],'r')
    fill(mmm,[0 (Rhm+Rdm+ESNm+Einm+Lkm+Qim)./TT 0],'g')
    fill(mmm,[0 (Rhm+Rdm+ESNm+Einm)./TT 0],'c')
    fill(mmm,[0 (Rhm+Rdm)./TT 0],'b')
    grid on;
    axis([1 12 0 1])
    legend('Soil Evap.','Transp.','Rec. + Lat. Flow','Interc. Evap.','Runoff')
    ylabel('Fraction'); xlabel('Month')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    figure(1004)
    set(gca,'FontSize',9);
    mmm=[0,1:12,13];
    fill(mmm,[0 (Rhm+Rdm+Qim+ESNm+Einm+Lkm+Tm+EGm)./TT 0],'y')
    hold on ;
    fill(mmm,[0 (Rhm+Rdm+Lkm+Qim)./TT 0],'r')
    fill(mmm,[0 (Rhm+Rdm+Qim)./TT 0],'b')
    fill(mmm,[0 (Rhm+Rdm)./TT 0],'g')
    fill(mmm,[0 (Rhm)./TT 0],'c')
    grid on;
    axis([1 12 0 1])
    legend('Evapotransp.','Rec.','Lat. flow','Sat. Excess','Infilt. Excess')
    ylabel('Fraction'); xlabel('Month')
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
end



