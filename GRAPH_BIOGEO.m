%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


figure(106)
subplot(3,1,1)
plot(P(:,1),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[gC m^{-2}]')
plot(P(:,2),'g','LineWidth', 1.5);
plot(P(:,3),'r','LineWidth', 1.5);
legend('Ab. Met','Ab Str','Ab Str Lig')
subplot(3,1,2)
plot(P(:,6),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[gC m^{-2}]')
plot(P(:,7),'g','LineWidth', 1.5);
plot(P(:,8),'r','LineWidth', 1.5);
legend('Be. Met','Be Str','Be Str Lig')
subplot(3,1,3)
plot(P(:,4),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[gC m^{-2}]')
plot(P(:,5),'g','LineWidth', 1.5);
legend('Ab. Wood','Ab Wood Lign')
xlabel('Days');

figure(107)
subplot(2,2,1)
plot(P(:,9),'g','LineWidth', 1.5);
hold on; grid on;
plot(P(:,10),'k','LineWidth', 1.5);
plot(P(:,11),'y','LineWidth', 1.5);
legend('SOM POC Lign','SOM POC - Cell','SOM MOC')
title('Carbon Pool'); ylabel('[gC m^{-2}]')
subplot(2,2,2)
plot(P(:,12),'k','LineWidth', 1.5);
hold on; grid on;
plot(P(:,13),'r','LineWidth', 1.5);
legend('DOC-B','DOC-F')
title('Carbon Pool') ;ylabel('[gC m^{-2}]')
subplot(2,2,3)
plot(P(:,18),'g','LineWidth', 1.5);
hold on; grid on;
plot(P(:,19),'k','LineWidth', 1.5);
plot(P(:,20),'r','LineWidth', 1.5);
plot(P(:,21),'b','LineWidth', 1.5);
plot(P(:,22),'y','LineWidth', 1.5);
xlabel('Days'); ylabel('[gC m^{-2}]')
legend('Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza','Macrofauna')
subplot(2,2,4)
plot(P(:,14),'g','LineWidth', 1.5);
hold on; grid on;
plot(P(:,15),'k','LineWidth', 1.5);
plot(P(:,16),'r','LineWidth', 1.5);
plot(P(:,17),'b','LineWidth', 1.5);
xlabel('Days'); ylabel('[gC m^{-2}]')
legend('EM-POC-B','EM-POC-F','EM-MOC-B','EM-MOC-F')

figure(108)
subplot(2,2,1)
plot(P(:,23),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[gN m^{-2}]')
plot(P(:,24),'g','LineWidth', 1.5);
plot(P(:,25),'r','LineWidth', 1.5);
legend('Ab. Lit','Ab Wod','Be Lit')
title('Nitrogen Pool')
subplot(2,2,2)
plot(P(:,26),'g','LineWidth', 1.5);
hold on; grid on;
 ylabel('[gN m^{-2}]')
 title('Nitrogen Pool')
legend('SOM')
subplot(2,2,3)
plot(P(:,27),'g','LineWidth', 1.5);
hold on; grid on;
plot(P(:,28),'k','LineWidth', 1.5);
plot(P(:,29),'r','LineWidth', 1.5);
plot(P(:,30),'b','LineWidth', 1.5);
xlabel('Days'); ylabel('[gN m^{-2}]')
legend('Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza')
subplot(2,2,4)
plot(P(:,31),'g','LineWidth', 1.5);
hold on; grid on;
plot(P(:,32),'k','LineWidth', 1.5);
plot(P(:,33),'b','LineWidth', 1.5);
xlabel('Days');  ylabel('[gN m^{-2}]')
legend('NH4+ ','NO3-','DON')

figure(109)
subplot(3,2,1)
plot(P(:,35),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[gP m^{-2}]')
plot(P(:,36),'g','LineWidth', 1.5);
plot(P(:,37),'r','LineWidth', 1.5);
legend('Ab. Lit','Ab Wod','Be Lit')
title('Phosporus Pool')
subplot(3,2,2)
plot(P(:,38),'g','LineWidth', 1.5);
hold on; grid on; ylabel('[gP m^{-2}]')
legend('SOM')
title('Phosporus Pool')
subplot(3,2,3)
plot(P(:,39),'g','LineWidth', 1.5);
hold on; grid on;
plot(P(:,40),'k','LineWidth', 1.5);
plot(P(:,41),'r','LineWidth', 1.5);
plot(P(:,42),'b','LineWidth', 1.5);
ylabel('[gP m^{-2}]')
legend('Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza')
subplot(3,2,4)
plot(P(:,43),'g','LineWidth', 1.5);
hold on; grid on;
plot(P(:,47),'r','LineWidth', 1.5);
ylabel('[gP m^{-2}]')
legend('Mineral','DOP')
subplot(3,2,5)
plot(P(:,44),'k','LineWidth', 1.5);
hold on; grid on;
title('Phosporus Pool')
xlabel('Days'); ylabel('[gP m^{-2}]')
legend('Primary Material')
subplot(3,2,6)
plot(P(:,46),'k','LineWidth', 1.5);
hold on; grid on;
plot(P(:,45),'g','LineWidth', 1.5);
xlabel('Days');ylabel('[gP m^{-2}]')
legend('Occluded','Secondary')


figure(110)
subplot(3,2,1)
plot(P(:,48),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[gK m^{-2}]')
plot(P(:,49),'g','LineWidth', 1.5);
plot(P(:,50),'r','LineWidth', 1.5);
legend('Ab. Lit','Ab Wod','Be Lit')
title('Potassium Pool')
subplot(3,2,2)
plot(P(:,51),'g','LineWidth', 1.5);
hold on; grid on; ylabel('[gK m^{-2}]')
legend('SOM')
title('Potassium Pool')
subplot(3,2,3)
plot(P(:,52),'g','LineWidth', 1.5);
hold on; grid on;
ylabel('[gK m^{-2}]')
legend('Mineral Solution ')
subplot(3,2,4)
plot(P(:,53),'g','LineWidth', 1.5);
hold on; grid on;
plot(P(:,54),'k','LineWidth', 1.5);
xlabel('Days'); ylabel('[gK m^{-2}]')
legend('Excheangeable ','Non-Excheangeable')
subplot(3,2,5)
plot(P(:,55),'g','LineWidth', 1.5);
xlabel('Days');
hold on; grid on; ylabel('[gK m^{-2}]')
legend('Primary Minerals')



figure(111)
plot((P(:,2)+P(:,3)+P(:,1))./P(:,23),'g','LineWidth', 1.5);
hold on; grid on;
plot((P(:,5)+P(:,4))./P(:,24),'m','LineWidth', 1.5);
plot((P(:,6)+P(:,7)+P(:,8))./P(:,25),'k','LineWidth', 1.5);
plot((P(:,9)+P(:,10)+P(:,11))./P(:,26),'b','LineWidth', 1.5);
plot(P(:,18)./P(:,27),'r','LineWidth', 1.5);
plot(P(:,19)./P(:,28),'y','LineWidth', 1.5);
plot(P(:,20)./P(:,29),'c','LineWidth', 1.5);
plot(P(:,21)./P(:,30),'Color',[0.168 0.50586 0.3372],'LineWidth', 1.5);
plot(P(:,22)./P(:,34),'Color',[0.06 0.7 0.6],'LineWidth', 1.5);
title('C:N Ratio')
xlabel('Days'); ylabel('C:N')
legend('AG Litter','AG Wood','BG Litter','SOM','Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza','Macrofauna')
%%%%%%%%%%%%%%%%%%%%%%%%%


figure(112)
plot((P(:,2)+P(:,3)+P(:,1))./P(:,35),'g','LineWidth', 1.5);
hold on; grid on;
plot((P(:,5)+P(:,4))./P(:,36),'m','LineWidth', 1.5);
plot((P(:,6)+P(:,7)+P(:,8))./P(:,37),'k','LineWidth', 1.5);
plot((P(:,9)+P(:,10)+P(:,11))./P(:,38),'b','LineWidth', 1.5);
plot(P(:,18)./P(:,39),'r','LineWidth', 1.5);
plot(P(:,19)./P(:,40),'y','LineWidth', 1.5);
plot(P(:,20)./P(:,41),'c','LineWidth', 1.5);
plot(P(:,21)./P(:,42),'Color',[0.168 0.50586 0.3372],'LineWidth', 1.5);
title('C:P Ratio')
xlabel('Days'); ylabel('C:P')
legend('AG Litter','AG Wood','BG Litter','SOM','Bacteria','Fungi','AM-Mycorrhiza','EM-Mycorrhiza')
%%%%%%%%%%%%%%%%%%%%%%%%%


figure(113)
plot((P(:,2)+P(:,3)+P(:,1))./P(:,48),'g','LineWidth', 1.5);
hold on; grid on;
plot((P(:,5)+P(:,4))./P(:,50),'m','LineWidth', 1.5);
plot((P(:,6)+P(:,7)+P(:,8))./P(:,37),'k','LineWidth', 1.5);
plot((P(:,9)+P(:,10)+P(:,11))./P(:,51),'b','LineWidth', 1.5);
title('C:K Ratio')
xlabel('Days'); ylabel('C:K')
legend('AG Litter','AG Wood','BG Litter','SOM')
%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
figure(114)
subplot(2,2,1)
plot(R_litter_sur,'r','LineWidth', 1.5);
hold on; grid on;
plot(R_litter-R_litter_sur,'m','LineWidth', 1.5);
plot(R_microbe,'b','LineWidth', 1.5);
title('Respiration Het.')
plot(R_ew,'g','LineWidth', 1.5);
xlabel('Days'); ylabel('[gC m^{-2} day^{-1}]')
legend('Litter AG','Litter BG','Microbe','Macrofauna')
subplot(2,2,2)
plot(VOL,'r','LineWidth', 1.5);
hold on; grid on;
plot(N2flx,'b','LineWidth', 1.5);
title('N - Fluxes')
xlabel('Days'); ylabel('[gN m^{-2} day^{-1}]')
legend('NH_4 Vol.','N_2')
subplot(2,2,3)
plot(Min_N,'k','LineWidth', 1.5);
hold on ;  grid on 
title('N - Fluxes')
xlabel('Days'); ylabel('[gN m^{-2} day^{-1}]')
legend('Min-N')
subplot(2,2,4)
plot(Min_P,'k','LineWidth', 1.5);
hold on ;  grid on 
title('P - Fluxes')
xlabel('Days'); ylabel('[gP m^{-2} day^{-1}]')
legend('Min-P')
%%%%%%%%%%%%%%%%%%%%%%%%%


figure(115)
subplot(3,1,1)
plot(Nuptake_H+Nuptake_L,'r','LineWidth', 1.5);
hold on; grid on;
title('N Uptake')
ylabel('[gN m^{-2} day^{-1}]')
subplot(3,1,2)
plot(Puptake_H+Puptake_L,'g','LineWidth', 1.5);
hold on; grid on;
title('P Uptake')
ylabel('[gP m^{-2} day^{-1}]')
subplot(3,1,3)
plot(Kuptake_H+Kuptake_L,'m','LineWidth', 1.5);
hold on; grid on;
title('K Uptake')
xlabel('Days'); ylabel('[gK m^{-2} day^{-1}]')
%%%%%%%%%%%%%%%%%%%%%%%%%

figure(116)
subplot(2,2,1)
plot(LEAK_DOC,'b','LineWidth', 1.5);
hold on; grid on;
title('DOC Leaching')
xlabel('Days'); ylabel('[gC m^{-2} day^{-1}]')
subplot(2,2,2)
plot(LEAK_NH4,'r','LineWidth', 1.5);
hold on; grid on;
plot(LEAK_NO3,'b','LineWidth', 1.5);
hold on; grid on;
plot(LEAK_DON,'m','LineWidth', 1.5);
title('N Leaching')
xlabel('Days'); ylabel('[gN m^{-2} day^{-1}]')
legend('NH_4','NO_3','DON')
subplot(2,2,3)
plot(LEAK_P,'g','LineWidth', 1.5);
hold on; grid on;
plot(LEAK_DOP,'m','LineWidth', 1.5);
title('P Leaching')
legend('P','DOP')
xlabel('Days'); ylabel('[gP m^{-2} day^{-1}]')
subplot(2,2,4)
plot(LEAK_K,'m','LineWidth', 1.5);
hold on; grid on;
title('K Leaching')
xlabel('Days'); ylabel('[gK m^{-2} day^{-1}]')
%%%%%%%%%%%%%%%%%%%%%%%%%


figure(117)
subplot(2,3,1)
plot(Nreserve_H,'b','LineWidth', 1.5);
hold on; grid on;
plot(Nreserve_L,'r','LineWidth', 1.5);
title('N Reserve')
 ylabel('[gN m^{-2}]')
legend('H-Veg','L-Veg')
subplot(2,3,2)
plot(Preserve_H,'b','LineWidth', 1.5);
hold on; grid on;
plot(Preserve_L,'r','LineWidth', 1.5);
hold on; grid on;
legend('H-Veg','L-Veg')
title('P Reserve')
 ylabel('[gP m^{-2}]')
subplot(2,3,3)
plot(Kreserve_H,'b','LineWidth', 1.5);
hold on; grid on;
plot(Kreserve_L,'r','LineWidth', 1.5);
hold on; grid on;
legend('H-Veg','L-Veg')
title('K Reserve')
 ylabel('[gK m^{-2}]')
subplot(2,3,4)
plot(TNIT_H,'b','LineWidth', 1.5);
hold on; grid on;
plot(TNIT_L,'r','LineWidth', 1.5);
title('Total N')
xlabel('Days'); ylabel('[gN m^{-2}]')
legend('H-Veg','L-Veg')
subplot(2,3,5)
plot(TPHO_H,'b','LineWidth', 1.5);
hold on; grid on;
plot(TPHO_L,'r','LineWidth', 1.5);
hold on; grid on;
legend('H-Veg','L-Veg')
title('Total P')
xlabel('Days'); ylabel('[gP m^{-2}]')
subplot(2,3,6)
plot(TPOT_H,'b','LineWidth', 1.5);
hold on; grid on;
plot(TPOT_L,'r','LineWidth', 1.5);
hold on; grid on;
legend('H-Veg','L-Veg')
title('Total K')
xlabel('Days'); ylabel('[gK m^{-2}]')


figure(118)
subplot(3,1,1)
plot(rNc_H,'b','LineWidth', 1.5);
hold on; grid on;
plot(rNc_L,'r','LineWidth', 1.5);
title('N Conc')
xlabel('Days'); ylabel('[-]')
legend('H-Veg','L-Veg')
subplot(3,1,2)
plot(rPc_H,'b','LineWidth', 1.5);
hold on; grid on;
plot(rPc_L,'r','LineWidth', 1.5);
hold on; grid on;
legend('H-Veg','L-Veg')
title('P Conc.')
xlabel('Days'); ylabel('[-]')
subplot(3,1,3)
plot(rKc_H,'b','LineWidth', 1.5);
hold on; grid on;
plot(rKc_L,'r','LineWidth', 1.5);
hold on; grid on;
legend('H-Veg','L-Veg')
title('K Conc.')
xlabel('Days'); ylabel('[-]')

n=length(Lk); fr=24;
m=floor(n/fr);
Lkday=reshape(Lk(1:m*fr),fr,m);
Lkday=sum(Lkday);Lkday=Lkday';
ff=length(Lkday); 

figure(120)
subplot(2,2,1)
plot(LEAK_DOC(1:ff)./(Lkday)*1000,'b','LineWidth', 1.5);
hold on; grid on;
title('DOC Conc.')
xlabel('Days'); ylabel('[mg-C/l]')
subplot(2,2,2)
plot(LEAK_NH4(1:ff)./(Lkday)*1000,'r','LineWidth', 1.5);
hold on; grid on;
plot(LEAK_NO3(1:ff)./(Lkday)*1000,'b','LineWidth', 1.5);
hold on; grid on;
plot(LEAK_DON(1:ff)./(Lkday)*1000,'m','LineWidth', 1.5);
title('N Conc.')
xlabel('Days'); ylabel('[mg-N/l]')
legend('NH_4','NO_3','DON')
subplot(2,2,3)
plot(LEAK_P(1:ff)./(Lkday)*1000*1000,'g','LineWidth', 1.5);
hold on; grid on;
plot(LEAK_DOP(1:ff)./(Lkday)*1000*1000,'m','LineWidth', 1.5);
title('P Conc.')
legend('P','DOP')
xlabel('Days');  ylabel('[ug-P/l]')
subplot(2,2,4)
plot(LEAK_K(1:ff)./(Lkday)*1000,'m','LineWidth', 1.5);
hold on; grid on;
title('K Conc.')
xlabel('Days'); ylabel('[mg-K/l]')


figure(121)
subplot(3,2,1)
plot((P(:,18)+P(:,19)+P(:,20) + P(:,21))./(P(:,9)+P(:,10)+P(:,11)+P(:,12)+P(:,13)),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[-]')
title('Microbial/Substrate')
subplot(3,2,2)
plot((P(:,22))./(P(:,18)+P(:,19)+P(:,20) + P(:,21)),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[-]')
title('Macrofauna/Microbial')
subplot(3,2,3)
plot(sum(P(:,6:21),2)./(Zbio*sum(Bio_Zs.*rsd)),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[gC /g soil]'); title('SOC')
subplot(3,2,4)
plot(sum(P(:,6:21),2),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[gC /m^2]'); title('SOC')
subplot(3,2,5)
plot(((P(:,19) + P(:,20) + P(:,21))./P(:,18)),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[-]')
title('Fungi/Bacteria')
xlabel('Days');
subplot(3,2,6)
plot( P(:,19)./(P(:,20)+P(:,21)),'b','LineWidth', 1.5);
hold on; grid on;  ylabel('[-]')
title('Saprotrophic/Mycorrhiza')
xlabel('Days');

Rlitter_sub= R_litter-R_litter_sur; 
rb = R_bacteria./R_microbe; 
rf = 1 -rb; 
Rmr= ((Rmr_L+Rmr_H)*Ccrown');
TSR = Rlitter_sub + R_microbe + R_ew + Rmr; 
figure(122)
subplot(2,1,1)
plot((Rmr)./TSR,'r','LineWidth', 1.5);
hold on; grid on;
plot((R_bacteria + rb.*Rlitter_sub)./TSR,'g','LineWidth', 1.5);
hold on; grid on;
plot((rf.*R_microbe + rf.*Rlitter_sub)./TSR,'b','LineWidth', 1.5);
title('Respiration A+H Subsurface')
plot(R_ew./TSR,'y','LineWidth', 1.5);
xlabel('Days'); ylabel('[%]')
legend('Root','Bacteria','Fungi','Macrofauna')
subplot(2,1,2)
plot((1e+6/24)*(Rlitter_sub+R_microbe+R_ew+Rmr)./sum(P(:,6:21),2),'r','LineWidth', 1.5);
hold on; grid on;
xlabel('Days'); ylabel('[ug CO_{2}-C / h gC-SOC]')
title('Respiration Soil A+H')


%%%%%%%% SOIL BIOGEOCHEMISTRY BALANCE CHECK
P(P<0)=0; 
if R_microbe(end) == 0
    ed=length(R_microbe)-1;
else
    ed=length(R_microbe);
end
IS= Ccrown*squeeze(sum(ISOIL_L,1)) + Ccrown*squeeze(sum(ISOIL_H,1));
C_exp = sum(IS(1:9));
N_exp = sum(IS(10:12));
P_exp = sum(IS(13:15));
K_exp = sum(IS(16:18));
dP_soil = P(1,:) - P(ed,:);
C_out = sum(LEAK_DOC)+sum(R_litter)+sum(R_microbe)+sum(R_ew);
N_out = sum(sum(Nuptake_H + Nuptake_L)*Ccrown') + sum(LEAK_NH4) + sum(LEAK_NO3) + sum(LEAK_DON) + sum(VOL) +sum(N2flx);
P_out = sum(sum(Puptake_H + Puptake_L)*Ccrown') + sum(LEAK_DOP) + sum(LEAK_P);
K_out = sum(sum(Kuptake_H + Kuptake_L)*Ccrown') +sum(LEAK_K);

%%%% Is missing the contribution from fertilization 
nttsp = length(VOL); 

CkC_s = sum(dP_soil([1:22])) + C_exp - C_out;
CkN_s= sum(dP_soil([23:34]))+ N_exp - N_out +B_IO.DepN*nttsp;
CkP_s= sum(dP_soil([35:47]))+ P_exp - P_out+(B_IO.DepP+B_IO.Tup_P)*nttsp;
CkK_s= sum(dP_soil(48:55))+ K_exp - K_out+(B_IO.DepK+B_IO.Tup_K)*nttsp;





