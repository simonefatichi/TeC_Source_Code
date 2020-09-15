%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Richards Model Integration %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[dV]=SOIL_MOISTURES_RICH_COMP(t,V,...
    Lk,f,EG,T_H,T_L,...
    Qi_in,Slo_pot,IS,SPAR,Osat,Ohy,O33,dz,Ks_Zs,Dz,numn,L,Pe,aR,aT,alpVG,nVG,lVG,Ks_mac,Omac,alpVGM,nVGM,lVGM,Phy1,s_SVG,bVG,Zs,cosalp,sinalp,SN)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Need only - Oi
%%% --> numn=length(V);
O = V'./dz + Ohy ;
I1 = O >= Osat -1e-5 ; I2 = O <= Ohy + 1e-5 ;
O(I1)=Osat(I1); O(I2)=Ohy(I2)+ 1e-5;
%O(O>Osat)=Osat; O(O<Ohy)=Ohy;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numnm1=numn-1;
vap=Osat(1:numnm1);
%%% Soil Water
dV=V*0;%zeros(numn,1);
K= Osat*0;%zeros(1,numn);
P= Osat*0;%zeros(1,numn);
Khalf=vap*0;%zeros(1,numnm1);
q=vap*0;%zeros(1,numnm1);
qf=vap*0;%zeros(1,numnm1);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
switch SPAR
    case 1
        mVG= 1-1./nVG; Se = (O-Ohy)./(Osat-Ohy);
        %%% Van-Genuchten, 1980 Corrected
        Phy = Phy1*101.9368; %% [mm]
        P1 = (1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
        P2 = -Phy*exp(Se.*-abs(bVG)); % [mm]
        P = (Se<s_SVG).*P2 +(Se>=s_SVG).*P1; % [mm]
        P(Se==0)=P2(Se==0); 
        K= Ks_Zs.*((Se).^(lVG)).*(1-(1-(Se).^(1./mVG)).^mVG).^2; %%% [mm/h]
    case 2
        B=1./L;
        A= exp(log(33)+B.*log(O33)); % Coefficient of moisture tension
        %%%%%%%%%%%%%%%%%%%%%%
        for i = 1:numn
            %[K(i),P(i)]=Conductivity_Suction(Ks_Zs(i),Osat(i),L(i),Pe(i),O33(i),O(i));
            %%%%%%%%%%%%%%%%%%%%
            K(i) = Ks_Zs(i)*(O(i)/Osat(i))^(3+(2/L(i))); %%% [mm/h]
            if O(i) < O33(i)
                P(i) = A(i).*(O(i)^-B(i)); %% [kPa]
            else
                P(i) =33-((O(i)-O33(i))*(33-Pe(i))/(Osat(i)-O33(i)));%% [kPa]
            end
            %%%%%%%%%%%%%%%
        end
        P =-101.9368*P;%%[mm]%
    case  3
        mVGM= 1-1./nVGM; mVG= 1-1./nVG;
        if Omac(1)==0
            Se = (O-Ohy)./(Osat-Ohy);
            %%% Van Genuchten  1980
            P = (1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
        else
            Se = (O-Ohy)./(Osat-Omac-Ohy);
            Se(Se>1)=1;
            h1 = (1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
            h1(Se==1 & Omac>0)=NaN; h1max = nanmax(h1);
            %%%%%%%%%%%%%
            Se2 = (O-(Osat-Omac))./(Osat-(Osat-Omac)); Se2(Se2<0)=0;
            h2 = (1./alpVGM).*((Se2).^(-1./mVGM)-1).^(1./nVGM); %%% [mm]
            h2(Se2==0)=NaN;  h2(h2<h1max)=h1max;
            %%%%%%%%%
            h=[h1 ; h2]; h=nanmean(h);
            P = h;
            
            
            %%%% Water content as superposition of matrix and macropore
            for i = 1:numn
                if P(i) < - 10 && P(i) > - 500
                    hp=-10.^(0:0.1:3);  %% [mm]
                    Op = (Omac(i))./((1+abs(alpVGM(i)*hp).^nVGM(i)).^mVGM(i)) +  Ohy(i) + (Osat(i)-Omac(i)-Ohy(i))./((1+abs(alpVG(i)*hp).^nVG(i)).^mVG(i));
                    P(i)=interp1(Op,hp,O(i),'linear');
                end
            end
            
        end
        h=P;
        
        Kom=Ks_Zs.*(((1 - ((abs(alpVG.*h)).^(nVG-1)).*(1+(abs(alpVG.*h)).^nVG).^-mVG).^2)./((1+(abs(alpVG.*h)).^nVG).^(lVG.*mVG)));
        Kop = Ks_mac.*(((1 - ((abs(alpVGM.*h)).^(nVGM-1)).*(1+(abs(alpVGM.*h)).^nVGM).^-mVGM).^2)./((1+(abs(alpVGM.*h)).^nVGM).^(lVGM.*mVGM)));
        %%%%
        K =  Kom + Kop ;
    case 4
        mVG= 1-1./nVG; Se = (O-Ohy)./(Osat-Ohy);
        %%% Van Genuchten  1980
        P = (1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
        K=  Ks_Zs.*((Se).^(lVG)).*(1-(1-(Se).^(1./mVG)).^mVG).^2; %%% [mm/h]
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
To = K*aR.*dz;  %%% Trasmissivity  unsaturated/saturated [mm^2/h]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Slo_pot: Slope Total Hydraulic Head
if SN==1
    %Qi_out = IS.*(To/aT).*cosalp.*(Slo_pot); %%%% [mm/h]
    Qi_out = IS.*(To/aT).*(sinalp); %%%% [mm/h]
else
    %Qi_out = (To/aT).*cosalp.*(Slo_pot); %%%% [mm/h]
    Qi_out = (To/aT).*(sinalp); %%%% [mm/h]
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Qi_in [mm/h] Lateral Flow incoming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%     Richards Model      %%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Khalf(1:numnm1)= 0.5*(K(1:numnm1) + K(2:numn));
%Khalf(1:numn-1)=exp(0.5*(log(K(1:numn-1)) + log(K(2:numn))));
q=Khalf.*(1*cosalp - (P(2:numn)-P(1:numnm1))./Dz(2:end)); %% [mm/h] %%% Flux positive downward from layer i+1 (above) to i (below)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
qf=(1-I1(2:numn)).*q;
q(q>0)=qf(q>0);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%% SOIL WATER BALANCE without other terms
dV(1) =  f - q(1)  - T_H(1) - T_L(1) - EG(1) + Qi_in(1) - Qi_out(1);
for i =2:numn-1
    dV(i)= q(i-1) - q(i) - T_H(i) - T_L(i) - EG(i) + Qi_in(i) - Qi_out(i) ; %
end
dV(numn) = q(numnm1) - Lk  - T_H(numn) - T_L(numn) - EG(numn) + Qi_in(numn) - Qi_out(numn);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%dV=dV';
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return