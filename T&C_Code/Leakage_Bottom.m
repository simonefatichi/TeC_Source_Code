%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Leakage Bottom             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Lk]=Leakage_Bottom(O,Ks_Zs,Osat,Ohy,L,alpVG,nVG,lVG,Ks_mac,Omac,alpVGM,nVGM,lVGM,Kbot,ms,SPAR)
%%REFERENCES %%
%%%INPUTS
%%% OUTPUTS
%%% Lk %% [mm/h]
%%%%%%% Kbot [mm/h] hydraulic conductivty bedrock
if isnan(Kbot)
    switch SPAR
        case 1
            Se = (O(ms)-Ohy(ms))/(Osat(ms)-Ohy(ms)); mVG= 1-1/nVG(ms);
            Ko= Ks_Zs(ms)*((Se)^(lVG(ms)))*(1-(1-(Se)^(1/mVG))^mVG)^2; %%% [mm/h]
        case 2
            Ko = Ks_Zs(ms)*(O(ms)./Osat(ms))^(3+(2/L(ms))); %%% [mm/h] %%% Estimation Conductivty last layer
        case 3
            mVG= 1-1./nVG(ms);
            mVGM= 1-1./nVGM(ms);
            Ohy= Ohy(ms); Osat=Osat(ms); Omac=Omac(ms); alpVG=alpVG(ms) ;nVG=nVG(ms); Ks_Zs=Ks_Zs(ms); lVG=lVG(ms);
            nVGM=nVGM(ms); alpVGM=alpVGM(ms); Ks_mac=Ks_mac(ms); lVGM=lVGM(ms);
            O=O(ms);
            %%%
            if Omac == 0
                Se = (O-Ohy)/(Osat-Ohy);
                Po = -(1./alpVG).*((Se).^(-1./mVG)-1).^(1./nVG); %%% [mm]
                h=-Po;
            else
                Se = (O-Ohy)/(Osat-Omac-Ohy);
                Se(Se>1)=1;
                h1 = (1/ alpVG)*((Se).^(-1/mVG)-1).^(1/nVG); %%% [mm]
                h1(Se==1 & Omac>0)=NaN; h1max = nanmax(h1);
                %%%%%%%%%%%%
                Se2 = (O-(Osat-Omac))/(Osat-(Osat-Omac)); Se2(Se2<0)=0;
                h2 = (1/alpVGM)*((Se2).^(-1/mVGM)-1).^(1/nVGM); %%% [mm]
                h2(Se2==0)=NaN;  h2(h2<h1max)=h1max;
                %%%%%%%%%
                h=[h1 ; h2]; h=nanmean(h);
                
                %%%%
                if h < - 10 && h > - 500
                    hp=-10.^(0:0.1:3);  %% [mm]
                    Op = (Omac(i))./((1+abs(alpVGM(i)*hp).^nVGM(i)).^mVGM(i)) +  Ohy(i) + (Osat(i)-Omac(i)-Ohy(i))./((1+abs(alpVG(i)*hp).^nVG(i)).^mVG(i));
                    h=interp1(Op,hp,O,'linear');
                end
            end
            Kom=Ks_Zs*(((1 - ((abs(alpVG*h)).^(nVG-1)).*(1+(abs(alpVG*h)).^nVG).^-mVG).^2)./((1+(abs(alpVG*h)).^nVG).^(lVG*mVG)));
            Kop = Ks_mac*(((1 - ((abs(alpVGM*h)).^(nVGM-1)).*(1+(abs(alpVGM*h)).^nVGM).^-mVGM).^2)./((1+(abs(alpVGM*h)).^nVGM).^(lVGM*mVGM)));
            %%%%
            Ko =  Kom + Kop ;
        case 4
            Se = (O(ms)-Ohy(ms))/(Osat(ms)-Ohy(ms)); mVG= 1-1/nVG(ms);
            Ko= Ks_Zs(ms)*((Se)^(lVG(ms)))*(1-(1-(Se)^(1/mVG))^mVG)^2; %%% [mm/h]
    end
    Lk=Ko;
else
    %Lk = exp((log(Kbot) + log(Ko))/2);% [mm/h] leakage between layer n and bedrock
    if O(ms)> Osat(ms)-1e-5
        Lk = Kbot;
    else
        Lk = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end