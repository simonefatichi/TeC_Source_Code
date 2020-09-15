%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Compute soil parameters II  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Ofc,Oss,Owp,Ohy]=Soil_parametersII(ms,Osat,L,Pe,Ks,O33,nVG,alpVG,Kfc,Pss,Pwp,Phy,Ohy,SPAR)
%%%INPUTS
%%% Osat [] Saturation moisture 0 kPa 
%%% L % Slope of logaritimc tension-moisture curve 
%%% Pe % Tension at air antry (bubbling pressure) [kPa]
%%% Ks  % saturation conductivty [mm/h]
%%% O33 %% 33 kPa Moisture 
%%% OUTPUTS 
%%% Ofc [] Field Capacity Moisture Kns < 0.2 mm/h  
%%% Oss [] Stomatal closure begin moisture 30 kPa - 0.03 MPa - 3 m 
%%% Owp [] Stomatal closure end moisture -Wilting point 3000 kPa - 3 MPa -300 m
%%% Ohy [] Hygroscopic Moisture Evaporation cessation 10000 kPa - 10 MPa - 1000 m                                 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% COMPUTATION
%%%%  Scb wp hy 
%Psi = %[ 30 3000 10000 ]; %% [kPa]  
%Kfc = 0.2; %% [mm/h]
if nargin <= 12 
    SPAR=2; 
end 
switch SPAR 
    case 1 
         mVG= 1-1./nVG; 
         Pss = -101.9368*Pss; %%[mm]
         Se = 1./((1+abs(alpVG*Pss).^nVG).^mVG);
         Se(Se>1)=1;   O= Ohy+(Osat-Ohy).*Se;
         Oss(1,:)=O; 
         %%%%%%%%%%%
         Pwp = -101.9368*Pwp; %%[mm]
         Se = 1./((1+abs(alpVG*Pwp).^nVG).^mVG);
         Se(Se>1)=1;   O= Ohy+(Osat-Ohy).*Se;
         Owp(1,:)=O;
         %%%%%%%%%%
         Ofc = zeros(1,ms);
         Ohy = zeros(1,ms);
         %Ohy=Ohy;
    case 2 
        for i=1:ms
            B=1/L(i);
            A= exp(log(33)+B*log(O33(i))); % Coefficient of moisture tension
            %%%%%%%%%%%%%%%%%%%%
            if Pss < 33
                Oss(1,i)= O33(i)+(33-Pss)*(Osat(i)-O33(i))/(33-Pe(i));
            else
                Oss(1,i)=(Pss/A)^(-1/B);
            end
            %%%%%%%%%%%%%%%%%%%%
            if Pwp < 33
                Owp(1,i)= O33(i)+(33-Pwp)*(Osat(i)-O33(i))/(33-Pe(i));
            else
                Owp(1,i)=(Pwp/A)^(-1/B);
            end
            %%%%%%%%%%%%%%%%%%%%
            if Phy < 33
                Ohy(1,i)= O33(i)+(33-Phy)*(Osat(i)-O33(i))/(33-Pe(i));
            else
                Ohy(1,i)=(Phy/A)^(-1/B);
            end
            %%%%%%%%%%%%%%%%%%%%%%%
            Ofc(1,i)= Osat(i)*(Kfc/Ks(i))^(1/(3+(2/L(i))));
        end
end 
return 
