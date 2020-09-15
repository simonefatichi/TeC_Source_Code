%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Compute soil parameters II  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Ofc,Oss,Owp,Ohy]=Soil_parametersII_spatial(Osat,L,Pe,Ks,O33,Kfc,Pss,Pwp,Phy)
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
B=1./L;
A= exp(log(33)+B.*log(O33)); % Coefficient of moisture tension
%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%% COMPUTATION
%%%%  Scb wp hy
Psi = [Pss, Pwp, Phy]; %[ 30 3000 10000 ]; %% [kPa]
%Kfc = 0.2; %% [mm/h]
for i=1:3
    case_Psi=i;
    switch case_Psi
        case 1
            %%% enter Psi [kPa]
            if Psi(i) < 33
                Oss= O33+(33-Psi(i)).*(Osat-O33)./(33-Pe);
            else
                Oss=(Psi(i)./A).^(-1./B);
            end
        case 2
            %%% enter Psi [kPa]
            if Psi(i) < 33
                Owp= O33+(33-Psi(i)).*(Osat-O33)./(33-Pe);
            else
                Owp=(Psi(i)./A).^(-1./B);
            end
        case 3
            %%% enter Psi [kPa]
            if Psi(i) < 33
                Ohy= O33+(33-Psi(i)).*(Osat-O33)./(33-Pe);
            else
                Ohy=(Psi(i)./A).^(-1./B);
            end
    end
end
%for i=1:length(Psi)
%%% enter Psi [kPa]

%if Psi(i) < 33
%   O(i)= O33+(33-Psi(i)).*(Osat-O33)./(33-Pe);
%else
%    O(i)=(Psi(i)/A)^(-1/B);
%end
%end
%Oss=O(1); Owp=O(2); Ohy=O(3);
%%%%%%%%%%%%%%%%%%%%%%%
Ofc= Osat.*(Kfc./Ks).^(1./(3+(2./L)));
return
