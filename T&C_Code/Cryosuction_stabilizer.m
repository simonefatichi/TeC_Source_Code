%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Thermal Mod %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V,O,Rd] = Cryosuction_stabilizer(Oice,V,dz,Osat,Ohy)
%%%%%%%%%%
n=length(V);
WTP= zeros(1,n); 
O = zeros(1,n); 
%%%%%%%%%%%%
for i=n:-1:1
    if i == n
        O(i) = (V(i)/dz(i))+ Ohy(i); %%%% Liquid  Water Content [] All layers
        WTP(i) = (O(i) + Oice(i) - Osat(i))*dz(i)*((O(i) + Oice(i)) > Osat(i)); %%% [mm] Water pushed up
    else
        O(i) = (V(i)+ WTP(i+1))/dz(i) + Ohy(i) ; %%%% Liquid  Water Content [] All layers
        WTP(i) = (O(i) + Oice(i) - Osat(i))*dz(i)*((O(i) + Oice(i))> Osat(i)); %%% [mm] Excess From the Reservoir - Below
    end
    if O(i) < Ohy(i)
        O(i) = Ohy(i);
    end
    if O(i) > (Osat(i)-Oice(i))
        O(i)=(Osat(i)-Oice(i));
    end
end
Rd = WTP(1); %%% [mm]  Dunne Runoff
%%%%%%%%%%%%%%%  Volume Correction for Rd and WTP
V(1) = V(1) +  WTP(2) - Rd;
V(2:n-1)= V(2:end-1)+ (WTP(3:n) - WTP(2:n-1));
V(n)= V(n) - WTP(n);
%%%%%%%%%%%%%%%
return