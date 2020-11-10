%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Thermal Mod %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [V,O,Rd] = Cryosuction_stabilizer(Oice,V,dz,Osat,Ohy)
%%%%%%%%%%
n=length(V);
WTP= zeros(1,n); 
WTR= zeros(1,n); 
O = zeros(1,n); 
%%%%%%%%%%%%
%%%%%%%%%%%%
for i=1:1:n
    if i == 1
        O(i) = (V(i)/dz(i))+ Ohy(i); %%%% Liquid  Water Content [] All layers
        WTP(i) = (O(i) + Oice(i) - Osat(i))*dz(i)*((O(i) + Oice(i)) > Osat(i)); %%% [mm] Water pushed down 
    else
        O(i) = (V(i)+ WTP(i-1))/dz(i) + Ohy(i) ; %%%% Liquid  Water Content [] All layers
        WTP(i) = (O(i) + Oice(i) - Osat(i))*dz(i)*((O(i) + Oice(i))> Osat(i)); %%% [mm] Excess From the Reservoir - Below
    end
    if O(i) < Ohy(i)
        O(i) = Ohy(i);
    end
    if O(i) > (Osat(i)-Oice(i))
        O(i)=(Osat(i)-Oice(i));
    end
end
%%%%%%
V(1) = V(1) -  WTP(1) ;
V(2:n-1)= V(2:n-1)+ (WTP(1:n-2) - WTP(2:n-1));
V(n)= V(n) + WTP(n-1);
%%%%%%%%%%%%%%%%%% 
%%%% As the last layer could be oversaturated - we should look at water table rise 
for i=n:-1:1
    if i == n
        O(i) = (V(i)/dz(i))+ Ohy(i); %%%% Liquid  Water Content [] All layers
        WTR(i) = (O(i) + Oice(i) - Osat(i))*dz(i)*((O(i) + Oice(i)) > Osat(i)); %%% [mm] Water pushed up
    else
        O(i) = (V(i)+ WTR(i+1))/dz(i) + Ohy(i) ; %%%% Liquid  Water Content [] All layers
        WTR(i) = (O(i) + Oice(i) - Osat(i))*dz(i)*((O(i) + Oice(i))> Osat(i)); %%% [mm] Excess From the Reservoir - Below
    end
    if O(i) < Ohy(i)
        O(i) = Ohy(i);
    end
    if O(i) > (Osat(i)-Oice(i))
        O(i)=(Osat(i)-Oice(i));
    end
end
Rd = WTR(1); %%% [mm]  Dunne Runoff
%%%%%%%%%%%%%%%  Volume Correction for Rd and WTP
V(1) = V(1) +  WTR(2) - Rd;
V(2:n-1)= V(2:end-1)+ (WTR(3:n) - WTR(2:n-1));
V(n)= V(n) - WTR(n);
%%%%%%%%%%%%%%%
return
