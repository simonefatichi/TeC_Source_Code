%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  FRUIT HARVEST         %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[B,RBi]= Fruit_Harvest(B,dtd,jDay,jDay_harv,B_harv)
%%%% INPUT
Btm1=B;
%%%%%%%%%% Harvest
if  not(isempty(intersect(jDay,jDay_harv)))
    if B_harv == -1 %%% Full harvesting of current fruits
        B(5)=0; 
    else
        B(5)=B(5) - B_harv; %%% Harvesting [gC/ m^2 ]
    end
    B(B<0)=0;
else
    %%% No  harvest
    B=B;
end
RBi=(Btm1(5)-B(5))/dtd; %%  [gC/ m^2 day]
return