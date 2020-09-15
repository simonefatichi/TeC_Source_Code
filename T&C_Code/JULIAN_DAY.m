%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  JULIAN_DAY              %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[jDay]= JULIAN_DAY(D)  
%%%% INPUT 
%%% Datam %% [Yr, MO, DA, HR]
% Determine the julian day of the current time
days = [31 28 31 30 31 30 31 31 30 31 30 31];    
nowYR=D(1); nowMO=D(2); nowDA=D(3); %nowHR =D(4);    
if(nowMO==1)
    jDay = nowDA;
elseif(nowMO==2)
    jDay = days(1) + nowDA;
else
    jDay = sum(days(1:(nowMO-1))) + nowDA;
    if(mod(nowYR,4)==0)
        if(mod(nowYR,400)==0)
            jDay = jDay + 1;
        elseif(mod(nowYR,100)~=0)
            jDay = jDay + 1;
        end
    end
end
return 
