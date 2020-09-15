%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Leakage Rock               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Lk]=Leakage_Rock(Krock,In,dth)
%%$ REFERENCES %%%%
%%%%% INPUTS
%%% Krock [mm/h] hydraulic conductivity rock 
%%% Intercepted Water [mm] 
%%% OUTPUTS
%%% Lk %% [mm/h]
if isnan(Krock) %% everything leaks 
    Lk=In/dth ;%% [mm/h]
else
    %%% Leakage from rock 
    Lk=min(Krock,In/dth); %[mm/h]  
end
Lk = Lk*dth; %%[mm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end