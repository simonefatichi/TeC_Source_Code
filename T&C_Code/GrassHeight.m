%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Grass Height      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hc] = GrassHeight(LAI,LAIdead)
%%%% Allen et al 1989
%%% INPUT
LAI = LAI + LAIdead/3; 
%%%
if LAI  < 3.6
    hc = (LAI)/24;
else
    hc = (LAI)/24;
    %hc= exp(((LAI)-5.5)/1.5);
    hc(hc>1.2)=1.2;
    %hc= exp(((LAI)-3.4)/2.4);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
hc(hc<0.005)=0.005; 
end