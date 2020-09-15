%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Temperature Root Zone 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Tdp_H,Tdp_L]=RootZone_Temp(Tdp,RfH_Zs,RfL_Zs)
[cc,n]=size(RfH_Zs); Tdp_H=zeros(1,cc); Tdp_L=zeros(1,cc); 
for i=1:cc
    Tdp_H(i) = sum(RfH_Zs(i,:).*Tdp);
    Tdp_L(i) = sum(RfL_Zs(i,:).*Tdp);
end
end 