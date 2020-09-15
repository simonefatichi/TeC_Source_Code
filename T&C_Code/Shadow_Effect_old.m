function [ShF] = Shadow_Effect(DTM,h_S,zeta_S,HZ,Z)
%%% INPUT  
%%% DTM : matrix digital elevation model  
%%% h_S solar height [rad] 
%%% zeta_S solar azimuth [rad] 
%%% HZ Horizon angle  [angular degree] from N 
%%% Z Azimuth  [angular degree] from N  
%%%% OUTPUT 
%%% ShF Shadow Effect 1 --> no shadow 0 --> shadow 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% COMPUTATION ShF 
[m,n]=size(DTM);  
S_azi=pi/2-h_S; %%%[rad] 
%%%%%%%%%%%%%%%
HZr=HZ*pi/180; %%[rad] 
Z= Z*pi/180; %%% [rad] 
ShF=zeros(m,n);
%%%%%%%%%%%%%%%%%%%
for i=1:m
    for j=1:n
        if not(isnan(DTM(i,j)))
             HZz=squeeze(squeeze(HZr(i,j,:))); 
             HZi = interp1([Z 2*pi],[HZz ; HZz(1)],zeta_S);
             %%%%%%%%%%%%%%%%%%%%
             if HZi > S_azi
                 ShF(i,j)=1;  %% No Shadow
             else
                 ShF(i,j)=0;  %% Shadow 
             end 
             %%%%%%%%%%%%%%%%%%%%%%
             if sum(isnan(HZz))>0
                 disp('I break !!!') 
                 return 
             end
        else
             ShF(i,j)=NaN;
        end 
    end 
end 
%%%%%
return 