function [SvF,Ct] = Sky_View_Factor(DTM,SLO,ASP,HZ,Z)
%%% INPUT  
%%% DTM : matrix digital elevation model  
%%% SLO slope matrix [angualar degree ]
%%% ASP aspect matrix [angular degree] from N  
%%% HZ Horizon angle  [angular degree] from N  
%%% Z Azimuth  [angular degree]  
%%%% OUTPUT 
%%% SvF Sky View Factor [0-1]
%%% Ct Terrain configuration factor [0-1] 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(DTM);  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% COMPUTATION SVF and Ct 
N= 2*pi/length(Z); 
%%%%%%%%%%%%%%%
Sl=SLO*pi/180; %%[rad] 
As=ASP*pi/180; %%%[rad] 
HZr=HZ*pi/180; %%%[rad] 
Z= Z*pi/180;  %%% [rad] 
SvF=zeros(m,n);
%%%%%%%%%%%%%%%%%%%
%bau = waitbar(0,'Internal Loop');
for i=1:m
    for j=1:n
        %waitbar(i/m,bau);
        if not(isnan(DTM(i,j)))
             HZz=squeeze(squeeze(HZr(i,j,:))); 
             %%%%%%%%%%%%%%%%%%%%
             if sum(isnan(HZz))>0
                 disp('Error in the Horizon Angle !!!') 
                 return 
             end
             SvF(i,j) = 1/(2*pi)*nansum((N)*(cos(Sl(i,j)).*(sin(HZz)).^2 + sin(Sl(i,j)).*cos(Z'-As(i,j)).*((HZz)-sin(HZz).*cos(HZz)))); 
        else
             SvF(i,j)=NaN;
        end 
    end 
end 
%close(bau);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Ct = (1 + cos(Sl))/2 - SvF; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return 
%%%%%%%%%