function [HZ,Z] = Horizon_angle_polar(DTM,cellsize,dmax)
%%% INPUT  
%%% DTM : matrix digital elevation model [m]  
%%% cellsize dimension cell [m] 
%%% dmax [m] max search distance 
%%%% OUTPUT 
%%% HZ Horizon angle array [angular degree]   
%%% Z Azimuth directions  [angualr degree] from N  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(DTM);  
xllcorner=0; 
yllcorner=0; 
x=xllcorner:cellsize:(xllcorner+cellsize*(n-1));
y=yllcorner:cellsize:(yllcorner+cellsize*(m-1));
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin == 2 
    dmax = max(max(x),max(y)); 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
dthe=30; 
Z= 0:dthe:360-dthe; %% search angle  [angular degree]   
spass = sqrt(2)*cellsize;  %% pass of search 
SD = 0.1:spass:dmax; %%[m] search distance 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% COMPUTATION HORIZON ANGLE f(azimuth) 
HZ=zeros(m,n,length(Z)); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%bau = waitbar(0,'Internal Loop');
for i=1:m %% row
    for j=1:n %%% column
        %waitbar(i/m,bau);
        if not(isnan(DTM(i,j)))
            for k=1:length(Z) %% search directions
                [xp,yp] = pol2cart(((90-Z(k))*pi/180)*ones(1,length(SD)),SD);
                Ep=interp2(x,y,DTM,x(j)+xp,y(i)+yp,'nearest');
                dE=(Ep-Ep(1));
                mang= 90-atan(dE./SD)*180/pi;
                HZ(i,j,k)=nanmin(mang);
                if sum(isnan(HZ(i,j,:)))>0
                    disp('Error in the Horizon Angle !!!')
                end
            end
        else
            HZ(i,j,:)=NaN*ones(1,length(Z));
        end
    end 
end 
%close(bau); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return 


