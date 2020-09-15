function [HZ,Z] = Horizon_Angle(DTM,cellsize)
%%% INPUT  
%%% DTM : matrix digital elevation model  
%%% cellsize dimension cell m 
%%%% OUTPUT 
%%% HZ Horizon angle array [angular degree]   
%%% Z Azimuth directions  [angualr degree] from N  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(DTM);  
%%%% COMPUTATION HORIZON ANGLE f(azimuth) 
%bau = waitbar(0,'Internal Loop');
HZ=zeros(m,n,8); 
Z= [0 45 90 135 180 225 270 315];
for i=1:m
    for j=1:n
        %waitbar(i/m,bau);
        if not(isnan(DTM(i,j)))
            %%% Direction North
            Ep = DTM(i:m,j);
            dE=(Ep-Ep(1));
            dL=cellsize:cellsize:cellsize*length(Ep);
            mang= 90-atan(dE./dL')*180/pi;
            H_z(1)=nanmin(mang);
            %%% Direction North -East
            r=min(length(i:m),length(j:n));
            for k=0:r-1
                Ep(k+1) = DTM(i+k,j+k);
                %x2(k+1)=i+k;
                %y2(k+1)=j+k;
            end
            dE=(Ep-Ep(1));
            dL=cellsize:sqrt(2)*cellsize:sqrt(2)*cellsize*length(Ep);
            mang= 90-atan(dE./dL')*180/pi;
            H_z(2)=nanmin(mang);
            %%%% Direction East 
            Ep = DTM(i,j:n);
            dE=(Ep-Ep(1));
            dL=cellsize:cellsize:cellsize*length(Ep);
            mang= 90-atan(dE./dL)*180/pi;
            H_z(3)=nanmin(mang);
            %%%%%%% Direction South - East 
            r=min(length(1:i),length(j:n));
            for k=0:r-1
                Ep(k+1) = DTM(i-k,j+k);
            end
            dE=(Ep-Ep(1));
            dL=cellsize:sqrt(2)*cellsize:sqrt(2)*cellsize*length(Ep);
            mang= 90-atan(dE./dL)*180/pi;
            H_z(4)=nanmin(mang);
            %%%%%% Direction South 
            Ep = DTM([i:-1:1],j);
            dE=(Ep-Ep(1));
            dL=cellsize:cellsize:cellsize*length(Ep);
            mang= 90-atan(dE./dL')*180/pi;
            H_z(5)=nanmin(mang);
            %%%%%%%%% Direction South - West 
            r=min(length(1:i),length(1:j));
            for k=0:r-1
                Ep(k+1) = DTM(i-k,j-k);
            end
            dE=(Ep-Ep(1));
            dL=cellsize:sqrt(2)*cellsize:sqrt(2)*cellsize*length(Ep);
            mang= 90-atan(dE./dL')*180/pi;
            H_z(6)=nanmin(mang);
            %%%%%%%%% Direction West 
            Ep = DTM(i,[j:-1:1]);
            dE=(Ep-Ep(1));
            dL=cellsize:cellsize:cellsize*length(Ep);
            mang= 90-atan(dE./dL)*180/pi;
            H_z(7)=nanmin(mang);
            %%%%%%% Direction North-West 
            r=min(length(i:m),length(1:j));
            for k=0:r-1
                Ep(k+1) = DTM(i+k,j-k);
            end
            dE=(Ep-Ep(1));
            dL=cellsize:sqrt(2)*cellsize:sqrt(2)*cellsize*length(Ep);
            mang= 90-atan(dE./dL)*180/pi;
            H_z(8)=nanmin(mang);
        else
           H_z=NaN*ones(1,8);  
        end
        HZ(i,j,:)=H_z; 
        clear H_z 
    end
end
%close(bau); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
return 

