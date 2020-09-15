function [SLO,ASP]=Slope_Aspect_indexes(DTM,cellsize,ANSW)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT -->
%%% DTM in [m]
%%% cellsize [m]
%%% ANSW ['grad'] method of gradient  or ['mste'] maximum steepness method 
%%% OUTPUT -->
%%% SLO [fraction]
%%% ASP [angular degree from North]
%%% SLOPE from GIS in Angular degree [°] --> method of gradient -- NOTE !!!
%%%%%%%%%%%%%%%%%%%%%%%%%%
[m,n]=size(DTM);
%%%%%%%%%%%%%%%%%%%%%%%%%%
if strcmp(ANSW,'mste')
    SLO = zeros(m,n);  SLO=SLO.*DTM;
    d = [ sqrt(2)*cellsize cellsize sqrt(2)*cellsize ;
        cellsize 0 cellsize ;
        sqrt(2)*cellsize cellsize sqrt(2)*cellsize ];
    DTMf= [ NaN*ones(1,n+2); [NaN*ones(m,1), DTM, NaN*ones(m,1)] ; NaN*ones(1,n+2) ];
    SLO= [ NaN*ones(1,n+2); [NaN*ones(m,1), SLO, NaN*ones(m,1)] ; NaN*ones(1,n+2) ];
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    for i=2:(m+1)
        for j=2:(n+1)
            if not(isnan(DTMf(i,j)))
                Z= (DTMf(i-1:i+1,j-1:j+1)-DTMf(i,j))./d;
                %s=nanmin(nanmin(Z)); %%% maximum slope
                s=min(min(Z(not(isnan(Z))))); %% 
                SLO(i,j)=-s;
            end
        end
    end
    SLO=SLO(2:end-1,2:end-1);
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Calculate slope and aspect (deg) using GRADIENT function
[fx,fy] = gradient(DTM,cellsize,cellsize); % uses simple, unweighted gradient of immediate neighbours
[asp,grad]=cart2pol(fy,fx); % convert to carthesian coordinates
grad=atan(grad); % [rad] steepest slope
%asp=asp.*-1+pi; % convert asp 0 facing south
asp = asp +pi ; % [rad] convert asp 0 facing north
%%%%%%
for i=1:m
    for j=1:n
        if not(isnan(DTM(i,j))) && (isnan(asp(i,j)) || isnan(grad(i,j)))
            grad(i,j)=0;
            asp(i,j)= 0;
        end
    end
end
%%%%%%%%%%%%%%%
if strcmp(ANSW,'grad')
    %%SLO = 180*grad/pi; %%% Slope in angular degree [°]
    SLO = tan(grad);%% Slope in fraction []
end
ASP =  180*asp/pi; %%% Aspect in angular degree [°]
%%%%%%%%%%%
end