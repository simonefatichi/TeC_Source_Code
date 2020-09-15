%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction SetSunVariables   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h_S,delta_S,zeta_S,T_sunrise,T_sunset,L_day] = SetSunVariables(Datam,DeltaGMT,Lon,Lat,t_bef,t_aft)
%%% INPUT 
%%% Datam %% [Yr, MO, DA, HR]
%%% DeltaGMT [°] 
%%% Lon [°]
%%% Lat [°]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUT 
%%delta_S,  Solar declination
%%tau_S,  Hour angle of the sun
%%h_S, [rad] solar altitude
%%zeta_S, Sun's azimuth
%%T_sunrise, [h]  sunrise time, 
%%T_sunset,  [h]sunset time, 
%%L_day, [h] total day length
%%r_ES,[] ratio of the actual Earth-Sun to the mean Earth-Sun  distance
%%jDay,    Julian Day 
%%Delta_TSL [h] Time difference between standard and local meridian
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Determine the julian day of the current time
days = [31 28 31 30 31 30 31 31 30 31 30 31];    
nowYR=Datam(1); nowMO=Datam(2); nowDA=Datam(3); nowHR =Datam(4);    
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
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compute solar declination
delta_S = 23.45*pi/180*cos(2*pi/365*(172 - jDay));
%%% Compute time difference between standard and local meridian
if(Lon<0)
    Delta_TSL = -1/15*(15*abs(DeltaGMT) - abs(Lon));
else
    Delta_TSL = 1/15*(15*abs(DeltaGMT) - abs(Lon));
end
%t_bef= 0.5; % 0.5; %% Integration intervals  1
%t_aft= 0.5; % 0.5; %% Integration intervals  0
t= nowHR-t_bef:0.0166666:nowHR+t_aft;
tau_S=zeros(1,length(t)); 
for i=1:length(t)
    %%%  Compute hour angle of the sun
    if(t(i) < (12 + Delta_TSL))
        tau_S(i) = 15*pi/180*(t(i) + 12 - Delta_TSL);
    else
        tau_S(i) = 15*pi/180*(t(i) - 12 - Delta_TSL);
    end
end
%%%% Compute solar altitude
Lat_rad = Lat*pi/180;
sinh_S = sin(Lat_rad)*sin(delta_S) + cos(Lat_rad)*cos(delta_S)*cos(tau_S);
h_S = asin(sinh_S);
h_S=mean(h_S); 
%%%% Compute Sun's azimuth
zeta_S = atan(-sin(tau_S)./(tan(delta_S)*cos(Lat_rad) - sin(Lat_rad)*cos(tau_S)));
%%%%%%%%%%%%%%%%%%%
for i=1:length(t)
    if (tau_S(i) >0 && tau_S(i) <= pi)
        if (zeta_S(i) > 0.)
            zeta_S(i) =  zeta_S(i) + pi;
        else
            zeta_S(i) = zeta_S(i) + (2.*pi);
        end
    elseif (tau_S(i) >=pi && tau_S(i) <= 2*pi)
        if (zeta_S(i) < 0.)
            zeta_S(i) =  zeta_S(i) + pi;
        end
    end
end
zeta_S=mean(zeta_S);
%%%% Compute sunrise time, sunset time, and total day length
T_sunrise = 180/(15*pi)*(2*pi - acos(-tan(delta_S)*tan(Lat_rad))) - 12;
T_sunset  = 180/(15*pi)*acos(-tan(delta_S)*tan(Lat_rad)) + 12;
L_day     = 360/(15*pi)*acos(-tan(delta_S)*tan(Lat_rad));
T_sunrise=real(T_sunrise); T_sunset = real(T_sunset); L_day=real(L_day); 
end