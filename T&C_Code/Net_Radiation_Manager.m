%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Net_Radiation_Manager      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Rnet]=Net_Radiation_Manager(Ts,Latm,SvF,...
Csno,Crock,Curb,Cwat,Cbare,Cice,Ccrown,...
    hc_L,SnoDep,ydepth,ICE_D,Cdeb,LAI_H,LAI_L,SAI_H,SAI_L,...
    RabsbSun_vegH,RabsbShd_vegH,Rabsb_soiH,...
    RabsbSun_vegL,RabsbShd_vegL,Rabsb_soiL,FsunH,FshdH,...
    FsunL,FshdL,Rabsb_sno,Rabsb_bare,Rabsb_urb,Rabsb_wat,Rabsb_rock,Rabsb_ice,Rabsb_deb,...
    e_sno,e_gr,e_sur,Cicew,Csnow) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT
%%%%%---> Csno - Crock - Curb - Cwat - Cbare - Ccrown (1...n)  <---%%%%%%%%
% SvF [0-1] sky view factor
% Ta air temperature [°C]
% ea vapor pressure [Pa]
% N Cloudiness [0-1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Ts surface radiative temperature [°C]
%%% LAI_H SAI_H LAI_L SAI_L
%%% Csno
%%% hc_L
%%% D
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OUTPUT
%%% Rnet
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% HIP. ONLY ONE TEMPERATURE %%%%%%%%%%%%%%%%%%%%%%%%%
Twat=Ts; Tice = Ts;  Trock=Ts;  Turb=Ts;
Tg=Ts;    TvHsun=Ts;  TvHshd=Ts;
TvLsun=Ts; TvLshd=Ts; Tsno=Ts; Tdeb=Ts;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Incoming Longwave  %%%%%%%%%%%
%[Latm]=Incoming_Longwave(Ta,ea,N); % Latm Incoming LongWave Radiation [W/m^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Absorbed Radiation Not Vegetated Patches % -Other Surfaces - %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Cwat > 0
    [Labsb_wat]=LongwaveFluxesOS(Twat,Latm,SvF,e_sur.W); %% [W/m^2]
else
    Labsb_wat=0;
end
if Crock > 0
    [Labsb_rock]=LongwaveFluxesOS(Trock,Latm,SvF,e_sur.R); %% [W/m^2]
else
    Labsb_rock=0;
end
if Curb > 0
    [Labsb_urb]=LongwaveFluxesOS(Turb,Latm,SvF,e_sur.U); %% [W/m^2]
else
    Labsb_urb=0;
end
%%%%%%%%%%%%%%%%%%%%
if Cbare > 0
    [Labsb_bare]=LongwaveFluxesOS(Tg,Latm,SvF,e_gr); %% [W/m^2]
else
    Labsb_bare=0;
end
if Csno > 0
    [Labsb_sno]=LongwaveFluxesOS(Tsno,Latm,SvF,e_sno); %% [W/m^2]
else
    Labsb_sno=0;
end
if Cice > 0
    if Cdeb == 0
        [Labsb_ice]=LongwaveFluxesOS(Tice,Latm,SvF,e_sur.I); %% [W/m^2]
        Labsb_deb=0; 
    else
        [Labsb_deb]=LongwaveFluxesOS(Tdeb,Latm,SvF,e_sur.D); %% [W/m^2]
        Labsb_ice=0; 
    end
else
    Labsb_ice=0;   Labsb_deb=0;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc=length(Ccrown);
Labsb_vegH=zeros(1,cc);
Labsb_vegL=zeros(1,cc);
Labsb_soi=zeros(1,cc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:cc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LONGWAVE %%%%%%%%%%%%%%%%%%%
    [Labsb_vegH(i),Labsb_vegL(i),Labsb_soi(i)]=LongwaveFluxes(Tg,TvHsun,TvHshd,TvLsun,TvLshd,Tsno,...
        Latm,SvF,e_gr,LAI_H(i),SAI_H(i),LAI_L(i),SAI_L(i),FshdH(i),FsunH(i),FshdL(i),FsunL(i),hc_L(i),SnoDep,e_sno,Csno);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ydepth > 0 && Csno == 0
        [Labsb_vegH(i),Labsb_vegL(i),Labsb_soi(i)]=LongwaveFluxes(Tg,TvHsun,TvHshd,TvLsun,TvLshd,Twat,...
            Latm,SvF,e_gr,LAI_H(i),SAI_H(i),LAI_L(i),SAI_L(i),FshdH(i),FsunH(i),FshdL(i),FsunL(i),hc_L(i),ydepth,e_sur.W,1);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
    if Cice == 1 && Csno == 0
        [Labsb_vegH(i),Labsb_vegL(i),Labsb_soi(i)]=LongwaveFluxes(Tg,TvHsun,TvHshd,TvLsun,TvLshd,Tice,...
            Latm,SvF,e_gr,LAI_H(i),SAI_H(i),LAI_L(i),SAI_L(i),FshdH(i),FsunH(i),FshdL(i),FsunL(i),hc_L(i),ICE_D,e_sur.I,Cice);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    end
end
%%%%%%%%%%%%%%%%%%%%%% CELL -QUANTITY
Rnet_urb = Curb*(Rabsb_urb + Labsb_urb)*(1-Csno)*(1-Cice);
Rnet_wat = Cwat*(Rabsb_wat + Labsb_wat)*(1-Csnow)*(1-Cicew)*(1-Cdeb);
Rnet_rock = Crock*(Rabsb_rock +Labsb_rock)*(1-Csno)*(1-Cice);
Rnet_ice =  Cice*(Rabsb_ice + Labsb_ice)*(1-Csno)*(1-sum(Ccrown)-Cwat) + Cicew*(Rabsb_ice + Labsb_ice)*Cwat*(1-Csnow);
Rnet_snow = Csno*(Rabsb_sno + Labsb_sno)*(1-sum(Ccrown)-Cwat) +  Csnow*(Rabsb_sno + Labsb_sno)*Cwat;
Rnet_deb = Cdeb*(Rabsb_deb +Labsb_deb)*(1-Csno);
%%%%%%%%%%%%
Rnet_ground =  Cbare*(Rabsb_bare + Labsb_bare)*(1-Csno)*(1-Cice) + sum(Ccrown.*(Rabsb_soiH + Rabsb_soiL + Labsb_soi));
RnetSun_vegH = (1-Cice)*(1-Csno)*sum(Ccrown.*(RabsbSun_vegH + Labsb_vegH.*FsunH));
RnetShd_vegH = (1-Cice)*(1-Csno)*sum(Ccrown.*(RabsbShd_vegH + Labsb_vegH.*FshdH));
RnetSun_vegL = (1-Cice)*(1-Csno)*sum(Ccrown.*(RabsbSun_vegL + Labsb_vegL.*FsunL));
RnetShd_vegL = (1-Cice)*(1-Csno)*sum(Ccrown.*(RabsbShd_vegL + Labsb_vegL.*FshdL));
%%%%%%%%%%%%%%%%%%%%%%%
Rnet_vegH = RnetSun_vegH + RnetShd_vegH;
Rnet_vegL = RnetSun_vegL + RnetShd_vegL;
%%% NET RADIATION [W/m^2]
Rnet = Rnet_ground + Rnet_vegH + Rnet_vegL + Rnet_snow + Rnet_rock +Rnet_wat + Rnet_urb + Rnet_ice + Rnet_deb;
%%%%%%%%%%%%%%%%%%%%%%%
return