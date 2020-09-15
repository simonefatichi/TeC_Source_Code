%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Net_Radiation_Manager      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Rnet]=Net_Radiation_Manager_VegSnow(Ts,Ts2,Latm,SvF,...
Csno,Ccrown,...
    hc_L,SnoDep,LAI_H,LAI_L,SAI_H,SAI_L,...
    RabsbSun_vegH,RabsbShd_vegH,FsunH,FshdH,...
    FsunL,FshdL,...
    e_sno,e_gr) 
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
Tg=Ts2;    TvHsun=Ts;  TvHshd=Ts;
TvLsun=Ts; TvLshd=Ts; Tsno=Ts2;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Incoming Longwave  %%%%%%%%%%%
%[Latm]=Incoming_Longwave(Ta,ea,N); % Latm Incoming LongWave Radiation [W/m^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%% Absorbed Radiation Not Vegetated Patches % -Other Surfaces - %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cc=length(Ccrown);
Labsb_vegH=zeros(1,cc);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for i=1:cc
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% LONGWAVE %%%%%%%%%%%%%%%%%%%
    [Labsb_vegH(i)]=LongwaveFluxes(Tg,TvHsun,TvHshd,TvLsun,TvLshd,Tsno,...
        Latm,SvF,e_gr,LAI_H(i),SAI_H(i),LAI_L(i),SAI_L(i),FshdH(i),FsunH(i),FshdL(i),FsunL(i),hc_L(i),SnoDep,e_sno,Csno);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%%%%%%%%%%%%%%%%%%%%%% CELL -QUANTITY
%%%%%%%%%%%%
RnetSun_vegH = (Csno)*sum(Ccrown.*(RabsbSun_vegH + Labsb_vegH.*FsunH));
RnetShd_vegH = (Csno)*sum(Ccrown.*(RabsbShd_vegH + Labsb_vegH.*FshdH));
%%%%%%%%%%%%%%%%%%%%%%%
Rnet_vegH = RnetSun_vegH + RnetShd_vegH;
%%% NET RADIATION [W/m^2]
Rnet =Rnet_vegH ;
%%%%%%%%%%%%%%%%%%%%%%%
return