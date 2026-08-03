function [SAB1,SAB2,SAD1,SAD2,PARB,PARD] = Radiation_Partition(Rsw,Datam,DeltaGMT,Lon,Lat,Zbas,COEFF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Radiation_Partition                                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Partition measured incoming shortwave into the six radiation variables
%%% T&C requires, for sites where only total Rsw is observed (FLUXNET, ICOS).
%%%
%%% INPUT
%%% Rsw       [W/m^2]  total incoming shortwave, hourly
%%% Datam     [Yr MO DA HR]
%%% DeltaGMT  [-]      offset of the local meridian from Greenwich
%%% Lon, Lat  [deg]
%%% Zbas      [m]      elevation, for the pressure correction to airmass
%%% COEFF     struct: f_diff, fvis_dir, fvis_dif, par_dir, par_dif, hour_offset
%%%
%%% OUTPUT  SAB1 SAD1 direct/diffuse visible, SAB2 SAD2 direct/diffuse NIR,
%%%         PARB PARD direct/diffuse PAR. Closure SAB1+SAB2+SAD1+SAD2 = Rsw
%%%         holds to machine precision by construction.
%%%
%%% WHY THREE FRACTIONS AND NOT ONE
%%% The visible share is a property of the STREAM, not of the site. At Zurich
%%% Fluntern over 1981-2012 the measured bands give
%%%       SAB1/(SAB1+SAB2) = 0.389 +/- 0.151     (direct)
%%%       SAD1/(SAD1+SAD2) = 0.537 +/- 0.101     (diffuse)
%%% a 15 point gap, which is Rayleigh scattering: the blue removed from the beam
%%% is what reappears as diffuse. Applying one visible fraction to both streams
%%% cannot reproduce the four bands however good the diffuse fraction is --
%%% imposing the true fractions per stream reconstructs them to RMSE 1.4 W/m^2,
%%% imposing one shared true fraction leaves 29.6 W/m^2.
%%%
%%% PAR is not a fixed multiple of the visible band either: PARB/SAB1 ranges
%%% 0.917 to 0.982 and PARD/SAD1 ranges 0.772 to 0.958, so both are predicted.
%%%
%%% HOUR CONVENTION
%%% COEFF.hour_offset shifts the timestamp before solar geometry is computed.
%%% Getting this wrong is expensive and silent: the shipped Zurich forcing is
%%% stamped in UTC while DeltaGMT = 1, and treating the stamp as local time
%%% leaves 6722 hours with Rsw > 20 W/m^2 while the computed sun is below the
%%% horizon, costing 19 per cent of aggregate RMSE before any partitioning
%%% model is chosen. calibrate_hour_offset.m determines it from the data.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 7 || isempty(COEFF)
    COEFF = Radiation_Partition_Coeff_Zurich();
end
Rsw = Rsw(:); NN = length(Rsw);
Yr = Datam(:,1); Mo = Datam(:,2); Da = Datam(:,3); Hr = Datam(:,4) + COEFF.hour_offset;

days = [31 28 31 30 31 30 31 31 30 31 30 31];
cum  = [0 cumsum(days)];
jDay = cum(Mo)' + Da;
leap = (mod(Yr,4)==0 & (mod(Yr,100)~=0 | mod(Yr,400)==0)) & Mo>2;
jDay = jDay + double(leap);
gam  = 2*pi*(jDay-1)/365;
delta = 0.006918 - 0.399912*cos(gam) + 0.070257*sin(gam) ...
      - 0.006758*cos(2*gam) + 0.000907*sin(2*gam) ...
      - 0.002697*cos(3*gam) + 0.00148*sin(3*gam);
EoT = 229.18*(0.000075 + 0.001868*cos(gam) - 0.032077*sin(gam) ...
    - 0.014615*cos(2*gam) - 0.040849*sin(2*gam));
TST  = Hr*60 + 4*(Lon - 15*DeltaGMT) + EoT;
HA   = (TST/4 - 180)*pi/180;
LatR = Lat*pi/180;
sinh_S = max(min(sin(LatR)*sin(delta) + cos(LatR)*cos(delta).*cos(HA),1),-1);
h_S  = asin(sinh_S);
r_ES = 1.00011 + 0.034221*cos(gam) + 0.00128*sin(gam) ...
     + 0.000719*cos(2*gam) + 0.000077*sin(2*gam);

cosz = max(sin(h_S),0);
I0   = 1367*r_ES.*cosz;
kt   = zeros(NN,1); ii = I0 > 1;
kt(ii) = Rsw(ii)./I0(ii);
kt   = max(min(kt,1.2),0);

zdeg = max(h_S,1e-6)*180/pi;
m    = 1./(sin(zdeg*pi/180) + 0.50572*(zdeg + 6.07995).^(-1.6364));
m    = max(min(m,40),1) * exp(-Zbas/8434);
lm   = log(max(min(m,40),1));

kprev = [kt(1); kt(1:end-1)];
knext = [kt(2:end); kt(end)];
pers  = 0.5*(kprev + knext);

one = ones(NN,1);
Xfd = [one kt kt.^2 lm cosz pers kt.*lm];
Xfv = [one lm lm.^2 kt kt.*lm cosz];
lg  = @(X,b) 1./(1 + exp(-max(min(X*b(:),40),-40)));

fd  = lg(Xfd, COEFF.f_diff);  fd(cosz <= 0) = 1;
vdr = lg(Xfv, COEFF.fvis_dir);
vdf = lg(Xfv, COEFF.fvis_dif);

Rdif = Rsw.*fd;  Rdir = Rsw.*(1-fd);
SAB1 = (Rdir.*vdr)';        SAB2 = (Rdir.*(1-vdr))';
SAD1 = (Rdif.*vdf)';        SAD2 = (Rdif.*(1-vdf))';
PARB = SAB1.*lg(Xfv, COEFF.par_dir)';
PARD = SAD1.*lg(Xfv, COEFF.par_dif)';
end
