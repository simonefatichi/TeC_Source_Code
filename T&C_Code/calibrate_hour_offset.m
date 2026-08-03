function [best, tbl] = calibrate_hour_offset(Rsw,Datam,DeltaGMT,Lon,Lat)
%%% Determine the timestamp convention from the data itself.
%%% Scans candidate offsets and returns the one minimising physically
%%% impossible hours: measurable shortwave with the sun below the horizon, or a
%%% high sun with zero shortwave. At Zurich this returns +1.0 h, reducing
%%% inconsistent hours from 6722 to 22 out of 276096.
cands = -2:0.5:2; tbl = zeros(numel(cands),3);
for k = 1:numel(cands)
    C.hour_offset = cands(k);
    C.f_diff=[0;0;0;0;0;0;0]; C.fvis_dir=zeros(6,1); C.fvis_dif=zeros(6,1);
    C.par_dir=zeros(6,1); C.par_dif=zeros(6,1);
    [~,~,~,~,~,~] = deal(0,0,0,0,0,0);
    Hr = Datam(:,4) + cands(k);
    cosz = local_cosz(Datam(:,1),Datam(:,2),Datam(:,3),Hr,DeltaGMT,Lon,Lat);
    bad1 = sum(Rsw(:) > 20 & cosz <= 0);
    bad2 = sum(cosz > 0.3 & Rsw(:) <= 0);
    tbl(k,:) = [cands(k) bad1 bad2];
end
[~,i] = min(tbl(:,2) + tbl(:,3));
best = tbl(i,1);
end

function cosz = local_cosz(Yr,Mo,Da,Hr,DeltaGMT,Lon,Lat)
days = [31 28 31 30 31 30 31 31 30 31 30 31]; cum = [0 cumsum(days)];
jDay = cum(Mo)' + Da;
leap = (mod(Yr,4)==0 & (mod(Yr,100)~=0 | mod(Yr,400)==0)) & Mo>2;
jDay = jDay + double(leap);
gam = 2*pi*(jDay-1)/365;
delta = 0.006918 - 0.399912*cos(gam) + 0.070257*sin(gam) - 0.006758*cos(2*gam) ...
      + 0.000907*sin(2*gam) - 0.002697*cos(3*gam) + 0.00148*sin(3*gam);
EoT = 229.18*(0.000075 + 0.001868*cos(gam) - 0.032077*sin(gam) ...
    - 0.014615*cos(2*gam) - 0.040849*sin(2*gam));
HA = ((Hr*60 + 4*(Lon - 15*DeltaGMT) + EoT)/4 - 180)*pi/180;
LatR = Lat*pi/180;
cosz = max(sin(LatR)*sin(delta) + cos(LatR)*cos(delta).*cos(HA),0);
end
