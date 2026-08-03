function [SIF_obs,mask] = SIF_Output(SIF_toc,Datam,DeltaGMT,Lon,Lat,revisit_d,overpass_h,window_h)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  SIF_Output                                               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Resample continuous modelled top-of-canopy SIF onto a satellite's
%%% observing schedule, so model output sits in the sensor's observation space
%%% and can be compared without further processing.
%%%
%%% Defaults are FLEX/FLORIS: 27-day repeat, ~10:00 local solar time descending
%%% node in tandem with Sentinel-3, averaged over a 1 h window.
%%%
%%% INPUT
%%% SIF_toc     [W m-2 sr-1 um-1]  hourly, from SIF_Escape.m
%%% Datam       [Yr MO DA HR]
%%% revisit_d   [d]   repeat cycle, default 27
%%% overpass_h  [h]   local solar time of overpass, default 10.0
%%% window_h    [h]   averaging window, default 1.0
%%%
%%% OUTPUT
%%% SIF_obs  same length as SIF_toc, NaN except on overpass steps
%%% mask     logical, true on overpass steps
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 6 || isempty(revisit_d);  revisit_d  = 27;   end
if nargin < 7 || isempty(overpass_h); overpass_h = 10.0; end
if nargin < 8 || isempty(window_h);   window_h   = 1.0;  end

Yr = Datam(:,1); Mo = Datam(:,2); Da = Datam(:,3); Hr = Datam(:,4);
dn = datenum(Yr,Mo,Da);
day0 = min(dn);
on_cycle = mod(dn - day0, revisit_d) == 0;

%%% local solar time, so the comparison is at the true overpass geometry
days = [31 28 31 30 31 30 31 31 30 31 30 31]; cum = [0 cumsum(days)];
jDay = cum(Mo)' + Da;
leap = (mod(Yr,4)==0 & (mod(Yr,100)~=0 | mod(Yr,400)==0)) & Mo>2;
jDay = jDay + double(leap);
gam = 2*pi*(jDay-1)/365;
EoT = 229.18*(0.000075 + 0.001868*cos(gam) - 0.032077*sin(gam) ...
    - 0.014615*cos(2*gam) - 0.040849*sin(2*gam));
LST = Hr + (4*(Lon - 15*DeltaGMT) + EoT)/60;

mask = on_cycle & abs(LST - overpass_h) <= window_h/2;
SIF_obs = nan(size(SIF_toc));
SIF_obs(mask) = SIF_toc(mask);
end
