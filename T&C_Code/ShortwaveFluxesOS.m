%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction ShortwaveFluxesOtherSurface   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Rabsb,Rrfl_vis,Rrfl_nir]=ShortwaveFluxesOS(Rsw,SvF,alb)
%%%%% INPUT
%%% Variables:  
% Rsw  shortwave [W/m^2]
%Rsw.dir_vis, Rsw.dir_nir,
%Rsw.dif_vis, Rsw.dif_nir
%alb [0-1]  Albedo 
%alb.dir_vis, alb.dir_nir,
%alb.dif_vis, alb.dif_nir
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% ShF shadow factor
% SvF sky view factor
%%%%%%% OUTPUT 
%%% Rabs =[W/m^2] Absorbed Radiation 
Rabs_vis = Rsw.dir_vis*(1-SvF*alb.dir_vis)+ Rsw.dif_vis*(1-SvF*alb.dif_vis); %%% [W/m^2] 
Rabs_nir = Rsw.dir_nir*(1-SvF*alb.dir_nir)+ Rsw.dif_nir*(1-SvF*alb.dif_nir); %%% [W/m^2] 
Rabsb = Rabs_vis + Rabs_nir; %% [W/m^2] 
Rrfl_vis =  Rsw.dir_vis*(SvF*alb.dir_vis)+ Rsw.dif_vis*(SvF*alb.dif_vis); 
Rrfl_nir =  Rsw.dir_nir*(SvF*alb.dir_nir)+ Rsw.dif_nir*(SvF*alb.dif_nir);
%%%%%%%%%%%%%%%%%%%%%%%%%%
return 