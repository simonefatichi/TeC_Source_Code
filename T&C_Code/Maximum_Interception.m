%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Maximum_Interception             %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%
function[In_max_SWE,In_max_L,In_max_H]=Maximum_Interception(Ccrown,LAI_L,LAI_H,SAI_H,SAI_L,...
    Ta,Sp_SN_In,Sp_LAI_L_In,Sp_LAI_H_In) 
%%%INPUTS
% LAI_H  Leaf Area Index First Layer [] 
% LAI_L Leaf Area Index Second Layer []
% SAI_H Stem Area Index First Layer []
% SAI_L Stem Area Index Second Layer []
% Sp_SN_In = 5.9; %%% [mm/LAI] -- [kg/m^2 LAI] [Specific interception for snow ]
% Sp_LAI_L_In = 0.2; %%%[mm/LAI] [Specific interception first layer]
% Sp_LAI_H_In = 0.2; %%% [mm/LAI] [Specific interception second layer]
%%% OUTPUTS 
%In_max_H  maximum interception in First Layer [mm]
%In_max_L  maximum interception in Second Layer [mm]
%In_max_SWE  = Intercepted Snow maximum [mm] 
%%% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%ros_n = 67.92 + 51.25*exp(Ta/2.59); %%  [Hedstrom and Pomeray 1998]
ros_n = 1000*(0.05 + (((9/5)*Ta +32)> 0).*(((9/5)*Ta+32)./100).^2); %% [kg/m^3] new snow density [from Bras 1990]
%%% Plot scale 
In_max_SWE = Sp_SN_In.*sum((LAI_H+SAI_H).*Ccrown)*(0.27 + 46/ros_n); %%% [Hedstrom and Pomeroy 1998] ;
%%%% For unit of Ground 
In_max_L = Sp_LAI_L_In.*(LAI_L + SAI_L).*Ccrown ; %%% [Mahfouf and Jacquemin 1989] --- [Ivanov et al., 2008 ]
In_max_H = Sp_LAI_H_In.*(LAI_H + SAI_H).*Ccrown ; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end