%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Compute Soil Heat Flux     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[G0,T,O,Oice]=Soil_Heat_Profile_New(Tup,dt,Ttm1,ms,dz,Zs,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
    Phy,s_SVG,bVG,Osat,Ohy,Oicetm1,Otm1,G0,OPT_FR_SOIL,OPT_ST)
%%%INPUTS
%%%%%% Boundary Mixed Condition
%Tdown = NaN; %% [K]
%G0 = NaN; %% [W/m^2]
%%%%
% OPZ case 1
% %%%% G0 Gn
% %%%%%%%%%%%%%
% OPZ case 2
% %%% Tup Tdown
% %%%%%%%%%%%%%%
% OPZ  case 3
% %%% Tup Gn
% %%%%%% Mixed Condition
if isnan(G0)
    OPZ=3;
end
if isnan(Tup)
    OPZ=1;
end
%%%
%%%%
Gn = 0; %% [W/m^2]
[lanS,~,~,Oice,O]=Soil_Thermal_properties_FT(Ttm1,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
    Phy,s_SVG,bVG,Osat,Ohy,(Oicetm1+Otm1),OPT_FR_SOIL);

roi = 916.2; %% ice density [kg/m^3]
Lf= 1000*333.7; %% [J/Kg] %% Latent heat melt/freezing
%%% Additional energy from cryosuction previous step 
dE = (roi*Lf*(Oice-Oicetm1))*0.001.*dz/dt; %% [W/m2]

%%%%%%%
Tdown= Ttm1(ms) - Gn*(0.001*dz(ms)*0.5)/lanS(ms); %%[K]
if OPZ == 1
    Tup = Ttm1(1) + G0*(0.001*dz(1)*0.5)/lanS(1);%% [K]
end
%%%%%%%%%%
OPZ_SOLV = 1; 

if OPZ_SOLV == 1
    Ttm1= [ Tup Ttm1 Tdown];
    %%%%%%%
    nit=10; %%% Number of Internal time step
    %%%%
    %G0p = zeros(1,nit); 
    for kk=1:nit
        [Tout] = Heat_CN_FT_Function(Ttm1,dt/nit,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
            Phy,s_SVG,bVG,Osat,Ohy,Oicetm1,Otm1,Zs,G0,Gn,Tup,Tdown,dE,OPZ,OPT_FR_SOIL);
        Ttm1 = Tout;
        
        %%%%     %%%%%%%%%%%%%%%
        %if  (OPZ == 2) || (OPZ == 3)
        %    G0p(kk)=-lanS(1).*((Ttm1(2)-Tup)./(0.001*dz(1)*0.5));  % [W/m^2]
        %end
        
    end
    T=Tout(2:end-1);% Temperature Layer %% [°C ]
    T=T';   
    
else
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    T_SPAN = [0 dt];
    %OPT_ST=odeset('AbsTol',0.01,'MaxStep',dt);
    %%%%%%%
    Dz=0.5*(dz(1:ms-1)+dz(2:ms));
    
    %%%%% Differential Equation Solving the Heat Transfer
    [tout,Tout]=ode23s(@SOIL_THERMAL_COND_FT,T_SPAN,Ttm1,OPT_ST,...
        Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
        Phy,s_SVG,bVG,Osat,Ohy,Oicetm1,Otm1,dz,ms,Dz,G0,Gn,Tup,Tdown,dE,OPZ,OPT_FR_SOIL);
    
    T = Tout(end,:);
end

%%%%%%%%%%%%%%%%%%%%%%
if isnan(sum(T))
    %%%%%%%%%%%%%%%
    disp('NaN values in Soil Temperature')
    return
end

[lanS,~,~,Oice,O]=Soil_Thermal_properties_FT(T,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
    Phy,s_SVG,bVG,Osat,Ohy,(Otm1+Oicetm1),OPT_FR_SOIL);
%%%%%%%%%%%%%%%%%%%%%%
%V = (O-Ohy).*dz; %%% [mm]
%Vice = (Oice).*dz;%% [mm]
%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%
if  (OPZ == 2) || (OPZ == 3) 
    G0=-lanS(1).*((T(1)-Tup)./(0.001*dz(1)*0.5));  % [W/m^2]
    %G0=mean(G0p);
end
%Tdown= T(nz) - Gn*(0.001*dz(nz)*0.5)/lanS(nz); %%[°C][K]
return