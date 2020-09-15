%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Thermal Mod %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T] = Heat_CN_FT_Function(Ttm1,dts,Pre,rsd,lan_dry,lan_s,cv_s,SPAR,L,Pe,O33,alpVG,nVG,...
    Phy1,s_SVG,bVG,Osat,Ohy,Oicetm1,Otm1,Zs,G0,Gn,Tup,Tdown,AE,OPZ,OPT_FR_SOIL)
Ttm1=flip(Ttm1); 
Ttm1=reshape(Ttm1,length(Ttm1),1); 
%%%%%%%%%%%%%%%%%%%%%%%%%%
Zs=sort(-Zs');
nz=length(Zs)-1; %% number of layers 
nn = nz + 2; %%% number of nodes for temperature 
dz= diff(Zs); %%%% [mm]  Thickness of the Layers
dz=dz*0.001; %% [m]

%Dz=[dz(1)*0.5, 0.5*(dz(1:nz-1)+dz(2:nz))];
DzC = [dz(1)*0.5 ; 0.5*(dz(1:nz-1)+dz(2:nz)) ; dz(end)*0.5];
%%%%%%%%%%%%%%%%%%%%%%%%%%%

g=9.81; %% m/s2
roi = 916.2; %% ice density [kg/m^3]
Lf= 1000*333.7; %% [J/Kg] %% Latent heat melt/freezing
%%%%%%%%%%%
[lanS,cv_Soil,~,~,~,Cw]=Soil_Thermal_properties_FT(Ttm1(2:nn-1),Pre,flip(rsd),flip(lan_dry),flip(lan_s),flip(cv_s),SPAR,flip(L),flip(Pe),flip(O33),flip(alpVG),flip(nVG),...
    Phy1,flip(s_SVG),flip(bVG),flip(Osat),flip(Ohy),flip(Oicetm1+Otm1),OPT_FR_SOIL);
%%%% cv %  Volumetric heat capacity Soil  [J/m^3 K]
%capp = cv_Soil - roi*Lf*dOice_dT;%%% [J/m^3 K]
capp = cv_Soil' + roi*Lf^2*(1000*Cw')./(g*(Ttm1(2:nn-1)+273.15));%% [J/m^3 K]
%capp = cv_Soil' ;% %%% [J/m^3 K]
capp(isnan(capp))=cv_Soil(isnan(capp));
%%%%%%%%%%%%%%%
lanS=reshape(lanS,length(lanS),1); 
capp=reshape(capp,length(capp),1); 
%%%%%%%%%%%%%%%%%

%AE [W/m2] Additional Energy to the layer 
RF = 0; %% W/m2 Absorbed radiation or additional energy in each layer of soil
RF= RF+flip(AE'); 

%%% Harmonic Mean of lanS
%lanS_half(2:nz)= (dz(1:nz-1) + dz(2:nz))./(dz(1:nz-1)./(lanS(1:nz-1)) + dz(2:nz)./lanS(2:nz));
lanS_half(2:nz)= 0.5*(lanS(1:nz-1) + lanS(2:nz));  %%  [W/m K ] 
lanS_half(1)=lanS(1); lanS_half(nz+1)=lanS(nz);
lanS_half=reshape(lanS_half,length(lanS_half),1); 

% --- Coefficients of the tridiagonal system
alp = 0.5; %% Tridiagonal system of Equation

%%% i ...  2:nn-1
%a = % subdiagonal a: coefficients of T(i-1)  %[1/s]
%c = a; % superdiagonal c: coefficients of T(i+1) %[1/s]
a = -(1-alp)*lanS_half(1:nz)./(capp(1:nz).*dz(1:nz).*DzC(1:nz));
c = -(1-alp)*lanS_half(2:nz+1)./(capp(1:nz).*dz(1:nz).*DzC(2:nz+1));
b = (1/dts)*ones(nz,1) - (a+c); % diagonal b: coefficients of T(i) %[1/s]
% Right hand side includes time derivative and Crank-Nicolson terms
a=[ NaN ; a ; NaN] ; c=[ NaN ; c ; NaN] ;  b=[ NaN ; b ; NaN] ; 
d = Ttm1/dts - [NaN; a(2:nn-1).*Ttm1(1:nn-2); NaN] ...
    + [NaN; (a(2:nn-1)+c(2:nn-1)).*Ttm1(2:nn-1); NaN] ...
    - [NaN; c(2:nn-1).*Ttm1(3:nn); NaN] + [NaN ; RF./(capp.*dz) ; NaN]; %%[°C /s]
%%%%%
switch OPZ
    case 1
        %%%%%%%% Boundary Condition
        %%%%%% Neumann Condition
        %%%% Set upper boundary condition.
        b(nn,1)= 1;
        a(nn,1)= -1;
        d(nn,1)= (G0)*DzC(nz+1)/lanS_half(nz+1);
        %%% Set lower boundary condition.
        b(1,1)=-1;
        c(1,1)=1;
        d(1,1)= (Gn)*DzC(1)/lanS_half(1);
        
    case 2
        %%%%%%%% Boundary Condition
        %%%%%% Dirichlet Condition
        %%%% Set upper boundary condition.
        b(nn,1)= 1;
        a(nn,1)= 0;
        d(nn,1)= Tup;
        %%% Set lower boundary condition.
        b(1,1)=1;
        c(1,1)=0;
        d(1,1)= Tdown;
        
    case 3
        %%%%%%%%% Mixed Condition
        %%%% Set upper boundary condition.
        b(nn,1)= 1;
        a(nn,1)= 0;
        d(nn,1)= Tup;
        %%% Set lower boundary condition.
        b(1,1)=-1;
        c(1,1)=1;
        d(1,1)= (Gn)*DzC(1)/lanS_half(1);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% MATRIX COMPOSITION
% S = sparse(i,j,s,m,n,nzmax) uses vectors i, j, and s to
%generate an m-by-n sparse matrix such that S(i(k),j(k)) = s(k), with
%space allocated for nzmax nonzeros.
maindiag=sparse(1:nn,1:nn,b,nn,nn);
upper=sparse(1:(nn-1),2:nn,c(1:nn-1),nn,nn);
lower=sparse(2:nn,1:(nn-1),a(2:nn),nn,nn);
A=maindiag+upper+lower;
T=A\d; %%[°C]
T=flip(T); 
return