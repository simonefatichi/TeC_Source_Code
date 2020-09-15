%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Debris Thermal Mod %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [T] = Heat_CN_Normal_Function(Ttm1,dts,lanS,cv,Zs,Tup,Tdown,Gn,G0,OPZ)
Ttm1=flip(Ttm1); lanS=flip(lanS);  cv=flip(cv);
%%%%%%%%%%%%%%%%%%
Ttm1=reshape(Ttm1,length(Ttm1),1);
lanS=reshape(lanS,length(lanS),1);
cv=reshape(cv,length(cv),1);

Zs=sort(-Zs');
nz=length(Zs)-1; %% number of layers
nn = nz + 2; %%% number of nodes for temperature
dz= diff(Zs); %%%% [mm]  Thickness of the Layers
dz=dz*0.001; %% [m]

%Dz=[dz(1)*0.5, 0.5*(dz(1:nz-1)+dz(2:nz))];
DzC = [dz(1)*0.5 ; 0.5*(dz(1:nz-1)+dz(2:nz)) ; dz(end)*0.5];

RF = 0; %% W/m2 Absorbed radiation or additional energy in each layer of soil


%%% Harmonic Mean of lanS
%lanS_half(2:nz)= (dz(1:nz-1) + dz(2:nz))./(dz(1:nz-1)./(lanS(1:nz-1)) + dz(2:nz)./lanS(2:nz));
lanS_half(2:nz)= 0.5*(lanS(1:nz-1) + lanS(2:nz));
lanS_half(1)=lanS(1); lanS_half(nz+1)=lanS(nz);
lanS_half=reshape(lanS_half,length(lanS_half),1);

% --- Coefficients of the tridiagonal system
alp = 0.5; %% Tridiagonal system of Equation

%%% i ...  2:nn-1
%a = % subdiagonal a: coefficients of T(i-1)  %[1/s]
%c = a; % superdiagonal c: coefficients of T(i+1) %[1/s]
a = -(1-alp)*lanS_half(1:nz)./(cv(1:nz).*dz(1:nz).*DzC(1:nz));
c = -(1-alp)*lanS_half(2:nz+1)./(cv(1:nz).*dz(1:nz).*DzC(2:nz+1));
b = (1/dts)*ones(nz,1) - (a+c); % diagonal b: coefficients of T(i) %[1/s]
% Right hand side includes time derivative and Crank-Nicolson terms
a=[ NaN ; a ; NaN] ; c=[ NaN ; c ; NaN] ;  b=[ NaN ; b ; NaN] ;
d = Ttm1/dts - [NaN; a(2:nn-1).*Ttm1(1:nn-2); NaN] ...
    + [NaN; (a(2:nn-1)+c(2:nn-1)).*Ttm1(2:nn-1); NaN] ...
    - [NaN; c(2:nn-1).*Ttm1(3:nn); NaN] + [NaN ; RF./(cv.*dz) ; NaN]; %%[°C /s]
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
    case 4
        %%%%%%%%% Mixed Condition
        %%%% Set upper boundary condition.
        b(nn,1)= 1;
        a(nn,1)= -1;
        d(nn,1)= (G0)*DzC(nz+1)/lanS_half(nz+1);
        %%% Set lower boundary condition.
        b(1,1)=1;
        c(1,1)=0;
        d(1,1)= Tdown;
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