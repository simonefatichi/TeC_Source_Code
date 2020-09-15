%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Aerodynamic Resistence     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%REFERENCES %%   - Mascart et al., 1995 --- Noilhan and Mafhouf 1996 ---
%%Kot and Song 1998 --- Bertoldi et al., 2004 ---
%%% Liu et al., 2007
%[Van den Hurk Holstag 1997]
%[Louis (1979)] [Mascart et al., 1995]  [Dyer 1974]
%%%% [ Businger et al., 1971] [Beljaars and Holstag (1991)]
%[Abdella and McFarlane 1996]   [Garrat 1992]
% [Su 2002]
function[ra]=Aerodynamic_Resistence(Ta,Ts,Pre,zatm,disp_h,zom,zoh,Ws,ea)
%%%INPUTS
%Ta = %% air temperature [°C] -- 
%Ts =  ; %% surface temperature [°C] -- 
Pre = Pre*100; %pressure [Pa]-- 
%z= %  reference height [m] % --- 
%zoh = %% roughness  eddy diffusivities for heat  [m]  
%zom = roughness eddy diffusivities for momentum [m]
u = Ws; %% [m/s] wind speed  ---
z= zatm - disp_h; 
%%% OUTPUTS 
%%% ra % [s/m]  Aerodynamic resistence  Heat flux 
%%% PARAMETERS 
g= 9.81; %% [m/s^2] gravity acceleration 
P_Ref= 100000; %% [Pa] reference pressure 
Rd =287.05; %% [J/kg K] dry air gas constant 
Rv = 1.607*Rd; %% [J/kg K]
cwv = 1858 + 0.382*(Ta) + 4.22e-04*(Ta)^2 - 1.996e-07*(Ta)^3 ; % [J/Kg K] specific heat water vapor  
cp=1005 + ((Ta +23.15)^2)/3364; %% specific heat air  [J/kg K]
k=0.4; %% Von Karman costant 
%%%%%
esatTs=611*exp(17.27*Ts./(237.3+Ts)); %% [Pa] Vapor pressure saturation
qTa=0.622*ea/(Pre-0.378*ea);% Specifc humidity Air  [] 
cair = cp*(1-qTa) + qTa.*cwv; cp=cair; 
Rair = Rd*(1-qTa) + qTa*Rv; Rd=Rair; 
%%%%% Tv = (1+0.61*q)*T; 
Tv = Ta*Pre./(Pre - 0.378*ea); %% Virtual  temperature air [K]
Tvs = Ts*Pre./(Pre - 0.378*1*esatTs); %% Virtual temperature surface [K] 
%%%% Tp = T* Pre/P_Ref^(Rd*(1-0.23*q))/cp)
Ova=(Tv +273.15)*(Pre/P_Ref)^-(Rd/cp); %% Virtual Potential temperature air [K]
Ovs=(Tvs +273.15)*(Pre/P_Ref)^-(Rd/cp); %% Virtual Potential temperature surface [K]
%%%%%%%%%%%%%%%%%%
%%% Richardson Number calculation
f2 = ((1-zom/z)^2)/(1-zoh/z); %% adj factor [Kot and Song (1998)]
Ri = -f2*((g*z*(Ovs-Ova))/(0.5*(Ovs+Ova)*u^2)); %% Richardson Number []
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ANS = 1;  %%% Approximate Solution CH = CDN*f(Ri) %%%  Exact Solution Obukohv Lenght lan  
%%%%%%%%%%%%%%%%%
if ANS == 1
    % - Mascart et al., 1995 --- Noilhan and Mafhouf 1996 ---
    %%Kot and Song 1998 --- Bertoldi et al., 2004 ---
    %%% Liu et al., 2007
    %%%%%%%%%% Aerodynamic resistence parameterization to Heat and Momentum exchange
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    mu=log(zom/zoh);
    ph=0.5802-0.1571*mu+0.0327*(mu^2)-0.0026*(mu^3);
    %pm=0.5233-0.0815*mu+0.0135*(mu^2)-0.0010*(mu^3);
    Chs=3.2165+4.3431*mu+0.5360*(mu^2)-0.0781*(mu^3);
    %Cms=6.8741+2.6933*mu-0.3601*(mu^2)+0.0154*(mu^3);
    %%%
    Cdn=(k^2)/((log(z/zom))^2); %% Transport coefficient --> Neutral condition
    Ch=15*Chs*Cdn*((z/zoh)^ph)*(log(z/zom)/log(z/zoh));
    %Cm=10*Cms*Cdn*((z/zoh)^pm);
    %%%%% Stability factor evaluation
    if Ri <= 0
        FH=(1-(15*Ri/(1+Ch*sqrt(abs(Ri)))))*(log(z/zom)/log(z/zoh));
     %   FM=(1-(10*Ri/(1+Cm*sqrt(abs(Ri))))); %%
    else
        FH=(1/(1+15*Ri*sqrt(1+5*Ri)))*(log(z/zom)/log(z/zoh)); %
      %  FM=(1/(1+(10*Ri)/(sqrt(1+5*Ri)))); %
    end
    %%%%%%%%
    CH = FH*Cdn; %% Transport coefficient -- Sensible Heat
    %CM = FM*Cdn; %%  Transport coefficient --- Momentum
    %%%%%%%%%
    ra = 1/(CH*u); %%% [s/m]  Aerodynamic resistence  Heat flux
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% Free Convection --  Beljaars 1994 
    if (u <= 0.05) && (Ovs>Ova)
        Cs =0.15; 
        ni= 1.51e-5; %% Kin. viscosity air [m2/s]
        Pr=0.71;
        CHu= Cs*(g*ni/(0.5*(Ovs+Ova)*Pr^2))^(1/3)*(Ovs-Ova)^(1/3);
        ra = 1/(CHu); %%% [s/m]  Aerodynamic resistence  free convection
    end
else
    if Ta == Ts 
        ra = (1/((k^2)*u))*(log(z/zom))*(log(z/zoh)); %%% Neutral aerodynamic resistance  %%[s/m]
    else
    %% First value  [m] Obukhov Length 
    if Ri > 0  %%%% Os < Oa  H < 0 LAN > 0 
        LANp= 2; %%% Stable condition 
    else %% Ri < 0 Os > Oa  H > 0 LAN < 0 
        LANp=  -2; %%% Unstable condition 
    end
    %%% LANp=((us^2)*(Ta+273.15))/(k*g*Oss); % [m] Obukhov Length 
    %%% LANp=-(ro*cp*us^3*(Ta+273.15))/(k*g*H); % [m] Obukhov Length 
    [LAN,us]=fzero(@Obukhov_length,LANp,[],z,zom,zoh,u,Ta,Ova,Ovs,k,g); %% 
    %%%%%%%%
    [Fih_z,Fim_z]=Businger_stability_functions(z/LAN);
    [Fih_zom,Fim_zom]=Businger_stability_functions(zom/LAN);
    [Fih_zoh,Fim_zoh]=Businger_stability_functions(zoh/LAN);
    %%%%%%%%%%%%%%%%%%%%%%%%%%  [Garrat 1992]
    ra = (1/(u*k^2))*((log(z/zom)-Fim_z +Fim_zom)*(log(z/zoh)-Fih_z + Fih_zoh)) ; %%% [s/m]  Aerodynamic resistence  Heat flux
    %CH = 1/(ra*u);
    end 
end    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% References 
        %[Abdella and McFarlane 1996] [Van den Hurk Holstag 1997]
        %[Dyer 1974] [Liu et al 2007] [Su 2002]
    function[LAN_d,us]=Obukhov_length(LAN,z,zom,zoh,u,Ta,Oa,Os,k,g)
        %%%%%%%%%%%%%%%%%%%%%%%%%
        Pr = 1; %%% neutral turbolent Prandtl Number 
        [Fih_z,Fim_z]=Businger_stability_functions(z/LAN);
        [Fih_zom,Fim_zom]=Businger_stability_functions(zom/LAN);
        [Fih_zoh,Fim_zoh]=Businger_stability_functions(zoh/LAN);
        %%%%%%%%%%%%%%%%%%%%%
        us= (k*u)/((log(z/zom)-Fim_z + Fim_zom));  % [m/s] Friction velocity
        Oss= (k*(Pr^-1)*(Oa-Os))/((log(z/zoh)-Fih_z + Fih_zoh)); %%  [K] Scaling Temperature
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        LAN_2=((us^2)*(Ta+273.15))/(k*g*Oss); % [m] Obukhov Length 
        %%%%%%%%%%%%%
        LAN_d = LAN-LAN_2; 
        %%%%%  Iteration 
    end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    function[Fih,Fim]=Businger_stability_functions(y)
        %%%%%%%% References 
        %[Van den Hurk Holstag 1997]
        %[Louis (1979)] [Mascart et al., 1995]  [Dyer 1974]
        %%%% [ Businger et al., 1971] [Beljaars and Holstag (1991)]   
        if y > 0
            %%% Stable Condition 
            a=1; b=0.667; c=5; d =0.35;
            Fim = -(a*y +b*(y-c/d).*exp(-d*y)+b*c/d) ;
            Fih = -( ((1+2*a*y/3)^1.5) +b*(y-c/d)*exp(-d*y) +(b*c/d-1));
        else
            %%% Unstable Condition 
            G=16; %% Dyer (1974)
            x=(1-G*y)^(1/4);
            Fim = log((0.5*(1+x^2)).*((0.5*(1+x))^2)) - 2*atan(x) + pi/2;
            Fih = 2*log((0.5*(1+x^2)));
        end
    end
end