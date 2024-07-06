%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Undercanopy_Resistence     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%REFERENCES %%   - Mahat et al 2013
function[rap_H,rap_L,rb_H,rb_L,Ws_und]=Undercanopy_Leaf_Resistence2(Ws,Ta,Ts,Ccrown,hc_H,hc_L,LAI_H,LAI_L,d_leaf_H,d_leaf_L,...
    zatm,disp_h,zom,zom_under,SND,disp_h_H,zom_H,disp_h_L,zom_L)
%%%INPUTS
%Ta = %% air temperature [°C] --
%Ts =  ; %% surface temperature [°C] --
%z= %  reference height [m] % ---
%%% zatm = measurement height
%%% disp_h = displacement height
%zom
%zom_under = roughness eddy diffusivities for momentum [m] undercanopy land use
%zoh_under = %% roughness  eddy diffusivities for heat  [m]  undercanopy land use
u = Ws; %% [m/s] wind speed  ---
g=9.81; %%[m/s2]
k= 0.41; %% Von Karman Constant
%SND snow height
%%% vectors
%hc_h hc_L = canpy heights [m]
% LAI_H LAI_L = [-] [Leaf Area Index ]
% d_leaf [m] [Leaf width]
%%% disp_h_H  = displacement height High Veg.
%zom_H
%%% disp_h_L  = displacement height High Veg.
%zom_L
%%%%%%%%%
%%% Generic Plot --
%z = zatm;
%d = disp_h;
%zom
%zom_under
%%%%
%%%%%
d = disp_h ; %% Zero plane displacement [m]
z = zatm; %% Measurement Height [m]
%%% OUTPUTS
%%% rap_H rap_L Undercanopy resistence
rb_H=zeros(1,length(Ccrown)); rb_L=zeros(1,length(Ccrown));
rap_H=zeros(1,length(Ccrown)); rap_L=zeros(1,length(Ccrown));
Ws_und=zeros(1,length(Ccrown)); 
%%%%
%%% All expression hc =max([hc_H(i),d+zom+0.01,0.05]); are for numerical
%%% stability 
%%%%%%%%%%%%%%%%%%%
if zom > max(max(zom_H),max(zom_L))
    %%% It means some other elements is above the vegetation and wind needs to be transferred below
    h_other =zom/0.123;
    d_other = h_other*(2/3);
    us =  k*u/log((z-d_other)/zom); %%% Friction Velocity [m/s]
    u_above = (us/k)*log((h_other-d_other)/zom); %% Wind Speed top  [m/s]
    alpha  = log(u/u_above)/(z/h_other -1); %% Attenuation Coefficient
    hc =max(max(hc_H),max(hc_L))+2;
    u_below = u_above*exp(-alpha*(1-hc/h_other)); %% Wind at reference height
    %%%
    z=hc;
    u=u_below;
    zom=max(max(zom_H),max(zom_L));
    d= max(max(disp_h_H),max(disp_h_L));
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Undercanopy and Leaf  resitence
for i=1:length(Ccrown)
    %%%%%%%%
    if disp_h_H(i) > 0 && disp_h_H(i) < d
        %%% High Veg below High Veg
        [hc,pi] = max(hc_H); %%% height vegetation above
        hc = max([max(hc_H),d+zom+0.01,0.05]); %%% height vegetation above
        d2 = disp_h_H(i) ;
        zom2  = zom_H(i);
        hc2 =max([hc_H(i),d2+zom2+0.01,0.05]);
        LAI_UP = LAI_H(pi); %% LAI vegetation above
        LAI_DW = LAI_H(i) ;
        d_leaf_UP =d_leaf_H(pi);
        d_leaf_DW =d_leaf_H(i);
        OPT_VEG = 1;
    else
        %%%% High and Low vegetation stucked --
        if  (hc_L(i) >0) && (hc_H(i)>0)
            OPT_VEG = 3;
            %%%%%%%%%
            hc = max([max(hc_H),d+zom+0.01,0.05]); %%% height vegetation above
            d2 = disp_h_L(i) ;
            zom2  = zom_L(i);
            hc2 =max([hc_L(i),d2+zom2+0.01,0.05]);
            LAI_UP = LAI_H(i); %% LAI vegetation above
            LAI_DW = LAI_L(i) ;
            d_leaf_UP =d_leaf_H(i);
            d_leaf_DW =d_leaf_L(i);
        else
            if disp_h_L(i) > 0 && disp_h_L(i) < d
                %%% Low Veg below Low or High Veg
                [hc11,pi1] = max(hc_H); %%% height vegetation above
                [hc12,pi2] = max(hc_L); %%% height vegetation above
                [hc,pi3] = max([hc11,hc12]);
                if pi3==1
                    LAI_UP = LAI_H(pi1); %% LAI vegetation above
                    d_leaf_UP =d_leaf_H(pi1);
                    hc = max([max(hc_H),d+zom+0.01,0.05]); %%% height vegetation above
                else
                    LAI_UP = LAI_L(pi2); %% LAI vegetation above
                    d_leaf_UP =d_leaf_L(pi2);
                    hc = max([max(hc_L),d+zom+0.01,0.05]); %%% height vegetation above
                end
                d2 = disp_h_L(i) ;
                zom2  = zom_L(i);
                hc2 =max([hc_L(i),d2+zom2+0.01,0.05]);
                LAI_DW = LAI_L(i) ;
                d_leaf_DW =d_leaf_L(i);
                OPT_VEG = 2;
            else
                if  (hc_L(i) == 0) && (hc_H(i) == 0)
                    %%% No Vegetation
                    OPT_VEG = 5 ;
                else
                    %%% Single Vegetation
                    OPT_VEG = 4;
                end
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%
    switch OPT_VEG
        case {1,2,3}
            %%% Aerodynamic lower vegetation in case of two vegetations
            if SND >= (d+zom)
                rap_UP=0;
                rb_UP=0;
                rap_DW=0;
                rb_DW=0;
                ums=Ws; ums2=0;
            else
                %%%%%
                us =  k*u/log((z-d)/zom); %%% Friction Velocity Canopy above [m/s]
                u_hc = (us/k)*log((hc-d)/zom); %% Wind Speed top Canopy above [m/s]
                alpha  = log(u/u_hc)/(z/hc -1); %% Attenuation Coefficient
                %%% Yi et al 2008
                alpha  = alpha*LAI_UP*0.5/2;
                %%%%%
                zms = 2 + max(hc2,SND) ; %% reference height above the surface
                %%%%
                zms_b= zms; 
                ums_hb = u_hc*exp(-alpha*(1-zms_b/hc)); %% Wind at reference height for Low vegetation
                if zms > d+zom
                    zms=d+zom;
                end
                %%%
                ums = u_hc*exp(-alpha*(1-zms/hc)); %% Wind at reference height for Low vegetation after correction used for High-undercanopy rap only  
                %%%%
                %%%
                Kh = k^2*u*(hc-d)/log((z-d)/zom);
                rap_n = hc*exp(alpha)/(Kh*alpha)*(exp(-alpha*(zms/hc))-exp(-alpha*((d+zom)/hc))) + 1/(k^2*ums)*log(zms/zom2)^2;
                %%% Stability correction
                Ri = (g*(Ta-Ts)*zms)/(ums^2*(0.5*(Ta+Ts)+273.15)); %% [-]
                Ri(Ri>0.16)=0.16; %% Max. Stability
                if Ri < 0 %% unstable
                    rap = rap_n/((1-5*Ri)^(3/4));
                else %% Stable
                    rap = rap_n/((1-5*Ri)^2);
                end
                rap_UP=rap;
                rb_UP=Leaf_BR(u_hc,Ts,Ta,d_leaf_UP,alpha);
                %%%%%%%%%%%%%%%%%%%%%%%%
                if SND > (d2+zom2)
                    rap_DW=0;
                    rb_DW=0;
                    ums2=Ws;
                else
                    us2 =  k*ums_hb/log((zms_b-d2)/zom2); %%% Friction Velocity Canopy above [m/s]
                    u_hc2 = (us2/k)*log((hc2-d2)/zom2); %% Wind Speed top Canopy [m/s]
                    alpha2  = log(ums_hb/u_hc2)/(zms_b/hc2 -1); %% Attenuation Coefficient
                    %%% Yi et al 2008
                    alpha2  = alpha2*LAI_DW*0.5/2;
                    %%%%%
                    zms2 = 2 +SND ; %% reference height above the surface
                    if zms2 > d2+zom2
                        zms2=d2+zom2;
                    end
                    %%%
                    ums2 = u_hc2*exp(-alpha2*(1-zms2/hc2)); %% Wind at reference height for Low vegetation
                    %%%
                    Kh2 = k^2*ums_hb*(hc2-d2)/log((zms_b-d2)/zom2);
                    rap_n2 = hc2*exp(alpha2)/(Kh2*alpha2)*(exp(-alpha2*(zms2/hc2))-exp(-alpha2*((d2+zom2)/hc2))) + 1/(k^2*ums2)*log(zms2/zom_under)^2;
                    %%% Stability correction
                    Ri = (g*(Ta-Ts)*zms2)/(ums2^2*(0.5*(Ta+Ts)+273.15)); %% [-]
                    Ri(Ri>0.16)=0.16; %% Max. Stability
                    if Ri < 0 %% unstable
                        rap2 = rap_n2/((1-5*Ri)^(3/4));
                    else %% Stable
                        rap2 = rap_n2/((1-5*Ri)^2);
                    end
                    rap_DW=rap2;
                    rb_DW=Leaf_BR(u_hc2,Ts,Ta,d_leaf_DW,alpha2);
                end
            end
            %%%%%%%%%%
            if OPT_VEG == 1
                %%% High Veg below High Veg
                rap_H(i)=rap_DW;
                rap_L(i)=0;
                rb_H(i)= rb_DW;
                rb_L(i)= 0;
                Ws_und(i)= ums2;
            end
            if OPT_VEG == 2
                %%% Low Veg below Low Veg
                rap_H(i)=0;
                rap_L(i)=rap_DW;
                rb_H(i)= 0;
                rb_L(i)=rb_DW;
                Ws_und(i)= ums2;
            end
            if OPT_VEG == 3
                %%% Both Vegetations 
                rap_H(i)=rap_UP;
                rap_L(i)=rap_DW;
                rb_H(i)= rb_UP;
                rb_L(i)=rb_DW;
                Ws_und(i) = ums2;
            end
        case 4
            %%%%%%%%%%%%%%%% ONLY ONE VEGETATION LAYER %%%%%%%%
            if hc_H(i) > 0 %% High Vegetation
                hc =max([hc_H(i),d+zom+0.01,0.05]);
                LAI = LAI_H(i);
                d_leaf = d_leaf_H(i);
            else %%% Low vegetation
                hc = max([hc_L(i),d+zom+0.01,0.05]);
                LAI = LAI_L(i);
                d_leaf = d_leaf_L(i);
            end
            if SND >= (d+zom)
                rb=0;
                rap=0;
                ums=Ws;
            else
                us =  k*u/log((z-d)/zom); %%% Friction Velocity Canopy above [m/s]
                u_hc = (us/k)*log((hc-d)/zom); %% Wind Speed top Canopy [m/s]
                %%% alpha = 0.6 -1.5 -- alpha = LAI/2;
                %%% Match exponential and logarithmic
                alpha = log(u/u_hc)/(z/hc -1); %% Attenuation Coefficient
                %%% Yi et al 2008
                alpha = alpha*LAI*0.5/2;
                % zom_under
                zms = 2 + SND ; %% reference height above the surface
                if zms > d+zom
                    zms=d+zom;
                end
                %%%
                ums = u_hc*exp(-alpha*(1-zms/hc)); %% Wind at reference height below canopy
                %%%
                Kh = k^2*u*(hc-d)/log((z-d)/zom);
                rap_n = hc*exp(alpha)/(Kh*alpha)*(exp(-alpha*(zms/hc))-exp(-alpha*((d+zom)/hc))) + 1/(k^2*ums)*log(zms/zom_under)^2;
                %%% Stability correction
                Ri = (g*(Ta-Ts)*zms)/(ums^2*(0.5*(Ta+Ts)+273.15)); %% [-]
                Ri(Ri>0.16)=0.16; %% Max. Stability
                if Ri < 0 %% unstable
                    rap = rap_n/((1-5*Ri)^(3/4));
                else %% Stable
                    rap = rap_n/((1-5*Ri)^2);
                end
                [rb]=Leaf_BR(u_hc,Ts,Ta,d_leaf,alpha);
            end
            if hc_H(i) > 0 %% High Vegetation
                rap_H(i)=rap;
                rb_H(i)=rb;
                rap_L(i)=0;
                rb_L(i)=0;
                Ws_und(i)=ums;
            else %%% Low vegetation
                rap_L(i)=rap;
                rb_L(i)=rb;
                rap_H(i)=0;
                rb_H(i)=0;
                Ws_und(i)=ums;
            end
        case 5
            %%% NO VEGETATION
            rap_H(i) = 0; rap_L(i) = 0 ;
            rb_H(i)=0; rb_L(i)=0;
            Ws_und(i) =Ws ; 
    end
end
return

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Leaf_Boundary_Resistence   %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%REFERENCES %%   - Vesala 1998 -- Chodhury and Monteith 1988 - Ivanov 2008
%%% Leigh et al 2011 ; Leuning et al 1995
function[rb]=Leaf_BR(u_hc,Ts,Ta,d_leaf,alpha)
%%%INPUTS
% Ta = %% air temperature [°C] --
% Ts =  ; %% surface temperature [°C] --
% Pre = %pressure [Pa]--
% z= %  reference height [m] % ---
% zoh = %% roughness  eddy diffusivities for heat  [m]
% zom = roughness eddy diffusivities for momentum [m]
% u = Ws; %% [m/s] wind speed  ---
% hc = canpy height [m]
% d_leaf = [cm] Leaf Dimension
% LAI = [Leaf Area Index ]
%%% OUTPUTS
%%% rb % [s/m] Leaf Boundary Resistence
%%% PARAMETERS
d_leaf = d_leaf/100; %% [m]
%k= 0.4; %% Von Karman Constant
a = 0.005; %% [m/s^0.5] --  Chodhury and Monteith 1988
%d = disp_h; %% Zero plane displacement [m]
%z= zatm ; %% Measurement Height [m]
%%% Hypothesis Logaritmic distribution of wind speed
%us =  k*u/log((z-d)/zom); %%% Friction Velocity  [m/s]
%u_hc = (us/k)*log((hc-d)/zom); %% Wind Speed top Canopy [m/s]
%%%
%alpha = log(u/u_hc)/(z/hc -1); %% Attenuation Coefficient
%alpha = 0.5*alpha*LAI/2;
%%% u(z) = u_hc*exp(alpha*((z/hc-1));
%%%%%%%%% Expression of Leaf Boundary Layer Resistence
%%% gb(z)= a*(u(z)/d_leaf)^0.5; %% [Jones 1992]  Integral between 0-LAI_TOT with LAI(z) = LAI_TOT*z/Hc
gb = (2*a/alpha)*((u_hc/d_leaf)^0.5)*(1-exp(-alpha/2)); %% [m/s]
%rb= (1/(gb*LAI)); %%   Leaf Boundary Layer Resistnce [s/m] Plant
%rb=1/gb; %%   Leaf Boundary Layer Resistnce [s/m] one-sided for unit leaf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Expression for free convection  (Leuning 1995, Monteith 1973)
Dh= 1.9*1e-5; %% [m2/s]
Gr= 1.6*1e+8*(Ts-Ta)*(d_leaf^3).*(Ts>Ta); %% [-]
gb_free = 0.5*Dh*Gr^(0.25)/d_leaf; %%[m/s]
%%%%%
gb=gb+gb_free;
rb=1/gb; %%   Leaf Boundary Layer Resistnce [s/m] one-sided for unit leaf
return