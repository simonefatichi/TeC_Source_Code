%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Root Depth Dynamic %%%%%%%   %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[ZR95,Rf,ZR50,ZRmax]=Root_Depth_Dynamics(CASE_ROOT,B,Btm1,Rl,Zs,ZR95tm1,ZR50tm1,ZRmaxtm1,Bfac_day,Psan,Tdep,O,Soil_Param,a_root)
%%%%%%%%%%%%%%%%
%%%%% Potential update for rooting depth --
if  CASE_ROOT~= 1
    disp('ERROR: IN ORDER TO HAVE A VARIABLE ROOT DEPTH CASE_ROOT MUST BE 1')
    return
end
if  isnan(ZRmaxtm1)
    disp('ERROR: A maximum rooting depth must be specified when the dynamic root option is enabled')
    return
end
O= mean(O,1);
%%%%%%%%%%%%% Effective saturation in the deepest roots
row=1000; %%[kg/m^3]
Osat=Soil_Param.Osat;
Ohy =Soil_Param.Ohy;
rsd=Soil_Param.rsd; %%  density dry soil [kg/m^3]
n=length(Zs)-1;
for j=1:n
    if Zs(j+1) > ZR95tm1
        Se=(O(j)-Ohy(j))./(Osat(j)-Ohy(j));
        rsoil = rsd(j) + (O(j)-Ohy(j))*row; %% Soil Density [kg/m^3]
        break
    end
end
%%%%%%%
dtd= 1; %% [day]
BRoot=B(:,1,3); %%% Root biomass [gC m-2]
dBRoot=B(:,1,3)-Btm1(:,1,3); %%% Root biomass [gC m-2 day-1]
%Rl % %% root length index [m root / m^2 PFT]
SRL = Rl/BRoot; %%% [m root / gC ]
Tdep = mean(Tdep);
% % %%%%
% % %%%% Root depth Update
if BRoot > 0
    %%% --> ZR95, ZR50, ZRmax should be all updated
    DYR_ROOT_OPT=4;
    switch DYR_ROOT_OPT
        case 1
            % %ZR=0.5*(2*BRoot)^r;
            % ZR = 1283*BRoot^0.1713-1914; ZR(ZR<5)=5; %%[mm]
            % ZR95(cc)=ZR;
            %%%% Arora and Boer 2003
            b=10;
            alpha=0.8;
            ZR= 3/b*(BRoot*2/1000)^alpha;
        case 2
            %%%%  Gayler et al 2014
            rext_max= 20; %% mm/day
            ZRmaxP = 1000 ; %% mm
            Topt=25;
            Tmin=0;
            Tmax=35;
            alpha=log(2)*log((Tmax-Tmin)./(Topt-Tmin));
            fT= (2*(Tdep-Tmin)^alpha*(Topt-Tmin)^alpha-(Tdep-Tmin)^(2*alpha))./((Topt-Tmin)^(2*alpha));
            fT(Tdp<0)=0;
            %fO= min(1,4*OR);
            fO=Bfac_day;
            ZR = ZR95tm1 + dtd*(rext_max)*fT*fO*(1-ZR95tm1./ZRmaxP);
        case 3
            %%%% Lu et al 2019 adapted
            ZRmaxP = 1500 ; %% mm
            if dBRoot>0
                a= 0.001; %% vertical vs total fine root length growth
                rext=a*SRL*dBRoot*1000; %%% mm/day
            else
                a= 0.0; %% vertical vs total fine root length decrease
                rext=a*SRL*dBRoot*1000; %%% mm/day
            end
            ZR = ZR95tm1 + dtd*(rext)*(1-ZR95tm1./ZRmaxP);
        case 4
            %%%% Adapted based on Asseng et al 1997;  Lu et al 2019 ;
            %%%% 
            ZRmaxP = ZRmaxtm1 ; %% mm
            WST=1-Bfac_day; %% Water Deficit to increase elongation
            %%% Soil Resistance factor - Jones et al 1991 
            Bdx=1000*(1.6+0.4*Psan);
            Bdo=1000*(1.1+0.5*Psan);
            fres=(Bdx-rsoil)./(Bdx-Bdo);
            fres(fres<0)=0; fres(fres>1)=1;
            %%% Soil Aereation factor
            fAer=0.002729*Se.^-56.03;  fAer(fAer>1)=1;
            %%% Soil Temperature factor
            ftem=(-0.0002667*Tdep.^3 + 0.007143*Tdep.^2 + 0.01381*Tdep).*(Tdep<20).*(Tdep>0) + 1*(Tdep>=20).*(Tdep<=25) + (-0.01.*Tdep.^2 + 0.5*Tdep -5.25).*(Tdep>25).*(Tdep<=35) ;
            ftem(ftem<0)=0; ftem(ftem>1)=1;
            %%%
            fr=min([fAer,fres,ftem]);
            if dBRoot>0
                a= a_root; %% vertical vs total fine root length growth
                rext=a*SRL*dBRoot*1000; %%% mm/day
            else
                a= a_root/10; %% vertical vs total fine root length decrease
                rext=a*SRL*dBRoot*1000; %%% mm/day
            end
            ZR = ZR95tm1 + dtd*(rext)*(1+WST)*fr*(1-ZR95tm1./ZRmaxP);
    end
    %%%
    ZR(ZR<5)=5; %%[mm]
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ZR95=ZR;
    ZR50 = ZR50tm1;
    ZRmax = ZRmaxtm1;
    [Rf]=Root_Fraction_General(Zs,CASE_ROOT,ZR95,ZR50,0*ZR95,0*ZR50,ZRmax,0*ZRmax);
else
    ZR95=0;
    ZR50 = ZR50tm1;
    ZRmax = ZRmaxtm1;
    [Rf]=Root_Fraction_General(Zs,CASE_ROOT,ZR95,ZR50,0*ZR95,0*ZR50,ZRmax,0*ZRmax);
end
return
%%%%%%
