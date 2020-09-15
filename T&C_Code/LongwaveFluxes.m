%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction LongwaveFluxes               %%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Labs_vegH,Labs_vegL,Labs_soi]=LongwaveFluxes(Tg,TvHsun,TvHshd,TvLsun,TvLshd,Tsno,...
    Latm,SvF,e_gr,LAI_H,SAI_H,LAI_L,SAI_L,FshdH,FsunH,FshdL,FsunL,hc_L,SnoDep,e_sno,Csno)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% Bonan 1996; Ivanov 2008;
%%%%%%% INPUT
% Latm Incoming LongWave Radiation [W/m^2]
% Ts surface radiative temperature [°C]
% SvF sky view factor
% e_gr ground emissivity
%%% LAI_H
%%% SAI_H
%%% LAI_L
%%% SAI_L
%%% FshdH
%%% FsunH
%%% FshdL
%%% FsunL
%%% hc_L Height Vegetation L [m]
%%% SnoDep [m]
%%%%%%% OUTPUT
%%Labs_vegH,  Long Wave absorbed by the vegetation H [W/m^2]
%%Labs_vegL  % Long Wave absorbed by the vegetation L [W/m^2]
%%Labs_soi % % Long Wave absorbed by the ground [W/m^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Tg_k =  Tg +273.15; %% % Tg ground temperature [K]
Tsno_k =  Tsno +273.15; %% % Tg ground temperature [K]
TvLsun_k =  TvLsun +273.15; %%% TvL vegetation L temperature sunlit [K]
TvHsun_k =  TvHsun +273.15; %%% TvH vegetation H temperature sunlit [K]
TvLshd_k =  TvLshd +273.15; %%% TvL vegetation L temperature shadow [K]
TvHshd_k =  TvHshd +273.15; %%% TvH vegetation H temperature shadow [K]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
e_veg_H= 1-exp(-(LAI_H+SAI_H)); %% Vegetation emissivity
e_veg_L= 1-exp(-(LAI_L+SAI_L)); %% Vegetation emissivity
%%e_gr= 0.96; % Ground emissivity
%%%%%% Hip Absorpitivity = Emissivity
a_veg_H= e_veg_H;
a_veg_L= e_veg_L;
a_gr= e_gr;
a_sno= e_sno;
%%%%%%%%%%%%%%
SvF_H = (1-a_veg_H)*SvF + a_veg_H; 
SvF_L = (1-a_veg_L)*SvF + a_veg_L; 
SvF_HL = (1-a_veg_L)*SvF_H + a_veg_L; 
%%%%%%%%%%%%%%%%%%%%%%%%%%
sigmaSB = 5.6704e-8; % Stefan-Boltzmann constant [W/m^2 K4]  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Latm = Latm*SvF; %%% Incoming LongWave Radiation [W/m^2]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if Csno == 0
    if (LAI_H+SAI_H)>0 && (LAI_L+SAI_L)>0 %%%%%%% BOTH VEGETATION ARE PRESENT
        LvegH_dn = FsunH*e_veg_H*sigmaSB*(TvHsun_k)^4  + FshdH*e_veg_H*sigmaSB*(TvHshd_k)^4 + ...
            (1-a_veg_H)*Latm; %%% Long Wave from vegetation H down [W/m^2]
        LvegL_dn = FsunL*e_veg_L*sigmaSB*(TvLsun_k)^4  + FshdL*e_veg_L*sigmaSB*(TvLshd_k)^4+ ...
            (1-a_veg_L)*LvegH_dn; %%% Long Wave from vegetation L down [W/m^2]
        Lground =  e_gr*sigmaSB*(Tg_k)^4+(1-a_gr)*LvegL_dn; %% Long Wave from soil [W/m^2]
        LvegL_up = FsunL*e_veg_L*sigmaSB*(TvLsun_k)^4+  FshdL*e_veg_L*sigmaSB*(TvLshd_k)^4 + ...
            (1-a_veg_L)*Lground; %%% Long Wave from vegetation L up [W/m^2]
        LvegH_up = FsunH*e_veg_H*sigmaSB*(TvHsun_k)^4+  FshdH*e_veg_H*sigmaSB*(TvHshd_k)^4  + ...
            (1-a_veg_H)*LvegL_up; %%% Long Wave from vegetation H up [W/m^2]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Labs_vegH= Latm - LvegH_dn - SvF*LvegH_up + SvF_H*LvegL_up; %%% Long Wave absorbed by the vegetation H [W/m^2]
        Labs_vegL= LvegH_dn - LvegL_dn - SvF_H*LvegL_up + SvF_HL*Lground; %%% Long Wave absorbed by the vegetation L [W/m^2]
        Labs_soi= LvegL_dn - SvF_HL*Lground; %%% % Long Wave absorbed by the ground [W/m^2]
        %%%%%%%%%%%%%%%%%%
    else
        if (LAI_H+SAI_H)>0  %%% ONLY HIGH VEGETATION
            LvegH_dn = FsunH*e_veg_H*sigmaSB*(TvHsun_k)^4  + FshdH*e_veg_H*sigmaSB*(TvHshd_k)^4 + ...
                (1-a_veg_H)*Latm; %%% Long Wave from vegetation H down [W/m^2]
            Lground =  e_gr*sigmaSB*(Tg_k)^4+(1-a_gr)*LvegH_dn; %% Long Wave from soil [W/m^2]
            LvegH_up = FsunH*e_veg_H*sigmaSB*(TvHsun_k)^4+  FshdH*e_veg_H*sigmaSB*(TvHshd_k)^4  + ...
                (1-a_veg_H)*Lground ; %%% Long Wave from vegetation H up [W/m^2]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Labs_vegH= Latm - LvegH_dn - SvF*LvegH_up + SvF_H*Lground; %%% Long Wave absorbed by the vegetation H [W/m^2]
            Labs_vegL= 0;
            Labs_soi= LvegH_dn  - SvF_H*Lground; %%% % Long Wave absorbed by the ground [W/m^2]
            %%%%%%%%%%%%%%%%%%
        else  %%% ONLY LOW VEGETATION
            LvegL_dn = FsunL*e_veg_L*sigmaSB*(TvLsun_k)^4  + FshdL*e_veg_L*sigmaSB*(TvLshd_k)^4+ ...
                (1-a_veg_L)*Latm; %%% Long Wave from vegetation L down [W/m^2]
            Lground =  e_gr*sigmaSB*(Tg_k)^4+(1-a_gr)*LvegL_dn; %% Long Wave from soil [W/m^2]
            LvegL_up = FsunL*e_veg_L*sigmaSB*(TvLsun_k)^4+  FshdL*e_veg_L*sigmaSB*(TvLshd_k)^4 + ...
                (1-a_veg_L)*Lground; %%% Long Wave from vegetation L up [W/m^2]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Labs_vegH=0;
            Labs_vegL= Latm - LvegL_dn - SvF*LvegL_up + SvF_L*Lground; %%% Long Wave absorbed by the vegetation L [W/m^2]
            Labs_soi= LvegL_dn - SvF_L*Lground; %%% % Long Wave absorbed by the ground [W/m^2]
            %%%%%%%%%%%%%%%%%%
        end
    end
else
    %%%%%%%%%%%% SNOW COVER PRESENCE 
    if (LAI_H+SAI_H)>0 && (LAI_L+SAI_L)>0 && (hc_L > SnoDep) %%%%%%% BOTH VEGETATION ARE PRESENT
        LvegH_dn = FsunH*e_veg_H*sigmaSB*(TvHsun_k)^4  + FshdH*e_veg_H*sigmaSB*(TvHshd_k)^4 + ...
            (1-a_veg_H)*Latm; %%% Long Wave from vegetation H down [W/m^2]
        LvegL_dn = FsunL*e_veg_L*sigmaSB*(TvLsun_k)^4  + FshdL*e_veg_L*sigmaSB*(TvLshd_k)^4+ ...
            (1-a_veg_L)*LvegH_dn; %%% Long Wave from vegetation L down [W/m^2]
        Lsno =  e_sno*sigmaSB*(Tsno_k)^4+(1-a_sno)*LvegL_dn; %% Long Wave from soil [W/m^2]
        LvegL_up = FsunL*e_veg_L*sigmaSB*(TvLsun_k)^4+  FshdL*e_veg_L*sigmaSB*(TvLshd_k)^4 + ...
            (1-a_veg_L)*Lsno; %%% Long Wave from vegetation L up [W/m^2]
        LvegH_up = FsunH*e_veg_H*sigmaSB*(TvHsun_k)^4+  FshdH*e_veg_H*sigmaSB*(TvHshd_k)^4  + ...
            (1-a_veg_H)*LvegL_up; %%% Long Wave from vegetation H up [W/m^2]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        Labs_vegH= Latm - LvegH_dn - SvF*LvegH_up + SvF_H*LvegL_up; %%% Long Wave absorbed by the vegetation H [W/m^2]
        Labs_vegL= LvegH_dn - LvegL_dn - SvF_H*LvegL_up + SvF_HL*Lsno; %%% Long Wave absorbed by the vegetation L [W/m^2]
        Labs_soi= LvegL_dn - SvF_HL*Lsno; %%% % Long Wave absorbed by the ground [W/m^2]
        %%%%%%%%%%%%%%%%%%
    else
        if (LAI_H+SAI_H)>0  %%% ONLY HIGH VEGETATION
            LvegH_dn = FsunH*e_veg_H*sigmaSB*(TvHsun_k)^4  + FshdH*e_veg_H*sigmaSB*(TvHshd_k)^4 + ...
                (1-a_veg_H)*Latm; %%% Long Wave from vegetation H down [W/m^2]
            Lsno =  e_sno*sigmaSB*(Tsno_k)^4+(1-a_sno)*LvegH_dn; %% Long Wave from soil [W/m^2]
            LvegH_up = FsunH*e_veg_H*sigmaSB*(TvHsun_k)^4+  FshdH*e_veg_H*sigmaSB*(TvHshd_k)^4  + ...
                (1-a_veg_H)*Lsno ; %%% Long Wave from vegetation H up [W/m^2]
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            Labs_vegH= Latm - LvegH_dn - SvF*LvegH_up + SvF_H*Lsno; %%% Long Wave absorbed by the vegetation H [W/m^2]
            Labs_vegL= 0;
            Labs_soi= LvegH_dn  - SvF_H*Lsno; %%% % Long Wave absorbed by the ground [W/m^2]
            %%%%%%%%%%%%%%%%%%
        else  %%% ONLY LOW VEGETATION Not Cover by snow 
            if (hc_L > SnoDep)
                LvegL_dn = FsunL*e_veg_L*sigmaSB*(TvLsun_k)^4  + FshdL*e_veg_L*sigmaSB*(TvLshd_k)^4+ ...
                    (1-a_veg_L)*Latm; %%% Long Wave from vegetation L down [W/m^2]
                Lsno =  e_sno*sigmaSB*(Tsno_k)^4+(1-a_sno)*LvegL_dn; %% Long Wave from soil [W/m^2]
                LvegL_up = FsunL*e_veg_L*sigmaSB*(TvLsun_k)^4+  FshdL*e_veg_L*sigmaSB*(TvLshd_k)^4 + ...
                    (1-a_veg_L)*Lsno; %%% Long Wave from vegetation L up [W/m^2]
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                Labs_vegH=0;
                Labs_vegL= Latm - LvegL_dn - SvF*LvegL_up + SvF_L*Lsno; %%% Long Wave absorbed by the vegetation L [W/m^2]
                Labs_soi= LvegL_dn - SvF_L*Lsno; %%% % Long Wave absorbed by the ground [W/m^2]
                %%%%%%%%%%%%%%%%%%
            else %%%%%% ONLY LOW VEGETATION Cover by snow  
                Lsno =  e_sno*sigmaSB*(Tsno_k)^4+(1-a_sno)*Latm; %% Long Wave from soil [W/m^2]
                Labs_vegH=0;
                Labs_vegL=0;
                Labs_soi= Latm - SvF*Lsno; %%% % Long Wave absorbed by the ground [W/m^2]
            end
        end
    end
end


