%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Albedo_Snow_Properties     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[snow_alb,tau_sno,e_sno]=Albedo_Snow_Properties(dt,SWE,h_S,Ts,SWEtm1,tau_snotm1,snow_albtm1,Th_Pr_sno,Pr_sno_day,Aice,Deb_Par,Cdeb,Cice,Ta_day,Pr_sno)
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% References [Oleson, et al., 2004]
%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% INPUT
%%% dt = [s]
%%% ros = snow density [kg/m^3]
%%% D =  snow depth [m]
%%% h_S solar heigth [rad]
%%% Ts
%%% rostm1
%%% Dtm1
%%% tau_snotm1 [] Relative Age of snow (t-1)
%%% MTa = maximum air temperature for the day divided so hourly
%%% MTatm1 = maximum air temperature for the day divided so hourly at time
%%% before
%%% SWE = snow water equivalen [mm w.e.]
%%% Aice = ice albedo
%%% Deb_Par.alb = debris albedo
%%% Cdeb = debris-covered glacier or not
%%% Cice = glacier or not
%%% Pr_sno = snowfall
%%% OUTPUT
% tau_sno [] %% Relative Age of snow
%snow_alb.dir_vis
%snow_alb.dif_vis
%snow_alb.dir_nir
%snow_alb.dif_nir
%%AMTa_out = accumulated maximum air temperature as output in tau_sno 
%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%
e_sno = 0.97; %%% Snow emissivity
%%%%%%%%%%
%Choose method based on if glacier or not
% if Cice==1
%     ANS=4; %Use snow on glacier albedo
% else
%     ANS=2; %Use previous method for most surfaces
% end
ANS=2;
%%%%%%%
switch ANS
    case 1
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if h_S <= 0
            h_S = 0; %% night
        end
        mu = sin(h_S); %mu is the cosine of the zenith angle of the incident beam
        b=2.0; %% Dickinson et al., 1993
        if mu <=0.5
            f_mu =(1+1/b)/(1+mu*2*b) -1/b;
            f_mu(f_mu<0)=0;
        else
            f_mu=0;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%Reference [Oleson et al., 2004]
        %%% vis nir
        C=[0.2 0.5]; %%% Dickinson 1993
        asno_bas=[0.95 0.65];
        to = 1*(10^-6); %%[1/s]
        Tf = 273.15; %[K]
        Ts_k = Ts +  273.15; %% [K]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%
        r1=exp(5000*(1/Tf - 1/Ts_k));
        r2=(r1^10); r2(r2>1)=1;
        r3 = 0.3;
        %%%%%%%%%%%%%%%%%%%%%
        if SWE <= 800 %%% Snow Mass [kg/m^2] = Snow water equivalent
            dtsno = to*(r1+r2+r3)*dt;
        else
            dtsno = 0;
        end
        tau_sno = (tau_snotm1 + dtsno)*(1-0.1*(SWE - SWEtm1));
        if  (tau_sno <= 0) || ((SWE - SWEtm1) > 10);%% 10 mm to restore
            tau_sno=0;
        end
        if tau_sno > 1000;
            tau_sno=1000;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%
        Fage = 1 - 1/(1+tau_sno);
        %%%%%%%%%%%%%%%%%%%%%%%%%%%
        snow_alb.dir_vis = asno_bas(1)+0.4*f_mu*(1-asno_bas(1));
        snow_alb.dir_nir = asno_bas(2)+0.4*f_mu*(1-asno_bas(2));
        snow_alb.dif_vis = (1-C(1)*Fage)*asno_bas(1);
        snow_alb.dif_nir = (1-C(2)*Fage)*asno_bas(2);
    case 2
        %%%%%
        Asnotm1= snow_albtm1.dir_vis;
        tau_sno=0;
        %%%%% ----->>>>>
        %%%%%%% OLD METHOD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% ALBEDO UPDATING
        %%% Reference [Douville 1995; Strack et al., 2004 ; Essery et al., 1999]
        %%% Parametrs
        %%%%%%%%%%%%%%
        %Th_Pr_sno = ##; %% [mm]
        %%%%%%%%%%%%%%%
        Amin =0.5; tau_A = 0.008; tau_1 = 86400; %[s]
        tau_F = 0.24; %
        %tau_F = 0.48;
        Amax = 0.85;
        if Ts < -0.01 %  A  % albedo of snow
            %%% frozen
            Asno = Asnotm1 -tau_A*dt/tau_1;
        else
            % Melting Season
            Asno= (Asnotm1 - Amin)*exp(-tau_F*dt/tau_1) + Amin; %%
        end
        %%% new snowafall
        if (Pr_sno_day >= Th_Pr_sno) || (SWEtm1) == 0 %% Check
            Asno = Amax;
        end
        if Asno < Amin
            Asno = Amin ;
        end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        snow_alb.dir_vis = Asno;
        snow_alb.dir_nir = Asno;
        snow_alb.dif_vis = Asno;
        snow_alb.dif_nir = Asno;
        %%%%
    case 3
        if SWE > 45
            Asno =   max(0.35,0.79-0.0018*Ts-0.0045*Ws);
        else
            Asno  = max(0.585,0.616-0.001*Ts - 0.0076*Ws + 0.01*SND);
        end
        if SWE < 8
            Asno = 0.565;
        end
        snow_alb.dir_vis = Asno;
        snow_alb.dir_nir = Asno;
        snow_alb.dif_vis = Asno;
        snow_alb.dif_nir = Asno;
        tau_sno=0;
    case 4
        %   Snow albedo function as described in Brock et al. (2000)
        AMTatm1 = tau_snotm1;
        %%%%
        MTa=max(Ta_day)/24;
        MTa(MTa<0)=0; %Only accumulate positive maximum temperature
        %%%%%%%%
        Br_Param.a = 0.713; %In deep snow equation
        Br_Param.b = 0.112; %In deep snow equation
        Br_Param.c = 0.442; %In shallow snow equation
        Br_Param.d = -0.058; %In shallow snow equation
        Br_Param.e = 0.024; %Length scale SWE
        
        %%%%%%%%%%
        SWE_m=SWE/1000; %Convert SWE to snow depth m w.e.
        
        %Calculate maximum accumulated Ta
        if Pr_sno==0 %No snow this hour
            ATa=MTa + AMTatm1; %Accumulated max T this timestep
            AMTa_out=ATa; %to go back as output for next timestep
        else %Pr_sno>0
            ATa=MTa; %Reset to only be accumulated max T for this hour as new snow
            AMTa_out=ATa;
        end
        
        %Determine underlying surface albedo
        if Cdeb==0
            a_u=Aice; %clean ice
        else
            a_u=Deb_Par.alb; %debris-covered ice
        end
        
        %Deep snow albedo
        a_ds = Br_Param.a - Br_Param.b*log10(ATa);
        
        %Shallow snow albedo
        a_ss = a_u + Br_Param.c*exp(Br_Param.d*ATa);
        
        %Final snow albedo
        Asno = (1-exp(-SWE_m/Br_Param.e))*a_ds + exp(-SWE_m/Br_Param.e)*a_ss;
        
        if Asno>0.85
            Asno=0.85;
        end
        snow_alb.dir_vis = Asno;
        snow_alb.dir_nir = Asno;
        snow_alb.dif_vis = Asno;
        snow_alb.dif_nir = Asno;
        tau_sno=AMTa_out;
end
return