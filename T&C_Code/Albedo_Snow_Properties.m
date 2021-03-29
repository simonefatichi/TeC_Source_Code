%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Albedo_Snow_Properties     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[snow_alb,tau_sno,e_sno]=Albedo_Snow_Properties(dt,SWE,h_S,Ts,Ta,SWEtm1,tau_snotm1,snow_albtm1,Th_Pr_sno,Pr_sno_day,Aice,Deb_Par,Cdeb,Cice,Ta_day,Pr_sno,Pr_liq,ros,N)

%%%%%%%%%%%%%%%%%%%%%%%%%
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
Pr_sno_day= Pr_sno_day+Pr_sno*dt/3600; %% [mm] 
Ta_day= [reshape(Ta_day(2:length(Ta_day)),1,length(Ta_day)-1), Ta]; %%[°C]
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
ANS=5;
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
        if  (tau_sno <= 0) || ((SWE - SWEtm1) > 10) %% 10 mm to restore
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
        if (Pr_sno_day >= Th_Pr_sno) || (SWEtm1==0) %Snow greater than threshold or there was no snow on the timestep before
            ATa=MTa; %Reset to only be accumulated max T for this hour as new snow
            AMTa_out=ATa;
        else %Pr_sno<threshold
            ATa=MTa + AMTatm1; %Accumulated max T this timestep
            AMTa_out=ATa; %to go back as output for next timestep
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
        elseif isnan(Asno)
            Asno=a_u; %So if=NaN (only occurs under conditions of Csno=1, SWE=0 and ATa=0, so its cold but the snow has melted, therefore assign the underlying albedo)
        end 
        
        snow_alb.dir_vis = Asno;
        snow_alb.dir_nir = Asno;
        snow_alb.dif_vis = Asno;
        snow_alb.dif_nir = Asno;
        tau_sno=AMTa_out;
    case 5
        %%%%%%%% Ding et al 2017 Albedo Parameterization
        %   Snow albedo function as described in Ding et al. (2017)
        
        %%%%%
        Asnotm1= tau_snotm1;
        
        Amax = 0.85;
        Amin = 0.5;
        if Cice == 1
            if Cdeb == 1
                Amin =  Deb_Par.alb;
            else
                Amin = Aice;
            end
        end
        
        Asnotm1(Asnotm1<Amin)=Amin;
        Asnotm1(Asnotm1>Amax)=Amax;
        %%%%%%%%%%%%%
        if Pr_sno > 0.0
            %%%% Fresh snow falling 
            row = 1000; % water density [kg/m^3]
            Ta_day=mean(Ta_day); 
            ros_n = 1000*(0.05 + (((9/5)*Ta_day +32)> 0).*(((9/5)*Ta_day+32)./100).^2); %% [kg/m^3] new snow density [from Bras 1990]
            dsn = 0.001*Pr_sno_day*row/ros_n ; %%%[m] the new cumulative snowfall depth since the beginning of the snowfall event
            
            ros_nd = ros_n; %% density of dry snow  --- [kg/m^3] assumed to be the density of the new snow 
            
            % d_snowfall: snow particle diameter of new snowfall (m) Andersen 1976
            if ros_nd >= 400
                d_snowfall = 2.976e-3;
            else
                d_snowfall = 1.6e-4 + 1.1e-13*ros_nd.^4;
            end
            
            alb_fs = 0.61 + 0.21*exp(-d_snowfall/(0.001069));  % parameterization, fresh snow albedo
            K_alb = (2.4*exp(-d_snowfall/(0.000116)) + 0.33)*(Pr_sno_day^(-0.91*exp(-d_snowfall/(0.000106))-0.03)); % parameterization
            K_alb(Pr_sno_day==0)=0; 
            alb_0 = alb_fs*(1- exp(-K_alb*(Pr_sno_day))) +  Asnotm1*exp(-2*K_alb*(Pr_sno_day)); % parameterization
            
            % Sleet impact
            frac_s_Prec = Pr_sno/(Pr_liq+Pr_sno);    % solid fraction of total precipitation
            
            if frac_s_Prec > 0.5
                alb_s = (alb_0 - Asnotm1)*(frac_s_Prec -0.5)/0.5 + Asnotm1;
            else
                alb_s = Asnotm1;
            end
            
            
            % Shallow snow impact
            frac_snow = dsn/(dsn + 0.02);  %  snow cover fraction of the new snowfall [-]
            Asno = alb_s*frac_snow + Asnotm1*(1-frac_snow);
            
            
        else
            % Snow Aging - No Fresh snow 
            
            %%%% Baker et al. [1990] and Verseghy [1991]
            %Amin =0.5; %% Amin =0.2;
            tau_A = 0.008; tau_1 = 86400; %[s]
            tau_F = 0.24; %
            %tau_F = 0.48;
            %Amax = 0.85;
            if Ts < -0.01 %  A  % albedo of snow
                %%% frozen
                Asno = Asnotm1 -tau_A*dt/tau_1;
            else
                % Melting Season
                Asno= (Asnotm1 - Amin)*exp(-tau_F*dt/tau_1) + Amin; %%
            end
            Asno(Asno<Amin)=Amin;
        end
        

        % Cloud (Petzold, 1977)
        if N > 0
            dac =  0.00449 + 0.097*N.^3; % Correction for Cloud cover
        else
            dac = 0;
        end
        
        % Sun-angle (Petzold, 1977)
        if (h_S/pi*180 <=40) &&  (h_S>0)
            dab =  -0.019 + 0.248*exp(-(h_S/pi*180)/15.5); % Correction for Solar Elevation
        else
            dab=0;
        end
        
        tau_sno=Asno;
        Asno = Asno + dac  + dab ;
        
        Asno(Asno<Amin)=Amin;
        Asno(Asno>Amax)=Amax;
        
        snow_alb.dir_vis = Asno;
        snow_alb.dir_nir = Asno;
        snow_alb.dif_vis = Asno;
        snow_alb.dif_nir = Asno;
        
end
return