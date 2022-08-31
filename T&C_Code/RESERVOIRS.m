%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  RESERVOIR_MODULE           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[WAT_out,Q_channel_out,q_runon_out,Q_out_Res,H_Res,VOL_Res]= RESERVOIRS(DTM,cellsize,t,dth,dt,SN,WAT,Q_channel,q_runon,...
    RES_ID_List,RES_ID,RES_Outlet,Res_prop,Res_TS)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% INPUT
%%RES_ID       MAP with the ID of reservoirs
%RES_ID_List  Vector with the list of reservoirs
%RES_Outlet   Outlet coordinate of each reservoir [y x ; ... ]
%Res_TS       Time series or single values for  Res_TS(i).U  %% uptake / withdrawals [m3/s] Res_TS(i).R  %% release [m3/s]  or  Res_TS(i).y  %%% Target Level [m]
%Res_prop     h [m],V [m3], Q [m3] curve of the reservoir:  Res_prop(i).h  Res_prop(i).V Res_prop(i).Q // empty reservoir h=0;
%%% dth [h] time step
%%% dt [s] time step
%%%%%%%%%%% OUTPUT
%Q_out_Res %% [mm/h] Total Outflow from the reservoir
%H_Res %%[m]  Level in the reservoir
%VOL_Res %% [m3]  Volume of the reservoir
%%%%%%%%%%%%%%
[m_cell,n_cell]=size(DTM);
q_runon_out = reshape(q_runon*dth,m_cell,n_cell); %%[mm]
WAT_in=reshape(WAT,m_cell,n_cell); %[mm]
WAT_out = WAT_in;
Q_channel_out=Q_channel;
Q_out_Res = zeros(1,length(RES_ID_List));
H_Res = zeros(1,length(RES_ID_List));
VOL_Res = zeros(1,length(RES_ID_List));
%%%%%%%%%%
for i=1:length(RES_ID_List)
    I=RES_ID==RES_ID_List(i);
    nc=sum(sum(I));
    Ares=nc*cellsize^2; %% [m2] Area Reservoir
    %%%%
    Q_ch_to_res = sum(sum(Q_channel_out(I)))*0.001*cellsize^2;%% [m3]
    Q_channel_out(I)=0;
    %%%%
    V_res=sum(sum(WAT_in(I)*0.001*cellsize^2))+Q_ch_to_res ; %% [m3]
    %%%%
    if V_res>max(Res_prop(i).V) %%% exceeding reservoir max given level
        h_res=max(Res_prop(i).h); %% [m]
    elseif  V_res<min(Res_prop(i).V)
        h_res=0;
    else
        h_res=interp1(Res_prop(i).V,Res_prop(i).h,V_res,'linear'); %% [m]
    end
    %H_res_ori(i)=h_res; 
    %V_ori(i)=V_res; 
    %%%%%%%%%%%%%
    %%% Single value vs time series
    if length(Res_TS(i).y)==1
        y_tar=Res_TS(i).y;
        Upt = Res_TS(i).U;
        Rel = Res_TS(i).R;
    else
        y_tar=Res_TS(i).y(t);
        Upt = Res_TS(i).U(t);
        Rel = Res_TS(i).R(t);
    end
    %%%%%
    if isnan(y_tar)
        %%% Uptake and Release prescribed option
        V_res =  V_res - Upt*dt - Rel*dt ; %% [m3]
        %%%%%%%%%%
        if V_res>max(Res_prop(i).V) %%% exceeding reservoir max given level
            h_res=max(Res_prop(i).h); %% [m]
        elseif  V_res<min(Res_prop(i).V)
            h_res=0;
        else
            h_res=interp1(Res_prop(i).V,Res_prop(i).h,V_res,'linear'); %% [m]
        end
        %%%%%%%%%
        Rel2=Rel; 
    else
        %%% Target level option
        if h_res>y_tar
            V_res2=interp1(Res_prop(i).h,Res_prop(i).V,y_tar,'linear');
            Rel2 = (V_res-V_res2)/dt; %% [m3/s]
            Q_upt_cap=interp1(Res_prop(i).h,Res_prop(i).Q,y_tar,'linear'); %% [m3/s]
            if Rel2>max(Q_upt_cap) %% Release exceed uptake capacity at target level 
                Rel2=Q_upt_cap;
                V_res2= V_res - Rel2*dt;
            end
            V_res=V_res2;
            %%%%
            if V_res>max(Res_prop(i).V) %%% exceeding reservoir max given level
                h_res=max(Res_prop(i).h); %% [m]
            elseif  V_res<min(Res_prop(i).V)
                h_res=0;
            else
                h_res=interp1(Res_prop(i).V,Res_prop(i).h,V_res,'linear'); %% [m]
            end
        else
            Rel2=0;
        end
    end

    %%%% Spillway control
    if h_res>Res_prop(i).reg_lev
        if V_res>max(Res_prop(i).V) %%% exceeding reservoir max given level
            Q_out=max(Res_prop(i).Q); %% [m3/s]
        else
            Q_out=interp1(Res_prop(i).h,Res_prop(i).Q,h_res,'linear'); %% [m3/s]
        end
        %%%%%
        Q_out=Q_out-Rel2; Q_out(Q_out<0)=0;
        if Q_out*dt<V_res
            V_res =  V_res - Q_out*dt  ; %% [m3]
        else
            V_res=0; Q_out=V_res/dt;
        end
        %%%
        Q_out=Q_out+Rel2;
        %%%%
    else
        Q_out=Rel2; 
    end
    %%%%%%
    if V_res>max(Res_prop(i).V) %%% exceeding reservoir max given level
        h_res=max(Res_prop(i).h); %% [m]
    elseif  V_res<min(Res_prop(i).V)
        h_res=0;
    else
        h_res=interp1(Res_prop(i).V,Res_prop(i).h,V_res,'linear'); %% [m]
    end

    %%%%%%%%%%% Final h_res V_res and Upt Rel Q_out

    Q_to_out = (Q_out)*dt/(cellsize^2)*1000; %% [mm]

    WAT_out(I)=V_res/Ares*1000; %% [mm] water in the cell


    k=sub2ind([m_cell,n_cell],RES_Outlet(i,1),RES_Outlet(i,2));
    if SN(k)==1
        Q_channel_out(k)= Q_channel_out(k) + Q_to_out ;
    else
        q_runon_out(k)= q_runon_out(k)+Q_to_out;
    end

    Q_out_Res(i) = Q_to_out/dth; %% [mm/h]
    H_Res(i)  = h_res; %%[m] 
    VOL_Res(i) = V_res; %% [m3] 
    %drun = reshape(q_runon_out,numel(DTM),1)-q_runon;
    %Ck_res(i) = sum(sum(WAT_in))-sum(sum(WAT_out)) +sum(sum(Q_channel))- sum(sum(Q_channel_out)) -  sum(drun); %% Mass balance

end

WAT_out = reshape(WAT_out,numel(DTM),1); %%[mm]
q_runon_out = reshape(q_runon_out,numel(DTM),1); %%[mm]
q_runon_out= q_runon_out/dth; %%% [mm/h]
return






