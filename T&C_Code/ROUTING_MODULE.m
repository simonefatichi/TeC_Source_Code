%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  ROUTING_MODULE             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[q_runon,q_channel_out,Qi_in,Slo_pot,Q_exit,Qsub_exit,T_pot,...
    QpointH,QpointC,UpointH,UpointC,Q_ch_added,Utot_H,Utot_C]= ROUTING_MODULE(dt,dth,Rd,Rh,Qi_out,q_channel_in,...
    cellsize,Area,DTM,NMAN_H,NMAN_C,MRough,WC,SN,T_flow,T_potI,Slo_top,ms_max,POT,ZWT,OPT_HEAD,Xout,Yout)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%% OUTPUT
%%% q_runon  [mm]  %% Runon
%%% Qi_in, [mm] %% Subsurface Lateral flow
%%% Slo_pot [fraction] %% Slope of Hydraulic head
%%% Q_exit [mm] %% discharge from domain surface
%%% Qsub_exit [mm] %% discharge from domain subsurface
%%% q_channel_out [mm] %%% Water in channels
%%% T_pot -- new flow direction cell of matrixs for subsurface
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% INPUT
%%% dt [s] time step
%%% Rd Dunne Runoff  [mm]
%%% Rh Horton Runoff [mm]
%%% Qi_out  Lateral Subsurface [mm/h]
%%% q_channel_in [mm] Water in channels in
%%%% cellsize [m]
%%% Area [m^2] watershed area
%%% DTM
%%% NMAN_H  [s/(m^1/3)] Manning Coefficient hillslope
%%% NMAN_C  [s/(m^1/3)] Manning Coefficient channels
%%% WC [m] Channel width
%%% SN [-] stream network identifier
%%% T_flow [ Sparse mn x mn]  Flow matrix surface
%%% T_potI  ms cells [ Sparse mn x mn]  Flow matrixs for subsurface
%%% Slo_top [ fraction] Topographic Slope
%%% ms_max Soil layers
%%% POT Head in a cell [mm]
%%% Zwt [mm] water table depth %%%
%%% OPT_HEAD options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[m_cell,n_cell]=size(DTM);
Qsur=Rh+Rd; %% [mm] %% Hillslope surface water
Qsur(Qsur<0)=0; %%%% Numerical Instability Issue
Qi_out = Qi_out*dth;  %% Subsurface flow [mm]
%%%%%%%%%%%
Q_exit = 0; %% Surface Flow exits the domain
Qsub_exit = 0; % Subsurface Flow exits the domain
npoint = length(Xout);
QpointH = zeros(npoint,1);
QpointC = zeros(npoint,1);
UpointH = zeros(npoint,1);
UpointC = zeros(npoint,1);
Qi_outR=zeros(m_cell,n_cell,ms_max);
Qi_seep=zeros(m_cell,n_cell,ms_max);
Qi_fall=zeros(m_cell,n_cell);
tH_store = zeros(m_cell,n_cell);
tC_store = zeros(m_cell,n_cell);
%Qi_in=zeros(m_cell,n_cell,ms_max);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% SUBSURFACE ROUTING
for jk=1:ms_max
    %%%% No Seepage
    Qi_seep(:,:,jk)=Qi_out(:,:,jk).*(SN); %%[mm]
    Qi_out(:,:,jk)= Qi_out(:,:,jk).*(1-SN); %% [mm]
    %%% SUBSURFACE ROUTING
    Qi_outR(:,:,jk)=Flow_Routing_Step2(DTM,T_potI{jk},Qi_out(:,:,jk)); %% [mm]
    %%%%%%%%%%%%%%%%%%%%%%%%%%
    Qsub_exit = Qsub_exit + (sum(sum(Qi_out(:,:,jk)))- sum(sum(Qi_outR(:,:,jk)))); %%%  %%[mm]
    Qi_out(:,:,jk)=Qi_outR(:,:,jk);   %% [mm]
end
Qi_in = Qi_out ; %%% [mm] Lateral subsurface for next step
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Qi_seep=sum(Qi_seep,3); %% [mm] Seepage flow from soils to channels
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% OVERLAND FLOW ROUTING
dti= 60; %%[s] Internal Time step for Surface Routing
cdti = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(sum(Qsur))>0
   while cdti < dt
      % for jj=1:1:dt/dti;
        %%%%%% SURFACE VELOCITY %%%%%%%%%%%
        Y = (Qsur/1000); %% [m]
        %%%%% Ponding - Microroughness
        Y = Y - MRough ; Y(Y<0)=0;
        %%%%
        t= (cellsize.*NMAN_H)./(Y.^(2/3).*sin(atan(Slo_top)).^0.5); %%% [s]
        t(isnan(t))=0;
        tH_store = tH_store + t*dti/dt;
        kdt= dti./t; %[]
        kdt(kdt>1)=1;
        kdt(kdt<0)=1; %% For numerical instabilities 
        kdt(isnan(kdt))=0;
        %%%% Surface Routing
        [QsurM]=Flow_Routing_Step2(DTM,T_flow,kdt.*Qsur); %%[mm]
        QsurM2 = QsurM.*(1-SN); %%[mm]
        Qi_fall = Qi_fall + QsurM.*(SN); %% [mm]
        I_fall = sum(sum(QsurM.*(SN)));
        QsurM = QsurM2;
        QsurR = QsurM + (Qsur - kdt.*Qsur) ; %% [mm]
        %%%%%%%%%%%%%%
        Q_exit= Q_exit + (sum(sum(Qsur))- sum(sum(QsurR))- I_fall); %% [mm]
        for ipo=1:npoint
            QpointH(ipo)= QpointH(ipo) + kdt(Yout(ipo),Xout(ipo))*Qsur(Yout(ipo),Xout(ipo)); %%%%%[mm]
        end
        %%%
        Qsur=QsurR; %%[mm]
        cdti = cdti +dti;
        dti =  min(min(t(t>0)));
        if cdti+dti>dt; dti = dt-cdti; end
    end
else
    Qsur = zeros(m_cell,n_cell);
end
%%%%%%%%%%
q_runon = Qsur;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% CHANNEL FLOW ROUTING
q_channel = q_channel_in + Qi_seep + Qi_fall; %%% [mm]
Q_ch_added = Qi_seep + Qi_fall; %% [mm]
dti= 2; %%[s] Internal Time step for Surface Routing
cdti = 0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(sum(q_channel))>0
   % for jj=1:1:dt/dti;
    while cdti < dt
        %%%%%% SURFACE VELOCITY %%%%%%%%%%%
        Y = (q_channel/1000).*cellsize./WC ; %% [m]
        t= (cellsize.*NMAN_C)./(Y.^(2/3).*sin(atan(Slo_top)).^0.5); %%% [s]
        t(isnan(t))=0;
        tC_store = tC_store + t*dti/dt;
        kdt= dti./t; %[]
        kdt(kdt>1)=1; 
        kdt(kdt<0)=1; %% For numerical instabilities 
        kdt(isnan(kdt))=0;
        %%%% Surface Routing
        [QchM]=Flow_Routing_Step2(DTM,T_flow,kdt.*q_channel); %%[mm]
        QchR = QchM + (q_channel - kdt.*q_channel) ; %% [mm]
        %%%%%%%%%%%%%%
        Q_exit= Q_exit + (sum(sum(q_channel))- sum(sum(QchR))); %% [mm]
        for ipo=1:npoint
            QpointC(ipo)= QpointC(ipo) + kdt(Yout(ipo),Xout(ipo))*q_channel(Yout(ipo),Xout(ipo)); %%%%%[mm]
        end
        %%%
        q_channel=QchR; %%[mm]
       cdti = cdti +dti;
        dti =  min(min(t(t>0)));
        if cdti+dti>dt; dti = dt-cdti; end
    end
else
    q_channel = zeros(m_cell,n_cell);
end
q_channel_out = q_channel;%[mm]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
for ipo=1:npoint
    UpointH(ipo)= cellsize./tH_store(Yout(ipo),Xout(ipo)); %%% Surface Velocity [m/s]
    UpointC(ipo)= cellsize./tC_store(Yout(ipo),Xout(ipo)); %%% Surface Velocity [m/s]
end
Utot_H= cellsize./tH_store; %%% Surface Velocity [m/s]
Utot_C= cellsize./tC_store; %%% Surface Velocity [m/s]
%%%%%%%%%%%%%%%%%%%%%%%%
Q_exit = Q_exit*(cellsize^2)/Area; %%[mm]
Qsub_exit = Qsub_exit*(cellsize^2)/Area; %%[mm]
%%%%%%%%%%%%%%%%%%%%%%%%
%QpointH = (QpointH/dth)*(cellsize^2)/(3600000); %% [m^3/s]
%QpointC = (QpointC/dth)*(cellsize^2)/(3600000); %% [m^3/s]
%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% HYDRAULIC HEAD PART %%%
T_pot=cell(1,ms_max);
Slo_pot = zeros(m_cell,n_cell,ms_max);
met_fdir=1;
if OPT_HEAD == 1
    %%%% NEW FLOW DIRECTIONS %%%%
    %%%% Estimation of Energy Slopes
    for jk=1:ms_max
        H = DTM + 0.001*reshape(POT(:,jk),m_cell,n_cell).*cos(atan(Slo_top)); %%% Hydraulic Head [m]
        %H = 0.001*reshape(POT(:,jk),m_cell,n_cell).*cos(atan(Slo_top)); %%% Water Potential Head [m]
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if met_fdir == 1
            %%% D-Inf flow direction method
            [R_pot, Slo_pot(:,:,jk)] = dem_flow(H,cellsize,cellsize);
            T_pot{jk} = flow_matrix(H,R_pot,cellsize,cellsize); %% Flow Matrix %%%%
        else
            Slo_pot(:,:,jk)=Slope_Aspect_indexes(H,cellsize,'mste');
            x=0:cellsize:(0+cellsize*(n_cell-1));
            y=0:cellsize:(0+cellsize*(m_cell-1));
            [x,y]=meshgrid(x,y);
            [Mpot] = flowdir(x,y,H,'type','multi'); %% Multiple D-Inf Quinn et al., 1993
            %[Mpot] = flowdir(x,y,H,'type','single'); %% D8  O'Callaghan & Mark, 1984
            T_pot{jk}=speye(m_cell*n_cell,m_cell*n_cell)-Mpot';
        end
    end
    Slo_pot(Slo_pot<0)=0;
else
    for jk=1:ms_max;
        T_pot{jk}= T_flow;
        Slo_pot(:,:,jk)=Slo_top;  %%%
    end
end
return









