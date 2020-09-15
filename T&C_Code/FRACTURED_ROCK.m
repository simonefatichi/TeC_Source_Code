%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  FRACTURED_ROCK             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Q_channel,FROCK,Qflow_rock]= FRACTURED_ROCK(Q_channeltm1,FROCKtm1,SNn,dth,m_cell,n_cell,num_cell,Kres_Rock)
%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% Reference: 
%%%%%% OUTPUT
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% INPUT
%%% FROCK [mm] Volume in fractured rock vertical 
%%% Q_channel [mm] Volume in the channels 
%%% SNn [0/1] Stream network identifier 
%%% dth [h] time step 
%Kres_Rock [h] Aquifer constant 
%%%%% 
Q_channeltm1=reshape(Q_channeltm1,num_cell,1);
Qrock = dth*FROCKtm1/Kres_Rock; %%% [mm]   
FROCK = FROCKtm1 - Qrock; %%[mm]
Qrock_TOT = sum(Qrock); 
Q_channel = Q_channeltm1 + SNn*(Qrock_TOT/sum(SNn)); %%[mm]
Qflow_rock = Qrock_TOT;
Q_channel=reshape(Q_channel,m_cell,n_cell);
end 