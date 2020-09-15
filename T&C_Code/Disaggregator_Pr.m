%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction   Rainfall Disaggregation and Re-Interpolation  %%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[Prs]=Disaggregator_Pr(Pr,pow,a)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% REFERENCES: Molnar and Burlando 2005; Onof et al., 2005 ; 
%Molnar and Burlando 2008 ; Rupp et al 2009 
%Sivakumar and Sharma 2008
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Disaggregation from 1 hour to 5 minutes 
%%% RANDOM CASCADE GENERATOR 
%%% Microcanonical Model 
%DTi=[1/2 1/4 1/8 1/16]; 
%a = a0*DTi.^gam1 ; 
%pow = pow0 + k2*log(DTi); 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%% INPUT 
%%% Pr  precipitation [mm]
%%% dt Recording Intervals [min]
b=2;%% branching number [#]
n=4;%% 60->30->15->7.5->3.75;%%Levels of subdivisions
%%% OUTPUT 
%%% Prs disaggregated precipitation [mm] 
%%%%%%%%%%%%%%%%%%%%%%%%%
if length(pow)==1
    pow=pow*ones(1,n);
    a=a*ones(1,n);
end
%%%%%%%%%%%%%%%%%%%%%%%
for i=1:n
    Prs=zeros(1,b*length(Pr));
    if mod(Prs,2)>0
        Prs=[Prs , 0];
    end
    m2=length(Pr);
    Pr0 = [Pr; Pr];
    Pr0=reshape(Pr0,1,2*m2);
    %dt=dt/2;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %%% MICROCANONICAL MODEL
    R1=unifrnd(0,1,1,length(Prs)/2);
    R2=unifrnd(0,1,1,length(Prs)/2);
    R3=betarnd(a(i),a(i),1,length(Prs)/2);
    W1=zeros(1,length(Prs)/2);
    W1(R1< pow(i) & R2 >0.5) = 1;
    W1(R1< pow(i) & R2 <=0.5) = 0;
    W1(R1 >= pow(i))= R3(R1 >= pow(i));
    W2=1-W1;
    W=[W1 ; W2];
    W=reshape(W,1,length(Prs));
    Prs=Pr0.*W;
    Pr=Prs;
end 
dt = 3.75 ; % new intervals  [min]
dtN = 5 ;% final target intervals [min]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
TL=length(Pr)*dt;%%% Total Length [min]  
%%%%%%%%%%%%%%%%%%
Prs=zeros(1,TL/dtN); %% Interpolated Rainfall 
%%%%%%%%%%%%%%
fr=4; 
n=length(Pr);
m=floor(n/fr);
Pr1=reshape(Pr(1:m*fr),fr,m);
%%%%%%%% [3X4] Row = new steps Column old steps -- sum in column == 1 
A=[ 1 1.25/3.75 0  0 ; 
    0 2.5/3.75 2.5/3.75 0 ; 
    0 0 1.25/3.75 1 ]; 
Pr2=A*Pr1; 
Prs=reshape(Pr2,1,length(Prs));
%%%%%%%%%%%%%%%%%
end