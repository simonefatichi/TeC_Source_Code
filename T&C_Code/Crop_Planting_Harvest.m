%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  CROP HARVESTING or PLANTING   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[B,RB,LAI,LAIdead,ManI,PHE_S,dflo,AgeL,AgeDL]= Crop_Planting_Harvest(B,RBtm1,dtd,LAI,LAIdead,PHE_S,dflo,AgeL,AgeDL,Datam,Mpar)
%%%%%%
DD = datenum(Datam(1),Datam(2),Datam(3),Datam(4),0,0);
Date_sowing=Mpar.Date_sowing+10; %%% Emergence as Date of sowing + 10 days  
Date_harvesting=Mpar.Date_harvesting;
Crop_B= Mpar.Crop_B; 
%%%%%
%%%% INPUT
ManI=0; %%% Management Indicator (-2) Planting or Harvesting 
Btm1=B;
%%%%%%%%%% Planting
if  sum(abs(DD-Date_sowing)<=0.49)>=1
    ManI=find(abs(DD-Date_sowing)<=0.49); %%% identify which sowing season 
    %%%%%
    B(3) = Crop_B(1);
    B(4) = Crop_B(2);
    %%%%
    B(B<0)=0;
    %%%%
end
%%%%%%%%%% Harvesting
if  sum(abs(DD-Date_harvesting)<=0.49)>=1
    ManI=-2;
    %%%%% Completely harvested pools // crop removed
    B(1:7)=0; %%; %% [gC/ m^2 ]
    %%%%%%%
    B(B<0)=0;
    %%%%
    LAI=0; LAIdead=0;
end
%%%%%%%%%%%
%if ManI~=0
if (sum(find(abs(DD-Date_sowing)<=0.49))>0) ||  (ManI == -2) %% planting or harvest 
    RB(1)=(Btm1(1)-B(1))/dtd; %%  [gC/ m^2 day] %%% Leaves - Grass
    RB(2)=(Btm1(2)-B(2))/dtd; %%  [gC/ m^2 day] %%% Sapwood
    RB(3)=(Btm1(3)-B(3))/dtd; %%  [gC/ m^2 day] %%% Fine roots
    RB(4)=(Btm1(4)-B(4))/dtd; %%  [gC/ m^2 day] %%% Carbohydrate Reserve
    RB(5)=(Btm1(5)-B(5))/dtd; %%  [gC/ m^2 day] %%% Fruit and Flower
    RB(6)=(Btm1(6)-B(6))/dtd; %%  [gC/ m^2 day] %%% Heartwood - Dead Sapwood
    RB(7)=(Btm1(7)-B(7))/dtd; %%  [gC/ m^2 day] %%% Leaves - Grass -- Standing Dead
else
    %%% No Management 
    RB(1:7)=0; %%
end
%%%% If crop still not out - re-initialize phenology 
if (dflo==2) && (LAI==0) 
    PHE_S=1; dflo= 0; AgeL= 0; AgeDL=0; 
end 
%%%%
RB=RB+RBtm1;
end
