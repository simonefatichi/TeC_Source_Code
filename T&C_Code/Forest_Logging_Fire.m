%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  HARVEST  TIMBERING  LOGGING   %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[B,RB,LAI,LAIdead,ManI]= Forest_Logging_Fire(B,RBtm1,dtd,Sl,mSl,aSE,LAI,LAIdead,Datam,Mpar)
%%%%%%
DD = datenum(Datam(1),Datam(2),Datam(3),Datam(4),0,0);
Date_log=Mpar.Date_log; %%% Date of logging
fract_log=Mpar.fract_log;
Date_fire=Mpar.Date_fire; %%% Date of fire
fire_eff = Mpar.fire_eff;
fract_resprout=Mpar.fract_resprout; %% fraction resrpout or surviving belowground
%%%%%
%%% Seed bank - longevity 1yr or 20yr  
FsR = 0.5;%%% Fraction of seed bank for Fine roots 
FsC = 0.5;%%% Fraction of seed bank for Reserves 
%%%% INPUT
ManI=0; %%% Management Indicator (0) Nothing (1) Fire (-1) Logging 
Btm1=B;
if fract_log > 0
    %%%%%%%%%% Logging
    if  sum(abs(DD-Date_log)<=0.49)>=1 
        ManI=-1;
        %%%%% Completely harvested pools
        B = B - fract_log*B; %% [gC/ m^2 ]
        %%%%
        if fract_resprout >0
            %%%% Resrpouting plant
            B(3) =  B(3) + fract_log*Btm1(3)*fract_resprout + FsR*fract_log*Btm1(5)*(1-fract_resprout);
            B(4) =  B(4) + fract_log*Btm1(4)*fract_resprout + FsC*fract_log*Btm1(5)*(1-fract_resprout);
        else
            %%% no-resrouting (growing from seeds)
            B(3) =  B(3) + FsR*fract_log*Btm1(5);
            B(4) =  B(4) + FsC*fract_log*Btm1(5);
        end
        %%%%%%%
        B(B<0)=0;
        %%%%
        if mSl == 0
            LAI = Sl*B(1);
            LAIdead = Sl*B(7);
        else
            LAI = Sl*((exp(mSl*B(1))-1)/mSl);
            LAIdead = Sl*((exp(mSl*B(7))-1)/mSl);
        end
    else
        %B=B;
        %LAI=LAI;
        %LAIdead=LAIdead;
    end
end
%%%% Fire -->
if fire_eff > 0
    if  sum(abs(DD-Date_fire)<=0.49)>=1 
        ManI=-5;
        %%%%% Completely harvested pools
        B = B - fire_eff*B; %% [gC/ m^2 ]
        %%%%
        if fract_resprout >0 || (aSE==2) 
            if (aSE==2)
                %%%% Resrpouting plant (Grass)
                B(3) =  B(3) + fire_eff*Btm1(3);
                B(4) =  B(4) + fire_eff*Btm1(4);
            else
                %%%% Resrpouting plant
                B(3) =  B(3) + fire_eff*Btm1(3)*fract_resprout +  FsR*fire_eff*Btm1(5)*(1-fract_resprout);
                B(4) =  B(4) + fire_eff*Btm1(4)*fract_resprout +  FsC*fire_eff*Btm1(5)*(1-fract_resprout);
            end
        else
            %%% no-resrouting (growing from seeds)
            B(3) =  B(3) + FsR*fire_eff*Btm1(5);
            B(4) =  B(4) + FsC*fire_eff*Btm1(5);
        end
        %%%%%%%
        B(B<0)=0;
        %%%%
        if mSl == 0
            LAI = Sl*B(1);
            LAIdead = Sl*B(7);
        else
            LAI = Sl*((exp(mSl*B(1))-1)/mSl);
            LAIdead = Sl*((exp(mSl*B(7))-1)/mSl);
        end
    else
        %B=B;
        %LAI=LAI;
        %LAIdead=LAIdead;
    end
end
if (fract_log > 0) || (fire_eff > 0)
    RB(1)=(Btm1(1)-B(1))/dtd; %%  [gC/ m^2 day] %%% Leaves - Grass
    RB(2)=(Btm1(2)-B(2))/dtd; %%  [gC/ m^2 day] %%% Sapwood
    RB(3)=(Btm1(3)-B(3))/dtd; %%  [gC/ m^2 day] %%% Fine roots
    RB(4)=(Btm1(4)-B(4))/dtd; %%  [gC/ m^2 day] %%% Carbohydrate Reserve
    RB(5)=(Btm1(5)-B(5))/dtd; %%  [gC/ m^2 day] %%% Fruit and Flower
    RB(6)=(Btm1(6)-B(6))/dtd; %%  [gC/ m^2 day] %%% Heartwood - Dead Sapwood
    RB(7)=(Btm1(7)-B(7))/dtd; %%  [gC/ m^2 day] %%% Leaves - Grass -- Standing Dead
else
    %%% No Logging / No Fire
    RB(1:7)=0; %%
end
%%%%
RB=RB+RBtm1;
end
