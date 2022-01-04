%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  Allocation_Coefficients  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% References 
%%%% Wolf et al 2011  Ecol. Appl and GBC 
%%%% Enquist-Niklas 2002 Science  Niklas-Enquist 2002 Am. Naturalist 
%%%% ALLOCATION %%% Friedlingstein et al 1998
%%% Krinner et al 2015 
function[fs1,fr1,fl1]= Allocation_Coefficients(TBio,LAI,Bfac,Se,Ts,FNC,aSE,age_cr,dflo,soCrop,OPT_VCA)
%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%% Compute the preliminary allocation coefficients 
%%% fr allocation to root
%%% fs allocation to sapwood
%%% fl allocation to leaf
%%% All. Constraints on Allocation partions 
Tot_Biomass = max(0,TBio) ;  %25:1:500; %275; %% [MgDM / ha] %%%  [ton DM / ha ] % 
%%%%%
%%% Allocation Coefficients are weakly depedent on productivity NPP or GPP
%%% increase in fine roots with low produc. spurious correlation with water limitation 
%%% This depedence is not accounted for, only depedence on the total biomass are accounted for.  
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%5
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ro=0.3; so =0.3;
%%% To considr minor allocation to woods for young growing forest 
if aSE ~= 2 && OPT_VCA >= 1 
    %so = min(0.3,0.1002*Tot_Biomass.^0.3061);  %% Wolf et al 2011 Ecol. Appl 
    so = min(0.3,0.01775*Tot_Biomass.^0.6937); %%% Recomputed by me 
    %so= min(0.3,0.3*Tot_Biomass./12.2);
end
if (aSE == 5) 
    so=soCrop; 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ke=0.15; %% Light extinction 
AH = max(0.1,Bfac); %% Moisture Stress
Al = max(exp(-ke*(LAI)),0.1); %% Light Lack of
%%% Nutrient Release --  
ANH =1;% min(1,max(0.5,Se)); %% Nitrogen Moisture Stress
%To = 30; Tc= 10;%%[°C]
ANT = 1; % min(1,max(0.1,2^((Ts -To)/Tc))); %% Nitrogen Temp Stress
FNC= max(0.1,FNC); 
AN = ANH*ANT*FNC; %% Nitrogen Stress
Ab = min(AH,AN);
%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Parameter controlling allocation
rmin=0.15; amin=0.2; amax=0.5; % ro=0.3; so =0.3; %% From Krinner 2005
fr1 = max(rmin,ro*(3*Al/(Al+2*Ab)));
fs1 = min(0.75,so*(3*Ab./(2*Al + Ab)));
fl1= max(amin,min(amax,1-fr1-fs1));
fr1= 1 -fs1 -fl1;
%%%%%%%%%
if (aSE == 3)  %%% Tropical Forest 
    %fl1 = fl1; 
    fl1 =  1-dflo/age_cr ;
    fl1(fl1>1)=1; fl1(fl1<0.0)=0.0; %
    frem = 1-fl1;
    tmp= frem*(fs1/(fr1+fs1));
    fr1 = frem*(fr1/(fr1+fs1));
    fs1=tmp;
end 
%%%%%%%%%%%%%%%%%%%%%%%%
if (aSE == 2) 
    fsl1 = fs1*(fl1/(fl1+fr1));
    fl1 = fl1 + fsl1;
    fr1 = fr1 + (fs1-fsl1);
    fs1=0;
end
return 
