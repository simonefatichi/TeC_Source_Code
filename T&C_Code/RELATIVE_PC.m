%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  RELATIVE PHOTOSYNTETIC CAPACITY %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[e_rel]= RELATIVE_PC(AgeL,dflo,NBL_Im,BLeaf,age_cr,aSE,L_day,Lmax_day,jDay)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
OPT_PC=6;
switch OPT_PC
    case 0
        e_rel=1;
    case 1
        e_rel=1;
    case 2
        e_rel=1;
    case 3
        e_rel=1;
    case 5
        r_age = AgeL/age_cr;
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        if r_age < 0.005
            e_rel = r_age/0.005;
        end
        if r_age >= 0.005 && r_age < 0.5
            e_rel = 1;
        end
        if r_age >= 0.5 && r_age < 1
            e_rel = -1.4*r_age + 1.7;
        end
        if r_age >= 1
            e_rel =0.3;
        end
    case 6
        if  aSE == 1  % || aSE == 2 || aSE == 0
            %%%% Bauerle et al 2012 PNAS
            e_rel = (L_day./(Lmax_day)).^2;
            e_rel(e_rel>1)=1;
        else
            e_rel=1;
        end
        %%%%%
        if  aSE == 3
            fNL = 30*NBL_Im/BLeaf; %% [1/month] 
            AgeLm = AgeL/30; %%[month]
            e_rel = 1.6104 + -0.0601*(AgeLm) + -1.2007*(fNL); % GM
            %e_rel = 1.4876 + -0.0505*(AgeLm) + -1.0806*(fNL);
            %e_rel = 0.3158 + (NBL_Im^-0.2381)  -0.001494*(AgeL);
            e_rel(e_rel>1)=1; 
            e_rel(e_rel<0.4)=0.4;
        end
        %%%%
        if  aSE == 5
            r_age2 = dflo/age_cr;
            %%%
            if r_age2<1
                e_rel=1;
            else
                e_rel= r_age2.^-8;  e_rel(e_rel>1)=1;
            end
        end 
    case 7
        %%%% Bauerle et al 2012 PNAS -- adjusted of 15 min
        e_rel = (L_day./(Lmax_day-0.25)).^2;
        e_rel(e_rel>1)=1;
    case 8
        %%%% Wu et al 2015 Science
        if  aSE == 3
            if (jDay<=121.6)
                e_rel =0.955;
            elseif ((jDay>121.6)&&(jDay<=212.8))
                e_rel =-0.0567*(jDay/30.4)+1.1817;
            elseif ((jDay>212.8)&&(jDay<=243.2))
                e_rel =0.7848;
            elseif((jDay>243.2)&&(jDay<=334.5))
                e_rel =0.0717*(jDay/30.4)+0.2117;
            else
                e_rel =1;
            end
        end
end
e_rel(e_rel<0)=0;
end