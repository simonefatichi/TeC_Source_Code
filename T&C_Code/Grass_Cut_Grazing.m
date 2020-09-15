%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction  GRASS CUT AND GRAZING %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[B,LAI,LAIdead,RBi]= Grass_Cut_Grazing(B,LAI,LAIdead,dtd,aSE,Sl,mSl,jDay,jDay_cut,LAI_cut)  
%%%% INPUT
Btm1=B; 
%%%%%%%%%% Grass Cut and Grazing
if (aSE==2)  &&  not(isempty(intersect(jDay,jDay_cut)))
   if LAI_cut < 0 %%% Grazing 
        B(1)=B(1)+ LAI_cut*dtd; %%% Grazing [gC/ m^2]
        B(B<0)=0;
        if mSl == 0
            LAI = Sl*B(1);
        else
            LAI = Sl*((exp(mSl*B(1))-1)/mSl);
        end
    else
        if LAI>LAI_cut
            LAI=LAI_cut;
            if mSl == 0
                B(1)= LAI/Sl;
            else
                B(1)= (log((mSl*LAI + Sl))-log(Sl))/mSl;
            end
        end
        if LAIdead >LAI_cut
            LAIdead = LAI_cut;
            if mSl == 0
                B(7)= LAIdead/Sl;
            else
                B(7)= (log((mSl*LAIdead + Sl))-log(Sl))/mSl ;
            end
        end
    end
else
    %%% No cut 
    B=B; 
    LAI=LAI;
    LAIdead=LAIdead; 
end 
RBi(1)=(Btm1(1)-B(1))/dtd; %%  [gC/ m^2 day]
RBi(2)=(Btm1(7)-B(7))/dtd; %%  [gC/ m^2 day]
return
      