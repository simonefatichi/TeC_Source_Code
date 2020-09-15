%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%% VOLUME CORRECTION FOR NEGATIVE VALUES  %%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function[V,T_H,T_L,EG,Lk]=Volume_Correction(V,EvL_Zs,RfH_Zs,RfL_Zs,EG,T_H,T_L,Lk)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
lay_H = sum(RfH_Zs'>0); %%%  number  of active layer for Transp H
lay_L = sum(RfL_Zs'>0); %%% number  of active layer for  Transp L
lay_G = sum(EvL_Zs >0); %%% number  of active layer for Evaporation
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% Compensatory Mechanism on deeper layers
for i=1:lay_G
    if (V(i)< 0) && (i < lay_G)
        V(i+1)= V(i+1) + V(i); V(i)=0;
    end
end
for i=1:max(lay_L)
    if (V(i)< 0) && (i < max(lay_L))
        V(i+1)= V(i+1) + V(i); V(i)=0;
    end
end
for i=1:max(lay_H)
    if (V(i)< 0) && (i < max(lay_H))
        V(i+1)= V(i+1) + V(i); V(i)=0;
    end
end
%%% Compensatory Mechanism on shallow layers
if sum(V < 0) > 0
    for i=max(find(V<0)):-1:2
        if (V(i)<0) &&  (V(i-1) >= 0)
            V(i-1) = V(i-1) + V(i); V(i) = 0;
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%% Brutal Correction %%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if sum(V < 0) > 0
    for i=1:lay_G
        if V(i) < 0
            if abs(V(i)) > EG
                V(i) = V(i)+ EG; EG = 0;
            else
                EG = EG + V(i); V(i) = 0;
            end
        end
    end
    for i=1:lay_L(1)
        if V(i) < 0
            if abs(V(i)) > sum(T_L)
                V(i) = V(i) + sum(T_L);  T_L=T_L*0;
            else
                T_L= T_L + V(i)*T_L/sum(T_L); V(i)= 0;
            end
        end
    end
    for i=1:lay_H(1)
        if V(i) < 0
            if abs(V(i)) > sum(T_H)
                V(i) = V(i) + sum(T_H);  T_H=T_H*0;
            else
                T_H= T_H + V(i)*T_H/sum(T_H); V(i)= 0;
            end
        end
    end
    if (V(end) < 0)
        for i=length(V):-1:2
            if (V(i)<0)
                Lk = Lk + V(i);
                V(i)=0;
            end
        end
    end
end
%%%%%%% FINAL CORRECTION 
if sum(V < 0) > 0
    %disp('WARNING: UNCONTROLLED NEGATIVE VOLUME') 
    EG = EG + sum(V.*(V<0)); V(V<0)=0; 
end
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%