function [ShF] = Shadow_Effect(DTM,h_S,zeta_S,HZ,Z)
%%% INPUT
%%% DTM : matrix digital elevation model
%%% h_S solar height [rad]
%%% zeta_S solar azimuth [rad]
%%% HZ Horizon angle  [angular degree] from N
%%% Z Azimuth  [angular degree] from N
%%%% OUTPUT
%%% ShF Shadow Effect 1 --> no shadow 0 --> shadow
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%% COMPUTATION ShF
[m,n]=size(DTM);
S_zen=pi/2-h_S; %%%[rad]
%%%%%%%%%%%%%%%
HZr=HZ*pi/180; %%[rad]
Z= Z*pi/180; %%% [rad]
Z=[Z 2*pi]; %%[rad]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
vec_m=(abs(Z-zeta_S));
[vm,p1]=min(vec_m);  vec_m(p1)=Inf;
[vm,p2]=min(vec_m);
p3=p1; if p3==length(Z) ; p3=1; end
p4=p2; if p4==length(Z) ; p4=1; end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
ShF=zeros(m,n);
%%%%%%%%%%%%%%%%%%%
for i=1:m
    for j=1:n
        if not(isnan(DTM(i,j)))
            HZz1=HZr(i,j,p3); HZz2=HZr(i,j,p4);
            HZi= HZz1+(HZz2-HZz1)*((zeta_S - Z(p1))/(Z(p2)-Z(p1)));
            %%%%%%%%%%%%%%%%%%%%
            if HZi > S_zen
                ShF(i,j)=1;  %% No Shadow
            else
                ShF(i,j)=0;  %% Shadow
            end
        else
            ShF(i,j)=NaN;
        end
    end
end
%%%%%
return