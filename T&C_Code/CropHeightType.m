%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Crop Height and type     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hc,SAI,B,Ccrown,ZR95,Nreserve,Preserve,Kreserve,TdpI,Bfac_weekI,Rf] = CropHeightType(LAI,LAIdead,cc,ZR95,B,Zs,CASE_ROOT,Ccrown,...
    Nreserve,Preserve,Kreserve,TdpI,Bfac_weekI,ManI,Mpar,Veg_Param_Dyn)
%%%% 
BRoot=B(:,cc,3); 
%Crop_type=Mpar.Crop_type; 
MHcrop =  Veg_Param_Dyn.MHcrop(cc); 
%%% INPUT
LAI = LAI + LAIdead/3; 
LAI_max=6.5; 
%%% CLM5.0 / AgroIBIS 
if LAI/(LAI_max-1)  < 1.0
    hc =  MHcrop*(LAI/(LAI_max-1))^2;
else
    hc = MHcrop; 
    %%%%%%
end
hc(hc<0.05)=0.05; 
SAI=max(0.15*LAI,0.001); 
%%%%%%%%%%%%%%%%%
if  ManI(cc)>0 
    ccrop=find(Ccrown~=0); fcrop=Mpar(cc).Crop_type(ManI(cc));
    Ccrown=Ccrown*0; ZR95=ZR95*0;
    Ccrown(fcrop)=Mpar(cc).Crop_crown(fcrop);
    ZR95(fcrop)= Mpar(cc).Crop_root(fcrop);
    if fcrop ~= ccrop
        B(:,fcrop,:)=B(:,ccrop,:);  B(:,ccrop,:)=zeros(1,8);
        Nreserve(:,fcrop)=Nreserve(:,ccrop); Nreserve(:,ccrop)=0; 
        Preserve(:,fcrop)=Preserve(:,ccrop); Preserve(:,ccrop)=0; 
        Kreserve(:,fcrop)=Kreserve(:,ccrop); Kreserve(:,ccrop)=0; 
        TdpI(:,fcrop)=TdpI(:,ccrop); 
        Bfac_weekI(:,fcrop)=Bfac_weekI(:,ccrop); 
    end
    cc=fcrop; BRoot=B(:,cc,3);
end
%%%%%%%%%%%
%%%%% Potential update for rooting depth -- 
if  CASE_ROOT~= 1 
    disp('IN ORDER TO HAVE A VARIABLE ROOT DEPTH CASE_ROOT MUST BE 1')
    return
end
% %ZR=0.5*(2*BRoot)^r;
% ZR = 1283*BRoot^0.1713-1914; ZR(ZR<5)=5; %%[mm]
% ZR95(cc)=ZR;
% % %%%%
ZR50= NaN*ZR95;
ZRmax= NaN*ZR95;
% % %%%% Root depth Update
[Rf]=Root_Fraction_General(Zs,CASE_ROOT,ZR95,ZR50,0*ZR95,0*ZR50,ZRmax,0*ZRmax);
return
%%%%%%

