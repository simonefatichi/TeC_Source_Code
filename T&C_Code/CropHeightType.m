%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Subfunction Crop Height and type     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [hc,SAI,B,Ccrown,Nreserve,Preserve,Kreserve,AgrHarNut] = CropHeightType(LAI,LAIdead,cc,B,Ccrown,...
    Nreserve,Preserve,Kreserve,Stoich,ManI,Mpar,Veg_Param_Dyn,OPT_SoilBiogeochemistry)
%%%%
if numel(size(B)) == 3
    BRoot=B(:,cc,3);
else
    BRoot=B(cc,3);
end
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
AgrHarNut=[0 0 0];
if  ManI(cc)>0 %% Sowing

    Ccrown(cc)=Mpar(cc).Crop_crown(ManI(cc));

    if Ccrown(cc)==0
        if numel(size(B)) == 3
            B(:,cc,:)=zeros(1,8);
        else
            B(cc,:)=zeros(1,8);
        end
    end

    %%% Adding enough Nutrients in the seeds to flush leaves 
    if OPT_SoilBiogeochemistry==1
        if numel(size(B)) == 3
            Nreserve(1,cc)= B(:,cc,3)*(1/Stoich.Nl-1/Stoich.Nr)  +  B(:,cc,4)*(1/Stoich.Nl-1/Stoich.Nh);
            Preserve(1,cc)= B(:,cc,3)*(1/Stoich.Phol-1/Stoich.Phor)  +  B(:,cc,4)*(1/Stoich.Phol-1/Stoich.Phoh);
            Kreserve(1,cc)= B(:,cc,3)*(1/Stoich.Kpotl-1/Stoich.Kpotr)  +  B(:,cc,4)*(1/Stoich.Kpotl-1/Stoich.Kpoth);
        else
            Nreserve(1,cc)= B(cc,3)*(1/Stoich.Nl-1/Stoich.Nr)  +  B(cc,4)*(1/Stoich.Nl-1/Stoich.Nh);
            Preserve(1,cc)= B(cc,3)*(1/Stoich.Phol-1/Stoich.Phor)  +  B(cc,4)*(1/Stoich.Phol-1/Stoich.Phoh);
            Kreserve(1,cc)= B(cc,3)*(1/Stoich.Kpotl-1/Stoich.Kpotr)  +  B(cc,4)*(1/Stoich.Kpotl-1/Stoich.Kpoth);
        end
        AgrHarNut=-[Nreserve(1,cc) Preserve(1,cc) Kreserve(1,cc)];
    end

elseif ManI(cc) == -3 %%%%% One day after Harvest
    if OPT_SoilBiogeochemistry==1
        AgrHarNut=[Nreserve(1,cc) Preserve(1,cc) Kreserve(1,cc)];
        Nreserve(1,cc)=0;
        Preserve(1,cc)=0;
        Kreserve(1,cc)=0;
    end
    %%% Removing the crop
    Ccrown(cc)=0;
end
%%%%%%%%%%%
%%%% All these part as been outsourced to the Root_Depth_Dynamics Function
%%%%% Potential update for rooting depth --
% if  CASE_ROOT~= 1
%     disp('IN ORDER TO HAVE A VARIABLE ROOT DEPTH CASE_ROOT MUST BE 1')
%     return
% end
% % %ZR=0.5*(2*BRoot)^r;
% % ZR = 1283*BRoot^0.1713-1914; ZR(ZR<5)=5; %%[mm]
% % ZR95(cc)=ZR;
% % % %%%%
% ZR50= NaN*ZR95;
% ZRmax= NaN*ZR95;
% % % %%%% Root depth Update
% [Rf]=Root_Fraction_General(Zs,CASE_ROOT,ZR95,ZR50,0*ZR95,0*ZR50,ZRmax,0*ZRmax);
return
%%%%%%

