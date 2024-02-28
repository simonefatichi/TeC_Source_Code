%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%% This Function assigns parameters whenever they have not been
%%% specified in the MOD_PARAM file
%%!! Be careful as it could be dangerous to have parameters you did not decide

%%% Rainfall Disaggregation
if not(exist('a_dis','var'))
    
    a_dis = NaN ;
    pow_dis = NaN;
end

%%%% Terrain Properties
if not(exist('aR','var'))
    fpr=1;
    SvF=1; %% Sky View Factor
    SN=0; %% Stream Identifier
    Slo_top=0;  %% [fraction dy/dx]
    Slo_pot=zeros(1,ms); %% [fraction dy/dx]
    Asur = 1./cos(atan(Slo_top)); %% Real Area/Projected Area [m^2/m^2]
    Ared = 1;
    aR =1; %%% anisotropy ratio
    %Kh=Ks*aR;
    %%%%%%%%%%%%%%%%%%%%%%%%%%%
    cellsize=1; %%[m^2];
    aTop = 1000*cellsize^2./cellsize; %% [mm] Ratio betweeen Area/ContourLenght
end

%%% Soil Structural Parameters
if not(exist('Omac','var'))
    if SPAR == 3
        disp('NOT SPECIFIED SOIL STRUCTURE PARAMETERS')
        return
    end
    nVGM = 3.0*ones(1,ms); %  %% [-]
    alpVGM= 30*alpVG; %% [1/mm]
    Omac = 0.00*Osat; % 0.021*Osat;
    Ks_mac = 100*Ks_Zs; %%[mm/h]
    lVGM =0.5*ones(1,ms); %  %% [-]
end
%%% Soil VG modifications
if not(exist('bVG','var'))
    if SPAR == 1
        disp('NOT SPECIFIED CORRECTION FOR VG')
        return
    end
    s_SVG =zeros(1,ms);
    bVG = zeros(1,ms);
end
if not(exist('lVG','var'))
    lVG =0.5*ones(1,ms); %  %% [-]
end
%%%%%%% Soil_Parameters
Soil_Param.Ks_Zs=Ks_Zs;
Soil_Param.Osat=Osat;
Soil_Param.Ohy=Ohy;
Soil_Param.L=L;
Soil_Param.Pe=Pe;
Soil_Param.O33=O33;
Soil_Param.Ofc=Ofc;
Soil_Param.alpVG=alpVG;
Soil_Param.nVG=nVG;
Soil_Param.Ks_mac=Ks_mac;
Soil_Param.Omac=Omac;
Soil_Param.alpVGM=alpVGM;
Soil_Param.nVGM=nVGM;
Soil_Param.rsd=rsd;
Soil_Param.lan_dry=lan_dry ;
Soil_Param.lan_s=lan_s;
Soil_Param.cv_s=cv_s;
Soil_Param.s_SVG=s_SVG;
Soil_Param.bVG=bVG;
Soil_Param.lVG=lVG;
Soil_Param.lVGM=lVGM;
Soil_Param.Phy = Phy;
%%% ColorClass
if not(exist('Color_Class','var'))
    Color_Class = 0;
end
%%% Interception Others
if not(exist('In_max_urb','var'))
    In_max_urb = 5;
end
if not(exist('In_max_rock','var'))
    In_max_rock=2; %% [mm]
end

%%%%%%%%%% Interception Parameters
if not(exist('KcI','var'))
    Interc_Param.Sp_SN_In=5.9;
    Interc_Param.Sp_LAI_L_In=0.2;
    Interc_Param.Sp_LAI_H_In=0.2;
    Interc_Param.gcI=3.7;
    Interc_Param.KcI=0.06;
else
    
    Interc_Param.Sp_SN_In=Sp_SN_In;
    Interc_Param.Sp_LAI_L_In=Sp_LAI_L_In;
    Interc_Param.Sp_LAI_H_In=Sp_LAI_H_In;
    Interc_Param.gcI=gcI;
    Interc_Param.KcI=KcI;
end
if not(exist('Kct','var'))
    Kct=0.75; %%%
end
if not(exist('d_leaf_H','var'))
    d_leaf_H = [1]; %%[cm]
    d_leaf_L= [1];  %% [cm]
end
if not(exist('OM_H','var'))
    OM_H=1;
    OM_L=1;
end
%%%%%% Snow-Ice Parameters
if not(exist('Aice','var'))
    SnowIce_Param.TminS=-0.8;
    SnowIce_Param.TmaxS=2.8;
    SnowIce_Param.WatFreez_Th=-8;
    SnowIce_Param.dz_ice=0.54;
    SnowIce_Param.Th_Pr_sno=10;
    SnowIce_Param.ros_max1=550;
    SnowIce_Param.ros_max2=260;
    SnowIce_Param.Ice_wc_sp=0.01;
    SnowIce_Param.ros_Ice_thr=500;
    SnowIce_Param.Aice=0.35;
else
    %%%%%% Snow-Ice Parameters
    SnowIce_Param.TminS=TminS;
    SnowIce_Param.TmaxS=TmaxS;
    SnowIce_Param.WatFreez_Th=WatFreez_Th;
    SnowIce_Param.dz_ice=dz_ice;
    SnowIce_Param.Th_Pr_sno=Th_Pr_sno;
    SnowIce_Param.ros_max1=ros_max1;
    SnowIce_Param.ros_max2=ros_max2;
    SnowIce_Param.Ice_wc_sp=Ice_wc_sp;
    SnowIce_Param.ros_Ice_thr=ros_Ice_thr;
    SnowIce_Param.Aice=Aice;
end

%%%%%%%%% Debris-Cover properties
if not(exist('dbThick','var')) || dbThick == 0
    Deb_Par.alb= 0.13 ;
    Deb_Par.e_sur =  0.94;
    Deb_Par.lan = 0.94;
    Deb_Par.rho = 1496;  % [kg/m^3]
    Deb_Par.cs = 948;   % [J/kg K]
    Deb_Par.zom = 0.016 ;
    dbThick = 0;
    Zs_deb=[];
else
    if dbThick>5
        nst=10;
        k=(dbThick/5)^(1/(nst-1));
        Zs_deb = [0 5*k.^(1:nst-1)]; %% [mm]
        clear nst k
    else
        disp('Debris is not thick enough')
        return
    end
end
if not(exist('Sllit','var'))
    Sllit = 2 ; %%% Litter Specific Leaf area [m2 Litter / kg DM]
end
%%%%
if not(exist('Urb_Par','var'))
    Urb_Par.alb= 0.15;
    Urb_Par.e_sur= 0.92;
    Urb_Par.BuildH= 12;
end
%%%%%%%%%%
%%%%% Vegetation High Parameters
if not(exist('Vmax_H','var'))
    VegH_Param.KnitH=NaN*ones(1,cc);
    VegH_Param.mSl_H=NaN*ones(1,cc);
    VegH_Param.Sl_H=NaN*ones(1,cc);
    VegH_Param.CT_H=NaN*ones(1,cc);
    VegH_Param.Vmax_H=NaN*ones(1,cc);
    VegH_Param.FI_H=NaN*ones(1,cc);
    VegH_Param.a1_H=NaN*ones(1,cc);
    VegH_Param.go_H=NaN*ones(1,cc);
    VegH_Param.DSE_H=NaN*ones(1,cc);
    VegH_Param.Ha_H=NaN*ones(1,cc);
    VegH_Param.Do_H=NaN*ones(1,cc);
    VegH_Param.gmes_H=NaN*ones(1,cc);
    VegH_Param.rjv_H=NaN*ones(1,cc);
    VegH_Param.Psi_sto_50_H=NaN*ones(1,cc);
    VegH_Param.Psi_sto_00_H=NaN*ones(1,cc);
    VegH_Param.Axyl_H=NaN*ones(1,cc);
    VegH_Param.PsiL50_H=NaN*ones(1,cc);
    VegH_Param.PsiL00_H=NaN*ones(1,cc);
    VegH_Param.Kleaf_max_H=NaN*ones(1,cc);
    VegH_Param.Cl_H=NaN*ones(1,cc);
    VegH_Param.Kx_max_H=NaN*ones(1,cc);
    VegH_Param.PsiX50_H=NaN*ones(1,cc);
    VegH_Param.Cx_H=NaN*ones(1,cc);
else
    VegH_Param.KnitH=KnitH;
    VegH_Param.mSl_H=mSl_H;
    VegH_Param.Sl_H=Sl_H;
    VegH_Param.CT_H=CT_H ;
    VegH_Param.Vmax_H=Vmax_H;
    VegH_Param.FI_H=FI_H;
    VegH_Param.a1_H=a1_H;
    VegH_Param.go_H=go_H;
    VegH_Param.DSE_H=DSE_H;
    VegH_Param.Ha_H=Ha_H;
    VegH_Param.Do_H=Do_H;
    VegH_Param.gmes_H=gmes_H;
    VegH_Param.rjv_H=rjv_H;
    VegH_Param.Psi_sto_50_H=Psi_sto_50_H;
    VegH_Param.Psi_sto_00_H=Psi_sto_00_H;
    VegH_Param.Axyl_H=Axyl_H;
    VegH_Param.PsiL50_H=PsiL50_H;
    VegH_Param.PsiL00_H=PsiL00_H;
    VegH_Param.Kleaf_max_H=Kleaf_max_H;
    VegH_Param.Cl_H=Cl_H;
    VegH_Param.Kx_max_H=Kx_max_H;
    VegH_Param.PsiX50_H=PsiX50_H;
    VegH_Param.Cx_H=Cx_H;
end

if not(exist('Osm_reg_Max_H','var'))
    VegH_Param.Osm_reg_Max_H = 0*ones(1,cc);
    VegH_Param.eps_root_base_H=  0.9*ones(1,cc);
else
    VegH_Param.Osm_reg_Max_H = Osm_reg_Max_H;   
    VegH_Param.eps_root_base_H=  eps_root_base_H; 
end


%%%%% Vegetation Low Parameters
if not(exist('Vmax_L','var'))
    VegL_Param.KnitL=NaN*ones(1,cc);
    VegL_Param.mSl_L=NaN*ones(1,cc);
    VegL_Param.Sl_L=NaN*ones(1,cc);
    VegL_Param.CT_L=NaN*ones(1,cc);
    VegL_Param.Vmax_L=NaN*ones(1,cc);
    VegL_Param.FI_L=NaN*ones(1,cc);
    VegL_Param.a1_L=NaN*ones(1,cc);
    VegL_Param.go_L=NaN*ones(1,cc);
    VegL_Param.DSE_L=NaN*ones(1,cc);
    VegL_Param.Ha_L=NaN*ones(1,cc);
    VegL_Param.Do_L=NaN*ones(1,cc);
    VegL_Param.gmes_L=NaN*ones(1,cc);
    VegL_Param.rjv_L=NaN*ones(1,cc);
    VegL_Param.Psi_sto_50_L = NaN*ones(1,cc);
    VegL_Param.Psi_sto_00_L=NaN*ones(1,cc);
    VegL_Param.Axyl_L=NaN*ones(1,cc);
    VegL_Param.PsiL50_L=NaN*ones(1,cc);
    VegL_Param.PsiL00_L=NaN*ones(1,cc);
    VegL_Param.Kleaf_max_L= NaN*ones(1,cc);
    VegL_Param.Cl_L=NaN*ones(1,cc);
    VegL_Param.Kx_max_L=NaN*ones(1,cc);
    VegL_Param.PsiX50_L=NaN*ones(1,cc);
    VegL_Param.Cx_L=NaN*ones(1,cc);
else
    VegL_Param.KnitL=KnitL;
    VegL_Param.mSl_L=mSl_L;
    VegL_Param.Sl_L=Sl_L;
    VegL_Param.CT_L=CT_L;
    VegL_Param.Vmax_L=Vmax_L;
    VegL_Param.FI_L=FI_L;
    VegL_Param.a1_L=a1_L;
    VegL_Param.go_L=go_L;
    VegL_Param.DSE_L=DSE_L;
    VegL_Param.Ha_L=Ha_L;
    VegL_Param.Do_L=Do_L;
    VegL_Param.gmes_L=gmes_L;
    VegL_Param.rjv_L=rjv_L;
    VegL_Param.Psi_sto_50_L = Psi_sto_50_L;
    VegL_Param.Psi_sto_00_L=Psi_sto_00_L;
    VegL_Param.Axyl_L=Axyl_L;
    VegL_Param.PsiL50_L=PsiL50_L;
    VegL_Param.PsiL00_L=PsiL00_L;
    VegL_Param.Kleaf_max_L= Kleaf_max_L;
    VegL_Param.Cl_L=Cl_L;
    VegL_Param.Kx_max_L=Kx_max_L;
    VegL_Param.PsiX50_L=PsiX50_L;
    VegL_Param.Cx_L=Cx_L;
end
%%%%%%%%%%%%%%%%%%%%%%%%%

if not(exist('Osm_reg_Max_L','var'))
    VegL_Param.Osm_reg_Max_L = 0*ones(1,cc);
    VegL_Param.eps_root_base_L=  0.9*ones(1,cc);
else
    VegL_Param.Osm_reg_Max_L = Osm_reg_Max_L;   
    VegL_Param.eps_root_base_L=  eps_root_base_L; 
end

%%%%% Vegetation High Parameters for vegetation dynamics
if not(exist('aSE_H','var'))
    VegH_Param_Dyn.Sl = NaN*ones(1,cc);
    VegH_Param_Dyn.mSl = NaN*ones(1,cc);
    VegH_Param_Dyn.Stoich = Stoich_H;
    VegH_Param_Dyn.r =NaN*ones(1,cc);
    VegH_Param_Dyn.gR = NaN*ones(1,cc);
    VegH_Param_Dyn.LtR = NaN*ones(1,cc);
    VegH_Param_Dyn.eps_ac = NaN*ones(1,cc);
    VegH_Param_Dyn.aSE = NaN*ones(1,cc);
    VegH_Param_Dyn.Trr = NaN*ones(1,cc);
    VegH_Param_Dyn.dd_max = NaN*ones(1,cc);
    VegH_Param_Dyn.dc_C = NaN*ones(1,cc);
    VegH_Param_Dyn.Tcold = NaN*ones(1,cc);
    VegH_Param_Dyn.drn = NaN*ones(1,cc);
    VegH_Param_Dyn.dsn= NaN*ones(1,cc);
    VegH_Param_Dyn.age_cr = NaN*ones(1,cc);
    VegH_Param_Dyn.Bfac_lo = NaN*ones(1,cc);
    VegH_Param_Dyn.Bfac_ls = NaN*ones(1,cc);
    VegH_Param_Dyn.Tlo= NaN*ones(1,cc);
    VegH_Param_Dyn.Tls = NaN*ones(1,cc);
    VegH_Param_Dyn.mjDay = NaN*ones(1,cc);
    VegH_Param_Dyn.LDay_min =NaN*ones(1,cc);
    VegH_Param_Dyn.dmg= NaN*ones(1,cc);
    VegH_Param_Dyn.Mf = NaN*ones(1,cc);
    VegH_Param_Dyn.Wm =NaN*ones(1,cc);
    VegH_Param_Dyn.LAI_min =NaN*ones(1,cc);
    VegH_Param_Dyn.LDay_cr= NaN*ones(1,cc);
    VegH_Param_Dyn.PsiG50 = NaN*ones(1,cc);
    VegH_Param_Dyn.PsiG99 =NaN*ones(1,cc);
    VegH_Param_Dyn.gcoef =NaN*ones(1,cc);
    VegH_Param_Dyn.fab= NaN*ones(1,cc);
    VegH_Param_Dyn.fbe = NaN*ones(1,cc);
    VegH_Param_Dyn.Klf =NaN*ones(1,cc);
    VegH_Param_Dyn.ff_r = NaN*ones(1,cc);
    VegH_Param_Dyn.PAR_th =NaN*ones(1,cc);
    VegH_Param_Dyn.PsiL50 = NaN*ones(1,cc);
    VegH_Param_Dyn.PsiL00 = NaN*ones(1,cc);
else
    VegH_Param_Dyn.Sl = Sl_H;
    VegH_Param_Dyn.mSl = mSl_H;
    VegH_Param_Dyn.Stoich = Stoich_H;
    VegH_Param_Dyn.r = r_H;
    VegH_Param_Dyn.gR = gR_H;
    VegH_Param_Dyn.LtR = LtR_H;
    VegH_Param_Dyn.eps_ac = eps_ac_H;
    VegH_Param_Dyn.aSE = aSE_H;
    VegH_Param_Dyn.Trr = Trr_H;
    VegH_Param_Dyn.dd_max = dd_max_H;
    VegH_Param_Dyn.dc_C = dc_C_H;
    VegH_Param_Dyn.Tcold = Tcold_H;
    VegH_Param_Dyn.drn = drn_H;
    VegH_Param_Dyn.dsn= dsn_H;
    VegH_Param_Dyn.age_cr = age_cr_H;
    VegH_Param_Dyn.Bfac_lo = Bfac_lo_H;
    VegH_Param_Dyn.Bfac_ls = Bfac_ls_H;
    VegH_Param_Dyn.Tlo= Tlo_H;
    VegH_Param_Dyn.Tls = Tls_H;
    VegH_Param_Dyn.mjDay = mjDay_H;
    VegH_Param_Dyn.LDay_min =LDay_min_H;
    VegH_Param_Dyn.dmg= dmg_H;
    VegH_Param_Dyn.Mf = Mf_H;
    VegH_Param_Dyn.Wm = Wm_H;
    VegH_Param_Dyn.LAI_min =LAI_min_H;
    VegH_Param_Dyn.LDay_cr= LDay_cr_H;
    VegH_Param_Dyn.PsiG50 = PsiG50_H;
    VegH_Param_Dyn.PsiG99 =PsiG99_H;
    VegH_Param_Dyn.gcoef =gcoef_H;
    VegH_Param_Dyn.fab= fab_H;
    VegH_Param_Dyn.fbe = fbe_H;
    VegH_Param_Dyn.Klf =Klf_H;
    VegH_Param_Dyn.ff_r =ff_r_H;
    VegH_Param_Dyn.PAR_th =PAR_th_H;
    VegH_Param_Dyn.PsiL50 = PsiL50_H;
    VegH_Param_Dyn.PsiL00 = PsiL00_H;
end
if not(exist('soCrop_H','var'))
    VegH_Param_Dyn.soCrop = NaN*ones(1,cc);
    VegH_Param_Dyn.MHcrop=  NaN*ones(1,cc);
    VegH_Param_Dyn.Sl_emecrop=  NaN*ones(1,cc);
else
    VegH_Param_Dyn.soCrop = soCrop_H;  
    VegH_Param_Dyn.MHcrop=  MHcrop_H; 
    VegH_Param_Dyn.Sl_emecrop= Sl_emecrop_H;
end

if not(exist('a_root_H','var'))
    a_root_H=0.001; 
else
    a_root_H=a_root_H; 
end 

%%%%% Vegetation Low Parameters for vegetation dynamics
if not(exist('aSE_L','var'))
    VegL_Param_Dyn.Sl = NaN*ones(1,cc);
    VegL_Param_Dyn.mSl = NaN*ones(1,cc);
    VegL_Param_Dyn.Stoich = Stoich_L;
    VegL_Param_Dyn.r = NaN*ones(1,cc);
    VegL_Param_Dyn.gR = NaN*ones(1,cc);
    VegL_Param_Dyn.LtR = NaN*ones(1,cc);
    VegL_Param_Dyn.eps_ac = NaN*ones(1,cc);
    VegL_Param_Dyn.aSE = NaN*ones(1,cc);
    VegL_Param_Dyn.Trr = NaN*ones(1,cc);
    VegL_Param_Dyn.dd_max = NaN*ones(1,cc);
    VegL_Param_Dyn.dc_C = NaN*ones(1,cc);
    VegL_Param_Dyn.Tcold = NaN*ones(1,cc);
    VegL_Param_Dyn.drn = NaN*ones(1,cc);
    VegL_Param_Dyn.dsn = NaN*ones(1,cc);
    VegL_Param_Dyn.age_cr = NaN*ones(1,cc);
    VegL_Param_Dyn.Bfac_lo = NaN*ones(1,cc);
    VegL_Param_Dyn.Bfac_ls = NaN*ones(1,cc);
    VegL_Param_Dyn.Tlo= NaN*ones(1,cc);
    VegL_Param_Dyn.Tls = NaN*ones(1,cc);
    VegL_Param_Dyn.mjDay = NaN*ones(1,cc);
    VegL_Param_Dyn.LDay_min =NaN*ones(1,cc);
    VegL_Param_Dyn.dmg= NaN*ones(1,cc);
    VegL_Param_Dyn.Mf = NaN*ones(1,cc);
    VegL_Param_Dyn.Wm = NaN*ones(1,cc);
    VegL_Param_Dyn.LAI_min =NaN*ones(1,cc);
    VegL_Param_Dyn.LDay_cr= NaN*ones(1,cc);
    VegL_Param_Dyn.PsiG50 = NaN*ones(1,cc);
    VegL_Param_Dyn.PsiG99 =NaN*ones(1,cc);
    VegL_Param_Dyn.gcoef =NaN*ones(1,cc);
    VegL_Param_Dyn.fab= NaN*ones(1,cc);
    VegL_Param_Dyn.fbe = NaN*ones(1,cc);
    VegL_Param_Dyn.Klf =NaN*ones(1,cc);
    VegL_Param_Dyn.ff_r =NaN*ones(1,cc);
    VegL_Param_Dyn.PAR_th =NaN*ones(1,cc);
    VegL_Param_Dyn.PsiL50 = NaN*ones(1,cc);
    VegL_Param_Dyn.PsiL00 = NaN*ones(1,cc);
else
    VegL_Param_Dyn.Sl = Sl_L;
    VegL_Param_Dyn.mSl = mSl_L;
    VegL_Param_Dyn.Stoich = Stoich_L;
    VegL_Param_Dyn.r = r_L;
    VegL_Param_Dyn.gR = gR_L;
    VegL_Param_Dyn.LtR = LtR_L;
    VegL_Param_Dyn.eps_ac = eps_ac_L;
    VegL_Param_Dyn.aSE = aSE_L;
    VegL_Param_Dyn.Trr = Trr_L;
    VegL_Param_Dyn.dd_max = dd_max_L;
    VegL_Param_Dyn.dc_C = dc_C_L;
    VegL_Param_Dyn.Tcold = Tcold_L;
    VegL_Param_Dyn.drn = drn_L;
    VegL_Param_Dyn.dsn = dsn_L;
    VegL_Param_Dyn.age_cr = age_cr_L;
    VegL_Param_Dyn.Bfac_lo = Bfac_lo_L;
    VegL_Param_Dyn.Bfac_ls = Bfac_ls_L;
    VegL_Param_Dyn.Tlo= Tlo_L;
    VegL_Param_Dyn.Tls = Tls_L;
    VegL_Param_Dyn.mjDay = mjDay_L;
    VegL_Param_Dyn.LDay_min =LDay_min_L;
    VegL_Param_Dyn.dmg= dmg_L;
    VegL_Param_Dyn.Mf = Mf_L;
    VegL_Param_Dyn.Wm = Wm_L;
    VegL_Param_Dyn.LAI_min =LAI_min_L;
    VegL_Param_Dyn.LDay_cr= LDay_cr_L;
    VegL_Param_Dyn.PsiG50 = PsiG50_L;
    VegL_Param_Dyn.PsiG99 =PsiG99_L;
    VegL_Param_Dyn.gcoef =gcoef_L;
    VegL_Param_Dyn.fab= fab_L;
    VegL_Param_Dyn.fbe = fbe_L;
    VegL_Param_Dyn.Klf =Klf_L;
    VegL_Param_Dyn.ff_r =ff_r_L;
    VegL_Param_Dyn.PAR_th =PAR_th_L;
    VegL_Param_Dyn.PsiL50 = PsiL50_L;
    VegL_Param_Dyn.PsiL00 = PsiL00_L;
end
if not(exist('soCrop_L','var'))
    VegL_Param_Dyn.soCrop = NaN*ones(1,cc);
    VegL_Param_Dyn.MHcrop=  NaN*ones(1,cc);
    VegL_Param_Dyn.Sl_emecrop=  NaN*ones(1,cc);
else
    VegL_Param_Dyn.soCrop = soCrop_L;  
    VegL_Param_Dyn.MHcrop=  MHcrop_L; 
    VegL_Param_Dyn.Sl_emecrop= Sl_emecrop_L;
end

if not(exist('a_root_L','var'))
    a_root_L=0.001; 
else
    a_root_L=a_root_L; 
end 