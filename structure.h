/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

#include"setting.h"

/*****************************/
/* declaration of structures */
/*****************************/

/* --suffix-- */
/* f: foliage */
/* c: stem and branch */
/* r: root */
/* l: litter */
/* h: humus */

double sres_co2[111];

/* grid conditions, derived from data files ***************************************/
struct Grid{ 
	/* location ***********/
	long 	row;					/* grid order, row  in 0.5 degree grid */
	long 	col;					/* grid order, column in 0.5 degree grid */
	long 	gcm_row;				/* grid order, row in GCM's grid */
	long 	gcm_col;				/* grid order, column in GCM's grid */
	long	ncep_lat, ncep_lon;
	long 	chaser_row;				/* grid order, row in CHASER grid */
	long 	chaser_col;				/* grid order, column in CHASER grid */
	
	double 	lat;					/* latitude of the center of grid cell */
	double 	lon;					/* longitude of the center of grid cell */

	double 	area;					/* biome area for each cell */
	/* double 	cell_area;	*/		
	/* double 	land_area;	*/		/* total area and land area for each cell */

	long 	country;				/* country number */
	long 	region;					/* continental region number */
	long 	n_olson;				/* cell numbers from the origin: olson map */
	long 	n_sage;					/* cell numbers from the origin: SAGE map */
	long 	n_crop;					/* cell numbers from the origin: SAGE map */
	short	reg_g;					/* region ID by F.Giorgi */
	short	rank_nat;				/* 1: developing, 2: developed */
	
	/* Olson's actual biome ID */
	long 	veg_olson;			
	/* SAGE's potential biome ID */
	long 	veg_sage;	
	/* dominant crop type ID */
	long 	veg_crop;	
	
	/* time ***********/
	long 	phase;					/* simulation phase: 0-spinup, 1-past, 2-future */
	long 	y;						/* calculation time from the simulation onset, in year */
	long 	m;						/* month of the year, from Jan. to Dec., 0 to 11 */
    long    h;
	
	long 	simy;					/* year for simulation (AD) */
	long 	co2y;					/* year for CO2 level estimation */
	long 	climy;					/* year of climate data */
    long    lucy;                   /* year of land-use data: 2014/09/11 by A.Ito */
    long    niny;                   /* year of nitrogen input: 2015/11/19 by A.Ito */

	/* atmospheric condition ***********/
	double 	bco2[ASTEP];			/* background CO2 concentration, in ppmv */
	double 	d13c_bco2[ASTEP];		/* stable carbon isotope composition of background CO2, permille */
	double	d14c_bco2[ASTEP];		/* D14C of atmospheric CO2: added by A.Ito (2009/06/23) */
    double  bo3[ASTEP];             /* monthly O3, ppb */

	/* climate condition: *[] means the transitional value */
	double 	tmp_sfc[ASTEP];			/* ground surface temperature, degree Celcius */
	double 	tmp_2m[ASTEP];			/* 2m air temperature, degree Celcius */
	double 	tmp10_soil[ASTEP];		/* soil temperature at 10 cm depth, degree Celcius */
	double 	tmp200_soil[ASTEP];		/* soil temperature at 200 cm depth, degree Celcius */
	double 	tcdc_clm[ASTEP];		/* total cloudiness, fraction */
	double 	prate_sfc[ASTEP];		/* precipitation, mm mon-1 */
	double 	prate_sfc_ann;			/* precipitation, mm mon-1 */
	double 	spfh_2m[ASTEP];			/* specific humidity, kg kg-1 */
	double 	wnd_10m[ASTEP];			/* wind velocity, m s-1 */

	double 	tmp_sfc_am;				/* mean temperature, degree Celcius */
	double 	tmp_sfc_mx;				/* maximum temperature, degree Celcius */
	double 	tmp_sfc_mn;				/* minimum temperature, degree Celcius */
	double 	tmp_soil_mean;
	double 	gp_atem;				/* average temperature during the growing period, degree Celcius */
	double 	gp_tem;					/* average temperature during the growing period, degree Celcius */
	double 	gp_pre;					/* precipitation during the prowing period, mm */
	
	/* climate condition: *_a[] means the average during 1965 to 1998 */	
	double 	tmp_sfc_a[ASTEP];		/* ground surface temperature, degree Celcius */
	double 	tmp_2m_a[ASTEP];		/* 2m air temperature, degree Celcius */
	double 	tmp10_soil_a[ASTEP];	/* soil temperature at 10 cm depth, degree Celcius */
	double 	tmp200_soil_a[ASTEP];	/* soil temperature at 200 cm depth, degree Celcius */
	double 	tcdc_clm_a[ASTEP];		/* total cloudiness, fraction */
	double 	prate_sfc_a[ASTEP];		/* precipitation, mm mon-1 */
	double 	spfh_2m_a[ASTEP];		/* specific humidity, kg kg-1 */
	double 	ugrd_10m_a[ASTEP];		/* zonal wind velocity, m s-1 */
	double	vgrd_10m_a[ASTEP];		/* meridional wind velocity, m s-1 */
	double 	rad_a[ASTEP];			/* solar radiation, W m-2 */
	double 	par_a[ASTEP];			/* PAR, micro mol m-2 s-1 */
    
    double 	tmp_soil_am;
	
	double 	prec_sub_a[ASTEP];		/* precipitation from substitute data (UEA/CRU), mm mon-1 */

	double 	topo;					/* topography, orology, and altitude, m above MSL */
	/* double 	whc;				*/	/* soil water holding capacity */
	/* double 	whc30;				*/	/* soil water holding capacity of above 300mm soil */
	double 	sd;						/* soil rooting depth */
	double 	hyd_cond;				/* hydraulic conductivity, */
	
	double	sl_dec[ASTEP];			/* solar declination, degree */
	double	sl_hgt[ASTEP];			/* solar hight at midday, degree */
	double	dlen[ASTEP];			/* monthly day length, hour */
	double	top_rad[ASTEP];			/* downward solar radiation at the atmosphere-top, W m-2 */
	double	gl_rad[ASTEP];			/* global radiation at the canopy-top, W m-2 */
	double	par[ASTEP];				/* monthly photosynthetically active radiation, micro mol photon m-2 s-1 */

	double	par_be[ASTEP];			/* PAR, beam, in energy W/m2 */
	double	par_de[ASTEP];			/* PAR, diffuse, in energy W/m2 */
	double	par_bp[ASTEP];			/* PAR, beam, in photon micro mol photon /m2 /s */
	double	par_dp[ASTEP];			/* PAR, diffuse, in photon micro mol photon /m2 /s */
	
	/* GCM climate: year x month x row x column */
	/* 160 x 320   */
#if GCM_RUN==1
	float	proj_tmp2m[DL_GCM][ASTEP][160][320];		/* temperature */
	float	proj_prec[DL_GCM][ASTEP][160][320];			/* precipitation */
	float	proj_hum[DL_GCM][ASTEP][160][320];			/* specific humidity */
	float	proj_rad[DL_GCM][ASTEP][160][320];			/* surface downward solar radiation */
	/* 1970-1999 average: month x row x column [ASTEP][160][320] */
	float	proj_tmp2m_b[ASTEP][160][320];				/* temperature */
	float	proj_prec_b[ASTEP][160][320];				/* precipitation */
	float	proj_hum_b[ASTEP][160][320];				/* specific humidity */
	float	proj_rad_b[ASTEP][160][320];				/* surface downward solar radiation */
#else	
	float	proj_tmp2m[1][1][1][1];			/* temperature */
	float	proj_prec[1][1][1][1];			/* precipitation */
	float	proj_hum[1][1][1][1];			/* specific humidity */
	float	proj_rad[1][1][1][1];			/* surface downward solar radiation */
	/* 1970-1999 average: month x row x column */
	float	proj_tmp2m_b[1][1][1];			/* temperature */
	float	proj_prec_b[1][1][1];			/* precipitation */
	float	proj_hum_b[1][1][1];			/* specific humidity */
	float	proj_rad_b[1][1][1];			/* surface downward solar radiation */
#endif	
	
	/* NCEP/NCAR 1948-2011 */
#if NCEP_RUN==1	
	float	ncep_tmp2m[DL_NCEP][ASTEP][94][192];		/* temperature */
	float	ncep_prate[DL_NCEP][ASTEP][94][192];		/* precipitation */
	float	ncep_tcdc[DL_NCEP][ASTEP][94][192];			/* total cloudiness */
	float	ncep_vpres[DL_NCEP][ASTEP][94][192];		/* vapor pressure */
	/* average */
	float	ncep_tmp2m_b[ASTEP][94][192];
	float	ncep_prate_b[ASTEP][94][192];
	float	ncep_tcdc_b[ASTEP][94][192];
	float	ncep_vpres_b[ASTEP][94][192]; 
#else	
	float	ncep_tmp2m[1][1][1][1];
	float	ncep_prate[1][1][1][1];
	float	ncep_tcdc[1][1][1][1];
	float	ncep_vpres[1][1][1][1]; 
	/* average */
	float	ncep_tmp2m_b[1][1][1];
	float	ncep_prate_b[1][1][1];
	float	ncep_tcdc_b[1][1][1];
	float	ncep_vpres_b[1][1][1]; 
#endif

    /* albedo perturbation: 2012/12/29 by A.Ito */
#if EX_ALBEDO==0
    float   albedo_av[1][1][1];
    float   albedo_sd[1][1][1];
    float   albedo_max[1][1][1];
    float   albedo_min[1][1][1];
#else
    float   albedo_av[ASTEP][36][72];
    float   albedo_sd[ASTEP][36][72];
    float   albedo_max[ASTEP][36][72];
    float   albedo_min[ASTEP][36][72];
#endif
    float   albedo_pert[12];
    float   glbalbedo[ASTEP];

	double	proj_prec_co;			/* carry-over of negative precipitation */

	/* historical (e.g., UEA/CRU TS2.1) data */
	long	flag_histdata;							/* flag of data availability */
	double	hist_tmp[DL_CRU][ASTEP];			/* temperature */
	double	hist_pre[DL_CRU][ASTEP];			/* precipitation */
	double	hist_cld[DL_CRU][ASTEP];			/* cloud cover */
	double	hist_vap[DL_CRU][ASTEP];			/* vapor pressure */
	/* historical average */
	double	hist_tmp_b[ASTEP];					/* temperature */
	double	hist_pre_b[ASTEP];					/* precipitation */
	double	hist_cld_b[ASTEP];					/* cloud cover */
	double	hist_vap_b[ASTEP];					/* vapor pressure */

	/* erosion */
	long 	rvbasin;						/* ID of river basin */
	double 	albedo_soil;					/* soil albedo */
	double 	pcnt_orgmat;					/* percent organic matter */
	double 	fls_slope;						/* slope factors LS for RUSLE erosion model */
	double 	fk_edodibility;					/* soil erodibility for RUSLE erosion model */
	
	/* land-use change */
	double 	fcrop_sage[293];				/* historical data, 1700-1992, Ramankutty & Foley */
	/* future, 1990-2100, IMAGE2, Wang et al. */
	double 	fcrop3_image[111];				/* C3 crop */
	double 	fcrop4_image[111];				/* C4 crop */
	double 	fgrass3_image[111];				/* C3 pasture grass */
	double 	fgrass4_image[111];				/* C4 pasture grass */
	
	double	fcrop;			/* fraction of cropland */
	double	frice;			/* fraction of rice */
	double	fwheat;			/* fraction of wheat */
	double	fmaize;			/* fraction of maize */
	double	fothers;		/* fraction of other crops */
	
	/* EOS-WEBSTER, 1700-2000/2005, Hurtt et al. */
	double	fcrop_unh_hmnzed[DL_LUH];		/* cropland fraction */
	double	fpast_unh_hmnzed[DL_LUH];		/* pasture fraction */
	double	fprim_unh_hmnzed[DL_LUH];		/* primary land fraction */
	double	fsecd_unh_hmnzed[DL_LUH];		/* secondary land fraction */
	double	ssma_unh_hmnzed[DL_LUH];		/* secondary land property 1 */
	double	ssmb_unh_hmnzed[DL_LUH];		/* secondary land property 2 */
	double	t_cp_unh_hmnzed[DL_LUH];		/* conversion crop to pasture */
	double	t_cs_unh_hmnzed[DL_LUH];		/* conversion crop to secondary */
	double	t_pc_unh_hmnzed[DL_LUH];		/* conversion pasture to crop */
	double	t_ps_unh_hmnzed[DL_LUH];		/* conversion pasture to secondary */
	double	t_sc_unh_hmnzed[DL_LUH];		/* conversion secondary to crop */
	double	t_sp_unh_hmnzed[DL_LUH];		/* conversion secondary to pasture */
	double	t_ss1_unh_hmnzed[DL_LUH];		/*  */
	double	t_ss2_unh_hmnzed[DL_LUH];		/*  */
	double	t_ss3_unh_hmnzed[DL_LUH];		/*  */
	double	t_vc_unh_hmnzed[DL_LUH];		/* conversion primary to crop */
	double	t_vp_unh_hmnzed[DL_LUH];		/* conversion primary to pasture */
	double	t_vs1_unh_hmnzed[DL_LUH];		/*  */
	double	t_vs2_unh_hmnzed[DL_LUH];		/*  */
	
	double 	f_crop_con;					/* contemporary cropland fraction */
	double 	f_crop_p;					/* previous cropland fraction */
	double 	f_crop_trend;				/* trend of cropland conversion, fraction/yr */
	double 	f_pasture_con;				/* contemporary pasture fraction */
	double 	f_pasture_p;				/* previous pasture fraction */
	double 	f_pasture_trend;			/* trend of pasture conversion, fraction/yr */
	
	double	f_deforest;					/* deforestation rate */
	double	f_deforest_v;				/* in primary lands */
	double	f_deforest_s;				/* in secondary lands */
	
	double	fcrop_rk[308];				/* crop fraction by Ramankutty & Kimball (2010) */
	double	fpast_rk[308];				/* pasture fraction by Ramankutty & Kimball (2010) */
	
	double	f_crop_base;				/* base cropland fraction in 2000 */
	double	f_pasture_base;				/* base pasture fraction in 2000 */
	
	/* wood harvest */
	double	hvst_p1[DL_LUH];
	double	hvst_p2[DL_LUH];
	double	hvst_s1[DL_LUH];
	double	hvst_s2[DL_LUH];
	double	hvst_s3[DL_LUH];

	/* RUSLE erosion model coefficients */
	double 	f_erosion_r;				/* rain factor */
	double 	f_erosion_ls;				/* slope factors */
	double 	f_erosion_k;				/* erodibility factor */
	double 	f_erosion_c;				/* vegetation cover factor */
	double 	f_erosion_p;				/* protection factor */
	
	/* soil parameters for CH4 oxy */
	long 	soiltexture;				/* soil texture */
	double 	field_cap1;					/* soil field capacity, upper */
	double 	field_cap2;					/* soil field capacity, lower */
	double 	pore_cap1;					/* pore capacity */
	double 	a_sw;						/* Saxton soil physical parameter, a */
	double 	b_sw;						/* Saxton soil physical parameter, b */
	
	double 	fc_30;						/* field capacity, 0-30 cm */
	double 	fc_150;						/* field capacity, 30-150 cm */
	double 	pc_sand;					/* sand fraction */
	double 	pc_clay;					/* clay fraction */
	double 	bulkdens;					/* bulk density, kg/litre */
	double 	soil_ph;					/* soil pH */
	double 	total_n_1m;					/* total soil N in upper 1m, gN/m2 */
	
	double 	f_paddy;					/* paddy land fraction */
	double 	f_paddy_b;					/* paddy land fraction */
	double 	f_wetland;					/* wetland fraction */
	double 	f_upland;					/* upland (e.g. forest, grassland) fraction */
	double 	f_lake;						/* lake fraction */
	double 	f_wetland0;					/* wetland fraction (base) */
	
	/* nitrogen deposition */
	double 	ndepo[3];					/* N deposition by Galloway et al. (2004) */
	
	/* CHASER 2001 monthly, by A.Ito (2010/05/21) */
	double	ndepo_chaser_dnhx[ASTEP][64][128];		/* NHx, dry */
	double	ndepo_chaser_dnoy[ASTEP][64][128];		/* NOy, dry */
	double	ndepo_chaser_wnhx[ASTEP][64][128];		/* NHx, wet */
	double	ndepo_chaser_wnoy[ASTEP][64][128];		/* NOy, wet */

	/* CHASER4.0 monthly, by A.Ito (2014/11/19) */
	double	ndepo_chaser4_nhx_h[ASTEP][64][128];		/* NHx */
	double	ndepo_chaser4_noy_h[ASTEP][64][128];		/* NOy */
	double	ndepo_chaser4_ont_h[ASTEP][64][128];		/* Org NOx */
	double	ndepo_chaser4_nhx_p[ASTEP][64][128];		/* NHx */
	double	ndepo_chaser4_noy_p[ASTEP][64][128];		/* NOy */
	double	ndepo_chaser4_ont_p[ASTEP][64][128];		/* Org NOx */
	
	/* radiation conversion model using SRB data */
	double	srb_dif_aa;					/* linear regression a */
	double	srb_dif_bb;					/* linear regression b */
	double	srb_dif_rr;					/* linear regression r */
	double	srb_dif_min_x;
	double	srb_dif_max_x;
	double	srb_dif_min_y;
	double	srb_dif_max_y;
	
	double	inundation_ssmi[ASTEP];		/* inundation by SSM/I: added by A.Ito (2009/07/13) */
	double	inundation_ssmi_av;
	double	inundation_ssmi_max;
    
    double  inundation_gcp_av[ASTEP];
    double  inundation_gcp_ts[15][ASTEP];
	
	long	type_permaforst;			/* permafrost type by NSIDC */
    double  tmp_base_permaforst;        /* 2012/10/26 by A.Ito */
    
    double  f_biofuel[DL_BF];                /* biofuel scenario: 2015/8/21 by A.Ito */
    
    /* NMIP input: 2015/11/19 by A.Ito */
    double  nmip_nfert[DL_NMIP];                /* nitrogen fertilizer */
    double  nmip_ndep_noy[DL_NMIP];                /* NOy deposition */
    double  nmip_ndep_nh4[DL_NMIP];                /* NH4 fertilizer */
    double  nmip_manure[DL_NMIP];                /* manure */
};

/* grid conditions, derived from submodules *******************************************/
struct Loct{ 
	short	v_type;                     /* vegetation classification types */
                                        /* 1: Olson+SAGE natural vegetation */
                                        /* 2: agricultural vegetation */
	
	long 	time_hyd;					/* time to reach stabilization of water budget */
	long 	time;						/* time to reach stabilization of carbon budget */

	double	aco2[ASTEP];				/* ambient CO2 concentration, in ppmv */
	double	d13c_aco2[ASTEP];			/* stable carbon isotope composition of CO2, dimensionless */
	double	cnpy_co2_recyc;				/* within-canopy CO2 recycling ratio */
    double  ao3;                        /* ambient O3, ppb */
	
	double	c4ptn[ASTEP];				/* ground coverage of C4 plants, fraction */
	double	c3ptn[ASTEP];				/* ground coverage of C3 plants, fraction */
	long	gd[ASTEP], bbm;				/* vegetative growing period, days */
	double	gdd[ASTEP];				 	/* cumulative growth degree days, degC days */
    
    double  est_maxlai;                 /* estimated max.LAI: 2014/05/20 by A.Ito */
	
	double	albedo_sfc[ASTEP];			/* land-surface albedo */
	double	gl_rad_g[ASTEP];			/* global radiation under the canopy, W m-2 */
	double	rad_net_p[ASTEP];			/* net radiation, canopy, W m-2 */
	double	rad_net_g[ASTEP];			/* net radiation, soil surface, W m-2 */
	double	rad_net_long[ASTEP];		/* net long-wave radiation, W m-2 */
	double	rad_net_short[ASTEP];		/* net short-wave radiation, W m-2 */
	double	rad_net[ASTEP];				/* net radiation, W m-2 */
	double	rdi;						/* radiative dryness index by Budyko */
    
    /* added: 2013/01/10 by A.Ito */
    double  grad_d[ASTEP];              /* daily average downward SW radiation, W m-2 */
    double  nsw_d[ASTEP];               /* daily average net SW radiation, W m-2 */
 	
    double  ppfd_h[DSTEP];
    double  ppfdb_h[DSTEP];
    double  ppfdd_h[DSTEP];

    double  ippfd_g[ASTEP];
    double  appfd_g[ASTEP];
	double	fappfd_g[ASTEP];
	double	fapar_df[ASTEP];

	double	pet_prty[ASTEP];			/* Priestley-Taylor potential evapotranspiration, mm month-1 */
	double	pet_prty_ann;				/* annual Priestley-Taylor potential evapotranspiration, mm yr-1 */
	
	double	prsr[ASTEP];				/* air pressure, hPa */
	double	dnsa[ASTEP];				/* density of air, kg m-3 */

	double	vp[ASTEP];					/* vapour pressure, hPa */
	double	vps[ASTEP];					/* saturation vapour pressure, hPa */
	double	slope_vps[ASTEP];			/* slope of saturation vapour pressure related to tempertaure, hPa deg C-1 */
	double	r_aero[ASTEP];				/* aerodynamic resistance, s m-1 */
	double	vpd[ASTEP];					/* vapour pressure deficit, hPa */
	
	double	lai[ASTEP];					/* leaf area index, m2 m-2 */
	double	canopy_con[ASTEP];			/* canopy conductance, mmol H2O m-2 s-1 */
	double	f_vegcov[ASTEP];			/* fractional vegetation cover */
	
	/* water pools, mm (= kg/m2) */
	double	snwa;						/* water equivalent snow depth, mm */
	double	msnwa[ASTEP];				/* monthly */
	double	sw30;						/* soil water content of above 30cm soil, mm */
	double	msw30[ASTEP];				/* monthly */
	double	sww;						/* whole soil water content, mm */
	double	msww[ASTEP];				/* monthly */
    
    double  b_sw30[ASTEP];              /* baseline soil water, 0-30cm */
    double  b_sww[ASTEP];               /* baseline soil water, 30-cm */

	/* water fluxes, mm / month */
	double	pm_evp[ASTEP];				/* potential soil evaporation rate, mm */
	double	pm_trn[ASTEP];				/* potential transpiration rate, mm */
	double	pm_incep[ASTEP];			/* potential intercepted-water evaporation rate, mm */
	double	incep[ASTEP];			
	double	evpr[ASTEP];				/* actual evaporation rate, mm */
	double	trspr[ASTEP];				/* actual transpiration rate, mm */
	double	ro1[ASTEP];					/* runoff from upper soil water, mm */
	double	ro2[ASTEP];					/* runoff from lower soil water, mm */
	double	pntrt[ASTEP];				/* water penetration from upper to lower layer, mm */
	double	thaw[ASTEP];				/* snow thaw water, mm */
	double	vmc30[ASTEP];				/* volumatric moisture content of upper layer, fraction */
	double	vmc[ASTEP];					/* volumatric moisture content of lower layer, fraction */
	double	snp[ASTEP];					/* snow fraction of precipitation */
	double	soil_appr30;				/* soil aperture of upper layer, fraction */
	double	soil_apprw;					/* soil aperture of lower layer, fraction */
	
	double	n_frtlz_in;					/* N-fertilization input */
	double	depo_no3[ASTEP];			/* NO3- deposition */
	double	depo_nh4[ASTEP];			/* NH4+ deposition */
	
	/* CASA moisture **********************************/
	/* long	mday;			*/
	double	m_m[ASTEP];					/* soil moisture index */
	double	m_m_pre;					/* precipitation */
	double	m_rdr[ASTEP];				/* dryness */
	double	m_vmc[ASTEP];				/* volumetric soil moisture */
	double	m_pet[ASTEP];				/* potential evapotranspiration */
	double	m_e[ASTEP];					/* soil water change */
	double	m_sw[ASTEP];				/* soil wetness */
	double	i_w[ASTEP];					/* moisture scalar */
	/* double	d_tmp[31];			
	double	d_sw[31];			
	double	d_vmc[31];			*/
	double	wfps[ASTEP];				/* water-filled pore space */
	
	double	f_inund_wet_wh[ASTEP];	/* inundation area for Wlater & Heimann CH4 scheme */
	double	f_inund_pad_wh[ASTEP];	/* inundation area for Wlater & Heimann CH4 scheme */
	
	/* maximum GPP for Cao CH4 scheme */
	double	gpp_max;						/* maximum GPP */
    double  npp_av[ASTEP];
	
	/* CH4 emission by Walter & Heimann: added by A.Ito (2009/08/05) */
	double	water_table_depth;				/* current time-step */
	double	water_table_depth_pre;			/* previous time-step */
	double	npp_max;						/* maximum NPP */
	double	prof_ch4[SOIL_LAYER+2];			/* CH4 concentration profile */
	
	double	cum_dprec;					/* cumulative precipitation change */
	
	/* tentative variables for debugging and monitoring */
	double	xx1[ASTEP];
	double	xx2[ASTEP];
	double	xx3[ASTEP];
	double	xx4[ASTEP];
	double	xx5[ASTEP];
	double	xx6[ASTEP];
	double	xx7[ASTEP];
	double	xx8[ASTEP];
	double	xx9[ASTEP];
};

/* vegetation characteristics ****************************************************/
struct Pchar{ 
	short	v_type;				/* vegetation classification types */
								/* 1: Olson+SAGE natural vegetation */
								/* 2: agricultural vegetation */
								
	/* optics */
	double	albedo;					/* reflectivity, or albedo */
	double	apar_bp[ASTEP];			/* absorbed PAR photon, beam, micro mol photon m-2 s-1 */
	double	apar_dp[ASTEP];			/* absorbed PAR photon, diffuse, micro mol photon m-2 s-1 */
	double	fapar[ASTEP];			/* fraction of absorbed PAR */

    double  ppfd_db[ASTEP];
    double  appfd_db[ASTEP];

	/* allocation *********/
	double	opt_lai[ASTEP];			/* optimum leaf area index */
	double	alloc_ass;				/* allocation ration for assimilation organ, fraction */
	double	alloc_abg;				/* allocation ration for aboveground non-assimilation organ, fraction */
	double	malloc_f[ASTEP];		/* monthly allocation ratio to leaf */
	double	malloc_c[ASTEP];		/* monthly allocation ratio to stem */
	double	malloc_r[ASTEP];		/* monthly allocation ratio to root */
	
	/* N parameters *******/
	double	n_salvage;				/* N re-uptake ratio */
	double	cn0_strg;				/* base C/N of storage nitrogen */
	double	n_conc_larea;			/* leaf N concentration, mmol N m-2 leaf area */
	double	n_conc_larea_m[ASTEP];	/* monthly */
	double	n_conc_lmass;			/* leaf N concentration, mmol N g-1 leaf weight */
	double	kn_nphoto;				/* canopy N gradient, exponential coefficient */
	double	amax_nphoto;			/* max. photosynthesis rate, micro mol CO2 m-2 s-1 */
	double	rd_n;					/* dark respiration, micro mol CO2 m-2 s-1 */
	
	double	n_leaf_df97;			/* leaf N concentration, mmol N m-2 */
	
	/* phenology *********/
	double	gdd;					/* growing degree days */
	double	mgdd[ASTEP];			/* monthly growing degree days */
	double	grw_pd;					/* growing period */
	long 	season[ASTEP];			/* phenological stage as a function of season */
	 	/* 0: dormancy */
	 	/* 1: vegetative growth */
	 	/* 2: emergence of new leaf */
	 	/* 3: abandon of old leaf */
	long 	frag_emg;				/* flag of leaf emergence */
	double	frag_dcd;				/* flag of deciduou */
	
	double	fleaf_age[49];			/* leaf age distribution, 0(new) to 48-over(old) months */
	
	/* critical temperature condiction for bur burst and leaf shedding */
	double	crit_temp;
	double	crit_gdd;
    
	/* photosynthesis *******/
	short	phototype;				/** photosynthetic metabolic pathway, 3=C3, 4=C4, 5=CAM **/
	/*  veg->psat[grid->m] = veg->pmax*ftem*fstl*fnstl   */
	double	psat[ASTEP];			/* light-saturated rate, micro mol CO2 m-2 s-1 */
	double	pmax;					/* potential maximum rate, micro mol CO2 m-2 s-1 */
	double	ft[ASTEP];				/* temperature coefficient */
	double	fcd[ASTEP];				/* CO2 coefficient */
	double	fsw[ASTEP];				/* soil water coefficient */
    double  fo3[ASTEP];             /* O3 coefficient: 2013/02/25 by A.Ito */
	double	ptop;					/* canopy-top photosynthetic rate */
	double	sla;					/* specific leaf area, cm2 g dm-1 */
	double	eK0;					/* light attenuation coefficient, no dimension */
	double	eK[ASTEP];				/* light attenuation coefficient, no dimension */
	double	lue0;					/* control light dependence coefficient, mol CO2 mol photon-1 */
	double	lue[ASTEP];				/* monthly quantum yield , mol CO2 mol photon-1 */ 
	double	topt;					/* optimum temperature, deg C */
	double	topt0;					/* optimum temperature (co2 dependent), deg C */
	double	tmin;					/* minimum temperature, deg C */
	double	tmax;					/* maximum temperature, deg C */
	double	ci[ASTEP];				/* monthly intercellular CO2 concentration, ppmv */
	double	kmci;					/* dependence of photosynthesis on intercellular CO2 concentration, ppmv */
	double	cmpcd0;					/* CO2 compensation point, ppmv */
	double	cmpcd[ASTEP];			/* CO2 compensation point, ppmv */
	/* parameters of conductance */
	double	gs[ASTEP];				/* monthly stomatal conductance, mmol H2O m-2 s-1 */
	double	gc[ASTEP];				/* monthly canopy conductance, mmol H2O m-2 s-1 */
	double	gs_b0;					/* Leuninig stomata model parameters */
	double	gs_b1;					/* Leuninig stomata model parameters */
	double	gs_b2;					/* Leuninig stomata model parameters */
	double	km_nstl;				/* maximum stomatal conductance */
	
	double	psat_df[ASTEP];
	double	lue_df[ASTEP];
    
    /* O3 effect: 2013/02/23 by A.Ito */
    double  f_o3;
	
	/*** photosynthesis: de Pury and Farquhar (1997) ***/
	
	/*** respiration ***/
	/* specific growth respiration rate, g C g C-1 alloc */
	double	rgf;					/* leaf */
	double	rgc;					/* stem */
	double	rgr;					/* root */
	/* specific maintenance respiration rate at 15 degC, mg C g C-1 day-1 */
	double	rmf;					/* leaf, ecosystem */
	double	rmc;					/* stem, ecosystem */
	double	rmr;					/* root, ecosystem */
	double	rmf0;					/* leaf */
	double	rmc_s;					/* sap wood */
	double	rmr_s;					/* fine root */
	double	rmc_h;					/* heart wood */
	double	rmr_h;					/* coarse root */
	/* temperature dependence Q10, dimensionless *******/
	double	qTf[ASTEP];				/* leaf */
	double	qTc[ASTEP];				/* stem */
	double	qTr[ASTEP];				/* root */
	/* temperature dependence Q10 at 15 deg C *****/
	double	qTf0;					/* leaf */
	double	qTc0;					/* stem */
	double	qTr0;					/* root */
	
	/* litter fall ********/
	/* specific litter fall rate, fraction */
	double	lf[ASTEP];				/* leaf */
	double	lc[ASTEP];				/* stem */
	double	lr[ASTEP];				/* root */
	/* base rate *********/
	double	lf0;					/* leaf */
	double	lc0;					/* stem */
	double	lr0;					/* root */
	/* deciduous leaf fraction *****/
	double	dcd;				
	
	/* photosynthetgic 13C fractionation, permil *******/
	double	photo_13c_frac[ASTEP];
		
	/* specific constants *********************************/
	/* root stratification parameters by Zeng (2001) *******/
	double	root_dist_a;			/* root profile parameter a, m-1 */
	double	root_dist_b;			/* root profile parameter b, m-1 */
	double	root_depth;				/* rooting depth */
};			

/* soil characteristics *************************************************/
struct Schar{ 
	short	v_type;				/* vegetation classification types */
								/* 1: Olson+SAGE natural vegetation */
								/* 2: agricultural vegetation */

	double	albedo0;			/* reflectivity, or albedo */
	double	albedo[ASTEP];		/* monthly */
	/* specific respiration rate at 15degC, mg C g dm-1 day-1 */
	double	rl;					/* litter */
	double	rh;					/* humus */
	double	rl0;				/* litter, base */
	double	rh0;				/* humus, base */
	double	qTl;				/* temperature dependence Q10, litter */
	double	qTh;				/* temperature dependence Q10, humus */
	double	kml;				/* moisture dependence, litter */
	double	kmh;				/* moisture dependence, humus */
	double	kmsl;				/* aerobic dependence, litter */
	double	kmsh;				/* aerobic dependence, humus */
	double	me;					/* mineral soil formation ratio to litter decomposition, fraction */
	
	/* environmental scalars, for acclimation sensitivity analysis */
	double	ft_l[ASTEP];		/* temperature for litter */
	double	ft_h[ASTEP];		/* temperature for humus */
	double	fm_l[ASTEP];		/* moisture for litter */
	double	fm_h[ASTEP];		/* moisture for humus */
	/* with "0", average value for sensitivity analysis */
	double	ft0_l[ASTEP];			
	double	ft0_h[ASTEP];
	double	fm0_l[ASTEP];
	double	fm0_h[ASTEP];
};			

/* ecosystem characteristics *******************************************/
struct Echar{ 
	struct	Pchar c3;			/* for C3 plants */
	struct	Pchar c4;			/* for C4 plants */
	struct	Schar soil;			/* for soil organic matter */
};			

/* plant biomass *******************************************************/
struct Pmas{ 
	short	v_type;				/* vegetation classification types */
								/* 1: Olson+SAGE natural vegetation */
								/* 2: agricultural vegetation */

	/* leaf area index (LAI), m2 m-2  */
	double	lai[ASTEP];			/* monthly value */
	double	lai0[ASTEP];		/* 1990's average value for sensitivity analysis: 2009/05/06 by A.Ito */
	
	/* carbon mass, Mg C ha-1 *****/ 
	double	fol;				/* leaf */
	double	mfol[ASTEP];		/* monthly */
	double	stm;				/* stem and branch */
	double	mstm[ASTEP];		/* monthly */
	double	rot;				/* root */
	double	mrot[ASTEP];		/* monthly */
	
	double	plant[ASTEP];		/* monthly plant mass */
	double	fol_p;				/* previous leaf carbon */
	
	/* stable carbon isotope composition, permille */
	double	d13c_fol;			/* leaf */
	double	d13c_mfol[ASTEP];	/* monthly */
	double	d13c_stm;			/* stem */
	double	d13c_mstm[ASTEP];	/* monthly */
	double	d13c_rot;			/* root */
	double	d13c_mrot[ASTEP];	/* monthly */
	
	double	d13c_plant[ASTEP];		/* total plant */
	
	/* radio isotope 14C: added by A.Ito (2009/06/23) */
	double	d14c_fol;			/* leaf */
	double	d14c_mfol[ASTEP];	/* monthly */
	double	d14c_stm;			/* stem */
	double	d14c_mstm[ASTEP];	/* monthly */
	double	d14c_rot;			/* root */
	double	d14c_mrot[ASTEP];	/* monthly */

	/* plant N, g N ha-1 */
	double	n_cnpy;				/* N in canopy */
	double	n_cnpy_m[ASTEP];	/* monthly */	
	double	n_strg;				/* N in storage pool */
	double	n_strg_m[ASTEP];	/* monthly */
};			

/* soil carbon storage *************************************************/
struct Smas{ 
	short	v_type;				/* vegetation classification types */
								/* 1: Olson+SAGE natural vegetation */
								/* 2: agricultural vegetation */

	/* carbon mass, Mg C ha-1 */ 
	double	ltr;				/* litter mass */
	double	ltr_m[ASTEP];		/* monthly litter mass */
	double	msl;				/* mineral soil and humus mass */
	double	msl_m[ASTEP];		/* monthly mineral soil and humus mass */
	double	soil[ASTEP];		/* monthly soil mass */
	
	double	doc;				/* dissolved organic carbon discharge */
	double	doc_m[ASTEP];		/* monthly */
	
	/* stable carbon isotope composition, permille */
	double	d13c_ltr;			/* litter */
	double	d13c_ltr_m[ASTEP];	/* monthly */
	double	d13c_msl;			/* mineral soil */
	double	d13c_msl_m[ASTEP];	/* monthly */
	double	d13c_soil[ASTEP];	/* total soil */
	
	/* radio isotope 14C: added by A.Ito (2009/06/23) */
	double	d14c_ltr;			/* litter */
	double	d14c_ltr_m[ASTEP];	/* monthly */
	double	d14c_msl;			/* mineral soil */
	double	d14c_msl_m[ASTEP];	/* monthly */

	/* soil inorganic N, g N ha-1 */
	double	n_no3;				/* NO3- */
	double	n_no3_m[ASTEP];		/* monthly */
	double	n_nh4;				/* NH4+ */
	double	n_nh4_m[ASTEP];		/* monthly */
	
	/* soil organic N, g N ha-1 */
	double	n_mcrb;				/* microbe */
	double	n_mcrb_m[ASTEP];	/* monthly */
	double	n_lttr;				/* soil organic litter */
	double	n_lttr_m[ASTEP];	/* monthly */
	double	n_hums;				/* soil organic humus */
	double	n_hums_m[ASTEP];	/* monthly */
};

/* ecosystem carbon storage *************************************/
struct Mass{ 
	struct	Pmas c3;			/* C3 plant mass */
	struct	Pmas c4;			/* C4 plant mass */
	struct	Pmas plant;			/* all plant mass */
	struct	Smas soil;			/* soil mass */

	/* carbon mass */ 
	double	total[ASTEP];		/* ecosystem total carbon storage, Mg C ha-1 */

	double	lai_p;				/* previous LAI, m2 m-2 */
	
	/* stable carbon isotope composition */
	double	d13c_total[ASTEP];		/* total d13C, permille */
};			

/* plant carbon fluxes, all monthly *******************************/
struct Pflx{ 
	short	v_type;				/* vegetation classification types */
								/* 1: Olson+SAGE natural vegetation */
								/* 2: agricultural vegetation */

	/* carbon flux, in Mg C ha-1 mon-1 */ 
	double	gpp[ASTEP];			/* gross primary production */
	double	spp[ASTEP];			/* net surplus production */
	double	epp[ASTEP];			/* net effective production */
	double	npp[ASTEP];			/* net primary production */
	
	double	gpp_df97[ASTEP];	/* gross primary production by de Pury & Farquhar scheme */

	double	tpf[ASTEP];			/* translocation of photosynthate to foliage */
	double	tpc[ASTEP];			/* translocation of photosynthate to stem */
	double	tpr[ASTEP];			/* translocation of photosynthate to root */
	double	tpp[ASTEP];			/* translocation of photosynthate, total */
 	
	double	ar[ASTEP];			/* plant respiration, =rpm+rpg */
	double	arg[ASTEP];			/* plant growth respiration */
	double	arm[ASTEP];			/* plant maintenance respiration */
	double	rfg[ASTEP];			/* foliage growth respiration */
	double	rfm[ASTEP];			/* foliage maintenance respiration */
	double	rcg[ASTEP];			/* stem and branch growth respiration */
	double	rcm[ASTEP];			/* stem and branch maintenance respiration */
	double	rrg[ASTEP];			/* root growth respiration */
	double	rrm[ASTEP];			/* root maintenance respiration */
	
	double	lf[ASTEP];				/* foliage litterfall */
	double	lc[ASTEP];				/* stem and branch litterfall */
	double	lr[ASTEP];				/* root litterfall */
	double	lL[ASTEP];				/* total litterfall */
	double	lf_c[ASTEP];			/* leaf shedding in C3/C4 altyeration in grassland */
	
	double	hvst[ASTEP];			/* harvest of crops */
	
	double	emit_ch4_kirschbaum_mass[ASTEP];		/* plant CH4 emission, mass-based */
	double	emit_ch4_kirschbaum_photo[ASTEP];		/* plant CH4 emission, photosynthesis-based */
	
	/* stable carbon isotope composition, d13C, permille ******/
	/* production */
	double	d13c_gpp[ASTEP];		/* GPP */
	double	d13c_spp[ASTEP];		/* SPP */
	double	d13c_epp[ASTEP];		/* EPP */
	double	d13c_npp[ASTEP];		/* NPP */
	/* translocation */
	double	d13c_tpf[ASTEP];		/* leaf  */	
	double	d13c_tpc[ASTEP];		/* stem  */
	double	d13c_tpr[ASTEP];		/* root  */
	double	d13c_tpp[ASTEP];		/* total  */
	/* autotrophic respiration */
	double	d13c_ar[ASTEP];			/* total  */	
	double	d13c_arg[ASTEP];		/* growth  */	
	double	d13c_arm[ASTEP];		/* maintenance  */	
	double	d13c_rfg[ASTEP];		/* growth, leaf  */	
	double	d13c_rfm[ASTEP];		/* maintenance, leaf  */	
	double	d13c_rcg[ASTEP];		/* growth, stem  */
	double	d13c_rcm[ASTEP];		/* maintenance, stem  */
	double	d13c_rrg[ASTEP];		/* growth, root  */	
	double	d13c_rrm[ASTEP];		/* maintenance, root  */
	/* leaf litterfall */
	double	d13c_lf[ASTEP];			/* leaf  */	
	double	d13c_lc[ASTEP];			/* stem  */
	double	d13c_lr[ASTEP];			/* root  */
	double	d13c_lL[ASTEP];			/* total  */
	double	d13c_lf_c[ASTEP];		

	double	d13c_hvst[ASTEP];		/* harvest */

	double	d14c_gpp[ASTEP];		/* GPP */
	double	d14c_lL[ASTEP];			/* litter input */

	/* nitrogen flux, g N ha-1 mon-1 */
	double	n_biofix[ASTEP];		/* biological N2 fixation */
	double	uptake_no3[ASTEP];		/* root NO3- uptake */
	double	uptake_nh4[ASTEP];		/* root NH4+ uptake */
	double	n_alloc_cnpy[ASTEP];	/* N allocation to canopy */
	double	n_alloc_strg[ASTEP];	/* N allocation to storage */
	double	n_salvage[ASTEP];		/* N salvage from shedding leaves */
	double	n_realloc[ASTEP];		/* N reallocation from storage to leaves */
	double	n_abdn_cnpy[ASTEP];		/* N abandonment from canopy */
	double	n_abdn_strg[ASTEP];		/* N abandonment from storage */
};			

/* soil carbon fluxes  *******************************************************/
struct Sflx{	 
	short	v_type;				/* vegetation classification types */
								/* 1: Olson+SAGE natural vegetation */
								/* 2: agricultural vegetation */

	/* carbon flux, in Mg C ha-1 mon-1 */ 
	double	lL[ASTEP];			/* litter input */
	double	rl[ASTEP];			/* litter decomposition */
	double	rh[ASTEP];			/* mineral soil and humus decomposition */
	double	hr[ASTEP];			/* total decomposition */
	double	sf[ASTEP];			/* humus formation */
	
	/* dissolved organic carbon (DOC) discharge */
	double	doc_boyer[ASTEP];

	/* stable carbon isotope composition, permil */
	double	d13c_lL[ASTEP];			/* litter input */
	double	d13c_rl[ASTEP];			/* litter decomposition */
	double	d13c_rh[ASTEP];			/* humus decomposition */
	double	d13c_hr[ASTEP];			/* heterotrophic respiration */
	double	d13c_sf[ASTEP];			/* soil formation from litter to humus */
	
	double	d14c_lL[ASTEP];			/* litter input */
	
	/* CH4 oxydation, in mg CH4 m-2 month-1  */
	double	ch4oxy_ridg[ASTEP];				/* Ridgwell scheme */
	double	ch4oxy_casa[ASTEP];				/* CASA-Potter scheme */
	double	ch4oxy_delgrosso[ASTEP];		/* DayCnet-Del Grosso scheme */
	double	ch4oxy_curry[ASTEP];			/* Curry scheme */
	
	/* CH4 emission, mg CH4 m-2 month-1 */
	double	ch4prod_wetland_cao[ASTEP];		/* CH4 production, Cao scheme, at wetland */
	double	ch4oxy_wetland_cao[ASTEP];		/* CH4 oxidation, Cao scheme, at wetland */
	double	ch4flux_wetland_cao[ASTEP];		/* CH4 net flux, Cao scheme, at wetland */
	
	double	ch4prod_paddy_cao[ASTEP];		/* CH4 production, Cao scheme, at paddy field */
	double	ch4oxy_paddy_cao[ASTEP];		/* CH4 oxidation, Cao scheme, at paddy field */
	double	ch4flux_paddy_cao[ASTEP];		/* CH4 net flux, Cao scheme, at paddy field */
	
	/* Walter & Heimann *********/
	/* wetlands */
	double	ch4_wetland_wh_plant[ASTEP];	/* plant-mediated flux */
	double	ch4_wetland_wh_ebull[ASTEP];	/* ebullition flux */
	double	ch4_wetland_wh_diff[ASTEP];		/* diffusion flux */
	double	ch4_wetland_wh_release[ASTEP];	/* release */
	/* paddy fields */
	double	ch4_paddy_wh_plant[ASTEP];		/* plant-mediated flux */
	double	ch4_paddy_wh_ebull[ASTEP];		/* ebullition flux */
	double	ch4_paddy_wh_diff[ASTEP];		/* diffusion flux */
	double	ch4_paddy_wh_release[ASTEP];	/* release */
	
	/* N flow, g N ha-1 mon-1 **/
	double	n_abdn[ASTEP];					/* N abandonment as litter */
	
	double	d_n2o_ntr_ngas[ASTEP];			/* N2O from nitrification by NGAS */
	double	d_n2o_dnt_ngas[ASTEP];			/* N2O from denitrification by NGAS */
	double	d_n2o_ngas[ASTEP];				/* total N2O by NGAS */
	double	d_n2_ngas[ASTEP];				/* N2 by NGAS */
	
	double	d_no_casa[ASTEP];				/* NO by CASA */
	double	d_n2_casa[ASTEP];				/* N2 by CASA */
	double	d_n2o_casa[ASTEP];				/* N2O by CASA */
	
	double	n_nh3vlt[ASTEP];				/* NH3 volatilization */
	double	n_leach[ASTEP];					/* NO3- leaching */
	double	n_minerlz_lttr[ASTEP];			/* litter N mineralization */
	double	n_minerlz_hums[ASTEP];			/* humus N mineralization */
	double	n_nitrif[ASTEP];				/* nitrification in NGAS */
	double	n_immbl[ASTEP];					/* N immobilization */
	double	n_mcrb_abdn[ASTEP];				/* microbial abandonment */
	
	double	n_fertin[ASTEP];				/* fertilizer input */
};

/* ecosystem carbon fluxes *************************************************/
struct Flux{ 
	/* components */
	struct	Pflx c3;				/* C3 plant fluxes */
	struct	Pflx c4;				/* C4 plant fluxes */
	struct	Pflx plant;				/* total plant fluxes */
	struct	Sflx soil;				/* soil fluxes */
	
	/* empirical NPP estimation: Mg C ha-1 yr-1 */
	double	npp_miami;				/* NPP estimated by Miami model  */
	double	npp_montreal;			/* NPP estimated by Montreal model */
	double	npp_rosenzweig;			/* NPP estimated by Rosenzweig model */
	double	npp_schuur;				/* NPP estimated by Schuur model */
	double	npp_madison_parwsi;		/* NPP estimated by Zaks's Madison model 1: PAR & WSI */
	double	npp_madison_gddswsi;	/* NPP estimated by Zaks's Madison model 1: GDD & WSI */
	double	npp_madison_tp;			/* NPP estimated by Zaks's Madison model 1: temp & prec */
	double	npp_nceas;				/* NPP estimated by NCEAS */
	double	npp_chikugo;			/* NPP estimated by Chikugo */

	/* carbon flux, in  Mg C ha-1 yr-1 */ 
	double	nep[ASTEP];				/* net ecosystem production */
	double	nbp[ASTEP];				/* net biome production (added by A.Ito: 2010/01/20) */
	double	ncb[ASTEP];				/* net carbon balance of grid */
	double	lL0[ASTEP];				/* total litter fall */
	double	sr[ASTEP];				/* soil respiration */
	double	er[ASTEP];				/* ecosystem respiration */
	
	/* land-use change,  Mg C ha-1 yr-1 */
	double	lu_detr;				/* LUC emission from detritus */
	double	lu_conv;				/* LUC emission from 1-year pool */
	double	lu_ten;					/* LUC emission from 10-year pool */
	double	lu_hund;				/* LUC emission from 100-year pool */
	/* historical land-use-generated detritus production */
	double	detr_ten[10];			/* 10-year pool */
	double	detr_hund[100];			/* 100-year pool */
	
	/* biomass burning */
	double	f_burnt;				/* burnt fraction */
	double	day_fire[ASTEP];		/* days of fire */
	double	a_burnt[ASTEP];			/* area burnt */
	/* CO2 (g species) */
	double	bb_co2_litter[ASTEP];		/* from litter */
	double	bb_co2_leaf[ASTEP];			/* from leaf */
	double	bb_co2_wood[ASTEP];			/* from wood */
	double	bb_co2_root[ASTEP];			/* from root */
	/* CO (g species) */
	double	bb_co_litter[ASTEP];		/* from litter */	
	double	bb_co_leaf[ASTEP];			/* from leaf */
	double	bb_co_wood[ASTEP];			/* from wood */
	double	bb_co_root[ASTEP];			/* from root */
	/* CH4 (g species) */
	double	bb_ch4_litter[ASTEP];		/* from litter */
	double	bb_ch4_leaf[ASTEP];			/* from leaf */
	double	bb_ch4_wood[ASTEP];			/* from wood */
	double	bb_ch4_root[ASTEP];			/* from root */
	/* NMHC (g species) */
	double	bb_nmhc_litter[ASTEP];		/* from litter */
	double	bb_nmhc_leaf[ASTEP];		/* from leaf */
	double	bb_nmhc_wood[ASTEP];		/* from wood */
	double	bb_nmhc_root[ASTEP];		/* from root */	
	/* OC (g species) */
	double	bb_oc_litter[ASTEP];		/* from litter */
	double	bb_oc_leaf[ASTEP];			/* from leaf */
	double	bb_oc_wood[ASTEP];			/* from wood */
	double	bb_oc_root[ASTEP];			/* from root */
	/* BC (g species) */
	double	bb_bc_litter[ASTEP];		/* from litter */
	double	bb_bc_leaf[ASTEP];			/* from leaf */
	double	bb_bc_wood[ASTEP];			/* from wood */
	double	bb_bc_root[ASTEP];			/* from root */
	/* NOx (g species) */
	double	bb_nox_litter[ASTEP];		/* from litter */
	double	bb_nox_leaf[ASTEP];			/* from leaf */
	double	bb_nox_wood[ASTEP];			/* from wood */
	double	bb_nox_root[ASTEP];			/* from root */
	/* SO2 (g species) */
	double	bb_so2_litter[ASTEP];		/* from litter */
	double	bb_so2_leaf[ASTEP];			/* from leaf */
	double	bb_so2_wood[ASTEP];			/* from wood */
	double	bb_so2_root[ASTEP];			/* from root */
	/* PM2.5 (g species) */
	double	bb_pm25_litter[ASTEP];		/* from litter */
	double	bb_pm25_leaf[ASTEP];		/* from leaf */
	double	bb_pm25_wood[ASTEP];		/* from wood */
	double	bb_pm25_root[ASTEP];		/* from root */	
	/* TPM (g species) */
	double	bb_tpm_litter[ASTEP];		/* from litter */
	double	bb_tpm_leaf[ASTEP];			/* from leaf */
	double	bb_tpm_wood[ASTEP];			/* from wood */
	double	bb_tpm_root[ASTEP];			/* from root */
	/* TEC (g species) */
	double	bb_tec_litter[ASTEP];		/* from litter */
	double	bb_tec_leaf[ASTEP];			/* from leaf */
	double	bb_tec_wood[ASTEP];			/* from wood */
	double	bb_tec_root[ASTEP];			/* from root */
	
	/* Biogenic VOC, in micro g C m-2 month-1 */
	double	voc_isopr_g97[ASTEP];				/* isoprene */
	double	voc_monotrp_g97[ASTEP];				/* monoterpene */
	double	voc_methanl_g97[ASTEP];				/* methanol */
	double	voc_acetone_g97[ASTEP];				/* acetone */
	double	voc_actaldhd_g97[ASTEP];			/* acetoaldehyde */
	double	voc_frmardhd_g97[ASTEP];			/* formaldehyde */
	double	voc_formacd_g97[ASTEP];				/* formacid */
	double	voc_acetacd_g97[ASTEP];				/* acetoacid */
	double	voc_co_g97[ASTEP];					/* CO */
    /* added 2014/09/11 by A.Ito */
	double	voc_afarnesene[ASTEP];			/* alpha-Farnesene */
	double	voc_bcaryophyllene[ASTEP];		/* beta-Caryophyllene */
	double	voc_othersesqui[ASTEP];			/* other sesquiterpenes */
	
	/* stable carbon isotope composition, d13C, permille */
	double	d13c_nep[ASTEP];					/* NEP */
	double	d13c_ncb[ASTEP];					/* NCB */
	double	d13c_sr[ASTEP];						/* SR, soil respiration */
	double	d13c_er[ASTEP];						/* ER, ecosystem respiration */

	double	efflux_p;							/* total CO2 efflux */
	double	d13c_efflux_p;						/* d13C */
	
	/* d14C: added by A.Ito (2009/07/12) */
	double	d14c_sr[ASTEP];						/* d14C of soil respiration */
	double	d14c_er[ASTEP];						/* d14C of ecosystem respiration */
	
	/* erosion, Mg ha-1 yr-1 ******/
	/* total */
	double	erod_soil;							/* erosion of mineral soil */
	double	erod_orgmat;						/* erosion of organic matter */
	double	erod_carbon;						/* erosion of carbon */
	
	/* wood harvest: 2010/11/09 */
	double	hvst_wood;
};
