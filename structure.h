/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
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

/******* grid conditions, derived from data files *********************************/
struct Grid{ 
	/* time ***********/
	long 	mm[12];					/* number of days in each month */

	long 	phase;					/* simulation phase: 0-spinup, 1-past, 2-future */
	long 	y;						/* calculation time from the simulation onset, in year */
	long 	m;						/* month of the year, from Jan. to Dec., 0 to 11 */
	long 	time_hyd;				/* time to reach stabilization of water budget */
	long 	time;					/* time to reach stabilization of carbon budget */
	
	/* location ***********/
	long 	row;					/* grid order, row  in 0.5 degree grid */
	long 	col;					/* grid order, column in 0.5 degree grid */
	long 	gcm_row;				/* grid order, row in GCM's grid */
	long 	gcm_col;				/* grid order, column in GCM's grid */

	double 	lat;					/* latitude of the center of grid cell */
	double 	lon;					/* longitude of the center of grid cell */

	double 	area;					/* biome area for each cell */
	/* double 	cell_area;	*/		
	/* double 	land_area;	*/		/* total area and land area for each cell */

	long 	country;				/* country number */
	long 	region;					/* continental region number */
	long 	nnn;					/* cell numbers from the origin */
	
	/* Olson's actual biome number */
	long 	veg_olson;			
	/* SAGE's potential biome number */
	long 	veg_sage;	
	
	/* atmospheric condition ***********/
	long 	CO2y;					/* year for CO2 level estimation */
	double 	bCO2[12];				/* background CO2 concentration, in ppmv */
	double 	d13C_bCO2[12];			/* stable carbon isotope composition of background CO2, permille */

	/* climate condition: *[] means the transitional value */
	long 	climy;					/* year of climate data */
	double 	tmp_sfc[12];			/* ground surface temperature, degree Celcius */
	double 	tmp_sfc_am;				/* mean temperature, degree Celcius */
	double 	tmp_sfc_mx;				/* maximum temperature, degree Celcius */
	double 	tmp_sfc_mn;				/* minimum temperature, degree Celcius */
	double 	tmp_2m[12];				/* 2m air temperature, degree Celcius */
	double 	gp_atem;				/* average temperature during the growing period, degree Celcius */
	double 	gp_tem;					/* average temperature during the growing period, degree Celcius */
	double 	gp_pre;					/* precipitation during the prowing period, mm */
	double 	tmp10_soil[12];			/* soil temperature at 10 cm depth, degree Celcius */
	double 	tmp200_soil[12];		/* soil temperature at 200 cm depth, degree Celcius */
	double 	tcdc_clm[12];			/* total cloudiness, fraction */
	double 	prate_sfc[12];			/* precipitation, mm mon-1 */
	double 	prate_sfc_ann;			/* precipitation, mm mon-1 */
	double 	spfh_2m[12];			/* specific humidity, kg kg-1 */
	double 	wnd_10m[12];			/* wind velocity, m s-1 */

	/* climate condition: *_a[] means the average during 1965 to 1998 */	
	double 	tmp_sfc_a[12];			/* ground surface temperature, degree Celcius */
	double 	tmp_2m_a[12];			/* 2m air temperature, degree Celcius */
	double 	tmp10_soil_a[12];		/* soil temperature at 10 cm depth, degree Celcius */
	double 	tmp200_soil_a[12];		/* soil temperature at 200 cm depth, degree Celcius */
	double 	tcdc_clm_a[12];			/* total cloudiness, fraction */
	double 	prate_sfc_a[12];		/* precipitation, mm mon-1 */
	double 	spfh_2m_a[12];			/* specific humidity, kg kg-1 */
	double 	ugrd_10m_a[12];			/* zonal wind velocity, m s-1 */
	double	vgrd_10m_a[12];			/* meridional wind velocity, m s-1 */
	double 	rad_a[12];				/* solar radiation, W m-2 */
	double 	par_a[12];				/* PAR, micro mol m-2 s-1 */
	
	double 	prec_sub_a[12];			/* precipitation from substitute data (UEA/CRU), mm mon-1 */

	double 	topo;					/* topography, orology, and altitude, m above MSL */
	/* double 	whc;				*/	/* soil water holding capacity */
	/* double 	whc30;				*/	/* soil water holding capacity of above 300mm soil */
	double 	sd;						/* soil rooting depth */
	double 	hyd_cond;				/* hydraulic conductivity, */
	
	double	sl_dec[12];				/* solar declination, degree */
	double	sl_hgt[12];				/* solar hight at midday, degree */
	double	dlen[12];				/* monthly day length, hour */
	double	top_rad[12];				/* downward solar radiation at the atmosphere-top, W m-2 */
	double	gl_rad[12];				/* global radiation at the canopy-top, W m-2 */
	double	par[12];				/* monthly photosynthetically active radiation, micro mol photon m-2 s-1 */

	double	par_be[12];				/* PAR, beam, in energy W/m2 */
	double	par_de[12];				/* PAR, diffuse, in energy W/m2 */
	double	par_bp[12];				/* PAR, beam, in photon micro mol photon /m2 /s */
	double	par_dp[12];				/* PAR, diffuse, in photon micro mol photon /m2 /s */

	/* GCM climate: year x month x row x column */
	float	proj_tmp2m[131][12][160][320];			/* temperature */
	float	proj_prec[131][12][160][320];			/* precipitation */
	float	proj_shum[131][12][160][320];			/* specific humidity */
	float	proj_rad[131][12][160][320];			/* surface downward solar radiation */
	/* 1970-1999 average: month x row x column */
	float	proj_tmp2m_b[12][160][320];				/* temperature */
	float	proj_prec_b[12][160][320];				/* precipitation */
	float	proj_shum_b[12][160][320];				/* specific humidity */
	float	proj_rad_b[12][160][320];				/* surface downward solar radiation */

	double	proj_prec_co;			/* carry-over of negative precipitation */

	/* UEA/CRU TS2.0 data */
	long	cru_exist;							/* flag of CRU TS2.1 data availability */
	double	hist_tmp[CRU_PD+2][12];				/* temperature */
	double	hist_pre[CRU_PD+2][12];				/* precipitation */
	double	hist_cld[CRU_PD+2][12];				/* cloud cover */
	double	hist_vap[CRU_PD+2][12];				/* vapor pressure */
	/* historical average */
	double	hist_tmp_b[12];						/* temperature */
	double	hist_pre_b[12];						/* precipitation */
	double	hist_cld_b[12];						/* cloud cover */
	double	hist_vap_b[12];						/* vapor pressure */
	
	/* erosion */
	long 	rvbasin;							/* ID of river basin */
	double 	albedo_soil;						/* soil albedo */
	double 	pcnt_orgmat;						/* percent organic matter */
	double 	fls_slope;							/* slope factors LS for RUSLE erosion model */
	double 	fk_edodibility;						/* soil erodibility for RUSLE erosion model */
	
	/* land-use change */
	double 	fcrop_sage[293];				/* historical data, 1700-1992, Ramankutty & Foley */
	/* future, 1990-2100, IMAGE2, Wang et al. */
	double 	fcrop3_image[111];				/* C3 crop */
	double 	fcrop4_image[111];				/* C4 crop */
	double 	fgrass3_image[111];				/* C3 pasture grass */
	double 	fgrass4_image[111];				/* C4 pasture grass */
	
	/* EOS-WEBSTER, 1700-2000, Hurtt et al. */
	double	fcrop_eossagehyde[301];		/* cropland fraction */
	double	fpast_eossagehyde[301];		/* pasture fraction */
	double	fprim_eossagehyde[301];		/* primary land fraction */
	double	fsecd_eossagehyde[301];		/* secondary land fraction */
	double	ssma_eossagehyde[301];		/* secondary land property 1 */
	double	ssmb_eossagehyde[301];		/* secondary land property 2 */
	double	t_cp_eossagehyde[301];		/* conversion crop to pasture */
	double	t_cs_eossagehyde[301];		/* conversion crop to secondary */
	double	t_pc_eossagehyde[301];		/* conversion pasture to crop */
	double	t_ps_eossagehyde[301];		/* conversion pasture to secondary */
	double	t_sc_eossagehyde[301];		/* conversion secondary to crop */
	double	t_sp_eossagehyde[301];		/* conversion secondary to pasture */
	double	t_ss1_eossagehyde[301];		/*  */
	double	t_ss2_eossagehyde[301];		/*  */
	double	t_ss3_eossagehyde[301];		/*  */
	double	t_vc_eossagehyde[301];		/* conversion primary to crop */
	double	t_vp_eossagehyde[301];		/* conversion primary to pasture */
	double	t_vs1_eossagehyde[301];		/*  */
	double	t_vs2_eossagehyde[301];		/*  */

	double 	f_crop_con;					/* contemporary cropland fraction */
	double 	f_crop_p;					/* previous cropland fraction */
	double 	f_crop_trend;				/* trend of cropland conversion, fraction/yr */
	double 	f_pasture_con;				/* contemporary pasture fraction */
	double 	f_pasture_p;				/* previous pasture fraction */
	double 	f_pasture_trend;			/* trend of pasture conversion, fraction/yr */
	
	double	f_deforest;					/* deforestation rate */
	double	f_deforest_v;				/* in primary lands */
	double	f_deforest_s;				/* in secondary lands */
	
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
	double 	f_wetland;					/* wetland fraction */
	double 	f_upland;					/* upland (e.g. forest, grassland) fraction */
	
	/* nitrogen deposition */
	double 	ndepo[3];					/* N deposition by Galloway et al. (2004) */
	double	nfert_nh4;					/* NH4+ fertilization input */
	double	nfert_no3;					/* NO3- fertilization input */
};			

/******* grid conditions, derived from submodules *************************************/
struct Loct{ 
	double	aCO2[12];				 	/* ambient CO2 concentration, in ppmv */
	double	d13C_aCO2[12];				/* stable carbon isotope composition of CO2, dimensionless */
	double	cnpy_co2_recyc;				/* within-canopy CO2 recycling ratio */

	double	C4ptn[12];					/* ground coverage of C3 plants, fraction */
	double	C3ptn[12];					/* ground coverage of C4 plants, fraction */
	long	gd[12], bbm;				/* vegetative growing period, days */
	double	gdd[12];				 	/* cumulative growth degree days, degC days */

	double	albedo_sfc[12];				/* land-surface albedo */
	double	gl_rad_g[12];				/* global radiation under the canopy, W m-2 */
	double	rad_net_p[12];				/* net radiation, canopy, W m-2 */
	double	rad_net_g[12];				/* net radiation, soil surface, W m-2 */
	double	rad_net_long[12];			/* net long-wave radiation, W m-2 */
	double	rad_net_short[12];			/* net short-wave radiation, W m-2 */
	double	rad_net[12];				/* net radiation, W m-2 */
	double	rdi;						/* radiative dryness index by Budyko */
	
	double	fapar_mono[12];				

	double	apar_bp[12];				/* absorbed PAR photon, beam, micro mol photon m-2 s-1 */
	double	apar_dp[12];				/* absorbed PAR photon, diffuse, micro mol photon m-2 s-1 */
	double	pet_prty[12];				/* Priestley-Taylor potential evapotranspiration, mm month-1 */
	double	pet_prty_ann;				/* annual Priestley-Taylor potential evapotranspiration, mm yr-1 */
	
	double	prsr[12];				/* air pressure, hPa */
	double	dnsa[12];				/* density of air, kg m-3 */

	double	vp[12];					/* vapour pressure, hPa */
	double	vps[12];				/* saturation vapour pressure, hPa */
	double	slope_vps[12];			/* slope of saturation vapour pressure related to tempertaure, hPa deg C-1 */
	double	r_aero[12];				/* aerodynamic resistance, s m-1 */
	double	vpd[12];				/* vapour pressure deficit, hPa */
	
	double	lai[12];				/* leaf area index, m2 m-2 */
	double	canopy_con[12];			/* canopy conductance, mmol H2O m-2 s-1 */
	double	f_vegcov[12];			/* fractional vegetation cover */
	
	/* water pools, mm (= kg/m2) */
	double	snwa;					/* water equivalent snow depth, mm */
	double	msnwa[12];				/* monthly */
	double	sw30;					/* soil water content of above 30cm soil, mm */
	double	msw30[12];				/* monthly */
	double	sww;					/* whole soil water content, mm */
	double	msww[12];				/* monthly */

	/* water fluxes, mm / month */
	double	pm_evp[12];				/* potential soil evaporation rate, mm */
	double	pm_trn[12];				/* potential transpiration rate, mm */
	double	pm_incep[12];			/* potential intercepted-water evaporation rate, mm */
	double	incep[12];			
	double	evpr[12];				/* actual evaporation rate, mm */
	double	trspr[12];				/* actual transpiration rate, mm */
	double	ro1[12];				/* runoff from upper soil water, mm */
	double	ro2[12];				/* runoff from lower soil water, mm */
	double	pntrt[12];				/* water penetration from upper to lower layer, mm */
	double	thaw[12];				/* snow thaw water, mm */
	double	vmc30[12];				/* volumatric moisture content of upper layer, fraction */
	double	vmc[12];				/* volumatric moisture content of lower layer, fraction */
	double	snp[12];				/* snow fraction of precipitation */
	double	soil_appr30;			/* soil aperture of upper layer, fraction */
	double	soil_apprw;				/* soil aperture of lower layer, fraction */
	
	double	n_frtlz_in;				/* N-fertilization input */
	double	depo_no3[12];			/* NO3- deposition */
	double	depo_nh4[12];			/* NH4+ deposition */
	
	/* CASA moisture **********************************/
	/* long	mday;			*/
	double	m_m[12];				/* soil moisture index */
	double	m_m_pre;				/* precipitation */
	double	m_rdr[12];				/* dryness */
	double	m_vmc[12];				/* volumetric soil moisture */
	double	m_pet[12];				/* potential evapotranspiration */
	double	m_e[12];				/* soil water change */
	double	m_sw[12];				/* soil wetness */
	double	i_w[12];				/* moisture scalar */
	/* double	d_tmp[31];			
	double	d_sw[31];			
	double	d_vmc[31];			*/
	double	wfps[12];				/* water-filled pore space */
	
	/* monitoring variables for debugging */
	double	xx1[12];
	double	xx2[12];
	double	xx3[12];
	double	xx4[12];
	double	xx5[12];
};			

/******* vegetation characteristics *******/
struct Pchar{ 
	double	albedo;					/* reflectivity, or albedo */
	double	fapar[12];				/* fraction of absorbed PAR */

	/* allocation *********/
	double	opt_lai[12];			/* optimum leaf area index */
	double	alloc_ass;				/* allocation ration for assimilation organ, fraction */
	double	alloc_abg;				/* allocation ration for aboveground non-assimilation organ, fraction */
	double	malloc_f[12];			/* monthly allocation ratio to leaf */
	double	malloc_c[12];			/* monthly allocation ratio to stem */
	double	malloc_r[12];			/* monthly allocation ratio to root */
	
	/* N parameters *******/
	double	n_salvage;				/* N re-uptake ratio */
	double	cn0_strg;				/* base C/N of storage nitrogen */
	double	n_conc_larea;			/* leaf N concentration, mmol N m-2 leaf area */
	double	n_conc_larea_m[12];		/* monthly */
	double	n_conc_lmass;			/* leaf N concentration, mmol N g-1 leaf weight */
	double	kn_nphoto;				/* canopy N gradient, exponential coefficient */
	double	amax_nphoto;			/* max. photosynthesis rate, micro mol CO2 m-2 s-1 */
	double	rd_n;					/* dark respiration, micro mol CO2 m-2 s-1 */
	
	/* phenology *********/
	double	gdd;					/* growing degree days */
	double	mgdd[12];				/* monthly growing degree days */
	double	grw_pd;					/* growing period */
	long 	season[12];				/* phenological stage as a function of season */
	 	/* 0: dormancy */
	 	/* 1: vegetative growth */
	 	/* 2: emergence of new leaf */
	 	/* 3: abandon of old leaf */
	long 	frag_emg;				/* flag of leaf emergence */
	double	frag_dcd;				/* flag of deciduou */
	
	double	fleaf_age[49];			/* leaf age distribution, 0(new) to 48-over(old) months */
	
	/* photosynthesis *******/
	long	phototype;				/** photosynthetic metabolic pathway, 3=C3, 4=C4, 5=CAM **/
	/*  veg->psat[grid->m] = veg->pmax*ftem*fstl*fnstl   */
	double	psat[12];				/* light-saturated rate, micro mol CO2 m-2 s-1 */
	double	pmax;					/* potential maximum rate, micro mol CO2 m-2 s-1 */
	double	ft[12];					/* temperature coefficient */
	double	fcd[12];				/* CO2 coefficient */
	double	fsw[12];				/* soil water coefficient */
	double	ptop;					/* canopy-top photosynthetic rate */
	double	sla;					/* specific leaf area, cm2 g dm-1 */
	double	eK0;					/* light attenuation coefficient, no dimension */
	double	eK[12];					/* light attenuation coefficient, no dimension */
	double	lue0;					/* control light dependence coefficient, mol CO2 mol photon-1 */
	double	lue[12];				/* monthly quantum yield , mol CO2 mol photon-1 */ 
	double	topt;					/* optimum temperature, deg C */
	double	topt0;					/* optimum temperature (co2 dependent), deg C */
	double	tmin;					/* minimum temperature, deg C */
	double	tmax;					/* maximum temperature, deg C */
	double	ci[12];					/* monthly intercellular CO2 concentration, ppmv */
	double	kmci;					/* dependence of photosynthesis on intercellular CO2 concentration, ppmv */
	double	cmpcd0;					/* CO2 compensation point, ppmv */
	double	cmpcd[12];				/* CO2 compensation point, ppmv */
	/* parameters of conductance */
	double	gs[12];					/* monthly stomatal conductance, mmol H2O m-2 s-1 */
	double	gc[12];					/* monthly canopy conductance, mmol H2O m-2 s-1 */
	double	gs_b0;					/* Leuninig stomata model parameters */
	double	gs_b1;					/* Leuninig stomata model parameters */
	double	gs_b2;					/* Leuninig stomata model parameters */
	double	km_nstl;				/* maximum stomatal conductance */
	
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
	double	qTf[12];				/* leaf */
	double	qTc[12];				/* stem */
	double	qTr[12];				/* root */
	/* temperature dependence Q10 at 15 deg C *****/
	double	qTf0;					/* leaf */
	double	qTc0;					/* stem */
	double	qTr0;					/* root */
	
	/* litter fall ********/
	/* specific litter fall rate, fraction */
	double	lf[12];					/* leaf */
	double	lc[12];					/* stem */
	double	lr[12];					/* root */
	/* base rate *********/
	double	lf0;					/* leaf */
	double	lc0;					/* stem */
	double	lr0;					/* root */
	/* deciduous leaf fraction *****/
	double	dcd;				
	
	/* photosynthetgic 13C fractionation, permil *******/
	double	photo_13c_frac[12];
		
	/* specific constants *********************************/
	/* root stratification *******/
	double	root_strat;			/* root profile parameter */
};			

/******* soil characteristics *******/
struct Schar{ 
	double	albedo0;			/* reflectivity, or albedo */
	double	albedo[12];			/* monthly */
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
	double	ft_l[12];			/* temperature for litter */
	double	ft_h[12];			/* temperature for humus */
	double	fm_l[12];			/* moisture for litter */
	double	fm_h[12];			/* moisture for humus */
	/* with "0", average value for sensitivity analysis */
	double	ft0_l[12];			
	double	ft0_h[12];
	double	fm0_l[12];
	double	fm0_h[12];			
};			

/******* ecosystem characteristics *******/
struct Echar{ 
	short	v_type;				/* vegetation classification types */
								/* 1: Olson actual vegetation */
								/* 2: SAGE natural vegetation */

	struct	Pchar c3;			/* for C3 plants */
	struct	Pchar c4;			/* for C4 plants */
	struct	Schar soil;			/* for soil organic matter */
};			

/* plant biomass *******************************/
struct Pmas{ 
	double	lai[12];			/* monthly leaf area index(LAI), m2 m-2 */
	/* carbon mass, Mg C ha-1 *****/ 
	double	fol;				/* leaf */
	double	mfol[12];			/* monthly */
	double	stm;				/* stem and branch */
	double	mstm[12];			/* monthly */
	double	rot;				/* root */
	double	mrot[12];			/* monthly */
	
	double	plant[12];			/* monthly plant mass */
	double	fol_p;				/* previous leaf carbon */
	
	/* stable carbon isotope composition, permille */
	double	d13c_fol;			/* leaf */
	double	d13c_mfol[12];		/* monthly */
	double	d13c_stm;			/* stem */
	double	d13c_mstm[12];		/* monthly */
	double	d13c_rot;			/* root */
	double	d13c_mrot[12];		/* monthly */
	
	double	d13c_plant[12];		/* total plant */

	/* plant N, g N ha-1 */
	double	n_cnpy;				/* N in canopy */
	double	n_cnpy_m[12];		/* monthly */	
	double	n_strg;				/* N in storage pool */
	double	n_strg_m[12];		/* monthly */
};			

/* soil carbon storage *******************************/
struct Smas{ 
	/* carbon mass, Mg C ha-1 */ 
	double	ltr;				/* litter mass */
	double	ltr_m[12];			/* monthly litter mass */
	double	msl;				/* mineral soil and humus mass */
	double	msl_m[12];			/* monthly mineral soil and humus mass */
	double	soil[12];			/* monthly soil mass */
	
	double	doc;				/* dissolved organic carbon discharge */
	double	doc_m[12];			/* monthly */
	
	/* stable carbon isotope composition, permille */
	double	d13c_ltr;			/* litter */
	double	d13c_ltr_m[12];		/* monthly */
	double	d13c_msl;			/* mineral soil */
	double	d13c_msl_m[12];		/* monthly */
	double	d13c_soil[12];		/* total soil */

	/* soil inorganic N, g N ha-1 */
	double	n_no3;				/* NO3- */
	double	n_no3_m[12];		/* monthly */
	double	n_nh4;				/* NH4+ */
	double	n_nh4_m[12];		/* monthly */
	
	/* soil organic N, g N ha-1 */
	double	n_mcrb;				/* microbe */
	double	n_mcrb_m[12];		/* monthly */
	double	n_lttr;				/* soil organic litter */
	double	n_lttr_m[12];		/* monthly */
	double	n_hums;				/* soil organic humus */
	double	n_hums_m[12];		/* monthly */
};			

/* ecosystem carbon storage *************************/
struct Mass{ 
	struct	Pmas c3;			/* C3 plant mass */
	struct	Pmas c4;			/* C4 plant mass */
	struct	Pmas plant;			/* all plant mass */
	struct	Smas soil;			/* soil mass */

	/* carbon mass */ 
	double	total[12];			/* ecosystem total carbon storage, Mg C ha-1 */

	double	lai_p;				/* previous LAI, m2 m-2 */
	
	/* stable carbon isotope composition */
	double	d13c_total[12];		/* total d13C, permille */
};			

/* plant carbon fluxes, all monthly  *******************/
struct Pflx{ 
	/* carbon flux, in Mg C ha-1 mon-1 */ 
	double	gpp[12];			/* gross primary production */
	double	spp[12];			/* net primary production */
	double	epp[12];			/* net primary production */
	double	npp[12];			/* net primary production */
	
	double	tpf[12];			/* translocation of photosynthate to foliage */
	double	tpc[12];			/* translocation of photosynthate to stem */
	double	tpr[12];			/* translocation of photosynthate to root */
	double	tpp[12];			/* translocation of photosynthate to root */
 	
	double	rp[12];				/* plant respiration, =rpm+rpg */
	double	rpg[12];			/* plant growth respiration */
	double	rpm[12];			/* plant maintenance respiration */
	double	rfg[12];			/* foliage growth respiration */
	double	rfm[12];			/* foliage maintenance respiration */
	double	rcg[12];			/* stem and branch growth respiration */
	double	rcm[12];			/* stem and branch maintenance respiration */
	double	rrg[12];			/* root growth respiration */
	double	rrm[12];			/* root maintenance respiration */
	
	double	lf[12];				/* foliage litterfall */
	double	lc[12];				/* stem and branch litterfall */
	double	lr[12];				/* root litterfall */
	double	lL[12];				/* total litterfall */
	double	lf_c[12];			/* leaf shedding in C3/C4 altyeration in grassland */
	
	double	hvst[12];			/* harvest of crops */
	
	double	emit_ch4_kirschbaum_mass[12];		/* plant CH4 emission, mass-based */
	double	emit_ch4_kirschbaum_photo[12];		/* plant CH4 emission, photosynthesis-based */
	
	/* stable carbon isotope composition, d13C, permille */
	double	d13c_gpp[12];		/* GPP */
	double	d13c_spp[12];		/* SPP */
	double	d13c_epp[12];		/* EPP */
	double	d13c_npp[12];		/* NPP */
	
	double	d13c_tpf[12];			
	double	d13c_tpc[12];			
	double	d13c_tpr[12];			
	double	d13c_tpp[12];			
 	
	double	d13c_rp[12];			
	double	d13c_rpg[12];			
	double	d13c_rpm[12];			
	double	d13c_rfg[12];			
	double	d13c_rfm[12];			
	double	d13c_rcg[12];			
	double	d13c_rcm[12];
	double	d13c_rrg[12];			
	double	d13c_rrm[12];			
	
	double	d13c_lf[12];			
	double	d13c_lc[12];			
	double	d13c_lr[12];			
	double	d13c_lL[12];
	double	d13c_lf_c[12];			

	double	d13c_hvst[12];			/* harvest */

	/* nitrogen flux, g N ha-1 mon-1 */
	double	n_biofix[12];			/* biological N2 fixation */
	double	uptake_no3[12];			/* root NO3- uptake */
	double	uptake_nh4[12];			/* root NH4+ uptake */
	double	n_alloc_cnpy[12];		/* N allocation to canopy */
	double	n_alloc_strg[12];		/* N allocation to storage */
	double	n_salvage[12];			/* N salvage from shedding leaves */
	double	n_realloc[12];			/* N reallocation from storage to leaves */
	double	n_abdn_cnpy[12];		/* N abandonment from canopy */
	double	n_abdn_strg[12];		/* N abandonment from storage */
};			

/******* soil carbon fluxes, all in Mg C ha-1 mon-1 *******/
struct Sflx{	 
	/* carbon flux */ 
	double	lL[12];			/* litter input */
	double	rl[12];			/* litter decomposition */
	double	rh[12];			/* mineral soil and humus decomposition */
	double	rS[12];			/* total decomposition */
	double	sf[12];			/* humus formation */
	
	/* dissolved organic carbon (DOC) discharge */
	double	doc_boyer[12];

	/* stable carbon isotope composition */
	double	d13c_lL[12];			/* litter input */
	double	d13c_rl[12];			/* litter decomposition */
	double	d13c_rh[12];			/* humus decomposition */
	double	d13c_rS[12];			/* heterotrophic respiration */
	double	d13c_sf[12];			/* soil formation from litter to humus */
	
	/* CH4 oxydation ********/
	double	ch4oxy_ridg[12];				/* Ridgwell scheme */
	double	ch4oxy_casa[12];				/* CASA-Potter scheme */
	double	ch4oxy_delgrosso[12];			/* DayCnet-Del Grosso scheme */
	double	ch4oxy_curry[12];				/* Curry scheme */
	
	/* CH4 emission **********/
	double	ch4emit_wetland_cao[12];		/* CH4 emission, Cao scheme, from wetland */
	double	ch4emit_paddy_cao[12];			/* CH4 emission, Cao scheme, from paddy field */
	
	/* N flow, g N ha-1 mon-1 **/
	double	n_abdn[12];						/* N abandonment as litter */
	
	double	d_n2o_ntr_ngas[12];				/* N2O from nitrification by NGAS */
	double	d_n2o_dnt_ngas[12];				/* N2O from denitrification by NGAS */
	double	d_n2o_ngas[12];					/* total N2O by NGAS */
	double	d_n2_ngas[12];					/* N2 by NGAS */
	
	double	d_no_casa[12];					/* NO by CASA */
	double	d_n2_casa[12];					/* N2 by CASA */
	double	d_n2o_casa[12];					/* N2O by CASA */
	
	double	n_nh3vlt[12];					/* NH3 volatilization */
	double	n_leach[12];					/* NO3- leaching */
	double	n_minerlz_lttr[12];				/* litter N mineralization */
	double	n_minerlz_hums[12];				/* humus N mineralization */
	double	n_nitrif[12];					/* nitrification in NGAS */
	double	n_immbl[12];					/* N immobilization */
	double	n_mcrb_abdn[12];				/* microbial abandonment */
};

/******* ecosystem carbon fluxes *******/
struct Flux{ 
	struct	Pflx c3;				/* C3 plant fluxes */
	struct	Pflx c4;				/* C4 plant fluxes */
	struct	Pflx plant;				/* total plant fluxes */
	struct	Sflx soil;				/* soil fluxes */

	double	npp_miami;				/* NPP estimated by Miami model, Mg C ha-1 yr-1*/
	double	npp_montreal;			/* NPP estimated by Montreal model, Mg C ha-1 yr-1*/
	double	npp_rosenzweig;			/* NPP estimated by Montreal model, Mg C ha-1 yr-1*/
	double	npp_schuur;				/* NPP estimated by Montreal model, Mg C ha-1 yr-1*/

	/* carbon flux */ 
	double	nep[12];				/* net ecosystem production, Mg C ha-1 mon-1 */
	double	ncb[12];				/* net carbon balance of grid, Mg C ha-1 mon-1 */
	double	lL0[12];				/* total litter fall */
	double	sresp[12];				/* soil respiration */
	
	/* land-use change */
	double	lu_detr;				/* LUC emission from detritus */
	double	lu_conv;				/* LUC emission from 1-year pool */
	double	lu_ten;					/* LUC emission from 10-year pool */
	double	lu_hund;				/* LUC emission from 100-year pool */
	/* historical land-use-generated detritus production */
	double	detr_ten[10];			/* 10-year pool */
	double	detr_hund[100];			/* 100-year pool */
	
	/* biomass burning */
	double	f_burnt;				/* burnt fraction */
	double	day_fire[12];			/* days of fire */
	double	a_burnt[12];			/* area burnt */
	/* CO2 (g species) */
	double	bb_co2_litter[12];			/* from litter */
	double	bb_co2_leaf[12];			/* from leaf */
	double	bb_co2_wood[12];			/* from wood */
	double	bb_co2_root[12];			/* from root */
	/* CO (g species) */
	double	bb_co_litter[12];			
	double	bb_co_leaf[12];			
	double	bb_co_wood[12];			
	double	bb_co_root[12];			
	/* CH4 (g species) */
	double	bb_ch4_litter[12];			
	double	bb_ch4_leaf[12];			
	double	bb_ch4_wood[12];			
	double	bb_ch4_root[12];			
	/* NMHC (g species) */
	double	bb_nmhc_litter[12];			
	double	bb_nmhc_leaf[12];			
	double	bb_nmhc_wood[12];			
	double	bb_nmhc_root[12];			
	/* OC (g species) */
	double	bb_oc_litter[12];			
	double	bb_oc_leaf[12];			
	double	bb_oc_wood[12];			
	double	bb_oc_root[12];			
	/* BC (g species) */
	double	bb_bc_litter[12];			
	double	bb_bc_leaf[12];			
	double	bb_bc_wood[12];			
	double	bb_bc_root[12];			
	/* NOx (g species) */
	double	bb_nox_litter[12];			
	double	bb_nox_leaf[12];			
	double	bb_nox_wood[12];			
	double	bb_nox_root[12];			
	/* SO2 (g species) */
	double	bb_so2_litter[12];			
	double	bb_so2_leaf[12];			
	double	bb_so2_wood[12];			
	double	bb_so2_root[12];			
	/* PM2.5 (g species) */
	double	bb_pm25_litter[12];			
	double	bb_pm25_leaf[12];			
	double	bb_pm25_wood[12];			
	double	bb_pm25_root[12];			
	/* TPM (g species) */
	double	bb_tpm_litter[12];			
	double	bb_tpm_leaf[12];			
	double	bb_tpm_wood[12];			
	double	bb_tpm_root[12];			
	/* TEC (g species) */
	double	bb_tec_litter[12];			
	double	bb_tec_leaf[12];			
	double	bb_tec_wood[12];			
	double	bb_tec_root[12];			
	
	/* VOC */
	double	voc_isopr_g97[12];					/* isoprene */
	double	voc_monotrp_g97[12];				/* monoterpene */
	double	voc_methanl_g97[12];				/* methanol */
	double	voc_acetone_g97[12];				/* acetone */
	double	voc_actaldhd_g97[12];				/* acetoaldehyde */
	double	voc_frmardhd_g97[12];				/* formaldehyde */
	double	voc_formacd_g97[12];				/* formacid */
	double	voc_acetacd_g97[12];				/* acetoacid */
	double	voc_co_g97[12];						/* CO */
	
	/* stable carbon isotope composition, d13C, permille */
	double	d13c_nep[12];						/* NEP */
	double	d13c_ncb[12];						/* NCB */

	double	efflux_p;							/* total CO2 efflux */
	double	d13c_efflux_p;						/* d13C */
	
	/* erosion, Mg ha-1 yr-1 ******/
	/* total */
	double	erod_soil;							/* erosion of mineral soil */
	double	erod_orgmat;						/* erosion of organic matter */
	double	erod_carbon;						/* erosion of carbon */
	/* cropland */
	double	erod_soil_crop;						/* cropland erosion of mineral soil */
	double	erod_orgmat_crop;					/* cropland erosion of organic matter */
	double	erod_carbon_crop;					/* cropland erosion of carbon */
};
