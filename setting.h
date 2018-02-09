/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in November 29, 2007						*/

/* constants declaration */
#define dTr 0.0174533 /* angle conversion, from degree to radian */
#define rTd 57.29577951 /* angle conversion, from radian to degree */
#define PI 3.141592653 /** pai **/
#define cdTc 0.272727 /* from CO2-base to Crabon-base */
#define dmTc 2.2 /* from dry-matter-base to Carbon base */
#define cTdm 0.4545 /* from dry-matter-base to Carbon base */
#define lTs (3600.0*12.0/100000000.0) /* from micro-mol m-2 s-1 to Mg C ha-1 day-1 */
#define ZAT 273.15 /* zero degree centigrade in absolute temperature */
#define STCIR (0.0111/0.9889) /* standard stable carbon isotope ratio */
#define UGC 8.314 /* universal gas constant */
#define SLC 1367.0 /* solar constant */
#define SBC (5.6703 / 100000000.0) /* Stephan-Boltzman Constant, W m-2 K-4 */
#define GAC 9.8 /* gravity acceleration constant, m s-2 */

/***********************************************************/
/* year time-step, 12=monthly */
#define ASTEP 12
/* annual daily-step, 24=hourly */
#define DSTEP 24

/* grid configulation */
#define N_ROW 360
#define N_COL 720

/* number of geographical regions */
#define N_REG 23
/* defined in region_giorgi() in vegetdeal.c */

/* initial (minimal) carbon stock ***********/
#define INT_C 0.01

#define CALC_STEP 1
#define CALC_OFFET 0
/* 1: every grid */
/* 10: every 10 grid */

/***********************************************************/
#define ISIMIP_RUN 0
/* 0: normal (no ISI-MIP) */
/* 1: ISI-MIP 1st-phase runs + CD-LINKS (2016/11/17 by A.Ito ) */
/* 2: PLUME (ISI-MIP Phase 2) runs : 2014/07/31 by A.Ito */
/* 3: ISI-MIP 2a historical runs */
/* 4: ISI-MIP 2b 1.5 or 2.0 deg target */
#define ISIMIP2_FIXCD 0
/* ISI-MIP2b IIb (fixed CO2 after 2005) */

#define GEOMIP_RUN 0
/* 0: normal (no GEO-MIP) */
/* 1: GEO-MIP runs */

#define IMPRESSIONS_RUN 0
/* 0: off */
/* 1: sensitivity run (for IRS) */
/* 2: sensitivity run (for IRS) + seasonal change + RCP4.5-CO2 */
/* 3: sensitivity run (for IRS) + seasonal change + RCP8.5-CO2 */

/* S10-BECCS experiment: 2016/02/15 by A.Ito */
#define EX_BECCS 0
/* 0: off (inc. scenario S5) */
/* 1: on (scenario S3) */
/* 2: on (scenario S3) based on Kinoshita-san data 2017/02 */
/* 3: on (scenario S3) based on Kinoshita-san data 2017/10 */
/* sub-scenarios: 2017/10/31 by A.Ito */
#define EX_BECCS_SUB 0
/* 0: no sub-setting */
/* 1: no beccs */
/* 2: fix land-use after 2000 */

/* biofuel experiment: 2015/08/21 by A.Ito */
#define BIOFUEL_RUN 0
/* 0: off */
/* 1: current fertilization */
/* 2: low fertilization */
/* 3: middle fertilization */
#define DL_BF 91 /* biofuel data length */

/* NMIP: N2O model intercomparison runs */
#define NMIP_RUN 0
/* 0: off */
/* 1: climate + CO2 + LCLUC + Ndep + Nfer + manure (S1) */
/* 2: climate + CO2 + LCLUC + Ndep + Nfer (S2) */
/* 3: climate + CO2 + LCLUC + Ndep (S3) */
/* 4: climate + CO2 + LCLUC (S4) */
/* 5: climate + CO2 (S5) */
/* 6: climate (S6) */
/* 7: all 1860 (S0) */

/* old NMIP ID */
/* 1: on use NMIP data (S1) */
/* 2: all fix (S0) */
/* 3: all combined (S2) */
/* 4: climate only (S3) */
/* 5: temperatute only (S3a) */
/* 6: precipitation only (S3b) */
/* 7: radiation only (S3c) */
/* 8: CO2 only (S4) */
/* 9: N deposition only (S5) */
/* 10: N fertilizer only (S6) */
/* 11: land-conversion only (S7) */
/* 12: land-conversion + N fertilizer + irrigation (S8) */

#if ISIMIP_RUN==4
    /* ISI-MIP2b */
    #define DL_NINPUT 639
    #define FDY_NINY 1661
#else
    /* #define DL_NINPUT 156 */
    #define DL_NINPUT 157 /* updated: 2017/10/19 by A.Ito */
    #define FDY_NINY 1860
#endif

/***********************************************************/
/* output text files */
#define OUTPUT_CARBON1 1
#define OUTPUT_CARBON2 1
#define OUTPUT_ISOTOPE 1
#define OUTPUT_NITROGEN 1
#define OUTPUT_HYDMET 1
#define OUTPUT_EROSION 1
#define OUTPUT_GHG 1
#define OUTPUT_BB 1
#define OUTPUT_BVOC 1
/* output binary */
#define C13_GOUT 0
#define C14_GOUT 0
#define PHYS_GOUT 1

/***********************************************************/
/* total vegetation number */
#define NVEG_OLSON 34	/* Olson veg (modified) */
#define NVEG_SAGE 16		/* SAGE veg (modified) */
/* #define NVEG_CROP 3	*/	/* crop types */
#define NVEG_CROP 4		/* crop types: add biofuel crop: 2015/08/21 by A.Ito */

/* calculation for land covers */
#define CALC_OLSON 1    /* matural vegetation */
#define CALC_CROP 1     /* cropland */
/* 0:off 1:on */

/* Olson croplands replaced by SAGE natural vegetation */
#define REPLACE_OLSON_CROP 1
/* 0:off 1:on */

/* change crop types: 2015/04/24 by A.Ito */
#define EX_CROP 0
/* 0: off (control) */
/* 1: all C3 (e.g. wheat) */
/* 2: all rice */
/* 3: all C4 (e.g. maize) */

/* atmopsheric GHG data length */
#if ISIMIP_RUN==1
    /* ISI-MIP: 2012/06/27 by A.Ito */
    #define DL_AGHG 736
    #define FDY_AGHG 1765
#elif ISIMIP_RUN==2
    /* PLUME: 2014/07/31 by A.Ito */
    #define DL_AGHG 736
    #define FDY_AGHG 1765
#elif ISIMIP_RUN==3
    /* ISI-MIP2a: 2014/12/02 by A.Ito */
    #define DL_AGHG 736
    #define FDY_AGHG 1765
#elif ISIMIP_RUN==4
    /* ISI-MIP2b: 2016/12/22 by A.Ito */
    #define DL_AGHG 639
    #define FDY_AGHG 1661
#elif GEOMIP_RUN==1
    /* GEO-MIP: 2012/06/27 by A.Ito */
    #define DL_AGHG 736
    #define FDY_AGHG 1765
#elif IMPRESSIONS_RUN==2
    /* IMPRESSIONS: 2017/05/02 by A.Ito */
    #define DL_AGHG 736
    #define FDY_AGHG 1765
#elif IMPRESSIONS_RUN==3
    /* IMPRESSIONS: 2017/05/02 by A.Ito */
    #define DL_AGHG 736
    #define FDY_AGHG 1765
#else
    #define DL_AGHG 553 /* default */
    #define FDY_AGHG 1750
#endif

/*****************************************************************/
/* simulation framework duration (years) */
#if ISIMIP_RUN==4
    #define PD_SIM 640 /* */ /* ISI-MIP2b (1.5/2.0 deg): 2016/12/22 by A.Ito */
#else
    #define PD_SIM 201
#endif
/* for memory preparation; not always actual experimental length */

#if ISIMIP_RUN==1
    #define FSY_HIST 1950 /* ISI-MIP */
    #define LSY_HIST 2099
#elif ISIMIP_RUN==4
    #define FSY_HIST 1661 /* ISI-MIP2b */
    #define LSY_HIST 2299
    /* #define LSY_HIST 2099 */ /* ISI-MIP2b (rcp6.0) */
#elif IMPRESSIONS_RUN==2
    #define FSY_HIST 1901 /* IMPRESSIONS 2 */
    #define LSY_HIST 2015
#elif IMPRESSIONS_RUN==3
    #define FSY_HIST 1901 /* IMPRESSIONS 2 */
    #define LSY_HIST 2015
#elif NMIP_RUN>=1
    #define FSY_HIST 1861 /* NMIP */
    #define LSY_HIST 2016 /* NMIP */
#else
    #define FSY_HIST 1901 /* */
    /* #define LSY_HIST 2016 */ /* history */
    #define LSY_HIST 2017 /* history */
#endif

/* start year (AD) of CO2 time series */
/* #define BGY_CO2Y 1901 */
/* cru-init 1901 */
/* con 1990 */
/* dcd 2081 */
/* NCEP1 1948 */
#if ISIMIP_RUN==1
    #define BGY_CO2Y 1950  /* ISI-MIP: 2012/06/27 by A.Ito */
#elif ISIMIP_RUN==2
    #define BGY_CO2Y 1901  /* PLUME: 2014/07/31 by A.Ito */
#elif ISIMIP_RUN==3
    #define BGY_CO2Y 1901  /* ISI-MIP2 (historical): 2014/11/30 by A.Ito */
#elif ISIMIP_RUN==4
    #define BGY_CO2Y 1661  /* ISI-MIP2b (1.5/2.0 deg): 2016/12/22 by A.Ito */
#elif IMPRESSIONS_RUN==2
    #define BGY_CO2Y 1901
#elif IMPRESSIONS_RUN==3
    #define BGY_CO2Y 1901
#else
    #define BGY_CO2Y 1901
#endif

/* total historical run: using CRU, NCEP, etc. ***/
#if ISIMIP_RUN==1
    #define PD_HIST 150  /* AD 1950 - 2099 */ /* ISI-MIP: 2012/06/27 by A.Ito */
                         /* and BIOFUEL RUN */
#elif ISIMIP_RUN==2
    #define PD_HIST 105  /* AD 1901 - 2005 */ /* PLUME: 2014/07/31 by A.Ito */
#elif ISIMIP_RUN==3
    /* #define PD_HIST 110 */ /* AD 1901 - 2010 */ /* ISI-MIP2 (historical gswp3, pgfv2): 2014/11/30 by A.Ito */
    /* #define PD_HIST 101 */ /* AD 1901 - 2001 */ /* ISI-MIP2 (historical watch): 2015/01/06 by A.Ito */
    #define PD_HIST 112 /* */ /* AD 1901 - 2012 */ /* ISI-MIP2 (historical watch): 2016/05/06 by A.Ito */
#elif ISIMIP_RUN==4
    #define PD_HIST 639  /* AD 1661 - 2299 */ /* ISI-MIP21.b: 2016/12/22 by A.Ito */
    /* #define PD_HIST 439 */  /* AD 1661 - 2099 */ /* ISI-MIP21.b: 2016/12/22 by A.Ito */
#elif GEOMIP_RUN==1
    #define PD_HIST 105 /* */  /* AD 1901 - 2005 --GEOMIP */
#elif NMIP_RUN>=1
    /* #define PD_HIST 155	*/	/* AD 1861 - 2015 */
    #define PD_HIST 156    /* */    /* AD 1861 - 2016 */
#elif IMPRESSIONS_RUN==2
    #define PD_HIST 115
#elif IMPRESSIONS_RUN==3
    #define PD_HIST 115
#else
    /* non-ISI-MIP: case dependent */
    /* #define PD_HIST 100  */	/* AD 1901 - 2000 */
    /* #define PD_HIST 102	*/	/* AD 1901 - 2002 */
    /* #define PD_HIST 108	*/	/* AD 1901 - 2008 */
    /* #define PD_HIST 109	*/	/* AD 1901 - 2009 */
    /* #define PD_HIST 110	*/	/* AD 1901 - 2010 */
    /* #define PD_HIST 111	*/	/* AD 1901 - 2011 */
    /* #define PD_HIST 112	*/	/* AD 1901 - 2012 */
    /* #define PD_HIST 113	*/	/* AD 1901 - 2013 */
    /* #define PD_HIST 114  */	/* AD 1901 - 2014 */
    /* #define PD_HIST 115	*/	/* AD 1901 - 2015 */
    /* #define PD_HIST 116	*/	/* AD 1901 - 2016 */
    #define PD_HIST 117    /* */    /* AD 1901 - 2017 */
#endif

/* start year (AD) of climate ***/
/* #define BGY_CLIM 1901 */
/* 1901: CRU */
/* 1990: control */
/* 1948: control */
#if ISIMIP_RUN==1
    #define BGY_CLIM 1950  /* ISI-MIP: 2012/06/27 by A.Ito */
#elif ISIMIP_RUN==2
    #define BGY_CLIM 1901  /* PLUME: 2014/07/31 by A.Ito */
#elif ISIMIP_RUN==3
    #define BGY_CLIM 1901  /* ISI-MIP2a (historical): 2014/11/30 by A.Ito */
#elif ISIMIP_RUN==4
    #define BGY_CLIM 1661  /* ISI-MIP2b (1.5/2.0 deg): 2016/12/22 by A.Ito */
#elif IMPRESSIONS_RUN==2
    #define BGY_CLIM 1901
#elif IMPRESSIONS_RUN==3
    #define BGY_CLIM 1901
#else
    #define BGY_CLIM 1901
#endif

/* historical climate (e.g., CRU) data length: 2010/01/04 (A.Ito) ***/
#if ISIMIP_RUN==1
    #define DL_HCLIM 180  /* SU 30 + AD 1950 - 2009 */
    /* note that DL_HCLIM data is not used in PLUME runs */
#elif ISIMIP_RUN==2
    /* PLUME: 2014/07/31 by A.Ito */
    #define DL_HCLIM 135  /* SU 30 + AD 1901 - 2005 */
#elif ISIMIP_RUN==3
    /* ISI-MIP2a (historical): 2014/11/30 by A.Ito */
    /* #define DL_HCLIM 140 */  /* SU 30 + AD 1901 - 2010 */
    /* #define DL_HCLIM 131 */  /* SU 30 + AD 1901 - 2001 */
    #define DL_HCLIM 142  /* SU 30 + AD 1901 - 2012 */
#elif ISIMIP_RUN==4
    /* ISI-MIP2b (1.5/2.0deg): 2016/12/24 by A.Ito */
    #define DL_HCLIM 639  /* AD 1661 - 2299 */
    /* #define DL_HCLIM 439 */  /* AD 1661 - 2099 */
#elif IMPRESSIONS_RUN==2
    #define DL_HCLIM 115
#elif IMPRESSIONS_RUN==3
    #define DL_HCLIM 115
#else
    /* non-ISI-MIP: case dependent */
    /* #define DL_HCLIM 111 */  /* AD 1901 - 2011 */
    #define DL_HCLIM 116  /* CRU TS3.25: AD 1901 - 2016 */
    /* 102: TS2.1 */
    /* 106: TS3.0 */
    /* 109: TS3.1 */
    /* 111: TS3.2 */
    /* 112: TS3.21 */
    /* 113: TS3.22 */
    /* 114: TS3.23 */
    /* 115: TS3.24 */
    /* 116: TS3.25 */
#endif

/* Simulation using NCEP/NCAR reanalysis data */
#define NCEP_RUN 1
/* 0: no  1:yes */
/* year of data beginning (AD) */
#define FDY_NCEP 1948
/* data length (years) */
/* extension to 2011: 2012/01/26 by A.Ito */
/* extension to 2012: 2013/04/14 by A.Ito */
/* extension to 2013: 2014/01/02 by A.Ito */
/* #define DL_NCEP 66 */   /* 1948-2013 */
/* #define DL_NCEP 67 */   /* 1948-2014 */
/* #define DL_NCEP 68 */   /* 1948-2015 */
/* #define DL_NCEP 69 */   /* 1948-2016 */
#define DL_NCEP 70   /* 1948-2017 */

/* Simulation using ISI-MIP data (yr) */
/* spinup 1951-1980 */
/* historical 1951-2005 */
/* projection 2006-2099 */
#if ISIMIP_RUN==1
    #define DL_ISIMIP 180  /* SU 30 + AD 1950 - 2009 */
#elif ISIMIP_RUN==2
    #define DL_ISIMIP 135  /* SU 30 + AD 1901 - 2005 */
#elif ISIMIP_RUN==3
    /* #define DL_ISIMIP 140 */  /* SU 30 + AD 1901 - 2010 */
    /* #define DL_ISIMIP 131 */ /* SU 30 + AD 1901 - 2001 */
    #define DL_ISIMIP 142  /* SU 30 + AD 1901 - 2012 */
#elif ISIMIP_RUN==4
    /* ISI-MIP2b (historical): 2014/11/30 by A.Ito */
    #define DL_ISIMIP 639  /* AD 1661 - 2299 */
    /* #define DL_ISIMIP 439 */  /* AD 1661 - 2099 */
#else
    #define DL_ISIMIP 1
#endif

/* future projection *******************************************/
/* simulation suing GCM-derived projection scenarios */
#define GCM_RUN 0
/* 0: no  1:yes */

/* year of GCM data (AD) */
#if GEOMIP_RUN==1
    #define DL_GCM 131 /* 1970-2100 --GEOMIP */
    #define BGY_GCM 2006  /* --GEOMIP PLUME */
    #define ENY_GCM 2100
    /* start year of GCM data (AD) */
    #define FDY_GCM 1970  /* --GEOMIP */
#else
    #define DL_GCM 94 /* 2006-2099 --ISI-MIP2 */
    #define BGY_GCM 2006  /* --PLUME */
    #define ENY_GCM 2099
    #define FDY_GCM 2006   /* --PLUME */
#endif
/* #define DL_GCM 241 */ /* 1860-2100 */
/* #define FDY_GCM 2001 */
/* #define FDY_GCM 1860 */

/*********************************************************/
/* NECB: coupling carbon loss */
/* 0: uncoupled */
/* 1: coupled */
/* erosion */
#define NECB_ERSN 1
/* biomass burning */
#define NECB_BB 1
/* bvoc */
#define NECB_BVOC 1
/* doc */
#define NECB_DOC 1
/* CH4 */
#define NECB_CH4 1
/* Wood harvest: 2010/10/15 by A.Ito */
#define NECB_WHVST 1
/* land-use change */
#define NECB_LUC 1
/* crop harvest */
#define NECB_CROP 1

/* land use change setting ********************************/
#define LANDUSE 10
/* 0: natural vegetation */
/* 1: no land-use change since 1901 */
/* 2: no land-use change since 1990 */
/* 3: linear land-use change since 1990 */
/* 4: high future land-use change, coupling */
/* 5: IMAGE2 land-use change since 1990 */
/* 6: EOS-WEBSTER Hurtt land-use change, 1700-2000 */
/* 7: Ramankutty land-use change, 1700-2007 */
/* 8: Hurtt harmonized land-use change, 1700-2005 (added 2010/01/31) */
/* 9: fixed land-use at 2000 (LUH 1500-2005) --GEOMIP PLUME */
/* 10: LUH 1500-2005/2005-2100 (RCP4.5) --GEOMIP */
/* 11: LUH 1500-2005/2005-2100 (RCP2.6) */
/* 12: LUH 1500-2005/2005-2100 (RCP6.0) */
/* 13: LUH 1500-2005/2005-2100 (RCP8.5) */
/* 14: SSP1 (tentative) */
/* 15: SSP2 (tentative) */
/* 16: SSP3 (tentative) */
/* 17: LUH 1500-2005/2005-2100 (RCP2.6) + Biofuel */
/* 18: ICARUS SSPs version 2016/08 */
/* 19: SSP1 (ICARUS v2016/08, RCP4.5-IPSL) */
/* 20: SSP2 (ICARUS v2016/08, RCP4.5-IPSL) */
/* 21: SSP3 (ICARUS v2016/08, RCP4.5-IPSL) */
/* 22: SSP4 (ICARUS v2016/08, RCP4.5-IPSL) */
/* 23: SSP5 (ICARUS v2016/08, RCP4.5-IPSL) */
/* 24: ISI-MIP2b land-use data (2016/12/22 by A.Ito) */
/* 25: ISI-MIP2b 2005 data (2017/11/01 by A.Ito) */

#define DL_LUC 601 /* */  /* 1500-2100 */
/* #define DL_LUC 306 */ /* 1700-2000/2005 */
/* #define DL_LUC 639 */  /* 1661-2299: ISI-MIP2b (2016/12/22 by A.Ito) */

/* begin year of land-use DATA */
#define FDY_LUC 1500 /* */
/* #define FDY_LUC 1700 */
/* #define FDY_LUC 1661 */ /* ISI-MIP2b (2016/12/22 by A.Ito) */

/* begin year of land-use SIMULATION */
#if ISIMIP_RUN==1
    #define BGY_LUC 2000    /* ISI-MIP: 2012/06/27 by A.Ito */
#elif ISIMIP_RUN==2
    #define BGY_LUC 2000    /* PLUME: 2014/07/31 by A.Ito */
#elif ISIMIP_RUN==3
    #define BGY_LUC 2000    /* ISI-MIP2 (historical): 2014/11/30 by A.Ito */
#elif ISIMIP_RUN==4
    #define BGY_LUC 1661    /* ISI-MIP2b (2016/12/22 by A.Ito) */
#else
    #define BGY_LUC 1900
#endif

/* erosion: setting of soil conservation */
#define SOIL_CONSV 0
/* 0: OECD nations = 0.75, developing = 0.95 */
/* 1: all 0.5 */

/* erosion: parameterization of vegetation cover */
#define PARA_VEGCV 0
/* 0: conventional */
/* 1: lai based */

/* constraint on maximum LAI by meta-analysis: 2015/03/23 by A.Ito */
#define CONSTRAIN_LAIMAX 0
/* 0: off */
/* 1: on */

/*******************************************************/
/* albedo perturbation experiment: 2012/12/30 by A.Ito */
#define EX_ALBEDO 0
/* 0: off */
/* 1: perturbation */
/* 2: CMIP5-max */
/* 3: CMIP5-min */
/* 4: CMIP5-mean */
/* 5: GlobAlbedo */

#define EX_TVAR 0
/* 0: off */
/* 1: albedo-induced temperature change */

/* ozone impacts: 2013/02/25 by A.Ito *****************/
#define EX_OZONE 0
/* 0: off */
/* 1: on */

/* experiment: soil decomposition temperature dependence */
#define EX_SDTD 0
/* 0: default (Lloyd & Taylor) */
/* 1: 2.0 */
/* 2: 1.5 for litter, 2.5 for humus */

/***************************************************/
/* PAR conversion */
#define D_PAR 1
/* 0: constant conversion factor */
/* 1: variable conversion factor after Dye et al. (2003) */

/* light attenuation coefficient */
#define EFF_K 1

/* SRB-based diffuse radiation estimation */
#define DIF_SRB 1
/* 0:off, 1:0 */

/* temperature dependence of plant respiration */
#define EX_TMP_RESP 0
/* 0: default */
/* 1: Yokota & Hagihara */
/* 2: Atkin */
/* 3: fix 2.0 */
/* 4: fix 1.5 */
/* 5: fix 2.5 */
/* 6: Heskel et al. 2016: added by A.Ito 2017/12/7 */

/***************************************************/
/* CH4 emission by Walter-Heimann scheme */
#define CH4_WH 1
/* 0:off, 1:0n */
#define N_SLAYER 20
/* number of soil layers */ 

/* sensitivity run of W&H scheme */
#define SENS_WHCH4 0
/* 0: control */
/* 1: t_veg * 1.1 */
/* 2: t_veg * 0.9 */
/* 3: bubble 450 microM */
/* 4: bubble 550 microM */

/* Alternative land-cover data for CH4 */
#define ALT_FWETLAND 0
/* 0: not use alternative data */
/* 1: use data */
/* 2: use Peregon-san data: 2014/02/04 */

/* inundation data */
#define ALT_INUND 0
/* 0: default (SSMI) */
/* 1: GCP-CH4  */
/* 2: IIS satellite observation */

/* specific scheme on permaforst */
#define EX_PERFROST 0
/* 0:off, 1:0n */

/* change in wetland extent due to permafrost melting: 2012/10/26 by A.Ito */
#define VAR_PFMWET 0
/* 0:off, 1:0n */

/* variable water-table depth: 2014/12/08 by A.Ito */
#define VAR_WTD 1
/* 0:off, 1:0n */

#define FIX_STMP 0
/* 0:off, 1:0n */
/* fix NPP carbon input to wetland: 2015/03/23 by A.Ito */

#define FIX_CH4_NPP 0
/* 0:off, 1:0n */

/*****************************************************/
/* parameter perturbation */
/* climate perturbation */
#define PRT_CLIM 0
/* 0:off, 1:0n */

/* fixed CO2 for gs: separation of CO2 effects */
#define FIX_GSCO2 0
/* 0:off, 1:0n */

/* carbon-nitrogen coupling */
#define CN_COUPLE 0
/* 0: no coupling (for safety) */
/* 1: physiological coupling: Amax, Rd, Decomp. */
/* 2: ecological coupling: mortality, etc. */

/* sensitivity run of biomass burning emission factor */
#define SENS_FIRE 0
/* 0:off, 1:0n */

/* parameter ensemble experiment ****/
#define ENSEMBLE_RUN 0
/* 0: off */
/* 1: on */

/* parameter ensemble */
#define N_PARA_ENS 20

/********************************************************/
/* 2016/06/08 by A.Ito */
#define EX_NITROGEN 0
/* 0: off (control) */
/* 1: half biological N2 fixation */
/* 2: double microbial abandonment */
/* 3: revised WFPS for NGAS */

/* sensitivity run of N deposition */
#define EX_CHASER_NDEPO 0
/* 0:off, 1:0n */

/* experiment: changed nitrogen deposition: 2015/08/12 by A.Ito */
#define EX_NDEPO 0
/* 0: control */
/* 1: +10% */
/* 2: -10% */
/* 3: +30% */
/* 4: -30% */
/* 5: +50% */
/* 6: -50% */
/* 7: all ammonium */
/* 8: all nitrate */
/* 9: 50:50 ammonium and nitrate */

/* sensitivity to nitrification N2O fraction: 2016/11/7 by A.Ito */
#define EX_NITR_N2O 2
/* 0: off (control) */
/* 01: 1.0% (N20-driven) */
/* 02: 0.5% (N20-driven) */
/* 03: 2.0% (N20-driven) */
/* 04: 1.0% (nitrification-driven) */
/* 05: 0.5% (nitrification-driven) */
/* 06: 2.0% (nitrification-driven) */

/* 21: DNDC (N20-driven) */
/* 22: DLEM (N20-driven) */
/* 23: CLM-CN (N20-driven) */

/* 31: DNDC (nitrification-driven) */
/* 32: DLEM (nitrification-driven) */
/* 33: CLM-CN (nitrification-driven) */

/* meta-analysis: 2017/09/15 */
/* 51: 0.051% (N20-driven) */
/* 52: 0.086% (N20-driven) */
/* 53: 0.789% (N20-driven) */
/* 54: 2.974% (N20-driven) */
/* 55: 1.612% (N20-driven) */
/* 56: 0.923% (N20-driven) */
/* 57: 2.014% (N20-driven) */

/* 61: 0.051% (nitrification-driven) */
/* 62: 0.086% (nitrification-driven) */
/* 63: 0.789% (nitrification-driven) */
/* 64: 2.974% (nitrification-driven) */
/* 65: 1.612% (nitrification-driven) */
/* 66: 0.923% (nitrification-driven) */
/* 67: 2.014% (nitrification-driven) */

/* future nitrogen fertilizer: 2016/11/22 by A.Ito  */
#define EX_NFERT 0
/* 0: off */
/* 1: HadGEM + RCP2.6 + SSP1 */
/* 2: HadGEM + RCP2.6 + SSP2 */
/* 3: HadGEM + RCP2.6 + SSP3 */
/* 4: HadGEM + RCP4.5 + SSP1 */
/* 5: HadGEM + RCP4.5 + SSP2 */
/* 6: HadGEM + RCP4.5 + SSP3 */
/* 7: MIROC + RCP2.6 + SSP1 */
/* 8: MIROC + RCP2.6 + SSP2 */
/* 9: MIROC + RCP2.6 + SSP3 */
/* 10: MIROC + RCP4.5 + SSP1 */
/* 11: MIROC + RCP4.5 + SSP2 */
/* 12: MIROC + RCP4.5 + SSP3 */
/* 101: Nishina ESSD data: 2017/02/13 by A.Ito */
#define EX_NFERT_SA 0
/* 0: reference */
/* 1: fix fertilizer */
/* 2: fix manure */
/* 3: fix deposition */
/* 4: fix land use */
/* 5: fix fertilizer + manure */
/* 6: fix fertilizer + manure + land use */
/* 7: no manure */

/* N input seasonality: 2018/02/09 by A.Ito */
#define EX_NIN_SEASON 0
/* 0: no seasonality */
/* 1: once in April (NH) or October (SH) */
/* 101: all January */

/****************************************************/
/* sensitivity analysis *****************/
/* temperature */
#define TM 0
/* 0: as it is */
/* 1: uniform warming by 1degｰC */
/* 2: gradual rise */
/* 3: descrete warming */
/* 4: spike rise */
/* 5: uniform cooling by 1degｰC */
/* 6: descrete cooling */

/* precipitation */
#define PR 0 /* */
/* 0: as it is */
/* 1: uniform precipitation increase by 10% */
/* 2: gradual rise */
/* 3: descrete increase */
/* 4: spike rise */
/* 5: uniform precipitation decrease by 10% */
/* 6: descrete decrease */

/* solar radiation */
#define SC 0
/* 0: present solar constant, 1367 W m-2 */
/* 1: increased solar constant by 1% */
/* 2: decreased solar constant by 1% */
/* 3: descrete increase */
/* 4: descrete decrease */
/* 5: increase PPFD */

/* respiration temperature dependence */
#define T_R 0
/* 0: control Q10 */
/* 1: increased Q10 */
/* 2: decreased Q10 */
/* decomposition temperture dependence */

#define T_D 0
/* 0: control Lloyd & Taylor (1994) */
/* 1: increased E0 */
/* 2: decreased E0 */
/* 3: increased T0 */
/* 4: decreased T0 */

/* decomposition moisture dependence */
#define M_D 0
/* 0: control */
/* 1: increased kd */
/* 2: decreased kd */
/* 3: increased ka */
/* 4: decreased ka */

/* acclimation */
#define EX_ACCLM 0 /* */
/* 0: control */
/* 1: soil temperature */
/* 2: soil moisture */
/* 3: litter quantity */

/* parameter sensitivity analysis */
#define SENS_PARA 0
/* 0: control */
/* 1: +10% gsmax */
/* 2: +10% Pmax */
/* 3: +10% K */
/* 4: +10% LUE */
/* 5: +10% Albedo */
/* 6: +10% WHC30/WHC */
/* 7: fixed LAI (1990s av) in 2000-2100 */

/* climate change ********************************/
/* 0:off   1:on */
/* temperature */
#define CC_T 1
/* aboveground temperature */
#define CC_T_A 1
/* belowground temperature */
#define CC_T_B 1
/* precipitation */
#define CC_P 1
/* humidity */
#define CC_H 1
/* constant future CO2 level */
#define CC_CD 1
/* 1: actual CO2 rise */
/* 2: no CO2 rise */
/* 3: fix CO2 after 2020 for GeoMIP runs */

/* deforestation ************************************/
#define EX_DEFOREST 0
/* 0: as present */
/* 1: entire deforestation, replaced by 19 */
/* 2: entire deforestation, replaced by 13 */
/* 3: entire deforestation, replaced by 31 */
/* 4: entire deforestation and regrowth: 2016/12/28 by A.Ito */

/****************************************************/
/* 0: no GCM */

/*** AR3 ***/
/* 1: CCSR/NIES A1 */
/* 2: CCSR/NIES A1FI */
/* 3: CCSR/NIES A1T */
/* 4: CCSR/NIES A2 */
/* 5: CCSR/NIES B1 */
/* 6: CCSR/NIES B2 */

/* 11: CCCma A2 */
/* 12: CCCma B2 */

/* 13: CCCma A2a */
/* 14: CCCma A2b */
/* 15: CCCma A2c */
/* 16: CCCma B2a */
/* 17: CCCma B2b */
/* 18: CCCma B2c */

/* 21: HADCM3 A2 */
/* 22: HADCM3 A2b */
/* 23: HADCM3 A2c */
/* 24: HADCM3 B2 */

/* 31: ECHAM A2 */
/* 32: ECHAM B2 */

/* 41: CSIRO A1 */
/* 42: CSIRO A2 */
/* 43: CSIRO B1 */
/* 44: CSIRO B2 */

/* 51: GFDL A2 */
/* 52: GFDL B2 */

/* 61: NCAR-PCM A1 */
/* 62: NCAR-PCM A2 */
/* 63: NCAR-PCM B2 */

/* 71: NCAR-CSM A2 */

/*** CO2 only ***/
/* 101: A1 */
/* 102: A1FI */
/* 103: A1T */
/* 104: A2 */
/* 105: B1 */
/* 106: B2 */

/*** AR4 ***/
/* 1000: MIROC-HIGH + A1B */
/* 1001: MIROC-HIGH + B1 */

/* 1010: MIROC-MED + A1B-1 */
/* 1011: MIROC-MED + A1B-2 */
/* 1012: MIROC-MED + A1B-3 */
/* 1013: MIROC-MED + A2-1 */
/* 1014: MIROC-MED + A2-2 */
/* 1015: MIROC-MED + A2-3 */
/* 1016: MIROC-MED + B1-1 */
/* 1017: MIROC-MED + B1-2 */
/* 1018: MIROC-MED + B1-3 */

/* 1050: BCCR + A1B */
/* 1051: BCCR + A2 */
/* 1052: BCCR + B1 */

/* 1060: INM-CM3 + A1B */
/* 1061: INM-CM3 + A2 */
/* 1062: INM-CM3 + B1 */

/* 1070: IPSL + A1B */
/* 1071: IPSL + A2 */
/* 1072: IPSL + B1 */

/* 1080: GFDL_CM21 + A1B */
/* 1081: GFDL_CM21 + A2 */
/* 1082: GFDL_CM21 + B1 */

/* 1090: HadCM3 + A1B */
/* 1091: HadCM3 + A2 */
/* 1092: HadCM3 + B1 */

/* 1100: GISS AOM + A1B 1 */
/* 1101: GISS AOM + A1B 2 */
/* 1102: GISS AOM + B1 1 */
/* 1103: GISS AOM + B1 2 */

/* 1110: CCCmaT63 + A1B */
/* 1111: CCCmaT63 + A2 */
/* 1112: CCCmaT63 + B1 */

/* 1120: CSIRO + A1B */
/* 1121: CSIRO + A1B */
/* 1122: CSIRO + A1B */

/* 1130: MRI + A1B 1 */
/* 1131: MRI + A1B 2 */
/* 1132: MRI + A1B 3 */
/* 1133: MRI + A1B 4 */
/* 1134: MRI + A1B 5 */
/* 1135: MRI + A2 1 */
/* 1136: MRI + A2 2 */
/* 1137: MRI + A2 3 */
/* 1138: MRI + A2 4 */
/* 1139: MRI + A2 5 */
/* 1140: MRI + B1 1 */
/* 1141: MRI + B1 2 */
/* 1142: MRI + B1 3 */
/* 1143: MRI + B1 4 */
/* 1144: MRI + B1 5 */

/* 1150: MPI + A1B 1 */
/* 1151: MPI + A1B 2 */
/* 1152: MPI + A1B 3 */
/* 1153: MPI + A1B 4 */
/* 1154: MPI + A2 1 */
/* 1155: MPI + A2 2 */
/* 1156: MPI + A2 3 */
/* 1157: MPI + B1 1 */
/* 1158: MPI + B1 2 */
/* 1159: MPI + B1 3 */

/* 1160: IAP + A1B 1 */
/* 1161: IAP + A1B 2 */
/* 1162: IAP + A1B 3 */
/* 1163: IAP + B1 1 */
/* 1164: IAP + B1 2 */
/* 1165: IAP + B1 3 */

/* 1170: HADGEM + A1B */
/* 1171: HADGEM + A2 */

/* 1180: GFDL_CM20 + A1B */
/* 1181: GFDL_CM20 + A2 */
/* 1182: GFDL_CM20 + B1 */

/* 1190: CNRM + A1B */
/* 1191: CNRM + A2 */
/* 1192: CNRM + B1 */

/* 1200: CCC47 + A1B 1 */
/* 1201: CCC47 + A1B 2 */
/* 1202: CCC47 + A1B 3 */
/* 1203: CCC47 + A1B 4 */
/* 1204: CCC47 + A1B 5 */
/* 1205: CCC47 + A2 1 */
/* 1206: CCC47 + A2 2 */
/* 1207: CCC47 + A2 3 */
/* 1208: CCC47 + A2 4 */
/* 1209: CCC47 + A2 5 */
/* 1210: CCC47 + B1 1 */
/* 1211: CCC47 + B1 2 */
/* 1212: CCC47 + B1 3 */
/* 1213: CCC47 + B1 4 */
/* 1214: CCC47 + B1 5 */

/* 1220: CCSM + A1B 1 */
/* 1221: CCSM + A1B 2 */
/* 1222: CCSM + A1B 3 */
/* 1223: CCSM + A1B 5 */
/* 1224: CCSM + A1B 6 */
/* 1225: CCSM + A1B 7 */
/* 1226: CCSM + A1B 9 */
/* 1227: CCSM + A2 1 */
/* 1228: CCSM + A2 2 */
/* 1229: CCSM + A2 3 */
/* 1230: CCSM + A2 4 */
/* 1231: CCSM + A2 5 */
/* 1232: CCSM + B1 1 */
/* 1233: CCSM + B1 2 */
/* 1234: CCSM + B1 3 */
/* 1235: CCSM + B1 4 */
/* 1236: CCSM + B1 5 */
/* 1237: CCSM + B1 6 */
/* 1238: CCSM + B1 7 */
/* 1239: CCSM + B1 9 */

/* 1240: GISS E + A1B 1 */
/* 1241: GISS E + A1B 2 */
/* 1242: GISS E + A1B 3 */

/* 1250: GISS R + A1B 2 */
/* 1251: GISS R + A1B 4 */
/* 1252: GISS R + A2 3 */
/* 1253: GISS R + B1 3 */

/* 1260: NCAR PCM + A1B 1 */
/* 1261: NCAR PCM + A1B 2 */
/* 1262: NCAR PCM + A1B 3 */
/* 1263: NCAR PCM + A1B 4 */
/* 1264: NCAR PCM + A2 1 */
/* 1265: NCAR PCM + A2 2 */
/* 1266: NCAR PCM + A2 3 */
/* 1267: NCAR PCM + A2 4 */

/*** ISI-MIP Phase 1: 2012/06/27 by A.Ito ********/
/* CD-LINK 2016/11/17 */
/* 2001: HadGEM2-ES RCP 2.6 +co2 */
/* 2002: HadGEM2-ES RCP 8.5 +co2 */
/* 2003: HadGEM2-ES RCP 4.5 +co2 */
/* 2004: HadGEM2-ES RCP 6.0 +co2 */
/* 2005: HadGEM2-ES RCP 2.6 +noco2 */
/* 2006: HadGEM2-ES RCP 8.5 +noco2 */
/* 2007: HadGEM2-ES RCP 4.5 +noco2 */
/* 2008: HadGEM2-ES RCP 6.0 +noco2 */

/* 2011: IPSL-CM5A-LR RCP 2.6 +co2 */
/* 2012: IPSL-CM5A-LR RCP 8.5 +co2 */
/* 2013: IPSL-CM5A-LR RCP 4.5 +co2 */
/* 2014: IPSL-CM5A-LR RCP 6.0 +co2 */
/* 2015: IPSL-CM5A-LR RCP 2.6 +noco2 */
/* 2016: IPSL-CM5A-LR RCP 8.5 +noco2 */
/* 2017: IPSL-CM5A-LR RCP 4.5 +noco2 */
/* 2018: IPSL-CM5A-LR RCP 6.0 +noco2 */

/* 2021: GFDL-ESM2M RCP 2.6 +co2 */
/* 2022: GFDL-ESM2M RCP 8.5 +co2 */
/* 2023: GFDL-ESM2M RCP 4.5 +co2 */
/* 2024: GFDL-ESM2M RCP 6.0 +co2 */
/* 2025: GFDL-ESM2M RCP 2.6 +noco2 */
/* 2026: GFDL-ESM2M RCP 8.5 +noco2 */
/* 2027: GFDL-ESM2M RCP 4.5 +noco2 */
/* 2028: GFDL-ESM2M RCP 6.0 +noco2 */

/* 2031: MIROC-ESM-CHEM RCP 2.6 +co2 */
/* 2032: MIROC-ESM-CHEM RCP 8.5 +co2 */
/* 2033: MIROC-ESM-CHEM RCP 4.5 +co2 */
/* 2034: MIROC-ESM-CHEM RCP 6.0 +co2 */
/* 2035: MIROC-ESM-CHEM RCP 2.6 +noco2 */
/* 2036: MIROC-ESM-CHEM RCP 8.5 +noco2 */
/* 2037: MIROC-ESM-CHEM RCP 4.5 +noco2 */
/* 2038: MIROC-ESM-CHEM RCP 6.0 +noco2 */

/* 2041: NorESM1-M RCP 2.6 +co2 */
/* 2042: NorESM1-M RCP 8.5 +co2 */
/* 2043: NorESM1-M RCP 4.5 +co2 */
/* 2044: NorESM1-M RCP 6.0 +co2 */
/* 2045: NorESM1-M RCP 2.6 +noco2 */
/* 2046: NorESM1-M RCP 8.5 +noco2 */
/* 2047: NorESM1-M RCP 4.5 +noco2 */
/* 2048: NorESM1-M RCP 6.0 +noco2 */

/** GEO-MIP: 2013/11/26 by A.Ito ***********/
/* revised 2015/07/21 */
/* 3000: BNU-ESM RCP4.5 */
/* 3001: BNU-ESM G3 */
/* 3002: BNU-ESM G4 */

/* 3010: CSIRO-mk3L-1-2 RCP4.5 */
/* 3011: CSIRO-mk3L-1-2 G4 */
/* 3012: CSIRO-mk3L-1-2 G3S */

/* 3020: GISS-EL-R RCP4.5 */
/* 3021: GISS-EL-R G3 */
/* 3022: GISS-EL-R G4 */

/* 3030: HadGEM2-ES RCP4.5 */
/* 3031: HadGEM2-ES G3 */
/* 3032: HadGEM2-ES G4 */
/* 3033: HadGEM2-ES G3S */
/* 3034: HadGEM2-ES G4cdcn */
/* 3035: HadGEM2-ES G4seaSalt */

/* 3040: IPSL-CM54-LR RCP4.5 */
/* 3041: IPSL-CM54-LR G3 */
/* 3042: IPSL-CM54-LR G5 */

/* 3050: MIROC-ESM RCP4.5 */
/* 3051: MIROC-ESM G4 */
/* 3052: MIROC-ESM G4cdcn */

/* 3060: MIROC-ESM-CHEM RCP4.5 */
/* 3061: MIROC-ESM-CHEM G4 */

/* 3070: CCCma RCP4.5 */
/* 3071: CCCma G4 */
/* 3072: CCCma G4ndcn */

/* 3080: MPI-ESM-LR RCP4.5 */
/* 3081: MPI-ESM-LR G3 */

/* 3090: CCSM4 RCP4.5 */
/* 3091: CCSM4 G3S */

/* 3100: NorESM RCP4.5 */
/* 3101: NorESM G4cdcn */

/** PLUME: 2014/07/31 by A.Ito *************/
/* 4011: GFDL RCP 4.5 */
/* 4012: GFDL RCP 8.5 */

/* 4021: IPSL RCP 4.5 */
/* 4022: IPSL RCP 8.5 */
/* 4023: IPSL RCP 2.6 */
/* 4024: IPSL RCP 6.0 */

/** ISI-MIP2a: 2014/11/30 by A.Ito ***********/
/* 5001: historical GSWP3 */
/* 5002: historical PGFv2.1 (2016/05/06) */
/* 5003: historical WATCH */
/* 5004: historical WATCH+WFDEI */

/** ISI-MIP2b: 2016/12/22 and 2017/05/17 by A.Ito  *****/
/* 5010: GFDL piControl */
/* 5011: GFDL piControl + historical + rcp2.6 */
/* 5012: GFDL piControl + historical + rcp6.0 */
/* 5020: IPSL piControl */
/* 5021: IPSL piControl + historical + rcp2.6 */
/* 5022: IPSL piControl + historical + rcp6.0 */
/* 5030: MIROC5 piControl */
/* 5031: MIROC5 piControl + historical + rcp2.6 */
/* 5032: MIROC5 piControl + historical + rcp6.0 */
/* 5040: HadGEM2-ES piControl */
/* 5041: HadGEM2-ES piControl + historical + rcp2.6 */
/* 5042: HadGEM2-ES piControl + historical + rcp6.0 */

/** IMPRESSIONS: 2015/07/17 **/
/* 6001: phase 1 */
/* temperature */
/* IMPRESSIONS_DT */
/* 0: temperature 0 K */
/* 1: temperature -3 K */
/* 2: temperature -2 K */
/* 3: temperature -1 K */
/* 4: temperature  0 K */
/* 5: temperature +1 K */
/* 6: temperature +2 K */
/* 7: temperature +3 K */
/* 8: temperature +4 K */
/* 9: temperature +5 K */
/* 10: temperature +6 K */
/* 11: temperature +7 K */
/* 12: temperature +8 K */
/* 13: temperature +9 K */
/* 14: temperature +10 K */
/* 15: temperature +11 K */

/* IMPRESSIONS_DP */
/* 0: precipitation 0% */
/* 1: precipitation -60% */
/* 2: precipitation -50% */
/* 3: precipitation -40% */
/* 4: precipitation -30% */
/* 5: precipitation -20% */
/* 6: precipitation -10% */
/* 7: precipitation   0% */
/* 8: precipitation +10% */
/* 9: precipitation +20% */
/* 10: precipitation +30% */
/* 11: precipitation +40% */
/* 12: precipitation +50% */
/* 13: precipitation +60% */

/** IMPRESSIONS: 2017/04/19 **/
/* 6002: phase 2 */
/* #define IMPRESSIONS_DT */
/* 0: temperature 0 K */
/* 1: temperature -3 K */
/* 2: temperature -2 K */
/* 3: temperature -1 K */
/* 4: temperature  0 K */
/* 5: temperature +0.5 K */
/* 6: temperature +1 K */
/* 7: temperature +2 K */
/* 8: temperature +3 K */
/* 9: temperature +4 K */
/* 10: temperature +5 K */
/* 11: temperature +6 K */
/* 12: temperature +7 K */
/* 13: temperature +8 K */

/* #define IMPRESSIONS_DP */
/* 0: precipitation 0% */
/* 1: precipitation -42% */
/* 2: precipitation -36% */
/* 3: precipitation -30% */
/* 4: precipitation -25% */
/* 5: precipitation -24% */
/* 6: precipitation -20% */
/* 7: precipitation -18% */
/* 8: precipitation -15% */
/* 9: precipitation -12% */
/* 10: precipitation -10% */
/* 11: precipitation -9% */
/* 12: precipitation -6% */
/* 13: precipitation -5% */
/* 14: precipitation -3% */
/* 15: precipitation  0% */
/* 16: precipitation +3% */
/* 17: precipitation +5% */
/* 18: precipitation +6% */
/* 19: precipitation +9% */
/* 20: precipitation +10% */
/* 21: precipitation +12% */
/* 22: precipitation +15% */
/* 23: precipitation +18% */

/* ICARUS 2016/08/12 ************/
/* 2201 (2101): gfdl SSP1 rcp2.6  */
/* 2202 (2102): gfdl SSP1 rcp4.5  */
/* 2203 (2103): gfdl SSP2 rcp2.6  */
/* 2204 (2104): gfdl SSP2 rcp4.5  */
/* 2205 (2105): gfdl SSP2 rcp6.0  */
/* 2206 (2106): gfdl SSP3 rcp4.5  */
/* 2207 (2107): gfdl SSP3 rcp6.0  */
/* 2208 (2108): gfdl SSP4 rcp2.6  */
/* 2209 (2109): gfdl SSP4 rcp4.5  */
/* 2210 (2110): gfdl SSP5 rcp2.6  */
/* 2211 (2111): gfdl SSP5 rcp4.5  */

/* 2221 (2201): hadgem SSP1 rcp2.6  */
/* 2222 (2202): hadgem SSP1 rcp4.5  */
/* 2223 (2203): hadgem SSP2 rcp2.6  */
/* 2224 (2204): hadgem SSP2 rcp4.5  */
/* 2225 (2205): hadgem SSP2 rcp6.0  */
/* 2226 (2206): hadgem SSP3 rcp4.5  */
/* 2227 (2207): hadgem SSP3 rcp6.0  */
/* 2228 (2208): hadgem SSP4 rcp2.6  */
/* 2229 (2209): hadgem SSP4 rcp4.5  */
/* 2230 (2210): hadgem SSP5 rcp2.6  */
/* 2231 (2211): hadgem SSP5 rcp4.5  */

/* 2241 (2301): ipsl SSP1 rcp2.6  */
/* 2242 (2302): ipsl SSP1 rcp4.5  */
/* 2243 (2303): ipsl SSP2 rcp2.6  */
/* 2244 (2304): ipsl SSP2 rcp4.5  */
/* 2245 (2305): ipsl SSP2 rcp6.0  */
/* 2246 (2306): ipsl SSP3 rcp4.5  */
/* 2247 (2307): ipsl SSP3 rcp6.0  */
/* 2248 (2308): ipsl SSP4 rcp2.6  */
/* 2249 (2309): ipsl SSP4 rcp4.5  */
/* 2250 (2310): ipsl SSP5 rcp2.6  */
/* 2251 (2311): ipsl SSP5 rcp4.5  */

/* 2261 (2401): miroc SSP1 rcp2.6  */
/* 2262 (2402): miroc SSP1 rcp4.5  */
/* 2263 (2403): miroc SSP2 rcp2.6  */
/* 2264 (2404): miroc SSP2 rcp4.5  */
/* 2265 (2405): miroc SSP2 rcp6.0  */
/* 2266 (2406): miroc SSP3 rcp4.5  */
/* 2267 (2407): miroc SSP3 rcp6.0  */
/* 2268 (2408): miroc SSP4 rcp2.6  */
/* 2269 (2409): miroc SSP4 rcp4.5  */
/* 2270 (2410): miroc SSP5 rcp2.6  */
/* 2271 (2411): miroc SSP5 rcp4.5  */

/* 2281 (2501): noresm SSP1 rcp2.6  */
/* 2282 (2502): noresm SSP1 rcp4.5  */
/* 2283 (2503): noresm SSP2 rcp2.6  */
/* 2284 (2504): noresm SSP2 rcp4.5  */
/* 2285 (2505): noresm SSP2 rcp6.0  */
/* 2286 (2506): noresm SSP3 rcp4.5  */
/* 2287 (2507): noresm SSP3 rcp6.0  */
/* 2288 (2508): noresm SSP4 rcp2.6  */
/* 2289 (2509): noresm SSP4 rcp4.5  */
/* 2290 (2510): noresm SSP5 rcp2.6  */
/* 2291 (2511): noresm SSP5 rcp4.5  */

/* add: 2016/10/17 */
/* 2601: SSP1 no climate change  */
/* 2602: SSP2 no climate change  */
/* 2603: SSP3 no climate change  */

