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
#define lTs (3600.0*12.0/100000000.0) /* from micro-mol m-2 s-1 to Mg C ha-1 day-1*/
#define ZAT 273.15 /* zero degree centigrade in absolute temperature */
#define STCIR (0.0111/0.9889) /* standard stable carbon isotope ratio */
#define UGC 8.314 /* universal gas constant */
#define SLC 1367.0
#define SBC (5.6703 / 100000000.0) /* Stephan-Boltzman Constant, W m-2 K-4 */
#define GAC 9.8 /* gravity acceleration constant, m s-2 */

/***********************************************************/
/* annual time-step, 12=monthly */
#define ASTEP 12

/* grid configulation */
#define N_ROW 360
#define N_COL 720

/* initial (minimal) carbon stock ***********/
#define INT_C 0.01

/***********************************************************/
#define ISIMIP_RUN 3
/* 0: normal (no ISI-MIP) */
/* 1: ISI-MIP 1st-phase runs */
/* 2: PLUME (ISI-MIP Phase 2) runs : 2014/07/31 by A.Ito */
/* 3: ISI-MIP 2nd-phase historical runs */

#define GEOMIP_RUN 0
/* 0: normal (no GEO-MIP) */
/* 1: GEO-MIP runs */

/********************************************************/
/* output text files */
#define OUTPUT_CARBON1 1
#define OUTPUT_CARBON2 1
#define OUTPUT_ISOTOPE 0
#define OUTPUT_NITROGEN 0
#define OUTPUT_HYDMET 1
#define OUTPUT_EROSION 0
#define OUTPUT_GHG 0
#define OUTPUT_BB 0
#define OUTPUT_BVOC 0
/* output binary */
#define C13_GOUT 1
#define C14_GOUT 1
#define PHYS_GOUT 1

/***********************************************************/
/* total vegetation number */
#define NVEG_OLSON 34	/* Olson veg (modified) */
#define NVEG_SAGE 16		/* SAGE veg (modified) */
#define NVEG_CROP 3		/* crop types */

/* calculation for land covers */
#define CALC_OLSON 1    /* matural vegetation */
#define CALC_CROP 1     /* cropland */
/* 0:off 1:on */

/* Olson croplands replaced by SAGE natural vegetation */
#define REPL_OLSON_CROP 1
/* 0:off 1:on */

/* number of geographical regions */
#define N_REG 23
/* defined in region_giorgi() in vegetdeal.c */

/* atmopsheric GHG data length */
#if ISIMIP_RUN==1
    /* ISI-MIP: 2012/06/27 by A.Ito */
    #define DL_AGHG 736
    #define BGY_AGHG 1765
#elif ISIMIP_RUN==2
    /* PLUME: 2014/07/31 by A.Ito */
    #define DL_AGHG 736
    #define BGY_AGHG 1765
#elif ISIMIP_RUN==3
    /* ISI-MIP2: 2014/12/02 by A.Ito */
    #define DL_AGHG 736
    #define BGY_AGHG 1765
#elif GEOMIP_RUN==1
    /* GEO-MIP: 2012/06/27 by A.Ito */
    #define DL_AGHG 736
    #define BGY_AGHG 1765
#else
    #define DL_AGHG 553 /* default */
    #define BGY_AGHG 1750
#endif

/***********************************************************/
/* simulation framework duration (years) */
#define PD_SIM 201	
/* only for memory preparation; not actual period */

/* start year (AD) of CO2 time series */
/* #define PIVOT_CO2Y 1901 */
/* cru-init 1901 */
/* con 1990 */
/* dcd 2081 */
/* NCEP1 1948 */
#if ISIMIP_RUN==1
    #define PIVOT_CO2Y 1950  /* ISI-MIP: 2012/06/27 by A.Ito */
#elif ISIMIP_RUN==2
    #define PIVOT_CO2Y 1901  /* PLUME: 2014/07/31 by A.Ito */
#elif ISIMIP_RUN==3
    #define PIVOT_CO2Y 1901  /* ISI-MIP2 (historical): 2014/11/30 by A.Ito */
#else
    #define PIVOT_CO2Y 1901
#endif

/* total historical run: using CRU, NCEP, etc. ***/
#if ISIMIP_RUN==1
    #define PD_HIST 150  /* AD 1950 - 2099 */ /* ISI-MIP: 2012/06/27 by A.Ito */
#elif ISIMIP_RUN==2
    #define PD_HIST 105  /* AD 1901 - 2005 */ /* PLUME: 2014/07/31 by A.Ito */
#elif ISIMIP_RUN==3
    #define PD_HIST 110  /* AD 1901 - 2010 */ /* ISI-MIP2 (historical): 2014/11/30 by A.Ito */
#elif GEOMIP_RUN==1
    #define PD_HIST 105 /* */  /* AD 1901 - 2005 --GEOMIP */
#else
    /* non-ISI-MIP: case dependent */
    /* #define PD_HIST 100  */	/* AD 1901 - 2000 */
    /* #define PD_HIST 102	*/	/* AD 1901 - 2002 */
    /* #define PD_HIST 108	*/	/* AD 1901 - 2008 */
    /* #define PD_HIST 109	*/	/* AD 1901 - 2009 */
    /* #define PD_HIST 110	*/	/* AD 1901 - 2010 */
    /* #define PD_HIST 111	*/	/* AD 1901 - 2011 */
    /* #define PD_HIST 112	*/	/* AD 1901 - 2012 */
    #define PD_HIST 113	/* */	/* AD 1901 - 2013 */
#endif

/* start year (AD) of climate ***/
/* #define PIVOT_CLIMY 1901 */
/* 1901: CRU */
/* 1990: control */
/* 1948: control */
#if ISIMIP_RUN==1
    #define PIVOT_CLIMY 1950  /* ISI-MIP: 2012/06/27 by A.Ito */
#elif ISIMIP_RUN==2
    #define PIVOT_CLIMY 1901  /* PLUME: 2014/07/31 by A.Ito */
#elif ISIMIP_RUN==3
    #define PIVOT_CLIMY 1901  /* ISI-MIP2 (historical): 2014/11/30 by A.Ito */
#else
    #define PIVOT_CLIMY 1901
#endif

/* historical (e.g., CRU) data length: 2010/01/04 (A.Ito) ***/
#if ISIMIP_RUN==1
    #define DL_CRU 180  /* SU 30 + AD 1950 - 2009 */
    /* note that DL_CRU data is not used in PLUME runs */
#elif ISIMIP_RUN==2
    /* PLUME: 2014/07/31 by A.Ito */
    #define DL_CRU 135  /* SU 30 + AD 1901 - 2005 */
#elif ISIMIP_RUN==3
    /* ISI-MIP2 (historical): 2014/11/30 by A.Ito */
    #define DL_CRU 140  /* SU 30 + AD 1901 - 2010 */
#else
    /* non-ISI-MIP: case dependent */
    /* #define DL_CRU 111 */  /* AD 1901 - 2011 */
    #define DL_CRU 113  /* CRU TS3.21: AD 1901 - 2012 */
    /* 102: TS2.1 */
    /* 106: TS3.0 */
    /* 109: TS3.1 */
    /* 111: TS3.2 */
    /* 112: TS3.21 */
    /* 113: TS3.22 */
#endif

/* Simulation using NCEP/NCAR reanalysis data */
#define NCEP_RUN 0
/* 0: no  1:yes */
/* year of data beginning (AD) */
#define PIVOT_NCEP 1948
/* data length (years) */
/* extension to 2011: 2012/01/26 by A.Ito */
/* extension to 2012: 2013/04/14 by A.Ito */
/* extension to 2013: 2014/01/02 by A.Ito */
#define DL_NCEP 66   /* 1948-2013 */

/* Simulation using ISI-MIP data (yr) */
/* spinup 1951-1980 */
/* historical 1951-2005 */
/* projection 2006-2099 */
#if ISIMIP_RUN==1
    #define DL_ISIMIP 180  /* SU 30 + AD 1950 - 2009 */
#elif ISIMIP_RUN==2
    #define DL_ISIMIP 135  /* SU 30 + AD 1901 - 2005 */
#elif ISIMIP_RUN==3
    #define DL_ISIMIP 140  /* SU 30 + AD 1901 - 2010 */
#else
    #define DL_ISIMIP 0
#endif

/* future projection *****************************/
/* simulation suing GCM-derived projection scenarios */
#define GCM_RUN 0
/* 0: no  1:yes */
/* #define GCM_PD 100 */	/* 100 : 2001-2100 */
/* #define GCM_PD 99 */	/* 99 : 2001-2099 */
/* year of data beginning (AD) */
/* #define BGY_GCM 2001 */
#define BGY_GCM 2006  /* --GEOMIP PLUME */
/* #define ENY_GCM 2100 */
#define ENY_GCM 2099  /* --PLUME */

/* GCM data length */
/* #define DL_GCM 131 */ /* 1970-2100 --GEOMIP */
#define DL_GCM 94 /* */ /* 2006-2099 --GEOMIP */
/* #define DL_GCM 241 */ /* 1860-2100 */

/* start year of GCM data (AD) */
/* #define PIVOT_GCMY 2001 */
/* #define PIVOT_GCMY 1970 */  /* --GEOMIP */
#define PIVOT_GCMY 2006   /* --PLUME */
/* #define PIVOT_GCMY 1860 */

/***************************************************/
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

/* parameter ensemble *************************************/
#define NPERT 20

/* land use change setting ***********/
#define LANDUSE 9
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

//#define DL_LUH 306 /* 1700-2000/2005 */
#define DL_LUH 601 /* 1500-2100 */

/* begin year of land-use data */
#define PIVOT_LUC 1500
//#define PIVOT_LUC 1700

#if ISIMIP_RUN==1
    #define BGY_LUC 2000    /* ISI-MIP: 2012/06/27 by A.Ito */
#elif ISIMIP_RUN==2
    #define BGY_LUC 2000    /* PLUME: 2014/07/31 by A.Ito */
#elif ISIMIP_RUN==3
    #define BGY_LUC 2000    /* ISI-MIP2 (historical): 2014/11/30 by A.Ito */
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

/***************************************************/
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

/* ozone impacts: 2013/02/25 by A.Ito *************/
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

/***************************************************/
/* CH4 emission by Walter-Heimann scheme */
#define CH4_WH 0
/* 0:off, 1:0n */
#define SOIL_LAYER 20
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

/***************************************************/
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

/* sensitivity run of N deposition */
#define SENS_N 0
/* 0:off, 1:0n */

/********************************************************/
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
#define BACC 0 /* */
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

/*** ISI-MIP Phase 1: 2012/06/27 by A.Ito ***/
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
/* 3000: BNU-ESM RCP4.5 */
/* 3003: BNU-ESM G3 */
/* 3004: BNU-ESM G4 */

/* 3100: CSIRO-mk3L-1-2 RCP4.5 */
/* 3104: CSIRO-mk3L-1-2 G4 */

/* 3200: GISS-EL-R RCP4.5 */
/* 3203: GISS-EL-R G3 */
/* 3204: GISS-EL-R G4 */

/* 3300: HadGEM2-ES RCP4.5 */
/* 3303: HadGEM2-ES G3 */
/* 3304: HadGEM2-ES G4 */
/* 3313: HadGEM2-ES G3S */

/* 3400: IPSL-CM54-LR RCP4.5 */
/* 3403: IPSL-CM54-LR G3 */
/* 3405: IPSL-CM54-LR G5 */

/* 3500: MIROC-ESM RCP4.5 */
/* 3504: MIROC-ESM G4 */

/* 3600: MIROC-ESM-CHEM RCP4.5 */
/* 3604: MIROC-ESM-CHEM G4 */

/* 3700: CCCma RCP4.5 */
/* 3704: CCCma G3 */

/* 3800: MPI-ESM-LR RCP4.5 */
/* 3803: MPI-ESM-LR G3 */

/* 3900: CCSM4 RCP4.5 */
/* 3913: CCSM4 G3S */

/** PLUME: 2014/07/31 by A.Ito ***********/
/* 4011: GFDL RCP 4.5 */
/* 4012: GFDL RCP 8.5 */

/* 4021: IPSL RCP 4.5 */
/* 4022: IPSL RCP 8.5 */
/* 4023: IPSL RCP 2.6 */
/* 4024: IPSL RCP 6.0 */

/** ISI-MIP2: 2014/11/30 by A.Ito ***********/
/* 5001: historical GSWP3 */
/* 5002: historical PGFv2 */

