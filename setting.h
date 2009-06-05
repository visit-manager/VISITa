/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in November 29, 2007						*/

/* prototype declaration */
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

/***************************************************/
/* total vegetation number */
#define VEG_NUM_OLSON 34		/* Olson veg (modified) */
#define VEG_NUM_SAGE 16		/* SAGE veg (modified) */
#define VEG_NUM_CROP 3		/* crop types */

/* 0:off 1:on */
#define CALC_OLSON 1
#define CALC_SAGE 0
#define CALC_CROP 0

/***************************************************/
/* 20th century */
#define CRU_PD 100  /* */		/* 100: AD 1901 - 2000 */
/* #define CRU_PD 102	*/	/* 102: AD 1901 - 2002 */
/* #define CRU_PD 108	*/	/* 102: AD 1901 - 2008 */
/* 21st century */
#define GCM_PD 100		/* 100 : 2001-2100 */
/* #define GCM_PD 99 */	/* 99 : 2001-2099 */

/* start year of CO2 */
#define PIVOT_CO2Y 1901
/* cru-init 1901 */
/* con 1990 */
/* dcd 2081 */
/* NCEP1 1948 */
/* start year of climate */
#define PIVOT_CLIMY 1901 
/* 1901: CRU */
/* 1990: control */
/* 1948: control */
/* #define PIVOT_GCMY 2001 */
/* start year of GCM climate */
#define PIVOT_GCMY 1970

/***************************************************/
/* coupling erosion carbon loss */
#define ERSN_CC 0
/* 0: uncoupled */
/* 1: coupled */

/* land use setting */
#define LANDUSE 6
/* 0: natural vegetation */
/* 1: no land-use change since 1901 */
/* 2: no land-use change since 1990 */
/* 3: linear land-use change since 1990 */
/* 4: high future land-use change, coupling */
/* 5: IMAGE2 land-use change since 1990 */
/* 6: EOS-WEBSTER Hurtt land-use change, 1700-2000 */

/* setting of soil conservation */
#define SOIL_CONSV 0
/* 0: OECD nations = 0.75, developing = 0.95 */
/* 1: all 0.5 */

/* parameterization of vegetation cover */
#define VEGCOVER 0
/* 0: conventional */
/* 1: lai based */

/***************************************************/
/* PAR conversion */
#define D_PAR 1
/* 0: constant conversion factor */
/* 1: variable conversion factor after Dye et al. (2003) */

/* light attenuation coefficient */
#define EFF_K 1

/* SRB-based diffuse radiation estimation */
#define DIF_SRB 1
/* 0: off */
/* 1: on */

/***************************************************/
/* carbon-nitrogen coupling */
#define CN_COUPLE 0
/* 0: no coupling (for safety) */
/* 1: physiological coupling: Amax, Rd, Decomp. */
/* 2: ecological coupling: mortality, etc. */

/***************************************************/
/* sensitivity analysis *****************/
/* temperature */
#define TM 0
/* 0: as it is */
/* 1: uniform warming by 1ﾃ釆ｩﾂｰC */
/* 2: gradual rise */
/* 3: descrete warming */
/* 4: spike rise */
/* 5: uniform cooling by 1ﾃ釆ｩﾂｰC */
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
#define SENS 0
/* 0: control */
/* 1: +10% gsmax */
/* 2: +10% Pmax */
/* 3: +10% K */
/* 4: +10% LUE */
/* 5: +10% Albedo */
/* 6: +10% WHC30/WHC */

/* 7: constant LAI (1990s av) in 2001-2100 */

/* climate change ************************/
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

/* deforestation ***************************/
#define DEFOREST 0
/* 0: as present */
/* 1: entire deforestation, replaced by 19 */
/* 2: entire deforestation, replaced by 13 */
/* 3: entire deforestation, replaced by 31 */

/*******************************************
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
