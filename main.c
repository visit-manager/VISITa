/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* History */
/* Revised August 15, 2007				*/
/* Revised August 19, 2007				*/

/* VISIT a: global model	*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"structure.h"
#include"prototype.h"

/* structure *****************/
struct Grid grid;
struct Loct loct;
struct Echar echar;
struct Mass mass;
struct Flux flux;

struct Loct loct_agr;
struct Echar echar_agr;
struct Mass mass_agr;
struct Flux flux_agr;

struct Loct loct_nat;
struct Echar echar_nat;
struct Mass mass_nat;
struct Flux flux_nat;

/* GCM variables ************/
long GCM, CO2S, GCM_R, GCM_C; /* */

/* atm. GHG scenario */
double aco2_a1[553], aco2_a2[553], aco2_b1[553], aco2_b2[553];
double ach4_a1[553], ach4_a2[553], ach4_b1[553], ach4_b2[553];
double an2o_a1[553], an2o_a2[553], an2o_b1[553], an2o_b2[553];

/* historical results *************/
double h_tmp[201], h_pre[201], h_dswr[201], h_aet[201], h_rof[201];
double h_gpp[201], h_npp[201], h_nep[201], h_plant[201], h_soil[201];
double h_sr[201], h_ersn_c[201];
double h_agrarea[201], h_luc[201];

double h_burnt_area[201];
double h_bioburn_co2[201], h_bioburn_ch4[201], h_bioburn_co[201];
double h_bioburn_nmhc[201], h_bioburn_oc[201], h_bioburn_bc[201];
double h_bioburn_nox[201], h_bioburn_so2[201], h_bioburn_pm25[201];
double h_bioburn_tpm[201], h_bioburn_tec[201];

double h_ch4ox1[201], h_ch4ox2[201], h_ch4ox3[201];
double h_ch4emit_cao_paddy[201], h_ch4emit_cao_wetland[201];
double h_n2o_emit_ngas[201], h_n2_emit_ngas[201];
double h_n2o_emit_casa[201], h_no_emit_casa[201], h_n2_emit_casa[201];
double h_nh3_emit[201];

double h_voc_isopr_g97[201], h_voc_monotrp_g97[201], h_voc_methanl_g97[201];
double h_voc_acetone_g97[201], h_voc_actaldhd_g97[201], h_voc_frmardhd_g97[201];
double h_voc_formacd_g97[201], h_voc_acetacd_g97[201], h_voc_co_g97[201];

/* monthly results **********/
double m_ch4ox1[12], m_ch4ox2[12], m_ch4ox3[12];
double m_bioburn_co2[12], m_bioburn_ch4[12], m_bioburn_co[12];
double m_bioburn_nmhc[12], m_bioburn_oc[12], m_bioburn_bc[12];

/* vegetation (olson) results */
double v_area[34];
double v_gpp[34], v_npp[34], v_nep[34];
double v_lai[34], v_fol[34], v_stm[34], v_rot[34], v_ltr[34], v_msl[34];

/* 070725 convensional conservation factor 0.5 */
/* 070802 revised conservation factor developed - developing */

/*** main simulation roop **************************************/
int main(
	int argc, 
	char *argv[]
){
	long f, g, h, zone;
	double glandarea;
	char filename[100];
	char s_date[25];
	char num[4];
	char s_case[25];
	char s_carbon[25];
	char s_nitrogen[25];
	char s_erosion[25];
	char s_ghg[25];
	char s_bioburn[25];
	char s_voc[25];
		
	FILE *fp_s[26];
	FILE *fp_c[40];
	FILE *fp_carbon;
	FILE *fp_nitrogen;
	FILE *fp_ersn;
	FILE *fp_ghg;
	FILE *fp_bioburn;	/* biomass burning: 070805 */
	FILE *fp_voc;	/* voc emission: 070905 */
	FILE *fp_glsum;

	if(argc<3){
		printf("No scenario specified !!!\n");
		exit(1);
	}
	
	GCM = atol(argv[1]);
	if(GCM>=0 && GCM<=2000){
		;
	}else{
		printf("Bad scenario specified !!!\n");
		exit(1);
	}
	
	/* GCM = 1003; */
	
	/* 0: no GCM */

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

	/* 101: A1 */
	/* 102: A1FI */
	/* 103: A1T */
	/* 104: A2 */
	/* 105: B1 */
	/* 106: B2 */
	
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

	/***************************************************/
	if(GCM==0){
		CO2S = 4;
		fp_s[23]=fopen("./data/image_a2_fcrop.dat","rt");
	}else if(GCM==1 || GCM==41 || GCM==61){
		CO2S = 1;
		fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt");
	}else if(GCM==2){
		CO2S = 2;
		fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt");
	}else if(GCM==3){
		CO2S = 3;
		fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt");
	}else if(GCM==4 || GCM==11 || GCM==13 || GCM==14 || GCM==15 || GCM==21 || GCM==22 || 
			GCM==23 || GCM==31 || GCM==42 || GCM==51 || GCM==62 || GCM==71 ){
		CO2S = 4;
		fp_s[23]=fopen("./data/image_a2_fcrop.dat","rt");
	}else if(GCM==5 || GCM==43){
		CO2S = 5;
		fp_s[23]=fopen("./data/image_b1_fcrop.dat","rt");
	}else if(GCM==6 || GCM==12 || GCM==16 || GCM==17 || GCM==18 || GCM==24 || GCM==32 || 
			GCM==44 || GCM==52 || GCM==63){
		CO2S = 6;
		fp_s[23]=fopen("./data/image_b1_fcrop.dat","rt");
	}
	
	if(GCM==1000 ||GCM==1010 ||GCM==1011 ||GCM==1012 ||GCM==1050 ||GCM==1060 ||GCM==1070 ||
		GCM==1080 ||GCM==1090 ||GCM==1100 ||GCM==1101 ||GCM==1110 ||GCM==1120 ||
		GCM==1130 ||GCM==1131 ||GCM==1132 ||GCM==1133 ||GCM==1134 ||GCM==1150 ||
		GCM==1151 ||GCM==1152 ||GCM==1153 ||GCM==1160 ||GCM==1161 ||GCM==1162 ||
		GCM==1170 ||GCM==1180 ||GCM==1190 ||GCM==1200 ||GCM==1201 ||GCM==1202 ||
		GCM==1203 ||GCM==1204 ||GCM==1220 ||GCM==1221 ||GCM==1222 ||GCM==1223 ||
		GCM==1224 ||GCM==1225 ||GCM==1226 ||GCM==1240 ||GCM==1250 ||GCM==1251 ||
		GCM==1260 ||GCM==1261 ||GCM==1262 ||GCM==1263){ 
		CO2S = 1;
		
		if((fp_s[23]=fopen("./data/image_a1b_fcrop.dat","rt"))==NULL){
			printf("NO DATA !!\n");
			exit(1);
		}
	}else if(GCM==1013 ||GCM==1014 ||GCM==1015 || GCM==1051 ||GCM==1061 ||GCM==1071 ||
		GCM==1081 ||GCM==1091 ||GCM==1111 ||GCM==1121 ||GCM==1135 ||GCM==1136 ||GCM==1137 ||
		GCM==1138 ||GCM==1139 ||GCM==1154 ||GCM==1155 ||GCM==1156 ||GCM==1171 ||GCM==1181 ||
		GCM==1191 ||GCM==1205 ||GCM==1206 ||GCM==1207 ||GCM==1208 ||GCM==1209 ||GCM==1227 ||
		GCM==1228 ||GCM==1229 ||GCM==1230 ||GCM==1231 ||GCM==1241 ||GCM==1252 ||GCM==1264 ||
		GCM==1265 ||GCM==1266 ||GCM==1267){ 
		CO2S = 4;

		if((fp_s[23]=fopen("./data/image_a2_fcrop.dat","rt"))==NULL){
			printf("NO DATA !!\n");
			exit(1);
		}
	}else if(GCM==1001 ||GCM==1016 ||GCM==1017 ||GCM==1018 ||GCM==1052 ||GCM==1062 ||GCM==1072 ||
		GCM==1082 ||GCM==1092 ||GCM==1102 ||GCM==1103 ||GCM==1112 ||GCM==1122 ||GCM==1140 ||
		GCM==1141 ||GCM==1142 ||GCM==1143 ||GCM==1144 ||GCM==1157 ||GCM==1158 ||GCM==1159 ||
		GCM==1163 ||GCM==1164 ||GCM==1165 ||GCM==1182 ||GCM==1192 ||GCM==1210 ||GCM==1211 ||
		GCM==1212 ||GCM==1213 ||GCM==1214 ||GCM==1232 ||GCM==1233 ||GCM==1234 ||GCM==1235 ||
		GCM==1236 ||GCM==1237 ||GCM==1238 ||GCM==1239 ||GCM==1242 ||GCM==1253){ 
		CO2S = 5;
		
		if((fp_s[23]=fopen("./data/image_b1_fcrop.dat","rt"))==NULL){
			printf("NO DATA !!\n");
			exit(1);
		}
	}
	/* 0: stable */
	/* 1: gradual rise (SRES A1 or A1B) */
	/* 2: gradual rise (SRES A1FI) */
	/* 3: gradual rise (SRES A1T) */
	/* 4: gradual rise (SRES A2) */
	/* 5: gradual rise (SRES B1) */
	/* 6: gradual rise (SRES B2) */
	/* 7: uniform rise */

	/***************************************************/
	if(GCM==0){
		GCM_R = 1;
		GCM_C = 1;
	}else if(GCM==1 || GCM==2 || GCM==3 || GCM==4 || GCM==5 || GCM==6){
		GCM_R = 32;
		GCM_C = 64;
	}else if(GCM==11 || GCM==12 || GCM==13 || GCM==14 || GCM==15 || GCM==16 || GCM==17 || GCM==18){
		GCM_R = 48;
		GCM_C = 96;
	}else if(GCM==21 || GCM==22 || GCM==23 || GCM==24){
		GCM_R = 73;
		GCM_C = 96;
	}else if(GCM==31 || GCM==32){
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM==41 || GCM==42 || GCM==43 || GCM==44){
		GCM_R = 56;
		GCM_C = 64;
	}else if(GCM==51 || GCM==52){
		GCM_R = 80;
		GCM_C = 96;
	}else if(GCM==61 || GCM==62 || GCM==63){
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM==71){
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM==1000 || GCM==1001){ /* MIROC-HIGH */
		GCM_R = 160;
		GCM_C = 320;
	}else if(GCM==1010 || GCM==1011 || GCM==1012 || GCM==1013 || GCM==1014 || 
			GCM==1015 || GCM==1016 || GCM==1017 || GCM==1018){ /* MIROC-MED */
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM==1050 || GCM==1051 ||GCM==1052){ /* BCCR */
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM==1061 || GCM==1062 ||GCM==1063){ /* INM-CM3 */
		GCM_R = 45;
		GCM_C = 72;
	}else if(GCM==1070 || GCM==1071 || GCM==1072){ /* IPSL */
		GCM_R = 72;
		GCM_C = 96;
	}else if(GCM==1080 || GCM==1081 || GCM==1082){ /* GFDL CM2-1 */
		GCM_R = 90;
		GCM_C = 144;
	}else if(GCM==1090 || GCM==1091 || GCM==1092){ /* HadCM3 */
		GCM_R = 73;
		GCM_C = 96;
	}else if(GCM==1100 || GCM==1101 || GCM==1102 || GCM==1103){ /* GISS AOM */
		GCM_R = 60;
		GCM_C = 90;
	}else if(GCM==1110 || GCM==1111 || GCM==1112){ /* CCCma-T63 */
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM==1120 || GCM==1121 || GCM==1122){ /* CSIRO */
		GCM_R = 96;
		GCM_C = 192;
	}else if(GCM==1130 || GCM==1131 || GCM==1132 || GCM==1133 || GCM==1134
			 || GCM==1135 || GCM==1136 || GCM==1137 || GCM==1138 || GCM==1139
			  || GCM==1140 || GCM==1141 || GCM==1142 || GCM==1143 || GCM==1144){ /* MRI */
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM==1150 ||GCM==1151 ||GCM==1152 ||GCM==1153 ||GCM==1154 ||
			GCM==1155 ||GCM==1156 ||GCM==1157 ||GCM==1158 ||GCM==1159){ /* MPI ECHAM */
		GCM_R = 96;
		GCM_C = 192;
	}else if(GCM==1160 || GCM==1161 || GCM==1162 || GCM==1163 || GCM==1164 || GCM==1165){ /* IAP */
		GCM_R = 60;
		GCM_C = 128;
	}else if(GCM==1170 || GCM==1171){ /* HadGEM */
		GCM_R = 145;
		GCM_C = 192;
	}else if(GCM==1180 || GCM==1181 || GCM==1182){ /* GFDL CM 2.0 */
		GCM_R = 90;
		GCM_C = 144;
	}else if(GCM==1190 || GCM==1191 || GCM==1192){ /* CNRM */
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM==1200 || GCM==1201 || GCM==1202 || GCM==1203 || GCM==1204
			 || GCM==1205 || GCM==1206 || GCM==1207 || GCM==1208 || GCM==1209
			 || GCM==1210 || GCM==1211 || GCM==1212 || GCM==1213 || GCM==1214){ /* CCC47 */
		GCM_R = 48;
		GCM_C = 96;
	}else if(GCM==1220 || GCM==1221 || GCM==1222 || GCM==1223 || GCM==1224
			 || GCM==1225 || GCM==1226 || GCM==1227 || GCM==1228 || GCM==1229
			  || GCM==1230 || GCM==1231 || GCM==1232 || GCM==1233 || GCM==1234
			   || GCM==1235 || GCM==1236 || GCM==1237 || GCM==1238 || GCM==1239){ /* NCAR CCSM */
		GCM_R = 128;
		GCM_C = 256;
	}else if(GCM==1240 || GCM==1241 || GCM==1242){ /* GISS E */
		GCM_R = 46;
		GCM_C = 72;
	}else if(GCM==1250 || GCM==1251 || GCM==1252 || GCM==1253){ /* GISS R */
		GCM_R = 46;
		GCM_C = 72;
	}else if(GCM==1260 || GCM==1261 || GCM==1262 || GCM==1263 || GCM==1264
			 || GCM==1265 || GCM==1266 || GCM==1267){ /* NCAR PCM */
		GCM_R = 64;
		GCM_C = 128;
	}
	
	printf("Open files...");
	/** open source files ***************************************/
	if( (fp_s[0]=fopen("./data/nat_reg_05.dat","rt"))==NULL ){  printf("No nat_reg_05.dat\n");  exit(1); }
	if( (fp_s[1]=fopen("./data/olson_biome3_05.dat","rt"))==NULL ){  printf("No olson_biome3_05.dat\n");  exit(1); }
	if( (fp_s[2]=fopen("./data/climate_05.dat","rt"))==NULL ){  printf("No climate_05.dat\n");  exit(1); }
	if( (fp_s[3]=fopen("./data/geo_prop_05.dat","rt"))==NULL ){  printf("No geo_prop_05.dat\n");  exit(1); }
	if( (fp_s[4]=fopen("./data/prec_cru.dat","rt"))==NULL ){  printf("No prec_cru.dat\n");  exit(1); }
	if( (fp_s[5]=fopen("./data/rivnum_trip_05.dat","rt"))==NULL ){  printf("No rivnum_trip_05.dat\n");  exit(1); }
	if( (fp_s[6]=fopen("./data/Albedo_soil_EX.dat","rt"))==NULL ){  printf("No Albedo_soil_EX.dat\n");  exit(1); }
	if( (fp_s[7]=fopen("./data/pcnt_orgmatter.dat","rt"))==NULL ){  printf("No pcnt_orgmatter.dat\n");  exit(1); }
	if( (fp_s[8]=fopen("./data/fls_slope.dat","rt"))==NULL ){  printf("No fls_slope.dat\n");  exit(1); }
	if( (fp_s[9]=fopen("./data/fk_erodibility.dat","rt"))==NULL ){  printf("No fk_erodibility.dat\n");  exit(1); }
	if( (fp_s[10]=fopen("./data/potveg_sage_05.dat","rt"))==NULL ){  printf("No potveg_sage_05.dat\n");  exit(1); }
	if( (fp_s[11]=fopen("./data/luc_1700-1992.dat","rt"))==NULL ){  printf("No luc_1700-1992.dat\n");  exit(1); }
	if( (fp_s[12]=fopen("./data/FAO_soil.textur.1nnegl.dat","rt"))==NULL ){  printf("No FAO_soil.textur.1nnegl.dat\n");  exit(1); }
	if( (fp_s[13]=fopen("./data/bulkdens.dat","rt"))==NULL ){  printf("No bulkdens.dat\n");  exit(1); }
	if( (fp_s[14]=fopen("./data/fieldcap.dat","rt"))==NULL ){  printf("No fieldcap.dat\n");  exit(1); }
	if( (fp_s[15]=fopen("./data/soil_field_cap0-150_30min.dat","rt"))==NULL ){  printf("No soil_field_cap0-150_30min.dat\n");  exit(1); }
	if( (fp_s[16]=fopen("./data/soil_field_cap0-30_30min.dat","rt"))==NULL ){  printf("No soil_field_cap0-30_30min.dat\n");  exit(1); }
	if( (fp_s[17]=fopen("./data/soil_bulk_dens0-30_30min.dat","rt"))==NULL ){  printf("No soil_bulk_dens0-30_30min.dat\n");  exit(1); }
	if( (fp_s[18]=fopen("./data/soil_clay_perc0-30_30min.dat","rt"))==NULL ){  printf("No soil_clay_perc0-30_30min.dat\n");  exit(1); }
	if( (fp_s[19]=fopen("./data/soil_sand_perc0-30_30min.dat","rt"))==NULL ){  printf("No soil_sand_perc0-30_30min.dat\n");  exit(1); }
	if( (fp_s[20]=fopen("./data/rastert_soilph1.dat","rt"))==NULL ){  printf("No rastert_soilph1.dat\n");  exit(1); }
	if( (fp_s[21]=fopen("./data/wetland_glwd_05.dat","rt"))==NULL ){  printf("No wetland_glwd_05.dat\n");  exit(1); }
	if( (fp_s[22]=fopen("./data/rice_sage_05.dat","rt"))==NULL ){  printf("No rice_sage_05.dat\n");  exit(1); }
	if( (fp_s[24]=fopen("./data/soil_n_1m.dat","rt"))==NULL ){  printf("No soil_n_1m.dat\n");  exit(1); }
	if( (fp_s[25]=fopen("./data/N-deposition.dat","rt"))==NULL ){  printf("No N-deposition.dat\n");  exit(1); }

	fp_c[0]=fopen("./data/cru2_cld_0110.dat","rt");
	fp_c[1]=fopen("./data/cru2_cld_1120.dat","rt");
	fp_c[2]=fopen("./data/cru2_cld_2130.dat","rt");
	fp_c[3]=fopen("./data/cru2_cld_3140.dat","rt");
	fp_c[4]=fopen("./data/cru2_cld_4150.dat","rt");
	fp_c[5]=fopen("./data/cru2_cld_5160.dat","rt");
	fp_c[6]=fopen("./data/cru2_cld_6170.dat","rt");
	fp_c[7]=fopen("./data/cru2_cld_7180.dat","rt");
	fp_c[8]=fopen("./data/cru2_cld_8190.dat","rt");
	fp_c[9]=fopen("./data/cru2_cld_9100.dat","rt");

	fp_c[10]=fopen("./data/cru2_pre_0110.dat","rt");
	fp_c[11]=fopen("./data/cru2_pre_1120.dat","rt");
	fp_c[12]=fopen("./data/cru2_pre_2130.dat","rt");
	fp_c[13]=fopen("./data/cru2_pre_3140.dat","rt");
	fp_c[14]=fopen("./data/cru2_pre_4150.dat","rt");
	fp_c[15]=fopen("./data/cru2_pre_5160.dat","rt");
	fp_c[16]=fopen("./data/cru2_pre_6170.dat","rt");
	fp_c[17]=fopen("./data/cru2_pre_7180.dat","rt");
	fp_c[18]=fopen("./data/cru2_pre_8190.dat","rt");
	fp_c[19]=fopen("./data/cru2_pre_9100.dat","rt");

	fp_c[20]=fopen("./data/cru2_tmp_0110.dat","rt");
	fp_c[21]=fopen("./data/cru2_tmp_1120.dat","rt");
	fp_c[22]=fopen("./data/cru2_tmp_2130.dat","rt");
	fp_c[23]=fopen("./data/cru2_tmp_3140.dat","rt");
	fp_c[24]=fopen("./data/cru2_tmp_4150.dat","rt");
	fp_c[25]=fopen("./data/cru2_tmp_5160.dat","rt");
	fp_c[26]=fopen("./data/cru2_tmp_6170.dat","rt");
	fp_c[27]=fopen("./data/cru2_tmp_7180.dat","rt");
	fp_c[28]=fopen("./data/cru2_tmp_8190.dat","rt");
	fp_c[29]=fopen("./data/cru2_tmp_9100.dat","rt");

	fp_c[30]=fopen("./data/cru2_vap_0110.dat","rt");
	fp_c[31]=fopen("./data/cru2_vap_1120.dat","rt");
	fp_c[32]=fopen("./data/cru2_vap_2130.dat","rt");
	fp_c[33]=fopen("./data/cru2_vap_3140.dat","rt");
	fp_c[34]=fopen("./data/cru2_vap_4150.dat","rt");
	fp_c[35]=fopen("./data/cru2_vap_5160.dat","rt");
	fp_c[36]=fopen("./data/cru2_vap_6170.dat","rt");
	fp_c[37]=fopen("./data/cru2_vap_7180.dat","rt");
	fp_c[38]=fopen("./data/cru2_vap_8190.dat","rt");
	fp_c[39]=fopen("./data/cru2_vap_9100.dat","rt");

	/** open result file ******************************************/
	switch(GCM){
		case 1000:	strcpy(s_case,"MHA11_"); break;
		case 1001:	strcpy(s_case,"MHB11_"); break;

		case 1010:	strcpy(s_case,"MMA11_"); break;
		case 1011:	strcpy(s_case,"MMA12_"); break;
		case 1012:	strcpy(s_case,"MMA13_"); break;
		case 1013:	strcpy(s_case,"MMA21_"); break;
		case 1014:	strcpy(s_case,"MMA22_"); break;
		case 1015:	strcpy(s_case,"MMA23_"); break;
		case 1016:	strcpy(s_case,"MMA21_"); break;
		case 1017:	strcpy(s_case,"MMA22_"); break;
		case 1018:	strcpy(s_case,"MMA23_"); break;

		case 1050:	strcpy(s_case,"BCA11_"); break;
		case 1051:	strcpy(s_case,"BCA21_"); break;
		case 1052:	strcpy(s_case,"BCB11_"); break;

		case 1060:	strcpy(s_case,"INA11_"); break;
		case 1061:	strcpy(s_case,"INA21_"); break;
		case 1062:	strcpy(s_case,"INB11_"); break;

		case 1070:	strcpy(s_case,"IPA11_"); break;
		case 1071:	strcpy(s_case,"IPA21_"); break;
		case 1072:	strcpy(s_case,"IPB11_"); break;

		case 1080:	strcpy(s_case,"GF1A11_"); break;
		case 1081:	strcpy(s_case,"GF1A21_"); break;
		case 1082:	strcpy(s_case,"GF1B11_"); break;

		case 1090:	strcpy(s_case,"HCA11_"); break;
		case 1091:	strcpy(s_case,"HCA21_"); break;
		case 1092:	strcpy(s_case,"HCB11_"); break;

		case 1100:	strcpy(s_case,"GSAA11_"); break;
		case 1101:	strcpy(s_case,"GSAA12_"); break;
		case 1102:	strcpy(s_case,"GSAA21_"); break;
		case 1103:	strcpy(s_case,"GSAB11_"); break;

		case 1110:	strcpy(s_case,"CC6A11_"); break;
		case 1111:	strcpy(s_case,"CC6A21_"); break;
		case 1112:	strcpy(s_case,"CC6B11_"); break;

		case 1120:	strcpy(s_case,"CSA11_"); break;
		case 1121:	strcpy(s_case,"CSA21_"); break;
		case 1122:	strcpy(s_case,"CSB11_"); break;

		case 1130:	strcpy(s_case,"MRA11_"); break;
		case 1131:	strcpy(s_case,"MRA12_"); break;
		case 1132:	strcpy(s_case,"MRA13_"); break;
		case 1133:	strcpy(s_case,"MRA14_"); break;
		case 1134:	strcpy(s_case,"MRA15_"); break;
		case 1135:	strcpy(s_case,"MRA21_"); break;
		case 1136:	strcpy(s_case,"MRA22_"); break;
		case 1137:	strcpy(s_case,"MRA23_"); break;
		case 1138:	strcpy(s_case,"MRA24_"); break;
		case 1139:	strcpy(s_case,"MRA25_"); break;
		case 1140:	strcpy(s_case,"MRB11_"); break;
		case 1141:	strcpy(s_case,"MRB12_"); break;
		case 1142:	strcpy(s_case,"MRB13_"); break;
		case 1143:	strcpy(s_case,"MRB14_"); break;
		case 1144:	strcpy(s_case,"MRB15_"); break;

		case 1150:	strcpy(s_case,"MPA11_"); break;
		case 1151:	strcpy(s_case,"MPA12_"); break;
		case 1152:	strcpy(s_case,"MPA13_"); break;
		case 1153:	strcpy(s_case,"MPA14_"); break;
		case 1154:	strcpy(s_case,"MPA21_"); break;
		case 1155:	strcpy(s_case,"MPA22_"); break;
		case 1156:	strcpy(s_case,"MPA23_"); break;
		case 1157:	strcpy(s_case,"MPB11_"); break;
		case 1158:	strcpy(s_case,"MPB12_"); break;
		case 1159:	strcpy(s_case,"MPB13_"); break;

		case 1160:	strcpy(s_case,"IAA11_"); break;
		case 1161:	strcpy(s_case,"IAA12_"); break;
		case 1162:	strcpy(s_case,"IAA13_"); break;
		case 1163:	strcpy(s_case,"IAB11_"); break;
		case 1164:	strcpy(s_case,"IAB12_"); break;
		case 1165:	strcpy(s_case,"IAB13_"); break;

		case 1170:	strcpy(s_case,"HGA11_"); break;
		case 1171:	strcpy(s_case,"HGA21_"); break;

		case 1180:	strcpy(s_case,"GF0A11_"); break;
		case 1181:	strcpy(s_case,"GF0A21_"); break;
		case 1182:	strcpy(s_case,"GF0B11_"); break;

		case 1190:	strcpy(s_case,"CNA11_"); break;
		case 1191:	strcpy(s_case,"CNA21_"); break;
		case 1192:	strcpy(s_case,"CNB11_"); break;

		case 1200:	strcpy(s_case,"CC4A11_"); break;
		case 1201:	strcpy(s_case,"CC4A12_"); break;
		case 1202:	strcpy(s_case,"CC4A13_"); break;
		case 1203:	strcpy(s_case,"CC4A14_"); break;
		case 1204:	strcpy(s_case,"CC4A15_"); break;
		case 1205:	strcpy(s_case,"CC4A21_"); break;
		case 1206:	strcpy(s_case,"CC4A22_"); break;
		case 1207:	strcpy(s_case,"CC4A23_"); break;
		case 1208:	strcpy(s_case,"CC4A24_"); break;
		case 1209:	strcpy(s_case,"CC4A25_"); break;
		case 1210:	strcpy(s_case,"CC4B11_"); break;
		case 1211:	strcpy(s_case,"CC4B12_"); break;
		case 1212:	strcpy(s_case,"CC4B13_"); break;
		case 1213:	strcpy(s_case,"CC4B14_"); break;
		case 1214:	strcpy(s_case,"CC4B15_"); break;

		case 1220:	strcpy(s_case,"NCA11_"); break;
		case 1221:	strcpy(s_case,"NCA12_"); break;
		case 1222:	strcpy(s_case,"NCA13_"); break;
		case 1223:	strcpy(s_case,"NCA15_"); break;
		case 1224:	strcpy(s_case,"NCA16_"); break;
		case 1225:	strcpy(s_case,"NCA17_"); break;
		case 1226:	strcpy(s_case,"NCA19_"); break;
		case 1227:	strcpy(s_case,"NCA21_"); break;
		case 1228:	strcpy(s_case,"NCA22_"); break;
		case 1229:	strcpy(s_case,"NCA23_"); break;
		case 1230:	strcpy(s_case,"NCA24_"); break;
		case 1231:	strcpy(s_case,"NCA25_"); break;
		case 1232:	strcpy(s_case,"NCB11_"); break;
		case 1233:	strcpy(s_case,"NCB12_"); break;
		case 1234:	strcpy(s_case,"NCB13_"); break;
		case 1235:	strcpy(s_case,"NCB14_"); break;
		case 1236:	strcpy(s_case,"NCB15_"); break;
		case 1237:	strcpy(s_case,"NCB16_"); break;
		case 1238:	strcpy(s_case,"NCB17_"); break;
		case 1239:	strcpy(s_case,"NCB19_"); break;

		case 1240:	strcpy(s_case,"GSEA11_"); break;
		case 1241:	strcpy(s_case,"GSEA21_"); break;
		case 1242:	strcpy(s_case,"GSEB11_"); break;

		case 1250:	strcpy(s_case,"GSRA11_"); break;
		case 1251:	strcpy(s_case,"GSRA12_"); break;
		case 1252:	strcpy(s_case,"GSRA21_"); break;
		case 1253:	strcpy(s_case,"GSRB11_"); break;

		case 1260:	strcpy(s_case,"GSRA11_"); break;
		case 1261:	strcpy(s_case,"GSRA12_"); break;
		case 1262:	strcpy(s_case,"GSRA13_"); break;
		case 1263:	strcpy(s_case,"GSRA14_"); break;
		case 1264:	strcpy(s_case,"GSRA21_"); break;
		case 1265:	strcpy(s_case,"GSRA22_"); break;
		case 1266:	strcpy(s_case,"GSRA23_"); break;
		case 1267:	strcpy(s_case,"GSRA24_"); break;
		
		default:
			printf("NULL CASE !!!");
		break;
		/* fl1 */
		/* fl2 */
		/* fl3 */
		/* nl */
		/* rcf */
		/* apf */
	}
	
	strcpy(s_carbon, s_case);
	strcpy(s_nitrogen, s_case);
	strcpy(s_erosion, s_case);
	strcpy(s_ghg, s_case);
	strcpy(s_bioburn, s_case);
	strcpy(s_voc, s_case);
	
	strcat(s_carbon, "crbn_");
	strcat(s_nitrogen, "nitr_");
	strcat(s_erosion, "ersn_");
	strcat(s_ghg, "ghg_");
	strcat(s_bioburn, "bb_");
	strcat(s_voc, "voc_");
	
	/***********/
	strcpy(s_date, argv[2]);
	strcat(s_date, "_");
	
	for(f=1;f<=18;f++){
	}
			
	printf("done\n");
	
	/******************************************************************************/
	printf("Initialize simulation...");

	glandarea = 0.0;
	for(f=0;f<201;f++){
		h_tmp[f] = h_pre[f] = h_dswr[f] = h_aet[f] = h_rof[f] = 0.0;
		h_gpp[f] = h_npp[f] = h_nep[f] = h_plant[f] = h_soil[f] = 0.0;
		h_sr[f] = h_ersn_c[f] = 0.0;
		
		h_agrarea[f] = h_luc[f] = 0.0;
		h_burnt_area[f] = 0.0;
		h_bioburn_co2[f] = h_bioburn_co[f] = h_bioburn_ch4[f] = 0.0;
		h_bioburn_nmhc[f] = h_bioburn_oc[f] = h_bioburn_bc[f] = 0.0;
		h_bioburn_nox[f] = h_bioburn_so2[f] = h_bioburn_pm25[f] = 0.0;
		h_bioburn_tpm[f] = h_bioburn_tec[f] = 0.0;
		
		h_ch4ox1[f] = h_ch4ox2[f] = h_ch4ox3[f] = 0.0;
		h_ch4emit_cao_paddy[f] = h_ch4emit_cao_wetland[f] = 0.0;
		h_n2o_emit_ngas[f] = h_n2_emit_ngas[f] = 0.0;
		h_n2o_emit_casa[f] = h_no_emit_casa[f] = h_n2_emit_casa[f] = 0.0;
		h_nh3_emit[f] = 0.0;

		h_voc_isopr_g97[f] = h_voc_monotrp_g97[f] = h_voc_methanl_g97[f] = 0.0;
		h_voc_acetone_g97[f] = h_voc_actaldhd_g97[f] = h_voc_frmardhd_g97[f] = 0.0;
		h_voc_formacd_g97[f] = h_voc_acetacd_g97[f] = h_voc_co_g97[f] = 0.0;
	}
	for(f=0;f<12;f++){
		m_ch4ox1[f] = m_ch4ox2[f] = m_ch4ox3[f] = 0.0;
		m_bioburn_co2[f] = m_bioburn_ch4[f] = m_bioburn_co[f] = 0.0;
		m_bioburn_nmhc[f] = m_bioburn_oc[f] = m_bioburn_bc[f] = 0.0;
	}
	for(f=0;f<34;f++){
		v_area[f] = 0.0;
		v_gpp[f] = v_npp[f] = v_nep[f] = 0.0;
		v_lai[f] = v_fol[f] = v_stm[f] = v_rot[f] = v_ltr[f] = v_msl[f] = 0.0;
	}
	
	/** initialize configuration **/
	initSim(&grid);
	printf("done\n");
	
	/*** read GCM climate scenario ***/
	printf("Reading GCM climate projection...");
	read_gcm_clim(&grid); /* */
	printf("done\n");
	
	/***************************************************************************************/
	/* latitude roop */
	printf("Start simulation...\n");
	for(f=0;f<360;f++){
		zone = f/20 + 1;
		
		if(f%20 == 0){
			/* carbon */
			strcpy(filename, s_date);
			strcat(filename, s_carbon);
			if(zone<=9){
				strcat(filename, "0");
			}
			snprintf(num, 4, "%ld", zone);
			strcat(filename, num);
			strcat(filename, ".dat");
			fp_carbon = fopen(filename,"wt");
			
			/* nitrogen */
			strcpy(filename, s_date);
			strcat(filename, s_nitrogen);
			if(zone<=9){
				strcat(filename, "0");
			}
			snprintf(num, 4, "%ld", zone);
			strcat(filename, num);
			strcat(filename, ".dat");
			fp_nitrogen = fopen(filename,"wt");
			
			/* erosion */
			strcpy(filename, s_date);
			strcat(filename, s_erosion);
			if(zone<=9){
				strcat(filename, "0");
			}
			snprintf(num, 4, "%ld", zone);
			strcat(filename, num);
			strcat(filename, ".dat");
			fp_ersn = fopen(filename,"wt");
			
			/* GHG */
			strcpy(filename, s_date);
			strcat(filename, s_ghg);
			if(zone<=9){
				strcat(filename, "0");
			}
			snprintf(num, 4, "%ld", zone);
			strcat(filename, num);
			strcat(filename, ".dat");
			fp_ghg = fopen(filename,"wt");

			/* Biomass burning */
			strcpy(filename, s_date);
			strcat(filename, s_bioburn);
			if(zone<=9){
				strcat(filename, "0");
			}
			snprintf(num, 4, "%ld", zone);
			strcat(filename, num);
			strcat(filename, ".dat");
			fp_bioburn = fopen(filename,"wt");

			/* VOC */
			strcpy(filename, s_date);
			strcat(filename, s_voc);
			if(zone<=9){
				strcat(filename, "0");
			}
			snprintf(num, 4, "%ld", zone);
			strcat(filename, num);
			strcat(filename, ".dat");
			fp_voc = fopen(filename,"wt");
		}

		/* grid area */
		grid.area = grid_area(90.0-0.5*(double)f, 89.5-0.5*(double)f, 0.0, 0.5);
		
		/* longitude roop */
		for(g=0;g<720;g++){
			/** roop of longitudinal rows, from W to E at the 0.5deg resolution **/	
			grid.row = f; 
			grid.col = g;	
			
			/* initialize grid condition */
			initG(fp_s, &grid);
			
			/* read CRU climate data */
			read_cru_clim(fp_c, &grid);
			
			printf("%3ld %3ld: %7.2lf %7.2lf: %2ld: %1ld :area %.1lf\n", 
				grid.row, grid.col, grid.lat, grid.lon, grid.veg_olson, grid.cru_exist, grid.area); /* */
			
			fprintf(fp_carbon,"%ld %ld %ld\n", grid.row, grid.col, grid.veg_olson); 
			fprintf(fp_nitrogen,"%ld %ld %ld\n", grid.row, grid.col, grid.veg_olson); 
			fprintf(fp_ersn,"%ld %ld %ld %ld %ld %lf\n", 
				grid.row, grid.col, grid.veg_olson, grid.veg_sage, grid.rvbasin, grid.f_crop_trend); 
			fprintf(fp_ghg,"%ld %ld %ld\n", grid.row, grid.col, grid.veg_olson); 
			fprintf(fp_bioburn,"%ld %ld %ld\n", grid.row, grid.col, grid.veg_olson); 
			fprintf(fp_voc,"%ld %ld %ld\n", grid.row, grid.col, grid.veg_olson); 
			
			if(grid.veg_olson!=0 && grid.veg_olson!=33 && (g+0)%1==0){ /* */ /* && (g+0)%10==0 */
			/* if(f==62 && g==377){ */
				grid.nnn++;
				glandarea += grid.area;
				v_area[grid.veg_olson] += grid.area;
				
				/*************************************/
				/** clear all parameters **/
				clear(&grid, &loct, &echar, &mass, &flux);
				
				/** initialize vegatation and soil conditions **/ 
				initVS(&grid, &loct, &mass, &flux, &echar);
				
				/** initialize climate conditions **/
				initC(&grid, &loct);
				
				/** initialize location conditions **/
				initL(&grid, &loct, &mass, &flux, &echar);
				
				/** initialize stable carbon isotope **/
				init_d13c(&grid, &flux, &echar, &mass);
					
				/** stabilization roop **/
				cal_stable(&grid, &loct, &echar, &mass, &flux, 
						fp_carbon, fp_nitrogen, fp_ersn, fp_ghg, fp_bioburn, fp_voc);
				
				/* printf("Lan use change: %.3lf<-%.3lf Emit: %.3lf\n", grid.f_crop_p, grid.f_crop, (flux.lu_conv+flux.lu_ten+flux.lu_hund)); */
			
				/** dynamic roop **/
				/* 1901-2000 */
				cal_cruclim(&grid, &loct, &echar, &mass, &flux, 
						fp_carbon, fp_nitrogen, fp_ersn, fp_ghg, fp_bioburn, fp_voc); 
				/* printf("Lan use change: %.3lf<-%.3lf Emit: %.3lf\n", grid.f_crop_p, grid.f_crop, (flux.lu_conv+flux.lu_ten+flux.lu_hund)); */
			
				/* 2001-2100 */
				cal_gcmclim2(&grid, &loct, &echar, &mass, &flux, 
						fp_carbon, fp_nitrogen, fp_ersn, fp_ghg, fp_bioburn, fp_voc);
				
				/** outputs **/
				screenshow(&grid, &loct, &mass, &flux, &echar);

				/* printf("Lan use change: %.3lf<-%.3lf Emit: %.3lf\n", grid.f_crop_p, grid.f_crop, (flux.lu_conv+flux.lu_ten+flux.lu_hund)); */
				
				/* printf("CUM LUC: %.1lf %.1lf %.1lf\n\n", h_luc[0], h_luc[100], h_luc[200]); */
				
				printf("\n");
			}
		}
		
		strcpy(filename, s_date);
		strcat(filename, s_case);
		strcat(filename, "glsum.dat");
		
		fp_glsum = fopen(filename,"wt");
		fprintf(fp_glsum,"%lf\n", glandarea);
		for(h=0;h<201;h++){
			fprintf(fp_glsum,"%ld ", h+1900);
			fprintf(fp_glsum,"%lf %lf %lf %lf %lf ", h_tmp[h], h_pre[h], h_dswr[h], h_aet[h], h_rof[h]);
			fprintf(fp_glsum,"%lf %lf %lf ", h_gpp[h], h_npp[h], h_nep[h]);
			fprintf(fp_glsum,"%lf %lf %lf ", h_plant[h], h_soil[h], h_sr[h]);
			
			fprintf(fp_glsum,"%lf ", h_ersn_c[h]);
			
			fprintf(fp_glsum,"%lf %lf ", h_agrarea[h], h_luc[h]);
			
			fprintf(fp_glsum,"%lf ", h_burnt_area[h]);
			fprintf(fp_glsum,"%lf %lf %lf ", h_bioburn_co2[h], h_bioburn_co[h], h_bioburn_ch4[h]);
			fprintf(fp_glsum,"%lf %lf %lf ", h_bioburn_nmhc[h], h_bioburn_oc[h], h_bioburn_bc[h]);
			fprintf(fp_glsum,"%lf %lf %lf ", h_bioburn_nox[h], h_bioburn_so2[h], h_bioburn_pm25[h]);
			fprintf(fp_glsum,"%lf %lf ", h_bioburn_tpm[h], h_bioburn_tec[h]);
			
			fprintf(fp_glsum,"%lf %lf ", h_n2o_emit_ngas[h], h_n2_emit_ngas[h]);
			fprintf(fp_glsum,"%lf %lf %lf ", h_n2o_emit_casa[h], h_n2_emit_casa[h], h_no_emit_casa[h]);
			fprintf(fp_glsum,"%lf %lf %lf ", h_ch4ox1[h], h_ch4ox2[h], h_ch4ox3[h]);
			fprintf(fp_glsum,"%lf %lf ", h_ch4emit_cao_paddy[h], h_ch4emit_cao_wetland[h]);
			fprintf(fp_glsum,"%lf ", h_nh3_emit[h]);
			
			fprintf(fp_glsum,"%lf %lf %lf ", h_voc_isopr_g97[h], h_voc_monotrp_g97[h], h_voc_methanl_g97[h]);
			fprintf(fp_glsum,"%lf %lf %lf ", h_voc_acetone_g97[h], h_voc_actaldhd_g97[h], h_voc_frmardhd_g97[h]);
			fprintf(fp_glsum,"%lf %lf %lf ", h_voc_formacd_g97[h], h_voc_acetacd_g97[h], h_voc_co_g97[h]);
			fprintf(fp_glsum,"\n");
		}
		fprintf(fp_glsum,"\n");
		for(h=0;h<12;h++){
			fprintf(fp_glsum,"%ld %lf %lf %lf ", h, m_ch4ox1[h], m_ch4ox2[h], m_ch4ox3[h]);
			fprintf(fp_glsum,"%lf %lf %lf ", m_bioburn_co2[h], m_bioburn_co[h], m_bioburn_ch4[h]);
			fprintf(fp_glsum,"%lf %lf %lf ", m_bioburn_nmhc[h], m_bioburn_oc[h], m_bioburn_bc[h]);
			fprintf(fp_glsum,"\n");
		}
		fprintf(fp_glsum,"\n");
		for(h=0;h<34;h++){
			fprintf(fp_glsum,"%lf ", v_area[h]);
			fprintf(fp_glsum,"%lf %lf %lf ", v_gpp[h], v_npp[h], v_nep[h]);
			fprintf(fp_glsum,"%lf %lf %lf %lf %lf %lf\n", v_lai[h], v_fol[h], v_stm[h], v_rot[h], v_ltr[h], v_msl[h]);
		}
		fclose(fp_glsum);
		
		if(f%20 == 19){
			fclose(fp_carbon); 
			fclose(fp_nitrogen); 
			fclose(fp_ersn);
			fclose(fp_ghg);
			fclose(fp_bioburn);
		}
	}
	
	/* close files */
	for(h=0;h<26;h++){
		fclose(fp_s[h]); 
	}
	for(h=0;h<40;h++){
		fclose(fp_c[h]); 
	}
	
	printf("Simulation ended\n");
	
	return 0;
}
