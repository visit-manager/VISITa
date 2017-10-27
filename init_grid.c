/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* Revised August 15, 2007 by A.Ito */
/* Revised August 19, 2007 by A.Ito */
/* Revised July 1, 2008 by A.Ito */
/* Separated February 17, 2009 by A.Ito */

/* initialize environmental characteristics of each grid */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* initial carbon pool size */
#define INT_C 0.01

/* initialize grid conditions **********************************************/
void f_init_grid(
	FILE *fp_s[IFILEN], 
	struct Grid *grid
){
	long e, h, g, country, region, aaa;
	double tmp_sfc, tmp_2m,tmp10_soil, tmp200_soil, dswrf_toa, dswrf_sfc, tcdc_clm;
	double prate_sfc, spfh_2m, soilw10, soilw200, ugrd_10m, vgrd_10m;
	double geo_prop, crit_tension;
	double lat, lon, total, wetland, lake, paddy, ddummy;
    double nfert_m, nfert_r, nfert_sw, nfert_ww;
    float rfdat[ASTEP], is2bdat[DL_NINPUT];
	
	/* grid latitudes of CSIRO AOGCM */
	double csiro_lat[56]={	
		 87.5613, 84.4022, 81.2245, 78.0425, 74.8590, -71.6747, 68.4899, 65.3049,
		 62.1197, 58.9343, 55.7489, 52.5634, 49.3779, 46.1924, 43.0068, 39.8211,
		 36.6355, 33.4498, 30.2642, 27.0785, 23.8928, 20.7071, 17.5214, 14.3357,
		 11.1500, 7.9643, 4.7786, 1.5929, -1.5929, -4.7786, -7.9643, -11.1500,
		 -14.3357, -17.5214, -20.7071, -23.8928, -27.0785, -30.2642, -33.4498, -36.6355,
		 -39.8211, -43.0068, -46.1924, -49.3779, -52.5634, -55.7489, -58.9343, -62.1197,
		 -65.3049, -68.4899, -71.6747, -74.8590, -78.0425, -81.2245, -84.4022, -87.5613
	}; 
	
	/* average field capacity (150cm) for Olson's biomes */
	double olveg_fc150[34] = {494.0844151, 592.0171207, 570.3302919, 551.9298516, 496.7930014,
		503.2473362, 509.419137, 472.4316637, 449.5906278, 409.8050576, 
		476.5335304, 435.8760416, 457.3711791, 531.5688406, 542.1963318, 
		502.4405676, 563.1641823, 503.0499857, 470.2081472, 504.7536108, 
		323.9272001, 415.1649954, 414.8235868, 531.6733625, 512.1259238, 
		506.076172, 461.7199239, 416.3817974, 347.099887, 551.8461504, 
		513.0419042, 509.2255222, 480.3831019, 245.5416336};
	/* average field capacity (30cm) for Olson's biomes */
	double olveg_fc30[34] = {128.0402598, 132.4363793, 134.744799, 128.4586405, 126.8214346, 
		130.0745942, 122.2639456, 125.8317647, 121.7714142, 117.3167714, 
		135.7699352, 130.4578967, 134.224135, 127.981914, 130.0432158, 
		123.7658455, 128.9114392, 122.3661531, 119.9684813, 123.5856945, 
		97.54575867, 125.0874554, 127.117925, 127.9560129, 137.5965516, 
		128.8016844, 119.5692298, 111.7970624, 109.7038513, 130.2056391, 
		128.7846698, 127.5222297, 121.7155247, 72.91636457};
	/* average bulk denisty for Olson's biomes */
	double olveg_bd[34] = {1.319663772, 1.337524685, 1.320005138, 1.425940593, 1.34650667, 
		1.371084061, 1.521147103, 1.310233567, 1.253574467, 1.124569698, 
		1.261084471, 1.141404339, 1.251011779, 1.383010262, 1.377164801, 
		1.474410875, 1.428756115, 1.488748801, 1.36807761, 1.458160626, 
		1.330929888, 1.283214093, 1.18470691, 1.3479263, 1.001494032, 
		1.349477665, 1.454684714, 1.538445206, 1.55359586, 1.426448633, 
		1.353602298, 1.403267173, 1.476433092, 0.899680244};
	/* average percentage sand fraction for Olson's biomes */
	double olveg_psand[34] = {46.52418085, 44.63135532, 44.83252566, 49.83793993, 44.40018048, 
		43.08494079, 53.44414746, 44.75929143, 37.19541727, 41.43155107, 
		36.4234834, 44.66832619, 38.57833011, 43.28312248, 41.98256651, 
		53.82325761, 50.19433257, 50.74725112, 49.92285997, 50.09217327, 
		41.71033249, 50.45165556, 46.78736248, 44.00399316, 33.84734969, 
		47.25795129, 50.13222498, 55.09942128, 68.81224605, 42.13913353, 
		38.83132202, 42.81201712, 45.56163789, 82.0};
	/* average percentage clay fraction for Olson's biomes */
	double olveg_pclay[34] = {24.21658857, 32.97702758, 29.69877674, 28.76695918, 23.7320935, 
		25.41978957, 23.75750454, 21.6693192, 23.98205022, 23.58236503, 
		25.7041984, 21.80942512, 23.58494412, 24.89991201, 28.20079033, 
		23.55827864, 28.21654069, 21.03711365, 21.2024205, 23.88689401, 
		28.05936266, 19.19199764, 21.56596459, 29.24609601, 27.77211489, 
		23.92326765, 20.40315154, 17.29372406, 12.1077655, 29.1624041, 
		24.38314974, 26.76259603, 23.10797052, 8.0};
	/* average soil pH for Olson's biomes */
	double olveg_ph[34] = {6.112273367, 5.278037138, 5.728076586, 6.052224862, 5.925479094, 
		6.333681185, 7.035558476, 5.977314591, 5.665770821, 5.277109702, 
		5.683146729, 5.525097779, 5.745947716, 5.880433971, 6.278726573, 
		6.841840519, 5.934110876, 7.229348227, 6.176573791, 6.873838987, 
		6.184740765, 6.205232922, 5.915702753, 6.095381045, 5.481843759, 
		6.317546811, 7.379898174, 7.438948581, 7.340133509, 6.045485427, 
		6.722120259, 6.742288787, 7.267506174, 5.759382756};
	double bulkdens, field_cap;
	double x, xx, y, yy, xy;
	
    /* added: A. Ito (with Hamada-san's comment) 2012/01/30 */
    grid->phase = 0;
    
	/* grid latitude and longitude ******************/
	grid->lat = 89.75 - 0.5*(double)grid->row;
	grid->lon = -179.75 + 0.5*(double)grid->col;
	
	/* region ID by F.Giorgi: added by A.Ito (2009/07/12) */
	grid->reg_g = region_giorgi(grid->lat, grid->lon);
	
	/* grid latitude and longitude corresponding to NCEP/NCAR reanalysis data */
	grid->ncep_lat = (long)((89.75-grid->lat)/(180.0/94.0));
	if(grid->lon>0.0){
		grid->ncep_lon = (long)(grid->lon/(360.0/192.0));
	}else{
		grid->ncep_lon = (long)(96.0 + (grid->lon+180.0)/(360.0/192.0));
	}
	
	/* GCM grid *************************************/
	grid->gcm_row = 0;
	grid->gcm_col = 0;
    
	/* AR3 ****************/
	if(SCENARIO_ID==1 || SCENARIO_ID==2 || SCENARIO_ID==3 || SCENARIO_ID==4
        || SCENARIO_ID==5 || SCENARIO_ID==6){ /* CCSR/NIES */
		grid->gcm_row = grid->row/11.25;	
		grid->gcm_col = grid->col/11.25;
	}else if(SCENARIO_ID==11 || SCENARIO_ID==12 || SCENARIO_ID==13 || SCENARIO_ID==14 || SCENARIO_ID==15 ||
            SCENARIO_ID==16 || SCENARIO_ID==17 || SCENARIO_ID==18){ /* CCCma */
		grid->gcm_row = grid->row/7.5;
		grid->gcm_col = grid->col/7.5;
	}else if(SCENARIO_ID==21 || SCENARIO_ID==22 || SCENARIO_ID==23 || SCENARIO_ID==24){ /* HadCM3 */
		if(grid->lat>=88.75){
			grid->gcm_row = 0;
		}else if(grid->lat<=-88.75){
			grid->gcm_row = 72;
		}else{
			grid->gcm_row = (long)((88.75-grid->lat)/2.5);
		}
		grid->gcm_col = grid->col/7.5;
	}else if(SCENARIO_ID==31 || SCENARIO_ID==32){ /* ECHAM */
		grid->gcm_row = grid->row/5.625;	
		grid->gcm_col = grid->col/5.625;
	}else if(SCENARIO_ID==41 || SCENARIO_ID==42 || SCENARIO_ID==43 || SCENARIO_ID==44){ /* CSIRO */
		for(e=0;e<55;e++){
			if(grid->lat>(csiro_lat[e]+csiro_lat[e+1])/2.0){
				grid->gcm_row = e;
			}
		}
		grid->gcm_col = grid->col/11.25;
	}else if(SCENARIO_ID==51 || SCENARIO_ID==52){ /* GFDL */
		grid->gcm_row = grid->row/4.5;
		grid->gcm_col = grid->col/7.5;
	}else if(SCENARIO_ID==61 || SCENARIO_ID==62 || SCENARIO_ID==63){ /* NCAR-PCM */
		grid->gcm_row = grid->row/5.625;	
		grid->gcm_col = grid->col/5.625;
	}else if(SCENARIO_ID==71){ /* NCAR-CSM */
		grid->gcm_row = grid->row/5.625;	
		grid->gcm_col = grid->col/5.625;
	}
	
	/* AR4 ****************/
	if(SCENARIO_ID==1000 || SCENARIO_ID==1001){ /* MIROC-HIGH */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1010 || SCENARIO_ID==1011 || SCENARIO_ID==1012 ||
            SCENARIO_ID==1013 || SCENARIO_ID==1014 || SCENARIO_ID==1015 ||
            SCENARIO_ID==1016 || SCENARIO_ID==1017 || SCENARIO_ID==1018){ /* MIROC-MED */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1050 || SCENARIO_ID==1051 ||SCENARIO_ID==1052){ /* BCCR */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1060 || SCENARIO_ID==1061 ||SCENARIO_ID==1062){ /* INM-CM3 : error-fixed 2009/04/26 A.Ito */
		if(grid->lat>=88.0){
			grid->gcm_row = 0;
		}else if(grid->lat<=-88.0){
			grid->gcm_row = 44;
		}else{
			grid->gcm_row = (long)((88.0-grid->lat)/4.0);
		}
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1070 || SCENARIO_ID==1071 || SCENARIO_ID==1072){ /* IPSL */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1080 || SCENARIO_ID==1081 || SCENARIO_ID==1082){ /* GFDL CM2-1 */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1090 || SCENARIO_ID==1091 || SCENARIO_ID==1092){ /* HadCM3 */
		if(grid->lat>=88.75){
			grid->gcm_row = 0;
		}else if(grid->lat<=-88.75){
			grid->gcm_row = 72;
		}else{
			grid->gcm_row = (long)((88.75-grid->lat)/2.5);
		}
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1100 || SCENARIO_ID==1101 || SCENARIO_ID==1102 || SCENARIO_ID==1103){ /* GISS AOM */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1110 || SCENARIO_ID==1111 || SCENARIO_ID==1112){ /* CCCma-T63 */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1120 || SCENARIO_ID==1121 || SCENARIO_ID==1122 || SCENARIO_ID==1270 || SCENARIO_ID==1271 || SCENARIO_ID==1272){ /* CSIRO 3.0 & 3.5 */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1130 || SCENARIO_ID==1131 || SCENARIO_ID==1132 ||
            SCENARIO_ID==1133 || SCENARIO_ID==1134 || SCENARIO_ID==1135 ||
            SCENARIO_ID==1136 || SCENARIO_ID==1137 || SCENARIO_ID==1138 ||
            SCENARIO_ID==1139 || SCENARIO_ID==1140 || SCENARIO_ID==1141 ||
            SCENARIO_ID==1142 || SCENARIO_ID==1143 || SCENARIO_ID==1144){ /* MRI */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1150 ||SCENARIO_ID==1151 ||SCENARIO_ID==1152 ||SCENARIO_ID==1153 ||SCENARIO_ID==1154 ||
			SCENARIO_ID==1155 ||SCENARIO_ID==1156 ||SCENARIO_ID==1157 ||SCENARIO_ID==1158 ||SCENARIO_ID==1159){ /* MPI ECHAM */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1160 || SCENARIO_ID==1161 || SCENARIO_ID==1162 || SCENARIO_ID==1163 || SCENARIO_ID==1164 || SCENARIO_ID==1165){ /* IAP */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1170 || SCENARIO_ID==1171){ /* HadGEM */
		if(grid->lat>=89.379){
			grid->gcm_row = 0;
		}else if(grid->lat<=-89.379){
			grid->gcm_row = 144;
		}else{
			grid->gcm_row = (long)((89.379-grid->lat)/(180.0/(double)GCM_R));
		}
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1180 || SCENARIO_ID==1181 || SCENARIO_ID==1182){ /* GFDL CM 2.0 */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1190 || SCENARIO_ID==1191 || SCENARIO_ID==1192){ /* CNRM */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1200 || SCENARIO_ID==1201 || SCENARIO_ID==1202 || SCENARIO_ID==1203 || SCENARIO_ID==1204
			 || SCENARIO_ID==1205 || SCENARIO_ID==1206 || SCENARIO_ID==1207 || SCENARIO_ID==1208 || SCENARIO_ID==1209
			 || SCENARIO_ID==1210 || SCENARIO_ID==1211 || SCENARIO_ID==1212 || SCENARIO_ID==1213 || SCENARIO_ID==1214){ /* CCC47 */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1220 || SCENARIO_ID==1221 || SCENARIO_ID==1222 ||
            SCENARIO_ID==1223 || SCENARIO_ID==1224 || SCENARIO_ID==1225 ||
            SCENARIO_ID==1226 || SCENARIO_ID==1227 || SCENARIO_ID==1228 ||
            SCENARIO_ID==1229 || SCENARIO_ID==1230 || SCENARIO_ID==1231 ||
            SCENARIO_ID==1232 || SCENARIO_ID==1233 || SCENARIO_ID==1234 ||
            SCENARIO_ID==1235 || SCENARIO_ID==1236 || SCENARIO_ID==1237 ||
            SCENARIO_ID==1238 || SCENARIO_ID==1239){ /* NCAR CCSM */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1240 || SCENARIO_ID==1241 || SCENARIO_ID==1242){ /* GISS E */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1250 || SCENARIO_ID==1251 || SCENARIO_ID==1252 || SCENARIO_ID==1253){ /* GISS R */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(SCENARIO_ID==1260 || SCENARIO_ID==1261 || SCENARIO_ID==1262 || SCENARIO_ID==1263 || SCENARIO_ID==1264
			 || SCENARIO_ID==1265 || SCENARIO_ID==1266 || SCENARIO_ID==1267){ /* NCAR PCM */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}
    
    /* GEO-MIP: 2013/11/26 by A.Ito */
    if(SCENARIO_ID==3000 || SCENARIO_ID==3001 || SCENARIO_ID==3002){ /* GEO-MIP BNU-ESM */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}if(SCENARIO_ID==3010 || SCENARIO_ID==3011 || SCENARIO_ID==3012){ /* GEO-MIP CSIRO */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}if(SCENARIO_ID==3020 || SCENARIO_ID==3021 || SCENARIO_ID==3022){ /* GEO-MIP GISS */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}if(SCENARIO_ID==3030 || SCENARIO_ID==3031 || SCENARIO_ID==3032 || SCENARIO_ID==3033 ||
            SCENARIO_ID==3034 || SCENARIO_ID==3035){ /* GEO-MIP HadGEM */
		if(grid->lat>=89.379){
			grid->gcm_row = 0;
		}else if(grid->lat<=-89.379){
			grid->gcm_row = 144;
		}else{
			grid->gcm_row = (long)((89.379-grid->lat)/(180.0/(double)GCM_R));
		}
        grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}if(SCENARIO_ID==3040 || SCENARIO_ID==3041 || SCENARIO_ID==3042){ /* GEO-MIP IPSL */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}if(SCENARIO_ID==3050 || SCENARIO_ID==3051 || SCENARIO_ID==3052){ /* GEO-MIP MIROC-ESM */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}if(SCENARIO_ID==3060 || SCENARIO_ID==3061){ /* GEO-MIP MIROC-ESM-CHEM */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}if(SCENARIO_ID==3070 || SCENARIO_ID==3071 || SCENARIO_ID==3072){ /* GEO-MIP CanESM2 */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}if(SCENARIO_ID==3080 || SCENARIO_ID==3081){ /* GEO-MIP MPI-ESM-LR */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}if(SCENARIO_ID==3090 || SCENARIO_ID==3091){ /* GEO-MIP CCSM4 */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}if(SCENARIO_ID==3100 || SCENARIO_ID==3101){ /* GEO-MIP NorESM1-M */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}
    
	/* input geography in the grid *************/	
	fscanf(fp_s[0],"%ld %ld", &country, &region); 
	grid->country = country;	/* country code */
	grid->region = region;		/* region code : note recommended to use Giorgi region code */
	
	/* added: 2010/05/11 */
	switch(grid->country){
		/* OECD countries */
		case 840:	grid->rank_nat = 2;		break;	/* United States */
		case 826:	grid->rank_nat = 2;		break;	/* United Kingdom */
		case 276:	grid->rank_nat = 2;		break;	/* Germany */
		case 250:	grid->rank_nat = 2;		break;	/* France */
		case 380:	grid->rank_nat = 2;		break;	/* Italy */
		case 528:	grid->rank_nat = 2;		break;	/* Netherland */
		case 56:	grid->rank_nat = 2;		break;	/* Belgium */
		case 442:	grid->rank_nat = 2;		break;	/* Luxenburg */
		case 246:	grid->rank_nat = 2;		break;	/* Finland */
		case 752:	grid->rank_nat = 2;		break;	/* Sweden */
		case 40:	grid->rank_nat = 2;		break;	/* Austria */
		case 208:	grid->rank_nat = 2;		break;	/* Denmark */
		case 724:	grid->rank_nat = 2;		break;	/* Spain */
		case 620:	grid->rank_nat = 2;		break;	/* Portugal */
		case 300:	grid->rank_nat = 2;		break;	/* Greece */
		case 372:	grid->rank_nat = 2;		break;	/* Ireland */
		case 200:	grid->rank_nat = 2;		break;	/* Czechoslovakia */
		case 348:	grid->rank_nat = 2;		break;	/* Hungary */
		case 616:	grid->rank_nat = 2;		break;	/* Poland */
		case 392:	grid->rank_nat = 2;		break;	/* Japan */
		case 124:	grid->rank_nat = 2;		break;	/* Canada */
		case 484:	grid->rank_nat = 2;		break;	/* Mexico */
		case 36:	grid->rank_nat = 2;		break;	/* Australia */
		case 554:	grid->rank_nat = 2;		break;	/* New Zealand */
		case 756:	grid->rank_nat = 2;		break;	/* Swiss */
		case 578:	grid->rank_nat = 2;		break;	/* Norway */
		case 352:	grid->rank_nat = 2;		break;	/* Iceland */
		case 792:	grid->rank_nat = 2;		break;	/* Turkey */
		case 410:	grid->rank_nat = 2;		break;	/* South Korea */
			
		default:	grid->rank_nat = 1;	/** developing countries */
	}
	
	/* Olson's vegetation area *****************************/
	/* Olson, J. S., J. A. Watts, et al. (1983). 
	Carbon in live vegetation of major world ecosystems, 
	Oak Ridge National Laboratory.  */
	fscanf(fp_s[1],"%ld", &aaa); 
	grid->veg_olson = aaa;
	/*
	0	water
	1	tropical & subtropical evergreen forest
	2	tropical montane forest
	3	tropical & subtropical dry forest
	4	mid-latitude mixed forest
	5	mid-latitude broad-leaved forest
	6	semiarid wood or low forest
	7	coniferous evergreen forest
	8	southern taiga
	9	main evergreen taiga
	10	main deciduous taiga
	11	northern evergreen taiga
	12	northern deciduous taiga
	13	second growth woods
	14	second growth field
	15	succulent & thorn wood
	16	tropical savanna, woodland
	17	mediterranean-type dry wood
	18	heath & moorland
	19	warm or hot shrub & grassland
	20	tibetan meadow & siberian highland
	21	tundra
	22	wooded tundra
	23	warm or hot wetlands
	24	cool bog & mire
	25	shore & hinterland
	26	cool semi-desert scrub
	27	non-polar desert
	28	non-polar sand desert
	29	paddyland
	30	cool cropland
	31	warm cropland
	32	irrigated
	33	antarctica
	*/
	
	/* sensitivity analysis for biome change by deforestation */
	if(EX_DEFOREST == 1){
		if(aaa>=1 && aaa<=12){
			grid->veg_olson = 19;
		}
	}else if(EX_DEFOREST == 2){
		if(aaa>=1 && aaa<=12){
			grid->veg_olson = 13;
		}
	}else if(EX_DEFOREST == 3){
		if(aaa>=1 && aaa<=12){
			grid->veg_olson = 31;
		}
	}
	
	/* input soil properties *************************************/
	fscanf(fp_s[3],"%lf", &geo_prop);
		geo_prop = (geo_prop>0.0)?geo_prop:0.0;
		grid->topo = geo_prop; /* topography, m ASL */
	fscanf(fp_s[3],"%lf", &geo_prop); 
		geo_prop = (geo_prop>10.0)?geo_prop:10.0;
		/* grid->whc30 = geo_prop; */ /* whc 30cm */
	fscanf(fp_s[3],"%lf", &geo_prop); 
		geo_prop = (geo_prop>10.0)?geo_prop:10.0;
		/* grid->whc = geo_prop; */ /* whc whole */
	fscanf(fp_s[3],"%lf", &geo_prop); 
		geo_prop = (geo_prop>50.0)?geo_prop:50.0;
		grid->sd = geo_prop; /* rooting depth */
	fscanf(fp_s[3],"%lf", &geo_prop); /* soil depth */
	fscanf(fp_s[3],"%lf", &geo_prop); /* slope */
	fscanf(fp_s[3],"%lf", &geo_prop); /* n */
	fscanf(fp_s[3],"%lf", &geo_prop); /* psi_s */
	fscanf(fp_s[3],"%lf", &geo_prop); 
		geo_prop = (geo_prop>0.003)?geo_prop:0.003;
		grid->hyd_cond = geo_prop; /* k_s */
	fscanf(fp_s[3],"%lf",&geo_prop); /* b */
	
	/* input monthly climate *************************************/
	/* long-term average, from NCEP/NCAR reanalysis */
	for(e=0;e<ASTEP;e++){
		fscanf(fp_s[2],"%lf", &tmp_sfc); 
		grid->tmp_sfc_a[e] = tmp_sfc-ZAT;
	}	
	for(e=0;e<ASTEP;e++){
		fscanf(fp_s[2],"%lf", &tmp_2m); 
		grid->tmp_2m_a[e] = tmp_2m-ZAT;
	}	
	for(e=0;e<ASTEP;e++){
		fscanf(fp_s[2],"%lf", &tmp10_soil); 
		grid->tmp10_soil_a[e] = tmp10_soil - ZAT;
	}	
	for(e=0;e<ASTEP;e++){
		fscanf(fp_s[2],"%lf", &tmp200_soil); 
		grid->tmp200_soil_a[e] = tmp200_soil - ZAT;
	}	
	for(e=0;e<ASTEP;e++){
		fscanf(fp_s[2],"%lf", &dswrf_toa); 
	}	
	for(e=0;e<ASTEP;e++){
		fscanf(fp_s[2],"%lf", &dswrf_sfc); 
	}	
	for(e=0;e<ASTEP;e++){
		fscanf(fp_s[2],"%lf", &tcdc_clm); 
		grid->tcdc_clm_a[e] = tcdc_clm/100.0;
	}	
	for(e=0;e<ASTEP;e++){
		fscanf(fp_s[2],"%lf", &prate_sfc); 
		grid->prate_sfc_a[e] = prate_sfc * MDN[e];
	}	
	for(e=0;e<ASTEP;e++){
		fscanf(fp_s[2],"%lf", &spfh_2m); 
		grid->spfh_2m_a[e] = spfh_2m;
	}	
	for(e=0;e<ASTEP;e++){
		fscanf(fp_s[2],"%lf", &soilw10); 
	}	
	for(e=0;e<ASTEP;e++){
		fscanf(fp_s[2],"%lf", &soilw200); 
	}	
	for(e=0;e<ASTEP;e++){
		fscanf(fp_s[2],"%lf", &ugrd_10m); 
		grid->ugrd_10m_a[e] = ugrd_10m;
	}	
	for(e=0;e<ASTEP;e++){
		fscanf(fp_s[2],"%lf", &vgrd_10m); 
		grid->vgrd_10m_a[e] = vgrd_10m;
	}
	
	/* set CRU TS2.1 *******************************************/
	/* set CRU TS3.0 *************** 2010/01/04 (A.Ito) ********/
	/* New, M., D. Lister, et al. (2002). "A high-resolution data set of 
	surface climate over global land areas." Climate Research 21: 1-25. */
	for(h=0;h<DL_HCLIM;h++){
		for(g=0;g<ASTEP;g++){
			grid->hist_cld[h][g] = grid->tcdc_clm_a[g];
			grid->hist_pre[h][g] = grid->prate_sfc_a[g];
			grid->hist_tmp[h][g] = grid->tmp_2m_a[g];
		}
	}
	for(g=0;g<ASTEP;g++){
		grid->hist_cld_b[g] = grid->tcdc_clm_a[g];
		grid->hist_pre_b[g] = grid->prate_sfc_a[g];
		grid->hist_tmp_b[g] = grid->tmp_2m_a[g];
	}
	
	/* substituted precipitation data (UEA/CRU) *********/
	/* probably, no more necessary (2009/08/16) */
	for(e=0;e<ASTEP;e++){
		fscanf(fp_s[4],"%lf", &prate_sfc); 
		grid->prec_sub_a[e] = prate_sfc * MDN[e];
	}
	
	/* river-basin ID */
	fscanf(fp_s[5],"%ld", &(grid->rvbasin)); 
	
	/* bare soil albedo */
	fscanf(fp_s[6],"%lf", &(grid->albedo_soil)); 
	
	/* soil percent organic matter */
	fscanf(fp_s[7],"%lf", &(grid->pcnt_orgmat)); 
	
	/* RUSLE slope factor */
	fscanf(fp_s[8],"%lf", &(grid->fls_slope)); 

	/* RUSLE soil erodibility factor */
	fscanf(fp_s[9],"%lf", &(grid->fk_edodibility)); 

	/* potential vegetation by SAGE (Ramankutty & Foley) */
	fscanf(fp_s[10],"%ld", &(grid->veg_sage)); 
	/*
	0	Ocean
	1	Tropical Evergreen Forest/Woodland
	2	Tropical Deciduous Forest/Woodland
	3	Temperate Broadleaf Evergreen Forest/Woodland
	4	Temperate Needleleaf Evergreen Forest/Woodland
	5	Temperate Deciduous Forest/Woodland
	6	Boreal Evergreen Forest/Woodland
	7	Boreal Deciduous Forest/Woodland
	8	Evergreen/Deciduous Mixed Forest/Woodland
	9	Savanna
	10	Grassland/Steppe
	11	Dense Shrubland
	12	Open Shrubland
	13	Tundra
	14	Desert
	15	Polar Desert/Rock/Ice
	*/
	
	/* Olson croplands replaced by SAGE natural vegetation */
	/* 2010/04/27 by A.Ito */
	if(REPLACE_OLSON_CROP == 1){
		if(grid->veg_olson==29 ||grid->veg_olson==30 ||grid->veg_olson==31 
		   ||grid->veg_olson==32){
			
			switch(grid->veg_sage){
				case 1: grid->veg_olson = 1; break;
				case 2: grid->veg_olson = 3; break;
				case 3: grid->veg_olson = 5; break;
				case 4: grid->veg_olson = 7; break;
				case 5: grid->veg_olson = 4; break;
				case 6: grid->veg_olson = 9; break;
				case 7: grid->veg_olson = 10; break;
				case 8: grid->veg_olson = 4; break;
				case 9: grid->veg_olson = 16; break;
				case 10: grid->veg_olson = 19; break;
				case 11: grid->veg_olson = 19; break;
				case 12: grid->veg_olson = 26; break;
				case 13: grid->veg_olson = 21; break;
				case 14: grid->veg_olson = 27; break;
				case 15: grid->veg_olson = 21; break;
				default: grid->veg_olson = 0; break;
			}
		}
	}
	
	/* experiment for biodiversity: 2011/01/15 (A.Ito) **********************************/
	if(grid->veg_olson >= 1 && grid->veg_olson <= 12){
		/* reduced biodiversity of forests */
		/* grid->veg_olson = 1; */
	}
	
	/* Cropland coverage by Ramankutty & Foley (1999) */
	for(h=0;h<293;h++){
		fscanf(fp_s[11],"%lf", &(grid->fcrop_net[h])); 
	}
	
	/* current linear trend of land-use chnage */
	x = xx = y = yy = xy = 0.0;
	for(h=0;h<10;h++){
		x += (double)h;
		xx += (double)h * (double)h;
		y += grid->fcrop_net[280+h];
		yy += grid->fcrop_net[280+h]*grid->fcrop_net[280+h];
		xy += (double)h*grid->fcrop_net[280+h];
	}
	grid->f_crop_trend = (10.0*xy - x*y)/(10.0*xx - x*x);
	
	/* future land-use change scenario by IMAGE2 (Wang et al. 2006) */
	/*
	Wang, A., Price, D.T. and Arora, V., 2006. Estimating changes in global 
	vegetation cover (1850-2100) for use in climate models. 
	Global Biogeochemical Cycles 20, 10.1029/2005GB002514.
	*/
    
    if(LANDUSE == 19 || LANDUSE == 20 ||LANDUSE == 21 ||LANDUSE == 22 ||LANDUSE == 23){
        /* CD-LINK 2016/11/17 */
        for(h=0;h<111;h++){
            fscanf(fp_s[23],"%lf", &(grid->fcrop3_future[h]));
            grid->fcrop4_future[h] = 0.0;
        }
        for(h=0;h<111;h++){
            fscanf(fp_s[45],"%lf", &(grid->fgrass3_future[h]));
        }
        for(h=0;h<111;h++){
            fscanf(fp_s[45],"%lf", &(grid->fgrass4_future[h])); 
        }
    }else{
        if(SCENARIO_ID>=2100 && SCENARIO_ID<=2999){
            /* ICARUS 2016/08/12 */
            for(h=0;h<111;h++){
                fscanf(fp_s[23],"%lf", &(grid->fcrop3_future[h]));
                grid->fcrop4_future[h] = 0.0;
            }
            for(h=0;h<111;h++){
                fscanf(fp_s[45],"%lf", &(grid->fgrass3_future[h]));
            }
            for(h=0;h<111;h++){
                fscanf(fp_s[45],"%lf", &(grid->fgrass4_future[h])); 
            }
        }else{
            for(h=0;h<111;h++){
                fscanf(fp_s[23],"%lf", &(grid->fcrop3_future[h])); 
            }
            for(h=0;h<111;h++){
                fscanf(fp_s[23],"%lf", &(grid->fcrop4_future[h])); 
            }
            for(h=0;h<111;h++){
                fscanf(fp_s[45],"%lf", &(grid->fgrass3_future[h])); 
            }
            for(h=0;h<111;h++){
                fscanf(fp_s[45],"%lf", &(grid->fgrass4_future[h])); 
            }
        }
    }
	
	/* soil texture */
	fscanf(fp_s[12],"%ld", &(grid->soiltexture)); 
	/* 0: ocean */
	/* 1: coarse */
	/* 2: medium/coarse */
	/* 3: medium */
	/* 4: fine/medium */
	/* 5: fine */
	/* 6: ice */
	/* 7: organic */
	
	/* soil bulk density */
	fscanf(fp_s[13],"%lf", &bulkdens); 
	
	/* soil field capacity */
	fscanf(fp_s[14],"%lf", &field_cap); 
	
	fscanf(fp_s[15],"%lf", &(grid->fc_150)); 
	if(grid->fc_150 <= 0.1){
		grid->fc_150 = olveg_fc150[grid->veg_olson];
	}
	fscanf(fp_s[16],"%lf", &(grid->fc_30)); 
	if(grid->fc_30 <= 0.1){
		grid->fc_30 = olveg_fc30[grid->veg_olson];
	}
	fscanf(fp_s[17],"%lf", &(grid->bulkdens)); 
	if(grid->bulkdens <= 0.1){
		grid->bulkdens = olveg_bd[grid->veg_olson];
	}
	fscanf(fp_s[18],"%lf", &(grid->pc_clay)); 
	if(grid->pc_clay <= 0.1){
		grid->pc_clay = olveg_pclay[grid->veg_olson];
	}
	fscanf(fp_s[19],"%lf", &(grid->pc_sand)); 
	if(grid->pc_sand<= 0.1){
		grid->pc_sand = olveg_psand[grid->veg_olson];
	}
	
	fscanf(fp_s[20],"%lf", &(grid->soil_ph)); 
	if(grid->soil_ph <= 0.5 || grid->soil_ph > 14.0){
		grid->soil_ph = olveg_ph[grid->veg_olson];
	}
	
	switch(grid->soiltexture){
		case 0: crit_tension = 33.0;	break;
		case 1: crit_tension = 10.0;	break;
		case 2: crit_tension = 10.0;	break;
		case 3: crit_tension = 33.0;	break;
		case 4: crit_tension = 33.0;	break;
		case 5: crit_tension = 33.0;	break;
		case 6: crit_tension = 33.0;	break;
		case 7: crit_tension = 33.0;	break;
		default: crit_tension = 33.0;
	}
	
	/* soil physical properties by Saxton (1986) ************/
	/* Saxton, K. E., et al. (1986), Estimating generalized soil-water characteristics 
	from texture, Soil Science Society of America Journal, 50, 1031-1036. */
	grid->a_sw = exp(-4.396 -0.0715*grid->pc_clay -4.488*0.0001*grid->pc_clay*grid->pc_clay 
					- 4.285*0.00001*grid->pc_sand*grid->pc_sand*grid->pc_clay)*100.0;
	grid->b_sw = -3.14 -0.00222 * grid->pc_clay * grid->pc_clay
					-3.484 * 0.00001 * grid->pc_sand * grid->pc_sand * grid->pc_clay;
	
	/* field capacity */
	grid->field_cap1 = pow(crit_tension/grid->a_sw, 1.0/grid->b_sw);
	grid->field_cap1 *= 300.0;
	/* pore space capacity */
	grid->pore_cap1 = 0.332 - 7.251 * 0.0001 * grid->pc_sand + 0.1276 * log10(grid->pc_clay);
	grid->pore_cap1 *= 300.0;
	
	grid->field_cap2 = grid->fc_150 - grid->field_cap1;
	if(grid->field_cap2 < 100.0){
		grid->field_cap2 = 100.0;
	}
	
	/* wetland fraction: data by Global Lakes and Wetlands Database by WWF *****/
	if(ALT_FWETLAND == 1){
		/* Alternative data (NASA/GISS): 2011/03/30 by A.Ito */
		fscanf(fp_s[21],"%lf", &wetland);
		grid->f_wetland = wetland;
        /* corrected: A. Ito (with Hamada-san's comment) 2012/01/30 */
        lake = 0.0;
	}else{
		/* default */
		/* revised wetland data: by A.Ito (2009/07/14) */
		fscanf(fp_s[21],"%lf %lf %lf %lf %lf", &lat, &lon, &total, &lake, &wetland); 
		grid->f_wetland = wetland/grid->area;
        
        /* WSL-Perogon welnad map: 2014/2/4 by A.Ito */
        fscanf(fp_s[83],"%ld %lf", &aaa, &wetland); 
        if(ALT_FWETLAND==2 && aaa==1){
            grid->f_wetland = wetland;
        }
        
		if(grid->f_wetland > 1.0){
			grid->f_wetland = 1.0;
		}

	}
	grid->f_lake = lake/grid->area;
	if(grid->f_lake > 1.0){
		grid->f_lake = 1.0;
	}
    
    /* base wetland extent: 2012/10/26 by A.Ito */
    grid->f_wetland0 = grid->f_wetland;
	
	/* paddy fraction *****************/
	if(ALT_FWETLAND == 1){
		/* Alternative data (IIS-UT + SAGE): 2011/03/30 by A.Ito */
		fscanf(fp_s[22],"%lf", &paddy); 
		if(paddy>0.0){
			grid->f_paddy = paddy;
			grid->f_paddy_b = paddy;
		}else{
			grid->f_paddy = 0.0;
			grid->f_paddy_b = 0.0;
		}
	}else{
		/* default: data by U.Wisconsin SAGE (Leff et al.) */
		fscanf(fp_s[22],"%lf", &paddy); 
		if(paddy > 0.0){
			grid->f_paddy = paddy;
			grid->f_paddy_b = paddy;
		}else{
			grid->f_paddy = 0.0;
			grid->f_paddy_b = 0.0;
		}
	}
	
	grid->f_upland = 1.0 - grid->f_wetland - grid->f_lake - grid->f_paddy;
	if(grid->f_upland < 0.0){
		grid->f_upland = 0.0;
	}
	
	/* soil total nitrogen */
	fscanf(fp_s[24],"%lf", &grid->total_n_1m); 
	if(grid->total_n_1m < 0.0){
		grid->total_n_1m = 0.0;
	}
	
    if(ISIMIP_RUN==4 && (SCENARIO_ID>=5010 && SCENARIO_ID<=5100)){
        /* ISI-MIP2b: 2016/12/24 by A.Ito */
        fread(is2bdat,sizeof(float),DL_AGHG, fp_s[25]);
        for(e=0;e<DL_NINPUT;e++){
            grid->nmip_ndep_nh4[e] = is2bdat[e];
        }
        
        fread(is2bdat,sizeof(float),DL_AGHG, fp_s[25]);
        for(e=0;e<DL_NINPUT;e++){
            grid->nmip_ndep_noy[e] = is2bdat[e];
        }
        
    }else{
        /* N deposition: Gallway & Dentener */
        /* URL  daac.ornl.gov/CLIMATE/guides/global_N_deposition_maps.html */
        fscanf(fp_s[25],"%lf", &lat);
        fscanf(fp_s[25],"%lf", &lon); 
        fscanf(fp_s[25],"%lf", &grid->ndepo[0]); 
        fscanf(fp_s[25],"%lf", &grid->ndepo[1]); 
        fscanf(fp_s[25],"%lf", &grid->ndepo[2]);
    }
		
	/* EOS-WEBSTER Land-use change data *********************/
	/*
	 Hurtt, G. C., et al. (2006), The underpinnings of land-use history: three centuries of 
	 global gridded land-use transitions, wood-harvest activity, and resulting secondary lands, 
	 Global Change Biology, 12, 1-22.
	 */
    /* revised by A.Ito (2013/12/20) */
	if(LANDUSE==6 || LANDUSE==8 || LANDUSE==9 || LANDUSE==10 || LANDUSE==11 || LANDUSE==12 ||
         LANDUSE==13 || LANDUSE==17|| LANDUSE==18|| LANDUSE == 19 || LANDUSE == 20 || LANDUSE == 21
        || LANDUSE == 22 || LANDUSE == 23){
        
		for(h=0;h<DL_LUC;h++){
			/* fractional cover */
            if(h<(BGY_GCM - FDY_LUC -1)){  /* 1500-2004 */
                fscanf(fp_s[26],"%lf", &grid->fcrop_unh_hmnzed[h]);
                fscanf(fp_s[27],"%lf", &grid->fpast_unh_hmnzed[h]);
                fscanf(fp_s[28],"%lf", &grid->fprim_unh_hmnzed[h]);
                fscanf(fp_s[29],"%lf", &grid->fsecd_unh_hmnzed[h]);
                fscanf(fp_s[30],"%lf", &grid->ssma_unh_hmnzed[h]);
                fscanf(fp_s[31],"%lf", &grid->ssmb_unh_hmnzed[h]);
                fscanf(fp_s[32],"%lf", &grid->t_cp_unh_hmnzed[h]);
                fscanf(fp_s[33],"%lf", &grid->t_cs_unh_hmnzed[h]);
                fscanf(fp_s[34],"%lf", &grid->t_pc_unh_hmnzed[h]);
                fscanf(fp_s[35],"%lf", &grid->t_ps_unh_hmnzed[h]);
                fscanf(fp_s[36],"%lf", &grid->t_sc_unh_hmnzed[h]);
                fscanf(fp_s[37],"%lf", &grid->t_sp_unh_hmnzed[h]);
                fscanf(fp_s[38],"%lf", &grid->t_ss1_unh_hmnzed[h]);
                fscanf(fp_s[39],"%lf", &grid->t_ss2_unh_hmnzed[h]);
                fscanf(fp_s[40],"%lf", &grid->t_ss3_unh_hmnzed[h]);
                fscanf(fp_s[41],"%lf", &grid->t_vc_unh_hmnzed[h]);
                fscanf(fp_s[42],"%lf", &grid->t_vp_unh_hmnzed[h]);
                fscanf(fp_s[43],"%lf", &grid->t_vs1_unh_hmnzed[h]);
                fscanf(fp_s[44],"%lf", &grid->t_vs2_unh_hmnzed[h]);
            }else{  /* 2005-2100 */
                fscanf(fp_s[59],"%lf", &grid->fcrop_unh_hmnzed[h]);
                fscanf(fp_s[60],"%lf", &grid->fpast_unh_hmnzed[h]);
                fscanf(fp_s[61],"%lf", &grid->fprim_unh_hmnzed[h]);
                fscanf(fp_s[62],"%lf", &grid->fsecd_unh_hmnzed[h]);
                fscanf(fp_s[63],"%lf", &grid->ssma_unh_hmnzed[h]);
                fscanf(fp_s[64],"%lf", &grid->ssmb_unh_hmnzed[h]);
                fscanf(fp_s[65],"%lf", &grid->t_cp_unh_hmnzed[h]);
                fscanf(fp_s[66],"%lf", &grid->t_cs_unh_hmnzed[h]);
                fscanf(fp_s[67],"%lf", &grid->t_pc_unh_hmnzed[h]);
                fscanf(fp_s[68],"%lf", &grid->t_ps_unh_hmnzed[h]);
                fscanf(fp_s[69],"%lf", &grid->t_sc_unh_hmnzed[h]);
                fscanf(fp_s[70],"%lf", &grid->t_sp_unh_hmnzed[h]);
                fscanf(fp_s[71],"%lf", &grid->t_ss1_unh_hmnzed[h]);
                fscanf(fp_s[72],"%lf", &grid->t_ss2_unh_hmnzed[h]);
                fscanf(fp_s[73],"%lf", &grid->t_ss3_unh_hmnzed[h]);
                fscanf(fp_s[74],"%lf", &grid->t_vc_unh_hmnzed[h]);
                fscanf(fp_s[75],"%lf", &grid->t_vp_unh_hmnzed[h]);
                fscanf(fp_s[76],"%lf", &grid->t_vs1_unh_hmnzed[h]);
                fscanf(fp_s[77],"%lf", &grid->t_vs2_unh_hmnzed[h]);
            }
		}
        
        if(LANDUSE==9 || LANDUSE==10 || LANDUSE==11 || LANDUSE==12 || LANDUSE==13 ||
            LANDUSE==17 || LANDUSE==18|| LANDUSE == 19 || LANDUSE == 20 ||
            LANDUSE == 21 || LANDUSE == 22 || LANDUSE == 23){
            /* skip Hist - 2005 data */
            fscanf(fp_s[26],"%lf", &ddummy);
            fscanf(fp_s[27],"%lf", &ddummy);
            fscanf(fp_s[28],"%lf", &ddummy);
            fscanf(fp_s[29],"%lf", &ddummy);
            fscanf(fp_s[30],"%lf", &ddummy);
            fscanf(fp_s[31],"%lf", &ddummy);
            fscanf(fp_s[32],"%lf", &ddummy);
            fscanf(fp_s[33],"%lf", &ddummy);
            fscanf(fp_s[34],"%lf", &ddummy);
            fscanf(fp_s[35],"%lf", &ddummy);
            fscanf(fp_s[36],"%lf", &ddummy);
            fscanf(fp_s[37],"%lf", &ddummy);
            fscanf(fp_s[38],"%lf", &ddummy);
            fscanf(fp_s[39],"%lf", &ddummy);
            fscanf(fp_s[40],"%lf", &ddummy);
            fscanf(fp_s[41],"%lf", &ddummy);
            fscanf(fp_s[42],"%lf", &ddummy);
            fscanf(fp_s[43],"%lf", &ddummy);
            fscanf(fp_s[44],"%lf", &ddummy);
        }
    
	}else if(LANDUSE == 14 || LANDUSE == 15 || LANDUSE == 16){
    
        for(h=2001;h<=2004;h++){
            fscanf(fp_s[59],"%lf", &ddummy);
            fscanf(fp_s[60],"%lf", &ddummy);
            fscanf(fp_s[61],"%lf", &ddummy);
            fscanf(fp_s[62],"%lf", &ddummy);
            fscanf(fp_s[63],"%lf", &ddummy);
            fscanf(fp_s[64],"%lf", &ddummy);
            fscanf(fp_s[65],"%lf", &ddummy);
            fscanf(fp_s[66],"%lf", &ddummy);
            fscanf(fp_s[67],"%lf", &ddummy);
            fscanf(fp_s[68],"%lf", &ddummy);
            fscanf(fp_s[69],"%lf", &ddummy);
            fscanf(fp_s[70],"%lf", &ddummy);
            fscanf(fp_s[71],"%lf", &ddummy);
            fscanf(fp_s[72],"%lf", &ddummy);
            fscanf(fp_s[73],"%lf", &ddummy);
            fscanf(fp_s[74],"%lf", &ddummy);
            fscanf(fp_s[75],"%lf", &ddummy);
            fscanf(fp_s[76],"%lf", &ddummy);
            fscanf(fp_s[77],"%lf", &ddummy);
        }
    
        for(h=FDY_LUC;h<(FDY_LUC+DL_LUC);h++){
			/* fractional cover */
            if(h < (BGY_GCM-1)){  /* 1500-2005 */
                fscanf(fp_s[26],"%lf", &grid->fcrop_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[27],"%lf", &grid->fpast_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[28],"%lf", &grid->fprim_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[29],"%lf", &grid->fsecd_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[30],"%lf", &grid->ssma_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[31],"%lf", &grid->ssmb_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[32],"%lf", &grid->t_cp_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[33],"%lf", &grid->t_cs_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[34],"%lf", &grid->t_pc_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[35],"%lf", &grid->t_ps_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[36],"%lf", &grid->t_sc_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[37],"%lf", &grid->t_sp_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[38],"%lf", &grid->t_ss1_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[39],"%lf", &grid->t_ss2_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[40],"%lf", &grid->t_ss3_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[41],"%lf", &grid->t_vc_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[42],"%lf", &grid->t_vp_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[43],"%lf", &grid->t_vs1_unh_hmnzed[h-FDY_LUC]);
                fscanf(fp_s[44],"%lf", &grid->t_vs2_unh_hmnzed[h-FDY_LUC]);
            }
            
            if(h >= (BGY_GCM-1)){  /* 2005-2100 */
                fscanf(fp_s[59],"%lf", &grid->fcrop_unh_hmnzed[h-FDY_LUC]);
                if(grid->fcrop_unh_hmnzed[h-FDY_LUC] < 0.0){
                    grid->fcrop_unh_hmnzed[h-FDY_LUC] = 0.0;
                }
                fscanf(fp_s[60],"%lf", &grid->fpast_unh_hmnzed[h-FDY_LUC]);
                if(grid->fpast_unh_hmnzed[h-FDY_LUC] < 0.0){
                    grid->fpast_unh_hmnzed[h-FDY_LUC] = 0.0;
                }
                fscanf(fp_s[61],"%lf", &ddummy); /* other */
                fscanf(fp_s[62],"%lf", &grid->fsecd_unh_hmnzed[h-FDY_LUC]);
                if(grid->fsecd_unh_hmnzed[h-FDY_LUC] < 0.0){
                    grid->fsecd_unh_hmnzed[h-FDY_LUC] = 0.0;
                }
                fscanf(fp_s[63],"%lf", &ddummy); /* grass */
                fscanf(fp_s[64],"%lf", &grid->fprim_unh_hmnzed[h-FDY_LUC]);
                if(grid->fprim_unh_hmnzed[h-FDY_LUC] < 0.0){
                    grid->fprim_unh_hmnzed[h-FDY_LUC] = 0.0;
                }
                fscanf(fp_s[65],"%lf", &grid->t_cp_unh_hmnzed[h-FDY_LUC]);
                if(grid->t_cp_unh_hmnzed[h-FDY_LUC] < 0.0){
                    grid->t_cp_unh_hmnzed[h-FDY_LUC] = 0.0;
                }
                fscanf(fp_s[66],"%lf", &grid->t_cs_unh_hmnzed[h-FDY_LUC]);
                if(grid->t_cs_unh_hmnzed[h-FDY_LUC] < 0.0){
                    grid->t_cs_unh_hmnzed[h-FDY_LUC] = 0.0;
                }
                fscanf(fp_s[67],"%lf", &ddummy);
                fscanf(fp_s[68],"%lf", &grid->t_pc_unh_hmnzed[h-FDY_LUC]);
                if(grid->t_pc_unh_hmnzed[h-FDY_LUC] < 0.0){
                    grid->t_pc_unh_hmnzed[h-FDY_LUC] = 0.0;
                }
                fscanf(fp_s[69],"%lf", &grid->t_ps_unh_hmnzed[h-FDY_LUC]);
                if(grid->t_ps_unh_hmnzed[h-FDY_LUC] < 0.0){
                    grid->t_ps_unh_hmnzed[h-FDY_LUC] = 0.0;
                }
                fscanf(fp_s[70],"%lf", &ddummy);
                fscanf(fp_s[71],"%lf", &grid->t_sc_unh_hmnzed[h-FDY_LUC]);
                if(grid->t_sc_unh_hmnzed[h-FDY_LUC] < 0.0){
                    grid->t_sc_unh_hmnzed[h-FDY_LUC] = 0.0;
                }
                fscanf(fp_s[72],"%lf", &grid->t_sp_unh_hmnzed[h-FDY_LUC]);
                if(grid->t_sp_unh_hmnzed[h-FDY_LUC] < 0.0){
                    grid->t_sp_unh_hmnzed[h-FDY_LUC] = 0.0;
                }
                fscanf(fp_s[73],"%lf", &ddummy);
                fscanf(fp_s[74],"%lf", &grid->t_vc_unh_hmnzed[h-FDY_LUC]);
                if(grid->t_vc_unh_hmnzed[h-FDY_LUC] < 0.0){
                    grid->t_vc_unh_hmnzed[h-FDY_LUC] = 0.0;
                }
                fscanf(fp_s[75],"%lf", &grid->t_vp_unh_hmnzed[h-FDY_LUC]);
                if(grid->t_vp_unh_hmnzed[h-FDY_LUC] < 0.0){
                    grid->t_vp_unh_hmnzed[h-FDY_LUC] = 0.0;
                }
                fscanf(fp_s[76],"%lf", &grid->t_vs1_unh_hmnzed[h-FDY_LUC]);
                if(grid->t_vs1_unh_hmnzed[h-FDY_LUC] < 0.0){
                    grid->t_vs1_unh_hmnzed[h-FDY_LUC] = 0.0;
                }
                fscanf(fp_s[77],"%lf", &ddummy);
                
                grid->t_vs2_unh_hmnzed[h-FDY_LUC] = 0.0;
                grid->ssma_unh_hmnzed[h-FDY_LUC] = 0.0;
                grid->ssmb_unh_hmnzed[h-FDY_LUC] = 0.0;
                grid->t_ss1_unh_hmnzed[h-FDY_LUC] = 0.0;
                grid->t_ss2_unh_hmnzed[h-FDY_LUC] = 0.0;
                grid->t_ss3_unh_hmnzed[h-FDY_LUC] = 0.0;
            }
        }
        
        /* skip Hist - 2005 data */
        fscanf(fp_s[26],"%lf", &ddummy);
        fscanf(fp_s[27],"%lf", &ddummy);
        fscanf(fp_s[28],"%lf", &ddummy);
        fscanf(fp_s[29],"%lf", &ddummy);
        fscanf(fp_s[30],"%lf", &ddummy);
        fscanf(fp_s[31],"%lf", &ddummy);
        fscanf(fp_s[32],"%lf", &ddummy);
        fscanf(fp_s[33],"%lf", &ddummy);
        fscanf(fp_s[34],"%lf", &ddummy);
        fscanf(fp_s[35],"%lf", &ddummy);
        fscanf(fp_s[36],"%lf", &ddummy);
        fscanf(fp_s[37],"%lf", &ddummy);
        fscanf(fp_s[38],"%lf", &ddummy);
        fscanf(fp_s[39],"%lf", &ddummy);
        fscanf(fp_s[40],"%lf", &ddummy);
        fscanf(fp_s[41],"%lf", &ddummy);
        fscanf(fp_s[42],"%lf", &ddummy);
        fscanf(fp_s[43],"%lf", &ddummy);
        fscanf(fp_s[44],"%lf", &ddummy);
    }else{
		for(h=0;h<DL_LUC;h++){
			grid->fcrop_unh_hmnzed[h] = 0.0;
			grid->fpast_unh_hmnzed[h] = 0.0;
			grid->fprim_unh_hmnzed[h] = 0.0;
			grid->fsecd_unh_hmnzed[h] = 0.0;
			grid->ssma_unh_hmnzed[h] = 0.0;
			grid->ssmb_unh_hmnzed[h] = 0.0;
            grid->t_cp_unh_hmnzed[h] = 0.0;
            grid->t_cs_unh_hmnzed[h] = 0.0;
            grid->t_pc_unh_hmnzed[h] = 0.0;
            grid->t_ps_unh_hmnzed[h] = 0.0;
            grid->t_sc_unh_hmnzed[h] = 0.0;
            grid->t_sp_unh_hmnzed[h] = 0.0;
            grid->t_ss1_unh_hmnzed[h] = 0.0;
            grid->t_ss2_unh_hmnzed[h] = 0.0;
            grid->t_ss3_unh_hmnzed[h] = 0.0;
            grid->t_vc_unh_hmnzed[h] = 0.0;
            grid->t_vp_unh_hmnzed[h] = 0.0;
            grid->t_vs1_unh_hmnzed[h] = 0.0;
            grid->t_vs2_unh_hmnzed[h] = 0.0;
		}
	}
	
	/* wood harvest based on RCP-harmonized data: LUHa.v1 */
	/* added by A.Ito (2010/10/15) */
	/* revised by A.Ito (2013/12/20, 24) */
    
    if(LANDUSE == 6 || LANDUSE == 8 || LANDUSE == 9 || LANDUSE == 10 ||
                LANDUSE==11 || LANDUSE==12 || LANDUSE==13 || LANDUSE==17){
        
        for(h=0;h<DL_LUC;h++){
            if(h<(BGY_GCM - FDY_LUC - 1)){
                fscanf(fp_s[53],"%lf", &grid->hvst_p1[h]);
                fscanf(fp_s[54],"%lf", &grid->hvst_p2[h]);
                fscanf(fp_s[55],"%lf", &grid->hvst_s1[h]);
                fscanf(fp_s[56],"%lf", &grid->hvst_s2[h]);
                fscanf(fp_s[57],"%lf", &grid->hvst_s3[h]);
            }else{
                fscanf(fp_s[78],"%lf", &grid->hvst_p1[h]);
                fscanf(fp_s[79],"%lf", &grid->hvst_p2[h]);
                fscanf(fp_s[80],"%lf", &grid->hvst_s1[h]);
                fscanf(fp_s[81],"%lf", &grid->hvst_s2[h]);
                fscanf(fp_s[82],"%lf", &grid->hvst_s3[h]);
            }
        }
        
        if(LANDUSE==9 || LANDUSE==10 || LANDUSE==11 || LANDUSE==12 || LANDUSE==13 || LANDUSE==17){
            fscanf(fp_s[53],"%lf", &ddummy);
            fscanf(fp_s[54],"%lf", &ddummy);
            fscanf(fp_s[55],"%lf", &ddummy);
            fscanf(fp_s[56],"%lf", &ddummy);
            fscanf(fp_s[57],"%lf", &ddummy);
        }
        
    }else if(LANDUSE == 14){
        
        for(h=2001;h<=2004;h++){
            fscanf(fp_s[78],"%lf", &ddummy);
        }
        
        for(h=FDY_LUC;h<(FDY_LUC+DL_LUC);h++){

            if(h < (BGY_GCM-1)){ /* 1500-2004 */
                fscanf(fp_s[53],"%lf", &grid->hvst_p1[h-FDY_LUC]);
                fscanf(fp_s[54],"%lf", &grid->hvst_p2[h-FDY_LUC]);
                fscanf(fp_s[55],"%lf", &grid->hvst_s1[h-FDY_LUC]);
                fscanf(fp_s[56],"%lf", &grid->hvst_s2[h-FDY_LUC]);
                fscanf(fp_s[57],"%lf", &grid->hvst_s3[h-FDY_LUC]);
            }
        
            if(h >= (BGY_GCM-1)){ /* 2005-2100 */
                fscanf(fp_s[78],"%lf", &grid->hvst_p1[h-FDY_LUC]);
                if(grid->hvst_p1[h-FDY_LUC] < 0.0){
                    grid->hvst_p1[h-FDY_LUC] = 0.0;
                }
                
                /* m3/grid/yr => MgC/grid/yr */
                grid->hvst_p1[h-FDY_LUC] = grid->hvst_p1[h-FDY_LUC] / 2.0 / 2.0;
                
                grid->hvst_p2[h-FDY_LUC] = 0.0;
                grid->hvst_s1[h-FDY_LUC] = 0.0;
                grid->hvst_s2[h-FDY_LUC] = 0.0;
                grid->hvst_s3[h-FDY_LUC] = 0.0;
            }
        }
        
        fscanf(fp_s[53],"%lf", &ddummy);
        fscanf(fp_s[54],"%lf", &ddummy);
        fscanf(fp_s[55],"%lf", &ddummy);
        fscanf(fp_s[56],"%lf", &ddummy);
        fscanf(fp_s[57],"%lf", &ddummy);

    }else{
        for(h=0;h<DL_LUC;h++){
            grid->hvst_p1[h] = 0.0;
            grid->hvst_p2[h] = 0.0;
            grid->hvst_s1[h] = 0.0;
            grid->hvst_s2[h] = 0.0;
            grid->hvst_s3[h] = 0.0;
        }
    }
	
	/* crop and pasture fractions: 1700-2007 */
	/* Ramankutty & Kimball: added 2010/01/07 (A.Ito) */
	for(h=0;h<308;h++){
		/* crop */
		fscanf(fp_s[50],"%lf", &grid->fcrop_rk[h]);
		if(grid->fcrop_rk[h] < 0.0){
			grid->fcrop_rk[h] = 0.0;
		}
		if(grid->fcrop_rk[h] > 1.0){
			grid->fcrop_rk[h] = 1.0;
		}
		/* pasture */
		fscanf(fp_s[51],"%lf", &grid->fpast_rk[h]);
		if(grid->fpast_rk[h] < 0.0){
			grid->fpast_rk[h] = 0.0;
		}
		if(grid->fpast_rk[h] > 1.0){
			grid->fpast_rk[h] = 1.0;
		}
	}
	
	/* base cropland and pasture area at 2000 */
	grid->f_crop_base = grid->fcrop_unh_hmnzed[2000 - 1700];
	grid->f_pasture_base = grid->fpast_unh_hmnzed[2000 - 1700];
	
	/* crop type ******************************************/	
	fscanf(fp_s[46],"%lf", &grid->fcrop);
	fscanf(fp_s[46],"%lf", &grid->frice);
	fscanf(fp_s[46],"%lf", &grid->fwheat);
	fscanf(fp_s[46],"%lf", &grid->fmaize);
	
	grid->type_crop = 1;
	/* 0: no crop */
	/* 1: C3 crops (non-rice) */
	/* 2: paddy */
	/* 3: C4 crops, e.g. maize */
	if(grid->frice > grid->fwheat && grid->frice > grid->fmaize){
		grid->type_crop = 2;
	}
	if(grid->fmaize > grid->fwheat && grid->fmaize > grid->frice){
		grid->type_crop = 3;
	}
    
    /* force change crop types: 2015/04/24 by A.Ito */
    if(EX_CROP == 1){
        grid->type_crop = 1;  /* C3:wheat */
    }
    if(EX_CROP == 2){
        grid->type_crop = 2;  /* rice */
    }
    if(EX_CROP == 3){
        grid->type_crop = 3;  /* C4:maize */
    }
    
    if(BIOFUEL_RUN >= 1){
         grid->type_crop = 4;  /* biofuel: added 2015/08/21 by A.Ito */
    }
	
	/* diffuse radiation estimation using SRB data ************/
	/* intercept */
	fscanf(fp_s[47],"%lf", &grid->srb_dif_aa);
	/* slope */
	fscanf(fp_s[47],"%lf", &grid->srb_dif_bb);
	/* R2 */
	fscanf(fp_s[47],"%lf", &grid->srb_dif_rr);
	fscanf(fp_s[47],"%lf", &grid->srb_dif_min_x);
	fscanf(fp_s[47],"%lf", &grid->srb_dif_max_x);
	fscanf(fp_s[47],"%lf", &grid->srb_dif_min_y);
	fscanf(fp_s[47],"%lf", &grid->srb_dif_max_y);
	
	/* inundation ********************************************/
	/* SSMI */
    fscanf(fp_s[48],"%lf", &lat);
	fscanf(fp_s[48],"%lf", &lon);
	grid->inundation_ssmi_av = 0.0;
	grid->inundation_ssmi_max = 0.0;
	for(h=0;h<ASTEP;h++){
		fscanf(fp_s[48],"%ld", &aaa);
		grid->inundation_ssmi[h] = (double)aaa/8.0;
		
		grid->inundation_ssmi_av += grid->inundation_ssmi[h] * MDN[h]/365.0;
		
		if(grid->inundation_ssmi[h] > grid->inundation_ssmi_max){
			grid->inundation_ssmi_max = grid->inundation_ssmi[h];
		}
	}
    
    /* NASA-JPL for GCP-CH4: 2014/05/28 (rev 2014/08/17) by A.Ito */
	fscanf(fp_s[84],"%lf", &lat);
	fscanf(fp_s[84],"%lf", &lon);
    if(ALT_INUND == 6){
        /* 2000/01-2012/12 */
        for(h=0;h<ASTEP;h++){
            grid->inundation_gcp_av[h] = 0.0;
        }
        for(g=0;g<13;g++){
            for(h=0;h<ASTEP;h++){
                fscanf(fp_s[84],"%lf", &grid->inundation_gcp_ts[g][h]);
                
                if(grid->inundation_gcp_ts[g][h] < 0.0){
                    grid->inundation_gcp_ts[g][h] = 0.0;
                }
                
                grid->inundation_gcp_av[h] += grid->inundation_gcp_ts[g][h] / 13.0;
            }
        }
    }else{
        /* 1999/07-2013/03 */
        /* average */
        for(h=0;h<ASTEP;h++){
            fscanf(fp_s[84],"%lf", &grid->inundation_gcp_av[h]);
        }
        for(g=0;g<15;g++){
            for(h=0;h<ASTEP;h++){
                fscanf(fp_s[84],"%lf", &grid->inundation_gcp_ts[g][h]);
                
                if(grid->inundation_gcp_ts[g][h] < 0.0){
                    grid->inundation_gcp_ts[g][h] = 0.0;
                }
            }
        }
    }
	
	/* permafrost type **************************************/
	/* from National Snow and Ice Data Center 
	 http://nsidc.org/data/docs/fgdc/ggd318_map_circumarctic/index.html
	*/
	fscanf(fp_s[49],"%ld", &grid->type_permaforst);
	/*
	 0: No information
	 1: Continuous permafrost extent with high ground ice content and thick overburden
	 2: Discontinuous permafrost extent with high ground ice content and thick overburden
	 3: Sporadic permafrost extent with high ground ice content and thick overburden
	 4: Isolated patches of permafrost extent with high ground ice content and thick overburden
	 5: Continuous permafrost extent with medium ground ice content and thick overburden
	 6: Discontinuous permafrost extent with medium ground ice content and thick overburden
	 7: Sporadic permafrost extent with medium ground ice content and thick overburden
	 8: Isolated patches of permafrost extent with medium ground ice content and thick overburden
	 9: Continuous permafrost extent with low ground ice content and thick overburden
	 10: Discontinuous permafrost extent with low ground ice content and thick overburden
	 11: Sporadic permafrost extent with low ground ice content and thick overburden
	 12: Isolated patches of permafrost extent with low ground ice content and thick overburden
	 13: Continuous permafrost extent with high ground ice content and thin overburden and exposed bedrock
	 14: Discontinuous permafrost extent with high ground ice content and thin overburden and exposed bedrock
	 15: Sporadic permafrost extent with high ground ice content and thin overburden and exposed bedrock
	 16: Isolated patches of permafrost extent with high ground ice content and thin overburden and exposed bedrock
	 17: Continuous permafrost extent with low ground ice content and thin overburden and exposed bedrock
	 18: Discontinuous permafrost extent with low ground ice content and thin overburden and exposed bedrock
	 19: Sporadic permafrost extent with low ground ice content and thin overburden and exposed bedrock
	 20: Isolated patches of permafrost extent with low ground ice content and thin overburden and exposed bedrock
	 21: Glaciers
	 22: Relict permafrost
	 23: Inland lakes
	 24: Ocean/inland seas
	 25: Land
	*/
	
	/* CHASER 2001 monthly, by A.Ito (2010/05/21) ******************************/
	grid->chaser_row = grid->row/(360.0/64.0);
	if(grid->col >= 360){
		grid->chaser_col = grid->col/(720.0/128.0) - 64;
	}else if(grid->col < 360){
		grid->chaser_col = grid->col/(720.0/128.0) + 64;
	}
		
    /* GlobAlbedo: 2013/01/16 **************************/
    fread(rfdat,sizeof(float),12, fp_s[58]);
    for(e=0;e<ASTEP;e++){
        grid->glbalbedo[e] = rfdat[e];
    }
    
    /* N input *****************************************/
    for(e=0;e<DL_NINPUT;e++){
        grid->nmip_frcrop[e] = 0.0;
        grid->nmip_nfert[e] = 0.0;
        grid->nmip_ndep_noy[e] = 0.0;
        grid->nmip_ndep_nh4[e] = 0.0;
        grid->nmip_manure[e] = 0.0;
        
        for(h=0;h<12;h++){
            grid->nmip_ndep_ccmi_noy[e][h] = 0.0;
            grid->nmip_ndep_ccmi_nh4[e][h] = 0.0;
        }
    }

    if(ISIMIP_RUN == 4){
        /* ISI-MIP2b: 2016/12/24 by A.Ito */
        
        fread(is2bdat,sizeof(float),DL_NINPUT, fp_s[87]);
        for(e=0;e<DL_NINPUT;e++){
            grid->nmip_frcrop[e] = is2bdat[e];
        }
    }else if(EX_BECCS == 2){
        /* BECCS scenario: 2017/02/20 by A.Ito */
        fscanf(fp_s[87],"%lf", &grid->beccs_s2b);
        fscanf(fp_s[87],"%lf", &grid->beccs_v2b);
        fscanf(fp_s[87],"%lf", &grid->beccs_v2s);
    }else if(EX_BECCS == 3){
        /* BECCS scenario: 2017/02/20 by A.Ito */
        fscanf(fp_s[87],"%lf", &grid->beccs_s2b);
        grid->beccs_v2b = 0.0;
        grid->beccs_v2s = 0.0;
    }else{
        /* Bio Fuel scenario: 2015/08/21 by A.Ito ***********/
        for(e=0;e<DL_BF;e++){
            fscanf(fp_s[87],"%lf", &grid->f_biofuel[e]);
            
            if(grid->f_biofuel[e] < 0.0){
                grid->f_biofuel[e] = 0.0;
            }
        }
    }

    if(ISIMIP_RUN == 4){
        fread(is2bdat,sizeof(float),DL_NINPUT, fp_s[88]);
        for(e=0;e<DL_NINPUT;e++){
            /* kg N /ha / yr */
            grid->nmip_nfert[e] = is2bdat[e];
        }
        
    }else{
        /* NMIP input: 2015/11/19 by A.Ito *************/
        for(e=0;e<DL_NINPUT;e++){
            
             /* crop fraction */
            fscanf(fp_s[88],"%lf", &grid->nmip_frcrop[e]);
            if(grid->nmip_frcrop[e] < 0.0){
                grid->nmip_frcrop[e] = 0.0;
            }

            /* g N/ha/yr => kg N/ha/yr */
            fscanf(fp_s[88],"%lf", &grid->nmip_nfert[e]);
            grid->nmip_nfert[e] *= 0.001;
            if(grid->nmip_nfert[e] < 0.0){
                grid->nmip_nfert[e] = 0.0;
            }

            /* g N/ha/yr */
            fscanf(fp_s[88],"%lf", &ddummy);
            if(grid->nmip_nfert[e] < 0.0){
                grid->nmip_nfert[e] = 0.0;
            }

            fscanf(fp_s[88],"%lf", &grid->nmip_manure[e]);
            grid->nmip_manure[e] *= 0.001;
            if(grid->nmip_manure[e] < 0.0){
                grid->nmip_manure[e] = 0.0;
            }

            fscanf(fp_s[88],"%lf", &ddummy);
            if(grid->nmip_manure[e] < 0.0){
                grid->nmip_manure[e] = 0.0;
            }

            for(h=0;h<12;h++){
                fscanf(fp_s[88],"%lf", &grid->nmip_ndep_ccmi_nh4[e][h]);
                if(grid->nmip_ndep_ccmi_nh4[e][h] < 0.0){
                    grid->nmip_ndep_ccmi_nh4[e][h] = 0.0;
                }
            }

            for(h=0;h<12;h++){
                fscanf(fp_s[88],"%lf", &grid->nmip_ndep_ccmi_noy[e][h]);
                if(grid->nmip_ndep_ccmi_noy[e][h] < 0.0){
                    grid->nmip_ndep_ccmi_noy[e][h] = 0.0;
                }
            }
            
        }
    }
    
    if(EX_NFERT == 101){
        /* 101: Nishina ESSD data: 2017/02/13 by A.Ito */
        fscanf(fp_s[89],"%lf %lf", &ddummy, &ddummy);
        
        for(e=0;e<ASTEP;e++){
            fscanf(fp_s[89],"%lf", &ddummy);
            if(ddummy<0.0){ ddummy = 0.0; }
            grid->nin_date[e] = ddummy;
        }
        
        for(e=0;e<50*ASTEP;e++){
            fscanf(fp_s[89],"%lf", &ddummy);
            if(ddummy<0.0){ ddummy = 0.0; }
            grid->nin_no3[e/ASTEP][e%ASTEP] = ddummy;
        }
        
        for(e=0;e<50*ASTEP;e++){
            fscanf(fp_s[89],"%lf", &ddummy);
            if(ddummy<0.0){ ddummy = 0.0; }
            grid->nin_nh4[e/ASTEP][e%ASTEP] = ddummy;
        }
        
        for(e=0;e<90;e++){
            grid->est_nfert[e] = 0.0;
        }
    }else{
        /* future nitrogen fertilizer: 2016/11/22 by A.Ito  */
        for(e=0;e<90;e++){
            /* data: 2010–2099, kg N ha-1 yr-1 */
            fscanf(fp_s[89],"%lf", &nfert_m);
            fscanf(fp_s[89],"%lf", &nfert_r);
            fscanf(fp_s[89],"%lf", &nfert_sw);
            fscanf(fp_s[89],"%lf", &nfert_ww);
            
            if(grid->type_crop == 1){
                
                if(grid->lat >=45.0){
                    /* spring wheat */
                    grid->est_nfert[e] = nfert_sw;
                }else{
                    /* winter wheat */
                    grid->est_nfert[e] = nfert_ww;
                }
                
            }else if(grid->type_crop == 2){
                /* rice */
                grid->est_nfert[e] = nfert_r;
            }else if(grid->type_crop == 3){
                /* maize */
                grid->est_nfert[e] = nfert_m;
            }
            
            if(grid->est_nfert[e] < 0.0){
                grid->est_nfert[e] = 0.0;
            }
        }
        
        for(e=0;e<50*ASTEP;e++){
            grid->nin_no3[e/ASTEP][e%ASTEP] = grid->nin_nh4[e/ASTEP][e%ASTEP] = 0.0;
        }
    }
    
    /* IMPRESSIONS: 2017/05/02 by A.Ito */
    grid->impressions_mask = 0;
    fscanf(fp_s[90],"%ld", &grid->impressions_mask);

    /* Manure input based on Potter: 2017/05/02 by A.Ito */
    grid->impressions_mask = 0;
    fscanf(fp_s[91],"%lf", &grid->nfert_potter);
    fscanf(fp_s[91],"%lf", &grid->nmanure_potter);
}
