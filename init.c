/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* Revised August 15, 2007 by A.Ito */
/* Revised August 19, 2007 by A.Ito */

/* initialize environmental characteristics of each grid */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

#define INT_C 0.01

extern long GCM, CO2S, GCM_R, GCM_C;
extern double aco2_a1[553], aco2_a2[553], aco2_b1[553], aco2_b2[553];
extern double ach4_a1[553], ach4_a2[553], ach4_b1[553], ach4_b2[553];
extern double an2o_a1[553], an2o_a2[553], an2o_b1[553], an2o_b2[553];

extern double glandarea;
extern double h_tmp[201], h_pre[201], h_dswr[201], h_aet[201], h_rof[201];
extern double h_gpp[201], h_npp[201], h_nep[201], h_plant[201], h_soil[201];
extern double h_sr[201], h_ersn_c[201], h_agrersn_c[201], h_doc[201];
extern double h_agrarea[201], h_luc[201];

extern double h_burnt_area[201];
extern double h_bioburn_co2[201], h_bioburn_ch4[201], h_bioburn_co[201];
extern double h_bioburn_nmhc[201], h_bioburn_oc[201], h_bioburn_bc[201];
extern double h_bioburn_nox[201], h_bioburn_so2[201], h_bioburn_pm25[201];
extern double h_bioburn_tpm[201], h_bioburn_tec[201];

extern double h_ch4ox1[201], h_ch4ox2[201], h_ch4ox3[201];
extern double h_ch4emit_cao_paddy[201], h_ch4emit_cao_wetland[201];
extern double h_n2o_emit_ngas[201], h_n2_emit_ngas[201];
extern double h_n2o_emit_casa[201], h_no_emit_casa[201], h_n2_emit_casa[201];
extern double h_nh3_emit[201], h_n2_biofix[201];
extern double h_ch4_emit_mass[201], h_ch4_emit_photo[201];

extern double h_voc_isopr_g97[201], h_voc_monotrp_g97[201], h_voc_methanl_g97[201];
extern double h_voc_acetone_g97[201], h_voc_actaldhd_g97[201], h_voc_frmardhd_g97[201];
extern double h_voc_formacd_g97[201], h_voc_acetacd_g97[201], h_voc_co_g97[201];

/* monthly results **********/
extern double m_ch4ox1[12], m_ch4ox2[12], m_ch4ox3[12];
extern double m_bioburn_co2[12], m_bioburn_ch4[12], m_bioburn_co[12];
extern double m_bioburn_nmhc[12], m_bioburn_oc[12], m_bioburn_bc[12];

/* vegetation (olson) results */
extern double v_area[34];
extern double v_gpp[34], v_npp[34], v_nep[34];
extern double v_lai[34], v_fol[34], v_stm[34], v_rot[34], v_ltr[34], v_msl[34];

/* initialize general parameters used in the simulation ******************/
void initSim(
	struct Grid *grid
){
	long f, year;
	double data;
	FILE *fp_co2;
	
	/* grid ID, sequential number */
	grid->nnn = 0;
	
	/* Number of days for each month */
	grid->mm[0] = 31;
	grid->mm[1] = 28;
	grid->mm[2] = 31;
	grid->mm[3] = 30;
	grid->mm[4] = 31;
	grid->mm[5] = 30;
	grid->mm[6] = 31;
	grid->mm[7] = 31;
	grid->mm[8] = 30;
	grid->mm[9] = 31;
	grid->mm[10] = 30;
	grid->mm[11] = 31;
	
	/* atm. CO2 scenario ****************************/
	printf("reading CO2 data...");
	if(CO2S==1){
		if((fp_co2 = fopen("./data/SRES_A1.dat","rt"))==NULL){
			printf("No SRES_A1.dat\n");
			exit(1);
		}
	}else if(CO2S==2){
		if((fp_co2 = fopen("./data/SRES_A1FI.dat","rt"))==NULL){
			printf("No SRES_A1FI.dat\n");
			exit(1);
		}
	}else if(CO2S==3){
		if((fp_co2 = fopen("./data/SRES_A1T.dat","rt"))==NULL){
			printf("No SRES_A1T.dat\n");
			exit(1);
		}
	}else if(CO2S==4){
		if((fp_co2 = fopen("./data/SRES_A2.dat","rt"))==NULL){
			printf("No SRES_A2.dat\n");
			exit(1);
		}
	}else if(CO2S==5){
		if((fp_co2 = fopen("./data/SRES_B1.dat","rt"))==NULL){
			printf("No SRES_B1.dat\n");
			exit(1);
		}
	}else if(CO2S==6){
		if((fp_co2 = fopen("./data/SRES_B2.dat","rt"))==NULL){
			printf("No SRES_B2.dat\n");
			exit(1);
		}
	}else{
		if((fp_co2 = fopen("./data/SRES_A1.dat","rt"))==NULL){
			printf("No SRES_A1.dat\n");
			exit(1);
		}
	}
	
	/* SRES Scenario CO2 */
	for(f=0;f<111;f++){
		fscanf(fp_co2,"%ld %lf", &year, &data);
		sres_co2[f] = data;
	}
	
	fclose(fp_co2);
	
	/* source: http://crga.atmos.uiuc.edu/research/post-sres.html
		M.E.Schlesinger and S.Malyshev			*/
	if((fp_co2 = fopen("./data/AtmGHG_timeseries.dat","rt"))==NULL){
		printf("No AtmGHG_timeseries.dat\n");
		exit(1);
	}
	for(f=0;f<553;f++){
		fscanf(fp_co2,"%ld", &year);
		/* CO2, ppmv */
		fscanf(fp_co2,"%lf", &aco2_a1[f]);
		fscanf(fp_co2,"%lf", &aco2_a2[f]);
		fscanf(fp_co2,"%lf", &aco2_b1[f]);
		fscanf(fp_co2,"%lf", &aco2_b2[f]);
		/* CH4, ppbv */
		fscanf(fp_co2,"%lf", &ach4_a1[f]);
		fscanf(fp_co2,"%lf", &ach4_a2[f]);
		fscanf(fp_co2,"%lf", &ach4_b1[f]);
		fscanf(fp_co2,"%lf", &ach4_b2[f]);
		/* N2O, ppbv */
		fscanf(fp_co2,"%lf", &an2o_a1[f]);
		fscanf(fp_co2,"%lf", &an2o_a2[f]);
		fscanf(fp_co2,"%lf", &an2o_b1[f]);
		fscanf(fp_co2,"%lf", &an2o_b2[f]);
	}
	fclose(fp_co2);
	
	/* global analysis ********************************/
	glandarea = 0.0;
	for(f=0;f<201;f++){
		h_tmp[f] = h_pre[f] = h_dswr[f] = h_aet[f] = h_rof[f] = 0.0;
		h_gpp[f] = h_npp[f] = h_nep[f] = h_plant[f] = h_soil[f] = 0.0;
		h_sr[f] = h_ersn_c[f] = h_agrersn_c[f] = h_doc[f] = 0.0;
		
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
		h_nh3_emit[f] = h_n2_biofix[f] = 0.0;
		h_ch4_emit_mass[f] = h_ch4_emit_photo[f] = 0.0;

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
}

/********* initialize grid conditions **************************************/
void initG(
	FILE *fp_s[IFILEN], 
	struct Grid *grid
){
	/** initialize climate and soil parameters by using the prepared data files **/
	
	long e, h, g, country, region, aaa;
	double tmp_sfc, tmp_2m,tmp10_soil, tmp200_soil, dswrf_toa, dswrf_sfc, tcdc_clm;
	double prate_sfc, spfh_2m, soilw10, soilw200, ugrd_10m, vgrd_10m;
	double geo_prop, crit_tension;
	double lat, lon, total, wetland, paddy;
	double csiro_lat[56]={	
		 87.5613, 84.4022, 81.2245, 78.0425, 74.8590, -71.6747, 68.4899, 65.3049,
		 62.1197, 58.9343, 55.7489, 52.5634, 49.3779, 46.1924, 43.0068, 39.8211,
		 36.6355, 33.4498, 30.2642, 27.0785, 23.8928, 20.7071, 17.5214, 14.3357,
		 11.1500, 7.9643, 4.7786, 1.5929, -1.5929, -4.7786, -7.9643, -11.1500,
		 -14.3357, -17.5214, -20.7071, -23.8928, -27.0785, -30.2642, -33.4498, -36.6355,
		 -39.8211, -43.0068, -46.1924, -49.3779, -52.5634, -55.7489, -58.9343, -62.1197,
		 -65.3049, -68.4899, -71.6747, -74.8590, -78.0425, -81.2245, -84.4022, -87.5613
	}; 
	double olveg_fc150[34] = {494.0844151, 592.0171207, 570.3302919, 551.9298516, 496.7930014,
		503.2473362, 509.419137, 472.4316637, 449.5906278, 409.8050576, 
		476.5335304, 435.8760416, 457.3711791, 531.5688406, 542.1963318, 
		502.4405676, 563.1641823, 503.0499857, 470.2081472, 504.7536108, 
		323.9272001, 415.1649954, 414.8235868, 531.6733625, 512.1259238, 
		506.076172, 461.7199239, 416.3817974, 347.099887, 551.8461504, 
		513.0419042, 509.2255222, 480.3831019, 245.5416336};
	double olveg_fc30[34] = {128.0402598, 132.4363793, 134.744799, 128.4586405, 126.8214346, 
		130.0745942, 122.2639456, 125.8317647, 121.7714142, 117.3167714, 
		135.7699352, 130.4578967, 134.224135, 127.981914, 130.0432158, 
		123.7658455, 128.9114392, 122.3661531, 119.9684813, 123.5856945, 
		97.54575867, 125.0874554, 127.117925, 127.9560129, 137.5965516, 
		128.8016844, 119.5692298, 111.7970624, 109.7038513, 130.2056391, 
		128.7846698, 127.5222297, 121.7155247, 72.91636457};
	double olveg_bd[34] = {1.319663772, 1.337524685, 1.320005138, 1.425940593, 1.34650667, 
		1.371084061, 1.521147103, 1.310233567, 1.253574467, 1.124569698, 
		1.261084471, 1.141404339, 1.251011779, 1.383010262, 1.377164801, 
		1.474410875, 1.428756115, 1.488748801, 1.36807761, 1.458160626, 
		1.330929888, 1.283214093, 1.18470691, 1.3479263, 1.001494032, 
		1.349477665, 1.454684714, 1.538445206, 1.55359586, 1.426448633, 
		1.353602298, 1.403267173, 1.476433092, 0.899680244};
	double olveg_psand[34] = {46.52418085, 44.63135532, 44.83252566, 49.83793993, 44.40018048, 
		43.08494079, 53.44414746, 44.75929143, 37.19541727, 41.43155107, 
		36.4234834, 44.66832619, 38.57833011, 43.28312248, 41.98256651, 
		53.82325761, 50.19433257, 50.74725112, 49.92285997, 50.09217327, 
		41.71033249, 50.45165556, 46.78736248, 44.00399316, 33.84734969, 
		47.25795129, 50.13222498, 55.09942128, 68.81224605, 42.13913353, 
		38.83132202, 42.81201712, 45.56163789, 82.0};
	double olveg_pclay[34] = {24.21658857, 32.97702758, 29.69877674, 28.76695918, 23.7320935, 
		25.41978957, 23.75750454, 21.6693192, 23.98205022, 23.58236503, 
		25.7041984, 21.80942512, 23.58494412, 24.89991201, 28.20079033, 
		23.55827864, 28.21654069, 21.03711365, 21.2024205, 23.88689401, 
		28.05936266, 19.19199764, 21.56596459, 29.24609601, 27.77211489, 
		23.92326765, 20.40315154, 17.29372406, 12.1077655, 29.1624041, 
		24.38314974, 26.76259603, 23.10797052, 8.0};
	double olveg_ph[34] = {6.112273367, 5.278037138, 5.728076586, 6.052224862, 5.925479094, 
		6.333681185, 7.035558476, 5.977314591, 5.665770821, 5.277109702, 
		5.683146729, 5.525097779, 5.745947716, 5.880433971, 6.278726573, 
		6.841840519, 5.934110876, 7.229348227, 6.176573791, 6.873838987, 
		6.184740765, 6.205232922, 5.915702753, 6.095381045, 5.481843759, 
		6.317546811, 7.379898174, 7.438948581, 7.340133509, 6.045485427, 
		6.722120259, 6.742288787, 7.267506174, 5.759382756};
	double bulkdens, field_cap;
	double x, xx, y, yy, xy;
	
	grid->lat = 89.75-0.5*(double)grid->row;
	grid->lon = -179.75+0.5*(double)grid->col;
	
	/* GCM grid *************************************/
	grid->gcm_row = 0;
	grid->gcm_col = 0;
	/* AR3 */
	if(GCM==1 || GCM==2 || GCM==3 || GCM==4 || GCM==5 || GCM==6){ /* CCSR/NIES */
		grid->gcm_row = grid->row/11.25;	
		grid->gcm_col = grid->col/11.25;
	}else if(GCM==11 || GCM==12 || GCM==13 || GCM==14 || GCM==15 || GCM==16 || GCM==17 || GCM==18){ /* CCCma */
		grid->gcm_row = grid->row/7.5;
		grid->gcm_col = grid->col/7.5;
	}else if(GCM==21 || GCM==22 || GCM==23 || GCM==24){ /* HadCM3 */
		if(grid->lat>=88.75){
			grid->gcm_row = 0;
		}else if(grid->lat<=-88.75){
			grid->gcm_row = 72;
		}else{
			grid->gcm_row = (long)((88.75-grid->lat)/2.5);
		}
		grid->gcm_col = grid->col/7.5;
	}else if(GCM==31 || GCM==32){ /* ECHAM */
		grid->gcm_row = grid->row/5.625;	
		grid->gcm_col = grid->col/5.625;
	}else if(GCM==41 || GCM==42 || GCM==43 || GCM==44){ /* CSIRO */
		for(e=0;e<55;e++){
			if(grid->lat>(csiro_lat[e]+csiro_lat[e+1])/2.0){
				grid->gcm_row = e;
			}
		}
		grid->gcm_col = grid->col/11.25;
	}else if(GCM==51 || GCM==52){ /* GFDL */
		grid->gcm_row = grid->row/4.5;
		grid->gcm_col = grid->col/7.5;
	}else if(GCM==61 || GCM==62 || GCM==63){ /* NCAR-PCM */
		grid->gcm_row = grid->row/5.625;	
		grid->gcm_col = grid->col/5.625;
	}else if(GCM==71){ /* NCAR-CSM */
		grid->gcm_row = grid->row/5.625;	
		grid->gcm_col = grid->col/5.625;
	}
	
	/* AR4 */
	if(GCM==1000 || GCM==1001){ /* MIROC-HIGH */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1010 || GCM==1011 || GCM==1012 || GCM==1013 || GCM==1014 || 
			GCM==1015 || GCM==1016 || GCM==1017 || GCM==1018){ /* MIROC-MED */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1050 || GCM==1051 ||GCM==1052){ /* BCCR */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1061 || GCM==1062 ||GCM==1063){ /* INM-CM3 */
		if(grid->lat>=88.0){
			grid->gcm_row = 0;
		}else if(grid->lat<=-88.0){
			grid->gcm_row = 44;
		}else{
			grid->gcm_row = (long)((88.0-grid->lat)/4.0);
		}
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1070 || GCM==1071 || GCM==1072){ /* IPSL */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1080 || GCM==1081 || GCM==1082){ /* GFDL CM2-1 */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1090 || GCM==1091 || GCM==1092){ /* HadCM3 */
		if(grid->lat>=88.75){
			grid->gcm_row = 0;
		}else if(grid->lat<=-88.75){
			grid->gcm_row = 72;
		}else{
			grid->gcm_row = (long)((88.75-grid->lat)/2.5);
		}
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1100 || GCM==1101 || GCM==1102 || GCM==1103){ /* GISS AOM */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1110 || GCM==1111 || GCM==1112){ /* CCCma-T63 */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1120 || GCM==1121 || GCM==1122){ /* CSIRO */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1130 || GCM==1131 || GCM==1132 || GCM==1133 || GCM==1134
			 || GCM==1135 || GCM==1136 || GCM==1137 || GCM==1138 || GCM==1139
			  || GCM==1140 || GCM==1141 || GCM==1142 || GCM==1143 || GCM==1144){ /* MRI */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1150 ||GCM==1151 ||GCM==1152 ||GCM==1153 ||GCM==1154 ||
			GCM==1155 ||GCM==1156 ||GCM==1157 ||GCM==1158 ||GCM==1159){ /* MPI ECHAM */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1160 || GCM==1161 || GCM==1162 || GCM==1163 || GCM==1164 || GCM==1165){ /* IAP */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1170 || GCM==1171){ /* HadGEM */
		if(grid->lat>=89.379){
			grid->gcm_row = 0;
		}else if(grid->lat<=-89.379){
			grid->gcm_row = 144;
		}else{
			grid->gcm_row = (long)((89.379-grid->lat)/(180.0/(double)GCM_R));
		}
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1180 || GCM==1181 || GCM==1182){ /* GFDL CM 2.0 */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1190 || GCM==1191 || GCM==1192){ /* CNRM */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);	
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1200 || GCM==1201 || GCM==1202 || GCM==1203 || GCM==1204
			 || GCM==1205 || GCM==1206 || GCM==1207 || GCM==1208 || GCM==1209
			 || GCM==1210 || GCM==1211 || GCM==1212 || GCM==1213 || GCM==1214){ /* CCC47 */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1220 || GCM==1221 || GCM==1222 || GCM==1223 || GCM==1224
			 || GCM==1225 || GCM==1226 || GCM==1227 || GCM==1228 || GCM==1229
			  || GCM==1230 || GCM==1231 || GCM==1232 || GCM==1233 || GCM==1234
			   || GCM==1235 || GCM==1236 || GCM==1237 || GCM==1238 || GCM==1239){ /* NCAR CCSM */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1240 || GCM==1241 || GCM==1242){ /* GISS E */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1250 || GCM==1251 || GCM==1252 || GCM==1253){ /* GISS R */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}else if(GCM==1260 || GCM==1261 || GCM==1262 || GCM==1263 || GCM==1264
			 || GCM==1265 || GCM==1266 || GCM==1267){ /* NCAR PCM */
		grid->gcm_row = grid->row/(360.0/(double)GCM_R);
		grid->gcm_col = grid->col/(720.0/(double)GCM_C);
	}

	/** input geography in the grid **/	
	fscanf(fp_s[0],"%ld %ld", &country, &region); 
	grid->country = country; /* country code */
	grid->region = region; /* region code */
	
	/* printf("reading Veg data..."); */
	/** input vegetation area **/
	fscanf(fp_s[1],"%ld", &aaa); 
	grid->veg_olson = aaa;
	/* printf("done\n"); */
	
	/* sensitivity analysis for deforestation */
	if(DEFOREST==1){
		if(aaa>=1&&aaa<=12){
			grid->veg_olson = 19;
		}
	}else if(DEFOREST==2){
		if(aaa>=1&&aaa<=12){
			grid->veg_olson = 13;
		}
	}else if(DEFOREST==3){
		if(aaa>=1&&aaa<=12){
			grid->veg_olson = 31;
		}
	}
	
	/* printf("reading Soil data..."); */
	/** input soil properties *******************************/
	fscanf(fp_s[3],"%lf", &geo_prop); 
		geo_prop = (geo_prop>0.0)?geo_prop:0.0;
		grid->topo = geo_prop; /* topography */
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
	/* printf("done\n"); */
	
	/* printf("reading Base climate data..."); */
	/************ input monthly climate ************/	
	/* long-term average, from NCEP/NCAR reanalysis */
	for(e=0;e<12;e++){
		fscanf(fp_s[2],"%lf", &tmp_sfc); 
		grid->tmp_sfc_a[e] = tmp_sfc-ZAT;
	}	
	for(e=0;e<12;e++){
		fscanf(fp_s[2],"%lf", &tmp_2m); 
		grid->tmp_2m_a[e] = tmp_2m-ZAT;
	}	
	for(e=0;e<12;e++){
		fscanf(fp_s[2],"%lf", &tmp10_soil); 
		grid->tmp10_soil_a[e] = tmp10_soil-ZAT;
	}	
	for(e=0;e<12;e++){
		fscanf(fp_s[2],"%lf", &tmp200_soil); 
		grid->tmp200_soil_a[e] = tmp200_soil-ZAT;
	}	
	for(e=0;e<12;e++){
		fscanf(fp_s[2],"%lf", &dswrf_toa); 
	}	
	for(e=0;e<12;e++){
		fscanf(fp_s[2],"%lf", &dswrf_sfc); 
	}	
	for(e=0;e<12;e++){
		fscanf(fp_s[2],"%lf", &tcdc_clm); 
		grid->tcdc_clm_a[e] = tcdc_clm/100.0;
	}	
	for(e=0;e<12;e++){
		fscanf(fp_s[2],"%lf", &prate_sfc); 
		grid->prate_sfc_a[e] = prate_sfc*(double)(grid->mm[e]);
	}	
	for(e=0;e<12;e++){
		fscanf(fp_s[2],"%lf", &spfh_2m); 
		grid->spfh_2m_a[e] = spfh_2m;
	}	
	for(e=0;e<12;e++){
		fscanf(fp_s[2],"%lf", &soilw10); 
	}	
	for(e=0;e<12;e++){
		fscanf(fp_s[2],"%lf", &soilw200); 
	}	
	for(e=0;e<12;e++){
		fscanf(fp_s[2],"%lf", &ugrd_10m); 
		grid->ugrd_10m_a[e] = ugrd_10m;
	}	
	for(e=0;e<12;e++){
		fscanf(fp_s[2],"%lf", &vgrd_10m); 
		grid->vgrd_10m_a[e] = vgrd_10m;
	}
	/* printf("done\n"); */
	
	/* init CRU memory */
	for(h=0;h<CRU_PD;h++){
		for(g=0;g<12;g++){
			grid->hist_cld[h][g] = grid->tcdc_clm_a[g];
			grid->hist_pre[h][g] = grid->prate_sfc_a[g];
			grid->hist_tmp[h][g] = grid->tmp_2m_a[g];
		}
	}
	for(g=0;g<12;g++){
		grid->hist_cld_b[g] = grid->tcdc_clm_a[g];
		grid->hist_pre_b[g] = grid->prate_sfc_a[g];
		grid->hist_tmp_b[g] = grid->tmp_2m_a[g];
	}
	
	/* substituted precipitation data (UEA/CRU) */
	for(e=0;e<12;e++){
		fscanf(fp_s[4],"%lf", &prate_sfc); 
		grid->prec_sub_a[e] = prate_sfc*(double)(grid->mm[e]);
	}
	
	/* river-basin ID */
	fscanf(fp_s[5],"%ld", &(grid->rvbasin)); 
	
	/* albedo */
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
	
	/* Cropland coverage by Ramankutty & Foley (1999) */
	for(h=0;h<293;h++){
		fscanf(fp_s[11],"%lf", &(grid->fcrop_sage[h])); 
	}
	
	x = xx = y = yy = xy = 0.0;
	for(h=0;h<10;h++){
		x += (double)h;
		xx += (double)h * (double)h;
		y += grid->fcrop_sage[280+h];
		yy += grid->fcrop_sage[280+h]*grid->fcrop_sage[280+h];
		xy += (double)h*grid->fcrop_sage[280+h];
	}

	grid->f_crop_trend = (10.0*xy - x*y)/(10.0*xx - x*x);
	
	/* future land-use change scenario by IMAGE2 (Wang et al. 2006) */
	/*
	Wang, A., Price, D.T. and Arora, V., 2006. Estimating changes in global 
	vegetation cover (1850-2100) for use in climate models. 
	Global Biogeochemical Cycles 20, 10.1029/2005GB002514.
	*/
	for(h=0;h<111;h++){
		fscanf(fp_s[23],"%lf", &(grid->fcrop3_image[h])); 
	}
	for(h=0;h<111;h++){
		fscanf(fp_s[23],"%lf", &(grid->fcrop4_image[h])); 
	}
	for(h=0;h<111;h++){
		fscanf(fp_s[45],"%lf", &(grid->fgrass3_image[h])); 
	}
	for(h=0;h<111;h++){
		fscanf(fp_s[45],"%lf", &(grid->fgrass4_image[h])); 
	}
	
	fscanf(fp_s[12],"%ld", &(grid->soiltexture)); 
	/* 0: ocean */
	/* 1: coarse */
	/* 2: medium/coarse */
	/* 3: medium */
	/* 4: fine/medium */
	/* 5: fine */
	/* 6: ice */
	/* 7: organic */
	
	fscanf(fp_s[13],"%lf", &bulkdens); 
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
	grid->a_sw = exp(-4.396 -0.0715*grid->pc_clay -4.488*0.0001*grid->pc_clay*grid->pc_clay 
					- 4.285*0.00001*grid->pc_sand*grid->pc_sand*grid->pc_clay)*100;
	grid->b_sw = -3.14 -0.00222*grid->pc_clay*grid->pc_clay 
					-3.484*0.00001*grid->pc_sand*grid->pc_sand*grid->pc_clay;
	
	/* field capacity */
	grid->field_cap1 = pow(crit_tension/grid->a_sw, 1.0/grid->b_sw);
	/* pore space capacity */
	grid->pore_cap1 = 0.332 -7.251*0.0001*grid->pc_sand + 0.1276*log10(grid->pc_clay);
	
	grid->field_cap1 *= 300.0;
	grid->field_cap2 = grid->fc_150 - grid->field_cap1;
	if(grid->field_cap2<0.0){
		grid->field_cap2 = 100.0;
	}
	grid->pore_cap1 *= 300.0;
	
	/* wetland fraction: data by Global Lakes and Wetlands Database */
	fscanf(fp_s[21],"%lf %lf", &total, &wetland); 
	wetland *= 1.0/1000000.0;
	grid->f_wetland = wetland/grid->area;
	if(grid->f_wetland > 1.0){
		grid->f_wetland = 1.0;
	}
	
	/* paddy fraction: data by U.Wisconsin SAGE (Leff et al.) */
	fscanf(fp_s[22],"%lf", &paddy); 
	if(paddy>0.0){
		grid->f_paddy = paddy;
	}else{
		grid->f_paddy = 0.0;
	}
	
	grid->f_upland = 1.0 - grid->f_wetland - grid->f_paddy;
	if(grid->f_upland < 0.0){
		grid->f_upland = 0.0;
	}
	
	/* soil total nitrogen */
	fscanf(fp_s[24],"%lf", &grid->total_n_1m); 
	if(grid->total_n_1m < 0.0){
		grid->total_n_1m = 0.0;
	}

	fscanf(fp_s[25],"%lf", &lat); 
	fscanf(fp_s[25],"%lf", &lon); 
	fscanf(fp_s[25],"%lf", &grid->ndepo[0]); 
	fscanf(fp_s[25],"%lf", &grid->ndepo[1]); 
	fscanf(fp_s[25],"%lf", &grid->ndepo[2]); 
	
	/* N fertilization */
	grid->nfert_nh4 = 0.0;
	grid->nfert_no3 = 0.0;
	
	/* EOS-WEBSTER Land-use change data */
	for(h=0;h<301;h++){
		fscanf(fp_s[26],"%lf", &grid->fcrop_eossagehyde[h]);
		fscanf(fp_s[27],"%lf", &grid->fpast_eossagehyde[h]);
		fscanf(fp_s[28],"%lf", &grid->fprim_eossagehyde[h]);
		fscanf(fp_s[29],"%lf", &grid->fsecd_eossagehyde[h]);
		fscanf(fp_s[30],"%lf", &grid->ssma_eossagehyde[h]);
		fscanf(fp_s[31],"%lf", &grid->ssmb_eossagehyde[h]);
		fscanf(fp_s[32],"%lf", &grid->t_cp_eossagehyde[h]);
		fscanf(fp_s[33],"%lf", &grid->t_cs_eossagehyde[h]);
		fscanf(fp_s[34],"%lf", &grid->t_pc_eossagehyde[h]);
		fscanf(fp_s[35],"%lf", &grid->t_ps_eossagehyde[h]);
		fscanf(fp_s[36],"%lf", &grid->t_sc_eossagehyde[h]);
		fscanf(fp_s[37],"%lf", &grid->t_sp_eossagehyde[h]);
		fscanf(fp_s[38],"%lf", &grid->t_ss1_eossagehyde[h]);
		fscanf(fp_s[39],"%lf", &grid->t_ss2_eossagehyde[h]);
		fscanf(fp_s[40],"%lf", &grid->t_ss3_eossagehyde[h]);
		fscanf(fp_s[41],"%lf", &grid->t_vc_eossagehyde[h]);
		fscanf(fp_s[42],"%lf", &grid->t_vp_eossagehyde[h]);
		fscanf(fp_s[43],"%lf", &grid->t_vs1_eossagehyde[h]);
		fscanf(fp_s[44],"%lf", &grid->t_vs2_eossagehyde[h]);
	}
}
