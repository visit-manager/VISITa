/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* read AOGCM output ***********************************/
void read_gcm_clim(
	struct Grid  *grid
){
	long f, g, h, i;
	long yr, mon;
	double atmp, apres, shum, alt;
	extern double MDN[12];
	FILE *fp_t, *fp_p, *fp_h, *fp_r;
	
	if(GCM==1){
		if( (fp_t=fopen("CN_A1_T.dat","rt"))==NULL ){  printf("No CN_A1_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CN_A1_P.dat","rt"))==NULL ){  printf("No CN_A1_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CN_A1_H.dat","rt"))==NULL ){  printf("No CN_A1_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CN_A1_R.dat","rt"))==NULL ){  printf("No CN_A1_R.dat\n");  exit(1); }
	}else if(GCM==2){
		if( (fp_t=fopen("CN_A1FI_T.dat","rt"))==NULL ){  printf("No CN_A1FI_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CN_A1FI_P.dat","rt"))==NULL ){  printf("No CN_A1FI_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CN_A1FI_H.dat","rt"))==NULL ){  printf("No CN_A1FI_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CN_A1FI_R.dat","rt"))==NULL ){  printf("No CN_A1FI_R.dat\n");  exit(1); }
	}else if(GCM==3){
		if( (fp_t=fopen("CN_A1T_T.dat","rt"))==NULL ){  printf("No CN_A1T_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CN_A1T_P.dat","rt"))==NULL ){  printf("No CN_A1T_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CN_A1T_H.dat","rt"))==NULL ){  printf("No CN_A1T_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CN_A1T_R.dat","rt"))==NULL ){  printf("No CN_A1T_R.dat\n");  exit(1); }
	}else if(GCM==4){
		if( (fp_t=fopen("CN_A2_T.dat","rt"))==NULL ){  printf("No CN_A2_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CN_A2_P.dat","rt"))==NULL ){  printf("No CN_A2_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CN_A2_H.dat","rt"))==NULL ){  printf("No CN_A2_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CN_A2_R.dat","rt"))==NULL ){  printf("No CN_A2_R.dat\n");  exit(1); }
	}else if(GCM==5){
		if( (fp_t=fopen("CN_B1_T.dat","rt"))==NULL ){  printf("No CN_B1_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CN_B1_P.dat","rt"))==NULL ){  printf("No CN_B1_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CN_B1_H.dat","rt"))==NULL ){  printf("No CN_B1_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CN_B1_R.dat","rt"))==NULL ){  printf("No CN_B1_R.dat\n");  exit(1); }
	}else if(GCM==6){
		if( (fp_t=fopen("CN_B2_T.dat","rt"))==NULL ){  printf("No CN_B2_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CN_B2_P.dat","rt"))==NULL ){  printf("No CN_B2_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CN_B2_H.dat","rt"))==NULL ){  printf("No CN_B2_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CN_B2_R.dat","rt"))==NULL ){  printf("No CN_B2_R.dat\n");  exit(1); }
	}else if(GCM==11){
		if( (fp_t=fopen("CC_A2_T.dat","rt"))==NULL ){  printf("No CC_A2_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CC_A2_P.dat","rt"))==NULL ){  printf("No CC_A2_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CC_A2_H.dat","rt"))==NULL ){  printf("No CC_A2_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CC_A2_R.dat","rt"))==NULL ){  printf("No CC_A2_R.dat\n");  exit(1); }
	}else if(GCM==12){
		if( (fp_t=fopen("CC_B2_T.dat","rt"))==NULL ){  printf("No CC_B2_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CC_B2_P.dat","rt"))==NULL ){  printf("No CC_B2_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CC_B2_H.dat","rt"))==NULL ){  printf("No CC_B2_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CC_B2_R.dat","rt"))==NULL ){  printf("No CC_B2_R.dat\n");  exit(1); }
	}else if(GCM==13){
		if( (fp_t=fopen("CC_A2a_T.dat","rt"))==NULL ){  printf("No CC_A2a_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CC_A2a_P.dat","rt"))==NULL ){  printf("No CC_A2a_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CC_A2a_H.dat","rt"))==NULL ){  printf("No CC_A2a_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CC_A2a_R.dat","rt"))==NULL ){  printf("No CC_A2a_R.dat\n");  exit(1); }
	}else if(GCM==14){
		if( (fp_t=fopen("CC_A2b_T.dat","rt"))==NULL ){  printf("No CC_A2b_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CC_A2b_P.dat","rt"))==NULL ){  printf("No CC_A2b_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CC_A2b_H.dat","rt"))==NULL ){  printf("No CC_A2b_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CC_A2b_R.dat","rt"))==NULL ){  printf("No CC_A2b_R.dat\n");  exit(1); }
	}else if(GCM==15){
		if( (fp_t=fopen("CC_A2c_T.dat","rt"))==NULL ){  printf("No CC_A2c_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CC_A2c_P.dat","rt"))==NULL ){  printf("No CC_A2c_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CC_A2c_H.dat","rt"))==NULL ){  printf("No CC_A2c_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CC_A2c_R.dat","rt"))==NULL ){  printf("No CC_A2c_R.dat\n");  exit(1); }
	}else if(GCM==16){
		if( (fp_t=fopen("CC_B2a_T.dat","rt"))==NULL ){  printf("No CC_B2a_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CC_B2a_P.dat","rt"))==NULL ){  printf("No CC_B2a_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CC_B2a_H.dat","rt"))==NULL ){  printf("No CC_B2a_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CC_B2a_R.dat","rt"))==NULL ){  printf("No CC_B2a_R.dat\n");  exit(1); }
	}else if(GCM==17){
		if( (fp_t=fopen("CC_B2b_T.dat","rt"))==NULL ){  printf("No CC_B2b_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CC_B2b_P.dat","rt"))==NULL ){  printf("No CC_B2b_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CC_B2b_H.dat","rt"))==NULL ){  printf("No CC_B2b_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CC_B2b_R.dat","rt"))==NULL ){  printf("No CC_B2b_R.dat\n");  exit(1); }
	}else if(GCM==18){
		if( (fp_t=fopen("CC_B2c_T.dat","rt"))==NULL ){  printf("No CC_B2c_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CC_B2c_P.dat","rt"))==NULL ){  printf("No CC_B2c_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CC_B2c_H.dat","rt"))==NULL ){  printf("No CC_B2c_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CC_B2c_R.dat","rt"))==NULL ){  printf("No CC_B2c_R.dat\n");  exit(1); }
	}else if(GCM==21){
		if( (fp_t=fopen("HD_A2_T.dat","rt"))==NULL ){  printf("No HD_A2_T.dat\n");  exit(1); }
		if( (fp_p=fopen("HD_A2_P.dat","rt"))==NULL ){  printf("No HD_A2_P.dat\n");  exit(1); }
		if( (fp_h=fopen("HD_A2_H.dat","rt"))==NULL ){  printf("No HD_A2_H.dat\n");  exit(1); }
		if( (fp_r=fopen("HD_A2_R.dat","rt"))==NULL ){  printf("No HD_A2_R.dat\n");  exit(1); }
	}else if(GCM==22){
		if( (fp_t=fopen("HD_A2b_T.dat","rt"))==NULL ){  printf("No HD_A2b_T.dat\n");  exit(1); }
		if( (fp_p=fopen("HD_A2b_P.dat","rt"))==NULL ){  printf("No HD_A2b_P.dat\n");  exit(1); }
		if( (fp_h=fopen("HD_A2b_H.dat","rt"))==NULL ){  printf("No HD_A2b_H.dat\n");  exit(1); }
		if( (fp_r=fopen("HD_A2b_R.dat","rt"))==NULL ){  printf("No HD_A2b_R.dat\n");  exit(1); }
	}else if(GCM==23){
		if( (fp_t=fopen("HD_A2c_T.dat","rt"))==NULL ){  printf("No HD_A2c_T.dat\n");  exit(1); }
		if( (fp_p=fopen("HD_A2c_P.dat","rt"))==NULL ){  printf("No HD_A2c_P.dat\n");  exit(1); }
		if( (fp_h=fopen("HD_A2c_H.dat","rt"))==NULL ){  printf("No HD_A2c_H.dat\n");  exit(1); }
		if( (fp_r=fopen("HD_A2c_R.dat","rt"))==NULL ){  printf("No HD_A2c_R.dat\n");  exit(1); }
	}else if(GCM==24){
		if( (fp_t=fopen("HD_B2_T.dat","rt"))==NULL ){  printf("No HD_B2_T.dat\n");  exit(1); }
		if( (fp_p=fopen("HD_B2_P.dat","rt"))==NULL ){  printf("No HD_B2_P.dat\n");  exit(1); }
		if( (fp_h=fopen("HD_B2_H.dat","rt"))==NULL ){  printf("No HD_B2_H.dat\n");  exit(1); }
		if( (fp_r=fopen("HD_B2_R.dat","rt"))==NULL ){  printf("No HD_B2_R.dat\n");  exit(1); }
	}else if(GCM==31){
		if( (fp_t=fopen("EH_A2_T.dat","rt"))==NULL ){  printf("No EH_A2_T.dat\n");  exit(1); }
		if( (fp_p=fopen("EH_A2_P.dat","rt"))==NULL ){  printf("No EH_A2_P.dat\n");  exit(1); }
		if( (fp_h=fopen("EH_A2_H.dat","rt"))==NULL ){  printf("No EH_A2_H.dat\n");  exit(1); }
		if( (fp_r=fopen("EH_A2_R.dat","rt"))==NULL ){  printf("No EH_A2_R.dat\n");  exit(1); }
	}else if(GCM==32){
		if( (fp_t=fopen("EH_B2_T.dat","rt"))==NULL ){  printf("No EH_B2_T.dat\n");  exit(1); }
		if( (fp_p=fopen("EH_B2_P.dat","rt"))==NULL ){  printf("No EH_B2_P.dat\n");  exit(1); }
		if( (fp_h=fopen("EH_B2_H.dat","rt"))==NULL ){  printf("No EH_B2_H.dat\n");  exit(1); }
		if( (fp_r=fopen("EH_B2_R.dat","rt"))==NULL ){  printf("No EH_B2_R.dat\n");  exit(1); }
	}else if(GCM==41){
		if( (fp_t=fopen("CS_A1_T.dat","rt"))==NULL ){  printf("No CS_A1_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CS_A1_P.dat","rt"))==NULL ){  printf("No CS_A1_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CS_A1_H.dat","rt"))==NULL ){  printf("No CS_A1_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CS_A1_R.dat","rt"))==NULL ){  printf("No CS_A1_R.dat\n");  exit(1); }
	}else if(GCM==42){
		if( (fp_t=fopen("CS_A2_T.dat","rt"))==NULL ){  printf("No CS_A2_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CS_A2_P.dat","rt"))==NULL ){  printf("No CS_A2_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CS_A2_H.dat","rt"))==NULL ){  printf("No CS_A2_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CS_A2_R.dat","rt"))==NULL ){  printf("No CS_A2_R.dat\n");  exit(1); }
	}else if(GCM==43){
		if( (fp_t=fopen("CS_B1_T.dat","rt"))==NULL ){  printf("No CS_B1_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CS_B1_P.dat","rt"))==NULL ){  printf("No CS_B1_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CS_B1_H.dat","rt"))==NULL ){  printf("No CS_B1_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CS_B1_R.dat","rt"))==NULL ){  printf("No CS_B1_R.dat\n");  exit(1); }
	}else if(GCM==44){
		if( (fp_t=fopen("CS_B2_T.dat","rt"))==NULL ){  printf("No CS_B2_T.dat\n");  exit(1); }
		if( (fp_p=fopen("CS_B2_P.dat","rt"))==NULL ){  printf("No CS_B2_P.dat\n");  exit(1); }
		if( (fp_h=fopen("CS_B2_H.dat","rt"))==NULL ){  printf("No CS_B2_H.dat\n");  exit(1); }
		if( (fp_r=fopen("CS_B2_R.dat","rt"))==NULL ){  printf("No CS_B2_R.dat\n");  exit(1); }
	}else if(GCM==51){
		if( (fp_t=fopen("GF_A2_T.dat","rt"))==NULL ){  printf("No GF_A2_T.dat\n");  exit(1); }
		if( (fp_p=fopen("GF_A2_P.dat","rt"))==NULL ){  printf("No GF_A2_P.dat\n");  exit(1); }
		if( (fp_h=fopen("GF_A2_H.dat","rt"))==NULL ){  printf("No GF_A2_H.dat\n");  exit(1); }
		if( (fp_r=fopen("GF_A2_R.dat","rt"))==NULL ){  printf("No GF_A2_R.dat\n");  exit(1); }
	}else if(GCM==52){
		if( (fp_t=fopen("GF_B2_T.dat","rt"))==NULL ){  printf("No GF_B2_T.dat\n");  exit(1); }
		if( (fp_p=fopen("GF_B2_P.dat","rt"))==NULL ){  printf("No GF_B2_P.dat\n");  exit(1); }
		if( (fp_h=fopen("GF_B2_H.dat","rt"))==NULL ){  printf("No GF_B2_H.dat\n");  exit(1); }
		if( (fp_r=fopen("GF_B2_R.dat","rt"))==NULL ){  printf("No GF_B2_R.dat\n");  exit(1); }
	}else if(GCM==61){
		if( (fp_t=fopen("NP_A1B_T.dat","rt"))==NULL ){  printf("No NP_A1B_T.dat\n");  exit(1); }
		if( (fp_p=fopen("NP_A1B_P.dat","rt"))==NULL ){  printf("No NP_A1B_P.dat\n");  exit(1); }
		if( (fp_h=fopen("NP_A1B_H.dat","rt"))==NULL ){  printf("No NP_A1B_H.dat\n");  exit(1); }
		if( (fp_r=fopen("NP_A1B_R.dat","rt"))==NULL ){  printf("No NP_A1B_R.dat\n");  exit(1); }
	}else if(GCM==62){
		if( (fp_t=fopen("NP_A2_T.dat","rt"))==NULL ){  printf("No NP_A2_T.dat\n");  exit(1); }
		if( (fp_p=fopen("NP_A2_P.dat","rt"))==NULL ){  printf("No NP_A2_P.dat\n");  exit(1); }
		if( (fp_h=fopen("NP_A2_H.dat","rt"))==NULL ){  printf("No NP_A2_H.dat\n");  exit(1); }
		if( (fp_r=fopen("NP_A2_R.dat","rt"))==NULL ){  printf("No NP_A2_R.dat\n");  exit(1); }
	}else if(GCM==63){
		if( (fp_t=fopen("NP_B2_T.dat","rt"))==NULL ){  printf("No NP_B2_T.dat\n");  exit(1); }
		if( (fp_p=fopen("NP_B2_P.dat","rt"))==NULL ){  printf("No NP_B2_P.dat\n");  exit(1); }
		if( (fp_h=fopen("NP_B2_H.dat","rt"))==NULL ){  printf("No NP_B2_H.dat\n");  exit(1); }
		if( (fp_r=fopen("NP_B2_R.dat","rt"))==NULL ){  printf("No NP_B2_R.dat\n");  exit(1); }
	}else if(GCM==71){
		if( (fp_t=fopen("NC_A2_T.dat","rt"))==NULL ){  printf("No NC_A2_T.dat\n");  exit(1); }
		if( (fp_p=fopen("NC_A2_P.dat","rt"))==NULL ){  printf("No NC_A2_P.dat\n");  exit(1); }
		if( (fp_h=fopen("NP_A2_H.dat","rt"))==NULL ){  printf("No NP_A2_H.dat\n");  exit(1); }  /*  alternative data */
		if( (fp_r=fopen("NC_A2_R.dat","rt"))==NULL ){  printf("No NC_A2_R.dat\n");  exit(1); }
	}
	
	if(GCM==1000){  /* MIROC-HIGH + A1B  */
		if( (fp_t=fopen("./data/tas_MIROCH_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MIROCH_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MIROCH_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MIROCH_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1001){  /* MIROC-HIGH + B1  */
		if( (fp_t=fopen("./data/tas_MIROCH_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_tas_B1-2000-2099-H.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MIROCH_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_pr_B1-2000-2099-H.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MIROCH_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_huss_B1-2000-2099-H.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MIROCH_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_rsds_B1-2000-2099-H.dat\n");  exit(1); }
	}
	
	if(GCM==1010){  /* MIROC-MED + A1B Run1 */
		if( (fp_t=fopen("./data/tas_MIROCM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_MIROCM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MIROCM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_MIROCM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MIROCM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_MIROCM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MIROCM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_MIROCM_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1011){  /* MIROC-MED + A1B Run2 */
		if( (fp_t=fopen("./data/tas_MIROCM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No tas_MIROCM_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MIROCM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No pr_MIROCM_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MIROCM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No huss_MIROCM_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MIROCM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No rsds_MIROCM_20C-A1B_R2.dat\n");  exit(1); }
	}else if(GCM==1012){  /* MIROC-MED + A1B Run3 */
		if( (fp_t=fopen("./data/tas_MIROCM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No tas_MIROCM_20C-A1B_R3.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MIROCM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No pr_MIROCM_20C-A1B_R3.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MIROCM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No huss_MIROCM_20C-A1B_R3.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MIROCM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No rsds_MIROCM_20C-A1B_R3.dat\n");  exit(1); }
	}else if(GCM==1013){  /* MIROC-MED + A2  */
		if( (fp_t=fopen("./data/tas_MIROCM_20C-A2_R1.dat","rt"))==NULL ){  printf("No tas_MIROCM_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MIROCM_20C-A2_R1.dat","rt"))==NULL ){  printf("No pr_MIROCM_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MIROCM_20C-A2_R1.dat","rt"))==NULL ){  printf("No huss_MIROCM_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MIROCM_20C-A2_R1.dat","rt"))==NULL ){  printf("No rsds_MIROCM_20C-A2_R1.dat\n");  exit(1); }
	}else if(GCM==1014){  /* MIROC-MED + A2  */
		if( (fp_t=fopen("./data/tas_MIROCM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_tas_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MIROCM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_pr_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MIROCM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_huss_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MIROCM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_rsds_A1B-2000-2099-R2.dat\n");  exit(1); }
	}else if(GCM==1015){  /* MIROC-MED + A2  */
		if( (fp_t=fopen("./data/tas_MIROCM_20C-A2_R3.dat","rt"))==NULL ){  printf("No tas_MIROCM_20C-A2_R3.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MIROCM_20C-A2_R3.dat","rt"))==NULL ){  printf("No pr_MIROCM_20C-A2_R3.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MIROCM_20C-A2_R3.dat","rt"))==NULL ){  printf("No huss_MIROCM_20C-A2_R3.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MIROCM_20C-A2_R3.dat","rt"))==NULL ){  printf("No rsds_MIROCM_20C-A2_R3.dat\n");  exit(1); }
	}else if(GCM==1016){  /* MIROC-MED + B1  */
		if( (fp_t=fopen("./data/tas_MIROCM_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_MIROCM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MIROCM_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_MIROCM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MIROCM_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_MIROCM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MIROCM_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_MIROCM_20C-B1_R1.dat\n");  exit(1); }
	}else if(GCM==1017){  /* MIROC-MED + B1  */
		if( (fp_t=fopen("./data/tas_MIROCM_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_tas_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MIROCM_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_pr_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MIROCM_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_huss_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MIROCM_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_rsds_A1B-2000-2099-R2.dat\n");  exit(1); }
	}else if(GCM==1018){  /* MIROC-MED + B1  */
		if( (fp_t=fopen("./data/tas_MIROCM_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_tas_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MIROCM_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_pr_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MIROCM_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_huss_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MIROCM_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_rsds_A1B-2000-2099-R2.dat\n");  exit(1); }
	}
	
	if(GCM==1050){  /* BCCR + A1B  */
		if( (fp_t=fopen("./data/tas_BCCR_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_BCCR_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_BCCR_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_BCCR_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_BCCR_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_BCCR_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_BCCR_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_BCCR_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1051){  /* BCCR + A2  */
		if( (fp_t=fopen("./data/tas_BCCR_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_BRA1B-2001-2100\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_BCCR_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_BRA1B-2001-2100\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_BCCR_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_BRA1B-2001-2100\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_BCCR_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_BRA1B-2001-2100\n");  exit(1); }
	}else if(GCM==1052){  /* BCCR + A1B  */
		if( (fp_t=fopen("./data/tas_BCCR_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_tas_BRA1B-2001-2100\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_BCCR_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_pr_BRA1B-2001-2100\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_BCCR_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_huss_BRA1B-2001-2100\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_BCCR_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_rsds_BRA1B-2001-2100\n");  exit(1); }
	}
	
	if(GCM==1060){  /* INM-CM3 + A1B  */
		if( (fp_t=fopen("./data/tas_INM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_INM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_INM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_INM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_INM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_INM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_INM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_INM_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1061){  /* INM-CM3 + A1B  */
		if( (fp_t=fopen("./data/tas_INM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_INA1B-2001-2200\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_INM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_INA1B-2001-2200\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_INM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_INA1B-2001-2200\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_INM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_INA1B-2001-2200\n");  exit(1); }
	}else if(GCM==1062){  /* INM-CM3 + A1B  */
		if( (fp_t=fopen("./data/tas_INM_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_INM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_INM_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_INM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_INM_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_INM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_INM_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_INM_20C-B1_R1.dat\n");  exit(1); }
	}
	
	if(GCM==1070){  /* IPSL + A1B  */
		if( (fp_t=fopen("./data/tas_IPSL_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_IPSL_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_IPSL_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_IPSL_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_IPSL_20C-A2_R1.dat\n");  exit(1); }
	}else if(GCM==1071){  /* IPSL + A2  */
		if( (fp_t=fopen("./data/tas_IPSL_20C-A2_R1.dat","rt"))==NULL ){  printf("No tas_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_IPSL_20C-A2_R1.dat","rt"))==NULL ){  printf("No pr_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_IPSL_20C-A2_R1.dat","rt"))==NULL ){  printf("No huss_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_IPSL_20C-A2_R1.dat","rt"))==NULL ){  printf("No rsds_IPSL_20C-A2_R1.dat\n");  exit(1); }
	}else if(GCM==1072){  /* IPSL + B1  */
		if( (fp_t=fopen("./data/tas_IPSL_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_IPSL_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_IPSL_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_IPSL_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_IPSL_20C-A2_R1.dat\n");  exit(1); }
	}
	
	if(GCM==1080){  /* GFDL CM2-1 + A1B  */
		if( (fp_t=fopen("./data/tas_GFDL21_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_tas_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GFDL21_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_pr_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GFDL21_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_hus_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GFDL21_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_rsds_GFA1B-2001-2100\n");  exit(1); }
	}else if(GCM==1081){  /* GFDL CM2-1 + A2  */
		if( (fp_t=fopen("./data/tas_GFDL21_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GFDL21_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GFDL21_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_hus_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GFDL21_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_GFA1B-2001-2100\n");  exit(1); }
	}else if(GCM==1082){  /* GFDL CM2-1 + B1  */
		if( (fp_t=fopen("./data/tas_GFDL21_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_tas_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GFDL21_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_pr_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GFDL21_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_hus_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GFDL21_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_rsds_GFA1B-2001-2100\n");  exit(1); }
	}
	
	if(GCM==1090){  /* HadCM3 + A1B  */
		if( (fp_t=fopen("./data/tas_HADCM3_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_tas_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_HADCM3_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_pr_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_HADCM3_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_huss_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_HADCM3_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_rsds_HCA1B-2000-2200\n");  exit(1); }
	}else if(GCM==1091){  /* HadCM3 + A2  */
		if( (fp_t=fopen("./data/tas_HADCM3_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_HADCM3_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_HADCM3_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_HADCM3_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_HCA1B-2000-2200\n");  exit(1); }
	}else if(GCM==1092){  /* HadCM3 + B1  */
		if( (fp_t=fopen("./data/tas_HADCM3_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_tas_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_HADCM3_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_pr_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_HADCM3_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_huss_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_HADCM3_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_rsds_HCA1B-2000-2200\n");  exit(1); }
	}
	
	else if(GCM==1100){  /*  GISS AOM + A1B R1 */
		if( (fp_t=fopen("./data/tas_GISSAOM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_GISSAOM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GISSAOM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_GISSAOM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GISSAOM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_GISSAOM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GISSAOM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_GISSAOM_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1101){  /*  GISS AOM + A1B R2 */
		if( (fp_t=fopen("./data/tas_GISSAOM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No tas_GISSAOM_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GISSAOM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No pr_GISSAOM_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GISSAOM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No huss_GISSAOM_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GISSAOM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No rsds_GISSAOM_20C-A1B_R2.dat\n");  exit(1); }
	}else if(GCM==1102){  /*  GISS AOM + B1 R1 (corrected 2009/04/08) */
		if( (fp_t=fopen("./data/tas_GISSAOM_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_GISSAOM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GISSAOM_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_GISSAOM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GISSAOM_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_GISSAOM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GISSAOM_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_GISSAOM_20C-B1_R1.dat\n");  exit(1); }
	}else if(GCM==1103){  /*  GISS AOM + B1 R2 (corrected 2009/04/08) */
		if( (fp_t=fopen("./data/tas_GISSAOM_20C-B1_R2.dat","rt"))==NULL ){  printf("No tas_GISSAOM_20C-B2_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GISSAOM_20C-B1_R2.dat","rt"))==NULL ){  printf("No pr_GISSAOM_20C-B2_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GISSAOM_20C-B1_R2.dat","rt"))==NULL ){  printf("No huss_GISSAOM_20C-B2_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GISSAOM_20C-B1_R2.dat","rt"))==NULL ){  printf("No rsds_GISSAOM_20C-B2_R1.dat\n");  exit(1); }
	}
	
	if(GCM==1110){  /* CCCma-T63 + A1B  */
		if( (fp_t=fopen("./data/tas_CCC63_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_tas_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC63_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_pr_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC63_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_huss_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC63_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_rsds_CC63A1B-2001-2100\n");  exit(1); }
	}else if(GCM==1111){  /* CCCma-T63 + A2  */
		if( (fp_t=fopen("./data/tas_CCC63_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC63_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC63_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC63_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_CC63A1B-2001-2100\n");  exit(1); }
	}else if(GCM==1112){  /* CCCma-T63 + B1  */
		if( (fp_t=fopen("./data/tas_CCC63_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_tas_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC63_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_pr_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC63_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_huss_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC63_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_rsds_CC63A1B-2001-2100\n");  exit(1); }
	}
	
	if(GCM==1120){  /* CSIRO + A1B  */
		if( (fp_t=fopen("./data/tas_CSIRO_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CSIRO_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CSIRO_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CSIRO_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1121){  /* CSIRO + A2  */
		if( (fp_t=fopen("./data/tas_CSIRO_20C-A2_R1.dat","rt"))==NULL ){  printf("No tas_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CSIRO_20C-A2_R1.dat","rt"))==NULL ){  printf("No pr_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CSIRO_20C-A2_R1.dat","rt"))==NULL ){  printf("No huss_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CSIRO_20C-A2_R1.dat","rt"))==NULL ){  printf("No rsds_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1122){  /* CSIRO + B1  */
		if( (fp_t=fopen("./data/tas_CSIRO_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CSIRO_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CSIRO_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CSIRO_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
	}
	
	if(GCM==1130){  /*  MRI + A1B  */
		if( (fp_t=fopen("./data/tas_MRI_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_MRI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MRI_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_MRI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MRI_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_MRI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MRI_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_MRI_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1131){  /*  MRI + A1B  */
		if( (fp_t=fopen("./data/tas_MRI_20C-A1B_R2.dat","rt"))==NULL ){  printf("No tas_MRI_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MRI_20C-A1B_R2.dat","rt"))==NULL ){  printf("No pr_MRI_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MRI_20C-A1B_R2.dat","rt"))==NULL ){  printf("No huss_MRI_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MRI_20C-A1B_R2.dat","rt"))==NULL ){  printf("No rsds_MRI_20C-A1B_R2.dat\n");  exit(1); }
	}else if(GCM==1132){  /*  MRI + A1B  */
		if( (fp_t=fopen("./data/tas_MRI_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_tas_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MRI_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_pr_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MRI_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_huss_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MRI_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_rsds_MRIA1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1133){  /*  MRI + A1B  */
		if( (fp_t=fopen("./data/tas_MRI_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_tas_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MRI_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_pr_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MRI_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_huss_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MRI_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_rsds_MRIA1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1134){  /*  MRI + A1B  */
		if( (fp_t=fopen("./data/tas_MRI_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_tas_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MRI_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_pr_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MRI_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_huss_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MRI_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_rsds_MRIA1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1135){  /*  MRI + A2  */
		if( (fp_t=fopen("./data/tas_MRI_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MRI_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MRI_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MRI_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_MRIA2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1136){  /*  MRI + A2  */
		if( (fp_t=fopen("./data/tas_MRI_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_tas_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MRI_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_pr_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MRI_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_huss_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MRI_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_rsds_MRIA2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1137){  /*  MRI + A2  */
		if( (fp_t=fopen("./data/tas_MRI_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_tas_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MRI_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_pr_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MRI_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_huss_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MRI_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_rsds_MRIA2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1138){  /*  MRI + A2  */
		if( (fp_t=fopen("./data/tas_MRI_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_tas_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MRI_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_pr_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MRI_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_huss_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MRI_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_rsds_MRIA2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1139){  /*  MRI + A2  */
		if( (fp_t=fopen("./data/tas_MRI_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_tas_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MRI_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_pr_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MRI_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_huss_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MRI_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_rsds_MRIA2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1140){  /*  MRI + B1  */
		if( (fp_t=fopen("./data/tas_MRI_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_tas_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MRI_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_pr_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MRI_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_huss_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MRI_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_rsds_MRIB1-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1141){  /*  MRI + B1  */
		if( (fp_t=fopen("./data/tas_MRI_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_tas_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MRI_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_pr_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MRI_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_huss_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MRI_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_rsds_MRIB1-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1142){  /*  MRI + B1  */
		if( (fp_t=fopen("./data/tas_MRI_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_tas_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MRI_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_pr_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MRI_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_huss_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MRI_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_rsds_MRIB1-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1143){  /*  MRI + B1  */
		if( (fp_t=fopen("./data/tas_MRI_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_tas_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MRI_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_pr_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MRI_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_huss_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MRI_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_rsds_MRIB1-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1144){  /*  MRI + B1  */
		if( (fp_t=fopen("./data/tas_MRI_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_tas_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MRI_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_pr_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MRI_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_huss_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MRI_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_rsds_MRIB1-2001-2100.dat\n");  exit(1); }
	}
	
	if(GCM==1150){  /* 1201: MPI + A1B  */
		if( (fp_t=fopen("./data/tas_MPI_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MPI_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MPI_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MPI_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1151){  /* 1201: MPI + A1B  */
		if( (fp_t=fopen("./data/tas_MPI_20C-A1B_R2.dat","rt"))==NULL ){  printf("No tas_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MPI_20C-A1B_R2.dat","rt"))==NULL ){  printf("No pr_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MPI_20C-A1B_R2.dat","rt"))==NULL ){  printf("No huss_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MPI_20C-A1B_R2.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1152){  /* 1201: MPI + A1B  */
		if( (fp_t=fopen("./data/tas_MPI_20C-A1B_R3.dat","rt"))==NULL ){  printf("No tas_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MPI_20C-A1B_R3.dat","rt"))==NULL ){  printf("No pr_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MPI_20C-A1B_R3.dat","rt"))==NULL ){  printf("No huss_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MPI_20C-A1B_R3.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1153){  /* 1201: MPI + A1B  */
		if( (fp_t=fopen("./data/tas_MPI_20C-A1B_R4.dat","rt"))==NULL ){  printf("No tas_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MPI_20C-A1B_R4.dat","rt"))==NULL ){  printf("No pr_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MPI_20C-A1B_R4.dat","rt"))==NULL ){  printf("No huss_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MPI_20C-A1B_R4.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1154){  /* 1201: MPI + A2  */
		if( (fp_t=fopen("./data/tas_MPI_20C-A2_R1.dat","rt"))==NULL ){  printf("No tas_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MPI_20C-A2_R1.dat","rt"))==NULL ){  printf("No pr_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MPI_20C-A2_R1.dat","rt"))==NULL ){  printf("No huss_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MPI_20C-A2_R1.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-A2_R1.dat\n");  exit(1); }
	}else if(GCM==1155){  /* 1201: MPI + A2  */
		if( (fp_t=fopen("./data/tas_MPI_20C-A2_R2.dat","rt"))==NULL ){  printf("No tas_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MPI_20C-A2_R2.dat","rt"))==NULL ){  printf("No pr_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MPI_20C-A2_R2.dat","rt"))==NULL ){  printf("No huss_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MPI_20C-A2_R2.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-A2_R1.dat\n");  exit(1); }
	}else if(GCM==1156){  /* 1201: MPI + A2  */
		if( (fp_t=fopen("./data/tas_MPI_20C-A2_R3.dat","rt"))==NULL ){  printf("No tas_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MPI_20C-A2_R3.dat","rt"))==NULL ){  printf("No pr_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MPI_20C-A2_R3.dat","rt"))==NULL ){  printf("No huss_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MPI_20C-A2_R3.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-A2_R1.dat\n");  exit(1); }
	}else if(GCM==1157){  /* 1201: MPI + B1  */
		if( (fp_t=fopen("./data/tas_MPI_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MPI_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MPI_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MPI_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-B1_R1.dat\n");  exit(1); }
	}else if(GCM==1158){  /* 1201: MPI + B1  */
		if( (fp_t=fopen("./data/tas_MPI_20C-B1_R2.dat","rt"))==NULL ){  printf("No tas_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MPI_20C-B1_R2.dat","rt"))==NULL ){  printf("No pr_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MPI_20C-B1_R2.dat","rt"))==NULL ){  printf("No huss_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MPI_20C-B1_R2.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-B1_R1.dat\n");  exit(1); }
	}else if(GCM==1159){  /* 1201: MPI + B1  */
		if( (fp_t=fopen("./data/tas_MPI_20C-B1_R3.dat","rt"))==NULL ){  printf("No tas_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MPI_20C-B1_R3.dat","rt"))==NULL ){  printf("No pr_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MPI_20C-B1_R3.dat","rt"))==NULL ){  printf("No huss_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MPI_20C-B1_R3.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-B1_R1.dat\n");  exit(1); }
	}

	else if(GCM==1160){  /* IAP + A1B  */
		if( (fp_t=fopen("./data/tas_IAP_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_IAP_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_IAP_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_IAP_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1161){  /* IAP + A1B  */
		if( (fp_t=fopen("./data/tas_IAP_20C-A1B_R2.dat","rt"))==NULL ){  printf("No tas_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_IAP_20C-A1B_R2.dat","rt"))==NULL ){  printf("No pr_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_IAP_20C-A1B_R2.dat","rt"))==NULL ){  printf("No huss_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_IAP_20C-A1B_R2.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1162){  /* IAP + A1B  */
		if( (fp_t=fopen("./data/tas_IAP_20C-A1B_R3.dat","rt"))==NULL ){  printf("No tas_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_IAP_20C-A1B_R3.dat","rt"))==NULL ){  printf("No pr_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_IAP_20C-A1B_R3.dat","rt"))==NULL ){  printf("No huss_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_IAP_20C-A1B_R3.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1163){  /* IAP + B1  */
		if( (fp_t=fopen("./data/tas_IAP_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_IAP_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_IAP_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_IAP_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-B1_R1.dat\n");  exit(1); }
	}else if(GCM==1164){  /* IAP + B1  */
		if( (fp_t=fopen("./data/tas_IAP_20C-B1_R2.dat","rt"))==NULL ){  printf("No tas_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_IAP_20C-B1_R2.dat","rt"))==NULL ){  printf("No pr_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_IAP_20C-B1_R2.dat","rt"))==NULL ){  printf("No huss_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_IAP_20C-B1_R2.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-B1_R1.dat\n");  exit(1); }
	}else if(GCM==1165){  /* IAP + B1  */
		if( (fp_t=fopen("./data/tas_IAP_20C-B1_R3.dat","rt"))==NULL ){  printf("No tas_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_IAP_20C-B1_R3.dat","rt"))==NULL ){  printf("No pr_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_IAP_20C-B1_R3.dat","rt"))==NULL ){  printf("No huss_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_IAP_20C-B1_R3.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-B1_R1.dat\n");  exit(1); }
	}

	if(GCM==1170){  /* HADGEM + A1B  */
		if( (fp_t=fopen("./data/tas_HADGEM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_HADGEM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_HADGEM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_HADGEM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_HADGEM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_HADGEM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_HADGEM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1171){  /* HADGEM + A2  */
		if( (fp_t=fopen("./data/tas_HADGEM_20C-A2_R1.dat","rt"))==NULL ){  printf("No tas_HADGEM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_HADGEM_20C-A2_R1.dat","rt"))==NULL ){  printf("No pr_HADGEM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_HADGEM_20C-A2_R1.dat","rt"))==NULL ){  printf("No huss_HADGEM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_HADGEM_20C-A2_R1.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
	}
	
	if(GCM==1180){  /* GFDL20 + A1B  */
		if( (fp_t=fopen("./data/tas_GFDL20_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GFDL20_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GFDL20_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GFDL20_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1181){  /* GFDL20 + A2  */
		if( (fp_t=fopen("./data/tas_GFDL20_20C-A2_R1.dat","rt"))==NULL ){  printf("No tas_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GFDL20_20C-A2_R1.dat","rt"))==NULL ){  printf("No pr_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GFDL20_20C-A2_R1.dat","rt"))==NULL ){  printf("No huss_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GFDL20_20C-A2_R1.dat","rt"))==NULL ){  printf("No rsds_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1182){  /* GFDL20 + B1  */
		if( (fp_t=fopen("./data/tas_GFDL20_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GFDL20_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GFDL20_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GFDL20_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
	}
	
	if(GCM==1190){  /* CNRM + A1B  */
		if( (fp_t=fopen("./data/tas_CNRM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CNRM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CNRM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CNRM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_CNRM_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1191){  /* CNRM + A2  */
		if( (fp_t=fopen("./data/tas_CNRM_20C-A2_R1.dat","rt"))==NULL ){  printf("No tas_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CNRM_20C-A2_R1.dat","rt"))==NULL ){  printf("No pr_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CNRM_20C-A2_R1.dat","rt"))==NULL ){  printf("No huss_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CNRM_20C-A2_R1.dat","rt"))==NULL ){  printf("No rsds_CNRM_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1192){  /* CNRM + B1  */
		if( (fp_t=fopen("./data/tas_CNRM_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CNRM_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CNRM_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CNRM_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_CNRM_20C-A1B_R1.dat\n");  exit(1); }
	}
	
	if(GCM==1200){  /*  CCC47 + A1B  */
		if( (fp_t=fopen("./data/tas_CCC47_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC47_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC47_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC47_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1201){  /*  CCC47 + A1B  */
		if( (fp_t=fopen("./data/tas_CCC47_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC47_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC47_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC47_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1202){  /*  CCC47 + A1B  */
		if( (fp_t=fopen("./data/tas_CCC47_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC47_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC47_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC47_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1203){  /*  CCC47 + A1B  */
		if( (fp_t=fopen("./data/tas_CCC47_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC47_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC47_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC47_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1204){  /*  CCC47 + A1B  */
		if( (fp_t=fopen("./data/tas_CCC47_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC47_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC47_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC47_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1205){  /*  CCC47 + A2  */
		if( (fp_t=fopen("./data/tas_CCC47_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC47_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC47_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC47_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1206){  /*  CCC47 + A2  */
		if( (fp_t=fopen("./data/tas_CCC47_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC47_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC47_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC47_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1207){  /*  CCC47 + A2  */
		if( (fp_t=fopen("./data/tas_CCC47_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC47_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC47_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC47_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1208){  /*  CCC47 + A2  */
		if( (fp_t=fopen("./data/tas_CCC47_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC47_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC47_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC47_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1209){  /*  CCC47 + A2  */
		if( (fp_t=fopen("./data/tas_CCC47_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC47_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC47_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC47_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1210){  /*  CCC47 + B1  */
		if( (fp_t=fopen("./data/tas_CCC47_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_tas_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC47_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_pr_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC47_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_huss_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC47_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47B1-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1211){  /*  CCC47 + B1  */
		if( (fp_t=fopen("./data/tas_CCC47_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_tas_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC47_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_pr_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC47_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_huss_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC47_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47B1-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1212){  /*  CCC47 + B1  */
		if( (fp_t=fopen("./data/tas_CCC47_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_tas_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC47_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_pr_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC47_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_huss_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC47_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47B1-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1213){  /*  CCC47 + B1  */
		if( (fp_t=fopen("./data/tas_CCC47_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_tas_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC47_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_pr_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC47_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_huss_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC47_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47B1-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1214){  /*  CCC47 + B1  */
		if( (fp_t=fopen("./data/tas_CCC47_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_tas_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCC47_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_pr_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCC47_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_huss_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCC47_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47B1-2001-2100.dat\n");  exit(1); }
	}
	
	if(GCM==1220){  /*  CCSM + A1B  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_CCSM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_CCSM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_CCSM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_CCSM_20C-A1B_R1.dat\n");  exit(1); }
	}else if(GCM==1221){  /*  CCSM + A1B  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1222){  /*  CCSM + A1B  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1223){  /*  CCSM + A1B  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1224){  /*  CCSM + A1B  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-A1B_R6.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-A1B_R6.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-A1B_R6.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-A1B_R6.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1225){  /*  CCSM + A1B  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-A1B_R7.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-A1B_R7.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-A1B_R7.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-A1B_R7.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1226){  /*  CCSM + A1B  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-A1B_R9.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-A1B_R9.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-A1B_R9.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-A1B_R9.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1227){  /*  CCSM + A2  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1228){  /*  CCSM + A2  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1229){  /*  CCSM + A2  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1230){  /*  CCSM + A2  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1231){  /*  CCSM + A2  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1232){  /*  CCSM + B1  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_tas_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_pr_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_huss_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMB1-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1233){  /*  CCSM + B1  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_tas_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_pr_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_huss_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMB1-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1234){  /*  CCSM + B1  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_tas_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_pr_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_huss_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMB1-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1235){  /*  CCSM + B1  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_tas_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_pr_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_huss_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMB1-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1236){  /*  CCSM + B1  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_tas_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_pr_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_huss_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMB1-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1237){  /*  CCSM + B1  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-B1_R6.dat","rt"))==NULL ){  printf("No ext_tas_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-B1_R6.dat","rt"))==NULL ){  printf("No ext_pr_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-B1_R6.dat","rt"))==NULL ){  printf("No ext_huss_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-B1_R6.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMB1-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1238){  /*  CCSM + B1  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-B1_R7.dat","rt"))==NULL ){  printf("No ext_tas_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-B1_R7.dat","rt"))==NULL ){  printf("No ext_pr_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-B1_R7.dat","rt"))==NULL ){  printf("No ext_huss_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-B1_R7.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMB1-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1239){  /*  CCSM + B1  */
		if( (fp_t=fopen("./data/tas_CCSM_20C-B1_R9.dat","rt"))==NULL ){  printf("No ext_tas_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CCSM_20C-B1_R9.dat","rt"))==NULL ){  printf("No ext_pr_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CCSM_20C-B1_R9.dat","rt"))==NULL ){  printf("No ext_huss_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CCSM_20C-B1_R9.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMB1-2001-2100.dat\n");  exit(1); }
	}
	
	if(GCM==1240){  /*  GISSE + A1B  */
		if( (fp_t=fopen("./data/tas_GISSE_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_tas_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GISSE_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_pr_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GISSE_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_huss_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GISSE_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_rsds_GISSEA1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1241){  /*  GISSE + A1B  */
		if( (fp_t=fopen("./data/tas_GISSE_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_tas_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GISSE_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_pr_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GISSE_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_huss_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GISSE_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_rsds_GISSEA1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1242){  /*  GISSE + A1B  */
		if( (fp_t=fopen("./data/tas_GISSE_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_tas_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GISSE_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_pr_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GISSE_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_huss_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GISSE_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_rsds_GISSEA1B-2001-2100.dat\n");  exit(1); }
	}
	
	if(GCM==1250){  /*  GISSR + A1B  */
		/* changed (R1 -> R4) 2009/04/07 by A.Ito */
		if( (fp_t=fopen("./data/tas_GISSR_20C-A1B_R4.dat","rt"))==NULL ){  printf("No tas_GISSR_20C-A1B_R4.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GISSR_20C-A1B_R4.dat","rt"))==NULL ){  printf("No pr_GISSR_20C-A1B_R4.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GISSR_20C-A1B_R4.dat","rt"))==NULL ){  printf("No huss_GISSR_20C-A1B_R4.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GISSR_20C-A1B_R4.dat","rt"))==NULL ){  printf("No rsds_GISSR_20C-A1B_R4.dat\n");  exit(1); }
	}else if(GCM==1251){  /*  GISSR + A1B  */
		if( (fp_t=fopen("./data/tas_GISSR_20C-A1B_R2.dat","rt"))==NULL ){  printf("No tas_GISSR_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GISSR_20C-A1B_R2.dat","rt"))==NULL ){  printf("No pr_GISSR_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GISSR_20C-A1B_R2.dat","rt"))==NULL ){  printf("No huss_GISSR_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GISSR_20C-A1B_R2.dat","rt"))==NULL ){  printf("No rsds_GISSR_20C-A1B_R2.dat\n");  exit(1); }
	}else if(GCM==1252){  /*  GISSR + A2  */
		if( (fp_t=fopen("./data/tas_GISSR_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_GISSRA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GISSR_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_GISSRA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GISSR_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_GISSRA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GISSR_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_GISSRA1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1253){  /*  GISSR + B1  */
		/* changed (R2 -> R1) 2009/04/07 by A.Ito */
		if( (fp_t=fopen("./data/tas_GISSR_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_GISSR_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GISSR_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_GISSR_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GISSR_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_GISSR_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GISSR_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_GISSR_20C-B1_R1.dat\n");  exit(1); }
	}
	
	if(GCM==1260){  /*  PCM + A1B  */
		if( (fp_t=fopen("./data/tas_PCM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_tas_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_PCM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_pr_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_PCM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_huss_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_PCM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_rsds_PCMA1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1261){  /*  PCM + A1B  */
		if( (fp_t=fopen("./data/tas_PCM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_tas_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_PCM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_pr_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_PCM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_huss_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_PCM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_rsds_PCMA1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1262){  /*  PCM + A1B  */
		if( (fp_t=fopen("./data/tas_PCM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_tas_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_PCM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_pr_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_PCM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_huss_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_PCM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_rsds_PCMA1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1263){  /*  PCM + A1B  */
		if( (fp_t=fopen("./data/tas_PCM_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_tas_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_PCM_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_pr_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_PCM_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_huss_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_PCM_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_rsds_PCMA1B-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1264){  /*  PCM + A2  */
		if( (fp_t=fopen("./data/tas_PCM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_PCM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_PCM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_PCM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_PCMA2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1265){  /*  PCM + A2  */
		if( (fp_t=fopen("./data/tas_PCM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_tas_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_PCM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_pr_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_PCM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_huss_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_PCM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_rsds_PCMA2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1266){  /*  PCM + A2  */
		if( (fp_t=fopen("./data/tas_PCM_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_tas_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_PCM_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_pr_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_PCM_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_huss_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_PCM_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_rsds_PCMA2-2001-2100.dat\n");  exit(1); }
	}else if(GCM==1267){  /*  PCM + A2  */
		if( (fp_t=fopen("./data/tas_PCM_20C-A2_R4.dat","rt"))==NULL ){  printf("No tas_PCM_20C-A2_R4.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_PCM_20C-A2_R4.dat","rt"))==NULL ){  printf("No pr_PCM_20C-A2_R4.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_PCM_20C-A2_R4.dat","rt"))==NULL ){  printf("No huss_PCM_20C-A2_R4.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_PCM_20C-A2_R4.dat","rt"))==NULL ){  printf("No rsds_PCM_20C-A2_R4.dat\n");  exit(1); }
	}
	
	if(GCM==1270){  /*  CSIRO 3.5 + A1  */
		if( (fp_t=fopen("./data/tas_CSIRO35_20C-A1_R1.dat","rt"))==NULL ){  printf("No tas_CSIRO35_20C-A1_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CSIRO35_20C-A1_R1.dat","rt"))==NULL ){  printf("No pr_CSIRO35_20C-A1_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CSIRO35_20C-A1_R1.dat","rt"))==NULL ){  printf("No huss_CSIRO35_20C-A1_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CSIRO35_20C-A1_R1.dat","rt"))==NULL ){  printf("No rsds_CSIRO35_20C-A1_R1.dat\n");  exit(1); }
	}else if(GCM==1271){  /*  CSIRO 3.5 + A2  */
		if( (fp_t=fopen("./data/tas_CSIRO35_20C-A2_R1.dat","rt"))==NULL ){  printf("No tas_CSIRO35_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CSIRO35_20C-A2_R1.dat","rt"))==NULL ){  printf("No pr_CSIRO35_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CSIRO35_20C-A2_R1.dat","rt"))==NULL ){  printf("No huss_CSIRO35_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CSIRO35_20C-A2_R1.dat","rt"))==NULL ){  printf("No rsds_CSIRO35_20C-A2_R1.dat\n");  exit(1); }
	}else if(GCM==1272){  /*  CSIRO 3.5 + B1  */
		if( (fp_t=fopen("./data/tas_CSIRO35_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_CSIRO35_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CSIRO35_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_CSIRO35_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CSIRO35_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_CSIRO35_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CSIRO35_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_CSIRO35_20C-B1_R1.dat\n");  exit(1); }
	}
    
    if(GCM==3000){  /*  GEOMIP:BNU-ESM base (RCP4.5)  */
		if( (fp_t=fopen("./data/tas_BNU-ESM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_BNU-ESM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_BNU-ESM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_BNU-ESM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_BNU-ESM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No huss_BNU-ESM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_BNU-ESM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_BNU-ESM_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(GCM==3003){  /*  GEOMIP: BNU-ESM G3  */
		if( (fp_t=fopen("./data/tas_BNU-ESM_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No tas_BNU-ESM_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_BNU-ESM_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No pr_BNU-ESM_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_BNU-ESM_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No huss_BNU-ESM_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_BNU-ESM_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No rsds_BNU-ESM_197001-210012_rcp45_g3.txt\n");  exit(1); }
	}else if(GCM==3004){  /*  GEOMIP: BNU-ESM G4  */
		if( (fp_t=fopen("./data/tas_BNU-ESM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No tas_BNU-ESM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_BNU-ESM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No pr_BNU-ESM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_BNU-ESM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No huss_BNU-ESM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_BNU-ESM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No rsds_BNU-ESM_197001-210012_rcp45_g4.txt\n");  exit(1); }
	}else if(GCM==3200){  /*  GEOMIP: GISS-E2-R base (RCP4.5)  */
		if( (fp_t=fopen("./data/tas_GISS-E2-R_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_GISS-E2-R_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GISS-E2-R_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_GISS-E2-R_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GISS-E2-R_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No huss_GISS-E2-R_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GISS-E2-R_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_GISS-E2-R_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(GCM==3203){  /*  GEOMIP: GISS-E2-R G3  */
		if( (fp_t=fopen("./data/tas_GISS-E2-R_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No tas_GISS-E2-R_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GISS-E2-R_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No pr_GISS-E2-R_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GISS-E2-R_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No huss_GISS-E2-R_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GISS-E2-R_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No rsds_GISS-E2-R_197001-210012_rcp45_g3.txt\n");  exit(1); }
	}else if(GCM==3204){  /*  GEOMIP: GISS-E2-R G4  */
		if( (fp_t=fopen("./data/tas_GISS-E2-R_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No tas_GISS-E2-R_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_GISS-E2-R_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No pr_GISS-E2-R_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_GISS-E2-R_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No huss_GISS-E2-R_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_GISS-E2-R_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No rsds_GISS-E2-R_197001-210012_rcp45_g4.txt\n");  exit(1); }
	}else if(GCM==3300){  /*  GEOMIP: HadGEM2-ES base (RCP4.5) */
		if( (fp_t=fopen("./data/tas_HadGEM2-ES_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_HadGEM2-ES_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_HadGEM2-ES_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_HadGEM2-ES_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_HadGEM2-ES_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No huss_HadGEM2-ES_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_HadGEM2-ES_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_HadGEM2-ES_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(GCM==3303){  /*  GEOMIP: HadGEM2-ES G3 */
		if( (fp_t=fopen("./data/tas_HadGEM2-ES_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No tas_HadGEM2-ES_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_HadGEM2-ES_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No pr_HadGEM2-ES_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_HadGEM2-ES_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No huss_HadGEM2-ES_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_HadGEM2-ES_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No rsds_HadGEM2-ES_197001-210012_rcp45_g3.txt\n");  exit(1); }
	}else if(GCM==3304){  /*  GEOMIP: HadGEM2-ES G4 */
		if( (fp_t=fopen("./data/tas_HadGEM2-ES_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No tas_HadGEM2-ES_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_HadGEM2-ES_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No pr_HadGEM2-ES_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_HadGEM2-ES_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No huss_HadGEM2-ES_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_HadGEM2-ES_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No rsds_HadGEM2-ES_197001-210012_rcp45_g4.txt\n");  exit(1); }
	}else if(GCM==3313){  /*  GEOMIP: HadGEM2-ES G3S */
		if( (fp_t=fopen("./data/tas_HadGEM2-ES_197001-210012_rcp45_g3s.txt","rt"))==NULL ){  printf("No tas_HadGEM2-ES_197001-210012_rcp45_g3s.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_HadGEM2-ES_197001-210012_rcp45_g3s.txt","rt"))==NULL ){  printf("No pr_HadGEM2-ES_197001-210012_rcp45_g3s.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_HadGEM2-ES_197001-210012_rcp45_g3s.txt","rt"))==NULL ){  printf("No huss_HadGEM2-ES_197001-210012_rcp45_g3s.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_HadGEM2-ES_197001-210012_rcp45_g3s.txt","rt"))==NULL ){  printf("No rsds_HadGEM2-ES_197001-210012_rcp45_g3s.txt\n");  exit(1); }
	}else if(GCM==3400){  /*  GEOMIP: IPSL-CM5A-LR base (RCP4.5) */
		if( (fp_t=fopen("./data/tas_IPSL-CM5A-LR_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_IPSL-CM5A-LR_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_IPSL-CM5A-LR_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_IPSL-CM5A-LR_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_IPSL-CM5A-LR_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No huss_IPSL-CM5A-LR_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_IPSL-CM5A-LR_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_IPSL-CM5A-LR_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(GCM==3403){  /*  GEOMIP: IPSL-CM5A-LR G3 */
		if( (fp_t=fopen("./data/tas_IPSL-CM5A-LR_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No tas_IPSL-CM5A-LR_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_IPSL-CM5A-LR_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No pr_IPSL-CM5A-LR_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_IPSL-CM5A-LR_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No huss_IPSL-CM5A-LR_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_IPSL-CM5A-LR_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No rsds_IPSL-CM5A-LR_197001-210012_rcp45_g3.txt\n");  exit(1); }
	}else if(GCM==3405){  /*  GEOMIP: IPSL-CM5A-LR G5 */
		if( (fp_t=fopen("./data/tas_IPSL-CM5A-LR_197001-210012_rcp45_g5.txt","rt"))==NULL ){  printf("No tas_IPSL-CM5A-LR_197001-210012_rcp45_g5.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_IPSL-CM5A-LR_197001-210012_rcp45_g5.txt","rt"))==NULL ){  printf("No pr_IPSL-CM5A-LR_197001-210012_rcp45_g5.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_IPSL-CM5A-LR_197001-210012_rcp45_g5.txt","rt"))==NULL ){  printf("No huss_IPSL-CM5A-LR_197001-210012_rcp45_g5.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_IPSL-CM5A-LR_197001-210012_rcp45_g5.txt","rt"))==NULL ){  printf("No rsds_IPSL-CM5A-LR_197001-210012_rcp45_g5.txt\n");  exit(1); }
	}else if(GCM==3500){  /*  GEOMIP: MIROC-ESM base (RCP4.5)  */
		if( (fp_t=fopen("./data/tas_MIROC-ESM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_MIROC-ESM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MIROC-ESM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_MIROC-ESM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MIROC-ESM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No huss_MIROC-ESM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MIROC-ESM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_MIROC-ESM_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(GCM==3504){  /*  GEOMIP: MIROC-ESM G4  */
		if( (fp_t=fopen("./data/tas_MIROC-ESM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No tas_MIROC-ESM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MIROC-ESM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No pr_MIROC-ESM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MIROC-ESM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No huss_MIROC-ESM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MIROC-ESM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No rsds_MIROC-ESM_197001-210012_rcp45_g4.txt\n");  exit(1); }
	}else if(GCM==3600){  /*  GEOMIP: MIROC-ESM-CHEM base (RCP4.5)  */
		if( (fp_t=fopen("./data/tas_MIROC-ESM-CHEM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_MIROC-ESM-CHEM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MIROC-ESM-CHEM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_MIROC-ESM-CHEM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MIROC-ESM-CHEM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No huss_MIROC-ESM-CHEM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MIROC-ESM-CHEM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_MIROC-ESM-CHEM_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(GCM==3604){  /*  GEOMIP: MIROC-ESM-CHEM G4  */
		if( (fp_t=fopen("./data/tas_MIROC-ESM-CHEM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No tas_MIROC-ESM-CHEM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MIROC-ESM-CHEM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No pr_MIROC-ESM-CHEM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_MIROC-ESM-CHEM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No huss_MIROC-ESM-CHEM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MIROC-ESM-CHEM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No rsds_MIROC-ESM-CHEM_197001-210012_rcp45_g4.txt\n");  exit(1); }
	}else if(GCM==3700){  /*  GEOMIP: CanESM2 base (RCP4.5)  */
		if( (fp_t=fopen("./data/tas_CanESM2_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_CanESM2_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CanESM2_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_CanESM2_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CanESM2_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No huss_CanESM2_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CanESM2_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_CanESM2_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(GCM==3704){  /*  GEOMIP: CanESM2 G4  */
		if( (fp_t=fopen("./data/tas_CanESM2_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No tas_CanESM2_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_CanESM2_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No pr_CanESM2_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/huss_CanESM2_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No huss_CanESM2_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_CanESM2_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No rsds_CanESM2_197001-210012_rcp45_g4.txt\n");  exit(1); }
	}else if(GCM==3800){  /*  GEOMIP: MPI-ESM-LR base (RCP4.5)  */
		if( (fp_t=fopen("./data/tas_MPI-ESM-LR_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_MPI-ESM-LR_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MPI-ESM-LR_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_MPI-ESM-LR_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/hus_MPI-ESM-LR_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No hus_MPI-ESM-LR_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MPI-ESM-LR_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_MPI-ESM-LR_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(GCM==3803){  /*  GEOMIP: MPI-ESM-LR G3  */
		if( (fp_t=fopen("./data/tas_MPI-ESM-LR_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No tas_MPI-ESM-LR_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_p=fopen("./data/pr_MPI-ESM-LR_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No pr_MPI-ESM-LR_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_h=fopen("./data/hus_MPI-ESM-LR_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No hus_MPI-ESM-LR_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_r=fopen("./data/rsds_MPI-ESM-LR_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No rsds_MPI-ESM-LR_197001-210012_rcp45_g3.txt\n");  exit(1); }
	}
	
	if(GCM!=0){
        /* altitude */
		alt = (grid->topo>=0.0)?grid->topo:0.0;
        
		for(f=0;f<DL_GCM;f++){   /*  131 => 1970-2100 */
			for(g=0;g<ASTEP;g++){
				fscanf(fp_t,"%ld %ld", &yr, &mon);
				fscanf(fp_p,"%ld %ld", &yr, &mon);
				fscanf(fp_h,"%ld %ld", &yr, &mon);
				fscanf(fp_r,"%ld %ld", &yr, &mon);
				for(h=0;h<GCM_R;h++){
					for(i=0;i<GCM_C;i++){
						fscanf(fp_t,"%f", &(grid->proj_tmp2m[f][g][h][i]));
						fscanf(fp_p,"%f", &(grid->proj_prec[f][g][h][i]));
						fscanf(fp_h,"%f", &(grid->proj_hum[f][g][h][i]));
						fscanf(fp_r,"%f", &(grid->proj_rad[f][g][h][i]));
					}
				}
				for(h=0;h<GCM_R;h++){
					for(i=0;i<GCM_C;i++){
						if(grid->proj_tmp2m[f][g][h][i] > 1000.0 || grid->proj_tmp2m[f][g][h][i] < 200.0){
							grid->proj_tmp2m[f][g][h][i] = (grid->proj_tmp2m[f][g][h][i-1] + grid->proj_tmp2m[f][g][h][i+1])/2.0;
						}
						
						if(grid->proj_prec[f][g][h][i] > 10.0 || grid->proj_prec[f][g][h][i] < 0.0){
							grid->proj_prec[f][g][h][i] = (grid->proj_prec[f][g][h][i-1] + grid->proj_prec[f][g][h][i+1])/2.0;
						}
						
						if(grid->proj_hum[f][g][h][i] > 10.0 || grid->proj_hum[f][g][h][i] < 0.0){
							grid->proj_hum[f][g][h][i] = (grid->proj_hum[f][g][h][i-1] + grid->proj_hum[f][g][h][i+1])/2.0;
						}
						
						if(grid->proj_rad[f][g][h][i] > 1000.0 || grid->proj_rad[f][g][h][i] < 0.0){
							grid->proj_rad[f][g][h][i] = (grid->proj_rad[f][g][h][i-1] + grid->proj_rad[f][g][h][i+1])/2.0;
						}
						
                        /* K => deg C */
						grid->proj_tmp2m[f][g][h][i] -= ZAT;
                        
                        /* kg m-2 s-1 => mm month–1 */
						grid->proj_prec[f][g][h][i] *= (float)MDN[g];
						if(GCM >= 1000){
							grid->proj_prec[f][g][h][i] *= 3600.0*24.0;
						}
						
						/* specific humidity (kg kg-1) to vapor pressure (hPa) */
						/* revided by A.Ito (2009/08/17) */
						atmp = grid->proj_tmp2m[f][g][h][i];
						apres = 1013.25*exp(-1.0*(28.964*0.001)*9.8*alt/(8.3144*(atmp+ZAT))); 
						shum = grid->proj_hum[f][g][h][i];
						grid->proj_hum[f][g][h][i] = apres * grid->proj_hum[f][g][h][i]/(0.622 + 0.378*grid->proj_hum[f][g][h][i]);
					}
				}
			}
		}
	}else{
		for(f=0;f<DL_GCM;f++){
			for(g=0;g<ASTEP;g++){
				for(h=0;h<GCM_R;h++){
					for(i=0;i<GCM_C;i++){
						grid->proj_tmp2m[f][g][h][i] = 0.0;
						grid->proj_prec[f][g][h][i] = 0.0;
						grid->proj_hum[f][g][h][i] = 0.0;
						grid->proj_rad[f][g][h][i] = 0.0;
					}
				}
			}
		}
	}
	
	/**********************************************************/
	for(f=0;f<30;f++){ /* 1970-1999 */
		for(g=0;g<ASTEP;g++){
			for(h=0;h<GCM_R;h++){
				for(i=0;i<GCM_C;i++){
					if(f==0){
						grid->proj_tmp2m_b[g][h][i] = 0.0;
						grid->proj_prec_b[g][h][i] = 0.0;
						grid->proj_hum_b[g][h][i] = 0.0;
						grid->proj_rad_b[g][h][i] = 0.0;
					}
				
					grid->proj_tmp2m_b[g][h][i] += grid->proj_tmp2m[f][g][h][i]/30.0;
					grid->proj_prec_b[g][h][i] += grid->proj_prec[f][g][h][i]/30.0;
					grid->proj_hum_b[g][h][i] += grid->proj_hum[f][g][h][i]/30.0;
					grid->proj_rad_b[g][h][i] += grid->proj_rad[f][g][h][i]/30.0;
				}
			}
		}
	}
	
	if(GCM!=0 && GCM!=106){
		fclose(fp_t);
		fclose(fp_p);
		fclose(fp_h);
		fclose(fp_r);
	}
}
