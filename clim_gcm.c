/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed  in CGER/NIES & RIGC/JAMSTEC						*/
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
    FILE *fp_c2[4],
	struct Grid  *grid
){
	long f, g, h, i;
	long yr, mon;
	double atmp, apres, shum, alt;
	extern double MDN[12];
	
	if(SCENARIO_ID == 1){
		if( (fp_c2[0]=fopen("CN_A1_T.dat","rt"))==NULL ){  printf("No CN_A1_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CN_A1_P.dat","rt"))==NULL ){  printf("No CN_A1_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CN_A1_H.dat","rt"))==NULL ){  printf("No CN_A1_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CN_A1_R.dat","rt"))==NULL ){  printf("No CN_A1_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID == 2){
		if( (fp_c2[0]=fopen("CN_A1FI_T.dat","rt"))==NULL ){  printf("No CN_A1FI_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CN_A1FI_P.dat","rt"))==NULL ){  printf("No CN_A1FI_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CN_A1FI_H.dat","rt"))==NULL ){  printf("No CN_A1FI_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CN_A1FI_R.dat","rt"))==NULL ){  printf("No CN_A1FI_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==3){
		if( (fp_c2[0]=fopen("CN_A1T_T.dat","rt"))==NULL ){  printf("No CN_A1T_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CN_A1T_P.dat","rt"))==NULL ){  printf("No CN_A1T_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CN_A1T_H.dat","rt"))==NULL ){  printf("No CN_A1T_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CN_A1T_R.dat","rt"))==NULL ){  printf("No CN_A1T_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==4){
		if( (fp_c2[0]=fopen("CN_A2_T.dat","rt"))==NULL ){  printf("No CN_A2_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CN_A2_P.dat","rt"))==NULL ){  printf("No CN_A2_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CN_A2_H.dat","rt"))==NULL ){  printf("No CN_A2_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CN_A2_R.dat","rt"))==NULL ){  printf("No CN_A2_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==5){
		if( (fp_c2[0]=fopen("CN_B1_T.dat","rt"))==NULL ){  printf("No CN_B1_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CN_B1_P.dat","rt"))==NULL ){  printf("No CN_B1_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CN_B1_H.dat","rt"))==NULL ){  printf("No CN_B1_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CN_B1_R.dat","rt"))==NULL ){  printf("No CN_B1_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==6){
		if( (fp_c2[0]=fopen("CN_B2_T.dat","rt"))==NULL ){  printf("No CN_B2_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CN_B2_P.dat","rt"))==NULL ){  printf("No CN_B2_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CN_B2_H.dat","rt"))==NULL ){  printf("No CN_B2_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CN_B2_R.dat","rt"))==NULL ){  printf("No CN_B2_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==11){
		if( (fp_c2[0]=fopen("CC_A2_T.dat","rt"))==NULL ){  printf("No CC_A2_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CC_A2_P.dat","rt"))==NULL ){  printf("No CC_A2_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CC_A2_H.dat","rt"))==NULL ){  printf("No CC_A2_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CC_A2_R.dat","rt"))==NULL ){  printf("No CC_A2_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==12){
		if( (fp_c2[0]=fopen("CC_B2_T.dat","rt"))==NULL ){  printf("No CC_B2_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CC_B2_P.dat","rt"))==NULL ){  printf("No CC_B2_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CC_B2_H.dat","rt"))==NULL ){  printf("No CC_B2_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CC_B2_R.dat","rt"))==NULL ){  printf("No CC_B2_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==13){
		if( (fp_c2[0]=fopen("CC_A2a_T.dat","rt"))==NULL ){  printf("No CC_A2a_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CC_A2a_P.dat","rt"))==NULL ){  printf("No CC_A2a_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CC_A2a_H.dat","rt"))==NULL ){  printf("No CC_A2a_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CC_A2a_R.dat","rt"))==NULL ){  printf("No CC_A2a_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==14){
		if( (fp_c2[0]=fopen("CC_A2b_T.dat","rt"))==NULL ){  printf("No CC_A2b_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CC_A2b_P.dat","rt"))==NULL ){  printf("No CC_A2b_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CC_A2b_H.dat","rt"))==NULL ){  printf("No CC_A2b_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CC_A2b_R.dat","rt"))==NULL ){  printf("No CC_A2b_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==15){
		if( (fp_c2[0]=fopen("CC_A2c_T.dat","rt"))==NULL ){  printf("No CC_A2c_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CC_A2c_P.dat","rt"))==NULL ){  printf("No CC_A2c_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CC_A2c_H.dat","rt"))==NULL ){  printf("No CC_A2c_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CC_A2c_R.dat","rt"))==NULL ){  printf("No CC_A2c_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==16){
		if( (fp_c2[0]=fopen("CC_B2a_T.dat","rt"))==NULL ){  printf("No CC_B2a_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CC_B2a_P.dat","rt"))==NULL ){  printf("No CC_B2a_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CC_B2a_H.dat","rt"))==NULL ){  printf("No CC_B2a_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CC_B2a_R.dat","rt"))==NULL ){  printf("No CC_B2a_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==17){
		if( (fp_c2[0]=fopen("CC_B2b_T.dat","rt"))==NULL ){  printf("No CC_B2b_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CC_B2b_P.dat","rt"))==NULL ){  printf("No CC_B2b_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CC_B2b_H.dat","rt"))==NULL ){  printf("No CC_B2b_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CC_B2b_R.dat","rt"))==NULL ){  printf("No CC_B2b_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==18){
		if( (fp_c2[0]=fopen("CC_B2c_T.dat","rt"))==NULL ){  printf("No CC_B2c_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CC_B2c_P.dat","rt"))==NULL ){  printf("No CC_B2c_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CC_B2c_H.dat","rt"))==NULL ){  printf("No CC_B2c_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CC_B2c_R.dat","rt"))==NULL ){  printf("No CC_B2c_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==21){
		if( (fp_c2[0]=fopen("HD_A2_T.dat","rt"))==NULL ){  printf("No HD_A2_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("HD_A2_P.dat","rt"))==NULL ){  printf("No HD_A2_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("HD_A2_H.dat","rt"))==NULL ){  printf("No HD_A2_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("HD_A2_R.dat","rt"))==NULL ){  printf("No HD_A2_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==22){
		if( (fp_c2[0]=fopen("HD_A2b_T.dat","rt"))==NULL ){  printf("No HD_A2b_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("HD_A2b_P.dat","rt"))==NULL ){  printf("No HD_A2b_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("HD_A2b_H.dat","rt"))==NULL ){  printf("No HD_A2b_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("HD_A2b_R.dat","rt"))==NULL ){  printf("No HD_A2b_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==23){
		if( (fp_c2[0]=fopen("HD_A2c_T.dat","rt"))==NULL ){  printf("No HD_A2c_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("HD_A2c_P.dat","rt"))==NULL ){  printf("No HD_A2c_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("HD_A2c_H.dat","rt"))==NULL ){  printf("No HD_A2c_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("HD_A2c_R.dat","rt"))==NULL ){  printf("No HD_A2c_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==24){
		if( (fp_c2[0]=fopen("HD_B2_T.dat","rt"))==NULL ){  printf("No HD_B2_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("HD_B2_P.dat","rt"))==NULL ){  printf("No HD_B2_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("HD_B2_H.dat","rt"))==NULL ){  printf("No HD_B2_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("HD_B2_R.dat","rt"))==NULL ){  printf("No HD_B2_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==31){
		if( (fp_c2[0]=fopen("EH_A2_T.dat","rt"))==NULL ){  printf("No EH_A2_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("EH_A2_P.dat","rt"))==NULL ){  printf("No EH_A2_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("EH_A2_H.dat","rt"))==NULL ){  printf("No EH_A2_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("EH_A2_R.dat","rt"))==NULL ){  printf("No EH_A2_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==32){
		if( (fp_c2[0]=fopen("EH_B2_T.dat","rt"))==NULL ){  printf("No EH_B2_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("EH_B2_P.dat","rt"))==NULL ){  printf("No EH_B2_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("EH_B2_H.dat","rt"))==NULL ){  printf("No EH_B2_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("EH_B2_R.dat","rt"))==NULL ){  printf("No EH_B2_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==41){
		if( (fp_c2[0]=fopen("CS_A1_T.dat","rt"))==NULL ){  printf("No CS_A1_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CS_A1_P.dat","rt"))==NULL ){  printf("No CS_A1_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CS_A1_H.dat","rt"))==NULL ){  printf("No CS_A1_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CS_A1_R.dat","rt"))==NULL ){  printf("No CS_A1_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==42){
		if( (fp_c2[0]=fopen("CS_A2_T.dat","rt"))==NULL ){  printf("No CS_A2_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CS_A2_P.dat","rt"))==NULL ){  printf("No CS_A2_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CS_A2_H.dat","rt"))==NULL ){  printf("No CS_A2_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CS_A2_R.dat","rt"))==NULL ){  printf("No CS_A2_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==43){
		if( (fp_c2[0]=fopen("CS_B1_T.dat","rt"))==NULL ){  printf("No CS_B1_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CS_B1_P.dat","rt"))==NULL ){  printf("No CS_B1_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CS_B1_H.dat","rt"))==NULL ){  printf("No CS_B1_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CS_B1_R.dat","rt"))==NULL ){  printf("No CS_B1_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==44){
		if( (fp_c2[0]=fopen("CS_B2_T.dat","rt"))==NULL ){  printf("No CS_B2_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("CS_B2_P.dat","rt"))==NULL ){  printf("No CS_B2_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("CS_B2_H.dat","rt"))==NULL ){  printf("No CS_B2_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("CS_B2_R.dat","rt"))==NULL ){  printf("No CS_B2_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==51){
		if( (fp_c2[0]=fopen("GF_A2_T.dat","rt"))==NULL ){  printf("No GF_A2_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("GF_A2_P.dat","rt"))==NULL ){  printf("No GF_A2_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("GF_A2_H.dat","rt"))==NULL ){  printf("No GF_A2_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("GF_A2_R.dat","rt"))==NULL ){  printf("No GF_A2_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==52){
		if( (fp_c2[0]=fopen("GF_B2_T.dat","rt"))==NULL ){  printf("No GF_B2_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("GF_B2_P.dat","rt"))==NULL ){  printf("No GF_B2_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("GF_B2_H.dat","rt"))==NULL ){  printf("No GF_B2_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("GF_B2_R.dat","rt"))==NULL ){  printf("No GF_B2_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==61){
		if( (fp_c2[0]=fopen("NP_A1B_T.dat","rt"))==NULL ){  printf("No NP_A1B_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("NP_A1B_P.dat","rt"))==NULL ){  printf("No NP_A1B_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("NP_A1B_H.dat","rt"))==NULL ){  printf("No NP_A1B_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("NP_A1B_R.dat","rt"))==NULL ){  printf("No NP_A1B_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==62){
		if( (fp_c2[0]=fopen("NP_A2_T.dat","rt"))==NULL ){  printf("No NP_A2_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("NP_A2_P.dat","rt"))==NULL ){  printf("No NP_A2_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("NP_A2_H.dat","rt"))==NULL ){  printf("No NP_A2_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("NP_A2_R.dat","rt"))==NULL ){  printf("No NP_A2_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==63){
		if( (fp_c2[0]=fopen("NP_B2_T.dat","rt"))==NULL ){  printf("No NP_B2_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("NP_B2_P.dat","rt"))==NULL ){  printf("No NP_B2_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("NP_B2_H.dat","rt"))==NULL ){  printf("No NP_B2_H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("NP_B2_R.dat","rt"))==NULL ){  printf("No NP_B2_R.dat\n");  exit(1); }
	}else if(SCENARIO_ID==71){
		if( (fp_c2[0]=fopen("NC_A2_T.dat","rt"))==NULL ){  printf("No NC_A2_T.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("NC_A2_P.dat","rt"))==NULL ){  printf("No NC_A2_P.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("NP_A2_H.dat","rt"))==NULL ){  printf("No NP_A2_H.dat\n");  exit(1); }  /*  alternative data */
		if( (fp_c2[3]=fopen("NC_A2_R.dat","rt"))==NULL ){  printf("No NC_A2_R.dat\n");  exit(1); }
	}
	
	if(SCENARIO_ID == 1000){  /* MIROC-HIGH + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_MIROCH_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MIROCH_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MIROCH_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MIROCH_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID == 1001){  /* MIROC-HIGH + B1  */
		if( (fp_c2[0]=fopen("./data/tas_MIROCH_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_tas_B1-2000-2099-H.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MIROCH_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_pr_B1-2000-2099-H.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MIROCH_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_huss_B1-2000-2099-H.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MIROCH_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_rsds_B1-2000-2099-H.dat\n");  exit(1); }
	}
	
	if(SCENARIO_ID == 1010){  /* MIROC-MED + A1B Run1 */
		if( (fp_c2[0]=fopen("./data/tas_MIROCM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_MIROCM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MIROCM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_MIROCM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MIROCM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_MIROCM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MIROCM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_MIROCM_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID == 1011){  /* MIROC-MED + A1B Run2 */
		if( (fp_c2[0]=fopen("./data/tas_MIROCM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No tas_MIROCM_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MIROCM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No pr_MIROCM_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MIROCM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No huss_MIROCM_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MIROCM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No rsds_MIROCM_20C-A1B_R2.dat\n");  exit(1); }
	}else if(SCENARIO_ID == 1012){  /* MIROC-MED + A1B Run3 */
		if( (fp_c2[0]=fopen("./data/tas_MIROCM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No tas_MIROCM_20C-A1B_R3.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MIROCM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No pr_MIROCM_20C-A1B_R3.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MIROCM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No huss_MIROCM_20C-A1B_R3.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MIROCM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No rsds_MIROCM_20C-A1B_R3.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1013){  /* MIROC-MED + A2  */
		if( (fp_c2[0]=fopen("./data/tas_MIROCM_20C-A2_R1.dat","rt"))==NULL ){  printf("No tas_MIROCM_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MIROCM_20C-A2_R1.dat","rt"))==NULL ){  printf("No pr_MIROCM_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MIROCM_20C-A2_R1.dat","rt"))==NULL ){  printf("No huss_MIROCM_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MIROCM_20C-A2_R1.dat","rt"))==NULL ){  printf("No rsds_MIROCM_20C-A2_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1014){  /* MIROC-MED + A2  */
		if( (fp_c2[0]=fopen("./data/tas_MIROCM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_tas_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MIROCM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_pr_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MIROCM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_huss_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MIROCM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_rsds_A1B-2000-2099-R2.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1015){  /* MIROC-MED + A2  */
		if( (fp_c2[0]=fopen("./data/tas_MIROCM_20C-A2_R3.dat","rt"))==NULL ){  printf("No tas_MIROCM_20C-A2_R3.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MIROCM_20C-A2_R3.dat","rt"))==NULL ){  printf("No pr_MIROCM_20C-A2_R3.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MIROCM_20C-A2_R3.dat","rt"))==NULL ){  printf("No huss_MIROCM_20C-A2_R3.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MIROCM_20C-A2_R3.dat","rt"))==NULL ){  printf("No rsds_MIROCM_20C-A2_R3.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1016){  /* MIROC-MED + B1  */
		if( (fp_c2[0]=fopen("./data/tas_MIROCM_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_MIROCM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MIROCM_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_MIROCM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MIROCM_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_MIROCM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MIROCM_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_MIROCM_20C-B1_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1017){  /* MIROC-MED + B1  */
		if( (fp_c2[0]=fopen("./data/tas_MIROCM_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_tas_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MIROCM_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_pr_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MIROCM_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_huss_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MIROCM_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_rsds_A1B-2000-2099-R2.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1018){  /* MIROC-MED + B1  */
		if( (fp_c2[0]=fopen("./data/tas_MIROCM_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_tas_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MIROCM_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_pr_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MIROCM_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_huss_A1B-2000-2099-R2.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MIROCM_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_rsds_A1B-2000-2099-R2.dat\n");  exit(1); }
	}
	
	if(SCENARIO_ID == 1050){  /* BCCR + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_BCCR_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_BCCR_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_BCCR_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_BCCR_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_BCCR_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_BCCR_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_BCCR_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_BCCR_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID == 1051){  /* BCCR + A2  */
		if( (fp_c2[0]=fopen("./data/tas_BCCR_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_BRA1B-2001-2100\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_BCCR_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_BRA1B-2001-2100\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_BCCR_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_BRA1B-2001-2100\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_BCCR_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_BRA1B-2001-2100\n");  exit(1); }
	}else if(SCENARIO_ID == 1052){  /* BCCR + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_BCCR_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_tas_BRA1B-2001-2100\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_BCCR_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_pr_BRA1B-2001-2100\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_BCCR_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_huss_BRA1B-2001-2100\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_BCCR_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_rsds_BRA1B-2001-2100\n");  exit(1); }
	}
	
	if(SCENARIO_ID == 1060){  /* INM-CM3 + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_INM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_INM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_INM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_INM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_INM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_INM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_INM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_INM_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1061){  /* INM-CM3 + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_INM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_INA1B-2001-2200\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_INM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_INA1B-2001-2200\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_INM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_INA1B-2001-2200\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_INM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_INA1B-2001-2200\n");  exit(1); }
	}else if(SCENARIO_ID==1062){  /* INM-CM3 + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_INM_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_INM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_INM_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_INM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_INM_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_INM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_INM_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_INM_20C-B1_R1.dat\n");  exit(1); }
	}
	
	if(SCENARIO_ID==1070){  /* IPSL + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_IPSL_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_IPSL_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_IPSL_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_IPSL_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_IPSL_20C-A2_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1071){  /* IPSL + A2  */
		if( (fp_c2[0]=fopen("./data/tas_IPSL_20C-A2_R1.dat","rt"))==NULL ){  printf("No tas_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_IPSL_20C-A2_R1.dat","rt"))==NULL ){  printf("No pr_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_IPSL_20C-A2_R1.dat","rt"))==NULL ){  printf("No huss_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_IPSL_20C-A2_R1.dat","rt"))==NULL ){  printf("No rsds_IPSL_20C-A2_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1072){  /* IPSL + B1  */
		if( (fp_c2[0]=fopen("./data/tas_IPSL_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_IPSL_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_IPSL_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_IPSL_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_IPSL_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_IPSL_20C-A2_R1.dat\n");  exit(1); }
	}
	
	if(SCENARIO_ID == 1080){  /* GFDL CM2-1 + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_GFDL21_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_tas_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GFDL21_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_pr_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GFDL21_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_hus_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GFDL21_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_rsds_GFA1B-2001-2100\n");  exit(1); }
	}else if(SCENARIO_ID == 1081){  /* GFDL CM2-1 + A2  */
		if( (fp_c2[0]=fopen("./data/tas_GFDL21_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GFDL21_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GFDL21_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_hus_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GFDL21_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_GFA1B-2001-2100\n");  exit(1); }
	}else if(SCENARIO_ID == 1082){  /* GFDL CM2-1 + B1  */
		if( (fp_c2[0]=fopen("./data/tas_GFDL21_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_tas_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GFDL21_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_pr_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GFDL21_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_hus_GFA1B-2001-2100\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GFDL21_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_rsds_GFA1B-2001-2100\n");  exit(1); }
	}
	
	if(SCENARIO_ID==1090){  /* HadCM3 + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_HADCM3_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_tas_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_HADCM3_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_pr_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_HADCM3_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_huss_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_HADCM3_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_rsds_HCA1B-2000-2200\n");  exit(1); }
	}else if(SCENARIO_ID==1091){  /* HadCM3 + A2  */
		if( (fp_c2[0]=fopen("./data/tas_HADCM3_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_HADCM3_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_HADCM3_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_HADCM3_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_HCA1B-2000-2200\n");  exit(1); }
	}else if(SCENARIO_ID==1092){  /* HadCM3 + B1  */
		if( (fp_c2[0]=fopen("./data/tas_HADCM3_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_tas_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_HADCM3_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_pr_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_HADCM3_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_huss_HCA1B-2000-2200\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_HADCM3_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_rsds_HCA1B-2000-2200\n");  exit(1); }
	}
	
	else if(SCENARIO_ID==1100){  /*  GISS AOM + A1B R1 */
		if( (fp_c2[0]=fopen("./data/tas_GISSAOM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_GISSAOM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GISSAOM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_GISSAOM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GISSAOM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_GISSAOM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GISSAOM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_GISSAOM_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1101){  /*  GISS AOM + A1B R2 */
		if( (fp_c2[0]=fopen("./data/tas_GISSAOM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No tas_GISSAOM_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GISSAOM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No pr_GISSAOM_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GISSAOM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No huss_GISSAOM_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GISSAOM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No rsds_GISSAOM_20C-A1B_R2.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1102){  /*  GISS AOM + B1 R1 (corrected 2009/04/08) */
		if( (fp_c2[0]=fopen("./data/tas_GISSAOM_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_GISSAOM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GISSAOM_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_GISSAOM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GISSAOM_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_GISSAOM_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GISSAOM_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_GISSAOM_20C-B1_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1103){  /*  GISS AOM + B1 R2 (corrected 2009/04/08) */
		if( (fp_c2[0]=fopen("./data/tas_GISSAOM_20C-B1_R2.dat","rt"))==NULL ){  printf("No tas_GISSAOM_20C-B2_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GISSAOM_20C-B1_R2.dat","rt"))==NULL ){  printf("No pr_GISSAOM_20C-B2_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GISSAOM_20C-B1_R2.dat","rt"))==NULL ){  printf("No huss_GISSAOM_20C-B2_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GISSAOM_20C-B1_R2.dat","rt"))==NULL ){  printf("No rsds_GISSAOM_20C-B2_R1.dat\n");  exit(1); }
	}
	
	if(SCENARIO_ID==1110){  /* CCCma-T63 + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_CCC63_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_tas_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC63_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_pr_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC63_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_huss_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC63_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_rsds_CC63A1B-2001-2100\n");  exit(1); }
	}else if(SCENARIO_ID==1111){  /* CCCma-T63 + A2  */
		if( (fp_c2[0]=fopen("./data/tas_CCC63_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC63_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC63_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC63_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_CC63A1B-2001-2100\n");  exit(1); }
	}else if(SCENARIO_ID==1112){  /* CCCma-T63 + B1  */
		if( (fp_c2[0]=fopen("./data/tas_CCC63_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_tas_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC63_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_pr_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC63_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_huss_CC63A1B-2001-2100\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC63_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_rsds_CC63A1B-2001-2100\n");  exit(1); }
	}
	
	if(SCENARIO_ID==1120){  /* CSIRO + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_CSIRO_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CSIRO_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CSIRO_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CSIRO_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1121){  /* CSIRO + A2  */
		if( (fp_c2[0]=fopen("./data/tas_CSIRO_20C-A2_R1.dat","rt"))==NULL ){  printf("No tas_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CSIRO_20C-A2_R1.dat","rt"))==NULL ){  printf("No pr_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CSIRO_20C-A2_R1.dat","rt"))==NULL ){  printf("No huss_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CSIRO_20C-A2_R1.dat","rt"))==NULL ){  printf("No rsds_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1122){  /* CSIRO + B1  */
		if( (fp_c2[0]=fopen("./data/tas_CSIRO_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CSIRO_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CSIRO_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CSIRO_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_CSIRO_20C-A1B_R1.dat\n");  exit(1); }
	}
	
	if(SCENARIO_ID==1130){  /*  MRI + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_MRI_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_MRI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MRI_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_MRI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MRI_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_MRI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MRI_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_MRI_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1131){  /*  MRI + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_MRI_20C-A1B_R2.dat","rt"))==NULL ){  printf("No tas_MRI_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MRI_20C-A1B_R2.dat","rt"))==NULL ){  printf("No pr_MRI_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MRI_20C-A1B_R2.dat","rt"))==NULL ){  printf("No huss_MRI_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MRI_20C-A1B_R2.dat","rt"))==NULL ){  printf("No rsds_MRI_20C-A1B_R2.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1132){  /*  MRI + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_MRI_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_tas_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MRI_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_pr_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MRI_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_huss_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MRI_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_rsds_MRIA1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1133){  /*  MRI + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_MRI_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_tas_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MRI_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_pr_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MRI_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_huss_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MRI_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_rsds_MRIA1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1134){  /*  MRI + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_MRI_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_tas_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MRI_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_pr_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MRI_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_huss_MRIA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MRI_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_rsds_MRIA1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1135){  /*  MRI + A2  */
		if( (fp_c2[0]=fopen("./data/tas_MRI_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MRI_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MRI_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MRI_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_MRIA2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1136){  /*  MRI + A2  */
		if( (fp_c2[0]=fopen("./data/tas_MRI_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_tas_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MRI_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_pr_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MRI_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_huss_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MRI_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_rsds_MRIA2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1137){  /*  MRI + A2  */
		if( (fp_c2[0]=fopen("./data/tas_MRI_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_tas_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MRI_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_pr_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MRI_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_huss_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MRI_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_rsds_MRIA2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1138){  /*  MRI + A2  */
		if( (fp_c2[0]=fopen("./data/tas_MRI_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_tas_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MRI_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_pr_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MRI_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_huss_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MRI_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_rsds_MRIA2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1139){  /*  MRI + A2  */
		if( (fp_c2[0]=fopen("./data/tas_MRI_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_tas_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MRI_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_pr_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MRI_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_huss_MRIA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MRI_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_rsds_MRIA2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1140){  /*  MRI + B1  */
		if( (fp_c2[0]=fopen("./data/tas_MRI_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_tas_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MRI_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_pr_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MRI_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_huss_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MRI_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_rsds_MRIB1-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1141){  /*  MRI + B1  */
		if( (fp_c2[0]=fopen("./data/tas_MRI_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_tas_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MRI_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_pr_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MRI_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_huss_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MRI_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_rsds_MRIB1-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1142){  /*  MRI + B1  */
		if( (fp_c2[0]=fopen("./data/tas_MRI_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_tas_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MRI_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_pr_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MRI_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_huss_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MRI_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_rsds_MRIB1-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1143){  /*  MRI + B1  */
		if( (fp_c2[0]=fopen("./data/tas_MRI_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_tas_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MRI_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_pr_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MRI_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_huss_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MRI_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_rsds_MRIB1-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1144){  /*  MRI + B1  */
		if( (fp_c2[0]=fopen("./data/tas_MRI_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_tas_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MRI_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_pr_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MRI_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_huss_MRIB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MRI_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_rsds_MRIB1-2001-2100.dat\n");  exit(1); }
	}
	
	if(SCENARIO_ID==1150){  /* 1201: MPI + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_MPI_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MPI_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MPI_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MPI_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1151){  /* 1201: MPI + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_MPI_20C-A1B_R2.dat","rt"))==NULL ){  printf("No tas_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MPI_20C-A1B_R2.dat","rt"))==NULL ){  printf("No pr_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MPI_20C-A1B_R2.dat","rt"))==NULL ){  printf("No huss_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MPI_20C-A1B_R2.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1152){  /* 1201: MPI + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_MPI_20C-A1B_R3.dat","rt"))==NULL ){  printf("No tas_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MPI_20C-A1B_R3.dat","rt"))==NULL ){  printf("No pr_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MPI_20C-A1B_R3.dat","rt"))==NULL ){  printf("No huss_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MPI_20C-A1B_R3.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1153){  /* 1201: MPI + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_MPI_20C-A1B_R4.dat","rt"))==NULL ){  printf("No tas_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MPI_20C-A1B_R4.dat","rt"))==NULL ){  printf("No pr_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MPI_20C-A1B_R4.dat","rt"))==NULL ){  printf("No huss_MPI_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MPI_20C-A1B_R4.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1154){  /* 1201: MPI + A2  */
		if( (fp_c2[0]=fopen("./data/tas_MPI_20C-A2_R1.dat","rt"))==NULL ){  printf("No tas_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MPI_20C-A2_R1.dat","rt"))==NULL ){  printf("No pr_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MPI_20C-A2_R1.dat","rt"))==NULL ){  printf("No huss_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MPI_20C-A2_R1.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-A2_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1155){  /* 1201: MPI + A2  */
		if( (fp_c2[0]=fopen("./data/tas_MPI_20C-A2_R2.dat","rt"))==NULL ){  printf("No tas_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MPI_20C-A2_R2.dat","rt"))==NULL ){  printf("No pr_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MPI_20C-A2_R2.dat","rt"))==NULL ){  printf("No huss_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MPI_20C-A2_R2.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-A2_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1156){  /* 1201: MPI + A2  */
		if( (fp_c2[0]=fopen("./data/tas_MPI_20C-A2_R3.dat","rt"))==NULL ){  printf("No tas_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MPI_20C-A2_R3.dat","rt"))==NULL ){  printf("No pr_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MPI_20C-A2_R3.dat","rt"))==NULL ){  printf("No huss_MPI_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MPI_20C-A2_R3.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-A2_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1157){  /* 1201: MPI + B1  */
		if( (fp_c2[0]=fopen("./data/tas_MPI_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MPI_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MPI_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MPI_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-B1_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1158){  /* 1201: MPI + B1  */
		if( (fp_c2[0]=fopen("./data/tas_MPI_20C-B1_R2.dat","rt"))==NULL ){  printf("No tas_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MPI_20C-B1_R2.dat","rt"))==NULL ){  printf("No pr_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MPI_20C-B1_R2.dat","rt"))==NULL ){  printf("No huss_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MPI_20C-B1_R2.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-B1_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1159){  /* 1201: MPI + B1  */
		if( (fp_c2[0]=fopen("./data/tas_MPI_20C-B1_R3.dat","rt"))==NULL ){  printf("No tas_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MPI_20C-B1_R3.dat","rt"))==NULL ){  printf("No pr_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MPI_20C-B1_R3.dat","rt"))==NULL ){  printf("No huss_MPI_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MPI_20C-B1_R3.dat","rt"))==NULL ){  printf("No rsds_MPI_20C-B1_R1.dat\n");  exit(1); }
	}

	else if(SCENARIO_ID==1160){  /* IAP + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_IAP_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_IAP_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_IAP_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_IAP_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1161){  /* IAP + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_IAP_20C-A1B_R2.dat","rt"))==NULL ){  printf("No tas_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_IAP_20C-A1B_R2.dat","rt"))==NULL ){  printf("No pr_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_IAP_20C-A1B_R2.dat","rt"))==NULL ){  printf("No huss_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_IAP_20C-A1B_R2.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1162){  /* IAP + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_IAP_20C-A1B_R3.dat","rt"))==NULL ){  printf("No tas_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_IAP_20C-A1B_R3.dat","rt"))==NULL ){  printf("No pr_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_IAP_20C-A1B_R3.dat","rt"))==NULL ){  printf("No huss_IAP_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_IAP_20C-A1B_R3.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1163){  /* IAP + B1  */
		if( (fp_c2[0]=fopen("./data/tas_IAP_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_IAP_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_IAP_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_IAP_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-B1_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1164){  /* IAP + B1  */
		if( (fp_c2[0]=fopen("./data/tas_IAP_20C-B1_R2.dat","rt"))==NULL ){  printf("No tas_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_IAP_20C-B1_R2.dat","rt"))==NULL ){  printf("No pr_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_IAP_20C-B1_R2.dat","rt"))==NULL ){  printf("No huss_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_IAP_20C-B1_R2.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-B1_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1165){  /* IAP + B1  */
		if( (fp_c2[0]=fopen("./data/tas_IAP_20C-B1_R3.dat","rt"))==NULL ){  printf("No tas_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_IAP_20C-B1_R3.dat","rt"))==NULL ){  printf("No pr_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_IAP_20C-B1_R3.dat","rt"))==NULL ){  printf("No huss_IAP_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_IAP_20C-B1_R3.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-B1_R1.dat\n");  exit(1); }
	}

	if(SCENARIO_ID==1170){  /* HADGEM + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_HADGEM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_HADGEM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_HADGEM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_HADGEM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_HADGEM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_HADGEM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_HADGEM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1171){  /* HADGEM + A2  */
		if( (fp_c2[0]=fopen("./data/tas_HADGEM_20C-A2_R1.dat","rt"))==NULL ){  printf("No tas_HADGEM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_HADGEM_20C-A2_R1.dat","rt"))==NULL ){  printf("No pr_HADGEM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_HADGEM_20C-A2_R1.dat","rt"))==NULL ){  printf("No huss_HADGEM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_HADGEM_20C-A2_R1.dat","rt"))==NULL ){  printf("No rsds_MIROCH_20C-A1B_R1.dat\n");  exit(1); }
	}
	
	if(SCENARIO_ID==1180){  /* GFDL20 + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_GFDL20_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GFDL20_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GFDL20_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GFDL20_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1181){  /* GFDL20 + A2  */
		if( (fp_c2[0]=fopen("./data/tas_GFDL20_20C-A2_R1.dat","rt"))==NULL ){  printf("No tas_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GFDL20_20C-A2_R1.dat","rt"))==NULL ){  printf("No pr_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GFDL20_20C-A2_R1.dat","rt"))==NULL ){  printf("No huss_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GFDL20_20C-A2_R1.dat","rt"))==NULL ){  printf("No rsds_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1182){  /* GFDL20 + B1  */
		if( (fp_c2[0]=fopen("./data/tas_GFDL20_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GFDL20_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GFDL20_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GFDL20_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_GFDL20_20C-A1B_R1.dat\n");  exit(1); }
	}
	
	if(SCENARIO_ID==1190){  /* CNRM + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_CNRM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CNRM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CNRM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CNRM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_CNRM_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1191){  /* CNRM + A2  */
		if( (fp_c2[0]=fopen("./data/tas_CNRM_20C-A2_R1.dat","rt"))==NULL ){  printf("No tas_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CNRM_20C-A2_R1.dat","rt"))==NULL ){  printf("No pr_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CNRM_20C-A2_R1.dat","rt"))==NULL ){  printf("No huss_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CNRM_20C-A2_R1.dat","rt"))==NULL ){  printf("No rsds_CNRM_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1192){  /* CNRM + B1  */
		if( (fp_c2[0]=fopen("./data/tas_CNRM_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CNRM_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CNRM_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_CNRM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CNRM_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_CNRM_20C-A1B_R1.dat\n");  exit(1); }
	}
	
	if(SCENARIO_ID==1200){  /*  CCC47 + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_CCC47_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC47_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC47_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC47_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1201){  /*  CCC47 + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_CCC47_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC47_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC47_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC47_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1202){  /*  CCC47 + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_CCC47_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC47_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC47_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC47_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1203){  /*  CCC47 + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_CCC47_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC47_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC47_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC47_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1204){  /*  CCC47 + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_CCC47_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC47_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC47_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC47_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1205){  /*  CCC47 + A2  */
		if( (fp_c2[0]=fopen("./data/tas_CCC47_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC47_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC47_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC47_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1206){  /*  CCC47 + A2  */
		if( (fp_c2[0]=fopen("./data/tas_CCC47_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC47_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC47_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC47_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1207){  /*  CCC47 + A2  */
		if( (fp_c2[0]=fopen("./data/tas_CCC47_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC47_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC47_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC47_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1208){  /*  CCC47 + A2  */
		if( (fp_c2[0]=fopen("./data/tas_CCC47_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC47_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC47_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC47_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1209){  /*  CCC47 + A2  */
		if( (fp_c2[0]=fopen("./data/tas_CCC47_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_tas_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC47_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_pr_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC47_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_huss_CCC47A2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC47_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47A2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1210){  /*  CCC47 + B1  */
		if( (fp_c2[0]=fopen("./data/tas_CCC47_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_tas_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC47_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_pr_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC47_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_huss_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC47_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47B1-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1211){  /*  CCC47 + B1  */
		if( (fp_c2[0]=fopen("./data/tas_CCC47_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_tas_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC47_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_pr_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC47_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_huss_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC47_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47B1-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1212){  /*  CCC47 + B1  */
		if( (fp_c2[0]=fopen("./data/tas_CCC47_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_tas_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC47_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_pr_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC47_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_huss_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC47_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47B1-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1213){  /*  CCC47 + B1  */
		if( (fp_c2[0]=fopen("./data/tas_CCC47_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_tas_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC47_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_pr_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC47_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_huss_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC47_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47B1-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1214){  /*  CCC47 + B1  */
		if( (fp_c2[0]=fopen("./data/tas_CCC47_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_tas_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCC47_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_pr_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCC47_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_huss_CCC47B1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCC47_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_rsds_CCC47B1-2001-2100.dat\n");  exit(1); }
	}
	
	if(SCENARIO_ID==1220){  /*  CCSM + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No tas_CCSM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No pr_CCSM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No huss_CCSM_20C-A1B_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No rsds_CCSM_20C-A1B_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1221){  /*  CCSM + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1222){  /*  CCSM + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1223){  /*  CCSM + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-A1B_R5.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1224){  /*  CCSM + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-A1B_R6.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-A1B_R6.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-A1B_R6.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-A1B_R6.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1225){  /*  CCSM + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-A1B_R7.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-A1B_R7.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-A1B_R7.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-A1B_R7.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1226){  /*  CCSM + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-A1B_R9.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-A1B_R9.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-A1B_R9.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-A1B_R9.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1227){  /*  CCSM + A2  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1228){  /*  CCSM + A2  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1229){  /*  CCSM + A2  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1230){  /*  CCSM + A2  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-A2_R4.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1231){  /*  CCSM + A2  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_tas_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_pr_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_huss_CCSMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-A2_R5.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMA2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1232){  /*  CCSM + B1  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_tas_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_pr_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_huss_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-B1_R1.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMB1-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1233){  /*  CCSM + B1  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_tas_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_pr_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_huss_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-B1_R2.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMB1-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1234){  /*  CCSM + B1  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_tas_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_pr_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_huss_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-B1_R3.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMB1-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1235){  /*  CCSM + B1  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_tas_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_pr_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_huss_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-B1_R4.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMB1-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1236){  /*  CCSM + B1  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_tas_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_pr_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_huss_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-B1_R5.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMB1-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1237){  /*  CCSM + B1  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-B1_R6.dat","rt"))==NULL ){  printf("No ext_tas_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-B1_R6.dat","rt"))==NULL ){  printf("No ext_pr_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-B1_R6.dat","rt"))==NULL ){  printf("No ext_huss_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-B1_R6.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMB1-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1238){  /*  CCSM + B1  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-B1_R7.dat","rt"))==NULL ){  printf("No ext_tas_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-B1_R7.dat","rt"))==NULL ){  printf("No ext_pr_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-B1_R7.dat","rt"))==NULL ){  printf("No ext_huss_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-B1_R7.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMB1-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1239){  /*  CCSM + B1  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM_20C-B1_R9.dat","rt"))==NULL ){  printf("No ext_tas_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM_20C-B1_R9.dat","rt"))==NULL ){  printf("No ext_pr_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM_20C-B1_R9.dat","rt"))==NULL ){  printf("No ext_huss_CCSMB1-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM_20C-B1_R9.dat","rt"))==NULL ){  printf("No ext_rsds_CCSMB1-2001-2100.dat\n");  exit(1); }
	}
	
	if(SCENARIO_ID==1240){  /*  GISSE + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_GISSE_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_tas_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GISSE_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_pr_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GISSE_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_huss_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GISSE_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_rsds_GISSEA1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1241){  /*  GISSE + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_GISSE_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_tas_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GISSE_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_pr_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GISSE_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_huss_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GISSE_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_rsds_GISSEA1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1242){  /*  GISSE + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_GISSE_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_tas_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GISSE_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_pr_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GISSE_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_huss_GISSEA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GISSE_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_rsds_GISSEA1B-2001-2100.dat\n");  exit(1); }
	}
	
	if(SCENARIO_ID==1250){  /*  GISSR + A1B  */
		/* changed (R1 -> R4) 2009/04/07  */
		if( (fp_c2[0]=fopen("./data/tas_GISSR_20C-A1B_R4.dat","rt"))==NULL ){  printf("No tas_GISSR_20C-A1B_R4.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GISSR_20C-A1B_R4.dat","rt"))==NULL ){  printf("No pr_GISSR_20C-A1B_R4.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GISSR_20C-A1B_R4.dat","rt"))==NULL ){  printf("No huss_GISSR_20C-A1B_R4.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GISSR_20C-A1B_R4.dat","rt"))==NULL ){  printf("No rsds_GISSR_20C-A1B_R4.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1251){  /*  GISSR + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_GISSR_20C-A1B_R2.dat","rt"))==NULL ){  printf("No tas_GISSR_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GISSR_20C-A1B_R2.dat","rt"))==NULL ){  printf("No pr_GISSR_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GISSR_20C-A1B_R2.dat","rt"))==NULL ){  printf("No huss_GISSR_20C-A1B_R2.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GISSR_20C-A1B_R2.dat","rt"))==NULL ){  printf("No rsds_GISSR_20C-A1B_R2.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1252){  /*  GISSR + A2  */
		if( (fp_c2[0]=fopen("./data/tas_GISSR_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_GISSRA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GISSR_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_GISSRA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GISSR_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_GISSRA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GISSR_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_GISSRA1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1253){  /*  GISSR + B1  */
		/* changed (R2 -> R1) 2009/04/07  */
		if( (fp_c2[0]=fopen("./data/tas_GISSR_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_GISSR_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GISSR_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_GISSR_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GISSR_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_GISSR_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GISSR_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_GISSR_20C-B1_R1.dat\n");  exit(1); }
	}
	
	if(SCENARIO_ID==1260){  /*  PCM + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_PCM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_tas_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_PCM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_pr_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_PCM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_huss_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_PCM_20C-A1B_R1.dat","rt"))==NULL ){  printf("No ext_rsds_PCMA1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1261){  /*  PCM + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_PCM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_tas_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_PCM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_pr_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_PCM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_huss_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_PCM_20C-A1B_R2.dat","rt"))==NULL ){  printf("No ext_rsds_PCMA1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1262){  /*  PCM + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_PCM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_tas_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_PCM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_pr_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_PCM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_huss_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_PCM_20C-A1B_R3.dat","rt"))==NULL ){  printf("No ext_rsds_PCMA1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1263){  /*  PCM + A1B  */
		if( (fp_c2[0]=fopen("./data/tas_PCM_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_tas_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_PCM_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_pr_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_PCM_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_huss_PCMA1B-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_PCM_20C-A1B_R4.dat","rt"))==NULL ){  printf("No ext_rsds_PCMA1B-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1264){  /*  PCM + A2  */
		if( (fp_c2[0]=fopen("./data/tas_PCM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_tas_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_PCM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_pr_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_PCM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_huss_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_PCM_20C-A2_R1.dat","rt"))==NULL ){  printf("No ext_rsds_PCMA2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1265){  /*  PCM + A2  */
		if( (fp_c2[0]=fopen("./data/tas_PCM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_tas_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_PCM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_pr_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_PCM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_huss_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_PCM_20C-A2_R2.dat","rt"))==NULL ){  printf("No ext_rsds_PCMA2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1266){  /*  PCM + A2  */
		if( (fp_c2[0]=fopen("./data/tas_PCM_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_tas_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_PCM_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_pr_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_PCM_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_huss_PCMA2-2001-2100.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_PCM_20C-A2_R3.dat","rt"))==NULL ){  printf("No ext_rsds_PCMA2-2001-2100.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1267){  /*  PCM + A2  */
		if( (fp_c2[0]=fopen("./data/tas_PCM_20C-A2_R4.dat","rt"))==NULL ){  printf("No tas_PCM_20C-A2_R4.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_PCM_20C-A2_R4.dat","rt"))==NULL ){  printf("No pr_PCM_20C-A2_R4.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_PCM_20C-A2_R4.dat","rt"))==NULL ){  printf("No huss_PCM_20C-A2_R4.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_PCM_20C-A2_R4.dat","rt"))==NULL ){  printf("No rsds_PCM_20C-A2_R4.dat\n");  exit(1); }
	}
	
	if(SCENARIO_ID==1270){  /*  CSIRO 3.5 + A1  */
		if( (fp_c2[0]=fopen("./data/tas_CSIRO35_20C-A1_R1.dat","rt"))==NULL ){  printf("No tas_CSIRO35_20C-A1_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CSIRO35_20C-A1_R1.dat","rt"))==NULL ){  printf("No pr_CSIRO35_20C-A1_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CSIRO35_20C-A1_R1.dat","rt"))==NULL ){  printf("No huss_CSIRO35_20C-A1_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CSIRO35_20C-A1_R1.dat","rt"))==NULL ){  printf("No rsds_CSIRO35_20C-A1_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1271){  /*  CSIRO 3.5 + A2  */
		if( (fp_c2[0]=fopen("./data/tas_CSIRO35_20C-A2_R1.dat","rt"))==NULL ){  printf("No tas_CSIRO35_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CSIRO35_20C-A2_R1.dat","rt"))==NULL ){  printf("No pr_CSIRO35_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CSIRO35_20C-A2_R1.dat","rt"))==NULL ){  printf("No huss_CSIRO35_20C-A2_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CSIRO35_20C-A2_R1.dat","rt"))==NULL ){  printf("No rsds_CSIRO35_20C-A2_R1.dat\n");  exit(1); }
	}else if(SCENARIO_ID==1272){  /*  CSIRO 3.5 + B1  */
		if( (fp_c2[0]=fopen("./data/tas_CSIRO35_20C-B1_R1.dat","rt"))==NULL ){  printf("No tas_CSIRO35_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CSIRO35_20C-B1_R1.dat","rt"))==NULL ){  printf("No pr_CSIRO35_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CSIRO35_20C-B1_R1.dat","rt"))==NULL ){  printf("No huss_CSIRO35_20C-B1_R1.dat\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CSIRO35_20C-B1_R1.dat","rt"))==NULL ){  printf("No rsds_CSIRO35_20C-B1_R1.dat\n");  exit(1); }
	}
    
    if(SCENARIO_ID == 3000){  /*  GEOMIP:BNU-ESM base (RCP4.5)  */
		if( (fp_c2[0]=fopen("./data/tas_BNU-ESM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_BNU-ESM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_BNU-ESM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_BNU-ESM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_BNU-ESM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No huss_BNU-ESM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_BNU-ESM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_BNU-ESM_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(SCENARIO_ID==3001){  /*  GEOMIP: BNU-ESM G3  */
		if( (fp_c2[0]=fopen("./data/tas_BNU-ESM_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No tas_BNU-ESM_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_BNU-ESM_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No pr_BNU-ESM_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_BNU-ESM_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No huss_BNU-ESM_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_BNU-ESM_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No rsds_BNU-ESM_197001-210012_rcp45_g3.txt\n");  exit(1); }
	}else if(SCENARIO_ID==3002){  /*  GEOMIP: BNU-ESM G4  */
		if( (fp_c2[0]=fopen("./data/tas_BNU-ESM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No tas_BNU-ESM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_BNU-ESM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No pr_BNU-ESM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_BNU-ESM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No huss_BNU-ESM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_BNU-ESM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No rsds_BNU-ESM_197001-210012_rcp45_g4.txt\n");  exit(1); }
	}else if(SCENARIO_ID==3010){  /*  GEOMIP: CSIRO-Mk3L-1-2 base (RCP4.5)  */
		if( (fp_c2[0]=fopen("./data/tas_CSIRO-Mk3L-1-2_197001-210001_rcp45_base.txt","rt"))==NULL ){  printf("No tas_CSIRO-Mk3L-1-2_197001-210001_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CSIRO-Mk3L-1-2_197001-210001_rcp45_base.txt","rt"))==NULL ){  printf("No pr_CSIRO-Mk3L-1-2_197001-210001_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/hus_CSIRO-Mk3L-1-2_197001-210001_rcp45_base.txt","rt"))==NULL ){  printf("No hus_CSIRO-Mk3L-1-2_197001-210001_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CSIRO-Mk3L-1-2_197001-210001_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_CSIRO-Mk3L-1-2_197001-210001_rcp45_base.txt\n");  exit(1); }
	}else if(SCENARIO_ID==3011){  /*  GEOMIP: CSIRO-Mk3L-1-2 G4  */
		if( (fp_c2[0]=fopen("./data/tas_CSIRO-Mk3L-1-2_197001-210001_rcp45_g4.txt","rt"))==NULL ){  printf("No tas_CSIRO-Mk3L-1-2_197001-210001_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CSIRO-Mk3L-1-2_197001-210001_rcp45_g4.txt","rt"))==NULL ){  printf("No pr_CSIRO-Mk3L-1-2_197001-210001_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/hus_CSIRO-Mk3L-1-2_197001-210001_rcp45_g4.txt","rt"))==NULL ){  printf("No hus_CSIRO-Mk3L-1-2_197001-210001_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CSIRO-Mk3L-1-2_197001-210001_rcp45_g4.txt","rt"))==NULL ){  printf("No rsds_CSIRO-Mk3L-1-2_197001-210001_rcp45_g4.txt\n");  exit(1); }
	}else if(SCENARIO_ID==3012){  /*  GEOMIP: CSIRO-Mk3L-1-2 G3S  */
		if( (fp_c2[0]=fopen("./data/tas_CSIRO-Mk3L-1-2_197001-210001_rcp45_g3s.txt","rt"))==NULL ){  printf("No tas_CSIRO-Mk3L-1-2_197001-210001_rcp45_g3s.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CSIRO-Mk3L-1-2_197001-210001_rcp45_g3s.txt","rt"))==NULL ){  printf("No pr_CSIRO-Mk3L-1-2_197001-210001_rcp45_g3s.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/hus_CSIRO-Mk3L-1-2_197001-210001_rcp45_g3s.txt","rt"))==NULL ){  printf("No hus_CSIRO-Mk3L-1-2_197001-210001_rcp45_g3s.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CSIRO-Mk3L-1-2_197001-210001_rcp45_g3s.txt","rt"))==NULL ){  printf("No rsds_CSIRO-Mk3L-1-2_197001-210001_rcp45_g3s.txt\n");  exit(1); }
	}else if(SCENARIO_ID==3020){  /*  GEOMIP: GISS-E2-R base (RCP4.5)  */
		if( (fp_c2[0]=fopen("./data/tas_GISS-E2-R_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_GISS-E2-R_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GISS-E2-R_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_GISS-E2-R_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GISS-E2-R_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No huss_GISS-E2-R_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GISS-E2-R_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_GISS-E2-R_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(SCENARIO_ID==3021){  /*  GEOMIP: GISS-E2-R G3  */
		if( (fp_c2[0]=fopen("./data/tas_GISS-E2-R_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No tas_GISS-E2-R_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GISS-E2-R_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No pr_GISS-E2-R_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GISS-E2-R_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No huss_GISS-E2-R_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GISS-E2-R_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No rsds_GISS-E2-R_197001-210012_rcp45_g3.txt\n");  exit(1); }
	}else if(SCENARIO_ID==3022){  /*  GEOMIP: GISS-E2-R G4  */
		if( (fp_c2[0]=fopen("./data/tas_GISS-E2-R_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No tas_GISS-E2-R_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_GISS-E2-R_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No pr_GISS-E2-R_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_GISS-E2-R_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No huss_GISS-E2-R_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_GISS-E2-R_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No rsds_GISS-E2-R_197001-210012_rcp45_g4.txt\n");  exit(1); }
	}else if(SCENARIO_ID==3030){  /*  GEOMIP: HadGEM2-ES base (RCP4.5) */
		if( (fp_c2[0]=fopen("./data/tas_HadGEM2-ES_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_HadGEM2-ES_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_HadGEM2-ES_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_HadGEM2-ES_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_HadGEM2-ES_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No huss_HadGEM2-ES_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_HadGEM2-ES_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_HadGEM2-ES_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(SCENARIO_ID==3031){  /*  GEOMIP: HadGEM2-ES G3 */
		if( (fp_c2[0]=fopen("./data/tas_HadGEM2-ES_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No tas_HadGEM2-ES_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_HadGEM2-ES_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No pr_HadGEM2-ES_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_HadGEM2-ES_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No huss_HadGEM2-ES_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_HadGEM2-ES_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No rsds_HadGEM2-ES_197001-210012_rcp45_g3.txt\n");  exit(1); }
	}else if(SCENARIO_ID==3032){  /*  GEOMIP: HadGEM2-ES G4 */
		if( (fp_c2[0]=fopen("./data/tas_HadGEM2-ES_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No tas_HadGEM2-ES_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_HadGEM2-ES_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No pr_HadGEM2-ES_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_HadGEM2-ES_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No huss_HadGEM2-ES_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_HadGEM2-ES_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No rsds_HadGEM2-ES_197001-210012_rcp45_g4.txt\n");  exit(1); }
	}else if(SCENARIO_ID==3033){  /*  GEOMIP: HadGEM2-ES G3S */
		if( (fp_c2[0]=fopen("./data/tas_HadGEM2-ES_197001-210012_rcp45_g3s.txt","rt"))==NULL ){  printf("No tas_HadGEM2-ES_197001-210012_rcp45_g3s.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_HadGEM2-ES_197001-210012_rcp45_g3s.txt","rt"))==NULL ){  printf("No pr_HadGEM2-ES_197001-210012_rcp45_g3s.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_HadGEM2-ES_197001-210012_rcp45_g3s.txt","rt"))==NULL ){  printf("No huss_HadGEM2-ES_197001-210012_rcp45_g3s.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_HadGEM2-ES_197001-210012_rcp45_g3s.txt","rt"))==NULL ){  printf("No rsds_HadGEM2-ES_197001-210012_rcp45_g3s.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3034){  /*  GEOMIP: HadGEM2-ES G4cdnc */
		if( (fp_c2[0]=fopen("./data/tas_HadGEM2-ES_197001-210012_rcp45_G4cdnc.txt","rt"))==NULL ){  printf("No tas_HadGEM2-ES_197001-210012_rcp45_G4cdnc.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_HadGEM2-ES_197001-210012_rcp45_G4cdnc.txt","rt"))==NULL ){  printf("No pr_HadGEM2-ES_197001-210012_rcp45_G4cdnc.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_HadGEM2-ES_197001-210012_rcp45_G4cdnc.txt","rt"))==NULL ){  printf("No huss_HadGEM2-ES_197001-210012_rcp45_G4cdnc.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_HadGEM2-ES_197001-210012_rcp45_G4cdnc.txt","rt"))==NULL ){  printf("No rsds_HadGEM2-ES_197001-210012_rcp45_G4cdnc.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3035){  /*  GEOMIP: HadGEM2-ES G4seaSalt */
		if( (fp_c2[0]=fopen("./data/tas_HadGEM2-ES_197001-210012_rcp45_G4seaSalt.txt","rt"))==NULL ){  printf("No tas_HadGEM2-ES_197001-210012_rcp45_G4seaSalt.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_HadGEM2-ES_197001-210012_rcp45_G4seaSalt.txt","rt"))==NULL ){  printf("No pr_HadGEM2-ES_197001-210012_rcp45_G4seaSalt.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_HadGEM2-ES_197001-210012_rcp45_G4seaSalt.txt","rt"))==NULL ){  printf("No huss_HadGEM2-ES_197001-210012_rcp45_G4seaSalt.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_HadGEM2-ES_197001-210012_rcp45_G4seaSalt.txt","rt"))==NULL ){  printf("No rsds_HadGEM2-ES_197001-210012_rcp45_G4seaSalts.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3040){  /*  GEOMIP: IPSL-CM5A-LR base (RCP4.5) */
		if( (fp_c2[0]=fopen("./data/tas_IPSL-CM5A-LR_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_IPSL-CM5A-LR_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_IPSL-CM5A-LR_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_IPSL-CM5A-LR_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_IPSL-CM5A-LR_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No huss_IPSL-CM5A-LR_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_IPSL-CM5A-LR_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_IPSL-CM5A-LR_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3041){  /*  GEOMIP: IPSL-CM5A-LR G3 */
		if( (fp_c2[0]=fopen("./data/tas_IPSL-CM5A-LR_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No tas_IPSL-CM5A-LR_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_IPSL-CM5A-LR_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No pr_IPSL-CM5A-LR_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_IPSL-CM5A-LR_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No huss_IPSL-CM5A-LR_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_IPSL-CM5A-LR_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No rsds_IPSL-CM5A-LR_197001-210012_rcp45_g3.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3042){  /*  GEOMIP: IPSL-CM5A-LR G5 */
		if( (fp_c2[0]=fopen("./data/tas_IPSL-CM5A-LR_197001-210012_rcp45_g5.txt","rt"))==NULL ){  printf("No tas_IPSL-CM5A-LR_197001-210012_rcp45_g5.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_IPSL-CM5A-LR_197001-210012_rcp45_g5.txt","rt"))==NULL ){  printf("No pr_IPSL-CM5A-LR_197001-210012_rcp45_g5.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_IPSL-CM5A-LR_197001-210012_rcp45_g5.txt","rt"))==NULL ){  printf("No huss_IPSL-CM5A-LR_197001-210012_rcp45_g5.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_IPSL-CM5A-LR_197001-210012_rcp45_g5.txt","rt"))==NULL ){  printf("No rsds_IPSL-CM5A-LR_197001-210012_rcp45_g5.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3050){  /*  GEOMIP: MIROC-ESM base (RCP4.5)  */
		if( (fp_c2[0]=fopen("./data/tas_MIROC-ESM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_MIROC-ESM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MIROC-ESM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_MIROC-ESM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MIROC-ESM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No huss_MIROC-ESM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MIROC-ESM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_MIROC-ESM_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3051){  /*  GEOMIP: MIROC-ESM G4  */
		if( (fp_c2[0]=fopen("./data/tas_MIROC-ESM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No tas_MIROC-ESM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MIROC-ESM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No pr_MIROC-ESM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MIROC-ESM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No huss_MIROC-ESM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MIROC-ESM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No rsds_MIROC-ESM_197001-210012_rcp45_g4.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3052){  /*  GEOMIP: MIROC-ESM G4cdnc  */
		if( (fp_c2[0]=fopen("./data/tas_MIROC-ESM_197001-210001_rcp45_G4cdnc.txt","rt"))==NULL ){  printf("No tas_MIROC-ESM_197001-210001_rcp45_G4cdnc.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MIROC-ESM_197001-210001_rcp45_G4cdnc.txt","rt"))==NULL ){  printf("No pr_MIROC-ESM_197001-210001_rcp45_G4cdnc.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MIROC-ESM_197001-210001_rcp45_G4cdnc.txt","rt"))==NULL ){  printf("No huss_MIROC-ESM_197001-210001_rcp45_G4cdnc.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MIROC-ESM_197001-210001_rcp45_G4cdnc.txt","rt"))==NULL ){  printf("No rsds_MIROC-ESM_197001-210001_rcp45_G4cdnc.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3060){  /*  GEOMIP: MIROC-ESM-CHEM base (RCP4.5)  */
		if( (fp_c2[0]=fopen("./data/tas_MIROC-ESM-CHEM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_MIROC-ESM-CHEM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MIROC-ESM-CHEM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_MIROC-ESM-CHEM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MIROC-ESM-CHEM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No huss_MIROC-ESM-CHEM_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MIROC-ESM-CHEM_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_MIROC-ESM-CHEM_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3061){  /*  GEOMIP: MIROC-ESM-CHEM G4  */
		if( (fp_c2[0]=fopen("./data/tas_MIROC-ESM-CHEM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No tas_MIROC-ESM-CHEM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MIROC-ESM-CHEM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No pr_MIROC-ESM-CHEM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_MIROC-ESM-CHEM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No huss_MIROC-ESM-CHEM_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MIROC-ESM-CHEM_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No rsds_MIROC-ESM-CHEM_197001-210012_rcp45_g4.txt\n");  exit(1); }
	} else if(SCENARIO_ID == 3070){  /*  GEOMIP: CanESM2 base (RCP4.5)  */
		if( (fp_c2[0]=fopen("./data/tas_CanESM2_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_CanESM2_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CanESM2_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_CanESM2_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CanESM2_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No huss_CanESM2_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CanESM2_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_CanESM2_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3071){  /*  GEOMIP: CanESM2 G4  */
		if( (fp_c2[0]=fopen("./data/tas_CanESM2_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No tas_CanESM2_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CanESM2_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No pr_CanESM2_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CanESM2_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No huss_CanESM2_197001-210012_rcp45_g4.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CanESM2_197001-210012_rcp45_g4.txt","rt"))==NULL ){  printf("No rsds_CanESM2_197001-210012_rcp45_g4.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3072){  /*  GEOMIP: CanESM2 G4cdnc  */
		if( (fp_c2[0]=fopen("./data/tas_CanESM2_197001-210012_rcp45_g4cdnc.txt","rt"))==NULL ){  printf("No tas_CanESM2_197001-210012_rcp45_g4cdnc.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CanESM2_197001-210012_rcp45_g4cdnc.txt","rt"))==NULL ){  printf("No pr_CanESM2_197001-210012_rcp45_g4cdnc.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CanESM2_197001-210012_rcp45_g4cdnc.txt","rt"))==NULL ){  printf("No huss_CanESM2_197001-210012_rcp45_g4cdnc.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CanESM2_197001-210012_rcp45_g4cdnc.txt","rt"))==NULL ){  printf("No rsds_CanESM2_197001-210012_rcp45_g4cdnc.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3080){  /*  GEOMIP: MPI-ESM-LR base (RCP4.5)  */
		if( (fp_c2[0]=fopen("./data/tas_MPI-ESM-LR_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_MPI-ESM-LR_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MPI-ESM-LR_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_MPI-ESM-LR_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/hus_MPI-ESM-LR_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No hus_MPI-ESM-LR_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MPI-ESM-LR_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_MPI-ESM-LR_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3081){  /*  GEOMIP: MPI-ESM-LR G3  */
		if( (fp_c2[0]=fopen("./data/tas_MPI-ESM-LR_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No tas_MPI-ESM-LR_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_MPI-ESM-LR_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No pr_MPI-ESM-LR_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/hus_MPI-ESM-LR_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No hus_MPI-ESM-LR_197001-210012_rcp45_g3.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_MPI-ESM-LR_197001-210012_rcp45_g3.txt","rt"))==NULL ){  printf("No rsds_MPI-ESM-LR_197001-210012_rcp45_g3.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3090){  /*  GEOMIP: CCSM4 base (RCP4.5)  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM4_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_CCSM4_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM4_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_CCSM4_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM4_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No huss_CCSM4_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM4_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_CCSM4_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3091){  /*  GEOMIP: CCSM4 G3  */
		if( (fp_c2[0]=fopen("./data/tas_CCSM4_197001-210012_rcp45_g3s.txt","rt"))==NULL ){  printf("No tas_CCSM4_197001-210012_rcp45_g3s.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_CCSM4_197001-210012_rcp45_g3s.txt","rt"))==NULL ){  printf("No pr_CCSM4_197001-210012_rcp45_g3s.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_CCSM4_197001-210012_rcp45_g3s.txt","rt"))==NULL ){  printf("No huss_CCSM4_197001-210012_rcp45_g3s.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_CCSM4_197001-210012_rcp45_g3s.txt","rt"))==NULL ){  printf("No rsds_CCSM4_197001-210012_rcp45_g3s.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3100){  /*  GEOMIP: NorESM1-M base (RCP4.5)  */
		if( (fp_c2[0]=fopen("./data/tas_NorESM1-M_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No tas_NorESM1-M_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_NorESM1-M_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No pr_NorESM1-M_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_NorESM1-M_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No huss_NorESM1-M_197001-210012_rcp45_base.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_NorESM1-M_197001-210012_rcp45_base.txt","rt"))==NULL ){  printf("No rsds_NorESM1-M_197001-210012_rcp45_base.txt\n");  exit(1); }
	}else if(SCENARIO_ID == 3101){  /*  GEOMIP: NorESM1-M G4cdnc  */
		if( (fp_c2[0]=fopen("./data/tas_NorESM1-M_197001-210012_rcp45_g4cdnc.txt","rt"))==NULL ){  printf("No tas_NorESM1-M_197001-210012_rcp45_g4cdnc.txt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_NorESM1-M_197001-210012_rcp45_g4cdnc.txt","rt"))==NULL ){  printf("No pr_NorESM1-M_197001-210012_rcp45_g4cdnc.txt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/huss_NorESM1-M_197001-210012_rcp45_g4cdnc.txt","rt"))==NULL ){  printf("No huss_NorESM1-M_197001-210012_rcp45_g4cdnc.txt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_NorESM1-M_197001-210012_rcp45_g4cdnc.txt","rt"))==NULL ){  printf("No rsds_NorESM1-M_197001-210012_rcp45_g4cdnc.txt\n");  exit(1); }
	}
    
    if(SCENARIO_ID == 4011){  /*  PLUME: GFDL-ESM2M RCP4.5  */
		if( (fp_c2[0]=fopen("./data/tas_bced_1960_1999_gfdl-esm2m_rcp4p5_2006-2099.flt","rt"))==NULL ){  printf("No tas_bced_1960_1999_gfdl-esm2m_rcp4p5_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_bced_1960_1999_gfdl-esm2m_rcp4p5_2006-2099.flt","rt"))==NULL ){  printf("No pr_bced_1960_1999_gfdl-esm2m_rcp4p5_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/hurs_gfdl-esm2m_rcp4p5_2006-2099.flt","rt"))==NULL ){  printf("No hurs_gfdl-esm2m_rcp4p5_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_bced_1960_1999_gfdl-esm2m_rcp4p5_2006-2099.flt","rt"))==NULL ){  printf("No rsds_bced_1960_1999_gfdl-esm2m_rcp4p5_2006-2099.flt\n");  exit(1); }
	}else if(SCENARIO_ID == 4012){  /*  PLUME: IPSL-CM5A-LR RCP8.5  */
		if( (fp_c2[0]=fopen("./data/tas_bced_1960_1999_gfdl-esm2m_rcp8p5_2006-2099.flt","rt"))==NULL ){  printf("No tas_bced_1960_1999_gfdl-esm2m_rcp8p5_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_bced_1960_1999_gfdl-esm2m_rcp8p5_2006-2099.flt","rt"))==NULL ){  printf("No pr_bced_1960_1999_gfdl-esm2m_rcp8p5_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/hurs_gfdl-esm2m_rcp8p5_2006-2099.flt","rt"))==NULL ){  printf("No hurs_gfdl-esm2m_rcp8p5_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_bced_1960_1999_gfdl-esm2m_rcp8p5_2006-2099.flt","rt"))==NULL ){  printf("No rsds_bced_1960_1999_gfdl-esm2m_rcp8p5_2006-2099.flt\n");  exit(1); }
	}
    if(SCENARIO_ID == 4021){  /*  PLUME: IPSL-CM5A-LR RCP4.5  */
		if( (fp_c2[0]=fopen("./data/tas_bced_1960_1999_ipsl-cm5a-lr_rcp4p5_2006-2099.flt","rt"))==NULL ){  printf("No tas_bced_1960_1999_ipsl-cm5a-lr_rcp4p5_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_bced_1960_1999_ipsl-cm5a-lr_rcp4p5_2006-2099.flt","rt"))==NULL ){  printf("No pr_bced_1960_1999_ipsl-cm5a-lr_rcp4p5_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/hurs_ipsl-cm5a-lr_rcp4p5_2006-2099.flt","rt"))==NULL ){  printf("No hurs_ipsl-cm5a-lr_rcp4p5_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_bced_1960_1999_ipsl-cm5a-lr_rcp4p5_2006-2099.flt","rt"))==NULL ){  printf("No rsds_bced_1960_1999_ipsl-cm5a-lr_rcp4p5_2006-2099.flt\n");  exit(1); }
	}else if(SCENARIO_ID == 4022){  /*  PLUME: IPSL-CM5A-LR RCP8.5  */
		if( (fp_c2[0]=fopen("./data/tas_bced_1960_1999_ipsl-cm5a-lr_rcp8p5_2006-2099.flt","rt"))==NULL ){  printf("No tas_bced_1960_1999_ipsl-cm5a-lr_rcp8p5_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_bced_1960_1999_ipsl-cm5a-lr_rcp8p5_2006-2099.flt","rt"))==NULL ){  printf("No pr_bced_1960_1999_ipsl-cm5a-lr_rcp8p5_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/hurs_ipsl-cm5a-lr_rcp8p5_2006-2099.flt","rt"))==NULL ){  printf("No hurs_ipsl-cm5a-lr_rcp8p5_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_bced_1960_1999_ipsl-cm5a-lr_rcp8p5_2006-2099.flt","rt"))==NULL ){  printf("No rsds_bced_1960_1999_ipsl-cm5a-lr_rcp8p5_2006-2099.flt\n");  exit(1); }
	}else if(SCENARIO_ID == 4023){  /*  PLUME: IPSL-CM5A-LR RCP2.6  */
		if( (fp_c2[0]=fopen("./data/tas_bced_1960_1999_ipsl-cm5a-lr_rcp2p6_2006-2099.flt","rt"))==NULL ){  printf("No tas_bced_1960_1999_ipsl-cm5a-lr_rcp2p6_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_bced_1960_1999_ipsl-cm5a-lr_rcp2p6_2006-2099.flt","rt"))==NULL ){  printf("No pr_bced_1960_1999_ipsl-cm5a-lr_rcp2p6_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/hurs_ipsl-cm5a-lr_rcp2p6_2006-2099.flt","rt"))==NULL ){  printf("No hurs_ipsl-cm5a-lr_rcp2p6_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_bced_1960_1999_ipsl-cm5a-lr_rcp2p6_2006-2099.flt","rt"))==NULL ){  printf("No rsds_bced_1960_1999_ipsl-cm5a-lr_rcp2p6_2006-2099.flt\n");  exit(1); }
	}else if(SCENARIO_ID == 4024){  /*  PLUME: IPSL-CM5A-LR RCP6.0  */
		if( (fp_c2[0]=fopen("./data/tas_bced_1960_1999_ipsl-cm5a-lr_rcp6p0_2006-2099.flt","rt"))==NULL ){  printf("No tas_bced_1960_1999_ipsl-cm5a-lr_rcp6p0_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[1]=fopen("./data/pr_bced_1960_1999_ipsl-cm5a-lr_rcp6p0_2006-2099.flt","rt"))==NULL ){  printf("No pr_bced_1960_1999_ipsl-cm5a-lr_rcp6p0_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[2]=fopen("./data/hurs_ipsl-cm5a-lr_rcp6p0_2006-2099.flt","rt"))==NULL ){  printf("No hurs_ipsl-cm5a-lr_rcp6p0_2006-2099.flt\n");  exit(1); }
		if( (fp_c2[3]=fopen("./data/rsds_bced_1960_1999_ipsl-cm5a-lr_rcp6p0_2006-2099.flt","rt"))==NULL ){  printf("No rsds_bced_1960_1999_ipsl-cm5a-lr_rcp6p0_2006-2099.flt\n");  exit(1); }
	}
	
	if(SCENARIO_ID >= 1 && SCENARIO_ID <=3999){
        /* altitude */
		alt = (grid->topo>=0.0)?grid->topo:0.0;
        
		for(f=0;f<DL_FUTURE;f++){   /*  131 => 1970-2100 */
			for(g=0;g<ASTEP;g++){
				fscanf(fp_c2[0],"%ld %ld", &yr, &mon);
				fscanf(fp_c2[1],"%ld %ld", &yr, &mon);
				fscanf(fp_c2[2],"%ld %ld", &yr, &mon);
				fscanf(fp_c2[3],"%ld %ld", &yr, &mon);
				for(h=0;h<GCM_R;h++){
					for(i=0;i<GCM_C;i++){
						fscanf(fp_c2[0],"%f", &(grid->proj_tmp2m[f][g][h][i]));
						fscanf(fp_c2[1],"%f", &(grid->proj_prec[f][g][h][i]));
						fscanf(fp_c2[2],"%f", &(grid->proj_hum[f][g][h][i]));
						fscanf(fp_c2[3],"%f", &(grid->proj_rad[f][g][h][i]));
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
						if(SCENARIO_ID >= 1000){
							grid->proj_prec[f][g][h][i] *= HSN * DHN;
						}
						
						/* specific humidity (kg kg-1) to vapor pressure (hPa) */
						/* revided  (2009/08/17) */
						atmp = grid->proj_tmp2m[f][g][h][i];
						apres = 1013.25 * exp(-1.0*(28.964*0.001)*9.8*alt / (8.3144*(atmp+ZAT)));
						shum = grid->proj_hum[f][g][h][i];
						grid->proj_hum[f][g][h][i] = apres * grid->proj_hum[f][g][h][i]/(0.622 + 0.378*grid->proj_hum[f][g][h][i]);
					}
				}
			}
		}
	}else if(SCENARIO_ID >= 4000 && SCENARIO_ID <=4999){
        /* PLUME: 2014/07/31  */
        /* read sequentially for each grid in clim_cru.c */
    
        for(f=0;f<DL_FUTURE;f++){
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
    }else{
        /* no SCENARIO_ID data */
		for(f=0;f<DL_FUTURE;f++){
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
}

/****************************************/
void f_read_chaser_ndepo(
    FILE *fp_s[IFILEN],
	struct Grid  *grid
){
    long f,g,h;
    
	for(h=0;h<64;h++){
		for(g=0;g<128;g++){
			for(f=0;f<12;f++){
				fscanf(fp_s[52],"%lf", &(ndepo_chaser_dnhx[f][h][g]));
				fscanf(fp_s[52],"%lf", &(ndepo_chaser_dnoy[f][h][g]));
				fscanf(fp_s[52],"%lf", &(ndepo_chaser_wnhx[f][h][g]));
				fscanf(fp_s[52],"%lf", &(ndepo_chaser_wnoy[f][h][g]));
			}
            
            /* preindustrial */
            for(f=0;f<12;f++){
                fscanf(fp_s[85],"%lf", &(ndepo_chaser4_nhx_h[f][h][g]));
            }
            for(f=0;f<12;f++){
                fscanf(fp_s[85],"%lf", &(ndepo_chaser4_noy_h[f][h][g]));
            }
            for(f=0;f<12;f++){
                fscanf(fp_s[85],"%lf", &(ndepo_chaser4_ont_h[f][h][g]));
            }
            /* contemporary */
            for(f=0;f<12;f++){
                fscanf(fp_s[86],"%lf", &(ndepo_chaser4_nhx_p[f][h][g]));
            }
            for(f=0;f<12;f++){
                fscanf(fp_s[86],"%lf", &(ndepo_chaser4_noy_p[f][h][g]));
            }
            for(f=0;f<12;f++){
                fscanf(fp_s[86],"%lf", &(ndepo_chaser4_ont_p[f][h][g]));
            }
		}
	}	
    
}
