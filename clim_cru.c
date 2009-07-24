/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* UEA/CRU data  **********************************************/
/*
New, M., Lister, D., Hulme, M. and Makin, I., 2002. A high-resolution data set of 
surface climate over global land areas. Climate Research, 21:1-25.
*/
void read_cru_clim(
	FILE *fp_c[4], 
	struct Grid *grid
){
	long kk[4];
	long f, g, h;
	double data;
	
	/*  printf("reading CRU data...");  */
	
	/* read CRU TS2.1 Cloud data */
	fscanf(fp_c[0],"%ld", &kk[0]);
	if(kk[0]==11){
		for(h=0;h<102;h++){
			for(g=0;g<ASTEP;g++){
				fscanf(fp_c[0],"%lf", &data);
				
				grid->hist_cld[h][g] = data/100.0;
			}
		}
	}
	
	/* read CRU TS2.1 Precipitation data */
	fscanf(fp_c[1],"%ld", &kk[1]);
	if(kk[1]==11){
		for(h=0;h<102;h++){
			for(g=0;g<ASTEP;g++){
				fscanf(fp_c[1],"%lf", &data);
				
				grid->hist_pre[h][g] = data;
			}
		}
	}
	
	/* read CRU TS2.1 Temperature data */
	fscanf(fp_c[2],"%ld", &kk[2]);
	if(kk[2]==11){
		for(h=0;h<102;h++){
			for(g=0;g<ASTEP;g++){
				fscanf(fp_c[2],"%lf", &data);
				
				grid->hist_tmp[h][g] = data;
			}
		}
	}
	
	/* read CRU TS2.1 Vapor-pressure data */
	fscanf(fp_c[3],"%ld", &kk[3]);
	if(kk[3]==11){
		for(h=0;h<102;h++){
			for(g=0;g<ASTEP;g++){
				fscanf(fp_c[3],"%lf", &data);
				
				grid->hist_vap[h][g] = data;
			}
		}
	}
	
	/*  printf("********* %ld\n",cru_flag);  */
	
	/* if valid CRU climate data are all available **/
	if((kk[0]+kk[1]+kk[2]+kk[3])==44){
		/* data available */
		grid->cru_exist = 1;

		/* base climate (average 1971 - 2000) ******************/
		for(g=0;g<ASTEP;g++){
			grid->hist_cld_b[g] = grid->hist_pre_b[g] = 0.0;
			grid->hist_tmp_b[g] = grid->hist_vap_b[g] = 0.0;
		}
		for(f=0;f<30;f++){
			for(g=0;g<ASTEP;g++){
				grid->hist_cld_b[g] += grid->hist_cld[70+f][g]/30.0;
				grid->hist_pre_b[g] += grid->hist_pre[70+f][g]/30.0;
				grid->hist_vap_b[g] += grid->hist_vap[70+f][g]/30.0;
				grid->hist_tmp_b[g] += grid->hist_tmp[70+f][g]/30.0;
			}
		}
	}else{
		/* unavailable CRU TS2.1 data, for example on ocean */
		grid->cru_exist = 0;
	}
	
	/*  printf("done\n");  */
}
