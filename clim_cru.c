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
void read_cru_clim(
	FILE *fp_c[40], 
	struct Grid *grid
){
	long cru_flag, kk;
	long f, g, h;
	double data;
	
	/*  printf("reading CRU data...");  */
	
	cru_flag = 0;
	for(f=0;f<40;f++){
		fscanf(fp_c[f],"%ld", &kk);
		cru_flag += kk;
		/*  printf("%ld ",kk);  */
		
		if(kk==1){
			for(h=0;h<10;h++){
				for(g=0;g<12;g++){
					fscanf(fp_c[f],"%lf",&data);
					
					if(f>=0&&f<10){
						grid->hist_cld[f*10+h][g] = data/100.0;
					}else if(f>=10&&f<20){
						grid->hist_pre[(f-10)*10+h][g] = data;
					}else if(f>=20&&f<30){
						grid->hist_tmp[(f-20)*10+h][g] = data;
					}else if(f>=30&&f<40){
						grid->hist_vap[(f-30)*10+h][g] = data;
					}
				}
			}
		}
	}
	
	/*  printf("********* %ld\n",cru_flag);  */
	
	/* if valid CRU climate data are all available **/
	if(cru_flag==40){
		grid->cru_exist = 1;

		/******* base climate (average 1971 - 2000) ******/
		for(g=0;g<12;g++){
			grid->hist_cld_b[g] = grid->hist_pre_b[g] = 0.0;
			grid->hist_tmp_b[g] = grid->hist_vap_b[g] = 0.0;
		}
		for(f=0;f<30;f++){
			for(g=0;g<12;g++){
				grid->hist_cld_b[g] += grid->hist_cld[70+f][g]/30.0;
				grid->hist_pre_b[g] += grid->hist_pre[70+f][g]/30.0;
				grid->hist_vap_b[g] += grid->hist_vap[70+f][g]/30.0;
			}
		}
		for(f=0;f<30;f++){
			for(g=0;g<12;g++){
				grid->hist_tmp_b[g] += grid->hist_tmp[70+f][g]/30.0;
			}
		}
	}else{
		grid->cru_exist = 0;
	}
	
	/*  printf("done\n");  */
}
