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

/**************************************************************************************/
void read_ncep_clim(
	struct Grid  *grid
){
	long f, g, h, i;
	long yr, mon;
	FILE *fp_t, *fp_p, *fp_h, *fp_r;
	
	/* if( (fp_t=fopen("./data/air.2m.gauss.mon4809.dat","rt"))==NULL ){  
		printf("No air.2m.gauss.mon4809.dat\n");  exit(1); }
	if( (fp_p=fopen("./data/prate.sfc.gauss.mon4809.dat","rt"))==NULL ){  
		printf("No prate.sfc.gauss.mon4809.dat\n");  exit(1); }
	if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4809.dat","rt"))==NULL ){  
		printf("No vpres.sfc.gauss.mon4809.dat\n");  exit(1); }
	if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4809.dat","rt"))==NULL ){  
		printf("No tcdc.eatm.gauss.mon4809.dat\n");  exit(1); } */
	
	/* if( (fp_t=fopen("./data/air.2m.gauss.mon4810.dat","rt"))==NULL ){  
		printf("No air.2m.gauss.mon4810.dat\n");  exit(1); }
	if( (fp_p=fopen("./data/prate.sfc.gauss.mon4810.dat","rt"))==NULL ){  
		printf("No prate.sfc.gauss.mon4810.dat\n");  exit(1); }
	if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4810.dat","rt"))==NULL ){  
		printf("No vpres.sfc.gauss.mon4810.dat\n");  exit(1); }
	if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4810.dat","rt"))==NULL ){  
		printf("No tcdc.eatm.gauss.mon4810.dat\n");  exit(1); } */
	
    /* extension to 2011: 2012/01/26 by A.Ito */
	if( (fp_t=fopen("./data/air.2m.gauss.mon4811.dat","rt"))==NULL ){  
		printf("No air.2m.gauss.mon4811.dat\n");  exit(1); }
	if( (fp_p=fopen("./data/prate.sfc.gauss.mon4811.dat","rt"))==NULL ){  
		printf("No prate.sfc.gauss.mon4811.dat\n");  exit(1); }
	if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4811.dat","rt"))==NULL ){  
		printf("No vpres.sfc.gauss.mon4811.dat\n");  exit(1); }
	if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4811.dat","rt"))==NULL ){  
		printf("No tcdc.eatm.gauss.mon4811.dat\n");  exit(1); }
	
	for(f=0;f<NCEP_TS;f++){   /*  1948-2011 */
		for(g=0;g<ASTEP;g++){
			fscanf(fp_t,"%ld %ld", &yr, &mon);
			fscanf(fp_p,"%ld %ld", &yr, &mon);
			fscanf(fp_h,"%ld %ld", &yr, &mon);
			fscanf(fp_r,"%ld %ld", &yr, &mon);
			for(h=0;h<94;h++){
				for(i=0;i<192;i++){
					fscanf(fp_t,"%f", &(grid->ncep_tmp2m[f][g][h][i]));
					fscanf(fp_p,"%f", &(grid->ncep_prate[f][g][h][i]));
					fscanf(fp_h,"%f", &(grid->ncep_vpres[f][g][h][i]));
					fscanf(fp_r,"%f", &(grid->ncep_tcdc[f][g][h][i]));
					
					grid->ncep_tmp2m[f][g][h][i] -= ZAT;
					grid->ncep_tcdc[f][g][h][i] *= 0.01;
				}
			}
		}
	}
	
	/* average 1970-1999 */
	for(f=0;f<30;f++){
		for(g=0;g<ASTEP;g++){
			for(h=0;h<94;h++){
				for(i=0;i<192;i++){
					if(f==0){
						grid->ncep_tmp2m_b[g][h][i] = 0.0;
						grid->ncep_prate_b[g][h][i] = 0.0;
						grid->ncep_vpres_b[g][h][i] = 0.0;
						grid->ncep_tcdc_b[g][h][i] = 0.0;
					}
					
					grid->ncep_tmp2m_b[g][h][i] += grid->ncep_tmp2m[f+22][g][h][i]/30.0;
					grid->ncep_prate_b[g][h][i] += grid->ncep_prate[f+22][g][h][i]/30.0;
					grid->ncep_vpres_b[g][h][i] += grid->ncep_vpres[f+22][g][h][i]/30.0;
					grid->ncep_tcdc_b[g][h][i] += grid->ncep_tcdc[f+22][g][h][i]/30.0;
				}
			}
		}
	}
	
	fclose(fp_t);
	fclose(fp_p);
	fclose(fp_h);
	fclose(fp_r);
}
