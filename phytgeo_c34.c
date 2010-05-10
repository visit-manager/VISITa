/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/** C3/C4 geography **/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* portion of C3 and C4 plant *****************************************/
void c34composition(
	short v_type, 
	struct Grid *grid, 
	struct Loct *loct
){
	long f;
	double tmp_factor, pre_factor;
	
	if(v_type == 1){  /* Olson */
		for(f=0;f<ASTEP;f++){
			/* empirical tharmal and moisture gradient, annual */
			switch(grid->veg_olson){
				case 14: case 15: case 16: case 17: case 18: case 19: case 27: case 28:
					if(grid->tmp_sfc_am<-5.0){
						loct->c4ptn[f]=0.0;
						loct->c3ptn[f]=1.0;
					}else if(grid->tmp_sfc_am>=-5.0){
						tmp_factor = (grid->tmp_sfc_am+5.0)*0.032;
						pre_factor = 1.0/(grid->prate_sfc_ann*0.00025+0.65);
						
						if((tmp_factor*pre_factor)<=0.96){
							loct->c4ptn[f] = tmp_factor*pre_factor;
							loct->c3ptn[f] = 1.0-loct->c4ptn[f];
						}else if((tmp_factor*pre_factor)>0.96){
							loct->c4ptn[f] = 0.96;
							loct->c3ptn[f] = 0.04;
						}
					}
					break;
				case 30: case 31: case 32:
					loct->c4ptn[f] = 0.2;
					loct->c3ptn[f] = 0.8;
					break;
				default:
					loct->c4ptn[f] = 0.0;
					loct->c3ptn[f] = 1.0;
			}
		}
	}else if(v_type==2){ /* crop */
		for(f=0;f<ASTEP;f++){
			if(grid->veg_crop==1 || grid->veg_crop==2){
				loct->c4ptn[f] = 0.0;
				loct->c3ptn[f] = 1.0;
			}else if(grid->veg_crop==3){
				loct->c4ptn[f] = 1.0;
				loct->c3ptn[f] = 0.0;
			}
		}
	}
}

/* monthly C3 - C4 alternation ******************************************/
void phenol_tran_c34(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){
	
	flux->lf_c[grid->m] = 0.0;
}
