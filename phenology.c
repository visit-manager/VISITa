/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
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

/******* growth period and phenological process *********/
void growthperiod(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	/** cumulative temperature, degree days**/
	if(grid->lat>=0.0){
		if(grid->m==0){
			pchar->gdd=0.0;
			if(grid->y==0){
				pchar->frag_dcd = 1;
				pchar->frag_emg = 0;
			}
		}
	}else if(grid->lat<0.0){
		if(grid->m==6){
			pchar->gdd=0.0;
			if(grid->y==0){
				pchar->frag_dcd = 0;
				pchar->frag_emg = 1;
			}
		}
	}
	if(grid->tmp_sfc[grid->m]>5.0){
		pchar->gdd+=grid->tmp_sfc[grid->m]*(double)grid->mm[grid->m];
	}	
	pchar->mgdd[grid->m]=pchar->gdd;
	
	/** growing period, days **/
	if(grid->m==0) pchar->grw_pd = 0.0;
	if(grid->tmp_sfc[grid->m] > 5.0){
		pchar->grw_pd += (double)grid->mm[grid->m];
	}	

	/** phenology **/
	switch(grid->veg_olson){
		case 0: case 33:
			phenology_bareland(grid, pchar);
			break;
		case 1: case 2: case 6: case 7: case 8: case 9: case 11: case 13: case 20: 
		case 21: case 22: case 23: case 24: case 25: case 26: 
			phenology_evergreen(grid, pchar);
			break;			
		case 4: case 5: case 10: case 12: 
			phenology_colddeciduous(grid, loct, pchar);
			break;
		case 3: 
			phenology_ariddeciduous(grid, pchar);
			break;
		case 14: case 15: case 16: case 17: case 18: case 19: case 27: case 28: 
			phenology_grass(grid, loct, pchar);
			break;
		case 29: case 30: case 31: case 32:
			phenology_agriculture(grid, loct, pchar);
			break;
	}
}

/******** bare land *********/
void phenology_bareland(
	struct Grid *grid, 
	struct Pchar *pchar
){
	pchar->season[grid->m] = 0;
}

/******** evergreen biomes ********/
void phenology_evergreen(
	struct Grid *grid, 
	struct Pchar *pchar
){
	if(grid->tmp_sfc[grid->m]>=5.0){ 
		/* warmer months */
		pchar->season[grid->m] = 1;
	}else{
		/* cooler months */
		pchar->season[grid->m] = 0;
	}
}

/******** Deciduous biomes in higher latitudes *********/
void phenology_colddeciduous(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	
	double crit_tem,crit_gdd;
	
	if(grid->veg_olson==4||grid->veg_olson==5){
		crit_tem = 5.0;
		crit_gdd = 300.0; /* 200->300, 2000/10/30*/
	}else if(grid->veg_olson==10||grid->veg_olson==12){
		crit_tem = 3.0;
		crit_gdd = 200.0; /* 120->200, 2000/10/30*/
	}
	
	if(grid->tmp_sfc[grid->m]<crit_tem){
		/* dormancy */
		pchar->season[grid->m] = 0;
		/* leaf-shedding */
		if(pchar->frag_emg==1&&pchar->frag_dcd==0){
			pchar->season[grid->m] = 3;
			pchar->frag_dcd = 1;
			pchar->frag_emg = 0;
		}
	}else if(grid->tmp_sfc[grid->m]>=crit_tem){
		/* growing-period */
		pchar->season[grid->m] = 1;
		/* leaf-emergence */
		if(pchar->gdd>crit_gdd){
			if(pchar->frag_emg==0){
				pchar->season[grid->m] = 2;
				pchar->frag_emg = 1;
				pchar->frag_dcd = 0;
			}
		}
	}
}

/*** arid deciduous biomes in lower latitudes ***/
void phenology_ariddeciduous(
	struct Grid *grid, 
	struct Pchar *pchar
){
	if(grid->prate_sfc[grid->m] < 50.0){
		/* dormancy */
		pchar->season[grid->m] = 0;
		/* leaf-shedding */
		if(pchar->frag_emg==1 && pchar->frag_dcd==0){
			pchar->season[grid->m] = 3;
			pchar->frag_dcd = 1;
			pchar->frag_emg = 0;
		}
	}else if(grid->prate_sfc[grid->m] >= 50.0){
		/* growing-period */
		pchar->season[grid->m] = 1;
		/* leaf-emergence */
		if(pchar->frag_emg == 0){
			pchar->season[grid->m] = 2;
			pchar->frag_emg = 1;
			pchar->frag_dcd = 0;
		}
	}
}

/******* grassland *******/
void phenology_grass(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	if(pchar->phototype==3){
		if(loct->msww[grid->m]/grid->field_cap2<0.1||grid->tmp_sfc[grid->m]<5.0){
			/* dormancy */
			pchar->season[grid->m]=0;
			/* leaf-shedding */
			if(pchar->frag_emg==1&&pchar->frag_dcd==0){
				pchar->season[grid->m]=3;
				pchar->frag_dcd=1;
				pchar->frag_emg=0;
			}
		}else if(loct->msww[grid->m]/grid->field_cap2>=0.1&&grid->tmp_sfc[grid->m]>=5.0){
			/* growing-period */
			pchar->season[grid->m]=1;
			/* leaf-emergence */
			if(pchar->frag_emg==0){
				pchar->season[grid->m]=2;
				pchar->frag_emg=1;
				pchar->frag_dcd=0;
			}
		}
	}else if(pchar->phototype==4){
		if(loct->msww[grid->m]/grid->field_cap2<0.1||grid->tmp_sfc[grid->m]<8.0){
			/* dormancy */
			pchar->season[grid->m]=0;
			/* leaf-shedding */
			if(pchar->frag_emg==1&&pchar->frag_dcd==0){
				pchar->season[grid->m]=3;
				pchar->frag_dcd=1;
				pchar->frag_emg=0;
			}
		}else if(loct->msww[grid->m]/grid->field_cap2>=0.1&&grid->tmp_sfc[grid->m]>=8.0){
			/* growing-period */
			pchar->season[grid->m]=1;
			/* leaf-emergence */
			if(pchar->frag_emg==0){
				pchar->season[grid->m]=2;
				pchar->frag_emg=1;
				pchar->frag_dcd=0;
			}
		}
	}
}

/******* crop land *******/
void phenology_agriculture(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	
	double crit_tem;
	
	/* critical temperature limiting growing period*/
	switch(grid->veg_olson){
		case 29:
			crit_tem=12.0; /* rice */
			break;
		case 30:
			crit_tem=7.0; /* cool crop */ /* 5->7, 2000/10/30*/
			break;
		case 31:
			crit_tem=10.0; /* warm crop */ /* 8->10, 2000/10/30*/
			break;
		case 32:
			crit_tem=10.0; /* irrigated crop */
			break;
	}
	
	if(grid->tmp_sfc_am>20.0&&grid->tmp_sfc_mn>10.0){
		/* multiple cropping */
		pchar->season[grid->m]=1;	
	}else{	
		/* seasonal single cropping */	
		if(grid->tmp_sfc[grid->m]<crit_tem){
			/* dormancy */
			pchar->season[grid->m] = 0;
			/* leaf-shedding */
			if(pchar->frag_emg==1 && pchar->frag_dcd==0){
				pchar->season[grid->m] = 3;
				pchar->frag_dcd = 1;
				pchar->frag_emg = 0;
			}
		}else if(grid->tmp_sfc[grid->m] >= crit_tem){
			/* growing-period */
			pchar->season[grid->m] = 1;
			/* leaf-emergence */
			if(pchar->frag_emg==0){
				pchar->season[grid->m] = 2;
				pchar->frag_emg = 1;
				pchar->frag_dcd = 0;
			}
		}
	}
}
