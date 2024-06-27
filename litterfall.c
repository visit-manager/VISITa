/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed  in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* functions of litterfall */
#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* turnover of plant organs ***************************************/
void f_mortality(
	struct Grid *grid, 
    struct Loct *loct, 
	struct Pchar *pchar
){
	double ftmp, fdrt;
	
	/* larger at warm, and smaller at cool */
    ftmp = 1.0 + (grid->tmp_sfc[grid->m] - grid->tmp_sfc_am)/100.0;
    
    fdrt = 1.0;
    /* experiments for drought-induced mortality: 2024/06/27 */
    if(EX_DROUGHT_MOTAL == 1){
        if(loct->sww < 100.0){
            fdrt = 1.0 - 0.001 * (100.0 - loct->sww);
        }
    }
	
	pchar->lf[grid->m] = pchar->lf0 * ftmp * fdrt; /* leaf */
	pchar->lc[grid->m] = pchar->lc0 * ftmp * fdrt; /* stem */
	pchar->lr[grid->m] = pchar->lr0 * ftmp * fdrt; /* root */
}

/* from foliage *******************************************************/
double flf(
	struct Grid *grid, 
	struct Pchar *veg, 
	struct Pmas *mass
){
	double lf;
	
	/* constant fraction of senescence organ, abandoned as litter */
	lf = veg->lf[grid->m] * mass->fol;
	
	return(lf);
}

/* from stem and branch *******************************************/		
double flc(
	struct Grid *grid, 
	struct Pchar *veg, 
	struct Pmas *mass
){
	double lc;
	
	/* constant fraction of senescence organ, abandoned as litter */
	lc = veg->lc[grid->m] * mass->stm;
	
	return(lc);
}

/* from root system ***********************************************/
double flr(
	struct Grid *grid, 
	struct Pchar *veg, 
	struct Pmas *mass
){
	double lr;
	
	/* constant fraction of senescence organ, abandoned as litter */
	lr = veg->lr[grid->m] * mass->rot;
	
	return(lr);
}
