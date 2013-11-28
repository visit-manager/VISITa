/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/** This file gives description of functions releted to photosynthesis, or GPP **/
#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

extern short RAD_SENS;

/* formula of daily gross primary production ****************************/
double fgpp(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar,
	struct Pmas *mass
){
	double gpp;
	double cc1, cc2, cc3, bb;
	
	/* give daily GPP based on the Monsi-Saeki theory, or Kuroiwa's equation */
	/*
	Kuroiwa, S., 1966. Dry matter production of plants. Ecology and Evolution. 
	Iwanami Shoten, Tokyo, pp. 71-100.
	*/
	if(pchar->psat[grid->m] > 0.0){
		cc1 = 2.0 * pchar->psat[grid->m]*grid->dlen[grid->m]*lTs/pchar->eK[grid->m]; 
		bb = pchar->eK[grid->m]*pchar->lue[grid->m]*grid->par[grid->m]/pchar->psat[grid->m];
		cc2 = 1.0 + sqrt(1.0 + bb);
		cc3 = 1.0 + sqrt(1.0 + bb*exp(-1.0*pchar->eK[grid->m]*mass->lai[grid->m]));
		gpp = cc1 * log(cc2 / cc3);
	}else{
		gpp = 0.0;
	}
	
	return (gpp);
}

/* light-saturated photosynthetic rate ************************************/
void pc_sat(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	double ftem, fstl, fnstl, fozone;
	double aa1, aa2, aa3;
	double f_acclim;
	
	/** optimum photosynthesis temperature **/
	if(pchar->phototype==3){ /* C3 plants with change */
		pchar->topt = pchar->topt0 + 0.01*pchar->ci[grid->m];	
	}else if(pchar->phototype==4){ /* C4 plants without change */
		pchar->topt = pchar->topt0;
	}
	
	/** CO2 compensation point **/
	if(pchar->phototype==3){ /* C3 plants with change: Brooks&Farquhar(1985) */
		aa3 = 0.000347 * (grid->tmp_sfc[grid->m] - 20.0)*(grid->tmp_sfc[grid->m] - 20.0);
		aa1 = 1.0 + 0.0451 * (grid->tmp_sfc[grid->m] - 20.0) + aa3;
		aa1 = (aa1>0.0)?aa1:0.0;
		pchar->cmpcd[grid->m] = pchar->cmpcd0*aa1;
	}else if(pchar->phototype==4){ /* C4 plants without change */
		pchar->cmpcd[grid->m] = pchar->cmpcd0;
	}
		
	/** temperature effect ******************************************************/
	if(grid->tmp_sfc[grid->m]<=pchar->tmax && grid->tmp_sfc[grid->m]>=pchar->tmin){
		aa1 = (grid->tmp_sfc[grid->m] - pchar->tmax)*(grid->tmp_sfc[grid->m] - pchar->tmin);
		aa2 = (grid->tmp_sfc[grid->m] - pchar->topt)*(grid->tmp_sfc[grid->m] - pchar->topt);
		ftem = aa1/(aa1 - aa2);
	}else{
		ftem = 0.0;
	}
	ftem = (ftem<=1.0)?ftem:1.0; 
	ftem = (ftem>=0.0)?ftem:0.0;
	
	pchar->ft[grid->m] = ftem;
	
	/** CO2 effect ***********************************************************/
	/** stomatal limitation via intercellular CO2 concentration **/
	if(pchar->phototype==3){ /* C3 plants */
		fstl = 0.05 + 0.95*(pchar->ci[grid->m] - pchar->cmpcd[grid->m])/(pchar->kmci + pchar->ci[grid->m]); 
	}else if(pchar->phototype==4){ /* C4 plants */
		fstl = 0.60 + 0.40*(pchar->ci[grid->m] - pchar->cmpcd[grid->m])/(pchar->kmci + pchar->ci[grid->m]); 
	}
	fstl = (fstl<=1.0)?fstl:1.0; 
	fstl = (fstl>=0.0)?fstl:0.0;
	
	pchar->fcd[grid->m] = fstl;

	/** soil water effect ******************************************************/
	/** non-stomatal limitation **/
	if(pchar->phototype==3){ /* C3 plants */
		fnstl = 0.95*loct->sww / (loct->sww + grid->field_cap2*pchar->km_nstl) + 0.05;
	}else if(pchar->phototype==4){ /* C4 plants */
		fnstl = 0.86*loct->sww / (loct->sww + grid->field_cap2*pchar->km_nstl) + 0.14;
	}
	fnstl = (fnstl<=1.0)?fnstl:1.0; fnstl=(fnstl>=0.0)?fnstl:0.0;
		
	pchar->fsw[grid->m] = fnstl;
    
    /* Ozone impact on photosynthesis and stomata: 2013/02/25 by A.Ito **/
    /* parameterization by:
    Sitch, S., P. M. Cox, W. J. Collins, and C. Huntingford (2007), 
    Indirect radiative forcing of climate change through ozone effects on the 
    land-carbon sink, Nature, 448, 791-794. */
    
    fozone = (loct->ao3*1000.0/22.4) / (loct->r_aero[grid->m] + 1.67 / (pchar->gs[grid->m]*0.0224/1000.0));
    
    fozone -= 3.0; /* biome-specific */
    if(fozone < 0.0){
        fozone = 0.0;
    }
    fozone = 1.0 - 1.0*fozone; /* biome-specific */
    if(fozone < 0.0){
        fozone = 0.0;
    }
    if(fozone > 1.0){
        fozone = 1.0;
    }
    pchar->fo3[grid->m] = fozone;
		
	/* leaf N effect ***********************************************************/
	if(CN_COUPLE >= 1){		
		pchar->pmax = pchar->amax_nphoto * pchar->n_conc_larea / (pchar->kn_nphoto + pchar->n_conc_larea);
		if(pchar->pmax < 0.0){
			pchar->pmax = 0.0;
		}
	}
	
	/* acclimation ************************************************ 2009/04/29 A.Ito */
	if(loct->aco2[grid->m]>400.0 && RAD_SENS==11){
		f_acclim = (1.16 - (loct->aco2[grid->m] - 400.0)*0.00075)/1.16;
	}else{
		f_acclim = 1.0;
	}

	/** light-saturated photosynthesis rate **/
	pchar->psat[grid->m] = pchar->pmax * ftem * fstl * fnstl * f_acclim;
    
    if(EX_OZONE == 1){
        pchar->psat[grid->m] *= fozone;
    }
	
	/* printf("%ld %.2lf %.2lf %.2lf\n", pchar->phototype, ftem, fstl, fnstl);*/
}
