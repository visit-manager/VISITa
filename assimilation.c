/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/** This file gives description of functions releted to photosynthesis, or GPP **/
#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* formula of daily gross primary production *******************/
double fgpp(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *veg, 
	struct Pmas *mass
){
	double gpp;
	double cc1, cc2, cc3, bb;
	
	/* give daily GPP based on the Monsi-Saeki theory, or Kuroiwa's equation */
	/*
	Kuroiwa, S., 1966. Dry matter production of plants. Ecology and Evolution. 
	Iwanami Shoten, Tokyo, pp. 71-100.
	*/
	if(veg->psat[grid->m] > 0.0){
		cc1 = 2.0*veg->psat[grid->m]*grid->dlen[grid->m]*lTs/veg->eK[grid->m]; 
		bb = veg->eK[grid->m]*veg->lue[grid->m]*grid->par[grid->m]/veg->psat[grid->m];
		cc2 = 1.0 + sqrt(1.0+bb);
		cc3 = 1.0 + sqrt(1.0+bb*exp(-1.0*veg->eK[grid->m]*mass->lai[grid->m]));
		gpp = cc1*log(cc2/cc3);
	}else{
		gpp = 0.0;
	}
	
	return (gpp);
}

/* light-saturated photosynthetic rate *************************/
void pc_sat(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *veg
){
	double ftem, fstl, fnstl;
	double aa1, aa2, aa3;
	
	/** optimum photosynthesis temperature **/
	if(veg->phototype==3){ /* C3 plants with change */
		veg->topt = veg->topt0 + 0.01*veg->ci[grid->m];	
	}else if(veg->phototype==4){ /* C4 plants without change */
		veg->topt = veg->topt0;
	}
	
	/** CO2 compensation point **/
	if(veg->phototype==3){ /* C3 plants with change: Brooks&Farquhar(1985) */
		aa3 = 0.000347*(grid->tmp_sfc[grid->m]-20.0)*(grid->tmp_sfc[grid->m]-20.0);
		aa1 = 1.0+0.0451*(grid->tmp_sfc[grid->m]-20.0)+aa3;
		aa1 = (aa1>0.0)?aa1:0.0;
		veg->cmpcd[grid->m] = veg->cmpcd0*aa1;	
	}else if(veg->phototype==4){ /* C4 plants without change */
		veg->cmpcd[grid->m] = veg->cmpcd0;
	}
		
	/** temperature effect ******************************************************/
	if(grid->tmp_sfc[grid->m]<=veg->tmax && grid->tmp_sfc[grid->m]>=veg->tmin){
		aa1 = (grid->tmp_sfc[grid->m]-veg->tmax)*(grid->tmp_sfc[grid->m]-veg->tmin);
		aa2 = (grid->tmp_sfc[grid->m]-veg->topt)*(grid->tmp_sfc[grid->m]-veg->topt);
		ftem = aa1/(aa1 - aa2);
	}else{
		ftem = 0.0;
	}
	ftem = (ftem<=1.0)?ftem:1.0; 
	ftem = (ftem>=0.0)?ftem:0.0;
	
	veg->ft[grid->m] = ftem;
	
	/** CO2 effect ******************************************************/
	/** stomatal limitation via intercellular CO2 concentration **/
	if(veg->phototype==3){ /* C3 plants */
		fstl = 0.05 + 0.95*(veg->ci[grid->m] - veg->cmpcd[grid->m])/(veg->kmci + veg->ci[grid->m]); 
	}else if(veg->phototype==4){ /* C4 plants */
		fstl = 0.60 + 0.40*(veg->ci[grid->m] - veg->cmpcd[grid->m])/(veg->kmci + veg->ci[grid->m]); 
	}
	fstl = (fstl<=1.0)?fstl:1.0; 
	fstl = (fstl>=0.0)?fstl:0.0;
	
	veg->fcd[grid->m] = fstl;

	/** soil water effect ******************************************************/
	/** non-stomatal limitation **/
	if(veg->phototype==3){ /* C3 plants */
		fnstl = 0.95*loct->sww / (loct->sww + grid->field_cap2*veg->km_nstl) + 0.05;
	}else if(veg->phototype==4){ /* C4 plants */
		fnstl = 0.86*loct->sww / (loct->sww + grid->field_cap2*veg->km_nstl) + 0.14;
	}
	fnstl = (fnstl<=1.0)?fnstl:1.0; fnstl=(fnstl>=0.0)?fnstl:0.0;
		
	veg->fsw[grid->m] = fnstl;
		
	/* leaf N effect ***********************************************/
	if(CN_COUPLE >= 1){		
		veg->pmax = veg->amax_nphoto * veg->n_conc_larea / (veg->kn_nphoto + veg->n_conc_larea);
		if(veg->pmax < 0.0){
			veg->pmax = 0.0;
		}
	}

	/** light-saturated photosynthesis rate **/
	veg->psat[grid->m] = veg->pmax*ftem*fstl*fnstl; 
	
	/* printf("%ld %.2lf %.2lf %.2lf\n", veg->phototype, ftem, fstl, fnstl);*/
}
