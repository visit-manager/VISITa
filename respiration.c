/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* functions of plant autotrophic respiration */ 

#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* from foliage, maintenance respiration ********************************/
double frfm(
	struct Grid *grid, 
	struct Pchar *veg, 
	struct Pmas *mass
){
	double rfmt0, t0, rfm, ft;
	
	t0 = 15.0; /* specific rate, at 15 deg C */
	
	if(CN_COUPLE >= 1){	
		/* micro mol m-2 s-1 */
		/* g g-1 */
		rfmt0 = veg->n_conc_larea * veg->rd_n * 12.0 / 
				1000000.0 * DHN * HSN / (10000.0/veg->sla * dmTc * 2.0);
	}else{
		/* g g-1 day-1 */
		rfmt0 = veg->rmf/1000.0; 
	}
	
	/* temperature dependence, exponential */
	ft = exp(log(veg->qTf[grid->m]) / 10.0*(grid->tmp_sfc[grid->m] - t0));
	
	if(mass->fol>=0.0){	
		rfm = mass->fol * rfmt0*ft;
	}else{
		rfm = 0.0;
	}
	
	return (rfm);
}

/* from foliage, growth and construction respiration ****************/
double frfg(
	struct Grid *grid, 
	struct Pchar *veg, 
	struct Pflx *flux
){
	double rfg;
	
	/* construction cost is proportional to biomass growth */
	rfg = veg->rgf * flux->tpf[grid->m]; 
	
	return(rfg);
}

/* from stem and branch, maintenance respiration *******************/
double frcm(
	struct Grid *grid, 
	struct Pchar *veg, 
	struct Pmas *mass
){
	double rfmt0, t0, rfm, ft;
	
	t0 = 15.0; /* specific rate, at 15 deg C */
	rfmt0 = veg->rmc / 1000.0; 
	
	/* temperature dependence, exponential */
	ft = exp(log(veg->qTc[grid->m]) / 10.0*(grid->tmp_sfc[grid->m] - t0));
	if(mass->stm >= 0.0){	
		rfm = mass->stm * rfmt0 * ft;
	}else{
		rfm = 0.0;
	}
	
	return(rfm);
}

/* from stem and branch, growth and construction respiration ********/
double frcg(
	struct Grid *grid, 
	struct Pchar *veg, 
	struct Pflx *flux
){
	double rcg;
	
	/* construction cost is proportional to biomass growth */
	rcg = veg->rgc * flux->tpc[grid->m]; 
	
	return(rcg);
}

/* from root, maintenance respiration ******************************/
double frrm(
	struct Grid *grid, 
	struct Pchar *veg, 
	struct Pmas *mass
){
	double rfmt0, t0, rfm, ft;
	t0 = 15.0; /* specific rate, at 15 deg C */
	rfmt0 = veg->rmr/1000.0; 
	
	/* temperature dependence, exponential */
	ft = exp(log(veg->qTr[grid->m]) / 10.0*(grid->tmp_sfc[grid->m] - t0));
	
	if(mass->rot >= 0.0){	
		rfm = mass->rot * rfmt0 * ft;
	}else{
		rfm = 0.0;
	}
	
	return(rfm);
}

/* from root, growth and construction respiration ****************/
double frrg(
	struct Grid *grid, 
	struct Pchar *veg, 
	struct Pflx *flux
){
	double rrg;
	
	/* construction cost is proportional to biomass growth */
	rrg = veg->rgr * flux->tpr[grid->m]; 
	
	return(rrg);
}
