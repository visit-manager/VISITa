/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* Revised August 15, 2007 by A.Ito */
/* Revised August 19, 2007 by A.Ito */

/* initialize environmental characteristics of each grid */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

#define INT_C 0.01

/****************************************************************/
/* Initialize vegetation conditions ****************/
void initVS(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Flux *flux, 
	struct Echar *echar
){	
	long f, g;
	
	/* a priori parameters ***************************/
	/* C3 */
	parameterC3(grid, &(echar->c3));
	/* C4 */
	parameterC4(grid, &(echar->c4));
	/* soil */
	parameterSoil(grid, &(echar->soil));

	/* Sensitivity Analysis ***************************/
	if(SENS==1){
		(echar->c3).gs_b1 *= 1.1;
		(echar->c4).gs_b1*=1.1;
	}else if(SENS==2){
		(echar->c3).pmax*=1.1;
		(echar->c4).pmax*=1.1;
	}else if(SENS==3){
		(echar->c3).eK0*=1.1;
		(echar->c4).eK0*=1.1;
	}else if(SENS==4){
		(echar->c3).lue0*=1.1;
		(echar->c4).lue0*=1.1;
	}else if(SENS==5){
		(echar->c3).albedo*=1.1;
		(echar->c4).albedo*=1.1;
		(echar->soil).albedo0*=1.1;
	}else if(SENS==6){
		/* grid->whc30*=1.1;
		grid->whc*=1.1; */
	}
	
	if(M_D==1){
		(echar->soil).kml*=1.3;
		(echar->soil).kmh*=1.3;
	}else if(M_D==2){
		(echar->soil).kml*=0.7;
		(echar->soil).kmh*=0.7;
	}else if(M_D==3){
		(echar->soil).kmsl*=1.3;
		(echar->soil).kmsh*=1.3;
	}else if(M_D==4){
		(echar->soil).kmsl*=0.7;
		(echar->soil).kmsh*=0.7;
	}
	
	/* growing period **********************/
	for(g=0;g<2;g++){
		for(f=0;f<ASTEP;f++){
			grid->m = f;
			growthperiod(grid, loct, &(echar->c3));
			growthperiod(grid, loct, &(echar->c4));
		}
	}
	
	/* land-use change detritus **********************/
	for(g=0;g<10;g++){
		flux->detr_ten[g] = 0.0;
	}
	
	for(g=0;g<100;g++){
		flux->detr_hund[g] = 0.0;
	}
	
	/* clear vegetation pools and fluxes */
	vlzero(grid, &(mass->c3), &(flux->c3));
	vlzero(grid, &(mass->c4), &(flux->c4));
	vlzero(grid, &(mass->plant), &(flux->plant));
	
	/* set initial biomass *****************************/
	/* Mg C / ha */
	(mass->c3).fol = (mass->c3).fol_p = (mass->c3).stm = (mass->c3).rot = INT_C;
	for(f=0;f<ASTEP;f++){
		(mass->c3).mfol[f] = (mass->c3).mstm[f] = (mass->c3).mrot[f] = INT_C; 
	}

	(mass->c4).fol = (mass->c4).fol_p = (mass->c4).stm= (mass->c4).rot= INT_C;
	for(f=0;f<ASTEP;f++){
		(mass->c4).mfol[f] = (mass->c4).mstm[f] = (mass->c4).mrot[f] = INT_C; 
	}

	(mass->soil).ltr = (mass->soil).msl = INT_C;
	for(f=0;f<ASTEP;f++){
		(mass->soil).ltr_m[f] = (mass->soil).msl_m[f] = INT_C;
	}
	
	/* set initial soil mineral nitrogen */
	/* kg N ha-1 */
	(mass->c3).n_cnpy = (mass->c3).fol/10.0 * 1000.0;
	(mass->c3).n_strg = ((mass->c3).stm + (mass->c3).rot)/10.0 * 1000.0;

	(mass->c4).n_cnpy = (mass->c4).fol/10.0 * 1000.0;
	(mass->c4).n_strg = ((mass->c4).stm + (mass->c4).rot)/10.0 * 1000.0;

	(mass->soil).n_no3 = 10.0 * 1000.0;
	(mass->soil).n_nh4 = 10.0 * 1000.0;

	(mass->soil).n_mcrb = 10.0 * 1000.0;
	(mass->soil).n_lttr = 10.0 * 1000.0;
	(mass->soil).n_hums = 10.0 * 1000.0;
	
	
	for(f=0;f<ASTEP;f++){
		(mass->c3).n_cnpy_m[f] = (mass->c3).n_cnpy;
		(mass->c3).n_strg_m[f] = (mass->c3).n_strg;

		(mass->c4).n_cnpy_m[f] = (mass->c4).n_cnpy;
		(mass->c4).n_strg_m[f] = (mass->c4).n_strg;

		(mass->soil).n_no3_m[f] = (mass->soil).n_no3;
		(mass->soil).n_nh4_m[f] = (mass->soil).n_nh4;

		(mass->soil).n_mcrb_m[f] = (mass->soil).n_mcrb;
		(mass->soil).n_lttr_m[f] = (mass->soil).n_lttr;
		(mass->soil).n_hums_m[f] = (mass->soil).n_hums;
	}
}
