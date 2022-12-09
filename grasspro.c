/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed  in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* biome dependent plant carbon exchange processes ----- grassland with C3 and C4 plants */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* annual process of grassland ecosystems, where C3 and C4 plants coexist ***/
void grass_process(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pflx *flux, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	/* preparation */
	f_before_deal(grid, flux);
	
	/*** mass change by phenological C3/C4 alternation ***/
	/*	phenol_tran_c34(grid, loct, pchar, mass, flux); */
		
	switch(pchar->season[grid->m]){
		case 0: /* dormancy */
			noleafperiod(grid, loct, pchar, mass, flux); 
			break;
		case 1: /* vegetative growth */
			greenperiod(grid, loct, pchar, mass, flux); 
			reallocation_survival(grid, pchar, mass);
			break;
		case 2: /* leaf emergence */
			leafemergence(grid, loct, pchar, mass, flux); 
			reallocation_survival(grid, pchar, mass);
			break;
		case 3: /* leaf shedding */
			leaffall(grid, loct, pchar, mass, flux); 
			break;
	}
	
	/* post-processing */
	f_after_deal(grid, loct, pchar, mass, flux);
}

