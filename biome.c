/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed  in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* corrected: A.Ito and E.Kato (2009/08/16) */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* Biome-specific processes **************************************************/
void f_biome_processes(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux
){
	/* C3 ***************************************************************/
	switch(grid->veg_olson){
		case 0: /* WATER */
			vlzero(grid, &(mass->c3), &(flux->c3)); /* no C3 plant */
			break;
		case 1: /* TROPICAL & SUBTROPICAL EVERGREEN FOREST */
			growthperiod(grid, loct, &(echar->c3)); /* processes of C3 plant */
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 2: /* TROPICAL MONTANE FOREST */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 3: /* TROPICAL & SUBTROPICAL DRY FOREST */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 4: /* MID-LATITUDE MIXED FOERST */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 5: /* MID-LATITUDE BROAD-LEAVED FOREST */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 6: /* SEMIARID WOOD OR LOW FOREST */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 7: /* CONIFEROUS EVERGREEN FOREST */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 8: /* SOUTHERN TAIGA */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 9: /* MAIN EVERGREEN TAIGA */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 10: /* MAIN DECIDUOUS TAIGA */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 11: /* NORTHERN EVERGREEN TAIGA */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 12: /* NORTHERN DECIDUOUS TAIGA */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 13: /* SECOND GROWTH WOODS */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 14: /* SECOND GROWTH FIELD */
			growthperiod(grid, loct, &(echar->c3)); /* processes of C3 plant */
			grass_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 15: /* SUCCULENT & THORN WOOD */
			growthperiod(grid, loct, &(echar->c3));
			grass_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 16: /* TROPICAL SAVANNA, WOODLAND */
			growthperiod(grid, loct, &(echar->c3));
			grass_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 17: /* MEDITERRANEAN-TYPE DRY WOOD */
			growthperiod(grid, loct, &(echar->c3));
			grass_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 18: /* HEATH & MOORLAND */
			growthperiod(grid, loct, &(echar->c3));
			grass_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 19: /* WARM OR HOT SHRUB & GRASSLAND */
			growthperiod(grid, loct, &(echar->c3));
			grass_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 20: /* TIBETAN MEADOW & SIBERIAN HIGHLAND */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 21: /* TUNDRA */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 22: /* WOODED TUNDRA */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 23: /* WARM OR HOT WETLANDS */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 24: /* COOL BOG & MIRE */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 25: /* SHORE & HINTERLANDS */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 26: /* COOL SEMI-DESERT SCRUB */
			growthperiod(grid, loct, &(echar->c3));
			forest_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 27: /* NON-POLAR DESERT */
			growthperiod(grid, loct, &(echar->c3));
			grass_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 28: /* NON-POLAR SAND DESERT */
			growthperiod(grid, loct, &(echar->c3));
			grass_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 29: /* PADDYLAND */
			growthperiod(grid, loct, &(echar->c3));
			agri_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 30: /* COOL CROPLAND */
			growthperiod(grid, loct, &(echar->c3));
			agri_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 31: /* WARM CROPLAND */
			growthperiod(grid, loct, &(echar->c3));
			agri_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 32: /* IRRIGATED */
			growthperiod(grid, loct, &(echar->c3));
			agri_process(grid, loct, &(flux->c3), &(echar->c3), &(mass->c3));
			break;
		case 33: /* ANTARCTICA */
			vlzero(grid, &(mass->c3), &(flux->c3));
			break;
	} 

	/* C4 *****************************************************************/
	switch(grid->veg_olson){
		case 0: /* WATER */
			vlzero(grid, &(mass->c4), &(flux->c4)); /* no C4 plant */
			break;
		case 1: /* TROPICAL & SUBTROPICAL EVERGREEN FOREST */
			vlzero(grid, &(mass->c4), &(flux->c4)); /* no C4 plant */
			break;
		case 2: /* TROPICAL MONTANE FOREST */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 3: /* TROPICAL & SUBTROPICAL DRY FOREST */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 4: /* MID-LATITUDE MIXED FOERST */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 5: /* MID-LATITUDE BROAD-LEAVED FOREST */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 6: /* SEMIARID WOOD OR LOW FOREST */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 7: /* CONIFEROUS EVERGREEN FOREST */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 8: /* SOUTHERN TAIGA */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 9: /* MAIN EVERGREEN TAIGA */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 10: /* MAIN DECIDUOUS TAIGA */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 11: /* NORTHERN EVERGREEN TAIGA */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 12: /* NORTHERN DECIDUOUS TAIGA */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 13: /* SECOND GROWTH WOODS */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 14: /* SECOND GROWTH FIELD */
			growthperiod(grid, loct, &(echar->c4)); /* processes of C4 plant */
			grass_process(grid, loct, &(flux->c4), &(echar->c4), &(mass->c4));
			break;
		case 15: /* SUCCULENT & THORN WOOD */
			growthperiod(grid, loct, &(echar->c4));
			grass_process(grid, loct, &(flux->c4), &(echar->c4), &(mass->c4));
			break;
		case 16: /* TROPICAL SAVANNA, WOODLAND */
			growthperiod(grid, loct, &(echar->c4));
			grass_process(grid, loct, &(flux->c4), &(echar->c4), &(mass->c4));
			break;
		case 17: /* MEDITERRANEAN-TYPE DRY WOOD */
			growthperiod(grid, loct, &(echar->c4));
			grass_process(grid, loct, &(flux->c4), &(echar->c4), &(mass->c4));
			break;
		case 18: /* HEATH & MOORLAND */
			growthperiod(grid, loct, &(echar->c4));
			grass_process(grid, loct, &(flux->c4), &(echar->c4), &(mass->c4));
			break;
		case 19: /* WARM OR HOT SHRUB & GRASSLAND */
			growthperiod(grid, loct, &(echar->c4));
			grass_process(grid, loct, &(flux->c4), &(echar->c4), &(mass->c4));
			break;
		case 20: /* TIBETAN MEADOW & SIBERIAN HIGHLAND */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 21: /* TUNDRA */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 22: /* WOODED TUNDRA */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 23: /* WARM OR HOT WETLANDS */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 24: /* COOL BOG & MIRE */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 25: /* SHORE & HINTERLANDS */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 26: /* COOL SEMI-DESERT SCRUB */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 27: /* NON-POLAR DESERT */
			growthperiod(grid, loct, &(echar->c4));
			grass_process(grid, loct, &(flux->c4), &(echar->c4), &(mass->c4));
			break;
		case 28: /* NON-POLAR SAND DESERT */
			growthperiod(grid, loct, &(echar->c4));
			grass_process(grid, loct, &(flux->c4), &(echar->c4), &(mass->c4));
			break;
		case 29: /* PADDYLAND */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
		case 30: /* COOL CROPLAND */
			growthperiod(grid, loct, &(echar->c4));
			agri_process(grid, loct, &(flux->c4), &(echar->c4), &(mass->c4));
			break;
		case 31: /* WARM CROPLAND */
			growthperiod(grid, loct, &(echar->c4));
			agri_process(grid, loct, &(flux->c4), &(echar->c4), &(mass->c4));
			/* vlzero(grid, &(mass->c4), &(flux->c4)); */
			break;
		case 32: /* IRRIGATED */
			growthperiod(grid, loct, &(echar->c4));
			agri_process(grid, loct, &(flux->c4), &(echar->c4), &(mass->c4));
			break;
		case 33: /* ANTARCTICA */
			vlzero(grid, &(mass->c4), &(flux->c4));
			break;
	} 
	
	/* N dynamics *******************************************/
	/* biological N2 fixation */
	f_biolfix(grid, loct, mass, flux);
	
	/* plant root N uptake */
	f_n_uptake(grid, loct, mass, flux);
	
	/* allocation of uptake N within biomass */
	f_n_alloc(grid, loct, &(echar->c3), &(mass->c3), &(flux->c3));
	f_n_alloc(grid, loct, &(echar->c4), &(mass->c4), &(flux->c4));
	
	/* re-allocation of storage N */
	f_n_realloc(grid, loct, &(echar->c3), &(mass->c3), &(flux->c3));
	f_n_realloc(grid, loct, &(echar->c4), &(mass->c4), &(flux->c4));
	
	/* abandon and salvage of aged organ N */
	f_n_abandon_salvage(grid, loct, &(echar->c3), &(mass->c3), &(flux->c3));
	f_n_abandon_salvage(grid, loct, &(echar->c4), &(mass->c4), &(flux->c4));
}
