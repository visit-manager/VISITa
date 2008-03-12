/*	VISIT: Vegetation Integrative SImulation Tool						*/
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

#define TER_CON 0.001 /* criteria for determining the equilibrium, NEP value in Mg C ha-1 yr-1 */

/**************** EQUILIBRIUM ******************/
void cal_stable(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux, 
	FILE *fp_o[6]
){
	long f, nn, term_time;
	double plantmass, ann_nep;
	
	/** maximum simulation times **/
	grid->phase = 0; /* spin-up */
	term_time = 4000;	
	
	(echar->soil).rl = (echar->soil).rl0;
	(echar->soil).rh = (echar->soil).rh0;
	
	f_cult_luc(grid);
	grid->f_crop_p = grid->fcrop_sage[199];
	
	/********** roop to stable stage ***********/
	nn = 0; 
	ann_nep = 10.0;
	while(ann_nep>TER_CON){ /*** acnep>TER_CON nn<10 ***/
		grid->y = nn;
				
		/***** empirical model NPP*****/
		if(grid->y==0){ /* for the first year */
			npp_empirical(grid, loct, flux);
		}
		
		plantmass = ann_nep = 0.0;
		for(f=0;f<12;f++){
			grid->m = f;
			
			/* initialize N fluxes ************/
			n_flux_zero(f, flux);

			cd_trend(grid);
			co2_in_canopy(grid, loct, mass, flux);
			
			/*** environmental condition ***/
			dynmcL(grid, loct, mass, echar);
			
			/***** vegetation processes *****/
			biome_processes(grid, loct, echar, mass, flux);

			/* VOC emission *****************/
			f_voc_emit_guenther97(grid, loct, echar, mass, flux);
			/* Plant CH4 emission *****************/
			f_ch4_emit_veg(grid, loct, echar, mass, flux);

			plant_stand(grid, loct, mass, flux);
			
			(flux->soil).lL[f] = (flux->plant).lL[f];			
			(flux->soil).d13c_lL[f] = (flux->plant).d13c_lL[f];
			
			/***** soil processes *****/
			soil_processes(grid, loct, &(echar->soil), &(mass->soil), &(flux->soil));
			
			/* CH4 oxydation **************/
			f_ch4oxy_ridgewell(grid, loct, flux);
			f_ch4oxy_casa(grid, loct, flux);
			f_ch4oxy_delgrosso(grid, loct, flux);

			/* CH4 emission (wetlands) */
			f_ch4_emit_cao(grid, loct, flux);
			
			/* N2O emission */
			f_n2o_emit_ngas(grid, loct, mass, flux);
			f_n2o_emit_casa(grid, loct, mass, flux);
			
			/***** ecosystem mass balance *****/	
			/* net ecosystem production */
			flux->nep[f] = (flux->plant).npp[f] - (flux->soil).rS[f];
			/* total ecosystem carbon storage */
			mass->total[f] = (mass->c3).plant[f]*loct->C3ptn[f]+(mass->c4).plant[f]*loct->C4ptn[f]+(mass->soil).soil[f];
			/** net carbon balance taking crop harvest into account **/
			flux->ncb[f] = flux->nep[f] + (flux->plant).hvst[f];
			
			/* carbon isotope */
			d13c_efflux(grid, loct, flux);

			/* nitrogen budget */
			n_budget(grid, loct, mass, flux);
						
			/* annual average plant mass */
			plantmass += (mass->plant).plant[f]/12.0;
			
			if(grid->veg_olson!=29&&grid->veg_olson!=30&&grid->veg_olson!=31&&grid->veg_olson!=32){
				ann_nep += flux->nep[f]; /* annula NEP, flux->nep[f]; */
			}else{
				ann_nep += flux->ncb[f]; /* annula NCB, flux->nep[f]; */
			}
			
			/*** acclimation ***/
			(echar->soil).ft0_l[f] = (echar->soil).ft_l[f];
			(echar->soil).ft0_h[f] = (echar->soil).ft_h[f];
			(echar->soil).fm0_l[f] = (echar->soil).fm_l[f];
			(echar->soil).fm0_h[f] = (echar->soil).fm_h[f];
			flux->lL0[f] = (flux->plant).lL[f];
		}
		
		/* biomass burning */
		f_biomassburning(grid, loct, mass, flux);
		
		/* erosion */
		f_erosion(grid, loct, echar, mass, flux);
		
		if(ERSN_CC==1){
			(mass->soil).ltr -= flux->erod_carbon;
			if((mass->soil).ltr<0.0){
				(mass->soil).ltr = 0.0;
			}
		}
		
		/****** terminal conditions ********/
		if(nn<200){	
			/* continued */
			ann_nep = 10.0; 
		}else if(nn>=200 && nn<term_time){	
			ann_nep = fabs(ann_nep); /**** 1. sufficiently stabilized ****/	
		}else if(nn>=term_time){
			break; /**** 3. stop by 2000 years ****/	
		}
		
		if(plantmass<0.0 || plantmass>=500.0){
			printf("!!! BAD plant biomass: %lf\n", plantmass);
			vanish(mass, flux);  /*** 2. excluding abnormal estimates ***/
			break;
		}
		grid->time = nn; /* simulation time of carbon budget */
		nn++;
	}
	/****** end of stabilization roop *******/
	
	/* land use change */
	f_luc_emit(grid, mass, flux);

	/* history data */
	f_set_history_data(0, grid, loct, mass, flux);
		
	/** output initial stable state **/
	publish_cbud(grid, loct, echar, mass, flux, fp_o[0]);
	
	/* output */
	f_output_result(1900, grid, loct, echar, mass, flux, fp_o);
}

