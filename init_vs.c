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
	double p_scale;
	
	/* a priori parameters ***************************/
	/* C3 */
	parameterC3(grid, &(echar->c3));
	/* C4 */
	parameterC4(grid, &(echar->c4));
	/* soil */
	parameterSoil(grid, &(echar->soil));

	/* Sensitivity Analysis ***************************/
	if(SENS_PARA == 1){
		(echar->c3).gs_b1 *= 1.1;
		(echar->c4).gs_b1 *= 1.1;
	}else if(SENS_PARA == 2){
		(echar->c3).pmax *= 1.1;
		(echar->c4).pmax *= 1.1;
	}else if(SENS_PARA == 3){
		(echar->c3).eK0 *= 1.1;
		(echar->c4).eK0 *= 1.1;
	}else if(SENS_PARA == 4){
		(echar->c3).lue0 *= 1.1;
		(echar->c4).lue0 *= 1.1;
	}else if(SENS_PARA == 5){
		(echar->c3).albedo *= 1.1;
		(echar->c4).albedo *= 1.1;
		(echar->soil).albedo0 *= 1.1;
	}else if(SENS_PARA == 6){
		/* grid->whc30*=1.1;
		grid->whc*=1.1; */
	}
    
    /* albedo perturbation: 2012/12/29 by A.Ito ********/
    if(EX_ALBEDO == 1){
        p_scale = 0.0;
        for(f=0;f<12;f++){
            p_scale += (double)rand() / (double)RAND_MAX;
        }
        p_scale -= 6.0;
        
        if(p_scale > 3.0){
            p_scale = 3.0;
        }
        if(p_scale < -3.0){
            p_scale = -3.0;
        }
        
        for(f=0;f<ASTEP;f++){
            grid->albedo_pert[f] = p_scale * grid->albedo_sd[f][grid->row/10][grid->col/10];
        }
    }else{
        for(f=0;f<ASTEP;f++){
            grid->albedo_pert[f] = 0.0;
        }
    }
	
	if(M_D == 1){
		(echar->soil).kml *= 1.3;
		(echar->soil).kmh *= 1.3;
	}else if(M_D == 2){
		(echar->soil).kml *= 0.7;
		(echar->soil).kmh *= 0.7;
	}else if(M_D == 3){
		(echar->soil).kmsl *= 1.3;
		(echar->soil).kmsh *= 1.3;
	}else if(M_D == 4){
		(echar->soil).kmsl *= 0.7;
		(echar->soil).kmsh *= 0.7;
	}
	
	/* parameter perturbation: 2010/05/10 added by A.Ito */
	p_scale = 0.06667;
    if(PARAM_PTB == 1){
        /**/
        (echar->c3).pmax *= 1.0 + p_scale*f_pert[0];
        (echar->c4).pmax *= 1.0 + p_scale*f_pert[0];
        /**/
        (echar->c3).lue0 *= 1.0 + p_scale*f_pert[1];
        (echar->c4).lue0 *= 1.0 + p_scale*f_pert[1];
        /**/
        (echar->c3).sla *= 1.0 + p_scale*f_pert[2];
        (echar->c4).sla *= 1.0 + p_scale*f_pert[2];
        /**/
        (echar->c3).qTf0 *= 1.0 + p_scale*f_pert[3];
        (echar->c4).qTf0 *= 1.0 + p_scale*f_pert[3];
        (echar->c3).qTc0 *= 1.0 + p_scale*f_pert[3];
        (echar->c4).qTc0 *= 1.0 + p_scale*f_pert[3];
        (echar->c3).qTr0 *= 1.0 + p_scale*f_pert[3];
        (echar->c4).qTr0 *= 1.0 + p_scale*f_pert[3];
        /**/
        (echar->c3).lf0 *= 1.0 + p_scale*f_pert[4];
        (echar->c4).lf0 *= 1.0 + p_scale*f_pert[4];
        /**/
        (echar->c3).topt0 += 0.3333*f_pert[5];
        (echar->c4).topt0 += 0.3333*f_pert[5];
        /**/
        (echar->c3).tmin += 0.3333*f_pert[6];
        (echar->c4).tmin += 0.3333*f_pert[6];
        /**/
        (echar->c3).kmci *= 1.0 + p_scale*f_pert[7];
        (echar->c4).kmci *= 1.0 + p_scale*f_pert[7];
        /**/
        (echar->c3).km_nstl *= 1.0 + p_scale*f_pert[8];
        (echar->c4).km_nstl *= 1.0 + p_scale*f_pert[8];
        
        /**/
        (echar->soil).rl0 *= 1.0 + p_scale*f_pert[9];
        /**/
        (echar->soil).rh0 *= 1.0 + p_scale*f_pert[10];
    }
    
    /* parameter ensemble: 2014/11/19 by A.Ito */
    if(PARAM_PTB == 2){
        p_scale = 0.1;
        /**/
        (echar->c3).pmax *= 1.0 + p_scale*f_pert[0];
        (echar->c4).pmax *= 1.0 + p_scale*f_pert[0];
        /**/
        (echar->c3).lue0 *= 1.0 + p_scale*f_pert[1];
        (echar->c4).lue0 *= 1.0 + p_scale*f_pert[1];
        /**/
        (echar->c3).sla *= 1.0 + p_scale*f_pert[2];
        (echar->c4).sla *= 1.0 + p_scale*f_pert[2];
        /**/
        (echar->c3).topt0 += 0.3333*f_pert[3];
        (echar->c4).topt0 += 0.3333*f_pert[3];
        /**/
        (echar->c3).tmin += 0.3333*f_pert[4];
        (echar->c4).tmin += 0.3333*f_pert[4];
        /**/
        (echar->c3).kmci *= 1.0 + p_scale*f_pert[5];
        (echar->c4).kmci *= 1.0 + p_scale*f_pert[5];
        /**/
        (echar->c3).km_nstl *= 1.0 + p_scale*f_pert[6];
        (echar->c4).km_nstl *= 1.0 + p_scale*f_pert[6];
    }
    
    if(PARAM_PTB == 3){
        p_scale = 0.1;
        /**/
        (echar->c3).qTf0 *= 1.0 + p_scale*f_pert[7];
        (echar->c4).qTf0 *= 1.0 + p_scale*f_pert[7];
        (echar->c3).qTc0 *= 1.0 + p_scale*f_pert[7];
        (echar->c4).qTc0 *= 1.0 + p_scale*f_pert[7];
        (echar->c3).qTr0 *= 1.0 + p_scale*f_pert[7];
        (echar->c4).qTr0 *= 1.0 + p_scale*f_pert[7];
        /**/
        (echar->c3).rmf0 *= 1.0 + p_scale*f_pert[8];
        (echar->c4).rmf0 *= 1.0 + p_scale*f_pert[8];
        (echar->c3).rmc_s *= 1.0 + p_scale*f_pert[8];
        (echar->c4).rmc_s *= 1.0 + p_scale*f_pert[8];
        (echar->c3).rmc_h *= 1.0 + p_scale*f_pert[8];
        (echar->c4).rmc_h *= 1.0 + p_scale*f_pert[8];
        (echar->c3).rmr_s *= 1.0 + p_scale*f_pert[8];
        (echar->c4).rmr_s *= 1.0 + p_scale*f_pert[8];
        (echar->c3).rmr_h *= 1.0 + p_scale*f_pert[8];
        (echar->c4).rmr_h *= 1.0 + p_scale*f_pert[8];
        /**/
        (echar->c3).rgf *= 1.0 + p_scale*f_pert[9];
        (echar->c4).rgf *= 1.0 + p_scale*f_pert[9];
        (echar->c3).rgc *= 1.0 + p_scale*f_pert[9];
        (echar->c4).rgc *= 1.0 + p_scale*f_pert[9];
        (echar->c3).rgr *= 1.0 + p_scale*f_pert[9];
        (echar->c4).rgr *= 1.0 + p_scale*f_pert[9];
    }
    
    if(PARAM_PTB == 4){
        p_scale = 0.1;
        /**/
        (echar->soil).rl0 *= 1.0 + p_scale*f_pert[10];
        /**/
        (echar->soil).rh0 *= 1.0 + p_scale*f_pert[11];
        /**/
        (echar->soil).kml *= 1.0 + p_scale*f_pert[12];
        (echar->soil).kmh *= 1.0 + p_scale*f_pert[12];
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
	/* g N ha-1 */
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
