/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/*
	Updated 13 Auguste 2007	
	2008 / 01 / 15, revised	for paddy field
*/

#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* erosion ********************************************************/
void f_erosion(
	struct Grid *grid, 
	struct Loct *loct,  
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux
){
	long f;
	double consv, fveg_a;
	double c_factor_v[16] = {0.0, 
				0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 0.001, 
				0.025, 0.025, 0.10, 0.10, 0.025, 1.0, 1.0
			};
	double p_factor_v[16] = {0.0, 
				1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 
				1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0
			};
	double c_factor_paddy, c_factor_upcrop, f_paddy, f_upcrop, prm_ensen;
	extern double MDN[12];
	
	/* C and P factors for paddy field */
	/* Soil Erosion and Conservation by R.P.C. Morgan
	Morgan, R.P.C., 2005. Soil Erosion and Conservation Third Edition. 
	Blackwell, Oxford, JK, 304 p.
	*/
	c_factor_paddy = 0.15;		/* annual average 0.10 - 0.20 */
	c_factor_upcrop = 0.50;
	
	/* grid->f_crop = f_upcrop + f_paddy */
	if(grid->f_crop_con > 0.0){
		if(grid->f_crop_con >= grid->f_paddy){
			f_paddy = grid->f_paddy;
			f_upcrop = grid->f_crop_con - grid->f_paddy;
		}else{  /*  if(grid->f_paddy > grid->f_crop_con) */
			f_paddy = grid->f_crop_con;
			f_upcrop = 0.0;
		}
	}else{
		f_paddy = 0.0;
		f_upcrop = 0.0;
	}
	
	/* R: rain */
	if(grid->prate_sfc_ann <= 850.0){
		grid->f_erosion_r = 0.0483 * pow(grid->prate_sfc_ann, 1.610);
	}else{
		grid->f_erosion_r = 587.8 - 1.219 * grid->prate_sfc_ann
				+ 0.004105*grid->prate_sfc_ann * grid->prate_sfc_ann;
	}
	if(grid->f_erosion_r < 0.0){
		grid->f_erosion_r = 0.0;
	}
	
	/* K: soil erodibility */
	grid->f_erosion_k = grid->fk_edodibility;
	if(grid->f_erosion_k<0.0){
		grid->f_erosion_k = 0.0;
	}

	/* LS: slope */
	grid->f_erosion_ls = grid->fls_slope;
	if(grid->f_erosion_ls<0.0){
		grid->f_erosion_ls = 0.0;
	}

	/* C: land cover */
	if(PARA_VEGCV == 0){
		/* conventional */
		/* grid->f_erosion_c = (1.0-grid->f_crop)*c_factor_v[grid->veg_sage] + grid->f_crop*0.5; */
		/* grid->f_erosion_c = (1.0-grid->f_crop_con)*c_factor_v[grid->veg_sage] 
			+ f_paddy*c_factor_paddy + f_upcrop * c_factor_upcrop; */
		
		/* case natural: 2011/12/15 (A.Ito) */
		if(loct->v_type == 1){
			grid->f_erosion_c = c_factor_v[grid->veg_sage];
		}
		
		/* case cropland: 2011/12/15 (A.Ito) */
		if(loct->v_type == 2){
			if((f_paddy+f_upcrop) > 0.0){
				grid->f_erosion_c = (f_paddy * c_factor_paddy + f_upcrop * c_factor_upcrop) 
						/ (f_paddy+f_upcrop);
			}else{
				grid->f_erosion_c = 0.0;
			}
		}
		
	}else if(PARA_VEGCV == 1){
		/* parameterization */
		/* fveg_a = 0.0;
		for(f=0;f<12;f++){
			fveg_a += loct->f_vegcov[f] * MDN[f]/YDN;
		}
		
		grid->f_erosion_c = (1.0 - grid->f_crop_con)*c_factor_v[grid->veg_sage]
			* (1.6 - fveg_a) + grid->f_crop_con*0.5; */
		
		fveg_a = 0.0;
		for(f=0;f<12;f++){
			fveg_a += loct->f_vegcov[f] * MDN[f]/YDN;
		}
		
		/* case natural: 2011/12/15 (A.Ito) */
		if(loct->v_type == 1){
			grid->f_erosion_c = c_factor_v[grid->veg_sage] * (1.6 - fveg_a);
		}
		
		/* case cropland: 2011/12/15 (A.Ito) */
		if(loct->v_type == 2){
			grid->f_erosion_c = 0.5;
		}
	}
	if(grid->f_erosion_c < 0.0){
		grid->f_erosion_c = 0.0;
	}

	/* P: protection */
	/* 070725 convensional conservation factor 0.5 */
	/* 070802 revised conservation factor developed - developing */
	if(SOIL_CONSV==0){
		/* Revised: 070802 */
		switch(grid->country){
			/* OECD countries */
			case 840:	consv =	0.75;		break;	/* United States */
			case 826:	consv =	0.75;		break;	/* United Kingdom */
			case 276:	consv =	0.75;		break;	/* Germany */
			case 250:	consv =	0.75;		break;	/* France */
			case 380:	consv =	0.75;		break;	/* Italy */
			case 528:	consv =	0.75;		break;	/* Netherland */
			case 56:	consv =	0.75;		break;	/* Beigium */
			case 442:	consv =	0.75;		break;	/* Luxenburg */
			case 246:	consv =	0.75;		break;	/* Finland */
			case 752:	consv =	0.75;		break;	/* Sweden */
			case 40:	consv =	0.75;		break;	/* Austria */
			case 208:	consv =	0.75;		break;	/* Denmark */
			case 724:	consv =	0.75;		break;	/* Spain */
			case 620:	consv =	0.75;		break;	/* Portugal */
			case 300:	consv =	0.75;		break;	/* Greece */
			case 372:	consv =	0.75;		break;	/* Ireland */
			case 200:	consv =	0.75;		break;	/* Czechoslovakia */
			case 348:	consv =	0.75;		break;	/* Hunagry */
			case 616:	consv =	0.75;		break;	/* Poland */
			case 392:	consv =	0.75;		break;	/* Japan */
			case 124:	consv =	0.75;		break;	/* Canada */
			case 484:	consv =	0.75;		break;	/* Mexico */
			case 36:	consv =	0.75;		break;	/* Australia */
			case 554:	consv =	0.75;		break;	/* New Zealand */
			case 756:	consv =	0.75;		break;	/* Swiss */
			case 578:	consv =	0.75;		break;	/* Norway */
			case 352:	consv =	0.75;		break;	/* Iceland */
			case 792:	consv =	0.75;		break;	/* Turkey */
			case 410:	consv =	0.75;		break;	/* South Korea */

			default:	consv = 0.95;	/** developing countries */
		}
	}else if(SOIL_CONSV == 1){
		/* Conventional assumption : before 070725 */
		consv = 0.5;
	}
	
	/* case natural: 2011/12/15 (A.Ito) */
	if(loct->v_type == 1){
		grid->f_erosion_p = p_factor_v[grid->veg_sage];
	}
	
	/* case cropland: 2011/12/15 (A.Ito) */
	if(loct->v_type == 2){	
		grid->f_erosion_p = consv;
	}

	if(grid->f_erosion_p < 0.0){
		grid->f_erosion_p = 0.0;
	}
    
    /* parameter ensemble: 2014/11/19 by A.Ito */
    prm_ensen = 1.0;
    if(PARAM_PTB == 10){
        if(PARAM_ENS==1){
            prm_ensen *= 0.7;
        }
        if(PARAM_ENS==2){
            prm_ensen *= 0.8;
        }
        if(PARAM_ENS==3){
            prm_ensen *= 0.9;
        }
        if(PARAM_ENS==4){
            prm_ensen *= 1.1;
        }
        if(PARAM_ENS==5){
            prm_ensen *= 1.2;
        }
        if(PARAM_ENS==6){
            prm_ensen *= 1.3;
        }
    }
	
    /* C-budget parameter ensemble: 2018/06/05 by A.Ito */
    if(PARAM_PTB == 20){
        prm_ensen = 1.0 + 0.3 * f_pert[5];
    }
    
	/* Erosion *************************************************************/
	/* case natural: 2011/12/15 (A.Ito) */
	if(loct->v_type == 1){
		flux->erod_soil = grid->f_erosion_r * grid->f_erosion_k * grid->f_erosion_ls * 
				grid->f_erosion_c * grid->f_erosion_p * prm_ensen;  /* t/ha/yr */
		
		if(flux->erod_soil < 0.0){
			flux->erod_soil = 0.0;
		}else if(flux->erod_soil > 130.0*5.0){
			flux->erod_soil = 130.0 * 5.0;
		}
		
		flux->erod_orgmat = flux->erod_soil * grid->pcnt_orgmat / 100.0;
		flux->erod_carbon = flux->erod_orgmat / dmTc;
	}
	
	/* case cropland: 2011/12/15 (A.Ito) */
	if(loct->v_type == 2){	
		flux->erod_soil = grid->f_erosion_r * grid->f_erosion_k * grid->f_erosion_ls * 
					grid->f_erosion_c * grid->f_erosion_p * prm_ensen;  /* t/ha/yr */
		
		if(flux->erod_soil < 0.0){
			flux->erod_soil = 0.0;
		}else if(flux->erod_soil > 130.0*5.0){
			flux->erod_soil = 130.0*5.0;
		}
		
		flux->erod_orgmat = flux->erod_soil * grid->pcnt_orgmat/100.0;
		flux->erod_carbon = flux->erod_orgmat / dmTc;
	}
}
