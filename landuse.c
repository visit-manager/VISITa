/*	VISIT: Vegetation Integrative SImulation Tool						*/
/*  Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/*  Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/*  Carbon cycle, erosion, biomass burning, land-use change,			*/
/*  CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* revised:	18/02/2008 by A.Ito */
/* revised:	04/06/2009 by A.Ito based on E.Kato */

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* land use change ************************************************************/
void f_cult_luc(
	struct Grid *grid
){
	if(LANDUSE==0){
		/* Natural condition */
		grid->f_crop_con = 0.0;
		grid->f_pasture_con = 0.0;
	}else if(LANDUSE>=1 && LANDUSE<=5){
		if(grid->climy<=1990){
			/* SAGE, net land use, only cropland */
			grid->f_crop_con = grid->fcrop_sage[grid->climy - 1700];
			grid->f_pasture_con = 0.0;
		}else{
			switch(LANDUSE){
				case 1:
					/* fixed 1900 */
					grid->f_crop_con = grid->fcrop_sage[200]; /*  */
					break;
				case 2:
					/* fixed 1990 */
					grid->f_crop_con = grid->fcrop_sage[290]; /*  */
					break;
				case 3:
					/* conventional scenario : ex b, d */
					grid->f_crop_con = grid->fcrop_sage[290] + grid->f_crop_trend 
										*(double)(grid->climy - 1990);
					break;
				case 4:
					/* high scenario : ex e */
					if(grid->f_crop_trend>0.0){
						grid->f_crop_con = grid->fcrop_sage[290] + grid->f_crop_trend*(1.0+0.01*
								((double)(grid->climy - 1990))) *(double)(grid->climy - 1990);
					}else{
						grid->f_crop_con = grid->fcrop_sage[290];
					}
					break;
				case 5:
					grid->f_crop_con = grid->fcrop_sage[290] + ((grid->fcrop3_image[grid->climy - 1990] 
						+ grid->fcrop4_image[grid->climy - 1990]) 
						- (grid->fcrop3_image[0]+grid->fcrop4_image[0]))/100.0;
					break;
			}
			
			grid->f_pasture_con = 0.0;
		}
	}else if(LANDUSE==6){
		/* EOS-WEBSTER */
		/* Hurtt, G. C., S. Frolking, M. G. Fearon, B. Moore, E. Shevliakova, S. Malyshev, 
			S. W. Pacala, and R. A. Houghton. 2006. The underpinnings of land-use history: 
			three centuries of global gridded land-use transitions, wood-harvest activity, 
			and resulting secondary lands. Global Change Biology 12:1-22. */
		if(grid->climy<=1999){
			grid->f_crop_con = grid->fcrop_eossagehyde[grid->climy - 1700];
			grid->f_pasture_con = grid->fpast_eossagehyde[grid->climy - 1700];
		}else if(grid->climy>=2000){
			grid->f_crop_con = grid->fcrop_eossagehyde[299] + 
				((grid->fcrop3_image[grid->climy - 1990] + grid->fcrop4_image[grid->climy - 1990]) 
				- (grid->fcrop3_image[9]+grid->fcrop4_image[9]))/100.0;
			grid->f_pasture_con = grid->fpast_eossagehyde[299] + 
				((grid->fgrass3_image[grid->climy - 1990] + grid->fgrass4_image[grid->climy - 1990]) 
				- (grid->fgrass3_image[9]+grid->fgrass4_image[9]))/100.0;
		}
	}
	
	/*********************************/
	
	if(grid->f_crop_con < 0.0){
		grid->f_crop_con = 0.0;
	}
	if(grid->f_pasture_con < 0.0){
		grid->f_pasture_con = 0.0;
	}
	
	if((grid->f_crop_con + grid->f_pasture_con) > 0.99){
		grid->f_crop_con = 0.99 / (grid->f_crop_con + grid->f_pasture_con) * grid->f_crop_con;
		grid->f_pasture_con = 0.99 / (grid->f_crop_con + grid->f_pasture_con) * grid->f_pasture_con;
	}

	/*********************************/
	/* annual deforestation */
	if(grid->phase==0){
		/* spin-up */
		grid->f_deforest = grid->fcrop_sage[1901 - 1700] - grid->fcrop_sage[1900 - 1700];
		/* 2008/08/20 corrected by A.Ito (thanks to E.Kato) */
	
	}else{
		if(LANDUSE>=1 && LANDUSE<=5){
			grid->f_deforest = grid->f_crop_con - grid->f_crop_p;
		}else if(LANDUSE==6){
			if(grid->climy<=1999){
				grid->f_deforest = grid->t_vc_eossagehyde[grid->climy - 1700] 
									+ grid->t_vp_eossagehyde[grid->climy - 1700]
									+ grid->t_sc_eossagehyde[grid->climy - 1700] 
									+ grid->t_sp_eossagehyde[grid->climy - 1700];
				grid->f_deforest_v = grid->t_vc_eossagehyde[grid->climy - 1700] 
									+ grid->t_vp_eossagehyde[grid->climy - 1700];
				grid->f_deforest_s = grid->t_sc_eossagehyde[grid->climy - 1700] 
									+ grid->t_sp_eossagehyde[grid->climy - 1700];
			}else{
				grid->f_deforest = (grid->f_crop_con - grid->f_crop_p) 
									+ (grid->f_pasture_con - grid->f_pasture_p);
				grid->f_deforest_v = grid->f_deforest;
				grid->f_deforest_s = 0.0;
			}
		}
	}
	if(grid->f_deforest < 0.0){
		grid->f_deforest = 0.0;
	}
}

/* emission from land use change ****************************************************/
void f_luc_emit(
	struct Grid *grid, 
	struct Mass *mass, 
	struct Flux *flux
){
	long f;
	double f_luc;
	double fe_conv;		/* fraction of conversion flux (1year) */
	double fe_ten;		/* fraction of 10-year pool flux */
	double fe_hund;		/* fraction of 100-year pool flux */
	double fe_detr;		/* fraction of detritus flux */
	double mass_detr, mass_conv, mass_ten, mass_hund;	/* added by A.Ito based on E.Kato (2009/03/30) */
	
	switch(grid->veg_sage){
		/* detritus production by land-use change:
		 McGuire, A. D., S. Sitch, J. S. Clein, R. Dargaville, G. Esser, J. Foley, M. Heimann, F. Joos, 
		 J. Kaplan, D. W. Kicklighter, R. A. Meier, J. M. Melillo, B. I. Moore, L. J. Williams, and 
		 U. Wittenberg, 2001: Carbon balance of the terrestrial biosphere in the twentieth century: analysis 
		 of CO2, climate and land use effects with four process-based ecosystem models. 
		 Global Biogeochemical Cycles, 15, 183-206.
		*/
		/* tropical forests */
		case 1: case 2:  
			fe_detr = 0.33;
			fe_conv = 0.40;
			fe_ten = 0.27;
			fe_hund = 0.0;
			break;
		/* temperate & boreal forests */
		case 3: case 4: case 5: case 6: case 7: case 8:
			fe_detr = 0.33;
			fe_conv = 0.40;
			fe_ten = 0.20;
			fe_hund = 0.07;
			break;
		/* grasslands & tundra */
		case 10: case 13: case 14: case 15:
			fe_detr = 0.50;
			fe_conv = 0.50;
			fe_ten = 0.0;
			fe_hund = 0.0;
			break;
		/* deserts & shrubs */
		case 9: case 11: case 12:
			fe_detr = 0.50;
			fe_conv = 0.40;
			fe_ten = 0.10;
			fe_hund = 0.0;
			break;
		default:
			fe_detr = 0.50;
			fe_conv = 0.40;
			fe_ten = 0.10;
			fe_hund = 0.0;
			break;
	}
	
	if(grid->phase==0){
		/* spin-up: fluxes for 1801-1900 *******************************/
		/* modified by A.Ito (2009/06/05) */
		if(LANDUSE>=1 && LANDUSE<=5){
			f_luc = grid->fcrop_sage[1900-1700] - grid->fcrop_sage[1900-1700-1];
		}else if(LANDUSE==6){
			f_luc = (grid->t_vc_eossagehyde[1900 - 1700] + grid->t_vp_eossagehyde[1900 - 1700])
				+ (grid->t_sc_eossagehyde[1900 - 1700] + grid->t_sp_eossagehyde[1900 - 1700])*0.5;
		}
		
		flux->lu_detr = f_luc * 0.2*(mass->plant).rot;
		flux->lu_conv = f_luc * ((mass->plant).fol + (mass->plant).stm + 0.8*(mass->plant).rot) * fe_conv/(fe_conv+fe_ten+fe_hund);
		
		for(f=1891;f<=1900;f++){
			/* senstivity analysis */
			if(LANDUSE>=1 && LANDUSE<=5){
				f_luc = grid->fcrop_sage[f-1700] - grid->fcrop_sage[f-1700-1];
			}else if(LANDUSE==6){
				f_luc = (grid->t_vc_eossagehyde[f - 1700] + grid->t_vp_eossagehyde[f - 1700])
						+ (grid->t_sc_eossagehyde[f - 1700] + grid->t_sp_eossagehyde[f - 1700])*0.5;
			}
			
			/* modified by A.Ito based on E.Kato (2009/03/30) */
			if(f_luc > 0.0){
				mass_ten = f_luc * ((mass->plant).fol + (mass->plant).stm + 0.8*(mass->plant).rot) * fe_ten/(fe_conv+fe_ten+fe_hund);
				flux->detr_ten[1900-f] = mass_ten;
			}else{
				fe_ten = 0.0;
				flux->detr_ten[1900-f] = 0.0;
			}
			flux->lu_ten = mass_ten*0.1;
		}
		
		for(f=1801;f<=1900;f++){
			/* senstivity analysis */
			if(LANDUSE>=1 && LANDUSE<=5){
				f_luc = grid->fcrop_sage[f-1700] - grid->fcrop_sage[f-1700-1];
			}else if(LANDUSE==6){
				f_luc = (grid->t_vc_eossagehyde[f - 1700] + grid->t_vp_eossagehyde[f - 1700])
						+ (grid->t_sc_eossagehyde[f - 1700] + grid->t_sp_eossagehyde[f - 1700])*0.5;
			}
			
			/* modified by A.Ito based on E.Kato (2009/03/30) */
			if(f_luc > 0.0){
				mass_ten = f_luc * ((mass->plant).fol + (mass->plant).stm + 0.8*(mass->plant).rot) * fe_hund/(fe_conv+fe_ten+fe_hund);				
				flux->detr_hund[1900-f] = mass_ten;
			}else{
				fe_hund = 0.0;
				flux->detr_hund[1900-f] = 0.0;
			}
			flux->lu_hund = fe_hund*0.01;
		}
	}else{
		/* experiment: 1901 - 2000 - 2100 *****************************/
		flux->lu_ten = 0.0;
		/* corrected by A.Ito (2009/06/03) based on E.Kato (2008/11/21) */
		for(f=1;f<10;f++){
			flux->lu_ten += 0.1 * flux->detr_ten[f-1];
		}
		for(f=1;f<10;f++){
			flux->detr_ten[f] = flux->detr_ten[f-1];
		}
		
		flux->lu_hund = 0.0;
		/* corrected by A.Ito (2009/06/03) based on E.Kato (2008/11/21) */
		for(f=1;f<100;f++){
			flux->lu_hund += 0.01 * flux->detr_hund[f-1];
		}
		for(f=1;f<100;f++){
			flux->detr_hund[f] = flux->detr_hund[f-1];
		}
		
		/* annual land use change */
		if(LANDUSE>=1 && LANDUSE<=5){
			f_luc = grid->f_deforest;			/*  grid->f_crop_con - grid->f_crop_p;  */
		}else if(LANDUSE==6){
			f_luc = grid->f_deforest_v + grid->f_deforest_s * 0.5;
		}
		
		if(f_luc > 0.0){ /* deforested */
			/* modified by A.Ito based on E.Kato (2009/03/30) */
			mass_detr = f_luc * 0.2*(mass->plant).rot;
			mass_conv = f_luc * ((mass->plant).fol + (mass->plant).stm + 0.8*(mass->plant).rot) * fe_conv/(fe_conv+fe_ten+fe_hund);
			mass_ten = f_luc * ((mass->plant).fol + (mass->plant).stm + 0.8*(mass->plant).rot) * fe_ten/(fe_conv+fe_ten+fe_hund);
			mass_hund = f_luc * ((mass->plant).fol + (mass->plant).stm + 0.8*(mass->plant).rot) * fe_hund/(fe_conv+fe_ten+fe_hund);
			
			flux->lu_detr = mass_detr;
			flux->lu_conv = mass_conv;
			
			flux->detr_ten[0] = mass_ten;
			flux->lu_ten += 0.1 * mass_ten;

			flux->detr_hund[0] = mass_hund;
			flux->lu_hund += 0.01 * mass_hund;
																		
		}else{
			flux->lu_detr = 0.0;
			flux->lu_conv = 0.0;
			flux->detr_ten[0] = 0.0;
			flux->detr_hund[0] = 0.0;
		}
	}
}

