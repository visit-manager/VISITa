/*	VISIT: Vegetation Integrative SImulation Tool						*/
/*  Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
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
	if(LANDUSE == 0){
		/* Natural condition *************/
		grid->f_crop_con = 0.0;
		grid->f_pasture_con = 0.0;
	}else if(LANDUSE>=1 && LANDUSE<=5){
		/* SAGE land-use data:
		 Ramankutty, N., and J. A. Foley (1999), Estimating historical changes in global 
		 land cover: croplands from 1700 to 1992, Global Biogeochemical Cycles, 13(4), 997-1027.
		 */
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
					/* IMAGE-based scenario */
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
			grid->f_crop_con = grid->fcrop_unh_hmnzed[grid->climy - 1700];
			grid->f_pasture_con = grid->fpast_unh_hmnzed[grid->climy - 1700];
		}else if(grid->climy>=2000){
			grid->f_crop_con = grid->fcrop_unh_hmnzed[299] + 
				((grid->fcrop3_image[grid->climy - 1990] + grid->fcrop4_image[grid->climy - 1990]) 
				- (grid->fcrop3_image[9]+grid->fcrop4_image[9]))/100.0;
			grid->f_pasture_con = grid->fpast_unh_hmnzed[299] + 
				((grid->fgrass3_image[grid->climy - 1990] + grid->fgrass4_image[grid->climy - 1990]) 
				- (grid->fgrass3_image[9]+grid->fgrass4_image[9]))/100.0;
		}
	}else if(LANDUSE==7){
		/* SAGE land-use data 1700-2007 (Revised: 2010/01/07):
		 Ramankutty, N., and J. A. Foley (1999), Estimating historical changes in global 
		 land cover: croplands from 1700 to 1992, Global Biogeochemical Cycles, 13(4), 997-1027.
		*/
		if(grid->climy<=2007){
			grid->f_crop_con = grid->fcrop_rk[grid->climy - 1700];
			grid->f_pasture_con = grid->fpast_rk[grid->climy - 1700];
		}else if(grid->climy>=2008){
			grid->f_crop_con = grid->fcrop_rk[307] + 
					((grid->fcrop3_image[grid->climy - 1990] + grid->fcrop4_image[grid->climy - 1990]) 
					- (grid->fcrop3_image[17] + grid->fcrop4_image[17]))/100.0;
			grid->f_pasture_con = grid->fpast_rk[307] + 
					((grid->fgrass3_image[grid->climy - 1990] + grid->fgrass4_image[grid->climy - 1990]) 
					 - (grid->fgrass3_image[17] + grid->fgrass4_image[17]))/100.0;
		}
	}else if(LANDUSE==8){
		/* Hurtt harmonized land-use change, 1700-2005 (added 2010/01/31) */
		if(grid->climy<=2005){
			grid->f_crop_con = grid->fcrop_unh_hmnzed[grid->climy - 1700];
			grid->f_pasture_con = grid->fpast_unh_hmnzed[grid->climy - 1700];
		}else if(grid->climy>=2006){
			grid->f_crop_con = grid->fcrop_unh_hmnzed[305] + 
				((grid->fcrop3_image[grid->climy - 1990] + grid->fcrop4_image[grid->climy - 1990]) 
				 - (grid->fcrop3_image[15]+grid->fcrop4_image[15]))/100.0;
			grid->f_pasture_con = grid->fpast_unh_hmnzed[305] + 
				((grid->fgrass3_image[grid->climy - 1990] + grid->fgrass4_image[grid->climy - 1990]) 
				 - (grid->fgrass3_image[15]+grid->fgrass4_image[15]))/100.0;
		}
	}else{
		printf("Wrong land-use setting ID\n");
		exit(1);
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
		if(LANDUSE == 0){
			grid->f_deforest = 0.0;
		}else if(LANDUSE>=1 && LANDUSE<=5){
			grid->f_deforest = grid->f_crop_con - grid->f_crop_p;
		}else if(LANDUSE==6 || LANDUSE==8){
			if(grid->climy<=1999){
				grid->f_deforest = grid->t_vc_unh_hmnzed[grid->climy - 1700] 
									+ grid->t_vp_unh_hmnzed[grid->climy - 1700]
									+ grid->t_sc_unh_hmnzed[grid->climy - 1700] 
									+ grid->t_sp_unh_hmnzed[grid->climy - 1700];
				grid->f_deforest_v = grid->t_vc_unh_hmnzed[grid->climy - 1700] 
									+ grid->t_vp_unh_hmnzed[grid->climy - 1700];
				grid->f_deforest_s = grid->t_sc_unh_hmnzed[grid->climy - 1700] 
									+ grid->t_sp_unh_hmnzed[grid->climy - 1700];
			}else{
				grid->f_deforest = (grid->f_crop_con - grid->f_crop_p) 
									+ (grid->f_pasture_con - grid->f_pasture_p);
				grid->f_deforest_v = grid->f_deforest;
				grid->f_deforest_s = 0.0;
			}
		}else if(LANDUSE==7){
			grid->f_deforest = (grid->f_crop_con - grid->f_crop_p) 
								+ (grid->f_pasture_con - grid->f_pasture_p);
		}
	}
	
	/* abandonment */
	if(grid->f_deforest < 0.0){
		grid->f_deforest = 0.0;
	}
	
	/* historical change in paddy field area: added by A.Ito (2011/2/28) ********/
	if(LANDUSE==7){
		if(grid->climy<=2005){
			if(grid->f_paddy_b > 0.0 && grid->fcrop_rk[2000 - 1700] > 0.0){
				grid->f_paddy = grid->f_paddy_b * 
					(grid->fcrop_rk[grid->climy - 1700] / grid->fcrop_rk[2000 - 1700]);
			}
		}else if(grid->climy>=2006){
			if(grid->f_paddy_b > 0.0 && grid->fcrop_rk[2000 - 1700] > 0.0){
				grid->f_paddy = grid->f_paddy_b * 
					(grid->fcrop_rk[2005 - 1700] / grid->fcrop_rk[2000 - 1700]);
			}
		}
		
		if(grid->f_paddy > 1.0){
			grid->f_paddy = 1.0;
		}
		if(grid->f_paddy < 0.0){
			grid->f_paddy = 0.0;
		}
	}if(LANDUSE==8){
		if(grid->climy<=2005){
			if(grid->f_paddy_b > 0.0 && grid->fcrop_unh_hmnzed[2000 - 1700] > 0.0){
				grid->f_paddy = grid->f_paddy_b * 
				(grid->fcrop_unh_hmnzed[grid->climy - 1700] / grid->fcrop_unh_hmnzed[2000 - 1700]);
			}
		}else if(grid->climy>=2006){
			if(grid->f_paddy_b > 0.0 && grid->fcrop_unh_hmnzed[2000 - 1700] > 0.0){
				grid->f_paddy = grid->f_paddy_b * 
				(grid->fcrop_unh_hmnzed[2005 - 1700] / grid->fcrop_unh_hmnzed[2000 - 1700]);
			}
		}
		
		if(grid->f_paddy > 1.0){
			grid->f_paddy = 1.0;
		}
		if(grid->f_paddy < 0.0){
			grid->f_paddy = 0.0;
		}
	}else if(LANDUSE==0){
		grid->f_paddy = 0.0;
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
	double eff_mass, f_mass_secfor;
	
	/* mass fraction of secondary forest */
	/* Nelson et al. (2000) Bioscience 50:419-431 */
	f_mass_secfor = 0.1;
	
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
	
	/* effectice biomass */
	eff_mass = (mass->plant).fol + (mass->plant).stm + 0.8*(mass->plant).rot;
	
	if(grid->phase==0){
		/* spin-up: fluxes for 1801-1900 *******************************/
		/* modified by A.Ito (2009/06/05: 2010/01/07) */
		if(LANDUSE == 0){
			f_luc = 0.0;
		}else if(LANDUSE>=1 && LANDUSE<=5){
			f_luc = grid->fcrop_sage[1900-1700] - grid->fcrop_sage[1900-1700-1];
		}else if(LANDUSE==6){
			f_luc = (grid->t_vc_unh_hmnzed[1900 - 1700] + grid->t_vp_unh_hmnzed[1900 - 1700])
				+ (grid->t_sc_unh_hmnzed[1900 - 1700] + grid->t_sp_unh_hmnzed[1900 - 1700])*f_mass_secfor;
		}else if(LANDUSE==7){
			/* added 2010/01/07 (A.Ito) */
			f_luc = (grid->fcrop_rk[1900-1700] - grid->fcrop_rk[1900-1700-1])
					+(grid->fpast_rk[1900-1700] - grid->fpast_rk[1900-1700-1]);
		}
		
		/* modified by A.Ito (2009/08/19) */
		if(f_luc > 0.0){
			flux->lu_detr = f_luc * 0.2*(mass->plant).rot;
			flux->lu_conv = f_luc * eff_mass * fe_conv/(fe_conv + fe_ten + fe_hund);
		}else{
			flux->lu_detr = 0.0;
			flux->lu_conv = 0.0;
		}
		
		for(f=1891;f<=1900;f++){
			/* senstivity analysis */
			if(LANDUSE == 0){
				f_luc = 0.0;
			}else if(LANDUSE>=1 && LANDUSE<=5){
				f_luc = grid->fcrop_sage[f-1700] - grid->fcrop_sage[f-1700-1];
			}else if(LANDUSE==6 || LANDUSE==8){
				f_luc = (grid->t_vc_unh_hmnzed[f - 1700] + grid->t_vp_unh_hmnzed[f - 1700])
						+ (grid->t_sc_unh_hmnzed[f - 1700] + grid->t_sp_unh_hmnzed[f - 1700])*f_mass_secfor;
				/* 0.5: assumption by A.Ito for secondary forest stock */
			}else if(LANDUSE==7){
				/* added 2010/01/07 (A.Ito) */
				f_luc = (grid->fcrop_rk[f-1700] - grid->fcrop_rk[f-1700-1])
						+ (grid->fpast_rk[f-1700] - grid->fpast_rk[f-1700-1]);
			}
			
			/* modified by A.Ito based on E.Kato (2009/03/30) */
			if(f_luc > 0.0){
				mass_ten = f_luc * eff_mass * fe_ten/(fe_conv + fe_ten + fe_hund);
				flux->detr_ten[1900-f] = mass_ten;
			}else{
				fe_ten = 0.0;
				flux->detr_ten[1900-f] = 0.0;
			}
			/* corrected: A.Ito and E.Kato (2009/08/16) */
			flux->lu_ten += mass_ten*0.1;
		}
		
		for(f=1801;f<=1900;f++){
			/* senstivity analysis */
			if(LANDUSE == 0){
				f_luc = 0.0;
			}else if(LANDUSE>=1 && LANDUSE<=5){
				f_luc = grid->fcrop_sage[f-1700] - grid->fcrop_sage[f-1700-1];
			}else if(LANDUSE==6 || LANDUSE==8){
				f_luc = (grid->t_vc_unh_hmnzed[f - 1700] + grid->t_vp_unh_hmnzed[f - 1700])
						+ (grid->t_sc_unh_hmnzed[f - 1700] + grid->t_sp_unh_hmnzed[f - 1700])*f_mass_secfor;
			}else if(LANDUSE==7){
				/* added 2010/01/07 (A.Ito) */
				f_luc = (grid->fcrop_rk[f-1700] - grid->fcrop_rk[f-1700-1])
						+ (grid->fpast_rk[f-1700] - grid->fpast_rk[f-1700-1]);
			}
			
			/* modified by A.Ito based on E.Kato (2009/03/30) */
			/* corrected: A.Ito and E.Kato (2009/08/16) */
			if(f_luc > 0.0){
				mass_hund = f_luc * eff_mass * fe_hund/(fe_conv + fe_ten + fe_hund);				
				flux->detr_hund[1900-f] = mass_hund;
			}else{
				mass_hund = 0.0;
				flux->detr_hund[1900-f] = 0.0;
			}
			flux->lu_hund += mass_hund*0.01;
		}
	}else{
		/* experiment: 1901 - 2000 - 2100 *****************************/
		
		/* emission from 10-year pool */
		/* corrected by A.Ito (2009/06/03) based on E.Kato (2008/11/21) */
		/* corrected: A.Ito and E.Kato (2009/08/16) */
		flux->lu_ten = 0.0;
		flux->lu_ten += 0.1 * flux->detr_ten[0];
		for(f=9;f>0;f--){
			flux->lu_ten += 0.1 * flux->detr_ten[f];
			flux->detr_ten[f] = flux->detr_ten[f-1];
		}
		
		/* emission from 100-year pool */
		/* corrected by A.Ito (2009/06/03) based on E.Kato (2008/11/21) */
		/* corrected: A.Ito and E.Kato (2009/08/16) */
		flux->lu_hund = 0.0;
		flux->lu_hund += 0.01 * flux->detr_hund[0];
		for(f=99;f>0;f--){
			flux->lu_hund += 0.01 * flux->detr_hund[f];
			flux->detr_hund[f] = flux->detr_hund[f-1];
		}
		
		/* annual land use change */
		if(LANDUSE == 0){
			f_luc = 0.0;
		}else if(LANDUSE>=1 && LANDUSE<=5){
			f_luc = grid->f_deforest;
			/*  grid->f_crop_con - grid->f_crop_p;  */
		}else if(LANDUSE==6 || LANDUSE==8){
			/* assumption: biomass in secondary forest is half (0.5) of primary forest */
			f_luc = grid->f_deforest_v + grid->f_deforest_s * f_mass_secfor;
		}else if(LANDUSE==7){
			/* added 2010/01/07 (A.Ito) */
			f_luc = grid->f_deforest;
		}
		
		if(f_luc > 0.0){ /* deforested */
			/* modified by A.Ito based on E.Kato (2009/03/30) */
			mass_detr = f_luc * 0.2*(mass->plant).rot;
			mass_conv = f_luc * eff_mass * fe_conv/(fe_conv + fe_ten + fe_hund);
			mass_ten = f_luc * eff_mass * fe_ten/(fe_conv + fe_ten + fe_hund);
			mass_hund = f_luc * eff_mass * fe_hund/(fe_conv + fe_ten + fe_hund);
			
			/* emission from 1-yr or instantaneous pool */
			flux->lu_detr = mass_detr;
			flux->lu_conv = mass_conv;
			
			/* corrected: A.Ito and E.Kato (2009/08/16) */
			flux->detr_ten[0] = mass_ten;
			flux->detr_hund[0] = mass_hund;
		}else{
			flux->lu_detr = 0.0;
			flux->lu_conv = 0.0;
			flux->detr_ten[0] = 0.0;
			flux->detr_hund[0] = 0.0;
		}
	}
}
