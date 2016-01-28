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
    double prm_ensen;

	if(LANDUSE == 0){
		/* Natural condition *************/
		grid->f_crop_con = 0.0;
		grid->f_pasture_con = 0.0;
	}else if(LANDUSE==1 || LANDUSE==2 || LANDUSE==3 || LANDUSE==4 || LANDUSE==5){
		/* SAGE land-use data:
		 Ramankutty, N., and J. A. Foley (1999), Estimating historical changes in global 
		 land cover: croplands from 1700 to 1992, Global Biogeochemical Cycles, 13(4), 997-1027.
        */
		if(grid->lucy <= 1990){
			/* SAGE, net land use, only cropland */
			grid->f_crop_con = grid->fcrop_sage[grid->lucy - PIVOT_LUC];
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
										*(double)(grid->lucy - 1990);
					break;
				case 4:
					/* high scenario : ex e */
					if(grid->f_crop_trend > 0.0){
						grid->f_crop_con = grid->fcrop_sage[290] + grid->f_crop_trend*(1.0+0.01*
								((double)(grid->lucy - 1990))) *(double)(grid->lucy - 1990);
					}else{
						grid->f_crop_con = grid->fcrop_sage[290];
					}
					break;
				case 5:
					/* IMAGE-based scenario */
					grid->f_crop_con = grid->fcrop_sage[290] + ((grid->fcrop3_image[grid->lucy - 1990]
						+ grid->fcrop4_image[grid->lucy - 1990])
						- (grid->fcrop3_image[0] + grid->fcrop4_image[0]))/100.0;
					break;
                default:
                    break;
			}
			
			grid->f_pasture_con = 0.0;
		}
	}else if(LANDUSE == 6){
		/* EOS-WEBSTER */
		/* Hurtt, G. C., S. Frolking, M. G. Fearon, B. Moore, E. Shevliakova, S. Malyshev, 
			S. W. Pacala, and R. A. Houghton. 2006. The underpinnings of land-use history: 
			three centuries of global gridded land-use transitions, wood-harvest activity, 
			and resulting secondary lands. Global Change Biology 12:1-22. */
		if(grid->lucy <= 1999){
			grid->f_crop_con = grid->fcrop_unh_hmnzed[grid->lucy - PIVOT_LUC];
			grid->f_pasture_con = grid->fpast_unh_hmnzed[grid->lucy - PIVOT_LUC];
		}else if(grid->lucy >= 2000){
			grid->f_crop_con = grid->fcrop_unh_hmnzed[299] + 
				((grid->fcrop3_image[grid->lucy - 1990] + grid->fcrop4_image[grid->lucy - 1990])
				- (grid->fcrop3_image[9]+grid->fcrop4_image[9]))/100.0;
			grid->f_pasture_con = grid->fpast_unh_hmnzed[299] + 
				((grid->fgrass3_image[grid->lucy - 1990] + grid->fgrass4_image[grid->lucy - 1990])
				- (grid->fgrass3_image[9] + grid->fgrass4_image[9]))/100.0;
		}
	}else if(LANDUSE == 7){
		/* SAGE land-use data 1700-2007 (Revised: 2010/01/07):
		 Ramankutty, N., and J. A. Foley (1999), Estimating historical changes in global 
		 land cover: croplands from 1700 to 1992, Global Biogeochemical Cycles, 13(4), 997-1027.
		*/
		if(grid->lucy <= 2007){
			grid->f_crop_con = grid->fcrop_rk[grid->lucy - PIVOT_LUC];
			grid->f_pasture_con = grid->fpast_rk[grid->lucy - PIVOT_LUC];
		}else if(grid->lucy >= 2008){
			grid->f_crop_con = grid->fcrop_rk[307] + 
					((grid->fcrop3_image[grid->lucy - 1990] + grid->fcrop4_image[grid->lucy - 1990])
					- (grid->fcrop3_image[17] + grid->fcrop4_image[17]))/100.0;
			grid->f_pasture_con = grid->fpast_rk[307] + 
					((grid->fgrass3_image[grid->lucy - 1990] + grid->fgrass4_image[grid->lucy - 1990])
					 - (grid->fgrass3_image[17] + grid->fgrass4_image[17]))/100.0;
		}
	}else if(LANDUSE == 8){
		/* UNH harmonized land-use change, 1700-2005 (added 2010/01/31) */
		if(grid->lucy <= 2005){
			grid->f_crop_con = grid->fcrop_unh_hmnzed[grid->lucy - PIVOT_LUC];
			grid->f_pasture_con = grid->fpast_unh_hmnzed[grid->lucy - PIVOT_LUC];
		}else if(grid->lucy >= 2006){
			grid->f_crop_con = grid->fcrop_unh_hmnzed[305] + 
				((grid->fcrop3_image[grid->lucy - 1990] + grid->fcrop4_image[grid->lucy - 1990])
				 - (grid->fcrop3_image[15] + grid->fcrop4_image[15]))/100.0;
			grid->f_pasture_con = grid->fpast_unh_hmnzed[305] + 
				((grid->fgrass3_image[grid->lucy - 1990] + grid->fgrass4_image[grid->lucy - 1990])
				 - (grid->fgrass3_image[15] + grid->fgrass4_image[15]))/100.0;
		}
	}else if(LANDUSE == 9){
        /* 9: fixed land-use at 2000 --GEOMIP */
        grid->f_crop_con = grid->fcrop_unh_hmnzed[2000 - PIVOT_LUC];
        grid->f_pasture_con = grid->fpast_unh_hmnzed[2000 - PIVOT_LUC];
    }else if(LANDUSE==10 || LANDUSE==11 || LANDUSE==12 || LANDUSE==13
        || LANDUSE==14 || LANDUSE==15 || LANDUSE==16 || LANDUSE==17){
		/* UNH harmonized land-use change, 1500-2100 (added 2013/12/20) */
			grid->f_crop_con = grid->fcrop_unh_hmnzed[grid->lucy - PIVOT_LUC];
			grid->f_pasture_con = grid->fpast_unh_hmnzed[grid->lucy - PIVOT_LUC];
    }else{
		printf("Wrong land-use setting ID\n");
		exit(1);
	}
    
    /************************/
    if(NMIP_RUN >= 1){
        if(NMIP_RUN == 3){
            grid->f_crop_con = 0.0;
        }else{
            if(grid->lucy>=1900 && grid->lucy<=2012){
                grid->f_crop_con = grid->nmip_frcrop[grid->lucy - PIVOT_NINY];
            }else if(grid->lucy<1900){
                grid->f_crop_con = grid->nmip_frcrop[1900 - PIVOT_NINY];
            }else if(grid->lucy>2012){
                grid->f_crop_con = grid->nmip_frcrop[2012 - PIVOT_NINY];
            }
        }
    }
    
    /* decouple land-use emission (fix to 2000 land cover): 2014/11/20 by A.Ito */
    if(EX_CCPL == 3 || EX_CCPL == 8){
        grid->f_crop_con = grid->fcrop_unh_hmnzed[2000 - PIVOT_LUC];
        grid->f_pasture_con = grid->fpast_unh_hmnzed[2000 - PIVOT_LUC];
    }
	
    /* biofuel experiment: 2015/08/21 by A.Ito */
    if(LANDUSE==17 || BIOFUEL_RUN >= 1){
        if(grid->lucy <= 2010){
            grid->f_crop_con = grid->f_biofuel[0];
        }else{
            grid->f_crop_con = grid->f_biofuel[grid->lucy - 2010];
        }
        
        /* grid->f_pasture_con = grid->fpast_unh_hmnzed[2000 - PIVOT_LUC]; */
        grid->f_pasture_con = 0.0;
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
	if(grid->phase == 0){
		/* spin-up */
		if(LANDUSE == 0){
			grid->f_deforest = 0.0;
		}else if(LANDUSE>=1 && LANDUSE<=5){
            grid->f_deforest = grid->fcrop_sage[(BGY_LUC+1) - PIVOT_LUC]
                        - grid->fcrop_sage[BGY_LUC - PIVOT_LUC];
        }else if(LANDUSE==6 || LANDUSE==8 || LANDUSE==10 || LANDUSE==11 || LANDUSE==12 ||
                LANDUSE==13 || LANDUSE==14 || LANDUSE==15 || LANDUSE==16 || LANDUSE==17){
            grid->f_deforest = grid->t_vc_unh_hmnzed[BGY_LUC - PIVOT_LUC] 
                                + grid->t_vp_unh_hmnzed[BGY_LUC - PIVOT_LUC]
                                + grid->t_sc_unh_hmnzed[BGY_LUC - PIVOT_LUC] 
                                + grid->t_sp_unh_hmnzed[BGY_LUC - PIVOT_LUC];
            grid->f_deforest_v = grid->t_vc_unh_hmnzed[BGY_LUC - PIVOT_LUC] 
                                + grid->t_vp_unh_hmnzed[BGY_LUC - PIVOT_LUC];
            grid->f_deforest_s = grid->t_sc_unh_hmnzed[BGY_LUC - PIVOT_LUC] 
                                + grid->t_sp_unh_hmnzed[BGY_LUC - PIVOT_LUC];
        }else if(LANDUSE==7){
            grid->f_deforest = grid->fcrop_sage[(BGY_LUC+1) - PIVOT_LUC]
                        - grid->fcrop_sage[BGY_LUC - PIVOT_LUC];
        }else if(LANDUSE==9){
             grid->f_deforest = 0.0;
        }
		/* 2008/08/20 corrected by A.Ito (thanks to E.Kato) */
	
	}else if(grid->phase==1 || grid->phase==2){
		if(LANDUSE == 0){
			grid->f_deforest = 0.0;
		}else if(LANDUSE>=1 && LANDUSE<=5){
			grid->f_deforest = grid->f_crop_con - grid->f_crop_p;
		}else if(LANDUSE==6 || LANDUSE==8){
			if(grid->lucy <= 1999){
				grid->f_deforest = grid->t_vc_unh_hmnzed[grid->lucy - PIVOT_LUC]
									+ grid->t_vp_unh_hmnzed[grid->lucy - PIVOT_LUC]
									+ grid->t_sc_unh_hmnzed[grid->lucy - PIVOT_LUC]
									+ grid->t_sp_unh_hmnzed[grid->lucy - PIVOT_LUC];
				grid->f_deforest_v = grid->t_vc_unh_hmnzed[grid->lucy - PIVOT_LUC]
									+ grid->t_vp_unh_hmnzed[grid->lucy - PIVOT_LUC];
				grid->f_deforest_s = grid->t_sc_unh_hmnzed[grid->lucy - PIVOT_LUC]
									+ grid->t_sp_unh_hmnzed[grid->lucy - PIVOT_LUC];
			}else{
				grid->f_deforest = (grid->f_crop_con - grid->f_crop_p) 
									+ (grid->f_pasture_con - grid->f_pasture_p);
				grid->f_deforest_v = grid->f_deforest;
				grid->f_deforest_s = 0.0;
			}
		}else if(LANDUSE == 7){
			grid->f_deforest = (grid->f_crop_con - grid->f_crop_p) 
								+ (grid->f_pasture_con - grid->f_pasture_p);
		}else if(LANDUSE == 9){
            /* grid->f_deforest = grid->f_crop_con - grid->f_crop_p; */
            grid->f_deforest = 0.0;
        }else if(LANDUSE==10 || LANDUSE==11 || LANDUSE==12 || LANDUSE==13
             || LANDUSE==14 || LANDUSE==15 || LANDUSE==16 || LANDUSE==17){
            grid->f_deforest = grid->t_vc_unh_hmnzed[grid->lucy - PIVOT_LUC]
                                + grid->t_vp_unh_hmnzed[grid->lucy - PIVOT_LUC]
                                + grid->t_sc_unh_hmnzed[grid->lucy - PIVOT_LUC]
                                + grid->t_sp_unh_hmnzed[grid->lucy - PIVOT_LUC];
            grid->f_deforest_v = grid->t_vc_unh_hmnzed[grid->lucy - PIVOT_LUC]
                                + grid->t_vp_unh_hmnzed[grid->lucy - PIVOT_LUC];
            grid->f_deforest_s = grid->t_sc_unh_hmnzed[grid->lucy - PIVOT_LUC]
                                + grid->t_sp_unh_hmnzed[grid->lucy - PIVOT_LUC];
        }
	}
    
    /* biofuel experiment: 2015/08/27 by A.Ito */
    if(LANDUSE==17 || BIOFUEL_RUN >= 1){
        grid->f_deforest = 0.0;
    }
    
    if(NMIP_RUN >= 1){
        if(NMIP_RUN == 2 || NMIP_RUN == 3 || NMIP_RUN == 4 || NMIP_RUN == 5 || NMIP_RUN == 6
            || NMIP_RUN == 7 || NMIP_RUN == 8 || NMIP_RUN == 9 || NMIP_RUN == 10){
            grid->f_deforest = 0.0;
        }
    }

    /* parameter ensemble: 2014/11/19 by A.Ito */
    prm_ensen = 1.0;
    if(PARAM_PTB == 7){
        if(PARAM_PTB == 1){
            prm_ensen *= 0.7;
        }
        if(PARAM_PTB == 2){
            prm_ensen *= 0.8;
        }
        if(PARAM_PTB == 3){
            prm_ensen *= 0.9;
        }
        if(PARAM_PTB == 4){
            prm_ensen *= 1.1;
        }
        if(PARAM_PTB == 5){
            prm_ensen *= 1.2;
        }
        if(PARAM_PTB == 6){
            prm_ensen *= 1.3;
        }
        grid->f_deforest *= prm_ensen;
        grid->f_deforest_v *= prm_ensen;
        grid->f_deforest_s *= prm_ensen;
    }
	
	/* abandonment */
	if(grid->f_deforest < 0.0){
		grid->f_deforest = 0.0;
	}
	
	/* historical change in paddy field area: added by A.Ito (2011/2/28) ********/
	if(LANDUSE == 0 || NMIP_RUN == 3){
		grid->f_paddy = 0.0;
	}else if(LANDUSE == 7){
		if(grid->lucy <= 2005){
			if(grid->f_paddy_b > 0.0 && grid->fcrop_rk[2000 - PIVOT_LUC] > 0.0){
				grid->f_paddy = grid->f_paddy_b * 
					(grid->fcrop_rk[grid->lucy - PIVOT_LUC] / grid->fcrop_rk[2000 - PIVOT_LUC]);
			}
		}else if(grid->lucy >= 2006){
			if(grid->f_paddy_b > 0.0 && grid->fcrop_rk[2000 - PIVOT_LUC] > 0.0){
				grid->f_paddy = grid->f_paddy_b * 
					(grid->fcrop_rk[2005 - PIVOT_LUC] / grid->fcrop_rk[2000 - PIVOT_LUC]);
			}
		}
		
		if(grid->f_paddy > 1.0){
			grid->f_paddy = 1.0;
		}
		if(grid->f_paddy < 0.0){
			grid->f_paddy = 0.0;
		}
	}else if(LANDUSE == 8){
		if(grid->lucy <= 2005){
        
			if(grid->f_paddy_b > 0.0 && grid->fcrop_unh_hmnzed[2000 - PIVOT_LUC] > 0.0){
				grid->f_paddy = grid->f_paddy_b * 
                    (grid->fcrop_unh_hmnzed[grid->lucy - PIVOT_LUC] / grid->fcrop_unh_hmnzed[2000 - PIVOT_LUC]);
			}
		}else if(grid->lucy >= 2006){
			if(grid->f_paddy_b > 0.0 && grid->fcrop_unh_hmnzed[2000 - PIVOT_LUC] > 0.0){
				grid->f_paddy = grid->f_paddy_b * 
                    (grid->fcrop_unh_hmnzed[2005 - PIVOT_LUC] / grid->fcrop_unh_hmnzed[2000 - PIVOT_LUC]);
			}
		}
		
		if(grid->f_paddy > 1.0){
			grid->f_paddy = 1.0;
		}
		if(grid->f_paddy < 0.0){
			grid->f_paddy = 0.0;
		}
	}else if(LANDUSE == 9){
        grid->f_paddy = grid->f_paddy_b;
    }else if(LANDUSE==10 || LANDUSE==11 || LANDUSE==12 || LANDUSE==13
        || LANDUSE==14 || LANDUSE==15 || LANDUSE==16 || LANDUSE==17){
        
        if(grid->f_paddy_b > 0.0 && grid->fcrop_unh_hmnzed[2000 - PIVOT_LUC] > 0.0){
            grid->f_paddy = grid->f_paddy_b * 
                (grid->fcrop_unh_hmnzed[grid->lucy - PIVOT_LUC] / grid->fcrop_unh_hmnzed[2000 - PIVOT_LUC]);
        }
        
        if(grid->f_paddy > 1.0){
			grid->f_paddy = 1.0;
		}
		if(grid->f_paddy < 0.0){
			grid->f_paddy = 0.0;
		}
    }
    
    /* biofuel experiment: 2015/08/27 by A.Ito */
    if(LANDUSE==17 || BIOFUEL_RUN >= 1){
        grid->f_paddy = grid->f_paddy_b;
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
	
    /* added: A. Ito (with Hamada-san's comment) 2012/01/30 */
    f_luc = 0.0;
    
	if(grid->phase==0){
		/* spin-up: fluxes for 1801-1900 *******************************/
		/* modified by A.Ito (2009/06/05: 2010/01/07) */
		if(LANDUSE == 0 || LANDUSE==9){
			f_luc = 0.0;
		}else if(LANDUSE>=1 && LANDUSE<=5){
			f_luc = grid->fcrop_sage[BGY_LUC-PIVOT_LUC] - grid->fcrop_sage[BGY_LUC-PIVOT_LUC-1];
		}else if(LANDUSE==6 || LANDUSE==8 || LANDUSE==10 || LANDUSE==11 || LANDUSE==12 ||
                LANDUSE==13 || LANDUSE==14 || LANDUSE==15 || LANDUSE==16 || LANDUSE==17){
            /* revised (after comments by E.Kato): 2013/10/02 by A.Ito */
			f_luc = (grid->t_vc_unh_hmnzed[BGY_LUC - PIVOT_LUC] + grid->t_vp_unh_hmnzed[BGY_LUC - PIVOT_LUC])
				+ (grid->t_sc_unh_hmnzed[BGY_LUC - PIVOT_LUC] + grid->t_sp_unh_hmnzed[BGY_LUC - PIVOT_LUC])*f_mass_secfor;
		}else if(LANDUSE==7){
			/* added 2010/01/07 (A.Ito) */
			f_luc = (grid->fcrop_rk[BGY_LUC-PIVOT_LUC] - grid->fcrop_rk[BGY_LUC-PIVOT_LUC-1])
					+(grid->fpast_rk[BGY_LUC-PIVOT_LUC] - grid->fpast_rk[BGY_LUC-PIVOT_LUC-1]);
		}
		
        /* NMIP: fixed land-use */
        if(NMIP_RUN == 2 || NMIP_RUN == 3 || NMIP_RUN == 4 || NMIP_RUN == 5 || NMIP_RUN == 6
            || NMIP_RUN == 7 || NMIP_RUN == 8 || NMIP_RUN == 9 || NMIP_RUN == 10){
            f_luc = 0.0;
        }

		/* modified by A.Ito (2009/08/19) */
		if(f_luc > 0.0){
			flux->lu_detr = f_luc * 0.2*(mass->plant).rot;
			flux->lu_conv = f_luc * eff_mass * fe_conv/(fe_conv + fe_ten + fe_hund);
		}else{
			flux->lu_detr = 0.0;
			flux->lu_conv = 0.0;
		}
		
		for(f=(BGY_LUC-9);f<=BGY_LUC;f++){
			/* senstivity analysis */
			if(LANDUSE == 0 || LANDUSE==9){
				f_luc = 0.0;
			}else if(LANDUSE>=1 && LANDUSE<=5){
				f_luc = grid->fcrop_sage[f - PIVOT_LUC] - grid->fcrop_sage[f - PIVOT_LUC - 1];
			}else if(LANDUSE==6 || LANDUSE==8 || LANDUSE==10 || LANDUSE==11 || LANDUSE==12 ||
                    LANDUSE==13 || LANDUSE==14 || LANDUSE==15 || LANDUSE==16 || LANDUSE==17){
				f_luc = (grid->t_vc_unh_hmnzed[f - PIVOT_LUC] + grid->t_vp_unh_hmnzed[f - PIVOT_LUC])
						+ (grid->t_sc_unh_hmnzed[f - PIVOT_LUC] + grid->t_sp_unh_hmnzed[f - PIVOT_LUC]) * f_mass_secfor;
				/* 0.5: assumption by A.Ito for secondary forest stock */
			}else if(LANDUSE==7){
				/* added 2010/01/07 (A.Ito) */
				f_luc = (grid->fcrop_rk[f - PIVOT_LUC] - grid->fcrop_rk[f-PIVOT_LUC-1])
						+ (grid->fpast_rk[f - PIVOT_LUC] - grid->fpast_rk[f-PIVOT_LUC-1]);
			}
			
            /* NMIP: fixed land-use */
            if(NMIP_RUN == 2 || NMIP_RUN == 3 || NMIP_RUN == 4 || NMIP_RUN == 5 || NMIP_RUN == 6
                || NMIP_RUN == 7 || NMIP_RUN == 8 || NMIP_RUN == 9 || NMIP_RUN == 10){
                f_luc = 0.0;
            }

			/* modified by A.Ito based on E.Kato (2009/03/30) */
			if(f_luc > 0.0){
				mass_ten = f_luc * eff_mass * fe_ten/(fe_conv + fe_ten + fe_hund);
				flux->detr_ten[BGY_LUC-f] = mass_ten;
			}else{
                /* corrected: A. Ito (with Hamada-san's comment) 2012/01/30 */
				mass_ten = 0.0;
				flux->detr_ten[BGY_LUC-f] = 0.0;
			}
			/* corrected: A.Ito and E.Kato (2009/08/16) */
			flux->lu_ten += mass_ten*0.1;
		}
		
		for(f=(BGY_LUC-99);f<=BGY_LUC;f++){
			/* senstivity analysis */
			if(LANDUSE == 0 || LANDUSE==9){
				f_luc = 0.0;
			}else if(LANDUSE>=1 && LANDUSE<=5){
				f_luc = grid->fcrop_sage[f - PIVOT_LUC] - grid->fcrop_sage[f-PIVOT_LUC-1];
			}else if(LANDUSE==6 || LANDUSE==8 || LANDUSE==10 || LANDUSE==11 || LANDUSE==12 ||
                    LANDUSE==13 || LANDUSE==14 || LANDUSE==15 || LANDUSE==16 || LANDUSE==17){
				f_luc = (grid->t_vc_unh_hmnzed[f - PIVOT_LUC] + grid->t_vp_unh_hmnzed[f - PIVOT_LUC])
						+ (grid->t_sc_unh_hmnzed[f - PIVOT_LUC] + grid->t_sp_unh_hmnzed[f - PIVOT_LUC])*f_mass_secfor;
			}else if(LANDUSE==7){
				/* added 2010/01/07 (A.Ito) */
				f_luc = (grid->fcrop_rk[f - PIVOT_LUC] - grid->fcrop_rk[f-PIVOT_LUC-1])
						+ (grid->fpast_rk[f - PIVOT_LUC] - grid->fpast_rk[f-PIVOT_LUC-1]);
			}
			
            /* NMIP: fixed land-use */
            if(NMIP_RUN == 2 || NMIP_RUN == 3 || NMIP_RUN == 4 || NMIP_RUN == 5 || NMIP_RUN == 6
                || NMIP_RUN == 7 || NMIP_RUN == 8 || NMIP_RUN == 9 || NMIP_RUN == 10){
                f_luc = 0.0;
            }

			/* modified by A.Ito based on E.Kato (2009/03/30) */
			/* corrected: A.Ito and E.Kato (2009/08/16) */
			if(f_luc > 0.0){
				mass_hund = f_luc * eff_mass * fe_hund/(fe_conv + fe_ten + fe_hund);				
				flux->detr_hund[BGY_LUC - f] = mass_hund;
			}else{
				mass_hund = 0.0;
				flux->detr_hund[BGY_LUC - f] = 0.0;
			}
			flux->lu_hund += mass_hund * 0.01;
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
		if(LANDUSE == 0 || LANDUSE == 9){
			f_luc = 0.0;
		}else if(LANDUSE>=1 && LANDUSE<=5){
			f_luc = grid->f_deforest;
			/*  grid->f_crop_con - grid->f_crop_p;  */
		}else if(LANDUSE==6 || LANDUSE==8 || LANDUSE==10 || LANDUSE==11 || LANDUSE==12 ||
                LANDUSE==13 || LANDUSE==14 || LANDUSE==15 || LANDUSE==16 || LANDUSE==17){
			/* assumption: biomass in secondary forest is lower (0.1) than primary forest */
			f_luc = grid->f_deforest_v + grid->f_deforest_s * f_mass_secfor;
		}else if(LANDUSE == 7){
			/* added 2010/01/07 (A.Ito) */
			f_luc = grid->f_deforest;
		}
		
        /* NMIP: fixed land-use */
        if(NMIP_RUN == 2 || NMIP_RUN == 3 || NMIP_RUN == 4 || NMIP_RUN == 5 || NMIP_RUN == 6
            || NMIP_RUN == 7 || NMIP_RUN == 8 || NMIP_RUN == 9 || NMIP_RUN == 10){
            f_luc = 0.0;
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
