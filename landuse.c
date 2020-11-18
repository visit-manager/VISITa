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
	}else if(LANDUSE == 1 || LANDUSE == 2 || LANDUSE == 3 || LANDUSE == 4 || LANDUSE == 5){
		/* SAGE land-use data:
		 Ramankutty, N., and J. A. Foley (1999), Estimating historical changes in global 
		 land cover: croplands from 1700 to 1992, Global Biogeochemical Cycles, 13(4), 997-1027.
        */
		if(grid->lucy <= 1990){
			/* SAGE, net land use, only cropland */
			grid->f_crop_con = grid->fcrop_net[grid->lucy - FDY_LUC];
			grid->f_pasture_con = 0.0;
		}else{
			switch(LANDUSE){
				case 1:
					/* fixed 1900 */
					grid->f_crop_con = grid->fcrop_net[200]; /*  */
					break;
				case 2:
					/* fixed 1990 */
					grid->f_crop_con = grid->fcrop_net[290]; /*  */
					break;
				case 3:
					/* conventional scenario : ex b, d */
					grid->f_crop_con = grid->fcrop_net[290] + grid->f_crop_trend 
										*(double)(grid->lucy - 1990);
					break;
				case 4:
					/* high scenario : ex e */
					if(grid->f_crop_trend > 0.0){
						grid->f_crop_con = grid->fcrop_net[290] + grid->f_crop_trend*(1.0+0.01*
								((double)(grid->lucy - 1990))) *(double)(grid->lucy - 1990);
					}else{
						grid->f_crop_con = grid->fcrop_net[290];
					}
					break;
				case 5:
					/* IMAGE-based scenario */
					grid->f_crop_con = grid->fcrop_net[290] + ((grid->fcrop3_image[grid->lucy - 1990]
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
			grid->f_crop_con = grid->fcrop_luh[grid->lucy - FDY_LUC];
			grid->f_pasture_con = grid->fpast_luh[grid->lucy - FDY_LUC];
		}else if(grid->lucy >= 2000){
			grid->f_crop_con = grid->fcrop_luh[299] + 
				((grid->fcrop3_image[grid->lucy - 1990] + grid->fcrop4_image[grid->lucy - 1990])
				- (grid->fcrop3_image[9]+grid->fcrop4_image[9]))/100.0;
			grid->f_pasture_con = grid->fpast_luh[299] + 
				((grid->fgrass3_image[grid->lucy - 1990] + grid->fgrass4_image[grid->lucy - 1990])
				- (grid->fgrass3_image[9] + grid->fgrass4_image[9]))/100.0;
		}
	}else if(LANDUSE == 7){
		/* SAGE land-use data 1700-2007 (Revised: 2010/01/07):
		 Ramankutty, N., and J. A. Foley (1999), Estimating historical changes in global 
		 land cover: croplands from 1700 to 1992, Global Biogeochemical Cycles, 13(4), 997-1027.
		*/
		if(grid->lucy <= 2007){
			grid->f_crop_con = grid->fcrop_rk[grid->lucy - FDY_LUC];
			grid->f_pasture_con = grid->fpast_rk[grid->lucy - FDY_LUC];
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
			grid->f_crop_con = grid->fcrop_luh[grid->lucy - FDY_LUC];
			grid->f_pasture_con = grid->fpast_luh[grid->lucy - FDY_LUC];
		}else if(grid->lucy >= 2006){
			grid->f_crop_con = grid->fcrop_luh[305] + 
				((grid->fcrop3_image[grid->lucy - 1990] + grid->fcrop4_image[grid->lucy - 1990])
				 - (grid->fcrop3_image[15] + grid->fcrop4_image[15]))/100.0;
			grid->f_pasture_con = grid->fpast_luh[305] + 
				((grid->fgrass3_image[grid->lucy - 1990] + grid->fgrass4_image[grid->lucy - 1990])
				 - (grid->fgrass3_image[15] + grid->fgrass4_image[15]))/100.0;
		}
	}else if(LANDUSE == 9){
        /* 9: fixed land-use at 2000 --GEOMIP */
        grid->f_crop_con = grid->fcrop_luh[2000 - FDY_LUC];
        grid->f_pasture_con = grid->fpast_luh[2000 - FDY_LUC];
    }else if(LANDUSE == 10 || LANDUSE == 11 || LANDUSE == 12 || LANDUSE == 13
        || LANDUSE == 14 || LANDUSE == 15 || LANDUSE == 16 || LANDUSE == 17
        || LANDUSE == 26 || LANDUSE == 27 || LANDUSE == 28
        || LANDUSE == 30 || LANDUSE == 31 || LANDUSE == 32 || LANDUSE == 33
        || LANDUSE == 34 || LANDUSE == 35 || LANDUSE == 36 || LANDUSE == 37){
		/* UNH harmonized land-use change, 1500-2100 (added 2013/12/20) */
            if((grid->lucy - FDY_LUC) < DL_LUC){
                grid->f_crop_con = grid->fcrop_luh[grid->lucy - FDY_LUC];
                grid->f_pasture_con = grid->fpast_luh[grid->lucy - FDY_LUC];
            }else{
                grid->f_crop_con = grid->fcrop_luh[DL_LUC -1];
                grid->f_pasture_con = grid->fpast_luh[DL_LUC -1];
            }
    }else if(LANDUSE == 18){
        /* ICARUS SSPs: 2016/08/14 by A.Ito */
        if(grid->lucy <= 1990){
            grid->f_crop_con = grid->fcrop3_image[0];
        }else{
            grid->f_crop_con = grid->fcrop3_image[grid->lucy - 1990];
        }
        grid->f_pasture_con = 0.0;
    }else if(LANDUSE == 19 || LANDUSE == 20 || LANDUSE == 21 ||
                LANDUSE == 22 || LANDUSE == 23){
        /* CD-LINKS SSPs: 2016/11/20 by A.Ito */
        if(grid->lucy < 2000){
            grid->f_crop_con = grid->fcrop_luh[1999 - FDY_LUC];
        }else{
            grid->f_crop_con = grid->fcrop3_image[grid->lucy - 1990];
        }
        grid->f_pasture_con = 0.0;
    }else if(LANDUSE == 24 || LANDUSE == 25){
        /* ISI-MIP2b: 2016/12/24 by A.Ito */
        grid->f_crop_con = grid->mip_frcrop[grid->lucy - FDY_NINY];
        
    }else if(LANDUSE == 29){
        /* MIROC0INTEG TELUMO: 2019/02/21 by A.Ito */
        if(grid->lucy < 1950){
            grid->f_crop_con = grid->fcrop_luh[1950 - FDY_LUC];
        }else if(grid->lucy > 2099){
            grid->f_crop_con = grid->fcrop_luh[2099 - FDY_LUC];
        }else{
            grid->f_crop_con = grid->fcrop_luh[grid->lucy - FDY_LUC];
        }
        if(grid->lucy < 1950){
            grid->f_pasture_con = grid->fcrop_luh[1950 - FDY_LUC];
        }else if(grid->lucy > 2099){
            grid->f_pasture_con = grid->fcrop_luh[2099 - FDY_LUC];
        }else{
            grid->f_pasture_con = grid->fpast_luh[grid->lucy - FDY_LUC];
        }
    }else if(LANDUSE == 38 || LANDUSE == 39 || LANDUSE == 40 || LANDUSE == 41){
        /* AIM data, 2019/06/21 by A.Ito */
        if(grid->lucy < 1950){
            grid->f_crop_con = grid->fcrop_luh[1949 - FDY_LUC];
        }else if(grid->lucy >= 1950 && grid->lucy <= 1999){
            grid->f_crop_con = grid->fcrop_luh[grid->lucy - FDY_LUC];
        }else if(grid->lucy >= 2000 && grid->lucy <= 2100){
            grid->f_crop_con = grid->aim_luc_fcrop[grid->lucy-2000];
        }else{
            grid->f_crop_con = grid->aim_luc_fcrop[2100-2000];
        }
        grid->f_pasture_con = 0.0;
    }else if(LANDUSE == 45){
        /* ISIMIP3a: 2020/10/01 by A.Ito */
        if(grid->climy < 1850){
            grid->lucy = 1850;
        }else if(grid->climy >= 1850 && grid->climy <= 2018){
            grid->lucy = grid->climy;
        }else{
            grid->lucy = 2018;
        }
        grid->f_crop_con = grid->mip_frcrop[grid->lucy - FDY_NINY];
        
    }else if(LANDUSE == 46){
        /* ISIMIP3b: 2020/11/18 by A.Ito */
        if(grid->climy < 1601){
            grid->lucy = 1601;
        }else if(grid->climy >= 1601 && grid->climy <= 2100){
            grid->lucy = grid->climy;
        }else{
            grid->lucy = 2100;
        }
        grid->f_crop_con = grid->mip_frcrop[grid->lucy - FDY_NINY];
    
    }else if(LANDUSE == 47){
        /* 2-2002 S1: 2020/10/08 by A.Ito */
        /* FDY_LUC = 1866 */
        if(grid->lucy < 1950){
            grid->f_crop_con = grid->fcrop_luh[1949 - FDY_LUC];
            grid->f_pasture_con = grid->fpast_luh[1949 - FDY_LUC];
        }else if(grid->lucy >= 1950 && grid->lucy <= 2004){
            grid->f_crop_con = grid->fcrop_luh[grid->lucy - FDY_LUC];
            grid->f_pasture_con = grid->fpast_luh[grid->lucy - FDY_LUC];
        }else if(grid->lucy >= 2005 && grid->lucy <= 2100){
            grid->f_crop_con = grid->aim_luc_fcrop[grid->lucy - 2005];
            grid->f_pasture_con = grid->aim_luc_grass[grid->lucy - 2005];
        }else if(grid->lucy > 2100){
            grid->f_crop_con = grid->aim_luc_fcrop[2100 - 2005];
            grid->f_pasture_con = grid->aim_luc_grass[2100 - 2005];
        }
        
    }else{
		printf("Wrong land-use setting ID\n");
		exit(1);
	}
 
    grid->f_crop_ans = grid->f_crop_con;
    if(EXTRA_LU_FIX == 2){
        grid->f_crop_con = 0.0;
    }
    
    /*************************************************************/
    if(NMIP_RUN >= 1){
        if(grid->lucy >= FDY_NINY && grid->lucy <= 2015){
            grid->f_crop_con = grid->mip_frcrop[grid->lucy - FDY_NINY];
        }else if(grid->lucy < FDY_NINY){
            grid->f_crop_con = grid->mip_frcrop[FDY_NINY - FDY_NINY];
        }else if(grid->lucy > 2015){
            grid->f_crop_con = grid->mip_frcrop[2015 - FDY_NINY];
        }
    }
    
    /* S10-BECCS: 2016/02/15 by A.Ito ***********************/
    if(EX_BECCS == 1){
        if(grid->lucy <= 2000){
            grid->f_crop_con = grid->fcrop_luh[2000 - FDY_LUC];
            grid->f_pasture_con = grid->fpast_luh[2000 - FDY_LUC];
        }else if(grid->lucy >= 2001){
            if(grid->veg_sage>=1 && grid->veg_sage<=8){
                
                grid->f_crop_con = grid->fcrop_luh[2000 - FDY_LUC]
                        + (double)(grid->lucy - 2000) * 0.1 * (1.0 - grid->fcrop_luh[2000 - FDY_LUC])/100.0;
                
                grid->f_pasture_con = grid->fpast_luh[2000 - FDY_LUC];
            }else{
                grid->f_crop_con = grid->fcrop_luh[2000 - FDY_LUC];
                grid->f_pasture_con = grid->fpast_luh[2000 - FDY_LUC];
            }
        }
    }else if(EX_BECCS == 2 || EX_BECCS == 3){
        /* 2017/02/25, 2017/10/17 by A.Ito */
        if(grid->lucy <= 1999){
            grid->f_crop_con = grid->fcrop_luh[2000 - FDY_LUC];
            grid->f_pasture_con = grid->fpast_luh[2000 - FDY_LUC];
        }else if(grid->lucy >= 2000){
            /* if(grid->veg_sage>=1 && grid->veg_sage<=8){ */
            if(grid->veg_sage>=1 && grid->veg_sage<=15){

                grid->f_crop_con = grid->fcrop_luh[grid->lucy - FDY_LUC]
                        + (double)(grid->lucy - 1999) * (grid->beccs_s2b
                             + grid->beccs_v2b + grid->beccs_v2s)/100.0;
                
                if(EX_BECCS_SUB == 2){
                    grid->f_crop_con = grid->fcrop_luh[2000 - FDY_LUC];
                            /* + (double)(grid->lucy - 1999) * (grid->beccs_s2b
                                 + grid->beccs_v2b + grid->beccs_v2s)/100.0; */
                }
                
                if(grid->f_crop_con > 1.0){
                    grid->f_crop_con = 1.0;
                }
                
                grid->f_pasture_con = grid->fpast_luh[grid->lucy - FDY_LUC];
            }else{
                grid->f_crop_con = grid->fcrop_luh[grid->lucy - FDY_LUC];
                grid->f_pasture_con = grid->fpast_luh[grid->lucy - FDY_LUC];
            }
        }
    }
    
    /* decouple land-use emission (fix to 2000 land cover): 2014/11/20 by A.Ito */
    /* if(EX_CCPL == 3 || EX_CCPL == 8){
        grid->f_crop_con = grid->fcrop_luh[2000 - FDY_LUC];
        grid->f_pasture_con = grid->fpast_luh[2000 - FDY_LUC];
    } */
	
    /* biofuel experiment: 2015/08/21 by A.Ito */
    if(LANDUSE == 17 || BIOFUEL_RUN >= 1){
        if(grid->lucy <= 2010){
            grid->f_crop_con = grid->f_biofuel[0];
        }else{
            grid->f_crop_con = grid->f_biofuel[grid->lucy - 2010];
        }
        
        /* grid->f_pasture_con = grid->fpast_luh[2000 - FDY_LUC]; */
        grid->f_pasture_con = 0.0;
    }
    
    /* Forest management: 2019/10/15 by A.Ito *********/
    if(EX_FORMAN == 1){
        if(grid->lucy<=2000){
            grid->f_crop_con = grid->fcrop_luh[grid->lucy - FDY_LUC];
            grid->f_pasture_con = grid->fpast_luh[grid->lucy - FDY_LUC];
        }else{
            grid->f_crop_con = grid->fcrop_luh[2000 - FDY_LUC] * (1.0 - (double)(grid->lucy-2000)*0.005 );
            grid->f_pasture_con = grid->fpast_luh[2000 - FDY_LUC];
        }
    }
    
    /* forced afforestation: 2020/11/09 by A.Ito */
    if((EX_FORCED_AFFOREST_2 >=1 && EX_FORCED_AFFOREST_2 <=12) && grid->lucy >= EX_FORCED_AFFOREST_2_YR){
        grid->f_crop_con = 0.0;
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

	/*****************************************************************/
	/* annual deforestation */
    grid->f_deforest = grid->f_deforest_v = grid->f_deforest_s = 0.0;
	if(grid->phase == 0){
		/* spin-up */
		if(LANDUSE == 0){
			grid->f_deforest = 0.0;
		}else if(LANDUSE>=1 && LANDUSE<=5){
            grid->f_deforest = grid->fcrop_net[(BGY_LUC+1) - FDY_LUC]
                        - grid->fcrop_net[BGY_LUC - FDY_LUC];
        }else if(LANDUSE==6 || LANDUSE==8 || LANDUSE==10 || LANDUSE==11 || LANDUSE==12
                 || LANDUSE==13 || LANDUSE==14 || LANDUSE==15 || LANDUSE==16 || LANDUSE==17
                 || LANDUSE==26 || LANDUSE==27 || LANDUSE==28 || LANDUSE == 30
                 || LANDUSE == 31 || LANDUSE == 32 || LANDUSE == 33
                 || LANDUSE == 34 || LANDUSE == 35 || LANDUSE == 36 || LANDUSE == 37){
            grid->f_deforest = grid->t_vc_luh[BGY_LUC - FDY_LUC]
                                + grid->t_vp_luh[BGY_LUC - FDY_LUC]
                                + grid->t_sc_luh[BGY_LUC - FDY_LUC]
                                + grid->t_sp_luh[BGY_LUC - FDY_LUC];
            grid->f_deforest_v = grid->t_vc_luh[BGY_LUC - FDY_LUC]
                                + grid->t_vp_luh[BGY_LUC - FDY_LUC];
            grid->f_deforest_s = grid->t_sc_luh[BGY_LUC - FDY_LUC]
                                + grid->t_sp_luh[BGY_LUC - FDY_LUC];
        }else if(LANDUSE == 7){
            grid->f_deforest = grid->fcrop_net[(BGY_LUC+1) - FDY_LUC]
                        - grid->fcrop_net[BGY_LUC - FDY_LUC];
        }else if(LANDUSE==9){
            grid->f_deforest = 0.0;
        }else if(LANDUSE==18 || LANDUSE == 19 || LANDUSE == 20 ||
                LANDUSE == 21 || LANDUSE == 22 || LANDUSE == 23 ||
                LANDUSE == 24 || LANDUSE == 25){
            /* ICARUS SSPs: 2016/08/14 by A.Ito */
            grid->f_deforest = 0.0;
        }else if(LANDUSE == 29){
            /* MIROC-INTEG TELUMO 2019/02/21 by A.Ito */
            grid->f_deforest = grid->t_vc_luh[BGY_LUC - FDY_LUC]
                                + grid->t_vp_luh[BGY_LUC - FDY_LUC]
                                + grid->t_sc_luh[BGY_LUC - FDY_LUC]
                                + grid->t_sp_luh[BGY_LUC - FDY_LUC];
            grid->f_deforest_v = grid->t_vc_luh[BGY_LUC - FDY_LUC]
                                + grid->t_vp_luh[BGY_LUC - FDY_LUC];
            grid->f_deforest_s = grid->t_sc_luh[BGY_LUC - FDY_LUC]
                                + grid->t_sp_luh[BGY_LUC - FDY_LUC];
        }else if(LANDUSE == 38 || LANDUSE == 39 || LANDUSE == 40 || LANDUSE == 41){
            /* AIM land use, 2019/06/21 by A.Ito */
            grid->f_deforest = grid->fcrop_luh[1950 - FDY_LUC] - grid->fcrop_luh[1949 - FDY_LUC];
            if(EXTRA_LU_FIX == 2){
                grid->f_deforest = 0.0;
            }
        }else if(LANDUSE == 45){
            /* ISIMIP3a: 2020/10/01 by A.Ito */
            /* spin-up: 1801–1900 */
            grid->f_deforest = grid->f_crop_con - grid->f_crop_p;
        }else if(LANDUSE == 46){
            /* ISIMIP3b: 2020/11/18 by A.Ito */
            /* spin-up: 1601 */
            grid->f_deforest = grid->f_crop_con - grid->f_crop_p;
        }
		/* 2008/08/20 corrected by A.Ito (thanks to E.Kato) */
        
        if(EX_BECCS == 1 || EX_BECCS == 2 || EX_BECCS == 3){
            grid->f_deforest = 0.0;
        }
	
	}else if(grid->phase == 1 || grid->phase == 2){
		if(LANDUSE == 0){
			grid->f_deforest = 0.0;
		}else if(LANDUSE >= 1 && LANDUSE <= 5){
			grid->f_deforest = grid->f_crop_con - grid->f_crop_p;
		}else if(LANDUSE == 6 || LANDUSE == 8){
			if(grid->lucy <= 1999){
				grid->f_deforest = grid->t_vc_luh[grid->lucy - FDY_LUC]
									+ grid->t_vp_luh[grid->lucy - FDY_LUC]
									+ grid->t_sc_luh[grid->lucy - FDY_LUC]
									+ grid->t_sp_luh[grid->lucy - FDY_LUC];
				grid->f_deforest_v = grid->t_vc_luh[grid->lucy - FDY_LUC]
									+ grid->t_vp_luh[grid->lucy - FDY_LUC];
				grid->f_deforest_s = grid->t_sc_luh[grid->lucy - FDY_LUC]
									+ grid->t_sp_luh[grid->lucy - FDY_LUC];
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
        }else if(LANDUSE == 10 || LANDUSE == 11 || LANDUSE == 12 || LANDUSE == 13
             || LANDUSE == 14 || LANDUSE == 15 || LANDUSE == 16 || LANDUSE == 17
             || LANDUSE == 26 || LANDUSE == 27 || LANDUSE == 28 || LANDUSE == 30
             || LANDUSE == 31 || LANDUSE == 32 || LANDUSE == 33
             || LANDUSE == 34 || LANDUSE == 35 || LANDUSE == 36 || LANDUSE == 37){
            if((grid->lucy - FDY_LUC) < DL_LUC){
                grid->f_deforest = grid->t_vc_luh[grid->lucy - FDY_LUC]
                                    + grid->t_vp_luh[grid->lucy - FDY_LUC]
                                    + grid->t_sc_luh[grid->lucy - FDY_LUC]
                                    + grid->t_sp_luh[grid->lucy - FDY_LUC];
                grid->f_deforest_v = grid->t_vc_luh[grid->lucy - FDY_LUC]
                                    + grid->t_vp_luh[grid->lucy - FDY_LUC];
                grid->f_deforest_s = grid->t_sc_luh[grid->lucy - FDY_LUC]
                                    + grid->t_sp_luh[grid->lucy - FDY_LUC];
            }else{
                grid->f_deforest = grid->t_vc_luh[DL_LUC - 1]
                                    + grid->t_vp_luh[DL_LUC - 1]
                                    + grid->t_sc_luh[DL_LUC - 1]
                                    + grid->t_sp_luh[DL_LUC - 1];
                grid->f_deforest_v = grid->t_vc_luh[DL_LUC - 1]
                                    + grid->t_vp_luh[DL_LUC - 1];
                grid->f_deforest_s = grid->t_sc_luh[DL_LUC - 1]
                                    + grid->t_sp_luh[DL_LUC - 1];
            }
        }else if(LANDUSE == 18){
            /* ICARUS SSPs: 2016/08/14 by A.Ito */
            grid->f_deforest = 0.0;
        }else if (LANDUSE == 19 || LANDUSE == 20
                  || LANDUSE == 21 || LANDUSE == 22 || LANDUSE == 23){
            /* CD-LINKS: 2016/11/21 by A.Ito */
            if(grid->lucy < 1991){
                /* grid->f_deforest = 0.0; */
                grid->f_deforest = grid->t_vc_luh[grid->lucy - FDY_LUC]
                                    + grid->t_vp_luh[grid->lucy - FDY_LUC];
            }else{
                grid->f_deforest = grid->fcrop3_image[grid->lucy - 1990]
                            - grid->fcrop3_image[grid->lucy - 1990 -1];
            }
        }else if(LANDUSE == 24 || LANDUSE == 25){
			grid->f_deforest = grid->f_crop_con - grid->f_crop_p;
		}else if(LANDUSE == 29){
            /* MIROC-INTEG TELUMO 2019/02/21 by A.Ito */
            if((grid->lucy - FDY_LUC) < DL_LUC){
                grid->f_deforest = grid->t_vc_luh[grid->lucy - FDY_LUC]
                                    + grid->t_vp_luh[grid->lucy - FDY_LUC]
                                    + grid->t_sc_luh[grid->lucy - FDY_LUC]
                                    + grid->t_sp_luh[grid->lucy - FDY_LUC];
                grid->f_deforest_v = grid->t_vc_luh[grid->lucy - FDY_LUC]
                                    + grid->t_vp_luh[grid->lucy - FDY_LUC];
                grid->f_deforest_s = grid->t_sc_luh[grid->lucy - FDY_LUC]
                                    + grid->t_sp_luh[grid->lucy - FDY_LUC];
            }else{
                grid->f_deforest = grid->t_vc_luh[DL_LUC - 1]
                                    + grid->t_vp_luh[DL_LUC - 1]
                                    + grid->t_sc_luh[DL_LUC - 1]
                                    + grid->t_sp_luh[DL_LUC - 1];
                grid->f_deforest_v = grid->t_vc_luh[DL_LUC - 1]
                                    + grid->t_vp_luh[DL_LUC - 1];
                grid->f_deforest_s = grid->t_sc_luh[DL_LUC - 1]
                                    + grid->t_sp_luh[DL_LUC - 1];
            }
        }else if(LANDUSE == 38 || LANDUSE == 39 || LANDUSE == 40 || LANDUSE == 41){
            /* AIM land use, 2019/06/21 by A.Ito */
            grid->f_deforest = grid->f_crop_con - grid->f_crop_p;
        }else if(LANDUSE == 45){
            grid->f_deforest = grid->f_crop_con - grid->f_crop_p;
        }else if(LANDUSE == 46){
            /* ISIMIP3b: 2020/11/18 by A.Ito */
            grid->f_deforest = grid->f_crop_con - grid->f_crop_p;
        }
        
        /* BECCS S10 experiment: 2016/02/16 by A.Ito ******/
        if(EX_BECCS == 1){
            if(grid->lucy <= 2000){
                grid->f_deforest = grid->f_deforest_v = grid->f_deforest_s = 0.0;
            }else if(grid->lucy >= 2001){
                if(grid->veg_sage>=1 && grid->veg_sage<=8){
                    grid->f_deforest = 0.1*(1.0 - grid->fcrop_luh[2000 - FDY_LUC])/100.0;
                }else{
                    grid->f_deforest = grid->f_deforest_v = grid->f_deforest_s = 0.0;
                }
            }
        }else if(EX_BECCS == 2 || EX_BECCS == 3){
            /* 2017/02/25 by A.Ito */
            
            if(grid->lucy <= 2000){
                grid->f_deforest = grid->f_deforest_v = grid->f_deforest_s = 0.0;
            }else if(grid->lucy >= 2001){
                if(grid->veg_sage>=1 && grid->veg_sage<=8){
                
                    /* aa = 1.0 - grid->fcrop_luh[2000 - FDY_LUC];
                    bb = grid->beccs_v2s + grid->beccs_v2b + grid->beccs_s2b;
                    max_conv = (aa>bb)?bb:aa; */
                
                    grid->f_deforest = grid->f_crop_con - grid->f_crop_p;
                }else{
                    grid->f_deforest = grid->f_deforest_v = grid->f_deforest_s = 0.0;
                }
            }
        }
        
        /* Forest management: 2019/10/15 by A.Ito *********/
        if(EX_FORMAN == 1){
            grid->f_deforest = 0.0;
        }
	}
    
    /* biofuel experiment: 2015/08/27 by A.Ito */
    if(LANDUSE == 17 || BIOFUEL_RUN >= 1){
        grid->f_deforest = 0.0;
    }
    
    /* fixed land-use for NMIP */
    /* updated 2016/10/20 by A.Ito */
    if(NMIP_RUN == 5 || NMIP_RUN == 6 || NMIP_RUN == 7){
        grid->f_deforest = 0.0;
    }

    /* forced afforestation: 2020/11/09 by A.Ito */
    if((EX_FORCED_AFFOREST_2 >=1 && EX_FORCED_AFFOREST_2 <=12) && grid->lucy >= EX_FORCED_AFFOREST_2_YR){
        grid->f_deforest = 0.0;
        grid->f_deforest_v = 0.0;
        grid->f_deforest_s = 0.0;
    }

    /* parameter ensemble: 2014/11/19 by A.Ito */
    prm_ensen = 1.0;
    if(PARAM_PTB == 7){
        if(PARAM_ENS == 1){
            prm_ensen *= 0.7;
        }
        if(PARAM_ENS == 2){
            prm_ensen *= 0.8;
        }
        if(PARAM_ENS == 3){
            prm_ensen *= 0.9;
        }
        if(PARAM_ENS == 4){
            prm_ensen *= 1.1;
        }
        if(PARAM_ENS == 5){
            prm_ensen *= 1.2;
        }
        if(PARAM_ENS == 6){
            prm_ensen *= 1.3;
        }
        grid->f_deforest *= prm_ensen;
        grid->f_deforest_v *= prm_ensen;
        grid->f_deforest_s *= prm_ensen;
    }
	
    /* C-budget parameter ensemble: 2018/06/05 by A.Ito */
    if(PARAM_PTB == 20){
        prm_ensen = 1.0 + 0.3 * f_pert[6];

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
			if(grid->f_paddy_b > 0.0 && grid->fcrop_rk[2000 - FDY_LUC] > 0.0){
				grid->f_paddy = grid->f_paddy_b * 
					(grid->fcrop_rk[grid->lucy - FDY_LUC] / grid->fcrop_rk[2000 - FDY_LUC]);
			}
		}else if(grid->lucy >= 2006){
			if(grid->f_paddy_b > 0.0 && grid->fcrop_rk[2000 - FDY_LUC] > 0.0){
				grid->f_paddy = grid->f_paddy_b * 
					(grid->fcrop_rk[2005 - FDY_LUC] / grid->fcrop_rk[2000 - FDY_LUC]);
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
			if(grid->f_paddy_b > 0.0 && grid->fcrop_luh[2000 - FDY_LUC] > 0.0){
				grid->f_paddy = grid->f_paddy_b * 
                    (grid->fcrop_luh[grid->lucy - FDY_LUC] / grid->fcrop_luh[2000 - FDY_LUC]);
			}
		}else if(grid->lucy >= 2006){
			if(grid->f_paddy_b > 0.0 && grid->fcrop_luh[2000 - FDY_LUC] > 0.0){
				grid->f_paddy = grid->f_paddy_b * 
                    (grid->fcrop_luh[2005 - FDY_LUC] / grid->fcrop_luh[2000 - FDY_LUC]);
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
            || LANDUSE==14 || LANDUSE==15 || LANDUSE==16 || LANDUSE==17
            || LANDUSE==26 || LANDUSE==27 || LANDUSE==28|| LANDUSE == 29
            || LANDUSE == 30 || LANDUSE == 31 || LANDUSE == 32 || LANDUSE == 33
            || LANDUSE == 34 || LANDUSE == 35 || LANDUSE == 36 || LANDUSE == 37){
        
        if((grid->lucy - FDY_LUC) < DL_LUC){
            if(grid->f_paddy_b > 0.0 && grid->fcrop_luh[2000 - FDY_LUC] > 0.0){
                grid->f_paddy = grid->f_paddy_b *
                    (grid->fcrop_luh[grid->lucy - FDY_LUC] / grid->fcrop_luh[2000 - FDY_LUC]);
            }else{
                grid->f_paddy = 0.0;
            }
        }else{
            if(grid->f_paddy_b > 0.0 && grid->fcrop_luh[2000 - FDY_LUC] > 0.0){
                grid->f_paddy = grid->f_paddy_b *
                    (grid->fcrop_luh[DL_LUC - 1] / grid->fcrop_luh[2000 - FDY_LUC]);
            }else{
                grid->f_paddy = 0.0;
            }
        }
        
        if(grid->f_paddy > 1.0){
			grid->f_paddy = 1.0;
		}
		if(grid->f_paddy < 0.0){
			grid->f_paddy = 0.0;
		}
    }else if(LANDUSE == 18 || LANDUSE == 19 || LANDUSE == 20
            || LANDUSE == 21 || LANDUSE == 22 || LANDUSE == 23
            || LANDUSE == 24 || LANDUSE == 25 || LANDUSE == 45){
        grid->f_paddy = grid->f_paddy_b;
    }else if(LANDUSE == 38 || LANDUSE == 39 || LANDUSE == 40 || LANDUSE == 41){
        /* AIM land use, 2019/06/21 by A.Ito */
        grid->f_paddy = grid->f_paddy_b;
    }else if(LANDUSE == 46){
        /* ISIMIP3b: 2020/11/18 by A.Ito */
        grid->f_paddy = grid->f_paddy_b;
    }
    
    /* biofuel experiment: 2015/08/27 by A.Ito */
    if(LANDUSE == 17 || BIOFUEL_RUN >= 1){
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
	double fe_conv;		/* fraction of conversion flux (1year) */
	double fe_ten;		/* fraction of 10-year pool flux */
	double fe_hund;		/* fraction of 100-year pool flux */
	double fe_detr;		/* fraction of detritus flux */
	double mass_detr, mass_conv, mass_ten, mass_hund;	/* added by A.Ito based on E.Kato (2009/03/30) */
	double eff_mass, f_mass_secfor;
    double fluc_1 = 0.0, fluc_10 = 0.0, fluc_100 = 0.0;
	
	/* mass fraction of secondary forest */
	/* Nelson et al. (2000) Bioscience 50:419-431 */
	/* f_mass_secfor = 0.1; */
    f_mass_secfor = 0.3;

	switch(grid->veg_sage){
		/* detritus production by land-use change:
		 McGuire, A. D., et al., 2001: Carbon balance of the terrestrial biosphere 
         in the twentieth century: analysis of CO2, climate and land use effects 
         with four process-based ecosystem models.
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
	
	/* effective biomass */
	eff_mass = (mass->plant).fol + (mass->plant).stm + 0.8*(mass->plant).rot;
	
    /* added: A. Ito (with Hamada-san's comment) 2012/01/30 */
    grid->f_luc = 0.0;
    
	if(grid->phase == 0){
		/* spin-up: fluxes for 1801-1900 *******************************/
		/* modified by A.Ito (2009/06/05: 2010/01/07) */
		if(LANDUSE == 0 || LANDUSE==9){
			fluc_1 = 0.0;
		}else if(LANDUSE >= 1 && LANDUSE <= 5){
			fluc_1 = grid->fcrop_net[BGY_LUC-FDY_LUC] - grid->fcrop_net[BGY_LUC-FDY_LUC-1];
		}else if(LANDUSE == 6 || LANDUSE==8 || LANDUSE==10 || LANDUSE==11 || LANDUSE == 12
                || LANDUSE == 13 || LANDUSE==14 || LANDUSE==15 || LANDUSE==16 || LANDUSE == 17
                || LANDUSE == 26 || LANDUSE == 27 || LANDUSE == 28 || LANDUSE == 29
                || LANDUSE == 30 || LANDUSE == 31 || LANDUSE == 32 || LANDUSE == 33
                || LANDUSE == 34 || LANDUSE == 35 || LANDUSE == 36 || LANDUSE == 37){
            /* revised (after comments by E.Kato): 2013/10/02 by A.Ito */
			fluc_1 = (grid->t_vc_luh[BGY_LUC - FDY_LUC] + grid->t_vp_luh[BGY_LUC - FDY_LUC])
				+ (grid->t_sc_luh[BGY_LUC - FDY_LUC] + grid->t_sp_luh[BGY_LUC - FDY_LUC])
                * f_mass_secfor;
		}else if(LANDUSE == 7){
			/* added 2010/01/07 (A.Ito) */
			fluc_1 = (grid->fcrop_rk[BGY_LUC-FDY_LUC] - grid->fcrop_rk[BGY_LUC-FDY_LUC-1])
					+(grid->fpast_rk[BGY_LUC-FDY_LUC] - grid->fpast_rk[BGY_LUC-FDY_LUC-1]);
		}else if(LANDUSE == 18|| LANDUSE == 19 || LANDUSE == 20
                || LANDUSE == 21 || LANDUSE == 22 || LANDUSE == 23
                || LANDUSE == 24 || LANDUSE == 25 || LANDUSE == 45 || LANDUSE == 46){
			fluc_1 = 0.0;
		}else if(LANDUSE == 38 || LANDUSE == 39 || LANDUSE == 40 || LANDUSE == 41){
            /* AIM land use, 2019/06/21 by A.Ito */
            fluc_1 = grid->fcrop_luh[1950 - FDY_LUC] - grid->fcrop_luh[1949 - FDY_LUC];
        }
		
        /* NMIP: fixed land-use */
        /* updated 2016/10/20 by A.Ito */
        if(NMIP_RUN == 5 || NMIP_RUN == 6 || NMIP_RUN == 7){
            fluc_1 = 0.0;
        }

		/* modified by A.Ito (2009/08/19) */
		if(fluc_1 > 0.0){
			flux->lu_detr = fluc_1 * 0.2*(mass->plant).rot;
			flux->lu_conv = fluc_1 * eff_mass * fe_conv/(fe_conv + fe_ten + fe_hund);
		}else{
			flux->lu_detr = 0.0;
			flux->lu_conv = 0.0;
		}
		
		for(f=(BGY_LUC-9);f<=BGY_LUC;f++){
			/* senstivity analysis */
			if(LANDUSE == 0 || LANDUSE == 9){
				fluc_10 = 0.0;
			}else if(LANDUSE>=1 && LANDUSE<=5){
				fluc_10 = grid->fcrop_net[f - FDY_LUC] - grid->fcrop_net[f - FDY_LUC - 1];
			}else if(LANDUSE == 6 || LANDUSE == 8 || LANDUSE == 10 || LANDUSE == 11 || LANDUSE == 12
                    || LANDUSE==13 || LANDUSE==14 || LANDUSE==15 || LANDUSE==16 || LANDUSE==17
                    || LANDUSE==26 || LANDUSE==27 || LANDUSE==28 || LANDUSE == 29
                    || LANDUSE == 30 || LANDUSE == 31 || LANDUSE == 32 || LANDUSE == 33
                    || LANDUSE == 34 || LANDUSE == 35 || LANDUSE == 36 || LANDUSE == 37){
				fluc_10 = (grid->t_vc_luh[f - FDY_LUC] + grid->t_vp_luh[f - FDY_LUC])
						+ (grid->t_sc_luh[f - FDY_LUC] + grid->t_sp_luh[f - FDY_LUC]) * f_mass_secfor;
				/* 0.5: assumption by A.Ito for secondary forest stock */
			}else if(LANDUSE == 7){
				/* added 2010/01/07 (A.Ito) */
				fluc_10 = (grid->fcrop_rk[f - FDY_LUC] - grid->fcrop_rk[f-FDY_LUC-1])
						+ (grid->fpast_rk[f - FDY_LUC] - grid->fpast_rk[f-FDY_LUC-1]);
			}else if(LANDUSE == 18 || LANDUSE == 19 || LANDUSE == 20 ||
                LANDUSE == 21 || LANDUSE == 22 || LANDUSE == 23 || LANDUSE == 24 ||
                LANDUSE == 25 || LANDUSE == 45 || LANDUSE == 46){
                fluc_10 = 0.0;
            }else if(LANDUSE == 38 || LANDUSE == 39 || LANDUSE == 40 || LANDUSE == 41){
                /* AIM land use, 2019/06/21 by A.Ito */
                fluc_10 = grid->fcrop_luh[f - FDY_LUC] - grid->fcrop_luh[f - FDY_LUC -1];
            }
			
            /* NMIP: fixed land-use */
            /* updated 2016/10/20 by A.Ito */
            if(NMIP_RUN == 5 || NMIP_RUN == 6 || NMIP_RUN == 7){
                fluc_10 = 0.0;
            }

			/* modified by A.Ito based on E.Kato (2009/03/30) */
			if(fluc_10 > 0.0){
				mass_ten = fluc_10 * eff_mass * fe_ten/(fe_conv + fe_ten + fe_hund);
				flux->detr_ten[BGY_LUC-f] = mass_ten;
			}else{
                /* corrected: A. Ito (with Hamada-san's comment) 2012/01/30 */
				mass_ten = 0.0;
				flux->detr_ten[BGY_LUC-f] = 0.0;
			}
			/* corrected: A.Ito and E.Kato (2009/08/16) */
			flux->lu_ten += mass_ten * 0.1;
		}
		
		for(f=(BGY_LUC-99);f<=BGY_LUC;f++){
			/* senstivity analysis */
			if(LANDUSE == 0 || LANDUSE == 9){
				fluc_100 = 0.0;
			}else if(LANDUSE >= 1 && LANDUSE <= 5){
				fluc_100 = grid->fcrop_net[f - FDY_LUC] - grid->fcrop_net[f - FDY_LUC-1];
			}else if(LANDUSE == 6 || LANDUSE == 8 || LANDUSE == 10 || LANDUSE == 11 || LANDUSE == 12
                    || LANDUSE==13 || LANDUSE==14 || LANDUSE==15 || LANDUSE==16 || LANDUSE==17
                    || LANDUSE==26 || LANDUSE==27 || LANDUSE==28 || LANDUSE == 29
                    || LANDUSE == 30 || LANDUSE == 31 || LANDUSE == 32 || LANDUSE == 33
                    || LANDUSE == 34 || LANDUSE == 35 || LANDUSE == 36 || LANDUSE == 37){
				fluc_100 = (grid->t_vc_luh[f - FDY_LUC] + grid->t_vp_luh[f - FDY_LUC])
						+ (grid->t_sc_luh[f - FDY_LUC] + grid->t_sp_luh[f - FDY_LUC])*f_mass_secfor;
			}else if(LANDUSE == 7){
				/* added 2010/01/07 (A.Ito) */
				fluc_100 = (grid->fcrop_rk[f - FDY_LUC] - grid->fcrop_rk[f-FDY_LUC-1])
						+ (grid->fpast_rk[f - FDY_LUC] - grid->fpast_rk[f-FDY_LUC-1]);
			}else if(LANDUSE==18|| LANDUSE == 19 || LANDUSE == 20 ||
                LANDUSE == 21 || LANDUSE == 22 || LANDUSE == 23 ||
                LANDUSE == 24 || LANDUSE == 25 || LANDUSE == 45 || LANDUSE == 46){
                fluc_100 = 0.0;
            }else if(LANDUSE == 38 || LANDUSE == 39 || LANDUSE == 40 || LANDUSE == 41){
                /* AIM land use, 2019/06/21 by A.Ito */
                fluc_100 = grid->fcrop_luh[f - FDY_LUC] - grid->fcrop_luh[f - FDY_LUC -1];
            }
			
            /* NMIP: fixed land-use */
            /* updated 2016/10/20 by A.Ito */
            if(NMIP_RUN == 5 || NMIP_RUN == 6 || NMIP_RUN == 7){
                fluc_100 = 0.0;
            }

			/* modified by A.Ito based on E.Kato (2009/03/30) */
			/* corrected: A.Ito and E.Kato (2009/08/16) */
			if(fluc_100 > 0.0){
				mass_hund = fluc_100 * eff_mass * fe_hund/(fe_conv + fe_ten + fe_hund);
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
			fluc_1 = 0.0;
		}else if(LANDUSE >= 1 && LANDUSE <= 5){
			fluc_1 = grid->f_deforest;
			/*  grid->f_crop_con - grid->f_crop_p;  */
		}else if(LANDUSE==6 || LANDUSE==8 || LANDUSE==10 || LANDUSE==11 || LANDUSE==12
                || LANDUSE==13 || LANDUSE==14 || LANDUSE==15 || LANDUSE==16 || LANDUSE==17
                || LANDUSE==26 || LANDUSE==27 || LANDUSE==28 || LANDUSE == 29
                || LANDUSE == 30 || LANDUSE == 31 || LANDUSE == 32 || LANDUSE == 33
                || LANDUSE == 34 || LANDUSE == 35 || LANDUSE == 36 || LANDUSE == 37){
			/* assumption: biomass in secondary forest is lower (0.1) than primary forest */
			fluc_1 = grid->f_deforest_v + grid->f_deforest_s * f_mass_secfor;
		}else if(LANDUSE == 7){
			/* added 2010/01/07 (A.Ito) */
			fluc_1 = grid->f_deforest;
		}else if(LANDUSE == 18 || LANDUSE == 19 || LANDUSE == 20
                || LANDUSE == 21 || LANDUSE == 22 || LANDUSE == 23
                || LANDUSE == 24 || LANDUSE == 25 || LANDUSE == 45 || LANDUSE == 46){
            fluc_1 = grid->f_deforest;
        }else if(LANDUSE == 38 || LANDUSE == 39 || LANDUSE == 40 || LANDUSE == 41){
            /* AIM land use, 2019/06/21 by A.Ito */
            fluc_1 = grid->f_deforest;
        }
		
        /* NMIP: fixed land-use */
        /* updated 2016/10/20 by A.Ito */
        if(NMIP_RUN == 5 || NMIP_RUN == 6 || NMIP_RUN == 7){
            fluc_1 = 0.0;
        }
        
        /* BECCS S10 experiment: 2016/02/16 by A.Ito ******/
        if(EX_BECCS == 1 || EX_BECCS == 2|| EX_BECCS == 3){
            if(grid->lucy <= 2000){
                fluc_1 = 0.0;
            }else if(grid->lucy >= 2001){
                fluc_1 = grid->f_deforest;
            }
        }

		if(fluc_1 > 0.0){ /* deforested */
			/* modified by A.Ito based on E.Kato (2009/03/30) */
			mass_detr = fluc_1 * 0.2 * (mass->plant).rot;
			mass_conv = fluc_1 * eff_mass * fe_conv/(fe_conv + fe_ten + fe_hund);
			mass_ten = fluc_1 * eff_mass * fe_ten/(fe_conv + fe_ten + fe_hund);
			mass_hund = fluc_1 * eff_mass * fe_hund/(fe_conv + fe_ten + fe_hund);
			
			/* emission from 1-yr or instantaneous pool */
			flux->lu_detr = mass_detr;
			flux->lu_conv = mass_conv;
			
			/* corrected: A.Ito and E.Kato (2009/08/16) */
			flux->detr_ten[0] = mass_ten;
			flux->detr_hund[0] = mass_hund;
   
            /* carbon loss for each compartment: 2018/10/24 by A.Ito  */
   
            flux->lu_fol = fluc_1 * (mass->plant).fol * fe_conv/(fe_conv + fe_ten + fe_hund);
            flux->lu_stm = fluc_1 * (mass->plant).stm * fe_conv/(fe_conv + fe_ten + fe_hund);
            flux->lu_rot = fluc_1 * 0.8 * (mass->plant).rot * fe_conv/(fe_conv + fe_ten + fe_hund);
            flux->lu_ltr = fluc_1 * 0.2 * (mass->plant).rot;
            flux->lu_msl = 0.0;
   
		}else{
			flux->lu_detr = 0.0;
			flux->lu_conv = 0.0;
			flux->detr_ten[0] = 0.0;
			flux->detr_hund[0] = 0.0;

            flux->lu_fol = 0.0;
            flux->lu_stm = 0.0;
            flux->lu_rot = 0.0;
            flux->lu_ltr = 0.0;
            flux->lu_msl = 0.0;
		}
        
        grid->f_luc = fluc_1;
	}
}
