/*	VISIT: Vegetation Integrative SImulation Tool						*/
/*  Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed  in CGER/NIES & RIGC/JAMSTEC						*/
/*  Carbon cycle, erosion, biomass burning, land-use change,			*/
/*  CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

#define NOTICE 0

/* CH4 emission model by Cao **********************************************/
/*
Cao, M., Marshall, S. and Gregson, K., 1996. Global carbon exchange and 
methane emissions from natural wetlands: Application of a process-based model. 
Journal of Geophysical Research 101, 14399-14414.
*/
void f_ch4_emit_cao(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Flux *flux
){
	double fv_inund_wet, fv_inund_pad, wtable;	/* water table, cm */
	double fv_temp, fv_wtable, f_wtable_lake, fv_wetland, fv_lake;
	double hr_decomp, gpp_factor, diff_wtd, dtemp;
	
	/* soil decomposition rate. Mg C ha-1 month-1 */
	hr_decomp = (flux->soil).hr[grid->m];
	if(hr_decomp < 0.0){
		hr_decomp = 0.0;
	}
	
	/* temperature (deg C) coefficient */
	/* eq.7 */
    dtemp = 0.0;
    if(EX_PADDY == 4){
        dtemp = 1.0;
    }
	fv_temp = exp((grid->tmp10_soil[grid->m] + dtemp) * 0.0693) / 7.996;
	if(fv_temp < 0.0){
		fv_temp = 0.0;
	}
	
	/* CH4 emission, mg CH4 m-2 month-1 */
    
    if(ALT_INUND == 0){
        /* standard */
        fv_inund_wet = grid->inundation_ssmi[grid->m];
        fv_inund_pad = grid->inundation_ssmi[grid->m];
    }else if(ALT_INUND==1 || ALT_INUND==3 || ALT_INUND==4 || ALT_INUND==5){
        /* for GCP-CH4 analyses */
        if(grid->climy>=1999 && grid->climy<=2013){
        
            if(grid->f_wetland > 0.0){
                fv_inund_wet = grid->inundation_alt_ts[grid->climy - 1999][grid->m] / grid->f_wetland;
                
                if(fv_inund_wet > 1.0){
                    fv_inund_wet = 1.0;
                }
            }
            if(grid->f_paddy > 0.0){
                /* f_inund_pad = grid->inundation_alt_ts[grid->climy-1999][grid->m] / grid->f_paddy; */
                
                fv_inund_pad = grid->inundation_ssmi[grid->m];
                
                if(fv_inund_pad > 1.0){
                    fv_inund_pad = 1.0;
                }
            }
        
        }else{
            if(grid->f_wetland > 0.0){
                fv_inund_wet = grid->inundation_alt_av[grid->m] / grid->f_wetland;
                
                if(fv_inund_wet > 1.0){
                    fv_inund_wet = 1.0;
                }else{
                
                }
            }
            if(grid->f_paddy > 0.0){
                /* f_inund_pad = grid->inundation_alt_av[grid->m] / grid->f_paddy; */
                
                fv_inund_pad = grid->inundation_ssmi[grid->m];
                
                if(fv_inund_pad > 1.0){
                    fv_inund_pad = 1.0;
                }
            }
        }
    }else if(ALT_INUND == 6){
        /* for GCP-CH4 fw_swamps-glwd_2000-2012.txt */
        /* added: 2015/02/06  */
    
        if(grid->climy>=2000 && grid->climy<=2012){
        
            if(grid->f_wetland > 0.0){
                fv_inund_wet = grid->inundation_alt_ts[grid->climy - 2000][grid->m] / grid->f_wetland;
                
                if(fv_inund_wet > 1.0){
                    fv_inund_wet = 1.0;
                }
            }
            if(grid->f_paddy > 0.0){
                /* f_inund_pad = grid->inundation_alt_ts[grid->climy-2000][grid->m] / grid->f_paddy; */
                
                fv_inund_pad = grid->inundation_ssmi[grid->m];
                
                if(fv_inund_pad > 1.0){
                    fv_inund_pad = 1.0;
                }
            }
        
        }else{
            if(grid->f_wetland > 0.0){
                fv_inund_wet = grid->inundation_alt_av[grid->m] / grid->f_wetland;
                
                if(fv_inund_wet > 1.0){
                    fv_inund_wet = 1.0;
                }else{
                
                }
            }
            if(grid->f_paddy > 0.0){
                /* f_inund_pad = grid->inundation_alt_av[grid->m] / grid->f_paddy; */
                
                fv_inund_pad = grid->inundation_ssmi[grid->m];
                
                if(fv_inund_pad > 1.0){
                    fv_inund_pad = 1.0;
                }
            }
        }
    }else if(ALT_INUND == 7 || ALT_INUND == 8){
        /* for GCP-CH4 v2: fw_gcp-ch4_wetlands_2000-2017.txt */
        /* added: 2018/08/19  */
    
        if(grid->climy>=2000 && grid->climy<=2017){
        
            if(grid->f_wetland > 0.0){
                fv_inund_wet = grid->inundation_alt_ts[grid->climy - 2000][grid->m] / grid->f_wetland;
                
                if(fv_inund_wet > 1.0){
                    fv_inund_wet = 1.0;
                }
            }else{
                fv_inund_wet = 0.0;
            }
            if(grid->f_paddy > 0.0){
                /* f_inund_pad = grid->inundation_alt_ts[grid->climy-2000][grid->m] / grid->f_paddy; */
                
                fv_inund_pad = grid->inundation_ssmi[grid->m];
                
                if(fv_inund_pad > 1.0){
                    fv_inund_pad = 1.0;
                }
            }else{
                fv_inund_pad = 0.0;
            }
        
        }else{
            if(grid->f_wetland > 0.0){
                fv_inund_wet = grid->inundation_alt_av[grid->m] / grid->f_wetland;
                
                if(fv_inund_wet > 1.0){
                    fv_inund_wet = 1.0;
                }
            }else{
                fv_inund_wet = 0.0;
            }
            if(grid->f_paddy > 0.0){
                /* f_inund_pad = grid->inundation_alt_av[grid->m] / grid->f_paddy; */
                
                fv_inund_pad = grid->inundation_ssmi[grid->m];
                
                if(fv_inund_pad > 1.0){
                    fv_inund_pad = 1.0;
                }
            }else{
                fv_inund_pad = 0.0;
            }
        }
    }else if(ALT_INUND == 9){
        /* SWAMPS anomaly: 2021/06/25  */
        if(grid->climy>=1992 && grid->climy<=2020){
            fv_inund_wet = grid->inundation_ssmi[grid->m] + grid->inundation_alt_ts[grid->climy - 1992][grid->m];

            if(fv_inund_wet > 1.0){
                fv_inund_wet = 1.0;
            }
        }else{
            fv_inund_wet = grid->inundation_ssmi[grid->m];
        }
        
        fv_inund_pad = grid->inundation_ssmi[grid->m];
    }else if(ALT_INUND == 10){
        /* for WAD2M_wetlands_2000-2020_05deg_Ver2.0.flt */
        /* added: 2021/10/26  */
    
        if(grid->climy>=2000 && grid->climy<=2020){
        
            if(grid->f_wetland > 0.0){
                fv_inund_wet = grid->inundation_alt_ts[grid->climy - 2000][grid->m] / grid->f_wetland;
                
                if(fv_inund_wet > 1.0){
                    fv_inund_wet = 1.0;
                }
            }else{
                fv_inund_wet = 0.0;
            }
            if(grid->f_paddy > 0.0){
                /* f_inund_pad = grid->inundation_alt_ts[grid->climy-2000][grid->m] / grid->f_paddy; */
                
                fv_inund_pad = grid->inundation_ssmi[grid->m];
                
                if(fv_inund_pad > 1.0){
                    fv_inund_pad = 1.0;
                }
            }else{
                fv_inund_pad = 0.0;
            }
        
        }else{
            if(grid->f_wetland > 0.0){
                fv_inund_wet = grid->inundation_alt_av[grid->m] / grid->f_wetland;
                
                if(fv_inund_wet > 1.0){
                    fv_inund_wet = 1.0;
                }else{
                
                }
            }else{
                fv_inund_wet = 0.0;
            }
            if(grid->f_paddy > 0.0){
                /* f_inund_pad = grid->inundation_alt_av[grid->m] / grid->f_paddy; */
                
                fv_inund_pad = grid->inundation_ssmi[grid->m];
                
                if(fv_inund_pad > 1.0){
                    fv_inund_pad = 1.0;
                }
            }else{
                fv_inund_pad = 0.0;
            }
        }
    }else{
        fv_inund_wet = grid->inundation_ssmi[grid->m];
        fv_inund_pad = grid->inundation_ssmi[grid->m];
    }
    
    /* alternative paddy crop calendar: 2021/04/07  */
    if(EX_PADDY == 3){
        if(grid->iizumi_mon_paddy_end > grid->iizumi_mon_paddy_start){
            if(grid->m >= grid->iizumi_mon_paddy_start && grid->m <= grid->iizumi_mon_paddy_end){
                fv_inund_pad = 1.0;
            }else{
                fv_inund_pad = 0.0;
            }
        }else if(grid->iizumi_mon_paddy_end < grid->iizumi_mon_paddy_start){
            if(grid->m <= grid->iizumi_mon_paddy_start || grid->m >= grid->iizumi_mon_paddy_end){
                fv_inund_pad = 1.0;
            }else{
                fv_inund_pad = 0.0;
            }
        }else{
            fv_inund_pad = 0.0;
        }
    }

	/* water table (cm relative to surface) coefficient */
	/* wetland *****************************************************************************/
	/* eq.6 */
	/* wtable = 5.0; */
	/* f_wtable = 0.383 * (0.5*exp(0.096 * wtable) + 0.5*exp(0.096 * -5.0)); OLD */
	/* revised  (2009/07/13) **/
	/* wtable = 0.0;
	f_wtable = 0.383 * (grid->inundation_ssmi[grid->m]*exp(0.096 * wtable) 
						+ (1.0-grid->inundation_ssmi[grid->m])*exp(0.096 * -25.0)); */
	/* 2013/11/29  */
 	/* 2015/11/30  */
   wtable = 4.5;
    
    /* 2014/12/10  */
    if(VAR_WTD == 1){
        /* mm => cm */
        diff_wtd = ((loct->sw30 + loct->sww) - (loct->b_sw30[grid->m] + loct->b_sww[grid->m])) / 10.0;
        wtable += diff_wtd;
        if(wtable < -30.0){
            wtable = -30.0;
        }
        if(wtable > 12.0){
            wtable = 12.0;
        }
    }
    
	fv_wtable = 0.383 * (fv_inund_wet * exp(0.096 * wtable)
						+ (1.0 - fv_inund_wet) * exp(0.096 * (wtable - 12.0)));
	/* if(ALT_FWETLAND == 1){
		f_wtable = 0.383 * (grid->f_wetland*exp(0.096 * wtable) 
							+ (1.0-grid->f_wetland)*exp(0.096 * -25.0));
	} */
	if(fv_wtable < 0.0){
		fv_wtable = 0.0;
	}
     
	/* lake: added  (2009/07/14) */
	/* f_wtable_lake = 0.383 * exp(0.096 * 5.0); */
	/* lake: revised  (2013/11/29) */
	f_wtable_lake = 0.383 * exp(0.096 * 10.0);
 
    if(ALT_INUND == 8){
        if(grid->climy>=2000 && grid->climy<=2017){
            fv_wetland = grid->inundation_alt_ts[grid->climy - 2000][grid->m];
        }else{
            fv_wetland = grid->inundation_alt_av[grid->m];
        }
        if(fv_wetland < grid->f_wetland){
            fv_wetland = grid->f_wetland;
        }
        fv_lake = grid->f_lake;
    }else if(ALT_INUND == 9){
        if(grid->climy>=1992 && grid->climy<=2020){
            fv_wetland = grid->f_wetland + grid->inundation_alt_ts[grid->climy - 1992][grid->m];
            if(fv_wetland > 1.0){
                fv_wetland = 1.0;
            }
        }else{
            fv_wetland = grid->f_wetland;
        }
        fv_lake = grid->f_lake;
    }else if(ALT_INUND == 10){
        if(grid->climy>=2000 && grid->climy<=2020){
            fv_wetland = grid->inundation_alt_ts[grid->climy - 2000][grid->m];
        }else{
            fv_wetland = grid->inundation_alt_av[grid->m];
        }
        if(fv_wetland < grid->f_wetland){
            fv_wetland = grid->f_wetland;
        }
        fv_lake = grid->f_lake;
    }else{
        fv_wetland = grid->f_wetland;
        fv_lake = grid->f_lake;
    }
	
	/* Mg C ha-1 month-1 */
	(flux->soil).ch4prod_wetland_cao[grid->m] = hr_decomp * fv_temp *
					(fv_wtable * fv_wetland + f_wtable_lake * 0.2 * fv_lake); /* 0.2: 090717 */
	if(ALT_FWETLAND == 1){
		(flux->soil).ch4prod_wetland_cao[grid->m] = hr_decomp * fv_temp *
			        (fv_wtable * fv_wetland + f_wtable_lake * 0.2 * fv_lake);
	}
	
	/* mg CH4 m-2 month-1 */
	(flux->soil).ch4prod_wetland_cao[grid->m] *= 16.0/12.0 * 1000000000.0 / 10000.0;
	
	/* CH4 oxidation */
	if((flux->soil).ch4prod_wetland_cao[grid->m]>0.0 && loct->gpp_max>0.0){
		gpp_factor = (flux->plant).gpp[grid->m] / loct->gpp_max;
		gpp_factor = (gpp_factor>0.0)?gpp_factor:0.0;
		gpp_factor = (gpp_factor<1.0)?gpp_factor:1.0;
		
		(flux->soil).ch4oxy_wetland_cao[grid->m] = (flux->soil).ch4prod_wetland_cao[grid->m] *
            (0.60 + 0.30 * gpp_factor);
	}else{
		(flux->soil).ch4oxy_wetland_cao[grid->m] = 0.0;
	}
	
	/* net flux */
	(flux->soil).ch4flux_wetland_cao[grid->m] = (flux->soil).ch4prod_wetland_cao[grid->m] 
			- (flux->soil).ch4oxy_wetland_cao[grid->m];
	
    /* if((flux->soil).ch4flux_wetland_cao[grid->m] > 1000.0){
        printf("*****************************%ld %lf %lf\n", grid->n_olson, 
        (flux->soil).ch4flux_wetland_cao[grid->m], grid->tmp10_soil[grid->m]);
        exit(1);
    }*/
    
	/* paddy field *****************************************************************************/
	/* eq.6 */
	/* if(grid->tmp_2m[grid->m] > 15.0 && grid->prate_sfc[grid->m] > 50.0){
		wtable = 0.0;
	}else{
		wtable = -20.0;
	}
	f_wtable = 0.383 * exp(0.096 * wtable); */
	/* revised  (2009/07/13) */
	/* f_wtable = 0.383 * (f_inund * exp(0.096 * 3.0)
						+ (1.0 - f_inund)*exp(0.096 * -50.0));	 */
    /* revised  (2013/11/29) */

    /* 2014/12/10  */
    wtable = 4.0;
    
    if(EX_PADDY == 1){
        wtable = -2.0;
        //wtable = -8.0;
    }

	/* fv_wtable = 0.383 * (fv_inund_pad * exp(0.096 * wtable)
						+ (1.0 - fv_inund_pad)*exp(0.096 * (wtable - 14.0))); */
    fv_wtable = 0.383 * (fv_inund_pad * exp(0.096 * wtable)
                        + (1.0 - fv_inund_pad)*exp(0.096 * (wtable - 50.0)));

	if(fv_wtable < 0.0){
        fv_wtable = 0.0;
	}
	/* Mg C ha-1 month-1 */
	/* assuming low decomposition at paddy field: 0.4, 2008/06/11 */
	(flux->soil).ch4prod_paddy_cao[grid->m] = 0.5*hr_decomp * fv_temp * fv_wtable;
	/* mg CH4 m-2 month-1 */
	(flux->soil).ch4prod_paddy_cao[grid->m] *= 16.0/12.0 * 1000000000.0 / 10000.0;
	(flux->soil).ch4prod_paddy_cao[grid->m] *= grid->f_paddy;
	
	/* CH4 oxidation */
	if((flux->soil).ch4prod_paddy_cao[grid->m]>0.0 && loct->gpp_max>0.0){
		gpp_factor = (flux->plant).gpp[grid->m]/loct->gpp_max;
		gpp_factor = (gpp_factor>0.0)?gpp_factor:0.0;
		gpp_factor = (gpp_factor<1.0)?gpp_factor:1.0;

		(flux->soil).ch4oxy_paddy_cao[grid->m] = (flux->soil).ch4prod_paddy_cao[grid->m] *
			(0.60 + 0.30*gpp_factor);
	}else{
		(flux->soil).ch4oxy_paddy_cao[grid->m] = 0.0;
	}
	
	/* net flux */
	(flux->soil).ch4flux_paddy_cao[grid->m] = (flux->soil).ch4prod_paddy_cao[grid->m] 
			- (flux->soil).ch4oxy_paddy_cao[grid->m];
    
    /* 2014/01/22 for WSL output **********/
    /* loct->xx1[grid->m] = grid->f_wetland;
    loct->xx2[grid->m] = grid->f_wetland + 0.2*grid->f_lake + grid->f_paddy;
    loct->xx3[grid->m] = grid->inundation_ssmi[grid->m];
    loct->xx4[grid->m] = grid->inundation_ssmi[grid->m]*(-4.0) + (1.0-grid->inundation_ssmi[grid->m])*10.0;
    loct->xx5[grid->m] = 0.333*grid->tmp10_soil[grid->m] + 0.667*grid->tmp200_soil[grid->m]; */
}

/* CH4 emission by Walter & Heimann: added  (2009/08/05) *****************/
/*
	Walter, B. P., M. Heimann, et al. (2001). "Modeling modern methane emissions 
	from natural wetlands 1. Model description and results." 
	Journal of Geophysical Research 106(D24): 34189-34206.
*/
void f_ch4_emit_walter(
	short smode, 
	struct Grid *grid, 
	struct Loct *loct, 
	struct Flux *flux
){
	long f, g, cumtime;
	double ww[N_SLAYER+2];
	double df[N_SLAYER+2], dpth[N_SLAYER+2], t_mean;
	double ff[N_SLAYER+2], q_ebull[N_SLAYER+2], q_plant[N_SLAYER+2];
	double q_prod[N_SLAYER+2], q_oxid[N_SLAYER+2];
	double sdepth, wtdepth, b_thresh, tmp[N_SLAYER+2], rdepth, poro;
	double f_in, f_org[N_SLAYER+2], f_t[N_SLAYER+2], f_grow, t_gr, t_mat;
	double t_veg, flux_ebull, flux_plant, release, f_sand, f_clay;
	double hh, rr, kk, df_dry, fa_paddy, fa_wetland, day_produc, day_oxid;
	double r0, f_inundation, diff_wtd;
	double q10_ch4prod, dtemp;
    double efflux_ebul, efflux_plant, efflux_diffs, efflux_reles;
	
	/*
	EX_CH4_1
	 0: control
	 1: increased innundation: low sensitivity
	 2: increased innundation: high sensitivity
	*/
			
	/*
	 EX_CH4_2
	 0: control
	 1: decreased temperature dependence (Q10, CH4 production)
	 2: increased temperature dependence (Q10, CH4 production)
	 */
	
	/*
	 EX_CH4_3
	 0: control
	 1: increased wetland area by lake
	 */
    
    dtemp = 0.0;
    if(EX_PADDY == 4){
        dtemp = 1.0;
    }
	
	/************************************************************************/
	sdepth = 1.0;		/* soil depth, m */
	kk = 0.01;			/* time step, hour */
	hh = sdepth/N_SLAYER;	/* layer thickness, m */
	
    /* bubbling threshold, micro mol / liter */
	b_thresh = 500.0;
	
	if(SENS_WHCH4 == 3){
		b_thresh = 450.0;
	}else if(SENS_WHCH4 == 4){
		b_thresh = 550.0;
	}
	
	/* rdepth = loct->rtdpth; */		/* rooting depth, m */
	f_sand = grid->pc_sand / 100.0;
	f_clay = grid->pc_clay / 100.0;
	
	/* depth from the soil surface, m */
	dpth[0] = -0.04;
	for(f=1;f<=N_SLAYER;f++){
		dpth[f] = 0.0 + ((double)f - 0.5)*hh;
	}
	
	/* t_mean = grid->tmp_soil_mean; */
    t_mean = grid->tmp_soil_am;  /* lont-term mean soil temperature */
    /* if(FIX_STMP == 1){
        t_mean = grid->tmp_soil_am;
    } */
	
    /* vegetation factor: 2014/06/09  *********************************/
    /* last calibrated 2014/06/12 */
    /* last calibrated 2014/11/17 */
    /* last calibrated 2014/11/21 */
    /* last calibrated 2014/11/30 */
    /* last calibrated 2015/03/30 */
    /* last calibrated 2015/11/19 */
    /* last calibrated 2018/05/25 */
    /* last calibrated 2018/05/26 */
    /* last calibrated 2018/07/03 */
    /* last calibrated 2018/07/28 */
    switch(smode){
        case 1:
            if(grid->veg_olson==1 || grid->veg_olson==2 || grid->veg_olson==3){
                t_veg = 18.35;
                r0 = 1.835;
            }else if(grid->veg_olson==4 || grid->veg_olson==5 || grid->veg_olson==6){
                t_veg = 15.0;
                r0 = 1.50;
            }else if(grid->veg_olson==7 || grid->veg_olson==8){
                t_veg = 5.75;
                r0 = 0.575;
            }else if(grid->veg_olson==9 || grid->veg_olson==10){
                t_veg = 2.75;
                r0 = 0.29;
            }else if(grid->veg_olson==11 || grid->veg_olson==12){
                t_veg = 1.6;
                r0 = 0.39; /* 2019/04/17 */
            }else if(grid->veg_olson==21 || grid->veg_olson==22){ /* 2019/04/03  */
                t_veg = 2.0;
                r0 = 0.47; /* 2019/04/17 */
            }else if(grid->veg_olson==23){ /* 2019/04/03  */
                t_veg = 1.75;
                r0 = 0.175;
            }else if(grid->veg_olson==24){ /* 2019/04/03  */
                t_veg = 2.2;
                r0 = 0.46; /* 2019/04/17 */
            }else{
                t_veg = 1.45;
                r0 = 0.145;
            }
            break;
        case 2:
            if(grid->veg_olson==1 || grid->veg_olson==2 || grid->veg_olson==3){
                t_veg = 7.65;
                r0 = 0.765;
            }else if(grid->veg_olson==4 || grid->veg_olson==5 || grid->veg_olson==6){
                t_veg = 5.8;
                r0 = 0.58;
            }else if(grid->veg_olson==7 || grid->veg_olson==8){
                t_veg = 1.2;
                r0 = 0.12;
            }else if(grid->veg_olson==9 || grid->veg_olson==10){
                t_veg = 0.475;
                r0 = 0.049;
            }else if(grid->veg_olson==11 || grid->veg_olson==12){
                t_veg = 0.40;
                r0 = 0.074; /* 2019/04/17 */
            }else if(grid->veg_olson==21 || grid->veg_olson==22){ /* 2019/04/03  */
                t_veg = 0.40;
                r0 = 0.085; /* 2019/04/17 */
            }else if(grid->veg_olson==23){ /* 2019/04/03  */
                t_veg = 0.38;
                r0 = 0.038;
            }else if(grid->veg_olson==24){ /* 2019/04/03  */
                t_veg = 0.45;
                r0 = 0.085; /* 2019/04/17 */
            }else{
                t_veg = 0.43;
                r0 = 0.043;
            }
            break;
        case 3: /*  */
            t_veg = 7.4;
            r0 = 0.7;
            break;
        case 4:
            t_veg = 4.4;
            r0 = 0.4;
            break;
        default:
             t_veg = 1.0;
             r0 = 0.1;
    }
		
	/* TIME */
	cumtime = 600;
	
    /* C-budget parameter ensemble: 2018/06/05  */
    if(PARAM_PTB == 20){
        t_veg *= 1.0 + 0.3 * f_pert[2];
    }
    
	/* characteristics ***************************************/
    wtdepth = 0.0;
    rdepth = 0.1;
	if(smode == 1){	/* water-logged wetlands */
		/* t_veg = 6.0; */  /* vegetation factor */
		/* rdepth = 0.20; */		/* rooting depth, m */
		/* rdepth = 0.30; */		/* rooting depth, m */ /* revised 2013/11/29  */
        rdepth = 0.25;        /* rooting depth, m */ /* revised 2018/05/29  */

        if(grid->veg_olson==9 || grid->veg_olson==10){
            rdepth = 0.20;
        }else if(grid->veg_olson==11 || grid->veg_olson==12 || grid->veg_olson==21 ||
            grid->veg_olson==22 || grid->veg_olson==23){
            rdepth = 0.15;
        }
		
		/* water-table depth, m from surface */
		/* loct->water_table_depth = 0.00; */
		/* loct->water_table_depth = -0.02; */
		/* loct->water_table_depth = -0.02; */ /* revised 2013/11/29  */
		loct->water_table_depth = 0.025; /* revised 2015/02/23  */
		if(EX_CH4_1 == 1){
			loct->water_table_depth = 0.0 - loct->cum_dprec*0.0002;
		}else if(EX_CH4_1 == 2){
			loct->water_table_depth = 0.0 - loct->cum_dprec*0.001;
		}
        /* 2014/12/08  */
        if(VAR_WTD == 1){
            diff_wtd = (loct->sw30 + loct->sww) - (loct->b_sw30[grid->m] + loct->b_sww[grid->m]);
            loct->water_table_depth -= diff_wtd/1000.0;
            if(loct->water_table_depth > 0.3){
                loct->water_table_depth = 0.3;
            }
        }else{
            if(loct->water_table_depth > 0.0){
                loct->water_table_depth = 0.0;
            }
        }
		wtdepth = loct->water_table_depth;
        /* loct->xx1[grid->m] = wtdepth; */
		
		/* tuning parameter (cf. Table 2) */
		/* r0 = 0.4; */  /* 1.0 => 0.7: 2009/08/20 */
		/* r0 = 0.45; */  /* revised 2013/11/29  */
	}else if(smode == 2){	/* drainage wetlands */
		/* t_veg = 4.0; */ /* vegetation factor */
		/* rdepth = 0.15; */  /* rooting depth, m */
		/* rdepth = 0.10; */		/* rooting depth, m */ /* revised 2013/11/29  */
		rdepth = 0.15;		/* rooting depth, m */ /* revised 2015/11/19  */
		
        /* revised 2015/04/02  */
        if(grid->veg_olson==9 || grid->veg_olson==10){
            rdepth = 0.10;
        }else if(grid->veg_olson==11 || grid->veg_olson==12 || grid->veg_olson==21 || grid->veg_olson==22
             || grid->veg_olson==23){
            rdepth = 0.05;
        }

		/* water-table depth, m from surface */
		/* loct->water_table_depth = 0.25; */
		/* loct->water_table_depth = 0.20; */
		loct->water_table_depth = 0.20;  /* revised 2015/04/02  */
		if(EX_CH4_1 == 1){
			loct->water_table_depth = 0.25 - loct->cum_dprec * 0.0002;
		}else if(EX_CH4_1 == 2){
			loct->water_table_depth = 0.25 - loct->cum_dprec * 0.001;
		}
        /* 2014/12/08  */
        if(VAR_WTD == 1){
            diff_wtd = (loct->sw30 + loct->sww) - (loct->b_sw30[grid->m] + loct->b_sww[grid->m]);
            loct->water_table_depth -= diff_wtd/1000.0;
            if(loct->water_table_depth > 0.5){
                loct->water_table_depth = 0.5;
            }
        }else{
            if(loct->water_table_depth > 0.3){
                loct->water_table_depth = 0.3;
            }
        }
		wtdepth = loct->water_table_depth;
        /* loct->xx2[grid->m] = wtdepth; */
		
		/* tuning parameter (cf. Table 2) */
		/* r0 = 0.25; */  /* 1.0 => 0.7: 2009/08/20 */
		/* r0 = 0.35; */  /* revised 2013/11/29  */
	}else if(smode == 3){	/* water-logged paddy fields */
		/* t_veg = 10.0; */  /* vegetation factor */
		rdepth = 0.20;		/* rooting depth, m */ /* 0.3 => 0.2: 2009/08/20 */
		
		/* water-table depth, m from surface */
		loct->water_table_depth = -0.03;
        if(EX_PADDY == 1){
            //loct->water_table_depth = 0.15;
            loct->water_table_depth = 0.03;
            //loct->water_table_depth = 0.07;
        }
		wtdepth = loct->water_table_depth;
        /* loct->xx3[grid->m] = wtdepth; */
		
		/* tuning parameter (cf. Table 2) */
		/* r0 = 0.6;  */  /* 1.0 => 0.7: 2009/08/20 */
	}else if(smode == 4){	/* drainage paddy fields */
		/* t_veg = 10.0; */  /* vegetation factor */
		rdepth = 0.1;		/* rooting depth, m */
		
		/* water-table depth, m from surface */
		loct->water_table_depth = 0.5;
        if(EX_PADDY == 1){
            //loct->water_table_depth = 0.9;
            loct->water_table_depth = 0.56;
            //loct->water_table_depth = 0.6;
        }
		wtdepth = loct->water_table_depth;
        /*  loct->xx4[grid->m] = wtdepth; */
		
		/* tuning parameter (cf. Table 2) */
		/* r0 = 0.4;  */  /* 1.0 => 0.7: 2009/08/20 */
 	}
	
	/* porosity */
	poro = 0.45*f_sand + 0.20*(1.0 - f_sand - f_clay) + 0.14*f_clay + 0.45*0.03;	
	
	/* diffusion coefficient, m2 s-1 */
    df_dry = 0.2 * 0.66 * poro;
	for(f=1;f<=N_SLAYER;f++){
		if(dpth[f] < wtdepth){
            /* air-dominated */
			df[f] = 0.2 * 0.66 * poro;  /* pow(poro, 2.3) */		
		}else{
            /* water-dominated */
			df[f] = 0.00002 * 0.66 * poro;  
		}
	}
    df[0] = df[N_SLAYER+1] = df_dry;
    
    /********************/
	if(SENS_WHCH4 == 1){
		t_veg *= 1.1;
	}else if(SENS_WHCH4 == 2){
		t_veg = 0.9;
	}
		
	/* stability index, should be < 0.5 */
	rr = kk/hh/hh;
	if(rr*df[1] > 0.5){
		printf("BAD calculation condition: r = %lf\n", rr);
		exit(1);
	}
	
	/* plant growth state parameters **************/
	/* W&H2000 pp.763 */
	if(t_mean < 5.0){
		t_gr = 2.0;
	}else{
		t_gr = 7.0;
	}
	t_mat = t_gr + 10.0;
	for(f=1;f<=N_SLAYER;f++){
        if(FIX_STMP == 1){
            tmp[f] = grid->tmp10_soil_a[grid->m] * (double)(N_SLAYER - f)/(double)N_SLAYER
                    + (double)f/(double)N_SLAYER * grid->tmp200_soil_a[grid->m] +dtemp;
        }else{
            tmp[f] = grid->tmp10_soil[grid->m] * (double)(N_SLAYER - f)/(double)N_SLAYER
                    + (double)f/(double)N_SLAYER * grid->tmp200_soil[grid->m] +dtemp;
        }
	}
	/* fgow: Eq. 20 */
	if(tmp[5] < t_gr){
		/* f_grow = 0.0; */
        f_grow = 0.5;
	}else if(tmp[5] >= t_gr && tmp[5] <= t_mat){
		/* f_grow = 0.0 + 4.0 * (1.0 - pow((t_mat - tmp[5])/(t_mat - t_gr), 2.0)); */
        f_grow = 0.5 + 3.5 * (1.0 - pow((t_mat - tmp[5])/(t_mat - t_gr), 2.0));
	}else if(tmp[5] > t_mat){
		f_grow = 4.0;
	}else{
        f_grow = 1.0;
    }
	
	/* organic matter factor **************************/
	f_org[0] = 1.0;
	for(f=1;f<=N_SLAYER;f++){
		/* vegetated soil */
		if(dpth[f] < rdepth){
			f_org[f] = 1.0;
		}else{
			f_org[f] = exp(-fabs(dpth[f]*100.0 - rdepth*100.0)/10.0);
		}
	}
	
	/* soil-temperature coefficient profile *******************/
	tmp[0] = grid->tmp_sfc[grid->m];
	for(f=1;f<=N_SLAYER;f++){		
		if(tmp[f] > 0.0){
			f_t[f] = 1.0;
		}else{
			f_t[f] = 0.0;
		}
	}
	
	/* carbon input factor: Eq.2 */
	if(loct->npp_max > 0.0){
        if(FIX_CH4_NPP == 0){
            f_in = 1.0 + (flux->plant).npp[grid->m]/loct->npp_max;
        }else if(FIX_CH4_NPP == 1){
            f_in = 1.0 + loct->npp_av[grid->m]/loct->npp_max;
        }
	}else{
		f_in = 0.0;
	}
	
	/* sensitivity of Q10 of CH4 production: A.Ito (2010/08/02) */
	/* q10_ch4prod = 6.0; */
	/* q10_ch4prod = 4.0; */
	/* q10_ch4prod = 3.0; */
	/* q10_ch4prod = 3.2; */
	/* q10_ch4prod = 3.85; */
    /* q10_ch4prod = 2.5; */ /* 2018/05/25  */
    /* q10_ch4prod = 2.4; */ /* 2018/05/26  */
    q10_ch4prod = 2.0; /* 2018/05/28  */
	if(EX_CH4_2 == 1){
		q10_ch4prod = 1.5;
	}else if(EX_CH4_2 == 2){
		q10_ch4prod = 2.5;
	}else if(EX_CH4_2 == 3){
        /* original Walter-Heimann (2001) value */
        q10_ch4prod = 6.0;
    }else if(EX_CH4_2 == 4){
        /* 2014/12/10 
        Yvon-Durocher, G., A. P. Allen, D. Bastviken, R. Conrad, C. Gudasz, A. St-Pierre,
        N. Thanh-Duc, and P. A. del Giorgio (2014),
        Methane fluxes show consistent temperature dependence across microbial to
        ecosystem scales, Nature, 507, 488–491, doi:10.1038/nature13164.  */
        q10_ch4prod = 3.85;
    }
    
    /* parameter ensemble: 2014/11/19  */
    if(PARAM_PTB == 5){
        if(PARAM_ENS==1){
            r0 *= 0.7;
        }
        if(PARAM_ENS==2){
            r0 *= 0.8;
        }
        if(PARAM_ENS==3){
            r0 *= 0.9;
        }
        if(PARAM_ENS==4){
            r0 *= 1.1;
        }
        if(PARAM_ENS==5){
            r0 *= 1.2;
        }
        if(PARAM_ENS==6){
            r0 *= 1.3;
        }
    }
		
	release = 0.0;
    flux_ebull = flux_plant = 0.0;
	for(g=1;g<=cumtime;g++){ 
		
		/* DEPTH **************************************/
		flux_ebull = flux_plant = 0.0;
		day_produc = day_oxid = 0.0;
		ff[0] = 0.0;
		for(f=1;f<=N_SLAYER;f++){
			/* when water-table falls down */
			/* if(g==1 && loct->water_table_depth>loct->water_table_depth_pre){
				if(dpth[f]<=loct->water_table_depth && dpth[f]>=loct->water_table_depth_pre){
					release += df_dry * (loct->prof_ch4[f] - loct->prof_ch4[0]) * (1000.0*hh);
				}else{
					release += 0.0;
				}
			} */
			release = 0.0;
			
			/* ebullition at high CH4 concentration */
			if( loct->prof_ch4[f] > b_thresh){
				q_ebull[f] = -1.0 * (loct->prof_ch4[f] - b_thresh);
				flux_ebull += -1000.0 * hh * q_ebull[f];
			}else{
				q_ebull[f] = 0.0;
			}
			
			/* plant flux above rooting depth */
			if(dpth[f] <= rdepth){
				q_plant[f] = -0.01 * t_veg * 2.0*(rdepth - dpth[f])/rdepth * f_grow * loct->prof_ch4[f];
				flux_plant += -1000.0*hh * q_plant[f] * (1.0 - 0.5);
			}else{
				q_plant[f] = 0.0;
			}
			
			/* CH4 production and oxidation */
			if(dpth[f] >= wtdepth){
				/* below water table: eq.5 */
				q_prod[f] = r0 * f_org[f] * f_in * f_t[f] * pow(q10_ch4prod, (tmp[f] - t_mean)/10.0);
				day_produc += q_prod[f];
				q_oxid[f] = 0.0;
			}else{
				/* above water table */
				q_prod[f] = 0.0;
				q_oxid[f] = -1.0 * (20.0 * loct->prof_ch4[f])/(3.0 + loct->prof_ch4[f]) * pow(2.0, (tmp[f] - t_mean)/10.0);
				day_oxid += q_oxid[f];
			}
			
			/* total flux, micro M / hour */
			ff[f] = (q_oxid[f] + q_prod[f] + q_ebull[f] + q_plant[f]);
		}
		
		/* diffusion equation solved by explicit method */
		for(f=1;f<=(N_SLAYER-1);f++){
			/* ww[f] = df[f]*rr*(loct->prof_ch4[f+1] + loct->prof_ch4[f-1]) 
					+ (1.0 - 2.0*df[f]*rr)*loct->prof_ch4[f] + kk*ff[f]; */
			
			ww[f] = loct->prof_ch4[f] + rr * (df[f] * (loct->prof_ch4[f+1]
						- 2.0*loct->prof_ch4[f] + loct->prof_ch4[f-1]) + ff[f]);
			
			if(ww[f]<=0.0){
				ww[f] = 0.0;
			}
			if(ww[f] > b_thresh*1.2){
				ww[f] = b_thresh*1.2;
			}
		}
		/* boundary condition at the bottom: no downward flux */
		ww[N_SLAYER] = loct->prof_ch4[N_SLAYER] + rr*(df[N_SLAYER]*
								(2.0*loct->prof_ch4[N_SLAYER-1] - 2.0*loct->prof_ch4[N_SLAYER]));
		if(ww[N_SLAYER]<=0.0){
			ww[N_SLAYER] = 0.0;
		}
		
		/* ww[LAYER] = 500.0*16.0/(8.3144*(tmp[0]+273.15))/1000000.0*apress; */
		
		/* concentration update */
		for(f=1;f<=N_SLAYER;f++){
			loct->prof_ch4[f] = ww[f];
			
			/* loct->prof_ch4[f] = 100.0; */		/* debug ***************/
		}
	}
	
	/* loct->xx1[grid->m] = day_produc;
	loct->xx2[grid->m] = day_oxid;
	loct->xx3[grid->m] = f_grow;
	loct->xx4[grid->m] = loct->prof_ch4[5];
	loct->xx5[grid->m] = loct->prof_ch4[15]; */
	
	/* CH4 emission area fraction **************************************/
    f_inundation = 0.0;
    fa_wetland = fa_paddy = 0.0;

	/* sensitivity experiments: 2010/07/02  */
	if(EX_CH4_1 >= 1){
		/* West Siberia */
		switch(smode){
			case 1: case 2:
				/* f_inundation = grid->inundation_ssmi_av; */
				
				f_inundation = 0.333 * (grid->inundation_ssmi[4] + grid->inundation_ssmi[5] + 
								grid->inundation_ssmi[6] + grid->inundation_ssmi[7])/4.0;
				
				/* f_inundation = grid->inundation_ssmi_max * 0.5; */
				break;
			case 3: case 4:
				f_inundation = grid->inundation_ssmi[grid->m];
				break;
		}
	}else{
		/* global */
        if(ALT_INUND == 0){
            f_inundation = grid->inundation_ssmi[grid->m]; /* */
        }else if(ALT_INUND == 1 || ALT_INUND == 3 || ALT_INUND == 4 || ALT_INUND == 5){
            if(grid->climy >= 1999 && grid->climy <= 2013){
        
                if(grid->f_wetland > 0.0){
                    f_inundation = grid->inundation_alt_ts[grid->climy-1999][grid->m] / grid->f_wetland;
                    /* f_inundation = grid->inundation_alt_ts[grid->climy-1999][grid->m] */;
                    
                    if(f_inundation > 1.0){
                        f_inundation = 1.0;
                    }
                }
                if(grid->f_paddy > 0.0){
                    /* f_inundation = grid->inundation_alt_ts[grid->climy-1999][grid->m] / grid->f_paddy; */
                    
                    f_inundation = grid->inundation_ssmi[grid->m];
                    
                    if(f_inundation > 1.0){
                        f_inundation = 1.0;
                    }
                }
            
            }else{
                if(grid->f_wetland > 0.0){
                    f_inundation = grid->inundation_alt_av[grid->m] / grid->f_wetland;
                    /* f_inundation = grid->inundation_alt_av[grid->m] */;
                    
                    if(f_inundation > 1.0){
                        f_inundation = 1.0;
                    }
                }
                if(grid->f_paddy > 0.0){
                    /* f_inundation = grid->inundation_alt_av[grid->m] / grid->f_paddy; */
                    
                    f_inundation = grid->inundation_ssmi[grid->m];
                    
                    if(f_inundation > 1.0){
                        f_inundation = 1.0;
                    }
                }
            }
        }else if(ALT_INUND == 6){
            /* use GCP-CH4 data */
            if(grid->climy >= 2000 && grid->climy <= 2012){
        
                if(grid->f_wetland > 0.0){
                    f_inundation = grid->inundation_alt_ts[grid->climy-2000][grid->m] / grid->f_wetland;
                    /* f_inundation = grid->inundation_alt_ts[grid->climy-2000][grid->m] */;
                    
                    if(f_inundation > 1.0){
                        f_inundation = 1.0;
                    }
                }
                if(grid->f_paddy > 0.0){
                    /* f_inundation = grid->inundation_alt_ts[grid->climy-2000][grid->m] / grid->f_paddy; */
                    
                    f_inundation = grid->inundation_ssmi[grid->m];
                    
                    if(f_inundation > 1.0){
                        f_inundation = 1.0;
                    }
                }
            
            }else{
                if(grid->f_wetland > 0.0){
                    f_inundation = grid->inundation_alt_av[grid->m] / grid->f_wetland;
                    /* f_inundation = grid->inundation_alt_av[grid->m] */;
                    
                    if(f_inundation > 1.0){
                        f_inundation = 1.0;
                    }
                }
                if(grid->f_paddy > 0.0){
                    /* f_inundation = grid->inundation_alt_av[grid->m] / grid->f_paddy; */
                    
                    f_inundation = grid->inundation_ssmi[grid->m];
                    
                    if(f_inundation > 1.0){
                        f_inundation = 1.0;
                    }
                }
            }
        }else if(ALT_INUND == 7 || ALT_INUND == 8){
            /* use GCP-CH4 v2 data: 2018/08/29  */
            if(grid->climy >= 2000 && grid->climy <= 2017){
        
                if(grid->f_wetland > 0.0){
                    f_inundation = grid->inundation_alt_ts[grid->climy-2000][grid->m] / grid->f_wetland;
                    /* f_inundation = grid->inundation_alt_ts[grid->climy-2000][grid->m] */;
                    
                    if(f_inundation > 1.0){
                        f_inundation = 1.0;
                    }
                }
                if(grid->f_paddy > 0.0){
                    /* f_inundation = grid->inundation_alt_ts[grid->climy-2000][grid->m] / grid->f_paddy; */
                    f_inundation = grid->inundation_ssmi[grid->m];
                    
                    if(f_inundation > 1.0){
                        f_inundation = 1.0;
                    }
                }
            
            }else{
                if(grid->f_wetland > 0.0){
                    f_inundation = grid->inundation_alt_av[grid->m] / grid->f_wetland;
                    /* f_inundation = grid->inundation_alt_av[grid->m] */;
                    
                    if(f_inundation > 1.0){
                        f_inundation = 1.0;
                    }
                }
                if(grid->f_paddy > 0.0){
                    /* f_inundation = grid->inundation_alt_av[grid->m] / grid->f_paddy; */
                    f_inundation = grid->inundation_ssmi[grid->m];
                    
                    if(f_inundation > 1.0){
                        f_inundation = 1.0;
                    }
                }
            }
        }else if(ALT_INUND == 9){
            if(grid->climy >= 1992 && grid->climy <= 2020){
                /* SWAMPS anomaly: 2021/06/25  */
                f_inundation = grid->inundation_ssmi[grid->m] + grid->inundation_alt_ts[grid->climy-1992][grid->m];

                if(f_inundation < 0.0){
                    f_inundation = 0.0;
                }
            }else{
                f_inundation = grid->inundation_ssmi[grid->m];
            }
        }else if(ALT_INUND == 10){
            /* use GCP-CH4 WAD2M data: 2021/10/26  */
            if(grid->climy >= 2000 && grid->climy <= 2020){
        
                if(grid->f_wetland > 0.0){
                    f_inundation = grid->inundation_alt_ts[grid->climy-2000][grid->m] / grid->f_wetland;
                    /* f_inundation = grid->inundation_alt_ts[grid->climy-2000][grid->m] */;
                    
                    if(f_inundation > 1.0){
                        f_inundation = 1.0;
                    }
                }
                if(grid->f_paddy > 0.0){
                    /* f_inundation = grid->inundation_alt_ts[grid->climy-2000][grid->m] / grid->f_paddy; */
                    f_inundation = grid->inundation_ssmi[grid->m];
                    
                    if(f_inundation > 1.0){
                        f_inundation = 1.0;
                    }
                }
            
            }else{
                if(grid->f_wetland > 0.0){
                    f_inundation = grid->inundation_alt_av[grid->m] / grid->f_wetland;
                    /* f_inundation = grid->inundation_alt_av[grid->m] */;
                    
                    if(f_inundation > 1.0){
                        f_inundation = 1.0;
                    }
                }
                if(grid->f_paddy > 0.0){
                    /* f_inundation = grid->inundation_alt_av[grid->m] / grid->f_paddy; */
                    f_inundation = grid->inundation_ssmi[grid->m];
                    
                    if(f_inundation > 1.0){
                        f_inundation = 1.0;
                    }
                }
            }
        }
		
		/* when using NASA/GISS wetland data: 2011/03/31  */
		if(ALT_FWETLAND == 1 && (smode == 1 || smode == 2)){
			/* to avoid double-counting of inundation fraction */
			f_inundation = 1.0;
		}
		
		/* for West Siberia */
		/* f_inundation = 0.333 * (grid->inundation_ssmi[4] + grid->inundation_ssmi[5] + 
								grid->inundation_ssmi[6] + grid->inundation_ssmi[7])/4.0; */
		/* f_inundation = grid->inundation_ssmi_av; */		
	}
	
    /* ***************************/
	if(EX_CH4_1 == 1){
		/* experiment for Sasakawa-san: low sensitivity */
		switch(smode){
			case 1:
				if(f_inundation < 0.05){
					f_inundation = 0.05;
				}else{
					f_inundation *= 1.0 + loct->cum_dprec * 0.001;
				}
				
				fa_wetland = f_inundation * grid->f_wetland;
				
				if(fa_wetland > 0.99){
					fa_wetland = 0.99;
				}
				break;
			case 2:
				if(f_inundation < 0.05){
					f_inundation = 0.05;
				}else{
					f_inundation *= 1.0 + loct->cum_dprec * 0.001;
				}

				fa_wetland = (1.0 - f_inundation) * grid->f_wetland;
				
				if(fa_wetland < 0.01){
					fa_wetland = 0.01;
				}
				break;
			case 3:
				fa_paddy = f_inundation * grid->f_paddy;
				break;
			case 4:
				fa_paddy = (1.0 - f_inundation) * grid->f_paddy;
				break;
		}
	}else if(EX_CH4_1 == 2){
		/* experiment for Sasakawa-san: high sensitivity */
		switch(smode){
			case 1:
				if(f_inundation < 0.05){
					f_inundation = 0.05;
				}else{
					f_inundation *= 1.0 + loct->cum_dprec * 0.003;
				}
				
				fa_wetland = f_inundation * grid->f_wetland;
				
				if(fa_wetland > 0.99){
					fa_wetland = 0.99;
				}
				break;
			case 2:
				if(f_inundation < 0.05){
					f_inundation = 0.05;
				}else{
					f_inundation *= 1.0 + loct->cum_dprec * 0.003;
				}
				
				fa_wetland = (1.0 - f_inundation) * grid->f_wetland;
				
				if(fa_wetland < 0.01){
					fa_wetland = 0.01;
				}
				break;
			case 3:
				fa_paddy = f_inundation * grid->f_paddy;
				break;
			case 4:
				fa_paddy = (1.0 - f_inundation) * grid->f_paddy;
				break;
		}
	}else{
		/* control */
		switch(smode){
			case 1:
				fa_wetland = f_inundation * grid->f_wetland;
				break;
			case 2:
				fa_wetland = (1.0 - f_inundation) * grid->f_wetland;
				break;
			case 3:
				fa_paddy = f_inundation * grid->f_paddy;
				break;
			case 4:
				fa_paddy = (1.0 - f_inundation) * grid->f_paddy;
				break;
		}
        
        /* no limit by GLWD: 2018/08/29  */
        if(ALT_INUND == 8 && smode == 1){
        
            if(grid->climy >= 2000 && grid->climy <= 2017){
                fa_wetland = grid->inundation_alt_ts[grid->climy-2000][grid->m];
            }else{
                fa_wetland = grid->inundation_alt_av[grid->m];
            }
            
        }else if(ALT_INUND == 8 && smode == 2){
        
            if(grid->climy >= 2000 && grid->climy <= 2017){
                if(grid->inundation_alt_ts[grid->climy-2000][grid->m] > grid->f_wetland){
                    fa_wetland = 0.0;
                }else{
                    if(f_inundation <= 1.0){
                        fa_wetland = (1.0 - f_inundation) * grid->f_wetland;
                    }else{
                        fa_wetland = 0.0;
                    }
                }
            }else{
                fa_wetland = grid->inundation_alt_av[grid->m];
            }
        }
        /* no limit by WAD2M: 2021/10/26  */
        if(ALT_INUND == 10 && smode == 1){
        
            if(grid->climy >= 2000 && grid->climy <= 2020){
                fa_wetland = grid->inundation_alt_ts[grid->climy-2000][grid->m];
            }else{
                fa_wetland = grid->inundation_alt_av[grid->m];
            }
            
        }else if(ALT_INUND == 10 && smode == 2){
        
            if(grid->climy >= 2000 && grid->climy <= 2020){
                if(grid->inundation_alt_ts[grid->climy-2000][grid->m] > grid->f_wetland){
                    fa_wetland = 0.0;
                }else{
                    if(f_inundation <= 1.0){
                        fa_wetland = (1.0 - f_inundation) * grid->f_wetland;
                    }else{
                        fa_wetland = 0.0;
                    }
                }
            }else{
                fa_wetland = grid->inundation_alt_av[grid->m];
            }
        }
        
        /* SWAMPS anomaly: 2021/06/25  */
        if(ALT_INUND == 9 && smode == 1){
            fa_wetland = f_inundation * grid->f_wetland;
            
            if(fa_wetland > 1.0){
                fa_wetland = 1.0;
            }
            if(fa_wetland < 0.0){
                fa_wetland = 0.0;
            }
        }else if(ALT_INUND == 9 && smode == 2){
            fa_wetland = (1.0 - f_inundation) * grid->f_wetland;
            
            if(fa_wetland > 1.0){
                fa_wetland = 1.0;
            }
            if(fa_wetland < 0.0){
                fa_wetland = 0.0;
            }
        }
	}
	
	/* sensitivity analysis for lake area: 2010/08/18 (A.Ito)  */
	if(EX_CH4_3 == 1){
		if(smode==1){
			fa_wetland += 1.0 * grid->f_lake;
		}
	}else{
		if(smode==1){
			fa_wetland += 0.2 * grid->f_lake;
		}
	}
 
    /* GCP-CH4: per unit wetland area */
    /* 2022/04/14  */
    if(PUA_WETLAND == 1){
		switch(smode){
			case 1:
				fa_wetland = f_inundation * 1.0;
				break;
			case 2:
				fa_wetland = (1.0 - f_inundation) * 1.0;
				break;
			case 3:
				fa_paddy = f_inundation * 1.0;
				break;
			case 4:
				fa_paddy = (1.0 - f_inundation) * 1.0;
				break;
		}
    }
	
    /* alternative paddy crop calendar: 2021/04/07  */
    if(EX_PADDY == 3 && (smode == 3 || smode == 4)){
        if(grid->iizumi_mon_paddy_end > grid->iizumi_mon_paddy_start){
            if(grid->m >= grid->iizumi_mon_paddy_start && grid->m <= grid->iizumi_mon_paddy_end){
                fa_paddy =  grid->f_paddy;
            }else{
                fa_paddy = 0.0;
            }
        }else if(grid->iizumi_mon_paddy_end < grid->iizumi_mon_paddy_start){
            if(grid->m <= grid->iizumi_mon_paddy_start || grid->m >= grid->iizumi_mon_paddy_end){
                fa_paddy = grid->f_paddy;
            }else{
                fa_paddy = 0.0;
            }
        }else{
            fa_paddy = 0.0;
        }
    }

	/* flux: mg CH4 m-2 month-1 ************************/
	if(smode == 1){
		loct->f_inund_wet_wh[grid->m] = fa_wetland;
        /* loct->xx5[grid->m] = fa_wetland; */
        
        /* if((grid->f_wetland + grid->f_paddy) > 0.0){
            loct->xx7[grid->m] += f_inundation * grid->f_wetland / (grid->f_wetland + grid->f_paddy) * 100.0*wtdepth;
        } */
        
        efflux_plant = fa_wetland * flux_plant / 1000.0 * DHN * 16.0 * MDN[grid->m];
        efflux_ebul = fa_wetland * flux_ebull * DHN * 16.0 / 1000.0 * MDN[grid->m];
        efflux_diffs = fa_wetland * df[1]/(dpth[1] - dpth[0])*
					(loct->prof_ch4[1] - loct->prof_ch4[0]) * DHN * 16.0 / 1000.0 * MDN[grid->m];
        efflux_reles = fa_wetland * release * DHN * 16.0 / 1000.0 * MDN[grid->m];
		
		(flux->soil).ch4_wetland_wh_plant[grid->m] += efflux_plant;
		(flux->soil).ch4_wetland_wh_ebull[grid->m] += efflux_ebul;
		(flux->soil).ch4_wetland_wh_diff[grid->m] += efflux_diffs;
		(flux->soil).ch4_wetland_wh_release[grid->m] += efflux_reles;
        
        /* loct->xx1[grid->m] = efflux_ebul + efflux_plant + efflux_diffs + efflux_reles; */
	}if(smode == 2){
        /* loct->xx6[grid->m] = fa_wetland; */
        
        /* if((grid->f_wetland + grid->f_paddy) > 0.0){
            loct->xx7[grid->m] += (1.0 - f_inundation) * grid->f_wetland / (grid->f_wetland + grid->f_paddy) * 100.0*wtdepth;
        } */
        
        efflux_plant = fa_wetland * flux_plant / 1000.0 * DHN * 16.0 * MDN[grid->m];
        efflux_ebul = fa_wetland * flux_ebull * DHN * 16.0 / 1000.0 * MDN[grid->m];
        efflux_diffs = fa_wetland * df[1]/(dpth[1] - dpth[0])*
					(loct->prof_ch4[1]-loct->prof_ch4[0]) * DHN * 16.0 / 1000.0 * MDN[grid->m];
        efflux_reles = fa_wetland * release * DHN * 16.0 / 1000.0 * MDN[grid->m];
        
		(flux->soil).ch4_wetland_wh_plant[grid->m] += efflux_plant;
		(flux->soil).ch4_wetland_wh_ebull[grid->m] += efflux_ebul;
		(flux->soil).ch4_wetland_wh_diff[grid->m] += efflux_diffs;
		(flux->soil).ch4_wetland_wh_release[grid->m] += efflux_reles;
        
        /* loct->xx2[grid->m] = efflux_ebul + efflux_plant + efflux_diffs + efflux_reles; */
	}else if(smode == 3){
		loct->f_inund_pad_wh[grid->m] = fa_paddy;
        /* loct->xx7[grid->m] = fa_paddy; */
        
        /* loct->xx6[grid->m] += fa_paddy;
        if((grid->f_wetland+grid->f_paddy) > 0.0){
            loct->xx7[grid->m] += f_inundation * grid->f_paddy / (grid->f_wetland + grid->f_paddy) * 100.0*wtdepth;
        } */
        
        efflux_plant = fa_paddy * flux_plant / 1000.0 * DHN * 16.0 * MDN[grid->m];
        efflux_ebul = fa_paddy * flux_ebull * DHN * 16.0 / 1000.0 * MDN[grid->m];
        efflux_diffs = fa_paddy * df[1]/(dpth[1] - dpth[0])*
					(loct->prof_ch4[1]-loct->prof_ch4[0]) * DHN * 16.0 / 1000.0 * MDN[grid->m];
        efflux_reles = fa_paddy * release * DHN * 16.0 / 1000.0 * MDN[grid->m];
        
		(flux->soil).ch4_paddy_wh_plant[grid->m] += efflux_plant;
		(flux->soil).ch4_paddy_wh_ebull[grid->m] += efflux_ebul;
		(flux->soil).ch4_paddy_wh_diff[grid->m] += efflux_diffs;
		(flux->soil).ch4_paddy_wh_release[grid->m] += efflux_reles;
        
        /* loct->xx3[grid->m] = efflux_ebul + efflux_plant + efflux_diffs + efflux_reles; */
	}else if(smode == 4){
        /* loct->xx8[grid->m] = fa_paddy; */
        
        /* if((grid->f_wetland + grid->f_paddy) > 0.0){
            loct->xx7[grid->m] += (1.0 - f_inundation) * grid->f_paddy / (grid->f_wetland + grid->f_paddy)  * 100.0*wtdepth;
        } */
        
        efflux_plant = fa_paddy * flux_plant / 1000.0 * DHN * 16.0 * MDN[grid->m];
        efflux_ebul = fa_paddy * flux_ebull * DHN * 16.0 / 1000.0 * MDN[grid->m];
        efflux_diffs = fa_paddy * df[1]/(dpth[1] - dpth[0])*
					(loct->prof_ch4[1]-loct->prof_ch4[0]) * DHN * 16.0 / 1000.0 * MDN[grid->m];
        efflux_reles = fa_paddy * release * DHN * 16.0 / 1000.0 * MDN[grid->m];
        
		(flux->soil).ch4_paddy_wh_plant[grid->m] += efflux_plant;
		(flux->soil).ch4_paddy_wh_ebull[grid->m] += efflux_ebul;
		(flux->soil).ch4_paddy_wh_diff[grid->m] += efflux_diffs;
		(flux->soil).ch4_paddy_wh_release[grid->m] += efflux_reles;
        
        /* loct->xx4[grid->m] = efflux_ebul + efflux_plant + efflux_diffs + efflux_reles; */
	}
    
    /* if((grid->f_wetland + grid->f_paddy) <= 0.0){
        loct->xx7[grid->m] = -9999.0;
    } */

	loct->water_table_depth_pre = loct->water_table_depth;
}

/* aerobic CH4 emission ***********************************************************/
/*
Keppler, F., Hamilton, J.T.G., Bra, M. and Rkmann, T., 2006. Methane emissions from 
terrestrial plants under aerobic conditions. Nature 439, 187-191.
*/
/*
Kirschbaum, M.U.F., Bruhn, D., Etheridge, D.M., Evans, J.R., Farquhar, G.D., Gifford, R.M., 
Paul, K.I. and Winters, A.J., 2006. A comment on the quantitative significance of aerobic 
methane releasse by plants. Functional Plant Biology 33, 521-530.
*/
void f_ch4_emit_veg(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux
){
	double femit_sun, femit_shade;
	double sunshine;
	extern double MDN[ASTEP];
	
	/* base emission rate */
	femit_sun = 374.0;		/* ng gdw-1 h-1 */
	femit_shade = 119.0;	/* ng gdw-1 h-1 */
	
	sunshine = grid->dlen[grid->m] * (1.0 - grid->tcdc_clm[grid->m]);
	
	/* MASS-based scaling up **/
	/* C3, g m-2 month-1 */
	if((echar->c3).season[grid->m] != 0){
		(flux->c3).emit_ch4_kirschbaum_mass[grid->m] = ( (mass->c3).mfol[grid->m] * dmTc * 100.0 ) *
			(sunshine*femit_sun + (DHN - sunshine) * femit_shade) * pow(10.0, -9.0) * MDN[grid->m];
	}else{
		(flux->c3).emit_ch4_kirschbaum_mass[grid->m] = 0.0;
	}

	/* C4, g m-2 month-1 */
	if((echar->c4).season[grid->m] != 0){
		(flux->c4).emit_ch4_kirschbaum_mass[grid->m] = ( (mass->c4).mfol[grid->m] * dmTc * 100.0 ) * 
			(sunshine*femit_sun + (DHN - sunshine)*femit_shade) * pow(10.0, -9.0) * MDN[grid->m];
	}else{
		(flux->c4).emit_ch4_kirschbaum_mass[grid->m] = 0.0;
	}
	
	/* PHOTO-based scaling up **/
	if(sunshine > 0.0){
		(flux->c3).emit_ch4_kirschbaum_photo[grid->m] = 2.0 * (16.0/12.0) * 
			((flux->c3).npp[grid->m]*1000.0) / 30000.0 * (1.0 + (DHN - sunshine)/sunshine * 
			femit_shade / femit_sun);
		(flux->c4).emit_ch4_kirschbaum_photo[grid->m] = 2.0 * (16.0/12.0) * 
			((flux->c4).npp[grid->m]*1000.0) / 30000.0 * (1.0 + (DHN - sunshine)/sunshine *
			femit_shade / femit_sun);
	}else{
		(flux->c3).emit_ch4_kirschbaum_photo[grid->m] = 0.0;
		(flux->c4).emit_ch4_kirschbaum_photo[grid->m] = 0.0;
	}
}

void f_ch4_emit_termite(
    struct Grid *grid,
    struct Loct *loct,
	struct Echar *echar,
	struct Mass *mass,
    struct Flux *flux
){
    short f;
    double pch4, trmden, gppann;
    
    trmden = pch4 = 0.0;
    
    if(loct->v_type == 1){
        switch(grid->veg_sage){
            case 1: /* Tropical Evergreen Forest/Woodland */
                //pch4 = (5.6 * (5.4 * 0.000001 * 24.0 * MDN[grid->m]) + 11.0 * (5.4 * 0.000001 * 24.0 * MDN[grid->m]))/2.0;
                trmden = 5.6;
                pch4 = trmden * (5.4 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 2: /* Tropical Deciduous Forest/Woodland */
                //pch4 = (5.6 * (5.4 * 0.000001 * 24.0 * MDN[grid->m]) + 11.0 * (5.4 * 0.000001 * 24.0 * MDN[grid->m]))/2.0;
                trmden = 5.6;
                pch4 = trmden * (5.4 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 3: /* Temperate Broadleaf Evergreen Forest/Woodland */
                trmden = 3.0;
                pch4 = trmden * (1.8 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 4: /* Temperate Needleleaf Evergreen Forest/Woodland */
                trmden = 3.0;
                pch4 = trmden * (1.8 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 5: /* Temperate Deciduous Forest/Woodland */
                trmden = 3.0;
                pch4 = trmden * (1.8 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 6: /* Boreal Evergreen Forest/Woodland */
                trmden = 0.0;
                pch4 = trmden * (5.4 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 7: /* Boreal Deciduous Forest/Woodland */
                trmden = 0.0;
                pch4 = trmden * (5.4 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 8: /* Evergreen/Deciduous Mixed Forest/Woodland */
                trmden = 3.0;
                pch4 = trmden * (1.8 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 9: /* Savanna */
                //pch4 = (4.5 * (8.0 * 0.000001 * 24.0 * MDN[grid->m]) + 11.1 * (8.0 * 0.000001 * 24.0 * MDN[grid->m]))/2.0;
                trmden = 4.5;
                pch4 = trmden * (8.0 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 10: /* Grassland/Steppe */
                //pch4 = (3.0 * (1.8 * 0.000001 * 24.0 * MDN[grid->m]) + 5.2 * (1.8 * 0.000001 * 24.0 * MDN[grid->m]))/2.0;
                trmden = 3.0;
                pch4 = trmden * (1.8 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 11: /* Dense Shrubland */
                //pch4 = (3.1 * (1.0 * 0.000001 * 24.0 * MDN[grid->m]) + 3.1 * (1.0 * 0.000001 * 24.0 * MDN[grid->m]))/2.0;
                trmden = 3.1;
                pch4 = trmden * (1.0 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 12: /* Open Shrubland */
                //pch4 = (3.1 * (1.0 * 0.000001 * 24.0 * MDN[grid->m]) + 0.8 * (1.0 * 0.000001 * 24.0 * MDN[grid->m]))/2.0;
                trmden = 3.1;
                pch4 = trmden * (1.0 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 13: /* Tundra */
                trmden = 0.0;
                pch4 = trmden * (5.4 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 14: /* Desert */
                //pch4 = (3.1 * (1.0 * 0.000001 * 24.0 * MDN[grid->m]) + 1.0 * (1.0 * 0.000001 * 24.0 * MDN[grid->m]))/2.0;
                trmden = 3.1;
                pch4 = trmden * (1.0 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 15: /* Polar Desert/Rock/Ice */
                trmden = 0.0;
                pch4 = trmden * (5.4 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            default:
                trmden = 0.0;
                pch4 = 0.0;
        }
    }
    
    /* Land-use based estimation */
    if(loct->v_type == 2){
        trmden = 7.8;
        pch4 = trmden * (1.8 * 0.000001 * 24.0 * MDN[grid->m]);
    }
    
    /* Temperature-limitation */
    /* if(grid->tmp_sfc_am < 0.0){ */
    if(grid->tmp_sfc_mn < -8.0){
        trmden = pch4 = 0.0;
    }
    
    flux->termite_ch4_lu[grid->m] = pch4;
    flux->termite_dens_lu[grid->m] = trmden;
    
    /* GPP-based estimation *********************************************/
    pch4 = gppann = 0.0;
    for(f=0;f<12;f++){
        gppann += (flux->plant).gpp[f];
    }
    
    if(gppann >= 0.0 && gppann <= 100.0){
        //trmden = 1.21 * exp(0.08 * gppann);
        trmden = 1.21 * exp(0.08 * flux->gpp_ann);
    }else{
        trmden = 0.0;
    }
    
    if(loct->v_type == 1){
        switch(grid->veg_sage){
            case 1: /* Tropical Evergreen Forest/Woodland */
            case 2: /* Tropical Deciduous Forest/Woodland */
                pch4 = trmden * (2.8 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 3: /* Temperate Broadleaf Evergreen Forest/Woodland */
            case 4: /* Temperate Needleleaf Evergreen Forest/Woodland */
            case 5: /* Temperate Deciduous Forest/Woodland */
            case 6: /* Boreal Evergreen Forest/Woodland */
            case 7: /* Boreal Deciduous Forest/Woodland */
            case 8: /* Evergreen/Deciduous Mixed Forest/Woodland */
                trmden = 3.0;
                pch4 = trmden * (1.7 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 9: /* Savanna */
            case 10: /* Grassland/Steppe */
                trmden = (5.2 + 10.6)/2.0;
                pch4 = trmden * (1.7 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 11: /* Dense Shrubland */
                trmden = (5.3 + 8.43)/2.0;
                pch4 = trmden * (1.7 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 12: /* Open Shrubland */
                trmden = (5.3 + 8.43)/2.0;
                pch4 = trmden * (1.7 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 13: /* Tundra */
                trmden = 0.0;
                pch4 = trmden * (1.7 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 14: /* Desert */
                trmden = (0.96 + 0.98)/2.0;
                pch4 = trmden * (1.7 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            case 15: /* Polar Desert/Rock/Ice */
                trmden = 0.0;
                pch4 = trmden * (1.7 * 0.000001 * 24.0 * MDN[grid->m]);
                break;
            default:
                trmden = 0.0;
                pch4 = 0.0;
        }
    }
    
    if(loct->v_type == 2){
        trmden = (5.38 + 2.25)/2.0;
        pch4 = trmden * (1.7 * 0.000001 * 24.0 * MDN[grid->m]);
    }
    
    /* if(grid->lat > 45.0  || grid->lat < -45.0){
        pch4 = 0.0;
    } */
    /* if(grid->tmp_sfc_am < 0.0){ */
    if(grid->tmp_sfc_mn < -8.0){
        pch4 = trmden = 0.0;
    }
    
     flux->termite_ch4_gpp[grid->m] = pch4;
     flux->termite_dens_gpp[grid->m] = trmden;
}
