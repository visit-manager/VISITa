/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

#define TER_CON 0.001 /* criteria for determining the equilibrium, NEP value in Mg C ha-1 yr-1 */

/* EQUILIBRIUM ***************************************************************/
void cal_spinup(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux, 
	FILE *fp_o[OFILEN]
){
	long f, g, nn, term_time, dyr, y_nin;
	double plantmass, ann_nep, f_fert, total_hvst, aa,bb;
    double f_nat, iweight, iweight3, avc3, prm_ensen, icrop, base_nin;
	
	/** maximum simulation times **/
	grid->phase = 0; /* spin-up */
	
	/* max. spin-up time, years *************/
	if((echar->c3).v_type == 1){
		/* corrected: A.Ito and E.Kato (2009/08/16) */
		if(grid->veg_olson==29 || grid->veg_olson==30 || grid->veg_olson==31 
				|| grid->veg_olson==32){
            /* short-term ecosystems */
			term_time = 300;
		}else{
            /* long-term ecosystems */
			term_time = 4000;
		}
	}else if((echar->c3).v_type==2){
		term_time = 300;
	}else{
        term_time = 300;
    }
	
	(echar->soil).rl = (echar->soil).rl0;
	(echar->soil).rh = (echar->soil).rh0;
	
	/* land-use change ***********/
	f_cult_luc(grid);
	if(LANDUSE == 0){
		grid->f_crop_p = 0.0;
		grid->f_pasture_p = 0.0;
	}else if(LANDUSE>=1 && LANDUSE<=5){
		grid->f_crop_p = grid->fcrop_net[199];
		grid->f_pasture_p = 0.0;
	}else if(LANDUSE == 6 || LANDUSE == 8){
		grid->f_crop_p = grid->fcrop_luh_hmnzed[199];
		grid->f_pasture_p = grid->fpast_luh_hmnzed[199];
	}else if(LANDUSE == 7){
		grid->f_crop_p = grid->fcrop_rk[199];
		grid->f_pasture_p = grid->fpast_rk[199];
	}else if(LANDUSE == 9){
		grid->f_crop_p = grid->fcrop_luh_hmnzed[2000-BGY_LUC];
		grid->f_pasture_p = grid->fpast_luh_hmnzed[2000-BGY_LUC];
	}else if(LANDUSE == 10 || LANDUSE == 11 || LANDUSE == 12 || LANDUSE == 13
         || LANDUSE == 14 || LANDUSE == 15 || LANDUSE == 16){
		grid->f_crop_p = grid->fcrop_luh_hmnzed[BGY_LUC - FDY_LUC];
		grid->f_pasture_p = grid->fpast_luh_hmnzed[BGY_LUC - FDY_LUC];
	}else if(LANDUSE == 18){
		grid->f_crop_p = grid->fcrop3_future[0];
		grid->f_pasture_p = 0.0;
    }else if(LANDUSE == 17 || BIOFUEL_RUN >= 1){
		grid->f_crop_p = grid->f_biofuel[0];
		grid->f_pasture_p = grid->fpast_luh_hmnzed[2000-BGY_LUC];
	}else if(LANDUSE == 19 || LANDUSE == 20 ||
            LANDUSE == 21 || LANDUSE == 22 || LANDUSE == 23){
		grid->f_crop_p = grid->fcrop_luh_hmnzed[1989 - FDY_LUC];
		grid->f_pasture_p = 0.0;
	}
	
    /* historical fertilizer */
    grid->niny = 1901;
    if(NMIP_RUN >= 1 || ISIMIP_RUN == 4){
        /* ISI-MIP2b: 2016/12/24 by A.Ito */
        grid->niny = FDY_NINY;

		grid->f_crop_p = grid->mip_frcrop[0];
		grid->f_pasture_p = 0.0;
    }
    if((echar->soil).v_type == 2){
        n_fertilizer_in(grid, loct);
    }else{
        loct->n_frtlz_in = 0.0;
        loct->n_manure_in = 0.0;
    }
    
    if(NMIP_RUN >= 1 || ISIMIP_RUN == 4){
        /* NMIP input: 2015/11/19 by A.Ito */
        /* ISI-MIP2b: 2016/12/24 by A.Ito */
        f_fert = 1.0; /* driven by data */
    }else{
        if(grid->rank_nat == 1){
            /* developing countries */
            f_fert = 2.0217112 / (1.0 + exp(0.049849599 * (2000.6575 - 1900.0)))+0.0014929171;
        }else if(grid->rank_nat == 2){
            /* developed countries */
            f_fert = 0.92939393 / (1.0 + exp(0.044112692 * (2000.0097 - 1900.0)))+0.53533202;
        }else{
            f_fert = 1.0;
        }
    }
    
    grid->simy = 1900;
    if(ISIMIP_RUN == 1 || (EX_BECCS==1 || EX_BECCS==2 || EX_BECCS==3)){
        grid->simy = 1949;
    }
    if(NMIP_RUN >= 1){
        grid->simy = FSY_HIST -1; /* 1860 */
    }
    if(ISIMIP_RUN == 4){
        grid->simy = FSY_HIST - 1; /* 1660 */
    }
	
	/* LOOP to stable stage ************************************************/
	nn = 0; 
	ann_nep = 10.0;
    /* corrected: A. Ito (with Hamada-san's comment) 2012/01/30 */
	loct->gpp_max = loct->npp_max = 0.0; 
	while(ann_nep > TER_CON){ /*** acnep>TER_CON nn<10 ***/
		
        /* AGE: years */
        grid->y = nn;
				
        /* for considering leap years: 2014/09/29 by A.Ito */
        if(grid->y%4 == 0){
            MDN[1] = 29.0;
            YDN = 366.0;
        }else{
            MDN[1] = 28.0;
            YDN = 365.0;
        }
			
        if(ISIMIP_RUN == 1 && grid->flag_histdata == 1){
            ann_nep = 10.0;
            grid->climy = grid->lucy = nn%30 +1950;
			set_hist_clim(grid);
		}else if(ISIMIP_RUN == 2 && grid->flag_histdata == 1){
            /* PLUME: 2014/07/31 by A.Ito */
            ann_nep = 10.0;
            grid->climy = grid->lucy = nn%30 +1901;
			set_hist_clim(grid);
        }else if(ISIMIP_RUN == 3 && grid->flag_histdata == 1){
            ann_nep = 10.0;
            grid->climy = grid->lucy = nn%30 +1901;
			set_hist_clim(grid);
        }else if(ISIMIP_RUN == 4 && grid->flag_histdata == 1){
            ann_nep = 10.0;
            grid->climy = grid->lucy = nn%30 + FSY_HIST;
			set_hist_clim(grid);
        }
        
        /* NMIP: 2015/11/19 by A.Ito **/
        /* updated: 2016/10/20 */
        if(ISIMIP_RUN == 1){
            grid->niny = FSY_HIST-1;
        }else if(ISIMIP_RUN == 4){
            grid->niny = FSY_HIST;
        }else{
            grid->niny = 1901;
        }
        
        if(NMIP_RUN >= 1){
            grid->climy = 1901;
            grid->niny = FSY_HIST-1;
            grid->co2y = FSY_HIST-1;
            grid->lucy = FSY_HIST-1;
            set_hist_clim(grid);
            n_fertilizer_in(grid, loct);
        }
		
        if((echar->soil).v_type == 2){
            if(EX_NFERT >= 1){
                /* with new niny */
                n_fertilizer_in(grid, loct);
                f_fert = 1.0; /* driven by data */
            }else{
                ;
            }
        }else{
            loct->n_frtlz_in = 0.0;
            loct->n_manure_in = 0.0;
        }

		plantmass = ann_nep = 0.0;
		for(f=0;f<ASTEP;f++){
			grid->m = f;
            
			/* initialize GHG fluxes ************/
			ghg_flux_zero(f, flux);
			
			/* atmospheric CO2 ****************/
			f_co2_trend(grid);
			co2_in_canopy(grid, loct, mass, flux);
			
			/* environmental condition ***/
			f_dyn_loct(grid, loct, mass, echar);
			
			/* vegetation processes *****/
			f_biome_processes(grid, loct, echar, mass, flux);

			/* VOC emission *****************/
			f_voc_emit_guenther97(grid, loct, echar, mass, flux);
			/* Plant CH4 emission *****************/
			f_ch4_emit_veg(grid, loct, echar, mass, flux);

			f_plant_stand_budget(grid, loct, mass, flux);
			
			/* GPP max for Cao's CH4 emission scheme */
			if((flux->plant).gpp[f] > loct->gpp_max){
				loct->gpp_max = (flux->plant).gpp[f];
			}
			
			/* maximum monthly npp, used in Walter & Heimann scheme */
			if((flux->plant).npp[f] > loct->npp_max){
				loct->npp_max = (flux->plant).npp[f];
			}
            loct->npp_av[f] = (flux->plant).npp[f];
			
			(flux->soil).lL[f] = (flux->plant).lL[f];			
			(flux->soil).d13c_lL[f] = (flux->plant).d13c_lL[f];
			(flux->soil).d14c_lL[f] = (flux->plant).d14c_lL[f];
			
			/* soil processes *************/
			soil_processes(grid, loct, &(echar->soil), &(mass->soil), &(flux->soil));
			flux->sr[f] = loct->c3ptn[f]*((flux->c3).rrm[f]+(flux->c3).rrg[f]) + 
						loct->c4ptn[f]*((flux->c4).rrm[f]+(flux->c4).rrg[f]) + 
						(flux->soil).hr[f];
			
			if(NECB_DOC == 1){
				(mass->soil).msl -= (flux->soil).doc_boyer[f]/1000000.0;
				if((mass->soil).msl < 0.0){
					(mass->soil).msl = 0.0;
				}
			}
			
			/* fertilizaer input for croplands: revised by A.Ito (2009/06/04) */
			/* NH4:NO3 ratio is based on inventories */
            /* this routine may not be activated when using REPLACE_OLSON_CROP option */
			if((echar->soil).v_type == 1 && (grid->veg_olson==29 || grid->veg_olson==30 ||
											 grid->veg_olson==31 || grid->veg_olson==32)){
				(flux->soil).n_fertin[grid->m] = loct->n_frtlz_in * 1000.0 * f_fert;
				(mass->soil).n_no3 += loct->n_frtlz_in * 0.2 * 1000.0 * f_fert;
				(mass->soil).n_nh4 += loct->n_frtlz_in * 0.8 * 1000.0 * f_fert;
			}else{
				(flux->soil).n_fertin[grid->m] = 0.0;
                (flux->soil).n_manurein[grid->m] = 0.0;
			}
            
            if(EX_NFERT == 101){
                if((echar->soil).v_type == 2){
                    /* 101: Nishina ESSD data: 2017/02/13 by A.Ito */
                    
                    if(grid->f_crop_con > 0.0){
                        icrop = 1.0 / grid->f_crop_con;
                    }else{
                        icrop = 0.0;
                    }
                    
                    y_nin = grid->niny;
                    if(EX_NFERT_SA == 1 || EX_NFERT_SA == 5 || EX_NFERT_SA == 6){
                        (flux->soil).n_fertin[grid->m] = icrop * (grid->nin_no3[0][grid->m]
                                                        + grid->nin_nh4[0][grid->m]) * 1000.0;
                        (mass->soil).n_no3 += icrop * grid->nin_no3[0][grid->m] * 1000.0;
                        (mass->soil).n_nh4 += icrop * grid->nin_nh4[0][grid->m] * 1000.0;
                    }else{
                        if(grid->niny < 1960){
                            (flux->soil).n_fertin[grid->m] = icrop * (grid->nin_no3[0][grid->m]
                                                        + grid->nin_nh4[0][grid->m]) * 1000.0;
                            (mass->soil).n_no3 += icrop * grid->nin_no3[0][grid->m] * 1000.0;
                            (mass->soil).n_nh4 += icrop * grid->nin_nh4[0][grid->m] * 1000.0;
                        }else if(grid->niny >= 1960 && grid->niny <= 2009){
                            (flux->soil).n_fertin[grid->m] = icrop * (grid->nin_no3[y_nin - 1960][grid->m]
                                                        + grid->nin_nh4[y_nin - 1960][grid->m]) * 1000.0;
                            (mass->soil).n_no3 += icrop * grid->nin_no3[y_nin - 1960][grid->m] * 1000.0;
                            (mass->soil).n_nh4 += icrop * grid->nin_nh4[y_nin - 1960][grid->m] * 1000.0;
                        }else if(grid->niny > 2009){
                            (flux->soil).n_fertin[grid->m] = icrop * (grid->nin_no3[49][grid->m]
                                                        + grid->nin_nh4[49][grid->m]) * 1000.0;
                            (mass->soil).n_no3 += icrop * grid->nin_no3[49][grid->m] * 1000.0;
                            (mass->soil).n_nh4 += icrop * grid->nin_nh4[49][grid->m] * 1000.0;
                        }
                    }
                    
                    /* sensitivity run: 2017/06/24 by A.Ito */
                    if(EX_NFERT_SA == 2 || EX_NFERT_SA == 5 || EX_NFERT_SA == 6){
                        /* fixed at 1960 */
                        base_nin = icrop * (grid->nin_no3[0][grid->m] + grid->nin_nh4[0][grid->m]) * 1000.0;
                    }else if(EX_NFERT_SA == 1){
                        /* change only for mature */
                        if(y_nin < 1960){
                            base_nin = icrop * (grid->nin_no3[0][grid->m]
                                            + grid->nin_nh4[0][grid->m]) * 1000.0;
                        }else if(y_nin >= 1960 && y_nin <= 2009){
                            base_nin = icrop * (grid->nin_no3[y_nin - 1960][grid->m]
                                            + grid->nin_nh4[y_nin - 1960][grid->m]) * 1000.0;
                        }else if(y_nin > 2009){
                            base_nin = icrop * (grid->nin_no3[49][grid->m]
                                            + grid->nin_nh4[49][grid->m]) * 1000.0;
                        }
                    }

                    if(grid->nfert_potter > 0.0 || grid->nmanure_potter > 0.0){
                        
                        if(grid->nfert_potter > 0.0){
                            aa = grid->nmanure_potter / grid->nfert_potter;
                            bb = aa * ((flux->soil).n_fertin[grid->m]/1000.0);
                            
                            if(bb > 20.0*(grid->nmanure_potter * MDN[grid->m] / YDN) ){
                                bb = 20.0*(grid->nmanure_potter * MDN[grid->m] / YDN);
                            }
                        }else{
                            bb = icrop * grid->nmanure_potter * MDN[grid->m] / YDN;
                        }
                        
                        if(bb < 0.0){
                            bb = 0.0;
                        }
                        
                        /* No manure: 2017/07/10 by A.Ito */
                        if(EX_NFERT_SA == 7){
                            bb = 0.0;
                        }

                        loct->n_manure_in = bb * 1000.0;
                        (flux->soil).n_manurein[grid->m] = loct->n_manure_in;
                        (mass->soil).n_lttr += loct->n_manure_in;
                    }else{
                        (flux->soil).n_manurein[grid->m] = loct->n_manure_in = 0.0;
                    }
                }
            }else{
                if((echar->soil).v_type == 2){
                    (flux->soil).n_fertin[grid->m] = loct->n_frtlz_in * 1000.0 * f_fert;
                    (mass->soil).n_no3 += loct->n_frtlz_in * 0.2 * 1000.0 * f_fert;
                    (mass->soil).n_nh4 += loct->n_frtlz_in * 0.8 * 1000.0 * f_fert;
                    
                    /* 2016/10/20 by A.Ito */
                    (flux->soil).n_manurein[grid->m] = loct->n_manure_in * 1000.0 * f_fert;
                    (mass->soil).n_lttr += loct->n_manure_in * 1000.0 * f_fert;
                }
            }
			
			/* CH4 oxydation **************/
			f_ch4oxy_ridgewell(grid, loct, flux);
			f_ch4oxy_casa(grid, loct, flux);
			f_ch4oxy_delgrosso(grid, loct, flux);
			f_ch4oxy_curry(grid, loct, flux);
			
			/* CH4 emission (wetlands) */
			f_ch4_emit_cao(grid, loct, flux);
			
			/* N2O emission */
			f_n2o_emit_ngas(grid, loct, mass, flux);
			f_n2o_emit_casa(grid, loct, mass, flux);
			
			/* ecosystem mass balance *****************/	
			/* net ecosystem production */
			flux->nep[f] = (flux->plant).npp[f] - (flux->soil).hr[f];
			flux->er[f] = (flux->plant).ar[f] + (flux->soil).hr[f];
			/* total ecosystem carbon storage */
			mass->total[f] = (mass->c3).plant[f]*loct->c3ptn[f] + 
							(mass->c4).plant[f]*loct->c4ptn[f] + (mass->soil).soil[f];
			/** net carbon balance taking crop harvest into account **/
			flux->ncb[f] = flux->nep[f] - (flux->plant).net_crop[f];
			
			/* carbon isotope */
			f_cisotope_efflux(grid, loct, mass, flux);
			
			/* nitrogen budget */
			n_budget(grid, loct, mass, flux);
						
			/* annual average plant mass */
			plantmass += (mass->plant).plant[f] * MDN[f]/YDN;
			
			if((echar->c3).v_type == 1){
				/* corrected: A.Ito and E.Kato (2009/08/16) */
				/* if(grid->veg_olson==29 || grid->veg_olson==30 || grid->veg_olson== 31
						|| grid->veg_olson==32){
					ann_nep += flux->ncb[f];
				}else{ */
					ann_nep += flux->nep[f];
				/* } */
			}else if((echar->c3).v_type == 2){
				ann_nep += flux->ncb[f];
			}
			
			/*** acclimation ***/
			(echar->soil).ft0_l[f] = (echar->soil).ft_l[f];
			(echar->soil).ft0_h[f] = (echar->soil).ft_h[f];
			(echar->soil).fm0_l[f] = (echar->soil).fm_l[f];
			(echar->soil).fm0_h[f] = (echar->soil).fm_h[f];
			flux->lL0[f] = (flux->plant).lL[f];
            
            /* baseline soil water: 2014/12/08 by A.Ito */
            loct->b_sw30[f] = loct->sw30;
            loct->b_sww[f] = loct->sww;
		}
		
		/* biomass burning */
		f_biomassburning(grid, loct, mass, flux);
		
		/* erosion */
		f_erosion(grid, loct, echar, mass, flux);
		
		/* empirical model NPP *****************/
		//if(grid->y==0){ /* for the first year */
			npp_empirical(grid, loct, flux);
		//}
        
		if(NECB_ERSN == 1){
            if(PARAM_PTB == 20){
                prm_ensen = 1.0 + 0.3 * f_pert[7];
            }else{
                prm_ensen = 1.0;
            }

			(mass->soil).ltr -= flux->erod_carbon * (prm_ensen * 0.20);
			if((mass->soil).ltr < 0.0){
				(mass->soil).ltr = 0.0;
			}
		}
		
		/* terminal conditions ****************************/
        if(ISIMIP_RUN == 0){
            if(nn < 200){	
                /* continued */
                ann_nep = 10.0; 
            }else if(nn>=200 && nn<term_time){	
                ann_nep = fabs(ann_nep); /**** 1. sufficiently stabilized ****/	
            }else{  /*  if(nn>=term_time) */
                break; /**** 3. stop by 2000 years ****/	
            }
        }else if(ISIMIP_RUN == 1 || ISIMIP_RUN == 2 || ISIMIP_RUN == 3 || ISIMIP_RUN == 4){
            /* spin-up 3000 years (30 x 100 times): 2012/07/02 by A.Ito */
            ann_nep = 10.0;
            if(nn == 3000){
                ann_nep = 0.0;
            }
        }
		
		if(plantmass<0.0 || plantmass>=500.0){
			printf("!!! BAD plant biomass: %lf\n", plantmass);
			vanish(mass, flux);  /*** 2. excluding abnormal estimates ***/
			break;
		}
		loct->time = nn; /* simulation time of carbon budget */
		nn++;
	}
	/* end of stabilization loop ***********************************************/
	
	/* Walter & Heimann (paddy) */
	if(CH4_WH==1 && grid->f_wetland>0.0 && loct->v_type == 1){
		for(g=0;g<4;g++){
			for(f=0;f<ASTEP;f++){
				(flux->soil).ch4_wetland_wh_plant[f] = 0.0;
				(flux->soil).ch4_wetland_wh_ebull[f] = 0.0;
				(flux->soil).ch4_wetland_wh_diff[f] = 0.0;
				(flux->soil).ch4_wetland_wh_release[f] = 0.0;
			}
			for(f=0;f<ASTEP;f++){
				grid->m = f;
				/* water-logged */
				f_ch4_emit_walter(1, grid, loct, flux);
				/* drainage */
				f_ch4_emit_walter(2, grid, loct, flux);
			}
		}
	}else{
		for(f=0;f<ASTEP;f++){
			(flux->soil).ch4_wetland_wh_plant[f] = 0.0;
			(flux->soil).ch4_wetland_wh_ebull[f] = 0.0;
			(flux->soil).ch4_wetland_wh_diff[f] = 0.0;
			(flux->soil).ch4_wetland_wh_release[f] = 0.0;
		}
	}
	if(CH4_WH == 1 && grid->f_paddy>0.0 && loct->v_type == 2){
		for(g=0;g<4;g++){
			for(f=0;f<ASTEP;f++){
				(flux->soil).ch4_paddy_wh_plant[f] = 0.0;
				(flux->soil).ch4_paddy_wh_ebull[f] = 0.0;
				(flux->soil).ch4_paddy_wh_diff[f] = 0.0;
				(flux->soil).ch4_paddy_wh_release[f] = 0.0;
			}
			for(f=0;f<ASTEP;f++){
				grid->m = f;
				/* water-logged */
				f_ch4_emit_walter(3, grid, loct, flux);
				/* drainage */
				f_ch4_emit_walter(4, grid, loct, flux);
			}
		}
	}else{
		for(f=0;f<ASTEP;f++){
			(flux->soil).ch4_paddy_wh_plant[f] = 0.0;
			(flux->soil).ch4_paddy_wh_ebull[f] = 0.0;
			(flux->soil).ch4_paddy_wh_diff[f] = 0.0;
			(flux->soil).ch4_paddy_wh_release[f] = 0.0;
		}
	}
		
	if(NECB_CH4 == 1){
		for(f=0;f<ASTEP;f++){
			(mass->soil).msl += grid->f_upland * (flux->soil).ch4oxy_curry[f] * 0.00001
					- grid->f_paddy * ((flux->soil).ch4_paddy_wh_plant[f] + (flux->soil).ch4_paddy_wh_ebull[f] + 
						(flux->soil).ch4_paddy_wh_diff[f] + (flux->soil).ch4_paddy_wh_release[f]) * 0.00001
					- grid->f_wetland * ((flux->soil).ch4_wetland_wh_plant[f] + (flux->soil).ch4_wetland_wh_ebull[f] + 
						(flux->soil).ch4_wetland_wh_diff[f] + (flux->soil).ch4_wetland_wh_release[f]) * 0.00001;
			
			if((mass->soil).msl < INT_C){
				(mass->soil).msl = INT_C;
			}
		}
	}
	
	/* land use change **********************************/
	if(loct->v_type == 1){
		f_luc_emit(grid, mass, flux);
	}else{
		flux->lu_detr = 0.0;
		flux->lu_conv = 0.0;
		flux->lu_ten = 0.0;
		flux->lu_hund = 0.0;

        flux->lu_fol = 0.0;
        flux->lu_stm = 0.0;
        flux->lu_rot = 0.0;
        flux->lu_ltr = 0.0;
        flux->lu_msl = 0.0;
	}
	
    f_nat = 1.0 - grid->f_crop_con;
    if(f_nat > 0.0){
        iweight = 1.0 / f_nat;
    }else{
        iweight = 1.0;
    }
    avc3 = 0.0;
    for(f=0;f<ASTEP;f++){
        avc3 += (loct->c3ptn[f] * MDN[f]/YDN);
    }
    if(avc3 > 0.0){
        iweight3 = iweight * (1.0 / avc3);
    }else{
        iweight3 = iweight;
    }

	/* wood harvest: 2010/10/15 by A.Ito *****************/
	total_hvst = 0.0;
	if((mass->c3).v_type == 1){
        /* revised (after comments by E.Kato): 2013/10/02 by A.Ito */
        
        if(LANDUSE ==9 || LANDUSE ==10 || LANDUSE ==11 || LANDUSE ==12 || LANDUSE ==13
                    || LANDUSE ==14){
            dyr = 1900 - FDY_LUC;
        }else{
            dyr = 1900 - FDY_LUC;
        }
		
        /* from total grid */
        /* total_hvst = grid->hvst_p1[dyr] + grid->hvst_p2[dyr] + grid->hvst_s1[dyr]
                    + grid->hvst_s2[dyr] + grid->hvst_s3[dyr]; */

        total_hvst = grid->hvst_p1[dyr] + grid->hvst_s1[dyr];

        /* parameter ensemble: 2014/11/19 by A.Ito */
        prm_ensen = 1.0;
        if(PARAM_PTB == 9){
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
            prm_ensen = 1.0 + 0.3 * f_pert[1];
        }
    
		total_hvst *= 1.0/1000.0 * 1.0/grid->area * prm_ensen;
        
        if(NECB_WHVST == 1){
            if((mass->c3).stm > (total_hvst + INT_C)){
            
                if((mass->c3).stm > (total_hvst * iweight3 + INT_C)){
                    (mass->c3).stm -= total_hvst * iweight3;
                    flux->hvst_wood = total_hvst * iweight3;
                }else{
                    (mass->c3).stm -= total_hvst;
                    flux->hvst_wood = total_hvst;
                }
                
            }else{
                flux->hvst_wood = total_hvst - INT_C;
                (mass->c3).stm = INT_C;
            }
            
            if((mass->c3).stm < INT_C){
                (mass->c3).stm = INT_C;
            }
        }
	}else{
		flux->hvst_wood = 0.0;
	}

	/* net biome production (added by A.Ito: 2010/01/20) */
	for(f=0;f<ASTEP;f++){
		flux->nbp[f] = flux->nep[f];
		
        /* altered: 2018/10/16 by A.Ito */
        if((mass->c3).v_type == 1 && NECB_LUC == 1){
            flux->nbp[f] -= iweight * (flux->lu_conv/(double)ASTEP + flux->lu_ten/(double)ASTEP + flux->lu_hund/(double)ASTEP);
        }

        if(NECB_WHVST == 1){
            flux->nbp[f] -= flux->hvst_wood /(double)ASTEP;
        }
		
		if(NECB_BB == 1){
            /* revised (after comments by E.Kato): 2013/10/02 by A.Ito */
            
			flux->nbp[f] -= (flux->bb_co2_litter[f] + flux->bb_co2_leaf[f] 
							 + flux->bb_co2_wood[f] + flux->bb_co2_root[f])/1000.0*12.0/44.0;

			flux->nbp[f] -= (flux->bb_co_litter[f] + flux->bb_co_leaf[f] 
							 + flux->bb_co_wood[f] + flux->bb_co_root[f])/1000.0*12.0/28.0;

			flux->nbp[f] -= (flux->bb_ch4_litter[f] + flux->bb_ch4_leaf[f]
							 + flux->bb_ch4_wood[f] + flux->bb_ch4_root[f])/1000.0*12.0/16.0;

			flux->nbp[f] -= (flux->bb_bc_litter[f] + flux->bb_bc_leaf[f]
							 + flux->bb_bc_wood[f] + flux->bb_bc_root[f])/1000.0;
		}
        
        if(NECB_DOC == 1){
            /* revised (after comments by E.Kato): 2013/10/02 by A.Ito */
            
            flux->nbp[f] -= (flux->soil).doc_boyer[f]/1000000.0;
        }
        
        if(NECB_CH4 == 1){
            /* revised (after comments by E.Kato): 2013/10/02 by A.Ito */
            
            flux->nbp[f] += 12.0/16.0 * (grid->f_upland * (flux->soil).ch4oxy_curry[f] * 0.00001
                    - grid->f_paddy * ((flux->soil).ch4_paddy_wh_plant[f] + (flux->soil).ch4_paddy_wh_ebull[f] + 
                                (flux->soil).ch4_paddy_wh_diff[f] + (flux->soil).ch4_paddy_wh_release[f]) * 0.00001
                    - grid->f_wetland * ((flux->soil).ch4_wetland_wh_plant[f] + (flux->soil).ch4_wetland_wh_ebull[f] + 
                                (flux->soil).ch4_wetland_wh_diff[f] + (flux->soil).ch4_wetland_wh_release[f]) * 0.00001);
        }
        
        if(NECB_ERSN == 1){
             if(PARAM_PTB == 20){
                prm_ensen = 1.0 + 0.3 * f_pert[7];
            }else{
                prm_ensen = 1.0;
            }

           /* revised (after comments by E.Kato): 2013/10/02 by A.Ito */
        
            flux->nbp[f] -= flux->erod_carbon * (prm_ensen * 0.20) /(double)ASTEP;
        }
        
        if(NECB_BVOC == 1){
            /* revised (after comments by E.Kato): 2013/10/02 by A.Ito */
        
            flux->nbp[f] -= (flux->voc_isopr_g97[f] + flux->voc_monotrp_g97[f] + flux->voc_methanl_g97[f] +
                    flux->voc_acetone_g97[f] + flux->voc_actaldhd_g97[f] + flux->voc_frmardhd_g97[f] +
                    flux->voc_formacd_g97[f] + flux->voc_acetacd_g97[f] + flux->voc_co_g97[f] +
                    flux->voc_afarnesene[f] + flux->voc_bcaryophyllene[f] + flux->voc_othersesqui[f])*10000.0/1000000.0/1000000.0;
        }
        
        if(NECB_CROP == 1){
            /* revised (after comments by E.Kato): 2013/10/02 by A.Ito */
            flux->nbp[f] -= 1.0 * (flux->plant).net_crop[f]; /* ! hvst is positive */
        }
	}
	
	/* history data */
	f_set_history_data(grid->simy - (FSY_HIST-1), grid, loct, mass, flux);
		
	/** output initial stable state **/
    if(OUTPUT_CARBON1 == 1){
	    publish_cbud(grid, loct, echar, mass, flux, fp_o[0]); /* */
    }
	
	/* output */
	f_output_result(FSY_HIST-1, grid, loct, echar, mass, flux, fp_o); /* */
}

