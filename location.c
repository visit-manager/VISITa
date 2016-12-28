/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* initialize environmental locations */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

#define TERM_HYD 0.1

/* initialization of climatic conditions (Primary data) ******************************/
void f_init_clim(
	struct Grid *grid
){
	short h, impex_p, impex_t;
	double aaa, bbb, tmp_var, pre_var;
	
	/* in 1950 :311 ppmv*/
	/* in 1990 : 352.7 ppmv*/
	grid->co2y = BGY_CO2Y; 
	if(CO2S == 7){
		grid->co2y = 2081; /* in 2081 : 700 ppmv*/
	}
    if(CC_CD == 5){
        grid->co2y = 2000;
    }
	grid->climy = BGY_CLIM;
	grid->lucy = BGY_CLIM;

	for(h=0;h<ASTEP;h++){
		grid->m = h;
		
		/* ambient CO2 condition *******/
		f_co2_trend(grid);

		/* climatic conditions */
		/* NCEP/NCAR mean climate data ***************************************/
		grid->tmp_sfc[h] = grid->tmp_sfc_a[h]; /* temperature, degree C */
		grid->tmp_2m[h] = grid->tmp_2m_a[h];
		grid->tmp10_soil[h] = grid->tmp10_soil_a[h];
		grid->tmp200_soil[h] = grid->tmp200_soil_a[h];
		grid->tcdc_clm[h] = grid->tcdc_clm_a[h]; /* cloudiness, 0-1 fraction */
		grid->prate_sfc[h] = grid->prate_sfc_a[h]; /* precipitation, mm / month */
		
		/*** UEA/CRU data ***/
		if(grid->flag_histdata == 1){ /* for grids data are available */
			/* temporary */
			grid->tmp_sfc[h] = grid->hist_tmp_b[h] + (grid->tmp_sfc_a[h] - grid->tmp_2m_a[h]);
			grid->tmp10_soil[h] = grid->hist_tmp_b[h] + (grid->tmp10_soil_a[h] - grid->tmp_2m_a[h]);
			grid->tmp200_soil[h] = grid->hist_tmp_b[h] + (grid->tmp200_soil_a[h] - grid->tmp_2m_a[h]);
			grid->tmp_2m[h] = grid->hist_tmp_b[h];
			grid->tcdc_clm[h] = grid->hist_cld_b[h];
			grid->prate_sfc[h] = grid->hist_pre_b[h];
			
			/* long-term average */
			grid->tmp_sfc_a[h] = grid->hist_tmp_b[h] + (grid->tmp_sfc_a[h] - grid->tmp_2m_a[h]);
			grid->tmp10_soil_a[h] = grid->hist_tmp_b[h] + (grid->tmp10_soil_a[h] - grid->tmp_2m_a[h]);
			grid->tmp200_soil_a[h] = grid->hist_tmp_b[h] + (grid->tmp200_soil_a[h] - grid->tmp_2m_a[h]);
			grid->tmp_2m_a[h] = grid->hist_tmp_b[h];
			grid->tcdc_clm_a[h] = grid->hist_cld_b[h];
			grid->prate_sfc_a[h] = grid->hist_pre_b[h];
		}
		
		/* alternative precipitation data *************************/
		if(grid->prec_sub_a[h] >= 0.0){
			grid->prate_sfc[h] = grid->prec_sub_a[h];
		}else if(grid->prec_sub_a[h] < 0.0){
			grid->prate_sfc[h] = grid->prate_sfc_a[h];
		}
		
		grid->spfh_2m[h] = grid->spfh_2m_a[h];
		
		aaa = grid->ugrd_10m_a[h] * grid->ugrd_10m_a[h];
		bbb = grid->vgrd_10m_a[h] * grid->vgrd_10m_a[h];
		grid->wnd_10m[h] = sqrt(aaa + bbb);
        
        /*******************************************/
        /* IMPRESSIONS IRS: 2015/07/17 by A.Ito */
        if(IMPRESSIONS_RUN == 1){
        
            if(SCENARIO_ID<6001 || SCENARIO_ID>7000){
                printf("BAD experimental ID\n");
                exit(1);
            }
            
            /*******/
            impex_t = (short)(SCENARIO_ID%50);
            
            switch(impex_t){
                case 1: tmp_var = -3.0; break;
                case 2: tmp_var = -2.0; break;
                case 3: tmp_var = -1.0; break;
                case 4: tmp_var = 0.0; break;
                case 5: tmp_var = 1.0; break;
                case 6: tmp_var = 2.0; break;
                case 7: tmp_var = 3.0; break;
                case 8: tmp_var = 4.0; break;
                case 9: tmp_var = 5.0; break;
                case 10: tmp_var = 6.0; break;
                case 11: tmp_var = 7.0; break;
                case 12: tmp_var = 8.0; break;
                case 13: tmp_var = 9.0; break;
                case 14: tmp_var = 10.0; break;
                case 15: tmp_var = 11.0; break;
                default: tmp_var = 0.0; break;
            }
            
            grid->tmp_sfc[h] += tmp_var;
            grid->tmp_2m[h] += tmp_var;
            grid->tmp10_soil[h] += tmp_var;
            grid->tmp200_soil[h] += tmp_var;
            
            /*******/
            impex_p = (short)((SCENARIO_ID - 6000)/50);
            
            switch(impex_p){
                case 0: pre_var = 0.4; break;
                case 1: pre_var = 0.5; break;
                case 2: pre_var = 0.6; break;
                case 3: pre_var = 0.7; break;
                case 4: pre_var = 0.8; break;
                case 5: pre_var = 0.9; break;
                case 6: pre_var = 1.0; break;
                case 7: pre_var = 1.1; break;
                case 8: pre_var = 1.2; break;
                case 9: pre_var = 1.3; break;
                case 10: pre_var = 1.4; break;
                case 11: pre_var = 1.5; break;
                case 12: pre_var = 1.6; break;
                default: pre_var = 1.0; break;
            }
            
            grid->prate_sfc[h] *= pre_var;
        }

        grid->tmp_soil_am += grid->tmp200_soil_a[h] * MDN[h] / 365.0;

		/* solar decrination and solar height */
		grid->sl_dec[h] = f_solar_decl(grid);
		grid->sl_hgt[h] = f_solar_hgt(grid);
		
		/* radiation fluxes and day-length */
		grid->dlen[h] = f_day_length(grid);		
		grid->top_rad[h] = f_top_rad(grid, 0); 	

		grid->gl_rad[h] = f_gl_rad(grid); 	
		grid->par[h] = f_par(grid); 
		
		grid->rad_a[h] = grid->gl_rad[h];
		grid->par_a[h] = grid->par[h];
        
		/* sensitivity analysis *************************/
		if(TM == 1){
			grid->tmp_sfc[h] += 1.0;
			grid->tmp_2m[h] += 1.0;
			grid->tmp10_soil[h] += 1.0;
			grid->tmp200_soil[h] += 1.0;
		}else if(TM == 5){
			grid->tmp_sfc[h] -= 1.0;
			grid->tmp_2m[h] -= 1.0;
			grid->tmp10_soil[h] -= 1.0;
			grid->tmp200_soil[h] -= 1.0;
		}
		if(PR == 1){
			grid->prate_sfc[h] *= 1.1;
		}else if(PR == 5){
			grid->prate_sfc[h] *= 0.9;
		}
	}
}

/* location conditions derived from the primary data (Secondary data1) *******************/
void f_init_loct(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Flux *flux, 
	struct Echar *echar
){
	long h, nn;
	double tem_grow, d_smc_a;
	double ftmp10b, ftmp200b, ftmp10, ftmp200;
	extern double MDN[12];
	
	/* long-term average annual climate condition */
	tem_grow = 0.0;
	grid->tmp_sfc_am = 0.0;
	grid->tmp_sfc_mx = -100.0;
	grid->tmp_sfc_mn = 100.0;
	grid->gp_atem = 0.0;
	grid->prate_sfc_ann = 0.0;
	nn = 0;
	for(h=0;h<ASTEP;h++){
		/* annual mean temperature */
		grid->tmp_sfc_am += grid->tmp_sfc[h]/12.0; 
		/* annual maximum */
		grid->tmp_sfc_mx = (grid->tmp_sfc[h]>grid->tmp_sfc_mx)?grid->tmp_sfc[h]:grid->tmp_sfc_mx; 
		/* annual minimum */
		grid->tmp_sfc_mn = (grid->tmp_sfc[h]<grid->tmp_sfc_mn)?grid->tmp_sfc[h]:grid->tmp_sfc_mn; 
		/* annual total precipitation */
		grid->prate_sfc_ann += grid->prate_sfc[h]; 
		
		if(grid->tmp_sfc[h]>5.0){
			nn += MDN[h];
			grid->gp_atem += (grid->tmp_sfc[h]-5.0)*MDN[h];
			tem_grow += grid->tmp_sfc[h]*MDN[h];
			grid->gp_pre += grid->prate_sfc[h]; 
		}
	}
	grid->gp_tem = (nn>0)?tem_grow/(double)nn:0.0;

	/* portion of C3 and C4 plant ***************/
	c34composition((echar->c3).v_type, grid, loct);
	
	/* sensitivity analysis */
	if(T_D==1 || T_D==2 || T_D==3 || T_D==4){
		ftmp10b = ftmp200b = ftmp10 = ftmp200 = 0.0;
		for(h=0;h<ASTEP;h++){
			if(grid->tmp10_soil[h] > -20.0){
				ftmp10b += 0.05 + 0.95*exp(308.56*(1.0/56.02-1.0/(grid->tmp10_soil[h]+46.02))); 
				if(T_D==0){
					ftmp10 += 0.05 + 0.95*exp(308.56*(1.0/56.02-1.0/(grid->tmp10_soil[h]+46.02)));
				}else if(T_D==1){
					ftmp10 += 0.05 + 0.95*exp(308.56*1.3*(1.0/56.02-1.0/(grid->tmp10_soil[h]+46.02)));
				}else if(T_D==2){
					ftmp10 += 0.05 + 0.95*exp(308.56*0.7*(1.0/56.02-1.0/(grid->tmp10_soil[h]+46.02)));
				}else if(T_D==3){
					ftmp10 += 0.05 + 0.95*exp(308.56*(1.0/56.02-1.0/(grid->tmp10_soil[h]+46.02*1.3)));
				}else if(T_D==4){
					ftmp10 += 0.05 + 0.95*exp(308.56*(1.0/56.02-1.0/(grid->tmp10_soil[h]+46.02*0.7)));
				}
			}else{
				ftmp10b += 0.05;
				ftmp10 += 0.05;
			}
			if(grid->tmp200_soil[h] > -20.0){
				ftmp200b += 0.05 + 0.95*exp(308.56*(1.0/56.02-1.0/(grid->tmp200_soil[h]+46.02)));
				if(T_D==0){
					ftmp200 += 0.05 + 0.95*exp(308.56*(1.0/56.02-1.0/(grid->tmp200_soil[h]+46.02)));
				}else if(T_D==1){
					ftmp200 += 0.05 + 0.95*exp(308.56*1.3*(1.0/56.02-1.0/(grid->tmp200_soil[h]+46.02)));
				}else if(T_D==2){
					ftmp200 += 0.05 + 0.95*exp(308.56*0.7*(1.0/56.02-1.0/(grid->tmp200_soil[h]+46.02)));
				}else if(T_D==3){
					ftmp200 += 0.05 + 0.95*exp(308.56*(1.0/56.02-1.0/(grid->tmp200_soil[h]+46.02*1.3)));
				}else if(T_D==4){
					ftmp200 += 0.05 + 0.95*exp(308.56*(1.0/56.02-1.0/(grid->tmp200_soil[h]+46.02*0.7)));
				}
			}else{
				ftmp200b += 0.05;
				ftmp200 += 0.05;
			}
		}
		(echar->soil).rl0 *= ftmp10b/ftmp10;
		(echar->soil).rh0 *= ftmp200b/ftmp200;
	}
	
	/* water condition - Sim-HYDRO *****************/
	loct->sw30 = 0.5*grid->field_cap1;
	loct->sww = 0.5*grid->field_cap2; 
	loct->snwa = 0.0;
	loct->m_m_pre = 0.24;
	nn = 0; 
	d_smc_a = 10.0;
	while(d_smc_a > TERM_HYD){
		d_smc_a = loct->sww;
		for(h=0;h<ASTEP;h++){
			grid->m = h;
		/*	(mass->c3).lai[grid->m]=loct->c3ptn[grid->m];
			(mass->c4).lai[grid->m]=loct->c4ptn[grid->m]; */
			
			co2_in_canopy(grid, loct, mass, flux);
			
			f_dyn_loct(grid, loct, mass, echar);
		}
		d_smc_a = fabs(loct->sww - d_smc_a);
		loct->time_hyd = nn; /* simulation time of carbon budget */
		nn++;
		
		if(nn < 6){
			d_smc_a = 10.0; /* at least 5 years */
		}
		if(nn > 50){
			break;
		}
	}
	
	/* average fertilizer-N input for each county, kg N ha-1 yr-1 */
	n_fertilizer_in(grid, loct);
    
    /* empirical NPP and meta-analysis-based optimal LAI */
    npp_empirical(grid, loct, flux);
}

/* dynamic estimation of environmnetal conditions (Secondary data2) *********************/
void f_dyn_loct(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Echar *echar
){
	long h, offset;
	double alt, k_c, vpres_var, tmp_ann, tmp_var, wet_var;
	
    /* initialization */
    if(grid->m == 0){
        for(h=0;h<ASTEP;h++){
            loct->xx1[h] = 0.0;
            loct->xx2[h] = 0.0;
            loct->xx3[h] = 0.0;
            loct->xx4[h] = 0.0;
            loct->xx5[h] = 0.0;
            loct->xx6[h] = 0.0;
            loct->xx7[h] = 0.0;
            loct->xx8[h] = 0.0;
            loct->xx9[h] = 0.0;
        }
        loct->est_maxlai = 0.0;
    }
    for(h=0;h<DSTEP;h++){
        loct->ppfd_h[h] = 0.0;
        loct->ppfdb_h[h] = 0.0;
        loct->ppfdd_h[h] = 0.0;
    }

	/* solar constant sensitivity */
	if(SC==3 || SC==4){
		grid->top_rad[grid->m] = f_top_rad(grid, 0); 	
		grid->gl_rad[grid->m] = f_gl_rad(grid); 	
		grid->par[grid->m] = f_par(grid); 
	}else if(SC == 5){
		grid->par[grid->m] = f_par(grid); 
		grid->par[grid->m] += 10.0;
	}
    
	/* radiatin for cal_historical: 1901-2000 */
	if(grid->flag_histdata == 1){
    
        /* daily-mean and hourly: added: 2013/01/10 by A.Ito */
        loct->grad_d[grid->m] = 0.0;
        for(h=0;h<DSTEP;h++){
            grid->h = h;
        
            /* atmosphere-top, hourly */
            grid->top_rad[grid->m] = f_top_rad(grid, -180 + h*15);
            /* surface, hourly */
            grid->gl_rad[grid->m] = f_gl_rad(grid);
            /* daily mean */
            loct->grad_d[grid->m] += f_gl_rad(grid)/(double)DSTEP;
            
            /* hourly PAR, micro mol photon m-2 s-1 */
            f_par_h(grid, loct);
        }
        
        /* midday */
        grid->top_rad[grid->m] = f_top_rad(grid, 0);
		grid->gl_rad[grid->m] = f_gl_rad(grid); 
		grid->par[grid->m] = f_par(grid);
	}
	
	if(grid->m == 0){
		grid->tmp_sfc_am = 0.0;
		grid->tmp_sfc_mx = -100.0;
		grid->tmp_sfc_mn = 100.0;
		grid->prate_sfc_ann = 0.0;
		grid->tmp_soil_mean = 0.0;
		for(h=0;h<ASTEP;h++){
			/* annual mean temperature */
			grid->tmp_sfc_am += grid->tmp_sfc[h] * MDN[h]/365.0;
			/* annual maximum */
			grid->tmp_sfc_mx = (grid->tmp_sfc[h]>grid->tmp_sfc_mx)?grid->tmp_sfc[h]:grid->tmp_sfc_mx; 
			/* annual minimum */
			grid->tmp_sfc_mn = (grid->tmp_sfc[h]<grid->tmp_sfc_mn)?grid->tmp_sfc[h]:grid->tmp_sfc_mn; 
			/* annual total precipitation */
			grid->prate_sfc_ann += grid->prate_sfc[h]; 
			/* annual mean soil temperature */
			grid->tmp_soil_mean += grid->tmp10_soil[h] * MDN[h]/365.0;
		}
		loct->cum_dprec = 0.0;
        
        /* change in wetland area due to permafrost melting: 2012/10/26 by A.Ito ****/
        if(VAR_PFMWET == 1){
            
            tmp_ann = 0.0;
            for(h=0;h<ASTEP;h++){
                tmp_ann = grid->tmp_2m[h]*MDN[h]/365.0;
            }
            
            tmp_var = tmp_ann - grid->tmp_base_permaforst;
            
            wet_var = 0.0;
            if(tmp_ann > -2.0 && tmp_var > 0.0 ){
                switch(grid->type_permaforst){
                    case 1: case 5: case 9: case 13: case 17:
                        /* continuous permafrost */
                        wet_var = 0.05 * tmp_var;
                        break;
                    case 2: case 6: case 10: case 14: case 18:
                        /* discontinuous permafrost */
                        wet_var = 0.03 * tmp_var;
                        break;
                    case 3: case 7: case 11: case 15: case 19:
                        /* sporadic permafrost */
                        wet_var = 0.01 * tmp_var;
                        break;
                    case 4: case 8: case 12: case 16: case 20:
                        /* isolated permafrost */
                        wet_var = 0.003 * tmp_var;
                        break;
                    default:
                        /* no permafrost */
                        wet_var = 0.0;
               }
            }
            if(wet_var < 0.0){
                wet_var = 0.0;
            }
            
            grid->f_wetland += wet_var;
            
            if(grid->f_wetland > 1.0){
                grid->f_wetland = 1.0;
            }
        }
	}
	/* cumulative precipitation anomaly: 2010/07/23 by A.Ito */
	/* if(grid->climy<2005){ */
		loct->cum_dprec += (grid->prate_sfc[grid->m] - grid->hist_pre_b[grid->m]);
	/* } */
	
	/** air conditions **/
	/* altitude */
	alt = (grid->topo>=0.0)?grid->topo:0.0; 
	/* air pressure */
	loct->prsr[grid->m] = 1013.25 * exp(-1.0 * (28.964 * 0.001) * GAC * alt / (UGC * (grid->tmp_2m[grid->m] + ZAT)));
	
	/* saturated vapour pressure, hPa */
	loct->vps[grid->m] = vap_pre_sat(grid); 
	/* slope of saturated vapour pressure as a function of temperature */
	loct->slope_vps[grid->m] = slope_vps(grid);	
	 /* aerodynamic resistance */ 
	loct->r_aero[grid->m] = r_aero(grid);	 

	/* initial soil CH4 concentration */
	for(h=0;h<=(N_SLAYER+1);h++){
		loct->prof_ch4[h] = ach4_a1[grid->co2y - 1750]/1000.0 
			* loct->prsr[grid->m] / (UGC * (grid->tmp10_soil[grid->m] + ZAT));
	}
    
    offset = 0;
    if(ISIMIP_RUN == 3){
        if(grid->phase == 0){
            offset = 0;
        }else if(grid->phase == 1){
            /* skip spin-up data */
            offset = 30;
        }
    }

	if(grid->flag_histdata == 1){
		/* vapour pressure, hPa */
		if(grid->phase == 0){
			/* spin-up */
			loct->vp[grid->m] = grid->hist_vap_b[grid->m];
		}else if(grid->phase == 1){
			if(grid->climy < (BGY_CLIM + DL_HCLIM)){
				/* based on UEA/CRU or ISI-MIP data */
				loct->vp[grid->m] = grid->hist_vap[grid->climy - BGY_CLIM + offset][grid->m];
			}else{
				/* based on NCEP/NCAR */
				vpres_var = grid->ncep_vpres[grid->climy - FDY_NCEP][grid->m][grid->ncep_lat][grid->ncep_lon] 
								- grid->ncep_vpres_b[grid->m][grid->ncep_lat][grid->ncep_lon];
				
				loct->vp[grid->m] = grid->hist_vap_b[grid->m] + vpres_var;
			}
		}else if(grid->phase == 2){
			/* prediction using AOGCM */
			/* loct->vp[grid->m] = loct->prsr[grid->m]*grid->spfh_2m[grid->m]/(0.622 + 0.378*grid->spfh_2m[grid->m]);  */
			
			/* revided by A.Ito (2009/08/17) */
			vpres_var = grid->proj_hum[grid->climy - FDY_GCM-1][grid->m][grid->gcm_row][grid->gcm_col] -
							grid->proj_hum_b[grid->m][grid->gcm_row][grid->gcm_col];
			
			loct->vp[grid->m] = grid->hist_vap_b[grid->m] + vpres_var;
            
            if(ISIMIP_RUN == 2){
                loct->vp[grid->m] = grid->proj_hum[grid->climy - FDY_GCM][grid->m][0][0];
            }
		}
		if(loct->vp[grid->m] < 0.0){
			loct->vp[grid->m] = 0.0;
		}
		
		if(loct->vp[grid->m]<=loct->vps[grid->m]){
			loct->vpd[grid->m] = loct->vps[grid->m] - loct->vp[grid->m];
		}else{
			loct->vpd[grid->m] = 0.0;
		}
	}else{
		/* vapour pressure, hPa */
		loct->vp[grid->m] = loct->prsr[grid->m] * grid->spfh_2m[grid->m]/(0.622 + 0.378*grid->spfh_2m[grid->m]);
		/* vapour pressure deficit */
		loct->vpd[grid->m] = (loct->vps[grid->m] >= loct->vp[grid->m])?loct->vps[grid->m] - loct->vp[grid->m]:0.0;
	}

	/* air density */
	loct->dnsa[grid->m] = air_density(grid, loct); 
	
	/** ecophysiology : ecophysiology.c **/
	/* C3 plants */
	f_ecophysiology(grid, loct, &(echar->c3), &(mass->c3));
	/* C4 plants */
	f_ecophysiology(grid, loct, &(echar->c4), &(mass->c4));
	
	/** net radiation **/
	f_net_rad(grid, loct, mass, echar);
    	
	(mass->plant).lai[grid->m] = (mass->c3).lai[grid->m] * loct->c3ptn[grid->m]
					+ (mass->c4).lai[grid->m] * loct->c4ptn[grid->m];
	loct->lai[grid->m] = (mass->plant).lai[grid->m];
    
	/** hydrological water budget **/
	f_waterbudget(grid, loct, echar);
	
	/* fractional vegetation cover */
	/* light attenuation by non-photosynthetic organs */
	/*
	Kurachi, N., A. Hagihara, and K. Hozumi. 1989. Effect of light interception by 
	non-photosynthetic organs on canopy photosynthetic production. 
	Ecological Research 4:187-197.
	*/
	if(grid->veg_sage >= 1 && grid->veg_sage <= 8){
		k_c = 0.221;	/* Kurachi and Hagihara */
	}else{
		k_c = 0.001;
	}

	loct->f_vegcov[grid->m] = 1.0 - loct->c3ptn[grid->m] * exp(-((echar->c3).eK0 + k_c)*(mass->c3).lai[grid->m]) 
								- loct->c4ptn[grid->m] * exp(-((echar->c4).eK0 + 0.001)*(mass->c4).lai[grid->m]);
	if(loct->f_vegcov[grid->m]<0.0){
		loct->f_vegcov[grid->m] = 0.0;
	}
	if(loct->f_vegcov[grid->m]>1.0){
		loct->f_vegcov[grid->m] = 1.0;
	}
	
	/* CASA water scheme **************/
	loct->m_sw[grid->m] = loct->sw30;
	loct->m_vmc[grid->m] = loct->sw30/300.0;
	loct->m_pet[grid->m] = loct->pm_incep[grid->m] + loct->pm_evp[grid->m] + loct->pm_trn[grid->m];
	
	loct->wfps[grid->m] = ((loct->m_vmc[grid->m]*100.0) / ((1.0 - grid->bulkdens/2.65)*100.0));
	
	if(loct->wfps[grid->m] > 1.5){
		loct->wfps[grid->m] = 1.5;
	}
	if(loct->wfps[grid->m] < 0.05){
		loct->wfps[grid->m] = 0.05;
	}
	/* soil moisture index */
	f_casa_moisture(grid, loct);
	
	/* nitrogen deposition ***********/
	f_n_deposit(grid, loct);
	
	/* decay of 14C: added by A.Ito (2009/06/24) ************************/
	(mass->c3).d14c_fol = f_decay_14c((mass->c3).d14c_fol);
	(mass->c3).d14c_stm = f_decay_14c((mass->c3).d14c_stm);
	(mass->c3).d14c_rot = f_decay_14c((mass->c3).d14c_rot);
	(mass->c4).d14c_fol = f_decay_14c((mass->c4).d14c_fol);
	(mass->c4).d14c_stm = f_decay_14c((mass->c4).d14c_stm);
	(mass->c4).d14c_rot = f_decay_14c((mass->c4).d14c_rot);
	(mass->soil).d14c_ltr = f_decay_14c((mass->soil).d14c_ltr);
	(mass->soil).d14c_msl = f_decay_14c((mass->soil).d14c_msl);
}
