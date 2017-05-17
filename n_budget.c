/*	VISIT: Vegetation Integrative SImulation Tool						*/
/*  Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/*  Carbon cycle, erosion, biomass burning, land-use change,			*/
/*  CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* revised 20 August 2007 */
/* revised 29 August 2007 */
/* revised 31 March 2010 by A.Ito */

#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* N budget **************************************************/
void n_budget(
	struct Grid *grid, 
	struct Loct *loct,
	struct Mass *mass, 
	struct Flux *flux
){
	/*  (flux->soil).d_n2o_ntr_ngas[grid->m] = 0.0;
	(flux->soil).d_n2o_dnt_ngas[grid->m] = 0.0;
	(flux->soil).d_n2_ngas[grid->m] = 0.0;  */
	
/*	(flux->c3).n_alloc_cnpy[grid->m] = (flux->c4).n_alloc_cnpy[grid->m] = 0.0;
	(flux->c3).n_alloc_strg[grid->m] = (flux->c4).n_alloc_strg[grid->m] = 0.0;
	(flux->c3).n_salvage[grid->m] = (flux->c4).n_salvage[grid->m] = 0.0;
	(flux->c3).n_abdn_cnpy[grid->m] = (flux->c4).n_abdn_cnpy[grid->m] = 0.0;
	(flux->c3).n_abdn_strg[grid->m] = (flux->c4).n_abdn_strg[grid->m] = 0.0; */
	
	/* (flux->c3).n_realloc[grid->m] = (flux->c4).n_realloc[grid->m] = 0.0; */
	
	/* plant N ****************************************/
	/* canopy */
	(mass->c3).n_cnpy += (flux->c3).n_alloc_cnpy[grid->m] 
						+ (flux->c3).n_realloc[grid->m] 
						- (flux->c3).n_salvage[grid->m] 
						- (flux->c3).n_abdn_cnpy[grid->m];
	if((mass->c3).n_cnpy < 0.0){
		(mass->c3).n_cnpy = 0.0;
	}
	(mass->c3).n_cnpy_m[grid->m] = (mass->c3).n_cnpy;

	(mass->c4).n_cnpy += (flux->c4).n_alloc_cnpy[grid->m] 
						+ (flux->c4).n_realloc[grid->m] 
						- (flux->c4).n_salvage[grid->m] 
						- (flux->c4).n_abdn_cnpy[grid->m];
	if((mass->c4).n_cnpy < 0.0){
		(mass->c4).n_cnpy = 0.0;
	}
	(mass->c4).n_cnpy_m[grid->m] = (mass->c4).n_cnpy;
	
	/* storage */
	(mass->c3).n_strg += (flux->c3).n_alloc_strg[grid->m] 
						- (flux->c3).n_realloc[grid->m] 
						+ (flux->c3).n_salvage[grid->m] 
						- (flux->c3).n_abdn_strg[grid->m];
	if((mass->c3).n_strg < 0.0){
		(mass->c3).n_strg = 0.0;
	}
	(mass->c3).n_strg_m[grid->m] = (mass->c3).n_strg;
	
	(mass->c4).n_strg += (flux->c4).n_alloc_strg[grid->m] 
						- (flux->c4).n_realloc[grid->m] 
						+ (flux->c4).n_salvage[grid->m] 
						- (flux->c4).n_abdn_strg[grid->m];
	if((mass->c4).n_strg < 0.0){
		(mass->c4).n_strg = 0.0;
	}
	(mass->c4).n_strg_m[grid->m] = (mass->c4).n_strg;
	
	/* soil N *************************************************/
	/* microbe */
	(mass->soil).n_mcrb += (flux->soil).n_immbl_nh4[grid->m]
                        + (flux->soil).n_immbl_no3[grid->m]
						- (flux->soil).n_mcrb_abdn[grid->m];
	if((mass->soil).n_mcrb < 0.0){
		(mass->soil).n_mcrb = 0.0;
	}
	(mass->soil).n_mcrb_m[grid->m] = (mass->soil).n_mcrb;
	
	/* litter */
	(mass->soil).n_lttr += (flux->plant).n_abdn_cnpy[grid->m]
						+ (flux->plant).n_abdn_strg[grid->m]
						- (flux->soil).n_minerlz_lttr[grid->m];
	if((mass->soil).n_lttr < 0.0){
		(mass->soil).n_lttr = 0.0;
	}
	(mass->soil).n_lttr_m[grid->m] = (mass->soil).n_lttr;

	/* humus */
	(mass->soil).n_hums += (flux->soil).n_mcrb_abdn[grid->m]
						- (flux->soil).n_minerlz_hums[grid->m];
	if((mass->soil).n_hums < 0.0){
		(mass->soil).n_hums = 0.0;
	}
	(mass->soil).n_hums_m[grid->m] = (mass->soil).n_hums;

	/* nitrate */
	(mass->soil).n_no3 += loct->depo_no3[grid->m] 
						+ (flux->soil).n_nitrif[grid->m]
                        - (flux->soil).n_immbl_no3[grid->m]
						- (flux->plant).uptake_no3[grid->m]
						- (flux->soil).d_n2o_dnt_ngas[grid->m] * 28.0/44.0
						- (flux->soil).d_n2_ngas[grid->m]
						- (flux->soil).n_leach[grid->m];
	if((mass->soil).n_no3 < 0.0){
		(mass->soil).n_no3 = 0.0;
	}
	(mass->soil).n_no3_m[grid->m] = (mass->soil).n_no3;

	/* ammonium */
	(mass->soil).n_nh4 += loct->depo_nh4[grid->m] 
						+ (flux->soil).n_minerlz_lttr[grid->m]
						+ (flux->soil).n_minerlz_hums[grid->m]
						- (flux->soil).n_immbl_nh4[grid->m]
						- (flux->plant).uptake_nh4[grid->m]
						- (flux->soil).n_nitrif[grid->m]
						- (flux->soil).d_n2o_ntr_ngas[grid->m] * 28.0/44.0
						- (flux->soil).n_nh3vlt[grid->m];
	if((mass->soil).n_nh4 < 0.0){
		(mass->soil).n_nh4 = 0.0;
	}
	(mass->soil).n_nh4_m[grid->m] = (mass->soil).n_nh4;
}

/* biological N2 fixation model by Cleveland **********************************/
/*
Cleveland, C.C., Townsend, A.R., Schimel, D.S., Fisher, H., Howarth, 
R.W., Hedin, L.O., Perakis, S.S., Latty, E.F., Von Fischer, J.C., 
Elseroad, A. and Wasson, M.F., 1999. Global patterns of terrestrial 
biological nitrogen (N2) fixation in natural ecosystems. 
Global Biogeochemical Cycles 13, 623-645.
*/
/**  biological N2 fixation, both symbiotic and non-symbiotic  **/
void f_biolfix(
	struct Grid *grid, 
	struct Loct *loct,
    struct Mass *mass,
	struct Flux *flux
){
	double aet, total_nbiofix, max_n_c3, max_n_c4, n_c3, n_c4;
    
    /* N saturation: 2016/08/03 by A.Ito */
    n_c3 = (mass->c3).n_cnpy + (mass->c3).n_strg;
    n_c4 = (mass->c4).n_cnpy + (mass->c4).n_strg;
	
    max_n_c3 = 1000000.0 * ((mass->c3).fol/12.0 + ((mass->c3).stm + (mass->c3).rot)/24.0);
    max_n_c4 = 1000000.0 * ((mass->c4).fol/12.0 + ((mass->c4).stm + (mass->c4).rot)/24.0);

	/* annual actual evapotranspiration */
	aet = loct->incep[grid->m] + loct->evpr[grid->m] + loct->trspr[grid->m];
	aet *= 0.1; /* cm month-1 */
	
	/* Cleveland et al. (1999) Figure 1, Central */
	total_nbiofix = 0.234 * aet - 0.172; /* kg N ha-1 yr-1 */
	total_nbiofix *= 1000.0; /* g N ha-1 yr-1 */
	if(total_nbiofix < 0.0){
		total_nbiofix = 0.0;
	}
    
    /* lower N2 fixation: 2015/09/03 by A.Ito */
    /* if(loct->v_type==2 && BIOFUEL_RUN >= 1){
        total_nbiofix *= 0.1;
    } */
    
    /* 2016/06/08 by A.Ito */
    if(EX_NITROGEN == 1){
        total_nbiofix *= 0.5;
    }
    
    /* 2016/07/25 by A.Ito  for debugging */
    /* total_nbiofix *= 0.01; */
    
    if(n_c3 > max_n_c3){
        (flux->c3).n_biofix[grid->m] = 0.0;
    }else{
        (flux->c3).n_biofix[grid->m] = total_nbiofix;
    }

    if(n_c4 > max_n_c4){
        (flux->c4).n_biofix[grid->m] = 0.0;
    }else{
        (flux->c4).n_biofix[grid->m] = total_nbiofix;
    }
}

/** NH3 volatilization **********************************************/
void f_nh3_volatilization(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Schar *schar, 
	struct Smas *mass, 
	struct Sflx *flux
){
	double f_ph, f_tmp, f_sw;
	double swp, nh4_soil, ph_soil, base_ph;
	
	/* g N ha-1 */
	/*  nh4_soil = grid->total_n_1m * 0.3 * 0.01*10000.0 *0.5;  */
	/*  nh4_soil = (mass->soil).n_nh4 /1000.0 /10000.0;  */
	nh4_soil = mass->n_nh4;
    
    /* nh4_soil = 10000.0; */
	
	/* pH: Lin et al. (2000) Eq.(16) */
	/* modified by A.Ito (2009/06/05) */
	/* f_ph = pow(10.0, grid->soil_ph - 10.0) / pow(10.0, 7.0 - 10.0); */
	/* if(schar->v_type==1 && (grid->veg_olson==29 || grid->veg_olson==30
							|| grid->veg_olson==31 || grid->veg_olson==32)){ */
    /* revised (after comments by E.Kato): 2013/10/02 by A.Ito */
    /* revised: 2014/11/27 by A.Ito */
	if(schar->v_type == 2){
        /* cropland */
		if((grid->soil_ph + 0.5) >= 6.5){ /* 6.0=>6.5 */
			ph_soil = grid->soil_ph + 1.5;
		}else{
			ph_soil = 6.5;
		}
	}else{
        /* natural */
		/* ph_soil = grid->soil_ph - 0.4; */
		ph_soil = grid->soil_ph;
	}
	
	/* base_ph = 6.5; */ /* 2010/03/28 (A.Ito) */
	/* base_ph = 5.5; */ /* 2010/03/30 (A.Ito) */
	/* base_ph = 5.4; */ /* 2014/12/02 (A.Ito) */
    /* revised: 2016/07/12 by A.Ito */
    base_ph = 6.0;
	f_ph = pow(10.0, ph_soil - 10.0) / pow(10.0, base_ph - 10.0);
	if(f_ph < 0.0){
		f_ph = 0.0;
	}
	
	/* temperature: Thornley (1998) Eq.(3.11a) */
	f_tmp = pow((grid->tmp10_soil[grid->m] - 0.0), 2.0)*(45.0 - grid->tmp10_soil[grid->m]) 
			/ (pow((20.0 - 0.0), 2.0)*(45.0 - 20.0));
	if(f_tmp < 0.0){
		f_tmp = 0.0;
	}
	
	/* soil water: Thornley (1998) Eq.(6.7a) */
	/* soil water potential Eq.(6.2g) */
    /* modified: 2016/07/03 by A.Ito */
	if(loct->sw30 > 1.0){
		/* modified by A.Ito (2009/06/05) */
		swp = -10.0 * pow(1.0 / (loct->sw30 / grid->field_cap1), 5.0);
		f_sw = exp((18.0 * swp) / (8314.0 * (grid->tmp10_soil[grid->m] + ZAT)));
        if(f_sw > 10.0){
            f_sw = 10.0;
        }
	}else{
		f_sw = 0.05; /* for dry land sublimation */
	}
    
    /* modified: 2016/07/11 by A.Ito */
    /* Thornley (1998) Eq.(5.11) */
    f_sw = pow(f_sw, 20.0);
	if(f_sw < 0.01){
		f_sw = 0.01;
	}
	
	/* Thornley (1998) Eq.(5.4i) */
	/* g NH3 ha-1 month-1 */
	flux->n_nh3vlt[grid->m] = nh4_soil * 0.02/30.0 * f_ph * f_tmp * f_sw *
			MDN[grid->m] * 17.0/14.0;
    
    /* safe guard: 2014/05/28 by A.Ito: 0.5 */
    /* safe guard: 2016/06/28 by A.Ito: 0.1 */
    /* safe guard: 2016/07/03 by A.Ito: 0.1 */
    /* safe guard: 2016/07/07 by A.Ito: 0.75 */
    if(flux->n_nh3vlt[grid->m] > (0.75 * mass->n_nh4)){
        flux->n_nh3vlt[grid->m] = 0.75 * mass->n_nh4;
    }
}

/* N deposition ********************************************************/
/* 
Dentener, F. J. 2006. Global Maps of Atmospheric Nitrogen Deposition, 
1860, 1993, and 2050. Data set. Available on-line from Oak Ridge National 
Laboratory Distributed Active Archive Center, Oak Ridge, Tennessee, U.S.A. 
*/
void f_n_deposit(
	struct Grid *grid, 
	struct Loct *loct
){
    long f, nyear;
	double pre_ann, ndepo_total, ndepo_dry, ndepo_wet, aa;
	double f_no3, f_nh4, ndepo_no3, ndepo_nh4;
	double f_wet, f_dry, uconv;
	extern double MDN[ASTEP];
	
	/* f_no3 = 0.75; */
	/*f_no3 = 0.47; */ /* revised by CHASER data: 2010/03/28 (A.Ito) */
	/* f_nh4 = 1.0 - f_no3; */
	
    if(EX_CHASER_NDEPO == 0 || EX_CHASER_NDEPO == 1){
    
        if(EX_CHASER_NDEPO == 0){
            /* CHASER-derived spatial and monthly NH4+/NO3- fraction */
            ndepo_total = ndepo_chaser_dnhx[grid->m][grid->chaser_row][grid->chaser_col]
                            + ndepo_chaser_dnoy[grid->m][grid->chaser_row][grid->chaser_col]
                            + ndepo_chaser_wnhx[grid->m][grid->chaser_row][grid->chaser_col]
                            + ndepo_chaser_wnoy[grid->m][grid->chaser_row][grid->chaser_col];
            
            if(ndepo_total <= 0.0){
                f_no3 = 0.5;
                f_nh4 = 0.5;
            }else{
                f_no3 = (ndepo_chaser_dnoy[grid->m][grid->chaser_row][grid->chaser_col]
                    + ndepo_chaser_wnoy[grid->m][grid->chaser_row][grid->chaser_col]) / ndepo_total;
                f_nh4 = (ndepo_chaser_dnhx[grid->m][grid->chaser_row][grid->chaser_col]
                    + ndepo_chaser_wnhx[grid->m][grid->chaser_row][grid->chaser_col]) / ndepo_total;
            }
            
            if(ndepo_total <= 0.0){
                f_wet = 0.5;
                f_dry = 0.5;
            }else{
                f_wet = (ndepo_chaser_wnhx[grid->m][grid->chaser_row][grid->chaser_col]
                    + ndepo_chaser_wnoy[grid->m][grid->chaser_row][grid->chaser_col]) / ndepo_total;
                f_dry = (ndepo_chaser_dnoy[grid->m][grid->chaser_row][grid->chaser_col]
                    + ndepo_chaser_dnhx[grid->m][grid->chaser_row][grid->chaser_col]) / ndepo_total;
            }
        }else if(EX_CHASER_NDEPO == 1){
            f_dry = 0.5;
            f_wet = 0.5;
        }
        
        /* annual precipitation */
        /* 1.0 : initial value: avoid zero */
        pre_ann = grid->prate_sfc_ann + 1.0;

        if(grid->climy<=1850){
            ndepo_dry = f_dry * grid->ndepo[0] * MDN[grid->m]/365.0;
            ndepo_wet = f_wet * grid->ndepo[0] * (grid->prate_sfc_a[grid->m] + 0.08333)/pre_ann;
        }else if(grid->climy>1850 && grid->climy<=1993){
            aa = grid->ndepo[0] + (grid->ndepo[1] - grid->ndepo[0])*(double)(grid->climy-1850)/143.0;
        
            ndepo_dry = f_dry * aa * MDN[grid->m]/365.0;
            ndepo_wet = f_wet * aa * (grid->prate_sfc_a[grid->m] + 0.08333)/pre_ann;
        }else if(grid->climy>1993 && grid->climy<=2050){
            aa = grid->ndepo[1] + (grid->ndepo[2] - grid->ndepo[1])*(double)(grid->climy-1993)/57.0;
        
            ndepo_dry = f_dry * aa * MDN[grid->m]/365.0;
            ndepo_wet = f_wet * aa * (grid->prate_sfc_a[grid->m] + 0.08333)/pre_ann;
        }else{ /*  if(grid->climy>2050) */
            ndepo_dry = grid->ndepo[2] * f_dry * MDN[grid->m]/365.0;
            ndepo_wet = grid->ndepo[2] * f_wet * (grid->prate_sfc_a[grid->m] + 0.08333)/pre_ann;
            /* 2008/08/20 corrected by A.Ito (thanks to E.Kato) */
        }
	
        if(ndepo_dry < 0.0){
            ndepo_dry = 0.0;
        }
        if(ndepo_wet < 0.0){
            ndepo_wet = 0.0;
        }

        /* original unit: mg N m-2 yr-1 */
        /* converted unit: g N ha-1 month-1 */
        loct->depo_no3[grid->m] = f_no3 * (ndepo_dry + ndepo_wet) *10.0;
        loct->depo_nh4[grid->m] = f_nh4 * (ndepo_dry + ndepo_wet) *10.0;
	}else if(EX_CHASER_NDEPO == 2){
        
        /* 2014/12/27 revised by A.Ito: add organic N deposition*/
        if(grid->climy<1850){
            ndepo_no3 = ndepo_chaser4_noy_h[grid->m][grid->chaser_row][grid->chaser_col]
                    + ndepo_chaser4_ont_h[grid->m][grid->chaser_row][grid->chaser_col];
            ndepo_nh4 = ndepo_chaser4_nhx_h[grid->m][grid->chaser_row][grid->chaser_col];
        }else if(grid->climy>=1850 && grid->climy<=2010){
            
            ndepo_no3 = ndepo_chaser4_noy_h[grid->m][grid->chaser_row][grid->chaser_col]
                + (ndepo_chaser4_noy_p[grid->m][grid->chaser_row][grid->chaser_col]
                        - ndepo_chaser4_noy_h[grid->m][grid->chaser_row][grid->chaser_col]) *
                    ((double)(grid->climy) - 1850.0)/160.0;
            ndepo_no3 += ndepo_chaser4_ont_h[grid->m][grid->chaser_row][grid->chaser_col]
                + (ndepo_chaser4_ont_p[grid->m][grid->chaser_row][grid->chaser_col]
                        - ndepo_chaser4_ont_h[grid->m][grid->chaser_row][grid->chaser_col]) *
                    ((double)(grid->climy) - 1850.0)/160.0;
            
            ndepo_nh4 = ndepo_chaser4_nhx_h[grid->m][grid->chaser_row][grid->chaser_col]
                + (ndepo_chaser4_nhx_p[grid->m][grid->chaser_row][grid->chaser_col]
                        - ndepo_chaser4_nhx_h[grid->m][grid->chaser_row][grid->chaser_col]) *
                    ((double)(grid->climy) - 1850.0)/160.0;
            
        }else if(grid->climy>2010){
            ndepo_no3 = ndepo_chaser4_noy_p[grid->m][grid->chaser_row][grid->chaser_col]
                + ndepo_chaser4_ont_p[grid->m][grid->chaser_row][grid->chaser_col];
            ndepo_nh4 = ndepo_chaser4_nhx_p[grid->m][grid->chaser_row][grid->chaser_col];
        }
        
        /* unit: g N ha-1 month-1 */
        loct->depo_no3[grid->m] = ndepo_no3;
        loct->depo_nh4[grid->m] = ndepo_nh4;
    }
    
    /* NMIP: 2015/11/19 by A.Ito *******/
    /* updated 2016/10/20 by A.Ito */
    /* ISI-MIP2b: 2016/12/24 by A.Ito */
    if(NMIP_RUN >= 1 || ISIMIP_RUN == 4){
    
        nyear = grid->niny;
        if(NMIP_RUN == 4 || NMIP_RUN == 5 || NMIP_RUN == 6){
            nyear = FDY_NINY+1;
        }
        
        /* to gN/ha */
        if(NMIP_RUN >= 1){
            uconv = 1000.0;
        }
        if(ISIMIP_RUN == 4){
            uconv = 10000.0;
        }

        /* seasonality based on CHASER */
        ndepo_no3 = ndepo_nh4 = 0.0;
        for(f=0;f<ASTEP;f++){
            ndepo_no3 += ndepo_chaser4_noy_h[f][grid->chaser_row][grid->chaser_col]
                    + ndepo_chaser4_ont_h[f][grid->chaser_row][grid->chaser_col];
            ndepo_nh4 += ndepo_chaser4_nhx_h[f][grid->chaser_row][grid->chaser_col];
        }
        
        if(ndepo_no3 > 0.0){
            loct->depo_no3[grid->m] = grid->nmip_ndep_noy[nyear - FDY_NINY]*uconv *
                    (ndepo_chaser4_noy_h[grid->m][grid->chaser_row][grid->chaser_col]
                    + ndepo_chaser4_ont_h[grid->m][grid->chaser_row][grid->chaser_col]) / ndepo_no3;
        }else{
            loct->depo_no3[grid->m] = grid->nmip_ndep_noy[grid->niny-FDY_NINY]*uconv / 12.0;
        }
        if(ndepo_nh4 > 0.0){
            loct->depo_nh4[grid->m] =  grid->nmip_ndep_nh4[nyear - FDY_NINY]*uconv *
                ndepo_chaser4_nhx_h[grid->m][grid->chaser_row][grid->chaser_col] / ndepo_nh4;
        }else{
            loct->depo_nh4[grid->m] = grid->nmip_ndep_nh4[grid->niny - FDY_NINY]*uconv / 12.0;
        }
    }
    
    /* experiment N deposition: 2015/08/12 by A.Ito *****/
    if(EX_NDEPO == 1){
        loct->depo_no3[grid->m] *= 1.1;
        loct->depo_nh4[grid->m] *= 1.1;
    }
    if(EX_NDEPO == 2){
        loct->depo_no3[grid->m] *= 0.9;
        loct->depo_nh4[grid->m] *= 0.9;
    }
    if(EX_NDEPO == 3){
        loct->depo_no3[grid->m] *= 1.3;
        loct->depo_nh4[grid->m] *= 1.3;
    }
    if(EX_NDEPO == 4){
        loct->depo_no3[grid->m] *= 0.7;
        loct->depo_nh4[grid->m] *= 0.7;
    }
    if(EX_NDEPO == 5){
        loct->depo_no3[grid->m] *= 1.5;
        loct->depo_nh4[grid->m] *= 1.5;
    }
    if(EX_NDEPO == 6){
        loct->depo_no3[grid->m] *= 0.5;
        loct->depo_nh4[grid->m] *= 0.5;
    }
    
    if(EX_NDEPO == 7){
        ndepo_total = loct->depo_no3[grid->m] + loct->depo_nh4[grid->m];
        loct->depo_nh4[grid->m] = ndepo_total;
        loct->depo_no3[grid->m] = 0.0;
    }
    if(EX_NDEPO == 8){
        ndepo_total = loct->depo_no3[grid->m] + loct->depo_nh4[grid->m];
        loct->depo_nh4[grid->m] = 0.0;
        loct->depo_no3[grid->m] = ndepo_total;
    }
    if(EX_NDEPO == 9){
        ndepo_total = loct->depo_no3[grid->m] + loct->depo_nh4[grid->m];
        loct->depo_nh4[grid->m] = 0.5 * ndepo_total;
        loct->depo_no3[grid->m] = 0.5 * ndepo_total;
    }
    
    /* 2016/07/25 by A.Ito  for debugging */
    /* loct->depo_no3[grid->m] *= 0.01;
    loct->depo_nh4[grid->m] *= 0.01; */
}

/*************************************************/
void f_n_mineralz(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Smas *mass, 
	struct Sflx *flux
){
	double f_c_min;
    double f_nmin_l, f_nmin_h;
    
    /* 2016/06/28 by A.Ito */
    /* f_nmin_l = 2.0;
    f_nmin_h = 20.0; */
	
    /* 2016/07/06 by A.Ito */
    /* f_nmin_l = 50.0;
    f_nmin_h = 80.0; */
	
    /* 2016/08/15 by A.Ito */
    f_nmin_l = 0.15;
    f_nmin_h = 0.7;
	
	/** litter **/
	if(mass->ltr > 0.01){
		f_c_min = flux->rl[grid->m] / mass->ltr;
	}else{
		f_c_min = 0.0;
	}
	
	/* g N ha-1 month-1 */
	flux->n_minerlz_lttr[grid->m] = f_nmin_l * mass->n_lttr * f_c_min;

	/** humus **/
	if(mass->msl > 0.01){
		f_c_min = flux->rh[grid->m] / mass->msl;
	}else{
		f_c_min = 0.0;
	}
	
	/* g N ha-1 month-1 */
	flux->n_minerlz_hums[grid->m] = f_nmin_h * mass->n_hums * f_c_min;
}

/* leaching of NO3 **************************************/
void f_n_leaching(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Smas *mass, 
	struct Sflx *flux
){
	double ntr_conc, aa;
	double fad_no3;
	
	/* adsorption of NO3- : 2010/03/28 by A.Ito */
	/* fad_no3 = 0.7; */
	/* fad_no3 = 0.3; */ /* 2010/03/29 by A.Ito */
	/* fad_no3 = 0.1; */ /* 2010/04/06 by A.Ito */
	/* fad_no3 = 0.5; */ /* 2016/06/01 by A.Ito */
	/* fad_no3 = 0.3; */ /* 2016/06/29 by A.Ito */
	fad_no3 = 0.70; /* 2016/08/15 by A.Ito */

	/* g N / ha */
	/* kg H2O / m2 */ /* 2016/07/05 by A.Ito */
	if((loct->sw30 + loct->sww + grid->prate_sfc[grid->m]) > 0.1){
		ntr_conc = (1.0 - fad_no3)*(mass->n_no3/10000.0) / 
					(loct->sw30 + loct->sww + grid->prate_sfc[grid->m]);
	}else{
		ntr_conc = 0.0;
	}
	/* g N / kg H2O */
	
	/* g NO3-N m-2 */
	aa = loct->ro2[grid->m] * ntr_conc * 10000.0;
	/* g NO3-N ha-1 */
    
    /* 2016/07/13 by A.Ito */
	/* if(aa > (mass->n_no3 * 0.95)){ */
	/* if(aa > (mass->n_no3 * 0.9)){ */
	if(aa > (mass->n_no3 * 0.75)){
		aa = mass->n_no3 * 0.75;
	}
    if(aa < 0.0){
        aa = 0.0;
    }
	
    /* 2016/07/07 by A.Ito */
    /* aa = mass->n_no3*0.1; */
    /* aa = 1.0; */
    
	/* g NO3-N ha-1 month-1 */
	flux->n_leach[grid->m] = aa;
}

/* plant N uptake *******************************************/
void f_n_uptake(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Flux *flux
){
	double f_temp, max_uptake, n_c3, n_c4;
	double n_max = 1.0;
	double ks, navil;
    double max_n_c3, max_n_c4, nsat_c3, nsat_c4, k_n;
	double uptake_no3, uptake_nh4;
    
    /* 2016/07/21 by A.Ito ************/
    /* N demand */
    
    n_c3 = (mass->c3).n_cnpy + (mass->c3).n_strg;
    n_c4 = (mass->c4).n_cnpy + (mass->c4).n_strg;
    
    /* g N ha-1 */
    max_n_c3 = 1000000.0 * ((mass->c3).fol/12.0 + ((mass->c3).stm + (mass->c3).rot)/24.0);
    max_n_c4 = 1000000.0 * ((mass->c4).fol/12.0 + ((mass->c4).stm + (mass->c4).rot)/24.0);
    
    k_n = 1.0 + max_n_c3 * 0.1;
    nsat_c3 = n_c3 / (k_n + n_c3);
    
    k_n = 1.0 + max_n_c4 * 0.1;
    nsat_c4 = n_c4 / (k_n + n_c4);

	/* Raich et al. (1991): Appendix  Eq.1.16 */
	
	/*****
	Effect of N allocatiom to root ?
	*****/
    /* n_max = 0.02; */ /* 2016/07/06 by A.Ito */
	
	/* temperature factor */
	f_temp = exp(0.0693 * grid->tmp10_soil[grid->m]);
	/* soil diffusion parameter */
	ks = 0.90 * pow(loct->sw30 / grid->field_cap1, 3.0) + 0.1;
	
	/* NO3 uptake */
    n_max = 0.14; /* 2016/08/15 by A.Ito */
	navil = (mass->soil).n_no3;
	/* C3 */
    max_uptake = (1.0 - nsat_c3) * navil * n_max * ks / (90.0 + ks*navil) * f_temp;
    if(max_uptake>=0.0 && max_uptake<navil){
        uptake_no3 = max_uptake;
    }else if(max_uptake < 0.0){
        uptake_no3 = 0.0;
    }else if(max_uptake >= navil){
        uptake_no3 = navil* n_max;
    }
	/* g N ha-1 month-1 */
	(flux->c3).uptake_no3[grid->m] = uptake_no3;
    
	/* C4 */
    max_uptake = (1.0 - nsat_c4) * navil * n_max * ks / (90.0 + ks*navil) * f_temp;
    if(max_uptake>=0.0 && max_uptake<navil){
        uptake_no3 = max_uptake;
    }else if(max_uptake < 0.0){
        uptake_no3 = 0.0;
    }else if(max_uptake >= navil){
        uptake_no3 = navil* n_max;
    }
	/* g N ha-1 month-1 */
	(flux->c4).uptake_no3[grid->m] = uptake_no3;
    
	
	/* NH4 uptake */
    n_max = 0.07; /* 2016/08/15 by A.Ito */
	navil = (mass->soil).n_nh4;
    /* C3 */
    max_uptake = (1.0 - nsat_c3) * navil * n_max * ks / (90.0 + ks*navil) * f_temp;
    if(max_uptake>0.0 && max_uptake<navil){
        uptake_nh4 = max_uptake;
    }else if(max_uptake < 0.0){
        uptake_nh4 = 0.0;
    }else if(max_uptake >= navil){
        uptake_nh4 = navil* n_max;
    }
	/* g N ha-1 month-1 */
	(flux->c3).uptake_nh4[grid->m] = uptake_nh4;
    
    /* C4 */
    max_uptake = (1.0 - nsat_c4) * navil * n_max * ks / (90.0 + ks*navil) * f_temp;
    if(max_uptake>0.0 && max_uptake<navil){
        uptake_nh4 = max_uptake;
    }else if(max_uptake < 0.0){
        uptake_nh4 = 0.0;
    }else if(max_uptake >= navil){
        uptake_nh4 = navil* n_max;
    }
	/* g N ha-1 month-1 */
	(flux->c4).uptake_nh4[grid->m] = uptake_nh4;
    	
	(flux->plant).uptake_no3[grid->m] = loct->c3ptn[grid->m] * (flux->c3).uptake_no3[grid->m]
				+ loct->c4ptn[grid->m] * (flux->c4).uptake_no3[grid->m];
	(flux->plant).uptake_nh4[grid->m] = loct->c3ptn[grid->m] * (flux->c3).uptake_nh4[grid->m]
				+ loct->c4ptn[grid->m] * (flux->c4).uptake_nh4[grid->m];
}

/* N abandoned as litter *********************************/
void f_n_abandon_salvage(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){
	/* canopy ****/
	/* abandon */
	flux->n_abdn_cnpy[grid->m] = (1.0 - pchar->n_salvage) * pchar->lf[grid->m] * 
					mass->n_cnpy;
	/* salvage (draw back to storage) */
	flux->n_salvage[grid->m] = pchar->n_salvage * pchar->lf[grid->m] * 
					mass->n_cnpy;

	/* storage ****/
	flux->n_abdn_strg[grid->m] = (flux->lc[grid->m] + flux->lr[grid->m]) *
				1000000.0 / pchar->cn0_strg;
	if(flux->n_abdn_strg[grid->m] > mass->n_strg*0.9){
		flux->n_abdn_strg[grid->m] = mass->n_strg*0.9;
	}
}

/* allocation of uptake N **********************************/
void f_n_alloc(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){
	double n_obtain;
	double n_leaf_conc;
	double amax;	/* maximum photosynthetic rate at high N */
	double kn;		/* half-saturation constant */
	double rd_nsp;	/* specific dark respiration rate for N */
	double n_opt;	/* optimal N for nitrogen utilization */
	double a_nlmt;
	double n_demand1, n_demand2;
	
	/* nitrogen obtained by root uptake and symbiotic fixation */
	n_obtain = flux->n_biofix[grid->m] + flux->uptake_no3[grid->m] + flux->uptake_nh4[grid->m];

	amax = pchar->amax_nphoto;
	kn = pchar->kn_nphoto;
	rd_nsp = pchar->rd_n;
	/* optimal leaf N */
	n_opt = sqrt(amax * kn/rd_nsp) - kn;
	
	/* leaf N concentration, mmol N m-2 */
	n_leaf_conc = pchar->n_conc_larea;
	a_nlmt = amax * n_leaf_conc / (kn + n_leaf_conc);
	
	if(n_opt > n_leaf_conc){
		n_demand1 = (n_opt - n_leaf_conc) * 14.0 / 1000.0 * 10000.0;
		
		if(n_obtain > n_demand1){
			flux->n_alloc_cnpy[grid->m] = n_demand1;
		}else{
			flux->n_alloc_cnpy[grid->m] = n_obtain;
		}
	}else{  /*  if(n_opt <= n_leaf_conc) */
        /* 2016/08/02 by A.Ito */
        if(n_leaf_conc < 1.1*n_opt){
            flux->n_realloc[grid->m] = 0.0;
        }else{
            n_demand2 = -0.1 * (n_leaf_conc - n_opt);
            flux->n_realloc[grid->m] = n_demand2;
        }
	}
	
	flux->n_alloc_strg[grid->m] = n_obtain - flux->n_alloc_cnpy[grid->m];
}

/* N reallocation from storage to canopy **********/
void f_n_realloc(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){
	double amax;	/* maximum photosynthetic rate at high N */
	double kn;		/* half-saturation constant */
	double rd_nsp;	/* specific dark respiration rate for N */
	double n_opt;	/* optimal N for nitrogen utilization */
	double n_leaf_conc;
	double n_demand1, n_demand2;
	double n_stock;
	double a_nlmt;
	
	if(mass->n_strg>0.0 && mass->n_cnpy>0.0){
		/* nitrogen stock for N reallocation */
		n_stock = mass->n_strg;
	
		amax = pchar->amax_nphoto;
		kn = pchar->kn_nphoto;
		rd_nsp = pchar->rd_n;
		/* optimal leaf N */
		n_opt = sqrt(amax * kn / rd_nsp) - kn;
		/* n_opt = 20.0; */
        
		/* leaf N concentration, mmol N m-2 */
		n_leaf_conc = pchar->n_conc_larea;
		a_nlmt = amax * n_leaf_conc / (kn + n_leaf_conc);
		/* a_nlmt = n_leaf_conc; */
		
        pchar->n_opt[grid->m] = n_opt;
        pchar->n_leaf[grid->m] = n_leaf_conc;

		if(n_opt > n_leaf_conc){
			/* g N ha-1 month-1 */
			n_demand1 = (n_opt - n_leaf_conc) * 14.0 / 1000.0 * 10000.0;
			
			if(n_stock > n_demand1){
				flux->n_realloc[grid->m] = n_demand1;
			}else{
				flux->n_realloc[grid->m] = n_stock;
			}
		}else{  /* if(n_opt <= n_leaf_conc) */
            /* 2016/08/02 by A.Ito */
            if(n_leaf_conc < 1.1*n_opt){
                flux->n_realloc[grid->m] = 0.0;
            }else{
                n_demand2 = -0.1 * (n_leaf_conc - n_opt);
                
                flux->n_realloc[grid->m] = n_demand2;
            }
		}
	}else{
		flux->n_realloc[grid->m] = 0.0;
	}
}

/* N immobilization by microbes *****************************/
void f_n_immoblz(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Schar *schar, 
	struct Smas *mass, 
	struct Sflx *flux
){
	extern double MDN[12];
	double f_immbl_no3, f_immbl_nh4;
	
	/* f_immbl_no3 = 0.002;
	f_immbl_nh4 = 0.001;
	f_immbl_no3 = 0.005;
	f_immbl_nh4 = 0.004; */
    
    /* 2016/06/01 by A.Ito *****/
	/* f_immbl_no3 = 0.006;
	f_immbl_nh4 = 0.003; */
    
    /* 2016/06/29 by A.Ito *****/
	/* f_immbl_no3 = 0.004;
	f_immbl_nh4 = 0.002; */
    
    /* 2016/07/06 by A.Ito *****/
    /* 2016/08/04 by A.Ito *****/
	f_immbl_no3 = 0.00025;
	f_immbl_nh4 = 0.00025;
    
    /* 2016/06/05 by A.Ito *****/
	/* flux->n_immbl[grid->m] = 0.05 * flux->n_minerlz_lttr[grid->m] +
		0.1 * flux->n_minerlz_hums[grid->m] +
		(f_immbl_no3 * mass->n_no3 + f_immbl_nh4 * mass->n_nh4) * MDN[grid->m]; */
    
    /* 2016/07/03 by A.Ito *****/
    /* 2016/07/06 by A.Ito *****/
	/* flux->n_immbl[grid->m] = 0.025 * flux->n_minerlz_lttr[grid->m] +
		0.025 * flux->n_minerlz_hums[grid->m] +
		(f_immbl_no3 * mass->n_no3 + f_immbl_nh4 * mass->n_nh4) * MDN[grid->m]; */

	flux->n_immbl_no3[grid->m] = f_immbl_no3 * mass->n_no3 * MDN[grid->m];
	flux->n_immbl_nh4[grid->m] = f_immbl_nh4 * mass->n_nh4 * MDN[grid->m];

	/* flux->n_immbl[grid->m] = 0.2 * flux->n_minerlz_lttr[grid->m] +
		0.4 * flux->n_minerlz_hums[grid->m] + 
		(f_immbl_no3 * mass->n_no3 + f_immbl_nh4 * mass->n_nh4) * MDN[grid->m]; */
    
    /* safe guard: 2014/05/28 by A.Ito */
    /* 2016/06/05 by A.Ito *****/
    /* 2016/07/12 by A.Ito */
    /* if((flux->n_immbl_no3[grid->m]+flux->n_immbl_nh4[grid->m]) > 0.1*mass->n_mcrb){
        flux->n_immbl_no3[grid->m] *= 0.1*mass->n_mcrb / (flux->n_immbl_no3[grid->m]+flux->n_immbl_nh4[grid->m]);
        flux->n_immbl_nh4[grid->m] *= 0.1*mass->n_mcrb / (flux->n_immbl_no3[grid->m]+flux->n_immbl_nh4[grid->m]);
    } */
}

/* N abandoned from microbes ******************************/
void f_n_mcrb_abdn(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Schar *schar, 
	struct Smas *mass, 
	struct Sflx *flux
){
	double f_temp;
	
	f_temp = exp(log(2.0)/10.0 * (grid->tmp10_soil[grid->m]-10.0));
    
    /* 2016/06/28 by A.Ito */
	/* flux->n_mcrb_abdn[grid->m] = 0.1 * f_temp * mass->n_mcrb; */
	/* flux->n_mcrb_abdn[grid->m] = 0.4 * f_temp * mass->n_mcrb; */
	/* flux->n_mcrb_abdn[grid->m] = 3.0 * f_temp * mass->n_mcrb; */
    /* 2016/08/04 by A.Ito */
	flux->n_mcrb_abdn[grid->m] = 0.025 * f_temp * mass->n_mcrb;
    
    /* 2016/06/08 by A.Ito */
    if(EX_NITROGEN == 2){
        flux->n_mcrb_abdn[grid->m] *= 2.0;
    }

    /* safe guard: 2014/05/28 by A.Ito */
    /* 2016/06/05 by A.Ito *****/
    if(flux->n_mcrb_abdn[grid->m] > (0.7 * mass->n_mcrb)){
        /* flux->n_mcrb_abdn[grid->m] = (0.5 * mass->n_mcrb); */
        flux->n_mcrb_abdn[grid->m] = (0.7 * mass->n_mcrb);
    }
}
