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
	(mass->soil).n_mcrb += (flux->soil).n_immbl[grid->m]
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
						- (flux->plant).uptake_no3[grid->m]
						- (flux->soil).d_n2o_dnt_ngas[grid->m]*28.0/44.0
						- (flux->soil).d_n2_ngas[grid->m]
						- (flux->soil).n_leach[grid->m];
	if((mass->soil).n_no3 < 0.0){
		(mass->soil).n_no3 = 0.0;
	}
	(mass->soil).n_no3_m[grid->m] = (mass->soil).n_no3;

	/* ammonium */
	(mass->soil).n_nh4 += loct->depo_nh4[grid->m] 
						+ ((flux->soil).n_minerlz_lttr[grid->m] 
							+(flux->soil).n_minerlz_hums[grid->m]
							- (flux->soil).n_immbl[grid->m])
						- (flux->plant).uptake_nh4[grid->m]
						- (flux->soil).n_nitrif[grid->m]
						- (flux->soil).d_n2o_ntr_ngas[grid->m]*28.0/44.0
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
	struct Flux *flux
){
	double aet, total_nbiofix;
	
	/* annual actual evapotranspiration */
	aet = loct->incep[grid->m] + loct->evpr[grid->m] + loct->trspr[grid->m];
	aet *= 0.1; /* cm month-1 */
	
	/* Cleveland et al. (1999) Figure 1, Central */
	total_nbiofix = 0.234 * aet - 0.172; /* kg N ha-1 yr-1 */
	total_nbiofix *= 1000.0; /* g N ha-1 yr-1 */
	if(total_nbiofix < 0.0){
		total_nbiofix = 0.0;
	}
	
	(flux->c3).n_biofix[grid->m] = total_nbiofix;
	(flux->c4).n_biofix[grid->m] = total_nbiofix;
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
	
	/* pH: Lin et al. (2000) Eq.(16) */
	/* modified by A.Ito (2009/06/05) */
	/* f_ph = pow(10.0, grid->soil_ph - 10.0) / pow(10.0, 7.0 - 10.0); */
	/* if(schar->v_type==1 && (grid->veg_olson==29 || grid->veg_olson==30
							|| grid->veg_olson==31 || grid->veg_olson==32)){ */
    /* revised (after comments by E.Kato): 2013/10/02 by A.Ito */
    /* revised: 2014/11/27 by A.Ito */
	if(schar->v_type == 2){
        /* cropland */
		if((grid->soil_ph+0.5) >= 6.5){ /* 6.0=>6.5 */
			ph_soil = grid->soil_ph + 1.5;
		}else{
			ph_soil = 6.5;
		}
	}else{
        /* natural */
		ph_soil = grid->soil_ph - 0.4;
	}
	
	/* base_ph = 6.5; */ /* 2010/03/28 (A.Ito) */
	/* base_ph = 5.5; */ /* 2010/03/30 (A.Ito) */
	base_ph = 5.4; /* 2014/12/02 (A.Ito) */
	f_ph = pow(10.0, ph_soil - 10.0) / pow(10.0, base_ph - 10.0);
	if(f_ph < 0.0){
		f_ph = 0.0;
	}
	
	/* temperature: Thornley (1998) Eq.(3.11a) */
	f_tmp = pow((grid->tmp10_soil[grid->m] - 0.0), 2.0)*(45.0 - grid->tmp10_soil[grid->m]) 
			/ (pow((20.0 - 0.0), 2.0)*(45.0 - 20.0));
	if(f_tmp<0.0){
		f_tmp = 0.0;
	}
	
	/* soil water: Thornley (1998) Eq.(6.7a) */
	/* soil water potential Eq.(6.2g) */
	if(loct->sw30 > 1.0){
		/* modified by A.Ito (2009/06/05) */
		swp = -10.0 * pow(1.0 / (loct->sw30 / grid->field_cap1), 5.0);
		f_sw = exp((18.0 * swp) / (8314.0 * (grid->tmp10_soil[grid->m] + ZAT)));
	}else{
		f_sw = 0.0;
	}
	if(f_sw<0.0){
		f_sw = 0.0;
	}
	
	/* Thornley (1998) Eq.(5.4i) */
	/* g NH3 ha-1 month-1 */
	flux->n_nh3vlt[grid->m] = nh4_soil * 0.02/30.0 * f_ph * f_tmp * pow(f_sw, 20.0) 
			* MDN[grid->m] * 17.0/14.0;
    
    /* safe guard: 2014/05/28 by A.Ito */
    if(flux->n_nh3vlt[grid->m] > (0.5*mass->n_nh4)){
        flux->n_nh3vlt[grid->m] = 0.5 * mass->n_nh4;
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
	double pre_ann, ndepo_total, ndepo_dry, ndepo_wet, aa;
	double f_no3, f_nh4, ndepo_no3, ndepo_nh4;
	double f_wet, f_dry;
	extern double MDN[12];
	
	/* f_no3 = 0.75; */
	/*f_no3 = 0.47; */ /* revised by CHASER data: 2010/03/28 (A.Ito) */
	/* f_nh4 = 1.0 - f_no3; */
	
    if(SENS_N == 0 || SENS_N == 1){
    
        if(SENS_N == 0){
            /* CHASER-derived spatial and monthly NH4+/NO3- fraction */
            ndepo_total = grid->ndepo_chaser_dnhx[grid->m][grid->chaser_row][grid->chaser_col]
                            + grid->ndepo_chaser_dnoy[grid->m][grid->chaser_row][grid->chaser_col]
                            + grid->ndepo_chaser_wnhx[grid->m][grid->chaser_row][grid->chaser_col]
                            + grid->ndepo_chaser_wnoy[grid->m][grid->chaser_row][grid->chaser_col];
            
            if(ndepo_total <= 0.0){
                f_no3 = 0.5;
                f_nh4 = 0.5;
            }else{
                f_no3 = (grid->ndepo_chaser_dnoy[grid->m][grid->chaser_row][grid->chaser_col]
                    + grid->ndepo_chaser_wnoy[grid->m][grid->chaser_row][grid->chaser_col]) / ndepo_total;
                f_nh4 = (grid->ndepo_chaser_dnhx[grid->m][grid->chaser_row][grid->chaser_col]
                    + grid->ndepo_chaser_wnhx[grid->m][grid->chaser_row][grid->chaser_col]) / ndepo_total;
            }
            
            if(ndepo_total <= 0.0){
                f_wet = 0.5;
                f_dry = 0.5;
            }else{
                f_wet = (grid->ndepo_chaser_wnhx[grid->m][grid->chaser_row][grid->chaser_col]
                    + grid->ndepo_chaser_wnoy[grid->m][grid->chaser_row][grid->chaser_col]) / ndepo_total;
                f_dry = (grid->ndepo_chaser_dnoy[grid->m][grid->chaser_row][grid->chaser_col]
                    + grid->ndepo_chaser_dnhx[grid->m][grid->chaser_row][grid->chaser_col]) / ndepo_total;
            }
        }else if(SENS_N == 1){
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
	}else if(SENS_N == 2){
        
        /* 2014/12/27 revised by A.Ito: add organic N deposition*/
        if(grid->climy<1850){
            ndepo_no3 = grid->ndepo_chaser4_noy_h[grid->m][grid->chaser_row][grid->chaser_col]
                    + grid->ndepo_chaser4_ont_h[grid->m][grid->chaser_row][grid->chaser_col];
            ndepo_nh4 = grid->ndepo_chaser4_nhx_h[grid->m][grid->chaser_row][grid->chaser_col];
        }else if(grid->climy>=1850 && grid->climy<=2010){
            
            ndepo_no3 = grid->ndepo_chaser4_noy_h[grid->m][grid->chaser_row][grid->chaser_col]
                + (grid->ndepo_chaser4_noy_p[grid->m][grid->chaser_row][grid->chaser_col]
                        - grid->ndepo_chaser4_noy_h[grid->m][grid->chaser_row][grid->chaser_col])
                * ((double)(grid->climy) - 1850.0)/160.0;
            ndepo_no3 += grid->ndepo_chaser4_ont_h[grid->m][grid->chaser_row][grid->chaser_col]
                + (grid->ndepo_chaser4_ont_p[grid->m][grid->chaser_row][grid->chaser_col]
                        - grid->ndepo_chaser4_ont_h[grid->m][grid->chaser_row][grid->chaser_col])
                * ((double)(grid->climy) - 1850.0)/160.0;
            
            ndepo_nh4 = grid->ndepo_chaser4_nhx_h[grid->m][grid->chaser_row][grid->chaser_col]
                + (grid->ndepo_chaser4_nhx_p[grid->m][grid->chaser_row][grid->chaser_col]
                        - grid->ndepo_chaser4_nhx_h[grid->m][grid->chaser_row][grid->chaser_col])
                * ((double)(grid->climy) - 1850.0)/160.0;
            
        }else if(grid->climy>2010){
            ndepo_no3 = grid->ndepo_chaser4_noy_p[grid->m][grid->chaser_row][grid->chaser_col]
                + grid->ndepo_chaser4_ont_p[grid->m][grid->chaser_row][grid->chaser_col];
            ndepo_nh4 = grid->ndepo_chaser4_nhx_p[grid->m][grid->chaser_row][grid->chaser_col];
        }
        
        /* unit: g N ha-1 month-1 */
        loct->depo_no3[grid->m] = ndepo_no3;
        loct->depo_nh4[grid->m] = ndepo_nh4;
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
}

/*************************************************/
void f_n_mineralz(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Smas *mass, 
	struct Sflx *flux
){
	double f_c_min;
	
	/** litter **/
	if(mass->ltr > 0.01){
		f_c_min = flux->rl[grid->m] / mass->ltr;
	}else{
		f_c_min = 0.0;
	}
	
	/* g N ha-1 month-1 */
	flux->n_minerlz_lttr[grid->m] = mass->n_lttr * f_c_min;

	/** humus **/
	if(mass->msl > 0.01){
		f_c_min = flux->rh[grid->m] / mass->msl;
	}else{
		f_c_min = 0.0;
	}
	
	/* g N ha-1 month-1 */
	flux->n_minerlz_hums[grid->m] = mass->n_hums * f_c_min;
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
	fad_no3 = 0.1; /* 2010/04/06 by A.Ito */

	/* g N / ha */
	/* kg H2O / m2 */
	if((loct->sw30 + grid->prate_sfc[grid->m]) > 0.1){
		ntr_conc = (1.0 - fad_no3)*(mass->n_no3/10000.0) / 
					(loct->sw30 + grid->prate_sfc[grid->m]);
	}else{
		ntr_conc = 0.0;
	}
	/* g N / kg H2O */
	
	aa = loct->ro2[grid->m] * ntr_conc;
	if(aa > (mass->n_no3 * 0.95)){
		aa = mass->n_no3 * 0.95;
	}
    if(aa < 0.0){
        aa = 0.0;
    }
	
	/* g NO3-N m-2 month-1 */
	flux->n_leach[grid->m] = aa*10000.0;
	/* g NO3-N ha-1 month-1 */
}

/* plant N uptake *******************************************/
void f_n_uptake(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Flux *flux
){
	double f_temp, aa;
	double n_max = 1.0;
	double ks, navil;
	double uptake_no3, uptake_nh4;

	/* Raich et al. (1991): Appendix  Eq.1.16 */
	
	/*****
	Effect of N allocatiom to root ?
	*****/
	
	/* temperature factor */
	f_temp = exp(0.0693 * grid->tmp10_soil[grid->m]);
	/* soil diffusion parameter */
	ks = 0.90 * pow(loct->sw30 / grid->field_cap1, 3.0) + 0.1;
	
	/* NO3 uptake */
	navil = (mass->soil).n_no3;
	aa = navil * n_max * ks / (90.0 + ks*navil) * f_temp;
    if(aa>0.0 && aa<navil){
        uptake_no3 = aa;
    }else if(aa<0.0){
        uptake_no3 = 0.0;
    }else if(aa>=navil){
        uptake_no3 = navil;
    }
    
	/* g N ha-1 month-1 */
	(flux->c3).uptake_no3[grid->m] = uptake_no3;
	(flux->c4).uptake_no3[grid->m] = uptake_no3;
	
	/* NH4 uptake */
	navil = (mass->soil).n_nh4;
    aa = navil * n_max * ks / (90.0 + ks*navil) * f_temp;
    if(aa>0.0 && aa<navil){
        uptake_nh4 = aa;
    }else if(aa<0.0){
        uptake_nh4 = 0.0;
    }else if(aa>=navil){
        uptake_nh4 = navil;
    }
	/* g N ha-1 month-1 */
	(flux->c3).uptake_nh4[grid->m] = uptake_nh4;
	(flux->c4).uptake_nh4[grid->m] = uptake_nh4;
	
	(flux->plant).uptake_no3[grid->m] = (flux->c3).uptake_no3[grid->m] 
				+ (flux->c4).uptake_no3[grid->m];
	(flux->plant).uptake_nh4[grid->m] = (flux->c3).uptake_nh4[grid->m] 
				+ (flux->c4).uptake_nh4[grid->m];
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
	/* salvage (drwa back to storage) */
	flux->n_salvage[grid->m] = pchar->n_salvage * pchar->lf[grid->m] * 
					mass->n_cnpy;

	/* storage ****/
	flux->n_abdn_strg[grid->m] = (flux->lc[grid->m] + flux->lr[grid->m])
				* 1000000.0 / pchar->cn0_strg;
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
	n_opt = sqrt(amax*kn/rd_nsp) - kn;
	
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
		n_demand2 = 1.0 * (amax - a_nlmt);
		
		if(n_obtain > n_demand2){
			flux->n_alloc_cnpy[grid->m] = n_demand2;
		}else{
			flux->n_alloc_cnpy[grid->m] = n_obtain;
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
		n_opt = sqrt(amax*kn/rd_nsp) - kn;
		/* n_opt = 20.0; */
		
		/* leaf N concentration, mmol N m-2 */
		n_leaf_conc = pchar->n_conc_larea;
		a_nlmt = amax * n_leaf_conc / (kn + n_leaf_conc);
		/* a_nlmt = n_leaf_conc; */
		
		if(n_opt > n_leaf_conc){
			/* g N ha-1 month-1 */
			n_demand1 = (n_opt - n_leaf_conc) * 14.0 / 1000.0 * 10000.0;
			
			if(n_stock > n_demand1){
				flux->n_realloc[grid->m] = n_demand1;
			}else{
				flux->n_realloc[grid->m] = n_stock;
			}
		}else{  /* if(n_opt <= n_leaf_conc) */
			n_demand2 = 1.0 * (amax - a_nlmt);
			
			if(n_stock > n_demand2){
				flux->n_realloc[grid->m] = n_demand2;
			}else{
				flux->n_realloc[grid->m] = n_stock;
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
	f_immbl_nh4 = 0.001; */

	f_immbl_no3 = 0.005;
	f_immbl_nh4 = 0.004;

	flux->n_immbl[grid->m] = 0.2 * flux->n_minerlz_lttr[grid->m] + 
		0.4 * flux->n_minerlz_hums[grid->m] + 
		(f_immbl_no3 * mass->n_no3 + f_immbl_nh4 * mass->n_nh4) * MDN[grid->m];
    
    /* safe guard: 2014/05/28 by A.Ito */
    if(flux->n_immbl[grid->m] > mass->n_mcrb){
        flux->n_immbl[grid->m] = mass->n_mcrb;
    }
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

	flux->n_mcrb_abdn[grid->m] = 0.1 * f_temp * mass->n_mcrb;
    
    /* safe guard: 2014/05/28 by A.Ito */
    if(flux->n_mcrb_abdn[grid->m] > (0.5 * mass->n_mcrb)){
        flux->n_mcrb_abdn[grid->m] = (0.5 * mass->n_mcrb);
    }
}
