/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
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

/* ecophysiological vegetation processes *******************************/
void ecophysiology(
	struct Grid *grid, 
	struct Loct *loct,
	struct Pchar *pchar, 
	struct Pmas *mass
){
	long g;
	double aaa, bbb;
	double sinb, ke_b1, ke_b2, irr_b, rfl_b, apar, fapar, eff_k;
	
	/* give leaf area index (LAI), m2 m-2*/
	mass->lai[grid->m] = lai_mass(grid, mass, pchar);
	
	/* canopy radiation absorption */
	sinb = sin(grid->lat*dTr)*sin(grid->sl_dec[grid->m]*dTr) 
			+ cos(grid->lat*dTr)*cos(grid->sl_dec[grid->m]*dTr)*1.0;
	sinb = (sinb<=1.0)?sinb:1.0; 
	sinb = (sinb>=-1.0)?sinb:-1.0;
	
	if(sinb>0.0 && grid->par[grid->m]>0.0 && mass->lai[grid->m]>0.0){
		ke_b1 = 0.5/sinb;
		ke_b2 = 0.46/sinb;
		irr_b = (1.0 - sqrt(1.0 - 0.15))/(1.0 + sqrt(1.0 - 0.15));
		rfl_b = 1.0-exp(-2.0*irr_b*ke_b1)/(1.0 + ke_b1);
		
		pchar->apar_bp[grid->m] = (1.0 - rfl_b)*grid->par_bp[grid->m]*(1.0 - exp(-ke_b2*mass->lai[grid->m]));
		pchar->apar_dp[grid->m] = (1.0 - 0.036)*grid->par_dp[grid->m]*(1.0 - exp(-0.719*mass->lai[grid->m]));
		
		apar = pchar->apar_bp[grid->m] + pchar->apar_dp[grid->m];
		fapar = apar/grid->par[grid->m];
				
		eff_k = -1.0*log(1.0 - fapar)/mass->lai[grid->m];
		eff_k = (eff_k>=0.46)?eff_k:0.1;
		eff_k = (eff_k<=5.0)?eff_k:10.0;
	}else{
		pchar->apar_bp[grid->m] = 0.0;
		pchar->apar_dp[grid->m] = 0.0;
		fapar = 0.0;
		eff_k = irr_attn(grid, loct, pchar);
	}
	
	/* give irradiance attenuation coefficient */
	if(EFF_K == 0){
		pchar->eK[grid->m] = irr_attn(grid, loct, pchar);
		pchar->fapar[grid->m] = (1.0 - pchar->albedo)*(1.0 - exp(-pchar->eK[grid->m]*mass->lai[grid->m]));
	}else if(EFF_K==1){
		pchar->eK[grid->m] = eff_k;
		pchar->fapar[grid->m] = fapar;
	}
	
	/* leaf N concentration */
	f_n_leaf_conc(grid, pchar, mass);
	
	/** initial ci: 70% of ambient level **/
	pchar->ci[grid->m] = loct->aco2[grid->m]*0.7; 
	
	/* stabilization of single-leaf processes */
	for(g=0;g<6;g++){			
		 /* give quantum yield , mol CO2 mon photon-1*/
		quantum_yield(grid, pchar);	
		
		/* give light-saturated photosynthesis rate, micromol m-2 s-1 */
		pc_sat(grid, loct, pchar);

		/* canopy-top photosynthetic rate*/
		aaa = pchar->psat[grid->m]*pchar->lue[grid->m]*grid->par[grid->m];
		bbb = pchar->psat[grid->m] + pchar->lue[grid->m]*grid->par[grid->m];
		if(bbb>0.0){
			pchar->ptop = aaa/bbb;
		}else{
			pchar->ptop = 0.0;
		}

		/* give stomatal conductance, mmol H2O m-2 s-1 */
		stom_cond(grid, loct, pchar); 

		/* give intercellular CO2 concentration, ppmv */
		incel_cdc(grid, loct, pchar); 
	}
		
	/** photosynthetic 13C discrimination **/
	photo_13c_frac(grid, loct, pchar);
	
	/** canopy conductance **/	
	pchar->gc[grid->m] = canopy_cond(grid, loct, pchar, mass);

	/** plant respiration **/
	qten_ar(grid, pchar); /*Q10*/
	spcfc_res_mass(pchar, mass); /* woody specific respiration rate */

	/** litterfall of plant respiration **/
	mortality(grid, pchar);
	
	/* GPP by de Pury & Farquhar scheme */
	if(DF97==1){
		f_df97_gpp(2, grid, loct, pchar, mass);
		f_df97_gpp(3, grid, loct, pchar, mass);
	}

	/** optimum leaf area index **/	
	opt_lai(grid,loct, pchar);
}

/* leaf area index **********************************************/
double lai_mass(
	struct Grid *grid, 
	struct Pmas *mass, 
	struct Pchar *pchar
){
	double sla, lai_est;
	
	/** specific leaf area as a function of... what? **/
	sla = pchar->sla;
	
	lai_est = sla*mass->fol*dmTc/100.0/2.0;
	lai_est = (lai_est>=0.0)?lai_est:0.0;
	
	/* dmTc: dry-matter to carbon */
	/* 100.0: cm2 g dm-1 to Mg ha-1 base */
	/* 2.0: single-sided leaf area */	
	
	/* sensitivity analysis: prescribed LAI **/
	if(SENS == 7 && (grid->climy>=2001) ){
		lai_est = mass->lai0[grid->m];
	}
	
	return(lai_est);
}

/* irradiance attenuation coefficient ********************************/
double irr_attn(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	double aaa, bbb;
	
	/* a function of solar hight angle */
	aaa = sin(grid->sl_hgt[grid->m]*dTr); 
	aaa = (aaa<=1.0)?aaa:1.0; 
	aaa = (aaa>=0.3)?aaa:0.3; /* to avoid extreme values*/
	bbb = pchar->eK0/aaa;
	
	return (bbb);
}

/* intercellular CO2 concentration ********************************/
void incel_cdc(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *plant
){
	double ci, gs_co2;
	/* give intercellular CO2 concentration, as a function of ambient 
		CO2 level and stomatal conductance */
	
	gs_co2 = plant->gs[grid->m]/1.56; 
	/* 1.56: conversion from H2O to CO2 conductance */
	
	ci = loct->aco2[grid->m]-(plant->ptop/(gs_co2/1000.0));
	/* 1000.0: conbert from mmol to micro mol */
	
	ci = (ci>=0.0)?ci:0.0;
	ci = (ci<=loct->aco2[grid->m])?ci:loct->aco2[grid->m];
		
	plant->ci[grid->m] = ci; 
}

/* quantum yield C3 and C4 *******************************************/
void quantum_yield(
	struct Grid *grid, 
	struct Pchar *pchar
){
	double eftem, efci;

	/** plant_type: 3=C3, 4=C4, (5=CAM) **/
	if(pchar->phototype == 3){
		/* temperature dependence */
		eftem = (52.0 - grid->tmp_sfc[grid->m])/(3.5 + 0.75*(52.0 - grid->tmp_sfc[grid->m])); 
		/* CO2 dependence */
		efci = pchar->ci[grid->m]/(90.0+0.6*pchar->ci[grid->m]); 
		/* 3.5, 52.0, etc.: empirical parameters */
	}else if(pchar->phototype == 4){ 
		/* insensitive QE of C4 species */
		eftem = 1.0;
		efci = 1.0;
	}
	/* give quantum yield */
	pchar->lue[grid->m] = pchar->lue0*eftem*efci;
}

/* stomatal conductance **********************************************/
void stom_cond(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	double b1d, cc;
	
	/* stomatal conductance model by Ball, Woodraw, and Berry (1987) */
	/*
	 Leuning, R. (1995), A critical appraisal of a combined stomatal-photosynthesis 
	 model for C3 plants, Plant, Cell and Environment, 18, 339-355.
	*/
	b1d = pchar->gs_b1/((loct->aco2[grid->m] - pchar->cmpcd[grid->m])*(1.0+loct->vpd[grid->m]/pchar->gs_b2)); /* */
	/* insensitive to CO2 */
	/* b1d=plant->gs_b1/(( 350.0 - 40.0 )*(1.0+loct->vpd[grid->m]/plant->gs_b2)); */

	/** add soil water factor **/
	cc = 1.0; /* not defined yet */

	if(pchar->psat[grid->m]>0.0){
		pchar->gs[grid->m] = pchar->gs_b0 + b1d*pchar->ptop/pchar->psat[grid->m]*cc;
	}else{
		pchar->gs[grid->m] = pchar->gs_b0; 
	}
}

/* canopy conductance ******************************************/
double canopy_cond(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	double aaa, sss, ttt, uuu, vvv, lue_gs, canopy_cond;
	
	/* NOTE: integrate leaf stomatal conductance with considering light attenuation in the canopy */
	/* aaa=plant->gs_b0+plant->gs_b1/(loct->aco2[grid->m]-plant->cmpcd[grid->m]); */
	aaa = pchar->gs_b0 + pchar->gs_b1/(350.0 - 40.0);
	lue_gs = pchar->lue[grid->m]*(aaa/pchar->pmax);
	
	if(mass->lai[grid->m]>0.0){
		sss = 2.0*pchar->gs[grid->m] / pchar->eK[grid->m]; 
		ttt = 1.0 + sqrt(1.0 + pchar->eK[grid->m]*lue_gs*grid->par[grid->m]/pchar->gs[grid->m]);
		vvv = -1.0*pchar->eK[grid->m]*mass->lai[grid->m];
		uuu = 1.0 + sqrt(1.0 + pchar->eK[grid->m]*lue_gs*grid->par[grid->m]*exp(vvv)/pchar->gs[grid->m]);
		canopy_cond = sss*log(ttt/uuu);
	}else{
		/* no leaf, no conductance*/
		canopy_cond = 0.0;
	} 
	
	return(canopy_cond);
}

/* optimum LAI by Kuroiwa (1966) ************************************/
void opt_lai(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar
){
	double aaa, bbb, ccc, cc4, ddd, eee;
	double arm, arg, ar;	
	double psat, lue;
	
	if(DF97==1){
		psat = pchar->psat_df[grid->m];
		lue = pchar->lue_df[grid->m];
	}else{
		psat = pchar->psat[grid->m];
		lue = pchar->lue[grid->m];
	}
	
	aaa = 1.0/pchar->eK[grid->m];
	bbb = pchar->eK[grid->m]*lue*grid->par[grid->m]; 
	
	/* daily respiratory cost */
	/* printf("%lf %lf\n", plant->qTc[grid->m], grid->tmp_sfc[grid->m]); */
	eee = log(pchar->qTc[grid->m])/10.0*(grid->tmp_sfc[grid->m] - 15.0);
	arm = pchar->rmf*exp(eee)/1000.0*dmTc*10000.0/(pchar->sla);
	arg = pchar->lf[grid->m]*dmTc*10000.0/(pchar->sla)*(1.0 + pchar->rgf);
	ar = arm + arg;

	cc4 = (psat*grid->dlen[grid->m])/(psat*grid->dlen[grid->m] - ar*24.0);
	ccc = psat*(cc4 - 1.0);
	
	if(ccc > 0.0){
		ddd = bbb/ccc;
		ddd = (ddd>1.0)?ddd:1.0;
		pchar->opt_lai[grid->m] = aaa*log(ddd);
	}else{
		pchar->opt_lai[grid->m] = 0.0;
	}
}

/* Q10 of autotrophic respiration **********************************************/
void qten_ar(
	struct Grid *grid, 
	struct Pchar *pchar
){
	double aaa;
	
	/* larger at cool and smaller at warm */
	aaa = exp(-0.009*(grid->tmp_sfc[grid->m]-15.0));
		
	pchar->qTf[grid->m] = pchar->qTf0*aaa;
	pchar->qTc[grid->m] = pchar->qTc0*aaa;
	pchar->qTr[grid->m] = pchar->qTr0*aaa;
	
	if(T_R==1){
		pchar->qTf[grid->m]*=0.9;
		pchar->qTc[grid->m]*=0.9;
		pchar->qTr[grid->m]*=0.9;
	}else if(T_R==2){
		pchar->qTf[grid->m]*=1.1;
		pchar->qTc[grid->m]*=1.1;
		pchar->qTr[grid->m]*=1.1;
	}
}

/* lowering specific respiration rate of massive biomass with heart wood tissue *******/
void spcfc_res_mass(
	struct Pchar *pchar, 
	struct Pmas *mass
){
	double powstm, powrot;
	double stm_sap, stm_hrt, rot_sap, rot_hrt;
	
	pchar->rmf = pchar->rmf0;
	
	/* specific respiration increasing in a power of 2/3 manner */
	powstm = 1.0 - 0.33334*mass->stm/(50.0 + mass->stm);
	powrot = 1.0 - 0.33334*mass->rot/(50.0 + mass->rot);
	
	stm_sap = pow(mass->stm, powstm); /* sapwood mass in stem */
	stm_hrt = mass->stm - stm_sap; /* heartwood mass in stem */
	rot_sap = pow(mass->rot, powrot); /* sapwood mass in root */
	rot_hrt = mass->rot - rot_sap; /* heartwood mass in root */
	
	pchar->rmc = (pchar->rmc_s*stm_sap + pchar->rmc_h*stm_hrt)/(mass->stm + 0.00001);
	pchar->rmr = (pchar->rmr_s*rot_sap + pchar->rmr_h*rot_hrt)/(mass->rot + 0.00001);
}

/* leaf N concentration ********************************************/
void f_n_leaf_conc(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	/* mmol N m-2 leaf area */
	if(mass->lai[grid->m] > 0.01){
		pchar->n_conc_larea = mass->n_cnpy / 10000.0 / 14.0 * 1000.0 / mass->lai[grid->m];
	}else{
		pchar->n_conc_larea = 1.0;
	}
	
	pchar->n_conc_larea_m[grid->m] = pchar->n_conc_larea;
	
	/* mmol N g-1 leaf weight */
	if(mass->mfol[grid->m] > 0.01){
		pchar->n_conc_lmass = mass->n_cnpy / 14.0 *1000.0 / (mass->mfol[grid->m] * 1000000.0);
	}else{
		pchar->n_conc_lmass = 1.0;
	}
}

/* leaf age **************************************************************/
/* used for BVOC emission */
void f_leaf_age(
	short update, 
	struct Pchar *pchar, 
	struct Pmas *mass, 
	double addshed
){
	long f;
	double aaa, bbb;
	
	/* age update, once per month */
	if(update == 1){
		/* update leaf age */
		pchar->fleaf_age[48] += pchar->fleaf_age[47];
		for(f=47;f>=1;f--){
			pchar->fleaf_age[f] = pchar->fleaf_age[f-1];
		}
		pchar->fleaf_age[0] = 0.0;
	}
	
	if(addshed > 0.0){
		/* new leaf */
		
		pchar->fleaf_age[0] += addshed;
	}else{
		/* senescent leaf */
		
		aaa = 0.0;
		for(f=0;f<=48;f++){
			aaa += pchar->fleaf_age[f];
		}
		
		if(aaa > 0.0 && mass->fol > 0.0){
			
			bbb = fabs(addshed) / aaa;
			
			for(f=0;f<=48;f++){
				/* modification by A.Ito (2009/06/03) based on E.Kato (2008/11/19) **/
				if((1.0 - bbb) > 0.0){
					pchar->fleaf_age[f] *= (1.0 - bbb);
				}else{
					pchar->fleaf_age[f] = 0.0;
				}
			}
			
			/* for(f=0;f<=48;f++){
				pchar->fleaf_age[f] *= mass->fol / aaa;
			} */
		}else{
			for(f=0;f<=48;f++){
				pchar->fleaf_age[f] = 0.0;
			}
		}
	}
}
