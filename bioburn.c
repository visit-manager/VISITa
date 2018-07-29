/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/*  Modified August 8, 2007						*/
/*  Modified November 3, 2007	by A.Ito			*/

#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

#define NOTICE 0

/* biomass burning *****************************************************/
void f_biomassburning(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Flux *flux
){
	short f;
	double aa, aad,  bb, cc, ss, n_fireseason;
	double fuel, fa_burnt;
	/* critical moisture */
	double me_crit[16] = {0.0, 
		0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3,
		0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2};
	double f_burnt_litter[16] = {0.0, 
		1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0,
		1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
	/* burnt fraction */
	/*
	Hoelzemann, J. J., M. G. Schultz, G. P. Brasseur, C. Granier, and M. Simon. 2004. 
	Global Wildland Fire Emission Model (GWEM): Evaluating the use of global area 
	burnt satellite data. Journal of Geophysical Research 109:10.1029/2003JD003666.
	*/
	double f_burnt_leaf[16] = {0.0, 
		0.6, 0.6, 0.3, 0.3, 0.3, 0.2, 0.2, 0.25,
		0.5, 0.5, 0.2, 0.2, 0.2, 0.5, 0.1};
	double f_burnt_wood[16] = {0.0, 
		0.25, 0.25, 0.1, 0.1, 0.1, 0.25, 0.25, 0.2,
		0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
	double f_burnt_root[16] = {0.0, 
		0.0, 0.0, 0.0, 0.0, 0.0, 0.05, 0.05, 0.0,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
	/* emission factor: g species/kg dry matter */
	double ef_co2[16] = {0.0, 
		1580.0, 1580.0, 1569.0, 1569.0, 1569.0, 1569.0, 1569.0, 1569.0,
		1663.0, 1663.0, 1663.0, 1663.0, 1663.0, 1663.0, 1663.0};
	double ef_co[16] = {0.0, 
		103.2, 103.2, 106.7, 106.7, 106.7, 106.7, 106.7, 106.7,
		61.6, 61.6, 61.6, 61.6, 61.6, 61.6, 61.6};
	
#if SENS_FIRE==1
	double ef_ch4[16] = {0.0, 
		9.0, 9.0, 4.70, 4.70, 4.70, 4.70, 4.70, 4.70,
		2.40, 2.40, 2.30, 2.30, 2.30, 2.30, 2.30};  /* alternative */
#else
	double ef_ch4[16] = {0.0, 
		6.80, 6.80, 4.70, 4.70, 4.70, 4.70, 4.70, 4.70,
		2.20, 2.20, 2.20, 2.20, 2.20, 2.20, 2.20}; /*  original */
#endif
	
	double ef_nmhc[16] = {0.0, 
		8.10, 8.10, 5.70, 5.70, 5.70, 5.70, 5.70, 5.70,
		3.40, 3.40, 3.40, 3.40, 3.40, 3.40, 3.40};
	double ef_oc[16] = {0.0, 
		5.20, 5.20, 9.14, 9.14, 9.14, 9.14, 9.14, 9.14,
		3.30, 3.30, 3.30, 3.30, 3.30, 3.30, 3.30};
	double ef_bc[16] = {0.0, 
		0.66, 0.66, 0.56, 0.56, 0.56, 0.56, 0.56, 0.56,
		0.47, 0.47, 0.47, 0.47, 0.47, 0.47, 0.47};
		
	double ef_nox[16] = {0.0, 
		1.85, 1.85, 3.00, 3.00, 3.00, 3.00, 3.00, 3.00,
		2.32, 2.32, 2.32, 2.32, 2.32, 2.32, 2.32};
	double ef_so2[16] = {0.0, 
		0.57, 0.57, 1.00, 1.00, 1.00, 1.00, 1.00, 1.00,
		0.71, 0.71, 0.71, 0.71, 0.71, 0.71, 0.71};
	double ef_pm25[16] = {0.0, 
		9.10, 9.10, 12.99, 12.99, 12.99, 12.99, 12.99, 12.99,
		4.90, 4.90, 4.90, 4.90, 4.90, 4.90, 4.90};
	double ef_tpm[16] = {0.0, 
		8.50, 8.50, 17.62, 17.62, 17.62, 17.62, 17.62, 17.62,
		9.20, 9.20, 9.20, 9.20, 9.20, 9.20, 9.20};
	double ef_tec[16] = {0.0, 
		6.60, 6.60, 8.28, 8.28, 8.28, 8.28, 8.28, 8.28,
		3.70, 3.70, 3.70, 3.70, 3.70, 3.70, 3.70};
		
    double ef_n2o[16] = {0.0,
        0.2, 0.2, 0.16, 0.16, 0.16, 0.41, 0.41, 38.0,
        0.2, 0.2, 0.2, 0.2, 0.2, 0.2, 0.2};
    
    /* ammonium: ref. Akagi et al. (2011): 2018/07/29 by A.Ito */
    double ef_nh4[16] = {0.0,
        0.00564, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
        0.0035, 0.00397, 0.0, 0.0, 0.0, 0.0, 0.0};

	/* burning efficiency */
	/* double burn_eff[16] = {0.0,
		0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
		0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; */
    double burn_eff[16] = {0.0,
        0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
        0.35, 0.35, 0.35, 0.35, 0.20, 0.25, 0.0}; /* */ /* test: 2018/07/23 by A.Ito */

	double closs_leaf, closs_wood, closs_root, closs_litter, prm_ensen;
	
	extern double MDN[ASTEP];
	
	/* veg_sage:
	0	Ocean
	1	Tropical Evergreen Forest/Woodland
	2	Tropical Deciduous Forest/Woodland
	3	Temperate Broadleaf Evergreen Forest/Woodland
	4	Temperate Needleleaf Evergreen Forest/Woodland
	5	Temperate Deciduous Forest/Woodland
	6	Boreal Evergreen Forest/Woodland
	7	Boreal Deciduous Forest/Woodland
	8	Evergreen/Deciduous Mixed Forest/Woodland
	9	Savanna
	10	Grassland/Steppe
	11	Dense Shrubland
	12	Open Shrubland
	13	Tundra
	14	Desert
	15	Polar Desert/Rock/Ice
	*/

	/* Estimation of burnt area
	   using LPJ-type paramaterization by 
	   Thonicke, K., S. Venevsky, S. Sitch, and W. Cramer. 2001. 
	   The role of fire disturbance for global vegetation dynamics: 
	   coupling fire into a Dynamic Global Vegetation Model. 
	   Global Ecology and Biogeography 10:661-677.
	*/
	n_fireseason = 0.0;
	for(f=0;f<ASTEP;f++){	
		/* fuel load: need 200 g dm/m2 */
		/* fuel = ((mass->soil).ltr_m[f] + (mass->plant).mfol[f] + (mass->plant).mstm[f] 
				+ (mass->plant).mrot[f]) * dmTc * 100.0; */ /* need 200 g dm/m2 */

		fuel = ((mass->soil).ltr_m[f] + (mass->plant).mstm[f]) * dmTc * 100.0; /* A.Ito (2010/06/17) */

		if(fuel >= 200.0){	/* fire threshold: 2007/11/03 */
			/* volumetric upper soil (litter-fuel) water content */
			aa = loct->msw30[f]/grid->field_cap1;
			
			/* Eq.1 in Thonicke  */
			aad = (0.4994 * (aa*100.0) + 1.02)/100.0;	
			
            /* corrected by A.Ito (2012/02/01) based on Hamada-san's comment */
            if(grid->veg_sage>=1 && grid->veg_sage<=15){
                bb = aad/me_crit[grid->veg_sage];
            }else{
                /* ocean */
                bb = 0.0;
            }
			/* Eq.2 in Thonicke  */
			cc = exp(-PI*bb*bb);
			
			if(cc<=0.0){
				cc = 0.0;
			}else if(cc >=1.0){
				cc = 1.0;
			}
			
			/* fire season length */
			flux->day_fire[f] = MDN[f] * cc;
			n_fireseason += MDN[f] * cc;
		}else{
			flux->day_fire[f] = 0.0;
		}
	}
	
	/* annual fraction of fire season */
	if(n_fireseason >= 0.05){
		ss = n_fireseason/YDN;	
		
		if(ss<=0.0){
			ss = 0.0;
		}else if(ss >=1.0){
			ss = 1.0;
		}
		
		aa = ss - 1.0;
		/* Eq.8 in Thonicke  */
		bb = 0.45 * pow(aa, 3.0) + 2.83 * pow(aa, 2.0) + 2.96 * aa + 1.04;
		
		/* fractional area burnt, Eq.6 */
		fa_burnt = ss * exp(aa/bb);
        
        /* constraint by GFED4s: 2018/05/18 by A.Ito */
        if(EX_FIRE_GFED >= 1){
            if(grid->simy <= 1997){
                loct->fb_base = fa_burnt;
            }
            if(grid->simy >=1998 && grid->simy <= 2016){
                fa_burnt = loct->fb_base * bf_gfed4s[grid->simy - 1997][grid->reg_g];
            }
            if(grid->simy >=2017){
                fa_burnt = loct->fb_base * bf_gfed4s[2016 - 1997][grid->reg_g];
            }
            
            if(EX_FIRE_GFED == 2){
                fa_burnt *= 0.73; /* adjust global total burnt area to GFED4s */
            }
        }

		if(fa_burnt<=0.0){
			fa_burnt = 0.0;
		}else if(fa_burnt >=1.0){
			fa_burnt = 1.0;
		}
	}else{
		fa_burnt = 0.0;
	}
    
    /* parameter ensemble: 2014/11/19 by A.Ito */
    prm_ensen = 1.0;
    if(PARAM_PTB == 6){
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
        prm_ensen = 1.0 + 0.3 * f_pert[0];
    }else{
        prm_ensen = 1.0;
    }
	
	/******************************/
	flux->f_burnt = fa_burnt;
	
	for(f=0;f<ASTEP;f++){
		/* seasonal change: fire-day length */
		if(n_fireseason>=0.05){
			aa = flux->day_fire[f]/n_fireseason;
		}else{
			aa = 0.0;
		}
		/* fractional area burnt */
		/* flux->a_burnt[f] = flux->f_burnt*aa*(1.0 - grid->f_crop_con); */
        /* revised (after comments by E.Kato): 2013/10/02 by A.Ito */
		flux->a_burnt[f] = flux->f_burnt * aa;
        
        /* burnt fractio for woods: 2017/11/30 */
        flux->wa_burnt[f] = flux->a_burnt[f] * f_burnt_wood[grid->veg_sage];
  
		/* g/kg = kg/Mg, per ha */
		/* CO2 emission */
		flux->bb_co2_litter[f] = flux->a_burnt[f] * (mass->soil).ltr_m[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_litter[grid->veg_sage] * ef_co2[grid->veg_sage] * prm_ensen;
		flux->bb_co2_leaf[f] = flux->a_burnt[f] * (mass->plant).mfol[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_leaf[grid->veg_sage] * ef_co2[grid->veg_sage] * prm_ensen;
		flux->bb_co2_wood[f] = flux->a_burnt[f] * (mass->plant).mstm[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_wood[grid->veg_sage] * ef_co2[grid->veg_sage] * prm_ensen;
		flux->bb_co2_root[f] = flux->a_burnt[f] * (mass->plant).mrot[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_root[grid->veg_sage] * ef_co2[grid->veg_sage] * prm_ensen;

		/* CO emission */
		flux->bb_co_litter[f] = flux->a_burnt[f] * (mass->soil).ltr_m[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_litter[grid->veg_sage] * ef_co[grid->veg_sage] * prm_ensen;
		flux->bb_co_leaf[f] = flux->a_burnt[f] * (mass->plant).mfol[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_leaf[grid->veg_sage] * ef_co[grid->veg_sage] * prm_ensen;
		flux->bb_co_wood[f] = flux->a_burnt[f] * (mass->plant).mstm[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_wood[grid->veg_sage] * ef_co[grid->veg_sage] * prm_ensen;
		flux->bb_co_root[f] = flux->a_burnt[f] * (mass->plant).mrot[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_root[grid->veg_sage] * ef_co[grid->veg_sage] * prm_ensen;

		/* CH4 emission */
		flux->bb_ch4_litter[f] = flux->a_burnt[f] * (mass->soil).ltr_m[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_litter[grid->veg_sage] * ef_ch4[grid->veg_sage] * prm_ensen;
		flux->bb_ch4_leaf[f] = flux->a_burnt[f] * (mass->plant).mfol[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_leaf[grid->veg_sage] * ef_ch4[grid->veg_sage] * prm_ensen;
		flux->bb_ch4_wood[f] = flux->a_burnt[f] * (mass->plant).mstm[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_wood[grid->veg_sage] * ef_ch4[grid->veg_sage] * prm_ensen;
		flux->bb_ch4_root[f] = flux->a_burnt[f] * (mass->plant).mrot[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_root[grid->veg_sage] * ef_ch4[grid->veg_sage] * prm_ensen;

		/* NMHC emission */
		flux->bb_nmhc_litter[f] = flux->a_burnt[f] * (mass->soil).ltr_m[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_litter[grid->veg_sage] * ef_nmhc[grid->veg_sage] * prm_ensen;
		flux->bb_nmhc_leaf[f] = flux->a_burnt[f] * (mass->plant).mfol[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_leaf[grid->veg_sage] * ef_nmhc[grid->veg_sage] * prm_ensen;
		flux->bb_nmhc_wood[f] = flux->a_burnt[f] * (mass->plant).mstm[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_wood[grid->veg_sage] * ef_nmhc[grid->veg_sage] * prm_ensen;
		flux->bb_nmhc_root[f] = flux->a_burnt[f] * (mass->plant).mrot[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_root[grid->veg_sage] * ef_nmhc[grid->veg_sage] * prm_ensen;

		/* OC emission */
		flux->bb_oc_litter[f] = flux->a_burnt[f] * (mass->soil).ltr_m[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_litter[grid->veg_sage] * ef_oc[grid->veg_sage] * prm_ensen;
		flux->bb_oc_leaf[f] = flux->a_burnt[f] * (mass->plant).mfol[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_leaf[grid->veg_sage] * ef_oc[grid->veg_sage] * prm_ensen;
		flux->bb_oc_wood[f] = flux->a_burnt[f] * (mass->plant).mstm[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_wood[grid->veg_sage] * ef_oc[grid->veg_sage] * prm_ensen;
		flux->bb_oc_root[f] = flux->a_burnt[f] * (mass->plant).mrot[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_root[grid->veg_sage] * ef_oc[grid->veg_sage] * prm_ensen;

		/* BC emission */
		flux->bb_bc_litter[f] = flux->a_burnt[f] * (mass->soil).ltr_m[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_litter[grid->veg_sage] * ef_bc[grid->veg_sage] * prm_ensen;
		flux->bb_bc_leaf[f] = flux->a_burnt[f] * (mass->plant).mfol[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_leaf[grid->veg_sage] * ef_bc[grid->veg_sage] * prm_ensen;
		flux->bb_bc_wood[f] = flux->a_burnt[f] * (mass->plant).mstm[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_wood[grid->veg_sage] * ef_bc[grid->veg_sage] * prm_ensen;
		flux->bb_bc_root[f] = flux->a_burnt[f] * (mass->plant).mrot[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_root[grid->veg_sage] * ef_bc[grid->veg_sage] * prm_ensen;

		/* NOx emission */
		flux->bb_nox_litter[f] = flux->a_burnt[f] * (mass->soil).ltr_m[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_litter[grid->veg_sage] * ef_nox[grid->veg_sage] * prm_ensen;
		flux->bb_nox_leaf[f] = flux->a_burnt[f] * (mass->plant).mfol[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_leaf[grid->veg_sage] * ef_nox[grid->veg_sage] * prm_ensen;
		flux->bb_nox_wood[f] = flux->a_burnt[f] * (mass->plant).mstm[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_wood[grid->veg_sage] * ef_nox[grid->veg_sage] * prm_ensen;
		flux->bb_nox_root[f] = flux->a_burnt[f] * (mass->plant).mrot[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_root[grid->veg_sage] * ef_nox[grid->veg_sage] * prm_ensen;

        /* N2O emission */
        flux->bb_n2o_litter[f] = flux->a_burnt[f] * (mass->soil).ltr_m[f]/cTdm * burn_eff[grid->veg_sage]
            * f_burnt_litter[grid->veg_sage] * ef_n2o[grid->veg_sage] * prm_ensen;
        flux->bb_n2o_leaf[f] = flux->a_burnt[f] * (mass->plant).mfol[f]/cTdm * burn_eff[grid->veg_sage]
            * f_burnt_leaf[grid->veg_sage] * ef_n2o[grid->veg_sage] * prm_ensen;
        flux->bb_n2o_wood[f] = flux->a_burnt[f] * (mass->plant).mstm[f]/cTdm * burn_eff[grid->veg_sage]
            * f_burnt_wood[grid->veg_sage] * ef_n2o[grid->veg_sage] * prm_ensen;
        flux->bb_n2o_root[f] = flux->a_burnt[f] * (mass->plant).mrot[f]/cTdm * burn_eff[grid->veg_sage]
            * f_burnt_root[grid->veg_sage] * ef_n2o[grid->veg_sage] * prm_ensen;

		/* SO2 emission */
		flux->bb_so2_litter[f] = flux->a_burnt[f] * (mass->soil).ltr_m[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_litter[grid->veg_sage] * ef_so2[grid->veg_sage] * prm_ensen;
		flux->bb_so2_leaf[f] = flux->a_burnt[f] * (mass->plant).mfol[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_leaf[grid->veg_sage] * ef_so2[grid->veg_sage] * prm_ensen;
		flux->bb_so2_wood[f] = flux->a_burnt[f] * (mass->plant).mstm[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_wood[grid->veg_sage] * ef_so2[grid->veg_sage] * prm_ensen;
		flux->bb_so2_root[f] = flux->a_burnt[f] * (mass->plant).mrot[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_root[grid->veg_sage] * ef_so2[grid->veg_sage] * prm_ensen;

		/* PM2.5 emission */
		flux->bb_pm25_litter[f] = flux->a_burnt[f] * (mass->soil).ltr_m[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_litter[grid->veg_sage] * ef_pm25[grid->veg_sage] * prm_ensen;
		flux->bb_pm25_leaf[f] = flux->a_burnt[f] * (mass->plant).mfol[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_leaf[grid->veg_sage] * ef_pm25[grid->veg_sage] * prm_ensen;
		flux->bb_pm25_wood[f] = flux->a_burnt[f] * (mass->plant).mstm[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_wood[grid->veg_sage] * ef_pm25[grid->veg_sage] * prm_ensen;
		flux->bb_pm25_root[f] = flux->a_burnt[f] * (mass->plant).mrot[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_root[grid->veg_sage] * ef_pm25[grid->veg_sage] * prm_ensen;

		/* TPM (total particulate matter) emission */
		flux->bb_tpm_litter[f] = flux->a_burnt[f] * (mass->soil).ltr_m[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_litter[grid->veg_sage] * ef_tpm[grid->veg_sage] * prm_ensen;
		flux->bb_tpm_leaf[f] = flux->a_burnt[f] * (mass->plant).mfol[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_leaf[grid->veg_sage] * ef_tpm[grid->veg_sage] * prm_ensen;
		flux->bb_tpm_wood[f] = flux->a_burnt[f] * (mass->plant).mstm[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_wood[grid->veg_sage] * ef_tpm[grid->veg_sage] * prm_ensen;
		flux->bb_tpm_root[f] = flux->a_burnt[f] * (mass->plant).mrot[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_root[grid->veg_sage] * ef_tpm[grid->veg_sage] * prm_ensen;

		/* TEC emission */
		flux->bb_tec_litter[f] = flux->a_burnt[f] * (mass->soil).ltr_m[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_litter[grid->veg_sage] * ef_tec[grid->veg_sage] * prm_ensen;
		flux->bb_tec_leaf[f] = flux->a_burnt[f] * (mass->plant).mfol[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_leaf[grid->veg_sage] * ef_tec[grid->veg_sage] * prm_ensen;
		flux->bb_tec_wood[f] = flux->a_burnt[f] * (mass->plant).mstm[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_wood[grid->veg_sage] * ef_tec[grid->veg_sage] * prm_ensen;
		flux->bb_tec_root[f] = flux->a_burnt[f] * (mass->plant).mrot[f]/cTdm * burn_eff[grid->veg_sage] 
			* f_burnt_root[grid->veg_sage] * ef_tec[grid->veg_sage] * prm_ensen;
		
		/* carbon budget ****************************************/
		if(NECB_BB == 1){
			closs_leaf = flux->bb_co2_leaf[f]*12.0/44.0/1000.0 + flux->bb_co_leaf[f]*12.0/28.0/1000.0 
				+ flux->bb_ch4_leaf[f]*12.0/16.0/1000.0 + flux->bb_bc_leaf[f]/1000.0;
			closs_wood = flux->bb_co2_wood[f]*12.0/44.0/1000.0 + flux->bb_co_wood[f]*12.0/28.0/1000.0 
				+ flux->bb_ch4_wood[f]*12.0/16.0/1000.0 + flux->bb_bc_wood[f]/1000.0;
			closs_root = flux->bb_co2_root[f]*12.0/44.0/1000.0 + flux->bb_co_root[f]*12.0/28.0/1000.0 
				+ flux->bb_ch4_leaf[f]*12.0/16.0/1000.0 + flux->bb_bc_leaf[f]/1000.0;
			closs_litter = flux->bb_co2_litter[f]*12.0/44.0/1000.0 + flux->bb_co_litter[f]*12.0/28.0/1000.0 
				+ flux->bb_ch4_litter[f]*12.0/16.0/1000.0 + flux->bb_bc_litter[f]/1000.0;
			
			/***************/
			(mass->c3).fol -= closs_leaf;
			if((mass->c3).fol < INT_C){
				(mass->c3).fol = INT_C;
			}
			(mass->c3).mfol[f] = (mass->c3).fol;
			
			(mass->c4).fol -= closs_leaf;
			if((mass->c4).fol < INT_C){
				(mass->c4).fol = INT_C;
			}
			(mass->c4).mfol[f] = (mass->c4).fol;
			
			(mass->plant).fol = (mass->c3).fol*loct->c3ptn[f] + (mass->c4).fol*loct->c4ptn[f];
			(mass->plant).mfol[f] = (mass->c3).mfol[f]*loct->c3ptn[f] + (mass->c4).mfol[f]*loct->c4ptn[f];
			
			/***************/
			(mass->c3).stm -= closs_wood;
			if((mass->c3).stm < INT_C){
				(mass->c3).stm = INT_C;
			}
			(mass->c3).mstm[f] = (mass->c3).stm;
			
			(mass->c4).stm -= closs_wood;
			if((mass->c4).stm < INT_C){
				(mass->c4).stm = INT_C;
			}
			(mass->c4).mstm[f] = (mass->c4).stm;
			
			(mass->plant).stm = (mass->c3).stm*loct->c3ptn[f] + (mass->c4).stm*loct->c4ptn[f];
			(mass->plant).mstm[f] = (mass->c3).mstm[f]*loct->c3ptn[f] + (mass->c4).mstm[f]*loct->c4ptn[f];

			/***************/
			(mass->c3).rot -= closs_root;
			if((mass->c3).rot < INT_C){
				(mass->c3).rot = INT_C;
			}
			(mass->c3).mrot[f] = (mass->c3).rot;
			
			(mass->c4).rot -= closs_root;
			if((mass->c4).rot < INT_C){
				(mass->c4).rot = INT_C;
			}
			(mass->c4).mrot[f] = (mass->c4).rot;
			
			(mass->plant).rot = (mass->c3).rot*loct->c3ptn[f] + (mass->c4).rot*loct->c4ptn[f];
			(mass->plant).mrot[f] = (mass->c3).mrot[f]*loct->c3ptn[f] + (mass->c4).mrot[f]*loct->c4ptn[f];	
			
			/***************/
			(mass->soil).ltr -= closs_litter;
			if((mass->soil).ltr < INT_C){
				(mass->soil).ltr = INT_C;
			}
			(mass->soil).ltr_m[f] = (mass->soil).ltr;
		}
	}
}
