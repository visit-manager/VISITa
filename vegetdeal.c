/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"structure.h"
#include"prototype.h"

/* aggregate C3 and C4 community *************************************************/
void f_plant_stand_budget(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Flux *flux
){
	long f;
	
	f = grid->m;
	
	/** mass **************************************************************************/
	mass->lai_p = (mass->c3).lai[f]*loct->c3ptn[f] + (mass->c4).lai[f]*loct->c4ptn[f];
	/* carbon */
	(mass->plant).fol = (mass->c3).fol*loct->c3ptn[f] + (mass->c4).fol*loct->c4ptn[f];
	(mass->plant).stm = (mass->c3).stm*loct->c3ptn[f] + (mass->c4).stm*loct->c4ptn[f];
	(mass->plant).rot = (mass->c3).rot*loct->c3ptn[f] + (mass->c4).rot*loct->c4ptn[f];
	(mass->plant).mfol[f] = (mass->c3).mfol[f]*loct->c3ptn[f] + (mass->c4).mfol[f]*loct->c4ptn[f];
	(mass->plant).mstm[f] = (mass->c3).mstm[f]*loct->c3ptn[f] + (mass->c4).mstm[f]*loct->c4ptn[f];
	(mass->plant).mrot[f] = (mass->c3).mrot[f]*loct->c3ptn[f] + (mass->c4).mrot[f]*loct->c4ptn[f];
	(mass->plant).lai[f] = (mass->c3).lai[f]*loct->c3ptn[f] + (mass->c4).lai[f]*loct->c4ptn[f];	
	
	/* nitrogen */
	(mass->plant).n_cnpy = (mass->c3).n_cnpy*loct->c3ptn[f] + (mass->c4).n_cnpy*loct->c4ptn[f];
	(mass->plant).n_strg = (mass->c3).n_strg*loct->c3ptn[f] + (mass->c4).n_strg*loct->c4ptn[f];
	(mass->plant).n_cnpy_m[f] = (mass->c3).n_cnpy_m[f]*loct->c3ptn[f] + (mass->c4).n_cnpy_m[f]*loct->c4ptn[f];
	(mass->plant).n_strg_m[f] = (mass->c3).n_strg_m[f]*loct->c3ptn[f] + (mass->c4).n_strg_m[f]*loct->c4ptn[f];

	/** flow **************************************************************************/
	/** production **/
	(flux->plant).gpp[f] = (flux->c3).gpp[f]*loct->c3ptn[f] + (flux->c4).gpp[f]*loct->c4ptn[f];
	(flux->plant).spp[f] = (flux->c3).spp[f]*loct->c3ptn[f] + (flux->c4).spp[f]*loct->c4ptn[f];
	(flux->plant).epp[f] = (flux->c3).epp[f]*loct->c3ptn[f] + (flux->c4).epp[f]*loct->c4ptn[f];
	(flux->plant).npp[f] = (flux->c3).npp[f]*loct->c3ptn[f] + (flux->c4).npp[f]*loct->c4ptn[f];
	
	(flux->plant).gpp_df97[f] = (flux->c3).gpp_df97[f]*loct->c3ptn[f] + (flux->c4).gpp_df97[f]*loct->c4ptn[f];

	/** maintenance respiration **/
	(flux->plant).rfm[f] = (flux->c3).rfm[f]*loct->c3ptn[f] + (flux->c4).rfm[f]*loct->c4ptn[f]; 
	(flux->plant).rcm[f] = (flux->c3).rcm[f]*loct->c3ptn[f] + (flux->c4).rcm[f]*loct->c4ptn[f]; 
	(flux->plant).rrm[f] = (flux->c3).rrm[f]*loct->c3ptn[f] + (flux->c4).rrm[f]*loct->c4ptn[f]; 
	(flux->plant).arm[f] = (flux->plant).rfm[f] + (flux->plant).rcm[f] + (flux->plant).rrm[f];
	
	/** growth respiration **/
	(flux->plant).rfg[f] = (flux->c3).rfg[f]*loct->c3ptn[f] + (flux->c4).rfg[f]*loct->c4ptn[f]; 
	(flux->plant).rcg[f] = (flux->c3).rcg[f]*loct->c3ptn[f] + (flux->c4).rcg[f]*loct->c4ptn[f]; 
	(flux->plant).rrg[f] = (flux->c3).rrg[f]*loct->c3ptn[f] + (flux->c4).rrg[f]*loct->c4ptn[f]; 
	(flux->plant).arg[f] = (flux->plant).rfg[f] + (flux->plant).rcg[f] + (flux->plant).rrg[f]; 
	
	/** total respiration **/
	(flux->plant).ar[f] = (flux->plant).arm[f] + (flux->plant).arg[f];			

	/** allocation **/
	(flux->plant).tpf[f] = (flux->c3).tpf[f]*loct->c3ptn[f] + (flux->c4).tpf[f]*loct->c4ptn[f]; 
	(flux->plant).tpc[f] = (flux->c3).tpc[f]*loct->c3ptn[f] + (flux->c4).tpc[f]*loct->c4ptn[f]; 
	(flux->plant).tpr[f] = (flux->c3).tpr[f]*loct->c3ptn[f] + (flux->c4).tpr[f]*loct->c4ptn[f]; 
	(flux->plant).tpp[f] = (flux->c3).tpp[f]*loct->c3ptn[f] + (flux->c4).tpp[f]*loct->c4ptn[f]; 
	
	/** litterfall **/
	(flux->plant).lf[f] = (flux->c3).lf[f]*loct->c3ptn[f] + (flux->c4).lf[f]*loct->c4ptn[f]; 
	(flux->plant).lc[f] = (flux->c3).lc[f]*loct->c3ptn[f] + (flux->c4).lc[f]*loct->c4ptn[f]; 
	(flux->plant).lr[f] = (flux->c3).lr[f]*loct->c3ptn[f] + (flux->c4).lr[f]*loct->c4ptn[f];
	(flux->plant).lL[f] = (flux->plant).lf[f] + (flux->plant).lc[f] + (flux->plant).lr[f];
	
	(flux->plant).lf_c[f] = (flux->c3).lf_c[f]*loct->c3ptn[f] + (flux->c4).lf_c[f]*loct->c4ptn[f]; 

	/** harvest **/
	(flux->plant).net_crop[f] = (flux->c3).net_crop[f]*loct->c3ptn[f] + (flux->c4).net_crop[f]*loct->c4ptn[f];
    (flux->plant).hvst_crop[f] = (flux->c3).hvst_crop[f]*loct->c3ptn[f] + (flux->c4).hvst_crop[f]*loct->c4ptn[f];

	(flux->plant).emit_ch4_kirschbaum_mass[f] = (flux->c3).emit_ch4_kirschbaum_mass[f]*loct->c3ptn[f] 
		+ (flux->c4).emit_ch4_kirschbaum_mass[f]*loct->c4ptn[f];
	(flux->plant).emit_ch4_kirschbaum_photo[f] = (flux->c3).emit_ch4_kirschbaum_photo[f]*loct->c3ptn[f] 
		+ (flux->c4).emit_ch4_kirschbaum_photo[f]*loct->c4ptn[f];

	/** stable carbon isotope **/
	/* litterfall */
	(flux->c3).d13c_lL[f] = d13c_addition3(
		(flux->c3).d13c_lf[f], (flux->c3).lf[f], (flux->c3).d13c_lc[f], (flux->c3).lc[f], (flux->c3).d13c_lr[f], (flux->c3).lr[f]);
	(flux->c4).d13c_lL[f] = d13c_addition3(
		(flux->c4).d13c_lf[f], (flux->c4).lf[f], (flux->c4).d13c_lc[f], (flux->c4).lc[f], (flux->c4).d13c_lr[f], (flux->c4).lr[f]);
	(flux->plant).d13c_lL[f] = d13c_addition(
		(flux->c3).d13c_lL[f], loct->c3ptn[f]*(flux->c3).lL[f], (flux->c4).d13c_lL[f], loct->c4ptn[f]*(flux->c4).lL[f]);
	
	(flux->plant).d13c_lf[f] = d13c_addition(
		(flux->c3).d13c_lf[f], loct->c3ptn[f]*(flux->c3).lf[f], (flux->c4).d13c_lf[f], loct->c4ptn[f]*(flux->c4).lf[f]);
	(flux->plant).d13c_lc[f] = d13c_addition(
		(flux->c3).d13c_lc[f], loct->c3ptn[f]*(flux->c3).lc[f], (flux->c4).d13c_lc[f], loct->c4ptn[f]*(flux->c4).lc[f]);
	(flux->plant).d13c_lr[f] = d13c_addition(
		(flux->c3).d13c_lr[f], loct->c3ptn[f]*(flux->c3).lr[f], (flux->c4).d13c_lr[f], loct->c4ptn[f]*(flux->c4).lr[f]);
	
	/* respiration */
	(flux->plant).d13c_rfm[f] = d13c_addition(
		(flux->c3).d13c_rfm[f], loct->c3ptn[f]*(flux->c3).rfm[f], (flux->c4).d13c_rfm[f], loct->c4ptn[f]*(flux->c4).rfm[f]);
	(flux->plant).d13c_rcm[f] = d13c_addition(
		(flux->c3).d13c_rcm[f], loct->c3ptn[f]*(flux->c3).rcm[f], (flux->c4).d13c_rcm[f], loct->c4ptn[f]*(flux->c4).rcm[f]);
	(flux->plant).d13c_rrm[f] = d13c_addition(
		(flux->c3).d13c_rrm[f], loct->c3ptn[f]*(flux->c3).rrm[f], (flux->c4).d13c_rrm[f], loct->c4ptn[f]*(flux->c4).rrm[f]);
	(flux->plant).d13c_arm[f] = d13c_addition3(
		(flux->plant).d13c_rfm[f], (flux->plant).rfm[f], (flux->plant).d13c_rcm[f], (flux->plant).rcm[f], (flux->plant).d13c_rrm[f], (flux->plant).rrm[f]);

	(flux->plant).d13c_rfg[f] = d13c_addition(
		(flux->c3).d13c_rfg[f], loct->c3ptn[f]*(flux->c3).rfg[f], (flux->c4).d13c_rfg[f], loct->c4ptn[f]*(flux->c4).rfg[f]);
	(flux->plant).d13c_rcg[f] = d13c_addition(
		(flux->c3).d13c_rcg[f], loct->c3ptn[f]*(flux->c3).rcg[f], (flux->c4).d13c_rcg[f], loct->c4ptn[f]*(flux->c4).rcg[f]);
	(flux->plant).d13c_rrg[f] = d13c_addition(
		(flux->c3).d13c_rrg[f], loct->c3ptn[f]*(flux->c3).rrg[f], (flux->c4).d13c_rrg[f], loct->c4ptn[f]*(flux->c4).rrg[f]);
	(flux->plant).d13c_arg[f] = d13c_addition3(
		(flux->plant).d13c_rfg[f], (flux->plant).rfg[f], (flux->plant).d13c_rcg[f], (flux->plant).rcg[f], (flux->plant).d13c_rrg[f], (flux->plant).rrg[f]);

	(flux->c3).d13c_arm[f] = d13c_addition3(
		(flux->c3).d13c_rfm[f], (flux->c3).rfm[f], (flux->c3).d13c_rcm[f], (flux->c3).rcm[f], (flux->c3).d13c_rrm[f], (flux->c3).rrm[f]);
	(flux->c4).d13c_arm[f] = d13c_addition3(
		(flux->c4).d13c_rfm[f], (flux->c4).rfm[f], (flux->c4).d13c_rcm[f], (flux->c4).rcm[f], (flux->c4).d13c_rrm[f], (flux->c4).rrm[f]);
	(flux->c3).d13c_arg[f] = d13c_addition3(
		(flux->c3).d13c_rfg[f], (flux->c3).rfg[f], (flux->c3).d13c_rcg[f], (flux->c3).rcg[f], (flux->c3).d13c_rrg[f], (flux->c3).rrg[f]);
	(flux->c4).d13c_arg[f] = d13c_addition3(
		(flux->c4).d13c_rfg[f], (flux->c4).rfg[f], (flux->c4).d13c_rcg[f], (flux->c4).rcg[f], (flux->c4).d13c_rrg[f], (flux->c4).rrg[f]);

	(flux->c3).d13c_ar[f] = d13c_addition(
		(flux->c3).d13c_arm[f], (flux->c3).arm[f], (flux->c3).d13c_arg[f], (flux->c3).arg[f]);
	(flux->c4).d13c_ar[f] = d13c_addition(
		(flux->c4).d13c_arm[f], (flux->c4).arm[f], (flux->c4).d13c_arg[f], (flux->c4).arg[f]);
	(flux->plant).d13c_ar[f] = d13c_addition(
		(flux->plant).d13c_arm[f],(flux->plant).arm[f], (flux->plant).d13c_arg[f], (flux->plant).arg[f]);
	
	/* d14C */
	(flux->plant).d14c_lL[f] = (loct->c3ptn[f]*(flux->c3).d14c_lL[f]*(flux->c3).lL[f] + loct->c4ptn[f]*(flux->c4).d14c_lL[f]*(flux->c4).lL[f]) 
						/ (loct->c3ptn[f]*(flux->c3).lL[f] + loct->c4ptn[f]*(flux->c4).lL[f]);
	
	/* nitrogen */
	(flux->plant).n_abdn_cnpy[f] = (flux->c3).n_abdn_cnpy[f]*loct->c3ptn[f] + (flux->c4).n_abdn_cnpy[f]*loct->c4ptn[f];
	(flux->plant).n_abdn_strg[f] = (flux->c3).n_abdn_strg[f]*loct->c3ptn[f] + (flux->c4).n_abdn_strg[f]*loct->c4ptn[f];
	(flux->plant).uptake_no3[f] = (flux->c3).uptake_no3[f]*loct->c3ptn[f] + (flux->c4).uptake_no3[f]*loct->c4ptn[f];
	(flux->plant).uptake_nh4[f] = (flux->c3).uptake_nh4[f]*loct->c3ptn[f] + (flux->c4).uptake_nh4[f]*loct->c4ptn[f];
	(flux->plant).n_biofix[f] = (flux->c3).n_biofix[f]*loct->c3ptn[f] + (flux->c4).n_biofix[f]*loct->c4ptn[f];

	(flux->plant).n_alloc_cnpy[f] = (flux->c3).n_alloc_cnpy[f]*loct->c3ptn[f] + (flux->c4).n_alloc_cnpy[f]*loct->c4ptn[f];
	(flux->plant).n_alloc_strg[f] = (flux->c3).n_alloc_strg[f]*loct->c3ptn[f] + (flux->c4).n_alloc_strg[f]*loct->c4ptn[f];
	(flux->plant).n_realloc[f] = (flux->c3).n_realloc[f]*loct->c3ptn[f] + (flux->c4).n_realloc[f]*loct->c4ptn[f];
	(flux->plant).n_salvage[f] = (flux->c3).n_salvage[f]*loct->c3ptn[f] + (flux->c4).n_salvage[f]*loct->c4ptn[f];
}

/* dealings before calculating carbon budget *******************************/
void f_before_deal(
	struct Grid *grid, 
	struct Pflx *flux
){
	plant_flux_zero(grid->m, flux);
}

/* dealings after calculating carbon budget *******************************/
void f_after_deal(
	struct Grid *grid,
    struct Loct *loct,
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){
    double aaa = 0.0;
	double in_d14c;
    
    /* constrain max LAI: 2015/04/06 */
    if(CONSTRAIN_LAIMAX == 1 && (grid->veg_olson >=1 && grid->veg_olson <= 30)){
        
        if(mass->lai[grid->m] > loct->est_maxlai && mass->lai[grid->m] > 1.0){
            aaa = (mass->lai[grid->m] - loct->est_maxlai) *100.0*2.0/2.2/pchar->sla;
            
            mass->fol -= aaa * 0.1;
            flux->lf[grid->m] += aaa * 0.1;
        }else{
            ;
        }
    }
	
	/* LAI update */
	mass->lai[grid->m] = lai_mass(grid, mass, pchar);
	
	/* total respiration */
	flux->arg[grid->m] = flux->rfg[grid->m] + flux->rcg[grid->m] + flux->rrg[grid->m];
	flux->ar[grid->m] = flux->arm[grid->m] + flux->arg[grid->m];
	
	if(DF97==1){
		flux->spp[grid->m] = flux->gpp_df97[grid->m] - flux->rfm[grid->m]-flux->rfg[grid->m];
		flux->npp[grid->m] = flux->gpp_df97[grid->m] - flux->arm[grid->m]-flux->arg[grid->m];
	}else{
		flux->spp[grid->m] = flux->gpp[grid->m] - flux->rfm[grid->m]-flux->rfg[grid->m];
		flux->npp[grid->m] = flux->gpp[grid->m] - flux->arm[grid->m]-flux->arg[grid->m];
	}
	
	/* total litterfall */
	flux->lL[grid->m] = flux->lf[grid->m] + flux->lc[grid->m] + flux->lr[grid->m];
		
	/* monthly mass ********************************************************/	
	mass->mfol[grid->m] = mass->fol;
	mass->mstm[grid->m] = mass->stm;
	mass->mrot[grid->m] = mass->rot;
	mass->plant[grid->m] = mass->fol + mass->stm + mass->rot;
	
	mass->fol_p = mass->fol;

	/* monthly d13c of mass */
	mass->d13c_mfol[grid->m] = mass->d13c_fol;
	mass->d13c_mstm[grid->m] = mass->d13c_stm;
	mass->d13c_mrot[grid->m] = mass->d13c_rot;
	
	/* d14C: added by A.Ito (2009/06/27) *****************************************/
	/* assumption: double discrimination */
	in_d14c = grid->d14c_bco2[grid->m] - 2.0*pchar->photo_13c_frac[grid->m];
	
	flux->d14c_gpp[grid->m] = in_d14c;
	
	/* litter */
	if((flux->lf[grid->m] + flux->lc[grid->m] + flux->lr[grid->m]) > 0.0){
		flux->d14c_lL[grid->m] = (mass->d14c_fol*flux->lf[grid->m] + mass->d14c_stm*flux->lc[grid->m] + mass->d14c_rot*flux->lr[grid->m]) / 
				(flux->lf[grid->m] + flux->lc[grid->m] + flux->lr[grid->m]);
	}else{
		/* no litterfall */
		flux->d14c_lL[grid->m] = 0.0;
	}
	
	/* leaves */
	if((flux->tpf[grid->m]-flux->rfg[grid->m]) > 0.0){
		mass->d14c_fol = (mass->d14c_fol*mass->fol + in_d14c*flux->tpf[grid->m]) 
						/ (mass->fol + flux->tpf[grid->m]);
	}else if((flux->tpf[grid->m]-flux->rfg[grid->m]) <= 0.0){
		mass->d14c_fol = (mass->d14c_fol*mass->fol + in_d14c*pchar->malloc_f[grid->m]*flux->gpp[grid->m]) 
						  / (mass->fol + pchar->malloc_f[grid->m]*flux->gpp[grid->m]);
	}
	if(mass->fol < 0.025){
		mass->d14c_fol = grid->d14c_bco2[grid->m];
	}
	mass->d14c_mfol[grid->m] = mass->d14c_fol;
	
	/* stems */
	if((flux->tpc[grid->m]-flux->rcg[grid->m]) > 0.0){
		mass->d14c_stm = (mass->d14c_stm*mass->stm + in_d14c*flux->tpc[grid->m]) 
						/ (mass->stm + flux->tpc[grid->m]);
	}else if((flux->tpc[grid->m]-flux->rcg[grid->m]) <= 0.0){
		mass->d14c_stm = (mass->d14c_stm*mass->stm +  in_d14c* pchar->malloc_c[grid->m]*flux->gpp[grid->m])
						 / (mass->stm + pchar->malloc_c[grid->m]*flux->gpp[grid->m]);
	}
	if(mass->stm < 0.025){
		mass->d14c_stm = grid->d14c_bco2[grid->m];
	}
	mass->d14c_mstm[grid->m] = mass->d14c_stm;
	
	/* roots */
	if((flux->tpr[grid->m]-flux->rrg[grid->m]) > 0.0){
		mass->d14c_rot = (mass->d14c_rot*mass->rot + in_d14c*flux->tpr[grid->m]) 
						/ (mass->rot + flux->tpr[grid->m]);
	}else if((flux->tpr[grid->m]-flux->rrg[grid->m]) <= 0.0){
		mass->d14c_rot = (mass->d14c_rot*mass->rot + in_d14c*pchar->malloc_r[grid->m]*flux->gpp[grid->m]) 
						  / (mass->rot + pchar->malloc_r[grid->m]*flux->gpp[grid->m]);
	}
	if(mass->rot < 0.025){
		mass->d14c_rot = grid->d14c_bco2[grid->m];
	}
	mass->d14c_mrot[grid->m] = mass->d14c_rot;
}

/**********************************************************/
/* area (ha) of grid lat1-lat2 x lon1-lon2 (in degree) */
double grid_area(
	double lat1, 
	double lat2, 
	double lon1, 
	double lon2
){
	double area;
	double e_rad = 6378137.0; /* m, GRS-80(revised) */
	double e_exc, e_flat = 298.257;
	double l_lat, l_lon1, l_lon2;
	double m_lat, aa;
	
	if(lat1>90.0) lat1 = 90.0;
	if(lat2<-90.0) lat2 = -90.0;
	
	m_lat = (lat1+lat2)/2.0 * PI/180.0;
	e_exc = sqrt(2.0/e_flat - 1.0/(e_flat*e_flat));
	
	aa = 1.0 - e_exc*e_exc*sin(m_lat)*sin(m_lat);
	l_lat = PI/180.0*e_rad*(1.0 - e_exc*e_exc)/pow(aa, 1.5) * fabs(lat1 - lat2);
	
	/* bug fixed 2007/12/25 E.Kato and A.Ito */
	aa = 1.0 - e_exc * e_exc * sin(lat1*PI/180.0)*sin(lat1*PI/180.0);
	l_lon1 = PI/180.0 * e_rad * cos(lat1*PI/180.0)/sqrt(aa) * fabs(lon1 - lon2);
	aa = 1.0 - e_exc * e_exc * sin(lat2*PI/180.0)*sin(lat2*PI/180.0);
	l_lon2 = PI/180.0 * e_rad * cos(lat2*PI/180.0)/sqrt(aa) * fabs(lon1 - lon2);
	
	area = (l_lon1 + l_lon2)*l_lat/2.0 / 10000.0;
	
	return (area); /* ha */
}

/******************************************************/
long f_basin_id_trip(
	long original
){
	long basin;

	switch(original){
		case 1: basin = 1; break;
		case 4: basin = 2; break;
		case 2: basin = 3; break;
		case 3: basin = 4; break;
		case 10: basin = 5; break;
		case 5: basin = 6; break;
		case 7: basin = 7; break;
		case 6: basin = 8; break;
		case 8: basin = 9; break;
		case 9: basin = 10; break;
		case 14: basin = 11; break;
		case 11: basin = 12; break;
		case 12: basin = 13; break;
		case 21: basin = 14; break;
		case 34: basin = 14; break;
		case 13: basin = 15; break;
		case 16: basin = 16; break;
		case 20: basin = 17; break;
		case 18: basin = 18; break;
		case 68: basin = 18; break;
		case 23: basin = 19; break;
		case 17: basin = 20; break;
		case 31: basin = 21; break;
		case 19: basin = 22; break;
		case 33: basin = 23; break;
		case 26: basin = 24; break;
		case 45: basin = 25; break;
		case 32: basin = 26; break;
		case 38: basin = 27; break;
		case 41: basin = 28; break;
		case 28: basin = 29; break;
		case 29: basin = 30; break;
		case 24: basin = 31; break;
		case 36: basin = 32; break;
		case 39: basin = 33; break;
		case 40: basin = 34; break;
		case 42: basin = 35; break;
		case 51: basin = 36; break;
		case 296: basin = 37; break;
		case 47: basin = 38; break;
		case 48: basin = 39; break;
		case 54: basin = 40; break;
		case 56: basin = 41; break;
		case 53: basin = 42; break;
		case 57: basin = 43; break;
		case 110: basin = 44; break;
		case 64: basin = 45; break;
		case 78: basin = 46; break;
		case 62: basin = 47; break;
		case 74: basin = 48; break;
		case 81: basin = 49; break;
		case 79: basin = 50; break;
		case 83: basin = 51; break;
		case 102: basin = 52; break;
		case 378: basin = 53; break;
		
		default:  basin=0;
	}
	
	return basin;
}

/******************************************************************************/
void set_rowcol_gcm(
	void
){
	if(SCENARIO_ID == 0){
		GCM_R = 1;
		GCM_C = 1;
	}else if(SCENARIO_ID==1 || SCENARIO_ID==2 || SCENARIO_ID==3 ||
                SCENARIO_ID==4 || SCENARIO_ID==5 || SCENARIO_ID==6){
		GCM_R = 32;
		GCM_C = 64;
	}else if(SCENARIO_ID==11 || SCENARIO_ID==12 || SCENARIO_ID==13 ||
                SCENARIO_ID==14 || SCENARIO_ID==15 || SCENARIO_ID==16 ||
                SCENARIO_ID==17 || SCENARIO_ID==18){
		GCM_R = 48;
		GCM_C = 96;
	}else if(SCENARIO_ID==21 || SCENARIO_ID==22 || SCENARIO_ID==23 || SCENARIO_ID==24){
		GCM_R = 73;
		GCM_C = 96;
	}else if(SCENARIO_ID==31 || SCENARIO_ID==32){
		GCM_R = 64;
		GCM_C = 128;
	}else if(SCENARIO_ID==41 || SCENARIO_ID==42 || SCENARIO_ID==43 || SCENARIO_ID==44){
		GCM_R = 56;
		GCM_C = 64;
	}else if(SCENARIO_ID==51 || SCENARIO_ID==52){
		GCM_R = 80;
		GCM_C = 96;
	}else if(SCENARIO_ID==61 || SCENARIO_ID==62 || SCENARIO_ID==63){
		GCM_R = 64;
		GCM_C = 128;
	}else if(SCENARIO_ID==71){
		GCM_R = 64;
		GCM_C = 128;
	}else if(SCENARIO_ID==1000 || SCENARIO_ID==1001){ /* MIROC-HIGH */
		GCM_R = 160;
		GCM_C = 320;
	}else if(SCENARIO_ID==1010 || SCENARIO_ID==1011 || SCENARIO_ID==1012 ||
            SCENARIO_ID==1013 || SCENARIO_ID==1014 || SCENARIO_ID==1015 ||
            SCENARIO_ID==1016 || SCENARIO_ID==1017 || SCENARIO_ID==1018){ /* MIROC-MED */
		GCM_R = 64;
		GCM_C = 128;
	}else if(SCENARIO_ID==1050 || SCENARIO_ID==1051 ||SCENARIO_ID==1052){ /* BCCR */
		GCM_R = 64;
		GCM_C = 128;
	}else if(SCENARIO_ID==1060 || SCENARIO_ID==1061 ||SCENARIO_ID==1062){ /* INM-CM3 : error-fixed 2009/04/26 A.Ito */
		GCM_R = 45;
		GCM_C = 72;
	}else if(SCENARIO_ID==1070 || SCENARIO_ID==1071 || SCENARIO_ID==1072){ /* IPSL */
		GCM_R = 72;
		GCM_C = 96;
	}else if(SCENARIO_ID==1080 || SCENARIO_ID==1081 || SCENARIO_ID==1082){ /* GFDL CM2-1 */
		GCM_R = 90;
		GCM_C = 144;
	}else if(SCENARIO_ID==1090 || SCENARIO_ID==1091 || SCENARIO_ID==1092){ /* HadCM3 */
		GCM_R = 73;
		GCM_C = 96;
	}else if(SCENARIO_ID==1100 || SCENARIO_ID==1101 || SCENARIO_ID==1102 || SCENARIO_ID==1103){ /* GISS AOM */
		GCM_R = 60;
		GCM_C = 90;
	}else if(SCENARIO_ID==1110 || SCENARIO_ID==1111 || SCENARIO_ID==1112){ /* CCCma-T63 */
		GCM_R = 64;
		GCM_C = 128;
	}else if(SCENARIO_ID==1120 || SCENARIO_ID==1121 || SCENARIO_ID==1122 || SCENARIO_ID==1270 || SCENARIO_ID==1271 || SCENARIO_ID==1272){ /* CSIRO */
		GCM_R = 96;
		GCM_C = 192;
	}else if(SCENARIO_ID==1130 || SCENARIO_ID==1131 || SCENARIO_ID==1132 || SCENARIO_ID==1133 || SCENARIO_ID==1134
			 || SCENARIO_ID==1135 || SCENARIO_ID==1136 || SCENARIO_ID==1137 || SCENARIO_ID==1138 || SCENARIO_ID==1139
			  || SCENARIO_ID==1140 || SCENARIO_ID==1141 || SCENARIO_ID==1142 || SCENARIO_ID==1143 || SCENARIO_ID==1144){ /* MRI */
		GCM_R = 64;
		GCM_C = 128;
	}else if(SCENARIO_ID==1150 ||SCENARIO_ID==1151 ||SCENARIO_ID==1152 ||SCENARIO_ID==1153 ||SCENARIO_ID==1154 ||
			SCENARIO_ID==1155 ||SCENARIO_ID==1156 ||SCENARIO_ID==1157 ||SCENARIO_ID==1158 ||SCENARIO_ID==1159){ /* MPI ECHAM */
		GCM_R = 96;
		GCM_C = 192;
	}else if(SCENARIO_ID==1160 || SCENARIO_ID==1161 || SCENARIO_ID==1162 || SCENARIO_ID==1163 || SCENARIO_ID==1164 || SCENARIO_ID==1165){ /* IAP */
		GCM_R = 60;
		GCM_C = 128;
	}else if(SCENARIO_ID==1170 || SCENARIO_ID==1171){ /* HadGEM */
		GCM_R = 145;
		GCM_C = 192;
	}else if(SCENARIO_ID==1180 || SCENARIO_ID==1181 || SCENARIO_ID==1182){ /* GFDL CM 2.0 */
		GCM_R = 90;
		GCM_C = 144;
	}else if(SCENARIO_ID==1190 || SCENARIO_ID==1191 || SCENARIO_ID==1192){ /* CNRM */
		GCM_R = 64;
		GCM_C = 128;
	}else if(SCENARIO_ID==1200 || SCENARIO_ID==1201 || SCENARIO_ID==1202 || SCENARIO_ID==1203 || SCENARIO_ID==1204
			 || SCENARIO_ID==1205 || SCENARIO_ID==1206 || SCENARIO_ID==1207 || SCENARIO_ID==1208 || SCENARIO_ID==1209
			 || SCENARIO_ID==1210 || SCENARIO_ID==1211 || SCENARIO_ID==1212 || SCENARIO_ID==1213 || SCENARIO_ID==1214){ /* CCC47 */
		GCM_R = 48;
		GCM_C = 96;
	}else if(SCENARIO_ID==1220 || SCENARIO_ID==1221 || SCENARIO_ID==1222 || SCENARIO_ID==1223 || SCENARIO_ID==1224
			 || SCENARIO_ID==1225 || SCENARIO_ID==1226 || SCENARIO_ID==1227 || SCENARIO_ID==1228 || SCENARIO_ID==1229
			  || SCENARIO_ID==1230 || SCENARIO_ID==1231 || SCENARIO_ID==1232 || SCENARIO_ID==1233 || SCENARIO_ID==1234
			   || SCENARIO_ID==1235 || SCENARIO_ID==1236 || SCENARIO_ID==1237 || SCENARIO_ID==1238 || SCENARIO_ID==1239){ /* NCAR CCSM */
		GCM_R = 128;
		GCM_C = 256;
	}else if(SCENARIO_ID==1240 || SCENARIO_ID==1241 || SCENARIO_ID==1242){ /* GISS E */
		GCM_R = 46;
		GCM_C = 72;
	}else if(SCENARIO_ID==1250 || SCENARIO_ID==1251 || SCENARIO_ID==1252 || SCENARIO_ID==1253){ /* GISS R */
		GCM_R = 46;
		GCM_C = 72;
	}else if(SCENARIO_ID==1260 || SCENARIO_ID==1261 || SCENARIO_ID==1262 || SCENARIO_ID==1263 || SCENARIO_ID==1264
			 || SCENARIO_ID==1265 || SCENARIO_ID==1266 || SCENARIO_ID==1267){ /* NCAR PCM */
		GCM_R = 64;
		GCM_C = 128;
	}else if(SCENARIO_ID==3000 || SCENARIO_ID==3001 || SCENARIO_ID==3002){
		GCM_R = 64;
		GCM_C = 128;
	}else if(SCENARIO_ID==3010 || SCENARIO_ID==3011 || SCENARIO_ID==3012){
		GCM_R = 56;
		GCM_C = 64;
	}else if(SCENARIO_ID==3020 || SCENARIO_ID==3021 || SCENARIO_ID==3022){
		GCM_R = 90;
		GCM_C = 144;
	}else if(SCENARIO_ID==3030 || SCENARIO_ID==3031 || SCENARIO_ID==3032 || SCENARIO_ID==3033 || SCENARIO_ID==3034 || SCENARIO_ID==3035){
		GCM_R = 145;
		GCM_C = 192;
	}else if(SCENARIO_ID==3040 || SCENARIO_ID==3041 || SCENARIO_ID==3042){
		GCM_R = 96;
		GCM_C = 96;
	}else if(SCENARIO_ID==3050 || SCENARIO_ID==3051 || SCENARIO_ID==3052){
		GCM_R = 64;
		GCM_C = 128;
	}else if(SCENARIO_ID==3060 || SCENARIO_ID==3061){
		GCM_R = 64;
		GCM_C = 128;
	}else if(SCENARIO_ID==3070 || SCENARIO_ID==3071 || SCENARIO_ID==3072){
		GCM_R = 64;
		GCM_C = 128;
	}else if(SCENARIO_ID==3080 || SCENARIO_ID==3081){
		GCM_R = 96;
		GCM_C = 192;
	}else if(SCENARIO_ID==3090 || SCENARIO_ID==3091){
		GCM_R = 192;
		GCM_C = 288;
	}else if(SCENARIO_ID==3100 || SCENARIO_ID==3101){
		GCM_R = 96;
		GCM_C = 144;
	}else{
        /* ISI-MIP, PLUME (no stored GCM data) */
		GCM_R = 1;
		GCM_C = 1;
	}
}

/****************************************************************/
void set_gcm_index(
	char s_case[]
){
	switch(SCENARIO_ID){
		case 1000:	strcpy(s_case,"MHA11_"); break;
		case 1001:	strcpy(s_case,"MHB11_"); break;

		case 1010:	strcpy(s_case,"MMA11_"); break;
		case 1011:	strcpy(s_case,"MMA12_"); break;
		case 1012:	strcpy(s_case,"MMA13_"); break;
		case 1013:	strcpy(s_case,"MMA21_"); break;
		case 1014:	strcpy(s_case,"MMA22_"); break;
		case 1015:	strcpy(s_case,"MMA23_"); break;
		case 1016:	strcpy(s_case,"MMB11_"); break;
		case 1017:	strcpy(s_case,"MMB12_"); break;
		case 1018:	strcpy(s_case,"MMB13_"); break;

		case 1050:	strcpy(s_case,"BCA11_"); break;
		case 1051:	strcpy(s_case,"BCA21_"); break;
		case 1052:	strcpy(s_case,"BCB11_"); break;

		case 1060:	strcpy(s_case,"INA11_"); break;
		case 1061:	strcpy(s_case,"INA21_"); break;
		case 1062:	strcpy(s_case,"INB11_"); break;

		case 1070:	strcpy(s_case,"IPA11_"); break;
		case 1071:	strcpy(s_case,"IPA21_"); break;
		case 1072:	strcpy(s_case,"IPB11_"); break;

		case 1080:	strcpy(s_case,"GF1A11_"); break;
		case 1081:	strcpy(s_case,"GF1A21_"); break;
		case 1082:	strcpy(s_case,"GF1B11_"); break;

		case 1090:	strcpy(s_case,"HCA11_"); break;
		case 1091:	strcpy(s_case,"HCA21_"); break;
		case 1092:	strcpy(s_case,"HCB11_"); break;
			
		case 1100:	strcpy(s_case,"GSAA11_"); break;
		case 1101:	strcpy(s_case,"GSAA12_"); break;
		case 1102:	strcpy(s_case,"GSAB11_"); break; /* corercted 2009/04/08 */
		case 1103:	strcpy(s_case,"GSAB12_"); break; /* corercted 2009/04/08 */

		case 1110:	strcpy(s_case,"CC6A11_"); break;
		case 1111:	strcpy(s_case,"CC6A21_"); break;
		case 1112:	strcpy(s_case,"CC6B11_"); break;

		case 1120:	strcpy(s_case,"CSA11_"); break;
		case 1121:	strcpy(s_case,"CSA21_"); break;
		case 1122:	strcpy(s_case,"CSB11_"); break;

		case 1130:	strcpy(s_case,"MRA11_"); break;
		case 1131:	strcpy(s_case,"MRA12_"); break;
		case 1132:	strcpy(s_case,"MRA13_"); break;
		case 1133:	strcpy(s_case,"MRA14_"); break;
		case 1134:	strcpy(s_case,"MRA15_"); break;
		case 1135:	strcpy(s_case,"MRA21_"); break;
		case 1136:	strcpy(s_case,"MRA22_"); break;
		case 1137:	strcpy(s_case,"MRA23_"); break;
		case 1138:	strcpy(s_case,"MRA24_"); break;
		case 1139:	strcpy(s_case,"MRA25_"); break;
		case 1140:	strcpy(s_case,"MRB11_"); break;
		case 1141:	strcpy(s_case,"MRB12_"); break;
		case 1142:	strcpy(s_case,"MRB13_"); break;
		case 1143:	strcpy(s_case,"MRB14_"); break;
		case 1144:	strcpy(s_case,"MRB15_"); break;

		case 1150:	strcpy(s_case,"MPA11_"); break;
		case 1151:	strcpy(s_case,"MPA12_"); break;
		case 1152:	strcpy(s_case,"MPA13_"); break;
		case 1153:	strcpy(s_case,"MPA14_"); break;
		case 1154:	strcpy(s_case,"MPA21_"); break;
		case 1155:	strcpy(s_case,"MPA22_"); break;
		case 1156:	strcpy(s_case,"MPA23_"); break;
		case 1157:	strcpy(s_case,"MPB11_"); break;
		case 1158:	strcpy(s_case,"MPB12_"); break;
		case 1159:	strcpy(s_case,"MPB13_"); break;

		case 1160:	strcpy(s_case,"IAA11_"); break;
		case 1161:	strcpy(s_case,"IAA12_"); break;
		case 1162:	strcpy(s_case,"IAA13_"); break;
		case 1163:	strcpy(s_case,"IAB11_"); break;
		case 1164:	strcpy(s_case,"IAB12_"); break;
		case 1165:	strcpy(s_case,"IAB13_"); break;

		case 1170:	strcpy(s_case,"HGA11_"); break;
		case 1171:	strcpy(s_case,"HGA21_"); break;

		case 1180:	strcpy(s_case,"GF0A11_"); break;
		case 1181:	strcpy(s_case,"GF0A21_"); break;
		case 1182:	strcpy(s_case,"GF0B11_"); break;

		case 1190:	strcpy(s_case,"CNA11_"); break;
		case 1191:	strcpy(s_case,"CNA21_"); break;
		case 1192:	strcpy(s_case,"CNB11_"); break;

		case 1200:	strcpy(s_case,"CC4A11_"); break;
		case 1201:	strcpy(s_case,"CC4A12_"); break;
		case 1202:	strcpy(s_case,"CC4A13_"); break;
		case 1203:	strcpy(s_case,"CC4A14_"); break;
		case 1204:	strcpy(s_case,"CC4A15_"); break;
		case 1205:	strcpy(s_case,"CC4A21_"); break;
		case 1206:	strcpy(s_case,"CC4A22_"); break;
		case 1207:	strcpy(s_case,"CC4A23_"); break;
		case 1208:	strcpy(s_case,"CC4A24_"); break;
		case 1209:	strcpy(s_case,"CC4A25_"); break;
		case 1210:	strcpy(s_case,"CC4B11_"); break;
		case 1211:	strcpy(s_case,"CC4B12_"); break;
		case 1212:	strcpy(s_case,"CC4B13_"); break;
		case 1213:	strcpy(s_case,"CC4B14_"); break;
		case 1214:	strcpy(s_case,"CC4B15_"); break;

		case 1220:	strcpy(s_case,"NCA11_"); break;
		case 1221:	strcpy(s_case,"NCA12_"); break;
		case 1222:	strcpy(s_case,"NCA13_"); break;
		case 1223:	strcpy(s_case,"NCA15_"); break;
		case 1224:	strcpy(s_case,"NCA16_"); break;
		case 1225:	strcpy(s_case,"NCA17_"); break;
		case 1226:	strcpy(s_case,"NCA19_"); break;
		case 1227:	strcpy(s_case,"NCA21_"); break;
		case 1228:	strcpy(s_case,"NCA22_"); break;
		case 1229:	strcpy(s_case,"NCA23_"); break;
		case 1230:	strcpy(s_case,"NCA24_"); break;
		case 1231:	strcpy(s_case,"NCA25_"); break;
		case 1232:	strcpy(s_case,"NCB11_"); break;
		case 1233:	strcpy(s_case,"NCB12_"); break;
		case 1234:	strcpy(s_case,"NCB13_"); break;
		case 1235:	strcpy(s_case,"NCB14_"); break;
		case 1236:	strcpy(s_case,"NCB15_"); break;
		case 1237:	strcpy(s_case,"NCB16_"); break;
		case 1238:	strcpy(s_case,"NCB17_"); break;
		case 1239:	strcpy(s_case,"NCB19_"); break;

		case 1240:	strcpy(s_case,"GSEA11_"); break;
		case 1241:	strcpy(s_case,"GSEA21_"); break;
		case 1242:	strcpy(s_case,"GSEB11_"); break;

		case 1250:	strcpy(s_case,"GSRA11_"); break;
		case 1251:	strcpy(s_case,"GSRA12_"); break;
		case 1252:	strcpy(s_case,"GSRA21_"); break;
		case 1253:	strcpy(s_case,"GSRB11_"); break;
		
		/*  corrected 2009/04/07 by A.Ito  */
		case 1260:	strcpy(s_case,"NPA11_"); break;
		case 1261:	strcpy(s_case,"NPA12_"); break;
		case 1262:	strcpy(s_case,"NPA13_"); break;
		case 1263:	strcpy(s_case,"NPA14_"); break;
		case 1264:	strcpy(s_case,"NPA21_"); break;
		case 1265:	strcpy(s_case,"NPA22_"); break;
		case 1266:	strcpy(s_case,"NPA23_"); break;
		case 1267:	strcpy(s_case,"NPA24_"); break;

		case 1270:	strcpy(s_case,"CS35A11_"); break;	/* added 2009/04/09 by A.Ito */
		case 1271:	strcpy(s_case,"CS35A21_"); break;
		case 1272:	strcpy(s_case,"CS35B11_"); break;

		case 2001:	strcpy(s_case,"ISIH26C_"); break;
		case 2002:	strcpy(s_case,"ISIH85C_"); break;
		case 2003:	strcpy(s_case,"ISIH45C_"); break;
		case 2004:	strcpy(s_case,"ISIH60C_"); break;
		case 2005:	strcpy(s_case,"ISIH26N_"); break;
		case 2006:	strcpy(s_case,"ISIH85N_"); break;
		case 2007:	strcpy(s_case,"ISIH45N_"); break;
		case 2008:	strcpy(s_case,"ISIH60N_"); break;
        
		case 2011:	strcpy(s_case,"ISII26C_"); break;
		case 2012:	strcpy(s_case,"ISII85C_"); break;
		case 2013:	strcpy(s_case,"ISII45C_"); break;
		case 2014:	strcpy(s_case,"ISII60C_"); break;
		case 2015:	strcpy(s_case,"ISII26N_"); break;
		case 2016:	strcpy(s_case,"ISII85N_"); break;
		case 2017:	strcpy(s_case,"ISII45N_"); break;
		case 2018:	strcpy(s_case,"ISII60N_"); break;
        
		case 2021:	strcpy(s_case,"ISIG26C_"); break;
		case 2022:	strcpy(s_case,"ISIG85C_"); break;
		case 2023:	strcpy(s_case,"ISIG45C_"); break;
		case 2024:	strcpy(s_case,"ISIG60C_"); break;
		case 2025:	strcpy(s_case,"ISIG26N_"); break;
		case 2026:	strcpy(s_case,"ISIG85N_"); break;
		case 2027:	strcpy(s_case,"ISIG45N_"); break;
		case 2028:	strcpy(s_case,"ISIG60N_"); break;
        
		case 2031:	strcpy(s_case,"ISIM26C_"); break;
		case 2032:	strcpy(s_case,"ISIM85C_"); break;
		case 2033:	strcpy(s_case,"ISIM45C_"); break;
		case 2034:	strcpy(s_case,"ISIM60C_"); break;
		case 2035:	strcpy(s_case,"ISIM26N_"); break;
		case 2036:	strcpy(s_case,"ISIM85N_"); break;
		case 2037:	strcpy(s_case,"ISIM45N_"); break;
		case 2038:	strcpy(s_case,"ISIM60N_"); break;
        
		case 2041:	strcpy(s_case,"ISIN26C_"); break;
		case 2042:	strcpy(s_case,"ISIN85C_"); break;
		case 2043:	strcpy(s_case,"ISIN45C_"); break;
		case 2044:	strcpy(s_case,"ISIN60C_"); break;
		case 2045:	strcpy(s_case,"ISIN26N_"); break;
		case 2046:	strcpy(s_case,"ISIN85N_"); break;
		case 2047:	strcpy(s_case,"ISIN45N_"); break;
		case 2048:	strcpy(s_case,"ISIN60N_"); break;
        
        /* ICARUS: 2016/08/15 by A.Ito */
        case 2201:	strcpy(s_case,"ICGFS1R26_"); break;
        case 2202:	strcpy(s_case,"ICGFS1R45_"); break;
        case 2203:	strcpy(s_case,"ICGFS2R26_"); break;
        case 2204:	strcpy(s_case,"ICGFS2R45_"); break;
        case 2205:	strcpy(s_case,"ICGFS2R60_"); break;
        case 2206:	strcpy(s_case,"ICGFS3R45_"); break;
        case 2207:	strcpy(s_case,"ICGFS3R60_"); break;
        case 2208:	strcpy(s_case,"ICGFS4R26_"); break;
        case 2209:	strcpy(s_case,"ICGFS4R45_"); break;
        case 2210:	strcpy(s_case,"ICGFS5R26_"); break;
        case 2211:	strcpy(s_case,"ICGFS5R45_"); break;
        
        case 2221:	strcpy(s_case,"ICHDS1R26_"); break;
        case 2222:	strcpy(s_case,"ICHDS1R45_"); break;
        case 2223:	strcpy(s_case,"ICHDS2R26_"); break;
        case 2224:	strcpy(s_case,"ICHDS2R45_"); break;
        case 2225:	strcpy(s_case,"ICHDS2R60_"); break;
        case 2226:	strcpy(s_case,"ICHDS3R45_"); break;
        case 2227:	strcpy(s_case,"ICHDS3R60_"); break;
        case 2228:	strcpy(s_case,"ICHDS4R26_"); break;
        case 2229:	strcpy(s_case,"ICHDS4R45_"); break;
        case 2230:	strcpy(s_case,"ICHDS5R26_"); break;
        case 2231:	strcpy(s_case,"ICHDS5R45_"); break;
        
        case 2241:	strcpy(s_case,"ICIPS1R26_"); break;
        case 2242:	strcpy(s_case,"ICIPS1R45_"); break;
        case 2243:	strcpy(s_case,"ICIPS2R26_"); break;
        case 2244:	strcpy(s_case,"ICIPS2R45_"); break;
        case 2245:	strcpy(s_case,"ICIPS2R60_"); break;
        case 2246:	strcpy(s_case,"ICIPS3R45_"); break;
        case 2247:	strcpy(s_case,"ICIPS3R60_"); break;
        case 2248:	strcpy(s_case,"ICIPS4R26_"); break;
        case 2249:	strcpy(s_case,"ICIPS4R45_"); break;
        case 2250:	strcpy(s_case,"ICIPS5R26_"); break;
        case 2251:	strcpy(s_case,"ICIPS5R45_"); break;
        
        case 2261:	strcpy(s_case,"ICMIS1R26_"); break;
        case 2262:	strcpy(s_case,"ICMIS1R45_"); break;
        case 2263:	strcpy(s_case,"ICMIS2R26_"); break;
        case 2264:	strcpy(s_case,"ICMIS2R45_"); break;
        case 2265:	strcpy(s_case,"ICMIS2R60_"); break;
        case 2266:	strcpy(s_case,"ICMIS3R45_"); break;
        case 2267:	strcpy(s_case,"ICMIS3R60_"); break;
        case 2268:	strcpy(s_case,"ICMIS4R26_"); break;
        case 2269:	strcpy(s_case,"ICMIS4R45_"); break;
        case 2270:	strcpy(s_case,"ICMIS5R26_"); break;
        case 2271:	strcpy(s_case,"ICMIS5R45_"); break;
 
        case 2281:	strcpy(s_case,"ICNOS1R26_"); break;
        case 2282:	strcpy(s_case,"ICNOS1R45_"); break;
        case 2283:	strcpy(s_case,"ICNOS2R26_"); break;
        case 2284:	strcpy(s_case,"ICNOS2R45_"); break;
        case 2285:	strcpy(s_case,"ICNOS2R60_"); break;
        case 2286:	strcpy(s_case,"ICNOS3R45_"); break;
        case 2287:	strcpy(s_case,"ICNOS3R60_"); break;
        case 2288:	strcpy(s_case,"ICNOS4R26_"); break;
        case 2289:	strcpy(s_case,"ICNOS4R45_"); break;
        case 2290:	strcpy(s_case,"ICNOS5R26_"); break;
        case 2291:	strcpy(s_case,"ICNOS5R45_"); break;

        case 2601:	strcpy(s_case,"ICS1_"); break;
        case 2602:	strcpy(s_case,"ICS2_"); break;
        case 2603:	strcpy(s_case,"ICS3_"); break;

        /* GEO-MIP runs: 2013/11/26 by A.Ito */
		case 3000:	strcpy(s_case,"GEBN45_"); break;
		case 3001:	strcpy(s_case,"GEBNG3_"); break;
		case 3002:	strcpy(s_case,"GEBNG4_"); break;
        
 		case 3010:	strcpy(s_case,"GECS45_"); break;
		case 3011:	strcpy(s_case,"GECSG4_"); break;
		case 3012:	strcpy(s_case,"GECSG3S_"); break;
       
		case 3020:	strcpy(s_case,"GEGI45_"); break;
		case 3021:	strcpy(s_case,"GEGIG3_"); break;
		case 3022:	strcpy(s_case,"GEGIG4_"); break;
        
		case 3030:	strcpy(s_case,"GEHD45_"); break;
		case 3031:	strcpy(s_case,"GEHDG3_"); break;
		case 3032:	strcpy(s_case,"GEHDG4_"); break;
		case 3033:	strcpy(s_case,"GEHDG3S_"); break;
		case 3034:	strcpy(s_case,"GEHDG4C_"); break;
		case 3035:	strcpy(s_case,"GEHDG4S_"); break;
        
		case 3040:	strcpy(s_case,"GEIP45_"); break;
		case 3041:	strcpy(s_case,"GEIPG3_"); break;
		case 3042:	strcpy(s_case,"GEIPG5_"); break;
        
		case 3050:	strcpy(s_case,"GEMR45_"); break;
		case 3051:	strcpy(s_case,"GEMRG4_"); break;
		case 3052:	strcpy(s_case,"GEMRG4C_"); break;
        
		case 3060:	strcpy(s_case,"GEMC45_"); break;
		case 3061:	strcpy(s_case,"GEMCG4_"); break;

		case 3070:	strcpy(s_case,"GECC45_"); break;
		case 3071:	strcpy(s_case,"GECCG4_"); break;
		case 3072:	strcpy(s_case,"GECCG4C_"); break;

		case 3080:	strcpy(s_case,"GEMP45_"); break;
		case 3081:	strcpy(s_case,"GEMPG3_"); break;

		case 3090:	strcpy(s_case,"GECM45_"); break;
		case 3091:	strcpy(s_case,"GECMG3S_"); break;
        
        case 3100:	strcpy(s_case,"GENE45_"); break;
        case 3101:	strcpy(s_case,"GENEG4C_"); break;
        
        /**/
		case 4011:	strcpy(s_case,"PLGF45_"); break;
		case 4012:	strcpy(s_case,"PLGF85_"); break;
		case 4021:	strcpy(s_case,"PLIP45_"); break;
		case 4022:	strcpy(s_case,"PLIP85_"); break;
		case 4023:	strcpy(s_case,"PLIP26_"); break;
		case 4024:	strcpy(s_case,"PLIP60_"); break;

		case 5001:	strcpy(s_case,"ISIM2H1_"); break;
		case 5002:	strcpy(s_case,"ISIM2H2_"); break;
		case 5003:	strcpy(s_case,"ISIM2H3_"); break;
		case 5004:	strcpy(s_case,"ISIM2H4_"); break;
        
        /* ISI-MIP2b: 2016/12/22 by A.Ito */
		case 5010:	strcpy(s_case,"ISM2bG-I_"); break;
		case 5011:	strcpy(s_case,"ISM2bG-II_"); break;
		case 5012:	strcpy(s_case,"ISM2bG-III_"); break;
		case 5020:	strcpy(s_case,"ISM2bI-I_"); break;
		case 5021:	strcpy(s_case,"ISM2bI-II_"); break;
		case 5022:	strcpy(s_case,"ISM2bI-III_"); break;
		case 5030:	strcpy(s_case,"ISM2bM-I_"); break;
		case 5031:	strcpy(s_case,"ISM2bM-II_"); break;
		case 5032:	strcpy(s_case,"ISM2bM-III_"); break;
		case 5040:	strcpy(s_case,"ISM2bH-I_"); break;
		case 5041:	strcpy(s_case,"ISM2bH-II_"); break;
		case 5042:	strcpy(s_case,"ISM2bH-III_"); break;
  
        /* ISI-MIP2b EX IV-VII: 2016/12/22 by A.Ito */
        case 5013:    strcpy(s_case,"ISM2bG-IV_"); break;
        case 5014:    strcpy(s_case,"ISM2bG-V_"); break;
        case 5015:    strcpy(s_case,"ISM2bG-VI_"); break;
        case 5016:    strcpy(s_case,"ISM2bG-VII_"); break;

        case 5023:    strcpy(s_case,"ISM2bI-IV_"); break;
        case 5024:    strcpy(s_case,"ISM2bI-V_"); break;
        case 5025:    strcpy(s_case,"ISM2bI-VI_"); break;
        case 5026:    strcpy(s_case,"ISM2bI-VII_"); break;

        case 5033:    strcpy(s_case,"ISM2bM-IV_"); break;
        case 5034:    strcpy(s_case,"ISM2bM-V_"); break;
        case 5035:    strcpy(s_case,"ISM2bM-VI_"); break;
        case 5036:    strcpy(s_case,"ISM2bM-VII_"); break;

        case 5043:    strcpy(s_case,"ISM2bH-IV_"); break;
        case 5044:    strcpy(s_case,"ISM2bH-V_"); break;
        case 5045:    strcpy(s_case,"ISM2bH-VI_"); break;
        case 5046:    strcpy(s_case,"ISM2bH-VII_"); break;
        
        /* ISI-MIP2b EX Ia-IIIb: 2018/03/01 by A.Ito */
        case 5060:    strcpy(s_case,"ISM2bG-Ia_"); break;
        case 5061:    strcpy(s_case,"ISM2bG-IIb_"); break;
        case 5062:    strcpy(s_case,"ISM2bG-IIIa_"); break;
        case 5063:    strcpy(s_case,"ISM2bG-IIIb_"); break;

        case 5070:    strcpy(s_case,"ISM2bI-Ia_"); break;
        case 5071:    strcpy(s_case,"ISM2bI-IIb_"); break;
        case 5072:    strcpy(s_case,"ISM2bI-IIIa_"); break;
        case 5073:    strcpy(s_case,"ISM2bI-IIIb_"); break;

        case 5080:    strcpy(s_case,"ISM2bM-Ia_"); break;
        case 5081:    strcpy(s_case,"ISM2bM-IIb_"); break;
        case 5082:    strcpy(s_case,"ISM2bM-IIIa_"); break;
        case 5083:    strcpy(s_case,"ISM2bM-IIIb_"); break;

        case 5090:    strcpy(s_case,"ISM2bH-Ia_"); break;
        case 5091:    strcpy(s_case,"ISM2bH-IIb_"); break;
        case 5092:    strcpy(s_case,"ISM2bH-IIIa_"); break;
        case 5093:    strcpy(s_case,"ISM2bH-IIIb_"); break;

        /* IMPRESSIONS *********/
		case 6001:	strcpy(s_case,"IMP1"); break;
		case 6002:	strcpy(s_case,"IMP2"); break;
        
		default:
			strcpy(s_case,"STCLIM_"); 
			printf("NULL CASE !!!");
		break;
		/* fl1 */
		/* fl2 */
		/* fl3 */
		/* nl */
		/* rcf */
		/* apf */
	}
 
    if(ISIMIP_RUN ==4 && ISIMIP2_FIXCD == 1){
        /* fixed CO2 level */
        strcat(s_case,"FXCD_");
    }
    
    /* IMPRESSIONS 1 */
    if(SCENARIO_ID == 6001){
        switch(IMP_DT){
            case 0: ; break;
            case 1: strcat(s_case,"TM03"); break;
            case 2: strcat(s_case,"TM02"); break;
            case 3: strcat(s_case,"TM01"); break;
            case 4: strcat(s_case,"T00"); break;
            case 5: strcat(s_case,"TP01"); break;
            case 6: strcat(s_case,"TP02"); break;
            case 7: strcat(s_case,"TP03"); break;
            case 8: strcat(s_case,"TP04"); break;
            case 9: strcat(s_case,"TP05"); break;
            case 10: strcat(s_case,"TP06"); break;
            case 11: strcat(s_case,"TP07"); break;
            case 12: strcat(s_case,"TP08"); break;
            case 13: strcat(s_case,"TP09"); break;
            case 14: strcat(s_case,"TP10"); break;
            case 15: strcat(s_case,"TP11"); break;
        }

        switch(IMP_DP){
            case 0: ; break;
            case 1: strcat(s_case,"PM60_"); break;
            case 2: strcat(s_case,"PM50_"); break;
            case 3: strcat(s_case,"PM40_"); break;
            case 4: strcat(s_case,"PM30_"); break;
            case 5: strcat(s_case,"PM20_"); break;
            case 6: strcat(s_case,"PM10_"); break;
            case 7: strcat(s_case,"P00_"); break;
            case 8: strcat(s_case,"PP10_"); break;
            case 9: strcat(s_case,"PP20_"); break;
            case 10: strcat(s_case,"PP30_"); break;
            case 11: strcat(s_case,"PP40_"); break;
            case 12: strcat(s_case,"PP50_"); break;
            case 13: strcat(s_case,"PP60_"); break;
        }
    }
    
    /* IMPRESSIONS 2 */
    if(SCENARIO_ID == 6002){
        switch(IMP_DT){
            case 0: ; break;
            case 1: strcat(s_case,"TM03"); break;
            case 2: strcat(s_case,"TM02"); break;
            case 3: strcat(s_case,"TM01"); break;
            case 4: strcat(s_case,"T00"); break;
            case 5: strcat(s_case,"TP005"); break;
            case 6: strcat(s_case,"TP01"); break;
            case 7: strcat(s_case,"TP02"); break;
            case 8: strcat(s_case,"TP03"); break;
            case 9: strcat(s_case,"TP04"); break;
            case 10: strcat(s_case,"TP05"); break;
            case 11: strcat(s_case,"TP06"); break;
            case 12: strcat(s_case,"TP07"); break;
            case 13: strcat(s_case,"TP08"); break;
        }

        switch(IMP_DP){
            case 0: ; break;
            case 1: strcat(s_case,"PM42"); break;
            case 2: strcat(s_case,"PM36"); break;
            case 3: strcat(s_case,"PM30"); break;
            case 4: strcat(s_case,"PM25"); break;
            case 5: strcat(s_case,"PM24"); break;
            case 6: strcat(s_case,"PM20"); break;
            case 7: strcat(s_case,"PM18"); break;
            case 8: strcat(s_case,"PM15"); break;
            case 9: strcat(s_case,"PM12"); break;
            case 10: strcat(s_case,"PM10"); break;
            case 11: strcat(s_case,"PM09"); break;
            case 12: strcat(s_case,"PM06"); break;
            case 13: strcat(s_case,"PM05"); break;
            case 14: strcat(s_case,"PM03"); break;
            case 15: strcat(s_case,"P00"); break;
            case 16: strcat(s_case,"PP03"); break;
            case 17: strcat(s_case,"PP05"); break;
            case 18: strcat(s_case,"PP06"); break;
            case 19: strcat(s_case,"PP09"); break;
            case 20: strcat(s_case,"PP10"); break;
            case 21: strcat(s_case,"PP12"); break;
            case 22: strcat(s_case,"PP15"); break;
            case 23: strcat(s_case,"PP18"); break;
        }

        switch(IMP_CD){
            case 0: ; break;
            case 1: strcat(s_case,"CD360_"); break;
            case 2: strcat(s_case,"CD423_"); break;
            case 3: strcat(s_case,"CD448_"); break;
            case 4: strcat(s_case,"CD474_"); break;
            case 5: strcat(s_case,"CD498_"); break;
            case 6: strcat(s_case,"CD517_"); break;
            case 7: strcat(s_case,"CD528_"); break;
            case 8: strcat(s_case,"CD532_"); break;
            case 9: strcat(s_case,"CD432_"); break;
            case 10: strcat(s_case,"CD469_"); break;
            case 11: strcat(s_case,"CD515_"); break;
            case 12: strcat(s_case,"CD572_"); break;
            case 13: strcat(s_case,"CD640_"); break;
            case 14: strcat(s_case,"CD718_"); break;
            case 15: strcat(s_case,"CD802_"); break;
        }
    }
}

/* REGIONS by F.Giorgi **********/
/*
 Giorgi, F., and R. Francisco (2000), Uncertainties in regional climate change prediction: 
 a regional analysis of ensemble simulations with the HADCM2 coupled AOGCM, 
 Climate Dynamics, 16, 169-182.
*/
long region_giorgi(
	double lat, 
	double lon
){
	long reg;
	
	reg = 0;
	if((lat>-48.0&&lat<-28.0) && (lon>110.0&&lon<180.0)){
		reg = 1;
	}else if((lat>-28.0&&lat<-11.0) && (lon>110.0&&lon<155.0)){
		reg = 2;
	}else if((lat>-20.0&&lat<12.0)&&(lon>-82.0&&lon<-34.0)){
		reg = 3;
	}else if((lat>-56.0&&lat<-20.0)&&(lon>-76.0&&lon<-40.0)){
		reg = 4;
	}else if((lat>10.0&&lat<30.0)&&(lon>-116.0&&lon<-83.0)){
		reg = 5;
	}else if((lat>30.0&&lat<60.0)&&(lon>-130.0&&lon<-103.0)){
		reg = 6;
	}else if((lat>30.0&&lat<50.0)&&(lon>-103.0&&lon<-85.0)){
		reg = 7;
	}else if((lat>25.0&&lat<50.0)&&(lon>-85.0&&lon<-60.0)){
		reg = 8;
	}else if((lat>60.0&&lat<72.0)&&(lon>-170.0&&lon<-103.0)){
		reg = 9;
	}else if((lat>50.0&&lat<85.0)&&(lon>-103.0&&lon<-10.0)){
		reg = 10;
	}else if((lat>30.0&&lat<48.0)&&(lon>-10.0&&lon<40.0)){
		reg = 11;
	}else if((lat>48.0&&lat<75.0)&&(lon>-10.0&&lon<40.0)){
		reg = 12;
	}else if((lat>-12.0&&lat<18.0)&&(lon>-20.0&&lon<22.0)){
		reg = 13;
	}else if((lat>-12.0&&lat<18.0)&&(lon>22.0&&lon<52.0)){
		reg = 14;
	}else if((lat>-35.0&&lat<-12.0)&&(lon>-10.0&&lon<52.0)){
		reg = 15;
	}else if((lat>18.0&&lat<30.0)&&(lon>-20.0&&lon<65.0)){
		reg = 16;
	}else if((lat>-11.0&&lat<20.0)&&(lon>95.0&&lon<155.0)){
		reg = 17;
	}else if((lat>20.0&&lat<50.0)&&(lon>100.0&&lon<145.0)){
		reg = 18;
	}else if((lat>5.0&&lat<30.0)&&(lon>65.0&&lon<100.0)){
		reg = 19;
	}else if((lat>30.0&&lat<50.0)&&(lon>40.0&&lon<75.0)){
		reg = 20;
	}else if((lat>30.0&&lat<50.0)&&(lon>75.0&&lon<100.0)){
		reg = 21;
	}else if((lat>50.0&&lat<70.0)&&(lon>40.0&&lon<180.0)){
		reg = 22;
	}
	
	return reg;
}

