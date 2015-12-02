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
	(flux->plant).hvst[f] = (flux->c3).hvst[f]*loct->c3ptn[f] + (flux->c4).hvst[f]*loct->c4ptn[f];

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
	if(GCM_ID == 0){
		GCM_R = 1;
		GCM_C = 1;
	}else if(GCM_ID==1 || GCM_ID==2 || GCM_ID==3 || GCM_ID==4 || GCM_ID==5 || GCM_ID==6){
		GCM_R = 32;
		GCM_C = 64;
	}else if(GCM_ID==11 || GCM_ID==12 || GCM_ID==13 || GCM_ID==14 || GCM_ID==15 || GCM_ID==16 || GCM_ID==17 || GCM_ID==18){
		GCM_R = 48;
		GCM_C = 96;
	}else if(GCM_ID==21 || GCM_ID==22 || GCM_ID==23 || GCM_ID==24){
		GCM_R = 73;
		GCM_C = 96;
	}else if(GCM_ID==31 || GCM_ID==32){
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM_ID==41 || GCM_ID==42 || GCM_ID==43 || GCM_ID==44){
		GCM_R = 56;
		GCM_C = 64;
	}else if(GCM_ID==51 || GCM_ID==52){
		GCM_R = 80;
		GCM_C = 96;
	}else if(GCM_ID==61 || GCM_ID==62 || GCM_ID==63){
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM_ID==71){
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM_ID==1000 || GCM_ID==1001){ /* MIROC-HIGH */
		GCM_R = 160;
		GCM_C = 320;
	}else if(GCM_ID==1010 || GCM_ID==1011 || GCM_ID==1012 || GCM_ID==1013 || GCM_ID==1014 || 
			GCM_ID==1015 || GCM_ID==1016 || GCM_ID==1017 || GCM_ID==1018){ /* MIROC-MED */
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM_ID==1050 || GCM_ID==1051 ||GCM_ID==1052){ /* BCCR */
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM_ID==1060 || GCM_ID==1061 ||GCM_ID==1062){ /* INM-CM3 : error-fixed 2009/04/26 A.Ito */
		GCM_R = 45;
		GCM_C = 72;
	}else if(GCM_ID==1070 || GCM_ID==1071 || GCM_ID==1072){ /* IPSL */
		GCM_R = 72;
		GCM_C = 96;
	}else if(GCM_ID==1080 || GCM_ID==1081 || GCM_ID==1082){ /* GFDL CM2-1 */
		GCM_R = 90;
		GCM_C = 144;
	}else if(GCM_ID==1090 || GCM_ID==1091 || GCM_ID==1092){ /* HadCM3 */
		GCM_R = 73;
		GCM_C = 96;
	}else if(GCM_ID==1100 || GCM_ID==1101 || GCM_ID==1102 || GCM_ID==1103){ /* GISS AOM */
		GCM_R = 60;
		GCM_C = 90;
	}else if(GCM_ID==1110 || GCM_ID==1111 || GCM_ID==1112){ /* CCCma-T63 */
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM_ID==1120 || GCM_ID==1121 || GCM_ID==1122 || GCM_ID==1270 || GCM_ID==1271 || GCM_ID==1272){ /* CSIRO */
		GCM_R = 96;
		GCM_C = 192;
	}else if(GCM_ID==1130 || GCM_ID==1131 || GCM_ID==1132 || GCM_ID==1133 || GCM_ID==1134
			 || GCM_ID==1135 || GCM_ID==1136 || GCM_ID==1137 || GCM_ID==1138 || GCM_ID==1139
			  || GCM_ID==1140 || GCM_ID==1141 || GCM_ID==1142 || GCM_ID==1143 || GCM_ID==1144){ /* MRI */
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM_ID==1150 ||GCM_ID==1151 ||GCM_ID==1152 ||GCM_ID==1153 ||GCM_ID==1154 ||
			GCM_ID==1155 ||GCM_ID==1156 ||GCM_ID==1157 ||GCM_ID==1158 ||GCM_ID==1159){ /* MPI ECHAM */
		GCM_R = 96;
		GCM_C = 192;
	}else if(GCM_ID==1160 || GCM_ID==1161 || GCM_ID==1162 || GCM_ID==1163 || GCM_ID==1164 || GCM_ID==1165){ /* IAP */
		GCM_R = 60;
		GCM_C = 128;
	}else if(GCM_ID==1170 || GCM_ID==1171){ /* HadGEM */
		GCM_R = 145;
		GCM_C = 192;
	}else if(GCM_ID==1180 || GCM_ID==1181 || GCM_ID==1182){ /* GFDL CM 2.0 */
		GCM_R = 90;
		GCM_C = 144;
	}else if(GCM_ID==1190 || GCM_ID==1191 || GCM_ID==1192){ /* CNRM */
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM_ID==1200 || GCM_ID==1201 || GCM_ID==1202 || GCM_ID==1203 || GCM_ID==1204
			 || GCM_ID==1205 || GCM_ID==1206 || GCM_ID==1207 || GCM_ID==1208 || GCM_ID==1209
			 || GCM_ID==1210 || GCM_ID==1211 || GCM_ID==1212 || GCM_ID==1213 || GCM_ID==1214){ /* CCC47 */
		GCM_R = 48;
		GCM_C = 96;
	}else if(GCM_ID==1220 || GCM_ID==1221 || GCM_ID==1222 || GCM_ID==1223 || GCM_ID==1224
			 || GCM_ID==1225 || GCM_ID==1226 || GCM_ID==1227 || GCM_ID==1228 || GCM_ID==1229
			  || GCM_ID==1230 || GCM_ID==1231 || GCM_ID==1232 || GCM_ID==1233 || GCM_ID==1234
			   || GCM_ID==1235 || GCM_ID==1236 || GCM_ID==1237 || GCM_ID==1238 || GCM_ID==1239){ /* NCAR CCSM */
		GCM_R = 128;
		GCM_C = 256;
	}else if(GCM_ID==1240 || GCM_ID==1241 || GCM_ID==1242){ /* GISS E */
		GCM_R = 46;
		GCM_C = 72;
	}else if(GCM_ID==1250 || GCM_ID==1251 || GCM_ID==1252 || GCM_ID==1253){ /* GISS R */
		GCM_R = 46;
		GCM_C = 72;
	}else if(GCM_ID==1260 || GCM_ID==1261 || GCM_ID==1262 || GCM_ID==1263 || GCM_ID==1264
			 || GCM_ID==1265 || GCM_ID==1266 || GCM_ID==1267){ /* NCAR PCM */
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM_ID==3000 || GCM_ID==3001 || GCM_ID==3002){
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM_ID==3010 || GCM_ID==3011 || GCM_ID==3012){
		GCM_R = 56;
		GCM_C = 64;
	}else if(GCM_ID==3020 || GCM_ID==3021 || GCM_ID==3022){
		GCM_R = 90;
		GCM_C = 144;
	}else if(GCM_ID==3030 || GCM_ID==3031 || GCM_ID==3032 || GCM_ID==3033 || GCM_ID==3034 || GCM_ID==3035){
		GCM_R = 145;
		GCM_C = 192;
	}else if(GCM_ID==3040 || GCM_ID==3041 || GCM_ID==3042){
		GCM_R = 96;
		GCM_C = 96;
	}else if(GCM_ID==3050 || GCM_ID==3051 || GCM_ID==3052){
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM_ID==3060 || GCM_ID==3061){
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM_ID==3070 || GCM_ID==3071 || GCM_ID==3072){
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM_ID==3080 || GCM_ID==3081){
		GCM_R = 96;
		GCM_C = 192;
	}else if(GCM_ID==3090 || GCM_ID==3091){
		GCM_R = 192;
		GCM_C = 288;
	}else if(GCM_ID==3100 || GCM_ID==3101){
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
	switch(GCM_ID){
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
        
        /* IMPRESSIONS *********/
		case 6001:	strcpy(s_case,"IMP_M03M60_"); break;
		case 6002:	strcpy(s_case,"IMP_M02M60_"); break;
		case 6003:	strcpy(s_case,"IMP_M01M60_"); break;
		case 6004:	strcpy(s_case,"IMP_000M60_"); break;
		case 6005:	strcpy(s_case,"IMP_P01M60_"); break;
		case 6006:	strcpy(s_case,"IMP_P02M60_"); break;
		case 6007:	strcpy(s_case,"IMP_P03M60_"); break;
		case 6008:	strcpy(s_case,"IMP_P04M60_"); break;
		case 6009:	strcpy(s_case,"IMP_P05M60_"); break;
		case 6010:	strcpy(s_case,"IMP_P06M60_"); break;
		case 6011:	strcpy(s_case,"IMP_P07M60_"); break;
		case 6012:	strcpy(s_case,"IMP_P08M60_"); break;
		case 6013:	strcpy(s_case,"IMP_P09M60_"); break;
		case 6014:	strcpy(s_case,"IMP_P10M60_"); break;
		case 6015:	strcpy(s_case,"IMP_P11M60_"); break;

		case 6051:	strcpy(s_case,"IMP_M03M50_"); break;
		case 6052:	strcpy(s_case,"IMP_M02M50_"); break;
		case 6053:	strcpy(s_case,"IMP_M01M50_"); break;
		case 6054:	strcpy(s_case,"IMP_000M50_"); break;
		case 6055:	strcpy(s_case,"IMP_P01M50_"); break;
		case 6056:	strcpy(s_case,"IMP_P02M50_"); break;
		case 6057:	strcpy(s_case,"IMP_P03M50_"); break;
		case 6058:	strcpy(s_case,"IMP_P04M50_"); break;
		case 6059:	strcpy(s_case,"IMP_P05M50_"); break;
		case 6060:	strcpy(s_case,"IMP_P06M50_"); break;
		case 6061:	strcpy(s_case,"IMP_P07M50_"); break;
		case 6062:	strcpy(s_case,"IMP_P08M50_"); break;
		case 6063:	strcpy(s_case,"IMP_P09M50_"); break;
		case 6064:	strcpy(s_case,"IMP_P10M50_"); break;
		case 6065:	strcpy(s_case,"IMP_P11M50_"); break;

		case 6101:	strcpy(s_case,"IMP_M03M40_"); break;
		case 6102:	strcpy(s_case,"IMP_M02M40_"); break;
		case 6103:	strcpy(s_case,"IMP_M01M40_"); break;
		case 6104:	strcpy(s_case,"IMP_000M40_"); break;
		case 6105:	strcpy(s_case,"IMP_P01M40_"); break;
		case 6106:	strcpy(s_case,"IMP_P02M40_"); break;
		case 6107:	strcpy(s_case,"IMP_P03M40_"); break;
		case 6108:	strcpy(s_case,"IMP_P04M40_"); break;
		case 6109:	strcpy(s_case,"IMP_P05M40_"); break;
		case 6110:	strcpy(s_case,"IMP_P06M40_"); break;
		case 6111:	strcpy(s_case,"IMP_P07M40_"); break;
		case 6112:	strcpy(s_case,"IMP_P08M40_"); break;
		case 6113:	strcpy(s_case,"IMP_P09M40_"); break;
		case 6114:	strcpy(s_case,"IMP_P10M40_"); break;
		case 6115:	strcpy(s_case,"IMP_P11M40_"); break;

		case 6151:	strcpy(s_case,"IMP_M03M30_"); break;
		case 6152:	strcpy(s_case,"IMP_M02M30_"); break;
		case 6153:	strcpy(s_case,"IMP_M01M30_"); break;
		case 6154:	strcpy(s_case,"IMP_000M30_"); break;
		case 6155:	strcpy(s_case,"IMP_P01M30_"); break;
		case 6156:	strcpy(s_case,"IMP_P02M30_"); break;
		case 6157:	strcpy(s_case,"IMP_P03M30_"); break;
		case 6158:	strcpy(s_case,"IMP_P04M30_"); break;
		case 6159:	strcpy(s_case,"IMP_P05M30_"); break;
		case 6160:	strcpy(s_case,"IMP_P06M30_"); break;
		case 6161:	strcpy(s_case,"IMP_P07M30_"); break;
		case 6162:	strcpy(s_case,"IMP_P08M30_"); break;
		case 6163:	strcpy(s_case,"IMP_P09M30_"); break;
		case 6164:	strcpy(s_case,"IMP_P10M30_"); break;
		case 6165:	strcpy(s_case,"IMP_P11M30_"); break;

		case 6201:	strcpy(s_case,"IMP_M03M20_"); break;
		case 6202:	strcpy(s_case,"IMP_M02M20_"); break;
		case 6203:	strcpy(s_case,"IMP_M01M20_"); break;
		case 6204:	strcpy(s_case,"IMP_000M20_"); break;
		case 6205:	strcpy(s_case,"IMP_P01M20_"); break;
		case 6206:	strcpy(s_case,"IMP_P02M20_"); break;
		case 6207:	strcpy(s_case,"IMP_P03M20_"); break;
		case 6208:	strcpy(s_case,"IMP_P04M20_"); break;
		case 6209:	strcpy(s_case,"IMP_P05M20_"); break;
		case 6210:	strcpy(s_case,"IMP_P06M20_"); break;
		case 6211:	strcpy(s_case,"IMP_P07M20_"); break;
		case 6212:	strcpy(s_case,"IMP_P08M20_"); break;
		case 6213:	strcpy(s_case,"IMP_P09M20_"); break;
		case 6214:	strcpy(s_case,"IMP_P10M20_"); break;
		case 6215:	strcpy(s_case,"IMP_P11M20_"); break;

		case 6251:	strcpy(s_case,"IMP_M03M10_"); break;
		case 6252:	strcpy(s_case,"IMP_M02M10_"); break;
		case 6253:	strcpy(s_case,"IMP_M01M10_"); break;
		case 6254:	strcpy(s_case,"IMP_000M10_"); break;
		case 6255:	strcpy(s_case,"IMP_P01M10_"); break;
		case 6256:	strcpy(s_case,"IMP_P02M10_"); break;
		case 6257:	strcpy(s_case,"IMP_P03M10_"); break;
		case 6258:	strcpy(s_case,"IMP_P04M10_"); break;
		case 6259:	strcpy(s_case,"IMP_P05M10_"); break;
		case 6260:	strcpy(s_case,"IMP_P06M10_"); break;
		case 6261:	strcpy(s_case,"IMP_P07M10_"); break;
		case 6262:	strcpy(s_case,"IMP_P08M10_"); break;
		case 6263:	strcpy(s_case,"IMP_P09M10_"); break;
		case 6264:	strcpy(s_case,"IMP_P10M10_"); break;
		case 6265:	strcpy(s_case,"IMP_P11M10_"); break;

		case 6301:	strcpy(s_case,"IMP_M03000_"); break;
		case 6302:	strcpy(s_case,"IMP_M02000_"); break;
		case 6303:	strcpy(s_case,"IMP_M01000_"); break;
		case 6304:	strcpy(s_case,"IMP_000000_"); break;
		case 6305:	strcpy(s_case,"IMP_P01000_"); break;
		case 6306:	strcpy(s_case,"IMP_P02000_"); break;
		case 6307:	strcpy(s_case,"IMP_P03000_"); break;
		case 6308:	strcpy(s_case,"IMP_P04000_"); break;
		case 6309:	strcpy(s_case,"IMP_P05000_"); break;
		case 6310:	strcpy(s_case,"IMP_P06000_"); break;
		case 6311:	strcpy(s_case,"IMP_P07000_"); break;
		case 6312:	strcpy(s_case,"IMP_P08000_"); break;
		case 6313:	strcpy(s_case,"IMP_P09000_"); break;
		case 6314:	strcpy(s_case,"IMP_P10000_"); break;
		case 6315:	strcpy(s_case,"IMP_P11000_"); break;

		case 6351:	strcpy(s_case,"IMP_M03P10_"); break;
		case 6352:	strcpy(s_case,"IMP_M02P10_"); break;
		case 6353:	strcpy(s_case,"IMP_M01P10_"); break;
		case 6354:	strcpy(s_case,"IMP_000P10_"); break;
		case 6355:	strcpy(s_case,"IMP_P01P10_"); break;
		case 6356:	strcpy(s_case,"IMP_P02P10_"); break;
		case 6357:	strcpy(s_case,"IMP_P03P10_"); break;
		case 6358:	strcpy(s_case,"IMP_P04P10_"); break;
		case 6359:	strcpy(s_case,"IMP_P05P10_"); break;
		case 6360:	strcpy(s_case,"IMP_P06P10_"); break;
		case 6361:	strcpy(s_case,"IMP_P07P10_"); break;
		case 6362:	strcpy(s_case,"IMP_P08P10_"); break;
		case 6363:	strcpy(s_case,"IMP_P09P10_"); break;
		case 6364:	strcpy(s_case,"IMP_P10P10_"); break;
		case 6365:	strcpy(s_case,"IMP_P11P10_"); break;

		case 6401:	strcpy(s_case,"IMP_M03P20_"); break;
		case 6402:	strcpy(s_case,"IMP_M02P20_"); break;
		case 6403:	strcpy(s_case,"IMP_M01P20_"); break;
		case 6404:	strcpy(s_case,"IMP_000P20_"); break;
		case 6405:	strcpy(s_case,"IMP_P01P20_"); break;
		case 6406:	strcpy(s_case,"IMP_P02P20_"); break;
		case 6407:	strcpy(s_case,"IMP_P03P20_"); break;
		case 6408:	strcpy(s_case,"IMP_P04P20_"); break;
		case 6409:	strcpy(s_case,"IMP_P05P20_"); break;
		case 6410:	strcpy(s_case,"IMP_P06P20_"); break;
		case 6411:	strcpy(s_case,"IMP_P07P20_"); break;
		case 6412:	strcpy(s_case,"IMP_P08P20_"); break;
		case 6413:	strcpy(s_case,"IMP_P09P20_"); break;
		case 6414:	strcpy(s_case,"IMP_P10P20_"); break;
		case 6415:	strcpy(s_case,"IMP_P11P20_"); break;

		case 6451:	strcpy(s_case,"IMP_M03P30_"); break;
		case 6452:	strcpy(s_case,"IMP_M02P30_"); break;
		case 6453:	strcpy(s_case,"IMP_M01P30_"); break;
		case 6454:	strcpy(s_case,"IMP_000P30_"); break;
		case 6455:	strcpy(s_case,"IMP_P01P30_"); break;
		case 6456:	strcpy(s_case,"IMP_P02P30_"); break;
		case 6457:	strcpy(s_case,"IMP_P03P30_"); break;
		case 6458:	strcpy(s_case,"IMP_P04P30_"); break;
		case 6459:	strcpy(s_case,"IMP_P05P30_"); break;
		case 6460:	strcpy(s_case,"IMP_P06P30_"); break;
		case 6461:	strcpy(s_case,"IMP_P07P30_"); break;
		case 6462:	strcpy(s_case,"IMP_P08P30_"); break;
		case 6463:	strcpy(s_case,"IMP_P09P30_"); break;
		case 6464:	strcpy(s_case,"IMP_P10P30_"); break;
		case 6465:	strcpy(s_case,"IMP_P11P30_"); break;

		case 6501:	strcpy(s_case,"IMP_M03P40_"); break;
		case 6502:	strcpy(s_case,"IMP_M02P40_"); break;
		case 6503:	strcpy(s_case,"IMP_M01P40_"); break;
		case 6504:	strcpy(s_case,"IMP_000P40_"); break;
		case 6505:	strcpy(s_case,"IMP_P01P40_"); break;
		case 6506:	strcpy(s_case,"IMP_P02P40_"); break;
		case 6507:	strcpy(s_case,"IMP_P03P40_"); break;
		case 6508:	strcpy(s_case,"IMP_P04P40_"); break;
		case 6509:	strcpy(s_case,"IMP_P05P40_"); break;
		case 6510:	strcpy(s_case,"IMP_P06P40_"); break;
		case 6511:	strcpy(s_case,"IMP_P07P40_"); break;
		case 6512:	strcpy(s_case,"IMP_P08P40_"); break;
		case 6513:	strcpy(s_case,"IMP_P09P40_"); break;
		case 6514:	strcpy(s_case,"IMP_P10P40_"); break;
		case 6515:	strcpy(s_case,"IMP_P11P40_"); break;

		case 6551:	strcpy(s_case,"IMP_M03P50_"); break;
		case 6552:	strcpy(s_case,"IMP_M02P50_"); break;
		case 6553:	strcpy(s_case,"IMP_M01P50_"); break;
		case 6554:	strcpy(s_case,"IMP_000P50_"); break;
		case 6555:	strcpy(s_case,"IMP_P01P50_"); break;
		case 6556:	strcpy(s_case,"IMP_P02P50_"); break;
		case 6557:	strcpy(s_case,"IMP_P03P50_"); break;
		case 6558:	strcpy(s_case,"IMP_P04P50_"); break;
		case 6559:	strcpy(s_case,"IMP_P05P50_"); break;
		case 6560:	strcpy(s_case,"IMP_P06P50_"); break;
		case 6561:	strcpy(s_case,"IMP_P07P50_"); break;
		case 6562:	strcpy(s_case,"IMP_P08P50_"); break;
		case 6563:	strcpy(s_case,"IMP_P09P50_"); break;
		case 6564:	strcpy(s_case,"IMP_P10P50_"); break;
		case 6565:	strcpy(s_case,"IMP_P11P50_"); break;

		case 6601:	strcpy(s_case,"IMP_M03P60_"); break;
		case 6602:	strcpy(s_case,"IMP_M02P60_"); break;
		case 6603:	strcpy(s_case,"IMP_M01P60_"); break;
		case 6604:	strcpy(s_case,"IMP_000P60_"); break;
		case 6605:	strcpy(s_case,"IMP_P01P60_"); break;
		case 6606:	strcpy(s_case,"IMP_P02P60_"); break;
		case 6607:	strcpy(s_case,"IMP_P03P60_"); break;
		case 6608:	strcpy(s_case,"IMP_P04P60_"); break;
		case 6609:	strcpy(s_case,"IMP_P05P60_"); break;
		case 6610:	strcpy(s_case,"IMP_P06P60_"); break;
		case 6611:	strcpy(s_case,"IMP_P07P60_"); break;
		case 6612:	strcpy(s_case,"IMP_P08P60_"); break;
		case 6613:	strcpy(s_case,"IMP_P09P60_"); break;
		case 6614:	strcpy(s_case,"IMP_P10P60_"); break;
		case 6615:	strcpy(s_case,"IMP_P11P60_"); break;

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
	if((lat>-48.0&&lat<-28.0)&&(lon>110.0&&lon<180.0)){
		reg = 1;
	}else if((lat>-28.0&&lat<-11.0)&&(lon>110.0&&lon<155.0)){
		reg = 2;
	}else if((lat>-20.0&&lat<12.0)&&(lon>-82.0&&lon<-34.0)){
		reg=3;
	}else if((lat>-56.0&&lat<-20.0)&&(lon>-76.0&&lon<-40.0)){
		reg=4;
	}else if((lat>10.0&&lat<30.0)&&(lon>-116.0&&lon<-83.0)){
		reg=5;
	}else if((lat>30.0&&lat<60.0)&&(lon>-130.0&&lon<-103.0)){
		reg=6;
	}else if((lat>30.0&&lat<50.0)&&(lon>-103.0&&lon<-85.0)){
		reg=7;
	}else if((lat>25.0&&lat<50.0)&&(lon>-85.0&&lon<-60.0)){
		reg=8;
	}else if((lat>60.0&&lat<72.0)&&(lon>-170.0&&lon<-103.0)){
		reg=9;
	}else if((lat>50.0&&lat<85.0)&&(lon>-103.0&&lon<-10.0)){
		reg=10;
	}else if((lat>30.0&&lat<48.0)&&(lon>-10.0&&lon<40.0)){
		reg=11;
	}else if((lat>48.0&&lat<75.0)&&(lon>-10.0&&lon<40.0)){
		reg=12;
	}else if((lat>-12.0&&lat<18.0)&&(lon>-20.0&&lon<22.0)){
		reg=13;
	}else if((lat>-12.0&&lat<18.0)&&(lon>22.0&&lon<52.0)){
		reg=14;
	}else if((lat>-35.0&&lat<-12.0)&&(lon>-10.0&&lon<52.0)){
		reg=15;
	}else if((lat>18.0&&lat<30.0)&&(lon>-20.0&&lon<65.0)){
		reg=16;
	}else if((lat>-11.0&&lat<20.0)&&(lon>95.0&&lon<155.0)){
		reg=17;
	}else if((lat>20.0&&lat<50.0)&&(lon>100.0&&lon<145.0)){
		reg=18;
	}else if((lat>5.0&&lat<30.0)&&(lon>65.0&&lon<100.0)){
		reg=19;
	}else if((lat>30.0&&lat<50.0)&&(lon>40.0&&lon<75.0)){
		reg=20;
	}else if((lat>30.0&&lat<50.0)&&(lon>75.0&&lon<100.0)){
		reg=21;
	}else if((lat>50.0&&lat<70.0)&&(lon>40.0&&lon<180.0)){
		reg=22;
	}
	
	return reg;
}

