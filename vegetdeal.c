/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/************ aggregate C3 and C4 community ****************/
void plant_stand(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Flux *flux
){
	long f;
	
	f = grid->m;
	
	/** mass **************************************************************************/
	mass->lai_p = (mass->c3).lai[f] + (mass->c4).lai[f];
	/* carbon */
	(mass->plant).fol = (mass->c3).fol*loct->C3ptn[f] + (mass->c4).fol*loct->C4ptn[f];
	(mass->plant).stm = (mass->c3).stm*loct->C3ptn[f] + (mass->c4).stm*loct->C4ptn[f];
	(mass->plant).rot = (mass->c3).rot*loct->C3ptn[f] + (mass->c4).rot*loct->C4ptn[f];
	(mass->plant).mfol[f] = (mass->c3).mfol[f]*loct->C3ptn[f] + (mass->c4).mfol[f]*loct->C4ptn[f];
	(mass->plant).mstm[f] = (mass->c3).mstm[f]*loct->C3ptn[f] + (mass->c4).mstm[f]*loct->C4ptn[f];
	(mass->plant).mrot[f] = (mass->c3).mrot[f]*loct->C3ptn[f] + (mass->c4).mrot[f]*loct->C4ptn[f];
	(mass->plant).lai[f] = (mass->c3).lai[f]*loct->C3ptn[f] + (mass->c4).lai[f]*loct->C4ptn[f];	
	
	/* nitrogen */
	(mass->plant).n_cnpy = (mass->c3).n_cnpy*loct->C3ptn[f] + (mass->c4).n_cnpy*loct->C4ptn[f];
	(mass->plant).n_strg = (mass->c3).n_strg*loct->C3ptn[f] + (mass->c4).n_strg*loct->C4ptn[f];
	(mass->plant).n_cnpy_m[f] = (mass->c3).n_cnpy_m[f]*loct->C3ptn[f] + (mass->c4).n_cnpy_m[f]*loct->C4ptn[f];
	(mass->plant).n_strg_m[f] = (mass->c3).n_strg_m[f]*loct->C3ptn[f] + (mass->c4).n_strg_m[f]*loct->C4ptn[f];

	/** flow **************************************************************************/
	/** production **/
	(flux->plant).gpp[f] = (flux->c3).gpp[f]*loct->C3ptn[f] + (flux->c4).gpp[f]*loct->C4ptn[f];
	(flux->plant).spp[f] = (flux->c3).spp[f]*loct->C3ptn[f] + (flux->c4).spp[f]*loct->C4ptn[f];
	(flux->plant).epp[f] = (flux->c3).epp[f]*loct->C3ptn[f] + (flux->c4).epp[f]*loct->C4ptn[f];
	(flux->plant).npp[f] = (flux->c3).npp[f]*loct->C3ptn[f] + (flux->c4).npp[f]*loct->C4ptn[f];
	
	/** maintenance respiration **/
	(flux->plant).rfm[f] = (flux->c3).rfm[f]*loct->C3ptn[f] + (flux->c4).rfm[f]*loct->C4ptn[f]; 
	(flux->plant).rcm[f] = (flux->c3).rcm[f]*loct->C3ptn[f] + (flux->c4).rcm[f]*loct->C4ptn[f]; 
	(flux->plant).rrm[f] = (flux->c3).rrm[f]*loct->C3ptn[f] + (flux->c4).rrm[f]*loct->C4ptn[f]; 
	(flux->plant).rpm[f] = (flux->plant).rfm[f] + (flux->plant).rcm[f] + (flux->plant).rrm[f];
	
	/** growth respiration **/
	(flux->plant).rfg[f] = (flux->c3).rfg[f]*loct->C3ptn[f] + (flux->c4).rfg[f]*loct->C4ptn[f]; 
	(flux->plant).rcg[f] = (flux->c3).rcg[f]*loct->C3ptn[f] + (flux->c4).rcg[f]*loct->C4ptn[f]; 
	(flux->plant).rrg[f] = (flux->c3).rrg[f]*loct->C3ptn[f] + (flux->c4).rrg[f]*loct->C4ptn[f]; 
	(flux->plant).rpg[f] = (flux->plant).rfg[f] + (flux->plant).rcg[f] + (flux->plant).rrg[f]; 
	
	/** total respiration **/
	(flux->plant).rp[f] = (flux->plant).rpm[f] + (flux->plant).rpg[f];			

	/** allocation **/
	(flux->plant).tpf[f] = (flux->c3).tpf[f]*loct->C3ptn[f] + (flux->c4).tpf[f]*loct->C4ptn[f]; 
	(flux->plant).tpc[f] = (flux->c3).tpc[f]*loct->C3ptn[f] + (flux->c4).tpc[f]*loct->C4ptn[f]; 
	(flux->plant).tpr[f] = (flux->c3).tpr[f]*loct->C3ptn[f] + (flux->c4).tpr[f]*loct->C4ptn[f]; 
	(flux->plant).tpp[f] = (flux->c3).tpp[f]*loct->C3ptn[f] + (flux->c4).tpp[f]*loct->C4ptn[f]; 
	
	/** litterfall **/
	(flux->plant).lf[f] = (flux->c3).lf[f]*loct->C3ptn[f] + (flux->c4).lf[f]*loct->C4ptn[f]; 
	(flux->plant).lc[f] = (flux->c3).lc[f]*loct->C3ptn[f] + (flux->c4).lc[f]*loct->C4ptn[f]; 
	(flux->plant).lr[f] = (flux->c3).lr[f]*loct->C3ptn[f] + (flux->c4).lr[f]*loct->C4ptn[f];
	(flux->plant).lL[f] = (flux->plant).lf[f] + (flux->plant).lc[f] + (flux->plant).lr[f];
	
	(flux->plant).lf_c[f] = (flux->c3).lf_c[f]*loct->C3ptn[f] + (flux->c4).lf_c[f]*loct->C4ptn[f]; 

	/** harvest **/
	(flux->plant).hvst[f] = (flux->c3).hvst[f]*loct->C3ptn[f] + (flux->c4).hvst[f]*loct->C4ptn[f];

	(flux->plant).emit_ch4_kirschbaum_mass[f] = (flux->c3).emit_ch4_kirschbaum_mass[f]*loct->C3ptn[f] 
		+ (flux->c4).emit_ch4_kirschbaum_mass[f]*loct->C4ptn[f];
	(flux->plant).emit_ch4_kirschbaum_photo[f] = (flux->c3).emit_ch4_kirschbaum_photo[f]*loct->C3ptn[f] 
		+ (flux->c4).emit_ch4_kirschbaum_photo[f]*loct->C4ptn[f];

	/** stable carbon isotope **/
	/* litterfall */
	(flux->c3).d13c_lL[f] = d13c_addition3(
		(flux->c3).d13c_lf[f], (flux->c3).lf[f], (flux->c3).d13c_lc[f], (flux->c3).lc[f], (flux->c3).d13c_lr[f], (flux->c3).lr[f]);
	(flux->c4).d13c_lL[f] = d13c_addition3(
		(flux->c4).d13c_lf[f], (flux->c4).lf[f], (flux->c4).d13c_lc[f], (flux->c4).lc[f], (flux->c4).d13c_lr[f], (flux->c4).lr[f]);
	(flux->plant).d13c_lL[f] = d13c_addition(
		(flux->c3).d13c_lL[f], loct->C3ptn[f]*(flux->c3).lL[f], (flux->c4).d13c_lL[f], loct->C4ptn[f]*(flux->c4).lL[f]);
	
	(flux->plant).d13c_lf[f] = d13c_addition(
		(flux->c3).d13c_lf[f], loct->C3ptn[f]*(flux->c3).lf[f], (flux->c4).d13c_lf[f], loct->C4ptn[f]*(flux->c4).lf[f]);
	(flux->plant).d13c_lc[f] = d13c_addition(
		(flux->c3).d13c_lc[f], loct->C3ptn[f]*(flux->c3).lc[f], (flux->c4).d13c_lc[f], loct->C4ptn[f]*(flux->c4).lc[f]);
	(flux->plant).d13c_lr[f] = d13c_addition(
		(flux->c3).d13c_lr[f], loct->C3ptn[f]*(flux->c3).lr[f], (flux->c4).d13c_lr[f], loct->C4ptn[f]*(flux->c4).lr[f]);
	
	/* respiration */
	(flux->plant).d13c_rfm[f] = d13c_addition(
		(flux->c3).d13c_rfm[f], loct->C3ptn[f]*(flux->c3).rfm[f], (flux->c4).d13c_rfm[f], loct->C4ptn[f]*(flux->c4).rfm[f]);
	(flux->plant).d13c_rcm[f] = d13c_addition(
		(flux->c3).d13c_rcm[f], loct->C3ptn[f]*(flux->c3).rcm[f], (flux->c4).d13c_rcm[f], loct->C4ptn[f]*(flux->c4).rcm[f]);
	(flux->plant).d13c_rrm[f] = d13c_addition(
		(flux->c3).d13c_rrm[f], loct->C3ptn[f]*(flux->c3).rrm[f], (flux->c4).d13c_rrm[f], loct->C4ptn[f]*(flux->c4).rrm[f]);
	(flux->plant).d13c_rpm[f] = d13c_addition3(
		(flux->plant).d13c_rfm[f], (flux->plant).rfm[f], (flux->plant).d13c_rcm[f], (flux->plant).rcm[f], (flux->plant).d13c_rrm[f], (flux->plant).rrm[f]);

	(flux->plant).d13c_rfg[f] = d13c_addition(
		(flux->c3).d13c_rfg[f], loct->C3ptn[f]*(flux->c3).rfg[f], (flux->c4).d13c_rfg[f], loct->C4ptn[f]*(flux->c4).rfg[f]);
	(flux->plant).d13c_rcg[f] = d13c_addition(
		(flux->c3).d13c_rcg[f], loct->C3ptn[f]*(flux->c3).rcg[f], (flux->c4).d13c_rcg[f], loct->C4ptn[f]*(flux->c4).rcg[f]);
	(flux->plant).d13c_rrg[f] = d13c_addition(
		(flux->c3).d13c_rrg[f], loct->C3ptn[f]*(flux->c3).rrg[f], (flux->c4).d13c_rrg[f], loct->C4ptn[f]*(flux->c4).rrg[f]);
	(flux->plant).d13c_rpg[f] = d13c_addition3(
		(flux->plant).d13c_rfg[f], (flux->plant).rfg[f], (flux->plant).d13c_rcg[f], (flux->plant).rcg[f], (flux->plant).d13c_rrg[f], (flux->plant).rrg[f]);

	(flux->c3).d13c_rpm[f] = d13c_addition3(
		(flux->c3).d13c_rfm[f], (flux->c3).rfm[f], (flux->c3).d13c_rcm[f], (flux->c3).rcm[f], (flux->c3).d13c_rrm[f], (flux->c3).rrm[f]);
	(flux->c4).d13c_rpm[f] = d13c_addition3(
		(flux->c4).d13c_rfm[f], (flux->c4).rfm[f], (flux->c4).d13c_rcm[f], (flux->c4).rcm[f], (flux->c4).d13c_rrm[f], (flux->c4).rrm[f]);
	(flux->c3).d13c_rpg[f] = d13c_addition3(
		(flux->c3).d13c_rfg[f], (flux->c3).rfg[f], (flux->c3).d13c_rcg[f], (flux->c3).rcg[f], (flux->c3).d13c_rrg[f], (flux->c3).rrg[f]);
	(flux->c4).d13c_rpg[f] = d13c_addition3(
		(flux->c4).d13c_rfg[f], (flux->c4).rfg[f], (flux->c4).d13c_rcg[f], (flux->c4).rcg[f], (flux->c4).d13c_rrg[f], (flux->c4).rrg[f]);

	(flux->c3).d13c_rp[f] = d13c_addition(
		(flux->c3).d13c_rpm[f], (flux->c3).rpm[f], (flux->c3).d13c_rpg[f], (flux->c3).rpg[f]);
	(flux->c4).d13c_rp[f] = d13c_addition(
		(flux->c4).d13c_rpm[f], (flux->c4).rpm[f], (flux->c4).d13c_rpg[f], (flux->c4).rpg[f]);
	(flux->plant).d13c_rp[f] = d13c_addition(
		(flux->plant).d13c_rpm[f],(flux->plant).rpm[f], (flux->plant).d13c_rpg[f], (flux->plant).rpg[f]);
	
	/* nitrogen */
	(flux->plant).n_abdn_cnpy[f] = (flux->c3).n_abdn_cnpy[f]*loct->C3ptn[f] + (flux->c4).n_abdn_cnpy[f]*loct->C4ptn[f];
	(flux->plant).n_abdn_strg[f] = (flux->c3).n_abdn_strg[f]*loct->C3ptn[f] + (flux->c4).n_abdn_strg[f]*loct->C4ptn[f];
	(flux->plant).uptake_no3[f] = (flux->c3).uptake_no3[f]*loct->C3ptn[f] + (flux->c4).uptake_no3[f]*loct->C4ptn[f];
	(flux->plant).uptake_nh4[f] = (flux->c3).uptake_nh4[f]*loct->C3ptn[f] + (flux->c4).uptake_nh4[f]*loct->C4ptn[f];
	(flux->plant).n_biofix[f] = (flux->c3).n_biofix[f]*loct->C3ptn[f] + (flux->c4).n_biofix[f]*loct->C4ptn[f];

	(flux->plant).n_alloc_cnpy[f] = (flux->c3).n_alloc_cnpy[f]*loct->C3ptn[f] + (flux->c4).n_alloc_cnpy[f]*loct->C4ptn[f];
	(flux->plant).n_alloc_strg[f] = (flux->c3).n_alloc_strg[f]*loct->C3ptn[f] + (flux->c4).n_alloc_strg[f]*loct->C4ptn[f];
	(flux->plant).n_realloc[f] = (flux->c3).n_realloc[f]*loct->C3ptn[f] + (flux->c4).n_realloc[f]*loct->C4ptn[f];
	(flux->plant).n_salvage[f] = (flux->c3).n_salvage[f]*loct->C3ptn[f] + (flux->c4).n_salvage[f]*loct->C4ptn[f];
}

/******** dealings before calculating carbon budget **********/
void beforedeal(
	struct Grid *grid, 
	struct Pflx *flux
){
	plant_flux_zero(grid->m, flux);
}

/******** dealings after calculating carbon budget **********/
void afterdeal(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){
	mass->lai[grid->m] = lai_mass(mass, pchar);

	flux->rpg[grid->m] = flux->rfg[grid->m] + flux->rcg[grid->m] + flux->rrg[grid->m];
	flux->rp[grid->m] = flux->rpm[grid->m] + flux->rpg[grid->m];
	flux->spp[grid->m] = flux->gpp[grid->m] - flux->rfm[grid->m]-flux->rfg[grid->m];
	flux->npp[grid->m] = flux->gpp[grid->m] - flux->rpm[grid->m]-flux->rpg[grid->m];

	flux->lL[grid->m] = flux->lf[grid->m] + flux->lc[grid->m] + flux->lr[grid->m] + flux->lf_c[grid->m];
	
	/* monthly mass */	
	mass->mfol[grid->m] = mass->fol;
	mass->mstm[grid->m] = mass->stm;
	mass->mrot[grid->m] = mass->rot;
	mass->plant[grid->m] = mass->fol + mass->stm + mass->rot;

	/* monthly d13c of mass */
	mass->d13c_mfol[grid->m] = mass->d13c_fol;
	mass->d13c_mstm[grid->m] = mass->d13c_stm;
	mass->d13c_mrot[grid->m] = mass->d13c_rot;
}

/**********************************************************/
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
	
	if(lat1>90.0) lat1=90.0;
	if(lat2<-90.0) lat2=-90.0;
	
	m_lat = (lat1+lat2)/2.0 * PI/180.0;
	e_exc = sqrt(2.0/e_flat - 1.0/(e_flat*e_flat));
	
	aa = 1.0 - e_exc*e_exc*sin(m_lat)*sin(m_lat);
	l_lat = PI/180.0*e_rad*(1.0 - e_exc*e_exc)/pow(aa, 1.5) * fabs(lat1 - lat2);
	
	aa = 1.0 - e_exc*e_exc*sin(lat1)*sin(lat1);
	l_lon1 = PI/180.0*e_rad*cos(lat1*PI/180.0)/sqrt(aa) * fabs(lon1-lon2);
	aa = 1.0 - e_exc*e_exc*sin(lat2)*sin(lat2);
	l_lon2 = PI/180.0*e_rad*cos(lat2*PI/180.0)/sqrt(aa) * fabs(lon1-lon2);
	
	area = (l_lon1+l_lon2)*l_lat/2.0 / 10000.0;
	
	return (area); /* ha */
}

/******************************************************/
long basin_id_trip(
	long original
){
	long basin;

	switch(original){
		case 1: basin=1; break;
		case 4: basin=2; break;
		case 2: basin=3; break;
		case 3: basin=4; break;
		case 10: basin=5; break;
		case 5: basin=6; break;
		case 7: basin=7; break;
		case 6: basin=8; break;
		case 8: basin=9; break;
		case 9: basin=10; break;
		case 14: basin=11; break;
		case 11: basin=12; break;
		case 12: basin=13; break;
		case 21: basin=14; break;
		case 34: basin=14; break;
		case 13: basin=15; break;
		case 16: basin=16; break;
		case 20: basin=17; break;
		case 18: basin=18; break;
		case 68: basin=18; break;
		case 23: basin=19; break;
		case 17: basin=20; break;
		case 31: basin=21; break;
		case 19: basin=22; break;
		case 33: basin=23; break;
		case 26: basin=24; break;
		case 45: basin=25; break;
		case 32: basin=26; break;
		case 38: basin=27; break;
		case 41: basin=28; break;
		case 28: basin=29; break;
		case 29: basin=30; break;
		case 24: basin=31; break;
		case 36: basin=32; break;
		case 39: basin=33; break;
		case 40: basin=34; break;
		case 42: basin=35; break;
		case 51: basin=36; break;
		case 296: basin=37; break;
		case 47: basin=38; break;
		case 48: basin=39; break;
		case 54: basin=40; break;
		case 56: basin=41; break;
		case 53: basin=42; break;
		case 57: basin=43; break;
		case 110: basin=44; break;
		case 64: basin=45; break;
		case 78: basin=46; break;
		case 62: basin=47; break;
		case 74: basin=48; break;
		case 81: basin=49; break;
		case 79: basin=50; break;
		case 83: basin=51; break;
		case 102: basin=52; break;
		case 378: basin=53; break;
		
		default:  basin=0;
	}
	
	return basin;
}
