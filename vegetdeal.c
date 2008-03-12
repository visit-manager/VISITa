/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
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

extern long GCM, GCM_R, GCM_C;

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
	/* LAI update */
	mass->lai[grid->m] = lai_mass(mass, pchar);
	
	/* total respiration */
	flux->rpg[grid->m] = flux->rfg[grid->m] + flux->rcg[grid->m] + flux->rrg[grid->m];
	flux->rp[grid->m] = flux->rpm[grid->m] + flux->rpg[grid->m];
	flux->spp[grid->m] = flux->gpp[grid->m] - flux->rfm[grid->m]-flux->rfg[grid->m];
	flux->npp[grid->m] = flux->gpp[grid->m] - flux->rpm[grid->m]-flux->rpg[grid->m];

	/* total litterfall */
	flux->lL[grid->m] = flux->lf[grid->m] + flux->lc[grid->m] + flux->lr[grid->m] + flux->lf_c[grid->m];
		
	/* monthly mass **********************************/	
	mass->mfol[grid->m] = mass->fol;
	mass->mstm[grid->m] = mass->stm;
	mass->mrot[grid->m] = mass->rot;
	mass->plant[grid->m] = mass->fol + mass->stm + mass->rot;
	
	mass->fol_p = mass->fol;

	/* monthly d13c of mass */
	mass->d13c_mfol[grid->m] = mass->d13c_fol;
	mass->d13c_mstm[grid->m] = mass->d13c_stm;
	mass->d13c_mrot[grid->m] = mass->d13c_rot;
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
	aa = 1.0 - e_exc*e_exc*sin(lat1*PI/180.0)*sin(lat1*PI/180.0);
	l_lon1 = PI/180.0*e_rad*cos(lat1*PI/180.0)/sqrt(aa) * fabs(lon1-lon2);
	aa = 1.0 - e_exc*e_exc*sin(lat2*PI/180.0)*sin(lat2*PI/180.0);
	l_lon2 = PI/180.0*e_rad*cos(lat2*PI/180.0)/sqrt(aa) * fabs(lon1-lon2);
	
	area = (l_lon1 + l_lon2)*l_lat/2.0 / 10000.0;
	
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

/******************************************************************************/
void set_rowcol_gcm(
	void
){
	if(GCM==0){
		GCM_R = 1;
		GCM_C = 1;
	}else if(GCM==1 || GCM==2 || GCM==3 || GCM==4 || GCM==5 || GCM==6){
		GCM_R = 32;
		GCM_C = 64;
	}else if(GCM==11 || GCM==12 || GCM==13 || GCM==14 || GCM==15 || GCM==16 || GCM==17 || GCM==18){
		GCM_R = 48;
		GCM_C = 96;
	}else if(GCM==21 || GCM==22 || GCM==23 || GCM==24){
		GCM_R = 73;
		GCM_C = 96;
	}else if(GCM==31 || GCM==32){
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM==41 || GCM==42 || GCM==43 || GCM==44){
		GCM_R = 56;
		GCM_C = 64;
	}else if(GCM==51 || GCM==52){
		GCM_R = 80;
		GCM_C = 96;
	}else if(GCM==61 || GCM==62 || GCM==63){
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM==71){
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM==1000 || GCM==1001){ /* MIROC-HIGH */
		GCM_R = 160;
		GCM_C = 320;
	}else if(GCM==1010 || GCM==1011 || GCM==1012 || GCM==1013 || GCM==1014 || 
			GCM==1015 || GCM==1016 || GCM==1017 || GCM==1018){ /* MIROC-MED */
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM==1050 || GCM==1051 ||GCM==1052){ /* BCCR */
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM==1061 || GCM==1062 ||GCM==1063){ /* INM-CM3 */
		GCM_R = 45;
		GCM_C = 72;
	}else if(GCM==1070 || GCM==1071 || GCM==1072){ /* IPSL */
		GCM_R = 72;
		GCM_C = 96;
	}else if(GCM==1080 || GCM==1081 || GCM==1082){ /* GFDL CM2-1 */
		GCM_R = 90;
		GCM_C = 144;
	}else if(GCM==1090 || GCM==1091 || GCM==1092){ /* HadCM3 */
		GCM_R = 73;
		GCM_C = 96;
	}else if(GCM==1100 || GCM==1101 || GCM==1102 || GCM==1103){ /* GISS AOM */
		GCM_R = 60;
		GCM_C = 90;
	}else if(GCM==1110 || GCM==1111 || GCM==1112){ /* CCCma-T63 */
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM==1120 || GCM==1121 || GCM==1122){ /* CSIRO */
		GCM_R = 96;
		GCM_C = 192;
	}else if(GCM==1130 || GCM==1131 || GCM==1132 || GCM==1133 || GCM==1134
			 || GCM==1135 || GCM==1136 || GCM==1137 || GCM==1138 || GCM==1139
			  || GCM==1140 || GCM==1141 || GCM==1142 || GCM==1143 || GCM==1144){ /* MRI */
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM==1150 ||GCM==1151 ||GCM==1152 ||GCM==1153 ||GCM==1154 ||
			GCM==1155 ||GCM==1156 ||GCM==1157 ||GCM==1158 ||GCM==1159){ /* MPI ECHAM */
		GCM_R = 96;
		GCM_C = 192;
	}else if(GCM==1160 || GCM==1161 || GCM==1162 || GCM==1163 || GCM==1164 || GCM==1165){ /* IAP */
		GCM_R = 60;
		GCM_C = 128;
	}else if(GCM==1170 || GCM==1171){ /* HadGEM */
		GCM_R = 145;
		GCM_C = 192;
	}else if(GCM==1180 || GCM==1181 || GCM==1182){ /* GFDL CM 2.0 */
		GCM_R = 90;
		GCM_C = 144;
	}else if(GCM==1190 || GCM==1191 || GCM==1192){ /* CNRM */
		GCM_R = 64;
		GCM_C = 128;
	}else if(GCM==1200 || GCM==1201 || GCM==1202 || GCM==1203 || GCM==1204
			 || GCM==1205 || GCM==1206 || GCM==1207 || GCM==1208 || GCM==1209
			 || GCM==1210 || GCM==1211 || GCM==1212 || GCM==1213 || GCM==1214){ /* CCC47 */
		GCM_R = 48;
		GCM_C = 96;
	}else if(GCM==1220 || GCM==1221 || GCM==1222 || GCM==1223 || GCM==1224
			 || GCM==1225 || GCM==1226 || GCM==1227 || GCM==1228 || GCM==1229
			  || GCM==1230 || GCM==1231 || GCM==1232 || GCM==1233 || GCM==1234
			   || GCM==1235 || GCM==1236 || GCM==1237 || GCM==1238 || GCM==1239){ /* NCAR CCSM */
		GCM_R = 128;
		GCM_C = 256;
	}else if(GCM==1240 || GCM==1241 || GCM==1242){ /* GISS E */
		GCM_R = 46;
		GCM_C = 72;
	}else if(GCM==1250 || GCM==1251 || GCM==1252 || GCM==1253){ /* GISS R */
		GCM_R = 46;
		GCM_C = 72;
	}else if(GCM==1260 || GCM==1261 || GCM==1262 || GCM==1263 || GCM==1264
			 || GCM==1265 || GCM==1266 || GCM==1267){ /* NCAR PCM */
		GCM_R = 64;
		GCM_C = 128;
	}
}

/****************************************************************/
void set_gcm_index(
	char s_case[]
){
	switch(GCM){
		case 1000:	strcpy(s_case,"MHA11_"); break;
		case 1001:	strcpy(s_case,"MHB11_"); break;

		case 1010:	strcpy(s_case,"MMA11_"); break;
		case 1011:	strcpy(s_case,"MMA12_"); break;
		case 1012:	strcpy(s_case,"MMA13_"); break;
		case 1013:	strcpy(s_case,"MMA21_"); break;
		case 1014:	strcpy(s_case,"MMA22_"); break;
		case 1015:	strcpy(s_case,"MMA23_"); break;
		case 1016:	strcpy(s_case,"MMA21_"); break;
		case 1017:	strcpy(s_case,"MMA22_"); break;
		case 1018:	strcpy(s_case,"MMA23_"); break;

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
		case 1102:	strcpy(s_case,"GSAA21_"); break;
		case 1103:	strcpy(s_case,"GSAB11_"); break;

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

		case 1260:	strcpy(s_case,"GSRA11_"); break;
		case 1261:	strcpy(s_case,"GSRA12_"); break;
		case 1262:	strcpy(s_case,"GSRA13_"); break;
		case 1263:	strcpy(s_case,"GSRA14_"); break;
		case 1264:	strcpy(s_case,"GSRA21_"); break;
		case 1265:	strcpy(s_case,"GSRA22_"); break;
		case 1266:	strcpy(s_case,"GSRA23_"); break;
		case 1267:	strcpy(s_case,"GSRA24_"); break;
		
		default:
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

