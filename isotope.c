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

/***** photosynthetic carbon isotope descrimination *****/
void photo_13c_frac(
	struct Grid *grid, 
	struct Loct *loct,
	struct Pchar *plant
){
	double aaa, bbb, bb3, bb4, fai;
	
	aaa = 4.4;
	bbb = 27.0;
	bb3 = 30.0;
	bb4 = -5.7;
	fai = 0.21;

	if(plant->phototype==3){
		plant->photo_13c_frac[grid->m] = aaa + (bbb-aaa)*plant->ci[grid->m]/loct->aCO2[grid->m];
	}else if(plant->phototype==4){
		plant->photo_13c_frac[grid->m] = aaa + (bb4+bb3*fai-aaa)*plant->ci[grid->m]/loct->aCO2[grid->m];
	}
}

/***** conversion from d13C to 12C/13C ratio *****/
double deltaTratio(
	double delta
){
	double ratio;
	
	ratio = (1.0 + delta/1000.0)*STCIR;

	return ratio;
}

/***** conversion from 12C/13C ratio to d13C *****/
double ratioTdelta(
	double ratio
){
	double delta;
	
	delta = (ratio/STCIR - 1.0)*1000.0;
	
	return delta;
}

/******* addition of two isotopically different substances, a and b *******/
double d13c_addition(
	double d13c_a, 
	double mass_a, 
	double d13c_b, 
	double mass_b
){
	double ratio_a, ratio_b;
	double c12_a, c13_a, c12_b, c13_b, c12_pro, c13_pro;
	double d13c_product;
	
	ratio_a = deltaTratio(d13c_a); 
	ratio_b = deltaTratio(d13c_b); 
	
	c12_a = mass_a/(1.0+ratio_a); c13_a = mass_a-c12_a;
	c12_b = mass_b/(1.0+ratio_b); c13_b = mass_b-c12_b;
	
	c12_pro = c12_a+c12_b;
	c13_pro = c13_a+c13_b;
	
	if(c12_pro>0.0){
		d13c_product = (c13_pro/c12_pro/STCIR - 1.0)*1000.0;
	}else{
		d13c_product = d13c_a;
	}
	
	return (d13c_product);
}

/******* addition of two isotopically different substances, a and b *******/
double d13c_addition3(
	double d13c_a, 
	double mass_a, 
	double d13c_b, 
	double mass_b, 
	double d13c_c, 
	double mass_c
){
	double d13c_ab, mass_ab, d13c_abc;

	mass_ab=mass_a+mass_b;
	d13c_ab=d13c_addition( d13c_a, mass_a, d13c_b, mass_b);
	
	d13c_abc=d13c_addition( d13c_ab, mass_ab, d13c_c, mass_c);
		
	return(d13c_abc);
}

/****** d13c of efflux CO2 from terrestrial ecosystems *******/
void d13c_efflux(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Flux *flux
){
	
	/* total ecosystem respiratory efflux of CO2 */
	flux->efflux_p=(flux->plant).rp[grid->m]+(flux->soil).rS[grid->m];

	(flux->c3).d13c_rpm[grid->m]=d13c_addition3((flux->c3).d13c_rfm[grid->m], (flux->c3).rfm[grid->m], 
					(flux->c3).d13c_rcm[grid->m], (flux->c3).rcm[grid->m], (flux->c3).d13c_rrm[grid->m], (flux->c3).rrm[grid->m]);
	(flux->c3).d13c_rpg[grid->m]=d13c_addition3((flux->c3).d13c_rfg[grid->m], (flux->c3).rfg[grid->m], 
					(flux->c3).d13c_rcg[grid->m], (flux->c3).rcg[grid->m], (flux->c3).d13c_rrg[grid->m], (flux->c3).rrg[grid->m]);
	(flux->c3).d13c_rp[grid->m]=d13c_addition((flux->c3).d13c_rpm[grid->m], (flux->c3).rpm[grid->m], (flux->c3).d13c_rpg[grid->m], (flux->c3).rpg[grid->m]);

	(flux->c4).d13c_rpm[grid->m]=d13c_addition3((flux->c4).d13c_rfm[grid->m], (flux->c4).rfm[grid->m], 
					(flux->c4).d13c_rcm[grid->m], (flux->c4).rcm[grid->m], (flux->c4).d13c_rrm[grid->m], (flux->c4).rrm[grid->m]);
	(flux->c4).d13c_rpg[grid->m]=d13c_addition3((flux->c4).d13c_rfg[grid->m], (flux->c4).rfg[grid->m], 
					(flux->c4).d13c_rcg[grid->m], (flux->c4).rcg[grid->m], (flux->c4).d13c_rrg[grid->m], (flux->c4).rrg[grid->m]);
	(flux->c4).d13c_rp[grid->m]=d13c_addition((flux->c4).d13c_rpm[grid->m], (flux->c4).rpm[grid->m], (flux->c4).d13c_rpg[grid->m], (flux->c4).rpg[grid->m]);
	
	(flux->plant).d13c_rp[grid->m]=d13c_addition((flux->c3).d13c_rp[grid->m], loct->C3ptn[grid->m]*(flux->c3).rp[grid->m], 
					(flux->c4).d13c_rp[grid->m], loct->C4ptn[grid->m]*(flux->c4).rp[grid->m]);

	(flux->soil).d13c_rS[grid->m]=d13c_addition((flux->soil).d13c_rl[grid->m], (flux->soil).rl[grid->m], (flux->soil).d13c_rh[grid->m], (flux->soil).rh[grid->m]);
	
	/* total ecosystem respiratory efflux of d13C - CO2 */
	flux->d13c_efflux_p=d13c_addition((flux->plant).d13c_rp[grid->m], (flux->plant).rp[grid->m], (flux->soil).d13c_rS[grid->m], (flux->soil).rS[grid->m]);
}


