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

/* photosynthetic carbon isotope descrimination *************************/
void photo_13c_frac(
	struct Grid *grid, 
	struct Loct *loct,
	struct Pchar *pchar
){
	double aaa, bbb, bb3, bb4, fai;
	
	/*
	 Farquhar GD, Ehleringer JR, Hubick KT (1989) 
	 Carbon isotope discrimination and photosynthesis. 
	 AnnuRevPlant PhysiolPlant Mol Biol 40:503-537
	*/
	aaa = 4.4;
	bbb = 27.0;
	bb3 = 30.0;
	bb4 = -5.7;
	fai = 0.21;

	if(pchar->phototype == 3){
		pchar->photo_13c_frac[grid->m] = aaa + (bbb-aaa)*pchar->ci[grid->m]
				/loct->aco2[grid->m];
	}else if(pchar->phototype == 4){
		pchar->photo_13c_frac[grid->m] = aaa + (bb4+bb3*fai-aaa)*pchar->ci[grid->m]
				/loct->aco2[grid->m];
	}
}

/* conversion from d13C to 12C/13C ratio ***************/
double deltaTratio(
	double delta
){
	double ratio;
	
	ratio = (1.0 + delta/1000.0)*STCIR;

	return ratio;
}

/* conversion from 12C/13C ratio to d13C ***************/
double ratioTdelta(
	double ratio
){
	double delta;
	
	delta = (ratio/STCIR - 1.0)*1000.0;
	
	return delta;
}

/* addition of two isotopically different substances, a and b ************/
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

/* addition of two isotopically different substances, a and b *************/
double d13c_addition3(
	double d13c_a, 
	double mass_a, 
	double d13c_b, 
	double mass_b, 
	double d13c_c, 
	double mass_c
){
	double d13c_ab, mass_ab, d13c_abc;

	mass_ab = mass_a+mass_b;
	d13c_ab = d13c_addition(d13c_a, mass_a,  d13c_b, mass_b);
	
	d13c_abc = d13c_addition(d13c_ab, mass_ab,  d13c_c, mass_c);
		
	return(d13c_abc);
}

/* d13c of efflux CO2 from terrestrial ecosystems *****************/
void f_cisotope_efflux(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Flux *flux
){
	double d13c_rr3, d13c_rr4, aaa;
	
	(flux->plant).d13c_gpp[grid->m] = d13c_addition((flux->c3).d13c_gpp[grid->m], loct->c3ptn[grid->m]*(flux->c3).gpp[grid->m], 
													(flux->c4).d13c_gpp[grid->m], loct->c4ptn[grid->m]*(flux->c4).gpp[grid->m]);
	aaa = loct->c3ptn[grid->m]*(flux->c3).gpp[grid->m] + loct->c4ptn[grid->m]*(flux->c4).gpp[grid->m];
	if(aaa > 0.0){
		(flux->plant).d14c_gpp[grid->m] = ((flux->c3).d14c_gpp[grid->m]*loct->c3ptn[grid->m]*(flux->c3).gpp[grid->m]+
									 (flux->c4).d14c_gpp[grid->m]*loct->c4ptn[grid->m]*(flux->c4).gpp[grid->m])/
									(loct->c3ptn[grid->m]*(flux->c3).gpp[grid->m] + loct->c4ptn[grid->m]*(flux->c4).gpp[grid->m]);
	}else{
		(flux->plant).d14c_gpp[grid->m] = 0.0;
	}
	
	(mass->plant).d13c_fol = d13c_addition((mass->c3).d13c_fol,loct->c3ptn[grid->m]*(mass->c3).fol, 
										   (mass->c4).d13c_fol,loct->c4ptn[grid->m]*(mass->c4).fol);
	(mass->plant).d13c_stm = d13c_addition((mass->c3).d13c_stm,loct->c3ptn[grid->m]*(mass->c3).stm, 
										   (mass->c4).d13c_stm,loct->c4ptn[grid->m]*(mass->c4).stm);
	(mass->plant).d13c_rot = d13c_addition((mass->c3).d13c_rot,loct->c3ptn[grid->m]*(mass->c3).rot, 
										   (mass->c4).d13c_rot,loct->c4ptn[grid->m]*(mass->c4).rot);
	
	(mass->plant).d14c_fol = ((mass->c3).d14c_fol*loct->c3ptn[grid->m]*(mass->c3).fol + 
								(mass->c4).d14c_fol*loct->c4ptn[grid->m]*(mass->c4).fol) / 
								((mass->c3).fol+(mass->c4).fol);
	(mass->plant).d14c_stm = ((mass->c3).d14c_stm*loct->c3ptn[grid->m]*(mass->c3).stm + 
							  (mass->c4).d14c_stm*loct->c4ptn[grid->m]*(mass->c4).stm) / 
							((mass->c3).stm+(mass->c4).stm);
	(mass->plant).d14c_rot = ((mass->c3).d14c_rot*loct->c3ptn[grid->m]*(mass->c3).rot + 
							  (mass->c4).d14c_rot*loct->c4ptn[grid->m]*(mass->c4).rot) / 
							((mass->c3).rot+(mass->c4).rot);
	
	(mass->plant).d13c_mfol[grid->m] = (mass->plant).d13c_fol;
	(mass->plant).d13c_mstm[grid->m] = (mass->plant).d13c_stm;
	(mass->plant).d13c_mrot[grid->m] = (mass->plant).d13c_rot;
	
	(mass->plant).d14c_mfol[grid->m] = (mass->plant).d14c_fol;
	(mass->plant).d14c_mstm[grid->m] = (mass->plant).d14c_stm;
	(mass->plant).d14c_mrot[grid->m] = (mass->plant).d14c_rot;
	
	/* total ecosystem respiratory efflux of CO2 **********************/
	flux->efflux_p = (flux->plant).ar[grid->m] + (flux->soil).hr[grid->m];

	/* total ecosystem respiratory efflux of d13C - CO2 */
	(flux->c3).d13c_arm[grid->m] = d13c_addition3((flux->c3).d13c_rfm[grid->m], (flux->c3).rfm[grid->m], 
												  (flux->c3).d13c_rcm[grid->m], (flux->c3).rcm[grid->m], 
												  (flux->c3).d13c_rrm[grid->m], (flux->c3).rrm[grid->m]);
	(flux->c3).d13c_arg[grid->m] = d13c_addition3((flux->c3).d13c_rfg[grid->m], (flux->c3).rfg[grid->m], 
												  (flux->c3).d13c_rcg[grid->m], (flux->c3).rcg[grid->m], 
												  (flux->c3).d13c_rrg[grid->m], (flux->c3).rrg[grid->m]);
	(flux->c3).d13c_ar[grid->m] = d13c_addition((flux->c3).d13c_arm[grid->m], (flux->c3).arm[grid->m], 
												(flux->c3).d13c_arg[grid->m], (flux->c3).arg[grid->m]);

	(flux->c4).d13c_arm[grid->m] = d13c_addition3((flux->c4).d13c_rfm[grid->m], (flux->c4).rfm[grid->m], 
												  (flux->c4).d13c_rcm[grid->m], (flux->c4).rcm[grid->m], 
												  (flux->c4).d13c_rrm[grid->m], (flux->c4).rrm[grid->m]);
	(flux->c4).d13c_arg[grid->m] = d13c_addition3((flux->c4).d13c_rfg[grid->m], (flux->c4).rfg[grid->m], 
												  (flux->c4).d13c_rcg[grid->m], (flux->c4).rcg[grid->m], 
												  (flux->c4).d13c_rrg[grid->m], (flux->c4).rrg[grid->m]);
	(flux->c4).d13c_ar[grid->m] = d13c_addition((flux->c4).d13c_arm[grid->m], (flux->c4).arm[grid->m], 
												(flux->c4).d13c_arg[grid->m], (flux->c4).arg[grid->m]);
	
	(flux->plant).d13c_arm[grid->m] = d13c_addition((flux->c3).d13c_arm[grid->m], loct->c3ptn[grid->m]*(flux->c3).arm[grid->m], 
													(flux->c4).d13c_arm[grid->m], loct->c4ptn[grid->m]*(flux->c4).arm[grid->m]);
	(flux->plant).d13c_arg[grid->m] = d13c_addition((flux->c3).d13c_arg[grid->m], loct->c3ptn[grid->m]*(flux->c3).arg[grid->m], 
													(flux->c4).d13c_arg[grid->m], loct->c4ptn[grid->m]*(flux->c4).arg[grid->m]);
	(flux->plant).d13c_ar[grid->m] = d13c_addition((flux->c3).d13c_ar[grid->m], loct->c3ptn[grid->m]*(flux->c3).ar[grid->m], 
												   (flux->c4).d13c_ar[grid->m], loct->c4ptn[grid->m]*(flux->c4).ar[grid->m]);
	
	(flux->soil).d13c_hr[grid->m] = d13c_addition((flux->soil).d13c_rl[grid->m], (flux->soil).rl[grid->m], 
												  (flux->soil).d13c_rh[grid->m], (flux->soil).rh[grid->m]);
	
	flux->d13c_er[grid->m] = flux->d13c_efflux_p = d13c_addition((flux->plant).d13c_ar[grid->m], (flux->plant).ar[grid->m], 
										(flux->soil).d13c_hr[grid->m], (flux->soil).hr[grid->m]);
	
	d13c_rr3 = d13c_addition( (flux->c3).d13c_rrm[grid->m], (flux->c3).rrm[grid->m], 
							 (flux->c3).d13c_rrg[grid->m], (flux->c3).rrg[grid->m]);
	d13c_rr4 = d13c_addition( (flux->c4).d13c_rrm[grid->m], (flux->c4).rrm[grid->m], 
							 (flux->c4).d13c_rrg[grid->m], (flux->c4).rrg[grid->m]);
	
	flux->d13c_sr[grid->m] = d13c_addition3(d13c_rr3, loct->c3ptn[grid->m]*((flux->c3).rrm[grid->m]+(flux->c3).rrg[grid->m]), 
											d13c_rr4, loct->c4ptn[grid->m]*((flux->c4).rrm[grid->m]+(flux->c4).rrg[grid->m]), 
											(flux->soil).d13c_hr[grid->m], (flux->soil).hr[grid->m]);
	
	flux->d14c_sr[grid->m] = ((mass->c3).d14c_rot*loct->c3ptn[grid->m]*((flux->c3).rrm[grid->m]+(flux->c3).rrg[grid->m])+
						(mass->c4).d14c_rot*loct->c4ptn[grid->m]*((flux->c4).rrm[grid->m]+(flux->c4).rrg[grid->m])+
						(mass->soil).d14c_ltr*(flux->soil).rl[grid->m] + (mass->soil).d14c_msl*(flux->soil).rh[grid->m]) / 
						(loct->c3ptn[grid->m]*((flux->c3).rrm[grid->m]+(flux->c3).rrg[grid->m])+
						loct->c4ptn[grid->m]*((flux->c4).rrm[grid->m]+(flux->c4).rrg[grid->m])+
						(flux->soil).rl[grid->m] + (flux->soil).rh[grid->m]);

	flux->d14c_er[grid->m] = ((mass->c3).d14c_fol*loct->c3ptn[grid->m]*((flux->c3).rfm[grid->m]+(flux->c3).rfg[grid->m])+
							  (mass->c3).d14c_stm*loct->c3ptn[grid->m]*((flux->c3).rcm[grid->m]+(flux->c3).rcg[grid->m])+
							  (mass->c3).d14c_rot*loct->c3ptn[grid->m]*((flux->c3).rrm[grid->m]+(flux->c3).rrg[grid->m])+
							  (mass->c4).d14c_fol*loct->c4ptn[grid->m]*((flux->c4).rfm[grid->m]+(flux->c4).rfg[grid->m])+
							  (mass->c4).d14c_stm*loct->c4ptn[grid->m]*((flux->c4).rcm[grid->m]+(flux->c4).rcg[grid->m])+
							  (mass->c4).d14c_rot*loct->c4ptn[grid->m]*((flux->c4).rrm[grid->m]+(flux->c4).rrg[grid->m])+
							  (mass->soil).d14c_ltr*(flux->soil).rl[grid->m] + (mass->soil).d14c_msl*(flux->soil).rh[grid->m]) / 
							(loct->c3ptn[grid->m]*((flux->c3).rfm[grid->m]+(flux->c3).rfg[grid->m])+
							 loct->c3ptn[grid->m]*((flux->c3).rcm[grid->m]+(flux->c3).rcg[grid->m])+
							 loct->c3ptn[grid->m]*((flux->c3).rrm[grid->m]+(flux->c3).rrg[grid->m])+
							 loct->c4ptn[grid->m]*((flux->c4).rfm[grid->m]+(flux->c4).rfg[grid->m])+
							 loct->c4ptn[grid->m]*((flux->c4).rcm[grid->m]+(flux->c4).rcg[grid->m])+
							 loct->c4ptn[grid->m]*((flux->c4).rrm[grid->m]+(flux->c4).rrg[grid->m])+
							 (flux->soil).rl[grid->m] + (flux->soil).rh[grid->m]);
}

/* d14C decay: added by A.Ito (2009/06/27) **********************/
double f_decay_14c(
	double init_val
){
	double end_val;
	double ratio, decayed;
	
	ratio = (init_val/1000.0 + 1.0) * 1.2;	/* denomination by 1.2 x 10^-10 */
	
	/* 14C decay for 1-month */
	decayed = ratio * exp(-0.0001245 * (1.0/12.0));
	
	end_val = (decayed/1.2 -1.0) * 1000.0;
	
	return end_val;
}

