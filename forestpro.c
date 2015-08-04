/*	VISIT: Vegetation Integrative SImulation Tool						*/
/*  Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/*  Carbon cycle, erosion, biomass burning, land-use change,			*/
/*  CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* biome dependent plant carbon exchange processes----evergreen and deciduous forest */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* forest processes *****************************************/
void forest_process(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pflx *flux, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	/* preparation */
	f_before_deal(grid, flux);

	switch(pchar->season[grid->m]){
		case 0: /* dormancy */
			noleafperiod(grid, loct, pchar, mass, flux); 
			break;
		case 1: /* vegetative growth */
			greenperiod(grid, loct, pchar, mass, flux); 
			reallocation_survival(grid, pchar, mass);
			break;
		case 2: /* new leaf emergence */
			leafemergence(grid, loct, pchar, mass, flux); 
			reallocation_survival(grid, pchar, mass);
			break;
		case 3: /* old leaf shedding */
			leaffall(grid, loct, pchar, mass, flux); 
			break;
	}
	
	/* post-processing */
	f_after_deal(grid, loct, pchar, mass, flux);
}

/* green period, while plants grow up actively *********************/
void greenperiod(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){						
	double nn;
	extern double MDN[ASTEP];
		
	nn = MDN[grid->m];

	/* litter */
	flux->lf[grid->m] = nn * flf(grid, pchar, mass);
	flux->lc[grid->m] = nn * flc(grid, pchar, mass);
	flux->lr[grid->m] = nn * flr(grid, pchar, mass);
	/* stable carbon isotope */
	flux->d13c_lf[grid->m] = mass->d13c_fol;
	flux->d13c_lc[grid->m] = mass->d13c_stm;
	flux->d13c_lr[grid->m] = mass->d13c_rot;

	/* litter fall */
	mass->fol -= flux->lf[grid->m];
	mass->lai[grid->m] = lai_mass(grid, mass, pchar);
	mass->stm -= flux->lc[grid->m];
	mass->rot -= flux->lr[grid->m];
	
	f_leaf_age(1, pchar, mass, -flux->lf[grid->m]);

	/* photosynthesis, gross primary production */
	flux->gpp[grid->m] = nn * fgpp(grid, loct, pchar, mass);
	/* stable carbon isotope */
	flux->d13c_gpp[grid->m] = loct->d13c_aco2[grid->m] - pchar->photo_13c_frac[grid->m];
	
	/* GPP by de Pury & Farquhar scheme */
	if(DF97==1){
		flux->gpp_df97[grid->m] = f_df97_gpp(1, grid, loct, pchar, mass);
	}
	
	/* maintenance respirations */
	flux->rfm[grid->m] = nn * frfm(grid, pchar, mass);
	flux->rcm[grid->m] = nn * frcm(grid, pchar, mass);
	flux->rrm[grid->m] = nn * frrm(grid, pchar, mass);
	flux->arm[grid->m] = flux->rfm[grid->m] + flux->rcm[grid->m] + flux->rrm[grid->m];
	/* stable carbon isotope */
	flux->d13c_rfm[grid->m] = mass->d13c_fol;
	flux->d13c_rcm[grid->m] = mass->d13c_stm;
	flux->d13c_rrm[grid->m] = mass->d13c_rot;
	
	/* tentative primary production */	
	if(DF97==1){
		flux->epp[grid->m] = flux->gpp_df97[grid->m] - flux->arm[grid->m];
	}else{
		flux->epp[grid->m] = flux->gpp[grid->m] - flux->arm[grid->m];
	}
	
	/* translocation of photosynthate */
	allocation(grid, loct, pchar, mass, flux);
	/* stable carbon isotope */
	flux->d13c_tpf[grid->m] = flux->d13c_gpp[grid->m];
	flux->d13c_tpc[grid->m] = flux->d13c_gpp[grid->m];
	flux->d13c_tpr[grid->m] = flux->d13c_gpp[grid->m];
	flux->d13c_tpp[grid->m] = flux->d13c_gpp[grid->m];
	
	if(flux->epp[grid->m]>0.0){
		/* growth construction respiration */
		flux->rfg[grid->m] = frfg(grid, pchar, flux);
		flux->rcg[grid->m] = frcg(grid, pchar, flux);
		flux->rrg[grid->m] = frrg(grid, pchar, flux);
	}else if(flux->epp[grid->m]<=0.0){
		flux->rfg[grid->m] = flux->rcg[grid->m] = flux->rrg[grid->m] = 0.0;
	}
	/* stable carbon isotope */
	flux->d13c_rfg[grid->m] = flux->d13c_tpf[grid->m];
	flux->d13c_rcg[grid->m] = flux->d13c_tpc[grid->m];
	flux->d13c_rrg[grid->m] = flux->d13c_tpr[grid->m];
	
	/* partitioning of photosynthate */
	mass->fol += (flux->tpf[grid->m] - flux->rfg[grid->m]);
	mass->stm += (flux->tpc[grid->m] - flux->rcg[grid->m]);
	mass->rot += (flux->tpr[grid->m] - flux->rrg[grid->m]);
	
	f_leaf_age(0, pchar, mass, (flux->tpf[grid->m] - flux->rfg[grid->m]));
	
	/* stable carbon isotope */
	if((flux->tpf[grid->m]-flux->rfg[grid->m]) > 0.0){
		mass->d13c_fol = d13c_addition(mass->d13c_fol, mass->fol, 
		flux->d13c_tpf[grid->m], (flux->tpf[grid->m]-flux->rfg[grid->m]));
	}else if((flux->tpf[grid->m]-flux->rfg[grid->m]) < 0.0){
		mass->d13c_fol = mass->d13c_fol;
	}
	if((flux->tpc[grid->m]-flux->rcg[grid->m]) > 0.0){
		mass->d13c_stm = d13c_addition(mass->d13c_stm, mass->stm, 
		flux->d13c_tpc[grid->m], (flux->tpc[grid->m]-flux->rcg[grid->m]));
	}else if((flux->tpf[grid->m]-flux->rfg[grid->m]) < 0.0){
		mass->d13c_stm = mass->d13c_stm;
	}
	if((flux->tpc[grid->m]-flux->rcg[grid->m]) > 0.0){
		mass->d13c_rot = d13c_addition(mass->d13c_rot, mass->rot, 
		flux->d13c_tpr[grid->m], (flux->tpr[grid->m]-flux->rrg[grid->m]));
	}else if((flux->tpr[grid->m]-flux->rrg[grid->m]) < 0.0){
		mass->d13c_rot = mass->d13c_rot;
	}
}

/* leaf shedding ****************************************************/
void leaffall(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){
	double nn;
	extern double MDN[12];
		
	nn = MDN[grid->m];

	/* litter */
	flux->lf[grid->m] = pchar->dcd * mass->fol; /* leaf-shedding */
	flux->lc[grid->m] = nn * flc(grid, pchar, mass);
	flux->lr[grid->m] = nn * flr(grid, pchar, mass);
	/* stable carbon isotope */
	flux->d13c_lf[grid->m] = mass->d13c_fol;
	flux->d13c_lc[grid->m] = mass->d13c_stm;
	flux->d13c_lr[grid->m] = mass->d13c_rot;

	/* litter fall */
	mass->fol -= flux->lf[grid->m];
	mass->lai[grid->m] = lai_mass(grid, mass, pchar);
	mass->stm -= flux->lc[grid->m];
	mass->rot -= flux->lr[grid->m];
	
	/* leaf age */
	f_leaf_age(1, pchar, mass, -flux->lf[grid->m]);

	/* photosynthesis, gross primary production */
	flux->gpp[grid->m] = nn * fgpp(grid, loct, pchar, mass);
	/* stable carbon isotope */
	flux->d13c_gpp[grid->m] = loct->d13c_aco2[grid->m] - pchar->photo_13c_frac[grid->m];
	
	/* GPP by de Pury & Farquhar scheme */
	if(DF97 == 1){
		flux->gpp_df97[grid->m] = f_df97_gpp(1, grid, loct, pchar, mass);
	}

	/* maintenance respirations */
	flux->rfm[grid->m] = nn*frfm(grid, pchar, mass);
	flux->rcm[grid->m] = nn*frcm(grid, pchar, mass);
	flux->rrm[grid->m] = nn*frrm(grid, pchar, mass);
	flux->arm[grid->m] = flux->rfm[grid->m] + flux->rcm[grid->m] + flux->rrm[grid->m];
	/* stable carbon isotope */
	flux->d13c_rfm[grid->m] = mass->d13c_fol;
	flux->d13c_rcm[grid->m] = mass->d13c_stm;
	flux->d13c_rrm[grid->m] = mass->d13c_rot;
	
	/* tentative primary production */	
	if(DF97==1){
		flux->epp[grid->m] = flux->gpp_df97[grid->m] - flux->arm[grid->m];
	}else{
		flux->epp[grid->m] = flux->gpp[grid->m] - flux->arm[grid->m];
	}
	
	/* translocation of photosynthate */
	allocation(grid, loct, pchar, mass, flux);
	/* stable carbon isotope */
	flux->d13c_tpf[grid->m] = flux->d13c_gpp[grid->m];
	flux->d13c_tpc[grid->m] = flux->d13c_gpp[grid->m];
	flux->d13c_tpr[grid->m] = flux->d13c_gpp[grid->m];
	flux->d13c_tpp[grid->m] = flux->d13c_gpp[grid->m];

	if(flux->epp[grid->m] > 0.0){
		/* growth construction respiration */
		flux->rfg[grid->m] = frfg(grid, pchar, flux);
		flux->rcg[grid->m] = frcg(grid, pchar, flux);
		flux->rrg[grid->m] = frrg(grid, pchar, flux);
	}else if(flux->epp[grid->m] <= 0.0){
		flux->rfg[grid->m] = flux->rcg[grid->m] = flux->rrg[grid->m] = 0.0;
	}
	/* stable carbon isotope */
	flux->d13c_rfg[grid->m] = flux->d13c_tpf[grid->m];
	flux->d13c_rcg[grid->m] = flux->d13c_tpc[grid->m];
	flux->d13c_rrg[grid->m] = flux->d13c_tpr[grid->m];
	
	/* partitioning of photosynthate */
	mass->fol += (flux->tpf[grid->m]-flux->rfg[grid->m]);
	mass->stm += (flux->tpc[grid->m]-flux->rcg[grid->m]);
	mass->rot += (flux->tpr[grid->m]-flux->rrg[grid->m]);
	
	f_leaf_age(0, pchar, mass, (flux->tpf[grid->m] - flux->rfg[grid->m]));
	
	/* stable carbon isotope */
	if((flux->tpf[grid->m]-flux->rfg[grid->m])>0.0){
		mass->d13c_fol = d13c_addition(mass->d13c_fol, mass->fol, 
					flux->d13c_tpf[grid->m], (flux->tpf[grid->m]-flux->rfg[grid->m]));
	}else if((flux->tpf[grid->m]-flux->rfg[grid->m])<0.0){
		mass->d13c_fol = mass->d13c_fol;
	}
	if((flux->tpc[grid->m]-flux->rcg[grid->m])>0.0){
		mass->d13c_stm = d13c_addition(mass->d13c_stm, mass->stm, 
					flux->d13c_tpc[grid->m], (flux->tpc[grid->m]-flux->rcg[grid->m]));
	}else if((flux->tpf[grid->m]-flux->rfg[grid->m])<0.0){
		mass->d13c_stm = mass->d13c_stm;
	}
	if((flux->tpc[grid->m]-flux->rcg[grid->m])>0.0){
		mass->d13c_rot = d13c_addition(mass->d13c_rot, mass->rot, 
					flux->d13c_tpr[grid->m], (flux->tpr[grid->m]-flux->rrg[grid->m]));
	}else if((flux->tpr[grid->m]-flux->rrg[grid->m])<0.0){
		mass->d13c_rot = mass->d13c_rot;
	}
}

/* new leaf emergence *****************************************/
void leafemergence(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){
	double nn, aaa, bbb, emerge;
	extern double MDN[12];
		
	nn = MDN[grid->m];

	/* leaf emergence */
	bbb = (pchar->opt_lai[grid->m]>2.0)?pchar->opt_lai[grid->m]:2.0;
	emerge = (bbb - mass->lai[grid->m])*100.0*2.0/2.2/pchar->sla;
	aaa = mass->stm + mass->rot;
	
	switch(grid->veg_olson){
		case 3: case 4: case 5: 
			emerge = (emerge>aaa*0.1)?aaa*0.1:emerge;
			break;
		case 10: case 12:
			emerge = (emerge>aaa*0.06)?aaa*0.06:emerge;
			break;
		case 14: case 15: case 16: case 17: case 18: case 19: case 27: case 28:
		case 29: case 30: case 31: case 32:
			emerge = (emerge>aaa*0.3)?aaa*0.3:emerge;
			break;
	}
	
	if(grid->veg_olson == 20) emerge=(emerge<1.0)?emerge:1.0;
	
	if(aaa>0.0){
		mass->fol += emerge;
		mass->stm -= emerge*mass->stm/aaa;
		mass->rot -= emerge*mass->rot/aaa;
	}
	
	f_leaf_age(1, pchar, mass, emerge);

	/* litter */
	flux->lf[grid->m] = nn*flf(grid, pchar, mass);
	flux->lc[grid->m] = nn*flc(grid, pchar, mass);
	flux->lr[grid->m] = nn*flr(grid, pchar, mass);
	/* stable carbon isotope */
	flux->d13c_lf[grid->m] = mass->d13c_fol;
	flux->d13c_lc[grid->m] = mass->d13c_stm;
	flux->d13c_lr[grid->m] = mass->d13c_rot;

	f_leaf_age(0, pchar, mass, -flux->lf[grid->m]);

	/* litter fall */
	mass->fol -= flux->lf[grid->m];
	mass->lai[grid->m] = lai_mass(grid, mass, pchar);
	mass->stm -= flux->lc[grid->m];
	mass->rot -= flux->lr[grid->m];

	/* photosynthesis, gross primary production */
	flux->gpp[grid->m] = nn*fgpp(grid, loct, pchar, mass);
	/* stable carbon isotope */
	flux->d13c_gpp[grid->m] = loct->d13c_aco2[grid->m]-pchar->photo_13c_frac[grid->m];
	
	/* GPP by de Pury & Farquhar scheme */
	if(DF97==1){
		flux->gpp_df97[grid->m] = f_df97_gpp(1, grid, loct, pchar, mass);
	}

	/* maintenance respirations */
	flux->rfm[grid->m] = nn * frfm(grid, pchar, mass);
	flux->rcm[grid->m] = nn * frcm(grid, pchar, mass);
	flux->rrm[grid->m] = nn * frrm(grid, pchar, mass);
	flux->arm[grid->m] = flux->rfm[grid->m]+flux->rcm[grid->m]+flux->rrm[grid->m];
	/* stable carbon isotope */
	flux->d13c_rfm[grid->m] = mass->d13c_fol;
	flux->d13c_rcm[grid->m] = mass->d13c_stm;
	flux->d13c_rrm[grid->m] = mass->d13c_rot;
	
	/* tentative primary production */	
	if(DF97==1){
		flux->epp[grid->m] = flux->gpp_df97[grid->m] - flux->arm[grid->m];
	}else{
		flux->epp[grid->m] = flux->gpp[grid->m] - flux->arm[grid->m];
	}
	
	/* translocation of photosynthate */
	allocation(grid, loct, pchar, mass, flux);
	/* stable carbon isotope */
	flux->d13c_tpf[grid->m] = flux->d13c_gpp[grid->m];
	flux->d13c_tpc[grid->m] = flux->d13c_gpp[grid->m];
	flux->d13c_tpr[grid->m] = flux->d13c_gpp[grid->m];
	flux->d13c_tpp[grid->m] = flux->d13c_gpp[grid->m];

	if(flux->epp[grid->m] > 0.0){
		/* growth construction respiration */
		flux->rfg[grid->m] = frfg(grid, pchar, flux);
		flux->rcg[grid->m] = frcg(grid, pchar, flux);
		flux->rrg[grid->m] = frrg(grid, pchar, flux);
	}else if(flux->epp[grid->m] <= 0.0){
		flux->rfg[grid->m] = flux->rcg[grid->m] = flux->rrg[grid->m] = 0.0;
	}
	/* stable carbon isotope */
	flux->d13c_rfg[grid->m] = flux->d13c_tpf[grid->m];
	flux->d13c_rcg[grid->m] = flux->d13c_tpc[grid->m];
	flux->d13c_rrg[grid->m] = flux->d13c_tpr[grid->m];
	
	/* partitioning of photosynthate */
	mass->fol += (flux->tpf[grid->m] - flux->rfg[grid->m]);
	mass->stm += (flux->tpc[grid->m] - flux->rcg[grid->m]);
	mass->rot += (flux->tpr[grid->m] - flux->rrg[grid->m]);
	
	f_leaf_age(0, pchar, mass, (flux->tpf[grid->m] - flux->rfg[grid->m]));
	
	/* stable carbon isotope */
	if((flux->tpf[grid->m] - flux->rfg[grid->m])>0.0){
		mass->d13c_fol = d13c_addition(mass->d13c_fol, mass->fol, flux->d13c_tpf[grid->m], 
					(flux->tpf[grid->m] - flux->rfg[grid->m]));
	}else if((flux->tpf[grid->m] - flux->rfg[grid->m])<0.0){
		mass->d13c_fol = mass->d13c_fol;
	}
	if((flux->tpc[grid->m] - flux->rcg[grid->m])>0.0){
		mass->d13c_stm = d13c_addition(mass->d13c_stm, mass->stm, flux->d13c_tpc[grid->m], 
					(flux->tpc[grid->m] - flux->rcg[grid->m]));
	}else if((flux->tpf[grid->m] - flux->rfg[grid->m])<0.0){
		mass->d13c_stm = mass->d13c_stm;
	}
	if((flux->tpc[grid->m] - flux->rcg[grid->m])>0.0){
		mass->d13c_rot = d13c_addition(mass->d13c_rot, mass->rot, flux->d13c_tpr[grid->m], 
					(flux->tpr[grid->m] - flux->rrg[grid->m]));
	}else if((flux->tpr[grid->m] - flux->rrg[grid->m])<0.0){
		mass->d13c_rot = mass->d13c_rot;
	}
}

/* dormancy period ***************************************************/
void noleafperiod(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){
	double nn;
	extern double MDN[ASTEP];
		
	nn = MDN[grid->m];

	/* litter */
	flux->lf[grid->m] = nn * flf(grid, pchar, mass);
	flux->lc[grid->m] = nn * flc(grid, pchar, mass);
	flux->lr[grid->m] = nn * flr(grid, pchar, mass);
	/* stable carbon isotope */
	flux->d13c_lf[grid->m] = mass->d13c_fol;
	flux->d13c_lc[grid->m] = mass->d13c_stm;
	flux->d13c_lr[grid->m] = mass->d13c_rot;

	/* litter fall */
	mass->fol -= flux->lf[grid->m];
	mass->lai[grid->m] = lai_mass(grid, mass, pchar);
	mass->stm -= flux->lc[grid->m];
	mass->rot -= flux->lr[grid->m];
	
	f_leaf_age(1, pchar, mass, -flux->lf[grid->m]);

	/* photosynthesis, gross primary production */
	flux->gpp[grid->m] = nn*fgpp(grid, loct, pchar, mass);
	/* stable carbon isotope: discrmination */
	flux->d13c_gpp[grid->m] = loct->d13c_aco2[grid->m] - pchar->photo_13c_frac[grid->m];
	
	/* GPP by de Pury & Farquhar scheme */
	if(DF97==1){
		flux->gpp_df97[grid->m] = f_df97_gpp(1, grid, loct, pchar, mass);
	}

	/* maintenance respirations */
	flux->rfm[grid->m] = nn * frfm(grid, pchar, mass);
	flux->rcm[grid->m] = nn * frcm(grid, pchar, mass);
	flux->rrm[grid->m] = nn * frrm(grid, pchar, mass);
	flux->arm[grid->m] = flux->rfm[grid->m] + flux->rcm[grid->m] + flux->rrm[grid->m];
	/* stable carbon isotope */
	flux->d13c_rfm[grid->m] = mass->d13c_fol;
	flux->d13c_rcm[grid->m] = mass->d13c_stm;
	flux->d13c_rrm[grid->m] = mass->d13c_rot;
	
	/* tentative primary production */	
	if(DF97==1){
		flux->epp[grid->m] = flux->gpp_df97[grid->m] - flux->arm[grid->m];
	}else{
		flux->epp[grid->m] = flux->gpp[grid->m] - flux->arm[grid->m];
	}
	
	/* translocation of photosynthate */
	allocation(grid, loct, pchar, mass, flux);
	/* stable carbon isotope */
	flux->d13c_tpf[grid->m] = flux->d13c_gpp[grid->m];
	flux->d13c_tpc[grid->m] = flux->d13c_gpp[grid->m];
	flux->d13c_tpr[grid->m] = flux->d13c_gpp[grid->m];
	flux->d13c_tpp[grid->m] = flux->d13c_gpp[grid->m];

	if(flux->epp[grid->m]>0.0){
		/* growth construction respiration */
		flux->rfg[grid->m] = frfg(grid, pchar, flux);
		flux->rcg[grid->m] = frcg(grid, pchar, flux);
		flux->rrg[grid->m] = frrg(grid, pchar, flux);
	}else if(flux->epp[grid->m]<=0.0){
		flux->rfg[grid->m] = flux->rcg[grid->m] = flux->rrg[grid->m] = 0.0;
	}
	/* stable carbon isotope */
	flux->d13c_rfg[grid->m] = flux->d13c_tpf[grid->m];
	flux->d13c_rcg[grid->m] = flux->d13c_tpc[grid->m];
	flux->d13c_rrg[grid->m] = flux->d13c_tpr[grid->m];
	
	/* partitioning of photosynthate */
	mass->fol += (flux->tpf[grid->m] - flux->rfg[grid->m]);
	mass->stm += (flux->tpc[grid->m] - flux->rcg[grid->m]);
	mass->rot += (flux->tpr[grid->m] - flux->rrg[grid->m]);
	
	f_leaf_age(0, pchar, mass, (flux->tpf[grid->m] - flux->rfg[grid->m]));
	
	/* stable carbon isotope */
	if((flux->tpf[grid->m] - flux->rfg[grid->m])>0.0){
		mass->d13c_fol = d13c_addition(mass->d13c_fol, mass->fol, flux->d13c_tpf[grid->m], 
					(flux->tpf[grid->m]-flux->rfg[grid->m]));
	}else if((flux->tpf[grid->m] - flux->rfg[grid->m])<0.0){
		mass->d13c_fol = mass->d13c_fol;
	}
	if((flux->tpc[grid->m] - flux->rcg[grid->m])>0.0){
		mass->d13c_stm = d13c_addition(mass->d13c_stm, mass->stm, flux->d13c_tpc[grid->m], 
					(flux->tpc[grid->m]-flux->rcg[grid->m]));
	}else if((flux->tpf[grid->m]-flux->rfg[grid->m])<0.0){
		mass->d13c_stm = mass->d13c_stm;
	}
	if((flux->tpc[grid->m] - flux->rcg[grid->m])>0.0){
		mass->d13c_rot = d13c_addition(mass->d13c_rot, mass->rot, flux->d13c_tpr[grid->m], 
					(flux->tpr[grid->m]-flux->rrg[grid->m]));
	}else if((flux->tpr[grid->m] - flux->rrg[grid->m])<0.0){
		mass->d13c_rot = mass->d13c_rot;
	}
}
