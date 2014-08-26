/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* biome dependent plant carbon exchange processes*/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* crop processes, e.g. grains, root crops, and pcharetables ***********/
void agri_process(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pflx *flux, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	/* preparation */	
	f_before_deal(grid, flux);
	
	switch(pchar->season[grid->m]){
		case 0: /* winter interval period */
			interval(grid, loct, pchar, mass, flux); 
			break;
		case 1: /* vegetative growth period */
			greenperiod(grid, loct, pchar, mass, flux);	 
			reallocation_survival(grid, pchar, mass);
			break;
		case 2: /* planting new crop */
			planting(grid, loct, pchar, mass, flux); 
			break;
		case 3: /* harvest */
			harvesting(grid, loct, pchar, mass, flux); 
			break;
	}
	
	/* settlement */
	f_after_deal(grid, pchar, mass, flux);
}

/* planting of new crops **********************************************/
void planting(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){
	double nn, aaa, bbb, emerge;
	extern double MDN[12];
	
	nn = MDN[grid->m];
	
	/* leaf emergence: perennial crops */
	aaa = mass->stm * 0.5;
	bbb = mass->rot * 0.5;
	emerge = aaa + bbb;
	
	if(emerge > 0.5){
		aaa = 0.5 / emerge * aaa;
		bbb = 0.5 / emerge * bbb;
		emerge = aaa + bbb;
	}
	mass->fol += emerge;
	mass->stm -= aaa;
	mass->rot -= bbb;
	
	f_leaf_age(1, pchar, mass, emerge);

	/* annual crops */
	mass->fol += 0.4;
	mass->stm += 0.3;
	mass->rot += 0.3;
	flux->hvst[grid->m] = 1.0;
	
	f_leaf_age(0, pchar, mass, 0.4);

	/* litter */
	flux->lf[grid->m] = nn*flf(grid, pchar, mass);
	flux->lc[grid->m] = nn*flc(grid, pchar, mass);
	flux->lr[grid->m] = nn*flr(grid, pchar, mass);
	/* stable carbon isotope */
	flux->d13c_lf[grid->m] = mass->d13c_fol;
	flux->d13c_lc[grid->m] = mass->d13c_stm;
	flux->d13c_lr[grid->m] = mass->d13c_rot;

	/* litter fall */
	mass->fol -= flux->lf[grid->m];
	mass->lai[grid->m] = lai_mass(grid, mass, pchar);
	mass->stm -= flux->lc[grid->m];
	mass->rot -= flux->lr[grid->m];
	
	f_leaf_age(0, pchar, mass, -flux->lf[grid->m]);

	/* photosynthesis, gross primary production */
	flux->gpp[grid->m] = nn * fgpp(grid, loct, pchar, mass);
	/* stable carbon isotope */
	flux->d13c_gpp[grid->m] = loct->d13c_aco2[grid->m] - pchar->photo_13c_frac[grid->m];
	
	/* GPP by de Pury & Farquhar scheme */
	if(DF97==1){
		flux->gpp_df97[grid->m] = f_df97_gpp(1, grid, loct, pchar, mass);
	}

	/* maintenance respirations */
	flux->rfm[grid->m] = nn*frfm(grid, pchar, mass);
	flux->rcm[grid->m] = nn*frcm(grid, pchar, mass);
	flux->rrm[grid->m] = nn*frrm(grid, pchar, mass);
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
	allocation(grid, pchar, mass, flux);
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
	
	f_leaf_age(0, pchar, mass, (flux->tpf[grid->m]-flux->rfg[grid->m]));
	
	/* stable carbon isotope */
	if((flux->tpf[grid->m] - flux->rfg[grid->m]) > 0.0){
		mass->d13c_fol = d13c_addition(mass->d13c_fol, mass->fol, 
			flux->d13c_tpf[grid->m], (flux->tpf[grid->m] - flux->rfg[grid->m]));
	}else if((flux->tpf[grid->m]-flux->rfg[grid->m])<0.0){
		mass->d13c_fol = mass->d13c_fol;
	}
	if((flux->tpc[grid->m] - flux->rcg[grid->m]) > 0.0){
		mass->d13c_stm = d13c_addition(mass->d13c_stm, mass->stm, 
			flux->d13c_tpc[grid->m], (flux->tpc[grid->m] - flux->rcg[grid->m]));
	}else if((flux->tpf[grid->m]-flux->rfg[grid->m])<0.0){
		mass->d13c_stm = mass->d13c_stm;
	}
	if((flux->tpc[grid->m] - flux->rcg[grid->m]) > 0.0){
		mass->d13c_rot = d13c_addition(mass->d13c_rot, mass->rot, 
			flux->d13c_tpr[grid->m], (flux->tpr[grid->m] - flux->rrg[grid->m]));
	}else if((flux->tpr[grid->m]-flux->rrg[grid->m])<0.0){
		mass->d13c_rot = mass->d13c_rot;
	}
}

/* harvest grain, leaves, and fruits ***************************************/
void harvesting(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){
	double hvst_index, clear, nn;
	extern double MDN[ASTEP];
	
	nn = MDN[grid->m];
		
	/** harvest of crops **/
	hvst_index = -0.45; /* harvest index -> 45% of biomass */
	flux->hvst[grid->m] = hvst_index*(mass->fol + mass->stm + mass->rot);

	mass->fol += (hvst_index * mass->fol);
	mass->stm += (hvst_index * mass->stm);
	mass->rot += (hvst_index * mass->rot);
	
	f_leaf_age(1, pchar, mass, hvst_index*mass->fol);

	/* litter */
	clear = 0.9 + hvst_index;
	flux->lf[grid->m] = clear * mass->fol;
	flux->lc[grid->m] = clear * mass->stm;
	flux->lr[grid->m] = clear * mass->rot;
	/* stable carbon isotope */
	flux->d13c_lf[grid->m] = mass->d13c_fol;
	flux->d13c_lc[grid->m] = mass->d13c_stm;
	flux->d13c_lr[grid->m] = mass->d13c_rot;

	/* litter fall */
	mass->fol -= flux->lf[grid->m];
	mass->lai[grid->m] = lai_mass(grid, mass, pchar);
	mass->stm -= flux->lc[grid->m];
	mass->rot -= flux->lr[grid->m];
	
	f_leaf_age(0, pchar, mass, -flux->lf[grid->m]);

	/* photosynthesis, gross primary production */
	flux->gpp[grid->m] = nn * fgpp(grid, loct, pchar, mass);
	/* stable carbon isotope */
	flux->d13c_gpp[grid->m] = loct->d13c_aco2[grid->m] - pchar->photo_13c_frac[grid->m];
	
	/* GPP by de Pury & Farquhar scheme */
	if(DF97==1){
		flux->gpp_df97[grid->m] = f_df97_gpp(1, grid, loct, pchar, mass);
	}

	/* maintenance respirations */
	flux->rfm[grid->m] = nn*frfm(grid, pchar, mass);
	flux->rcm[grid->m] = nn*frcm(grid, pchar, mass);
	flux->rrm[grid->m] = nn*frrm(grid, pchar, mass);
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
	allocation(grid, pchar, mass, flux);
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
	}else if(flux->epp[grid->m] <= 0.0){
		/* non productive period */
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
	
	f_leaf_age(0, pchar, mass, flux->tpf[grid->m]-flux->rfg[grid->m]);
	
	/* stable carbon isotope */
	if((flux->tpf[grid->m]-flux->rfg[grid->m])>0.0){
		mass->d13c_fol = d13c_addition(mass->d13c_fol, mass->fol, 
			flux->d13c_tpf[grid->m], (flux->tpf[grid->m] - flux->rfg[grid->m]));
	}else if((flux->tpf[grid->m]-flux->rfg[grid->m])<0.0){
		mass->d13c_fol = mass->d13c_fol;
	}
	if((flux->tpc[grid->m]-flux->rcg[grid->m])>0.0){
		mass->d13c_stm = d13c_addition(mass->d13c_stm, mass->stm, 
			flux->d13c_tpc[grid->m], (flux->tpc[grid->m] - flux->rcg[grid->m]));
	}else if((flux->tpf[grid->m]-flux->rfg[grid->m])<0.0){
		mass->d13c_stm = mass->d13c_stm;
	}
	if((flux->tpc[grid->m]-flux->rcg[grid->m])>0.0){
		mass->d13c_rot = d13c_addition(mass->d13c_rot, mass->rot, 
			flux->d13c_tpr[grid->m], (flux->tpr[grid->m] - flux->rrg[grid->m]));
	}else if((flux->tpr[grid->m]-flux->rrg[grid->m])<0.0){
		mass->d13c_rot = mass->d13c_rot;
	}
}

/* interval period, fallowing lands ********************************************/
void interval(
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
	flux->lf[grid->m] = nn*flf(grid, pchar, mass);
	flux->lc[grid->m] = nn*flc(grid, pchar, mass);
	flux->lr[grid->m] = nn*flr(grid, pchar, mass);
	/* stable carbon isotope */
	flux->d13c_lf[grid->m] = mass->d13c_fol;
	flux->d13c_lc[grid->m] = mass->d13c_stm;
	flux->d13c_lr[grid->m] = mass->d13c_rot;

	/* litter fall */
	mass->fol -= flux->lf[grid->m];
	mass->stm -= flux->lc[grid->m];
	mass->rot -= flux->lr[grid->m];
	
	f_leaf_age(1, pchar, mass, -flux->lf[grid->m]);

	/* photosynthesis, gross primary production */
	flux->gpp[grid->m] = nn*fgpp(grid, loct, pchar, mass);
	/* stable carbon isotope */
	flux->d13c_gpp[grid->m] = loct->d13c_aco2[grid->m] - pchar->photo_13c_frac[grid->m];
	
	/* GPP by de Pury & Farquhar scheme */
	if(DF97==1){
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
	allocation(grid, pchar, mass, flux);
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
	
	f_leaf_age(0, pchar, mass, (flux->tpf[grid->m]-flux->rfg[grid->m]));
	
	/* stable carbon isotope */
	if((flux->tpf[grid->m]-flux->rfg[grid->m])>0.0){
		mass->d13c_fol = d13c_addition(mass->d13c_fol, mass->fol, 
			flux->d13c_tpf[grid->m], (flux->tpf[grid->m]-flux->rfg[grid->m]));
	}else if((flux->tpf[grid->m]-flux->rfg[grid->m]) < 0.0){
		mass->d13c_fol = mass->d13c_fol;
	}
	if((flux->tpc[grid->m]-flux->rcg[grid->m])>0.0){
		mass->d13c_stm = d13c_addition(mass->d13c_stm, mass->stm, 
			flux->d13c_tpc[grid->m], (flux->tpc[grid->m]-flux->rcg[grid->m]));
	}else if((flux->tpf[grid->m]-flux->rfg[grid->m]) < 0.0){
		mass->d13c_stm = mass->d13c_stm;
	}
	if((flux->tpc[grid->m]-flux->rcg[grid->m])>0.0){
		mass->d13c_rot = d13c_addition(mass->d13c_rot, mass->rot, 
			flux->d13c_tpr[grid->m], (flux->tpr[grid->m]-flux->rrg[grid->m]));
	}else if((flux->tpr[grid->m]-flux->rrg[grid->m]) < 0.0){
		mass->d13c_rot = mass->d13c_rot;
	}
}


