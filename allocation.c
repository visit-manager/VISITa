/*	VISIT: Vegetation Integrative SImulation Tool						*/
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

extern short DF97;

/* allocation of photosynthate ********************************************/
void allocation(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pmas *mass, 
	struct Pflx *flux
){
						
	double aaa, bbb, cc1, ccc, ddd;
	double alloc_f, alloc_c, alloc_r; /* allocation ratios */
						
	if(flux->epp[grid->m] > 0.0){ /* during growing-period */
		/* allocation ratios of EPP */
		if(mass->lai[grid->m] > pchar->opt_lai[grid->m]){ 
			/* if holding LAI is greater than the optimam one */
			/**** avoid too high LAI under good condition: 2009/04/29 A.Ito ****/
			aaa = mass->lai[grid->m] - pchar->opt_lai[grid->m];
			bbb = 1.0 - 0.5*aaa;
			bbb = (bbb>0.0)?bbb:0.0;
			
			/* formar: bbb = 1.0; */
			
			ccc = pchar->alloc_ass*bbb;
			alloc_f = ccc;
			alloc_c = (1.0 - ccc)*pchar->alloc_abg;
			alloc_r = (1.0 - ccc)*(1.0-pchar->alloc_abg);
		}else{     /* if(mass->lai[grid->m] <= pchar->opt_lai[grid->m]) */
			/* if holding LAI is smaller than the optimam one */
			aaa = (pchar->opt_lai[grid->m]-mass->lai[grid->m])*100.0*2.0/2.2/pchar->sla;
			bbb = flux->epp[grid->m]*pchar->alloc_ass;
			
			/* allocate photosyntahte to foliage to attain the optimum one */
			if(aaa <= bbb){
				ccc = pchar->alloc_ass;
				
			}else if(aaa>bbb){
				cc1 = aaa/bbb;
				
				/* maximum allocation ratio to foliage is 40 %: 2008/08/25 by A.Ito */
				/* ccc = ((pchar->alloc_ass*cc1)<0.4)?(pchar->alloc_ass*cc1):0.4; */
				ccc = ((pchar->alloc_ass*cc1)<0.5)?(pchar->alloc_ass*cc1):0.5; 
			}

			/* ccc = pchar->alloc_ass*1.0; */
			
			if(pchar->season[grid->m] == 0){
			/*	ddd = 0.5; */ /* 2008/08/25 by A.Ito */
				ddd = 0.7; /* 2008/08/25 by A.Ito */
			}else{
				ddd = 1.0;
			}
			
			alloc_f = ccc * ddd;
			alloc_c = (1.0 - ccc*ddd) * pchar->alloc_abg;
			alloc_r = (1.0 - ccc*ddd) * (1.0 - pchar->alloc_abg);
		}
				
		/* monthly translocation fluxes */
		flux->tpp[grid->m] = (alloc_f + alloc_c + alloc_r) * flux->epp[grid->m];
		flux->tpf[grid->m] = alloc_f * flux->epp[grid->m];
		flux->tpc[grid->m] = alloc_c * flux->epp[grid->m];
		flux->tpr[grid->m] = alloc_r * flux->epp[grid->m];
		
	}else{ /* during NON growing-period */   /*  if(flux->epp[grid->m] <= 0.0) */
		/* allocation ratios of GPP, not EPP */
		alloc_f = pchar->alloc_ass;
		alloc_c = (1.0 - pchar->alloc_ass) * pchar->alloc_abg;
		alloc_r = (1.0 - pchar->alloc_ass) * (1.0 - pchar->alloc_abg);
				
		/* monthly translocation fluxes */
		flux->tpp[grid->m] = (alloc_f + alloc_c + alloc_r)*flux->epp[grid->m];
		
		if(DF97==1){
			flux->tpf[grid->m] = alloc_f*flux->gpp_df97[grid->m] - flux->rfm[grid->m];
			flux->tpc[grid->m] = alloc_c*flux->gpp_df97[grid->m] - flux->rcm[grid->m];
			flux->tpr[grid->m] = alloc_r*flux->gpp_df97[grid->m] - flux->rrm[grid->m];
		}else{
			flux->tpf[grid->m] = alloc_f*flux->gpp[grid->m] - flux->rfm[grid->m];
			flux->tpc[grid->m] = alloc_c*flux->gpp[grid->m] - flux->rcm[grid->m];
			flux->tpr[grid->m] = alloc_r*flux->gpp[grid->m] - flux->rrm[grid->m];
		}
	}
	
	/* monthly partitioning ratios */
	pchar->malloc_f[grid->m] = alloc_f;
	pchar->malloc_c[grid->m] = alloc_c;
	pchar->malloc_r[grid->m] = alloc_r;
}

/* plant-internal mass re-adjustment **********************************/
void reallocation_survival(
	struct Grid *grid, 
	struct Pchar *pchar, 
	struct Pmas *mass
){
	double ral_cap_stf, ral_cap_rtf;
	double aaa, bbb, ccc;
	double ral_stf, ral_rtf;
	
	if(pchar->v_type == 1){
		switch(grid->veg_olson){
			case 1: case 2: case 3: case 4: case 5: case 6: case 7: case 8: 
			case 9: case 11: 
				ral_cap_stf = 0.03;
				ral_cap_rtf = 0.065;
				break;
			case 10: case 12: 
				ral_cap_stf = 0.03;
				ral_cap_rtf = 0.03;
				break;
			case 13: case 14: case 15: case 16: 
			case 17: case 18: case 19: case 20: case 21: case 22: case 23: case 24: 
			case 25: case 26: case 27: case 28: case 29: case 30: case 31: case 32: case 33: 
				ral_cap_stf = 0.070;
				ral_cap_rtf = 0.12;
				break;
		}
	}else if(pchar->v_type == 2){
		/* switch(grid->veg_sage){
			case 1: case 2: case 3: case 4: case 5: case 6: case 8: 
				ral_cap_stf = 0.03;
				ral_cap_rtf = 0.065;
				break;
			case 7: 
				ral_cap_stf = 0.03;
				ral_cap_rtf = 0.03;
				break;
			case 9: case 10: case 11: case 12: case 13: case 14: case 15:  */
				ral_cap_stf = 0.070;
				ral_cap_rtf = 0.12;
				/* break; */
		/* } */
	}/* else if(pchar->v_type==3){
		ral_cap_stf = 0.070;
		ral_cap_rtf = 0.12;
	} */
		
	/* to fliage: critical reallocation for survival */
	aaa = 0.075 * 100.0 * 2.0 / 2.2 / pchar->sla;
	if(mass->fol < aaa){
		bbb = mass->stm * ral_cap_stf;
		ccc = mass->rot * ral_cap_rtf;
		
		ral_stf = aaa * pchar->alloc_abg * (bbb/aaa)/(0.5 + (bbb/aaa));
		ral_rtf = aaa * (1.0 - pchar->alloc_abg)*(ccc / aaa)/(0.5 + (ccc/aaa));

		mass->fol += ral_stf+ral_rtf;
		mass->stm -= ral_stf;
		mass->rot -= ral_rtf;
		
		f_leaf_age(0, pchar, mass, ral_stf+ral_rtf);
	}
	
	/* to fliage, improve production */
	aaa = (pchar->opt_lai[grid->m] - mass->lai[grid->m])*100.0*2.0/2.2/pchar->sla;
	if(aaa>0.0){
		bbb = mass->stm * ral_cap_stf;
		ccc = mass->rot * ral_cap_rtf;
		
		ral_stf = aaa*pchar->alloc_abg * (bbb/aaa)/(1.5 + (bbb/aaa));
		ral_rtf = aaa*(1.0 - pchar->alloc_abg)*(ccc/aaa)/(1.5 + (ccc/aaa));

		mass->fol += ral_stf + ral_rtf;
		mass->stm -= ral_stf;
		mass->rot -= ral_rtf;
		
		f_leaf_age(0, pchar, mass, ral_stf + ral_rtf);
	}
	
	mass->lai[grid->m] = lai_mass(grid, mass, pchar);
	
	/* to stem and branch */
	if(mass->stm <= 0.0){
		if(mass->rot > 0.0){
			aaa = mass->rot * pchar->alloc_abg;
			mass->rot -= aaa;
			mass->stm += aaa;
		}
	}
	
	/* to root */
	if(mass->rot <= 0.0){
		if(mass->stm>0.0){
			aaa = mass->stm*(1.0 - pchar->alloc_abg);
			mass->rot += aaa;
			mass->stm -= aaa;
		}
	}
}

