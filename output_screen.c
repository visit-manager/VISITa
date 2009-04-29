/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* make a result-output file and show values in the console port */
#include<stdio.h>
#include<stdlib.h>
#include"structure.h"
#include"prototype.h"

/******** show result on the console window ********/
void screenshow(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Flux *flux, 
	struct Echar *echar
){
	short h;
	double ann;
	double fol, stm, rot, ltr, msl;
	double gpp, ar, npp, hr, nep;
	
	printf("No:%ld A:%7.2lf B-olson:%ld B-sage:%ld HT:%ld CT:%ld WHC30:%.1lf WHC:%.1lf BD:%.3lf UL:%.3lf\n",
			grid->n_olson, grid->area, grid->veg_olson, grid->veg_sage, loct->time_hyd, loct->time, grid->field_cap1, 
			grid->field_cap2, grid->bulkdens, grid->f_upland);
	
	printf("CO2: ");
	for(h=0;h<12;h++){		printf("%6.0lf ",grid->bCO2[h]);		}		printf("\n"); 
	printf("T2M: ");
	for(h=0;h<12;h++){		printf("%6.0lf ",grid->tmp_2m[h]);		}		printf("\n"); 
	
	printf("PRE: ");	ann=0.0;
	for(h=0;h<12;h++){		printf("%6.0lf ",grid->prate_sfc[h]);		ann+=grid->prate_sfc[h];	}		printf(" %.1lf\n",ann); 
	printf("PAR: ");
	for(h=0;h<12;h++){		printf("%6.0lf ",grid->par[h]);		}			printf("\n"); 
	
	printf("SWU: ");
	for(h=0;h<12;h++){		printf("%6.0lf ",loct->msw30[h]);	}			printf("\n"); 
	printf("SWW: ");
	for(h=0;h<12;h++){		printf("%6.0lf ",loct->msww[h]);	}			printf("\n"); 
	
	printf("LA3: ");
	for(h=0;h<12;h++){		printf("%6.1lf ",(mass->c3).lai[h]);		}	printf("\n"); 
	printf("GD3: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%6.1lf ",(flux->c3).gpp_df97[h]);		ann+=(flux->c3).gpp_df97[h];	}		printf(" %.1lf\n",ann); 
	printf("NP3: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%6.1lf ",(flux->c3).npp[h]);		ann+=(flux->c3).npp[h];	}		printf(" %.1lf\n",ann); 
	printf("LA4: ");
	for(h=0;h<12;h++){		printf("%6.1lf ",(mass->c4).lai[h]);		}	printf("\n"); 
	printf("GD4: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%6.1lf ",(flux->c4).gpp_df97[h]);		ann+=(flux->c4).gpp_df97[h];	}		printf(" %.1lf\n",ann); 
	printf("NP4: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%6.1lf ",(flux->c4).npp[h]);		ann+=(flux->c4).npp[h];	}		printf(" %.1lf\n",ann); 
	printf("NEP: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%6.1lf ",flux->nep[h]);		ann+=flux->nep[h];	}		printf(" %.1lf\n",ann); 
	

/*	printf("GP3: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%6.1lf ",(flux->c3).gpp[h]);		ann+=(flux->c3).gpp[h];	}		printf(" %.1lf\n",ann); 
	printf("GP4: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%6.1lf ",(flux->c4).gpp[h]);		ann+=(flux->c4).gpp[h];	}		printf(" %.1lf\n",ann); 
	printf("XX1: ");
	for(h=0;h<12;h++){		printf("%6.3lf ",loct->xx1[h]);		}			printf("\n"); 
	printf("XX2: ");
	for(h=0;h<12;h++){		printf("%6.3lf ",loct->xx2[h]);		}			printf("\n"); 
	printf("XX3: ");
	for(h=0;h<12;h++){		printf("%6.3lf ",loct->xx3[h]);		}			printf("\n"); 
	printf("XX4: ");
	for(h=0;h<12;h++){		printf("%6.3lf ",loct->xx4[h]);		}			printf("\n"); 
	printf("XX5: ");
	for(h=0;h<12;h++){		printf("%6.3lf ",loct->xx5[h]);		}			printf("\n");

	printf("DS3: ");
	for(h=0;h<12;h++){		printf("%6.3lf ",(echar->c3).fsw[h]);		}	printf("\n"); 
	printf("DS4: ");
	for(h=0;h<12;h++){		printf("%6.3lf ",(echar->c4).fsw[h]);		}	printf("\n"); */

/*	printf("MO1: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.1lf ",(flux->soil).ch4oxy_ridg[h]);		ann+=(flux->soil).ch4oxy_ridg[h];	}		printf(" %.1lf\n",ann); 
	printf("MO2: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.1lf ",(flux->soil).ch4oxy_casa[h]);		ann+=(flux->soil).ch4oxy_casa[h];	}		printf(" %.1lf\n",ann); 
	printf("MO3: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.1lf ",(flux->soil).ch4oxy_delgrosso[h]);		ann+=(flux->soil).ch4oxy_delgrosso[h];	}		printf(" %.1lf\n",ann); 
	printf("MO4: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.1lf ",(flux->soil).ch4oxy_curry[h]);		ann+=(flux->soil).ch4oxy_curry[h];	}		printf(" %.1lf\n",ann); */

/*	printf("GS3: ");
	for(h=0;h<12;h++){		printf("%5.1lf ",(echar->c3).gs[h]);		}	printf("\n"); 
	printf("GS4: ");
	for(h=0;h<12;h++){		printf("%5.1lf ",(echar->c4).gs[h]);		}	printf("\n"); 
	printf("PS3: ");
	for(h=0;h<12;h++){		printf("%5.1lf ",(echar->c3).psat[h]);		}	printf("\n"); 
	printf("PS4: ");
	for(h=0;h<12;h++){		printf("%5.1lf ",(echar->c4).psat[h]);		}	printf("\n");  */
	
/*	printf("FL3: ");
	for(h=0;h<12;h++){		printf("%5.2lf ",(mass->c3).mfol[h]);		}	printf("\n");	*/
/*	printf("FVG: ");
	for(h=0;h<12;h++){		printf("%5.3lf ",loct->f_vegcov[h]);		}		printf("\n"); */

/*	printf("DOC: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.2lf ",(mass->soil).doc_m[h]);		ann+=(mass->soil).doc_m[h];	}		printf(" %.2lf\n",ann); 
	printf("DOC: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.2lf ",(flux->soil).doc_boyer[h]);		ann+=(flux->soil).doc_boyer[h];	}		printf(" %.2lf\n",ann); */

/*	printf("MTP: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.3lf ",(flux->soil).ch4flux_paddy_cao[h]);		ann+=(flux->soil).ch4flux_paddy_cao[h];	}		printf(" %.1lf\n",ann); 
	printf("MTW: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.3lf ",(flux->soil).ch4flux_wetland_cao[h]);		ann+=(flux->soil).ch4flux_wetland_cao[h];	}		printf(" %.1lf\n",ann); */

/*	printf("LAD\n");
	for(h=0;h<12;h++){		printf("%8.6lf ",(echar->c3).fleaf_age[h]);		}	printf("\n"); 
	for(h=12;h<24;h++){		printf("%8.6lf ",(echar->c3).fleaf_age[h]);		}	printf("\n"); 
	for(h=24;h<36;h++){		printf("%8.6lf ",(echar->c3).fleaf_age[h]);		}	printf("\n"); 
	for(h=36;h<=48;h++){		printf("%8.6lf ",(echar->c3).fleaf_age[h]);		}	printf("\n"); */

/*	printf("RSL: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.2lf ", (flux->soil).rS[h]);		ann+=(flux->soil).rS[h];	}		printf(" %.1lf\n",ann); */

/*	printf("WFP: ");
	for(h=0;h<12;h++){		printf("%5.2lf ",loct->wfps[h]);		}			printf("\n"); 
	printf("M_M: ");
	for(h=0;h<12;h++){		printf("%5.2lf ",loct->m_m[h]);		}			printf("\n"); 
	printf("M_E: ");
	for(h=0;h<12;h++){		printf("%5.2lf ",loct->m_e[h]);		}			printf("\n"); 
	printf("I_W: ");
	for(h=0;h<12;h++){		printf("%5.2lf ",loct->i_w[h]);		}			printf("\n"); */
	
/*	 */

/*	printf("NCN: ");
	for(h=0;h<12;h++){		printf("%5.0lf ",(mass->plant).n_cnpy_m[h]/1000.0);		}			printf("\n"); 
	printf("NST: ");
	for(h=0;h<12;h++){		printf("%5.0lf ",(mass->plant).n_strg_m[h]/1000.0);		}			printf("\n"); 
	printf("NH4: ");
	for(h=0;h<12;h++){		printf("%5.0lf ",(mass->soil).n_nh4_m[h]/1000.0);		}			printf("\n"); 
	printf("NO3: ");
	for(h=0;h<12;h++){		printf("%5.0lf ",(mass->soil).n_no3_m[h]/1000.0);		}			printf("\n"); 
	printf("MCB: ");
	for(h=0;h<12;h++){		printf("%5.0lf ",(mass->soil).n_mcrb_m[h]/1000.0);		}			printf("\n"); 
	printf("SLN: ");
	for(h=0;h<12;h++){		printf("%5.0lf ",(mass->soil).n_lttr_m[h]/1000.0);		}			printf("\n"); 
	printf("SHN: ");
	for(h=0;h<12;h++){		printf("%5.0lf ",(mass->soil).n_hums_m[h]/1000.0);		}			printf("\n"); 
	printf("LN3: ");
	for(h=0;h<12;h++){		printf("%5.1lf ",(echar->c3).n_conc_larea_m[h]/1000.0);		}			printf("\n"); 
	printf("LN4: ");
	for(h=0;h<12;h++){		printf("%5.1lf ",(echar->c4).n_conc_larea_m[h]/1000.0);		}			printf("\n"); */

/*	printf("DN3: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",loct->depo_no3[h]);		ann+=loct->depo_no3[h];	}		printf(" %7.1lf\n",ann); 
	printf("DN4: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",loct->depo_nh4[h]);		ann+=loct->depo_nh4[h];	}		printf(" %7.1lf\n",ann); */

/*	printf("BFX: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->plant).n_biofix[h]);		ann+=(flux->plant).n_biofix[h];	}		printf(" %7.1lf\n",ann); 
	printf("UN3: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->plant).uptake_no3[h]);		ann+=(flux->plant).uptake_no3[h];	}		printf(" %7.1lf\n",ann); 
	printf("UN4: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->plant).uptake_nh4[h]);		ann+=(flux->plant).uptake_nh4[h];	}		printf(" %7.1lf\n",ann); 
	printf("ALC: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->plant).n_alloc_cnpy[h]);		ann+=(flux->plant).n_alloc_cnpy[h];	}		printf(" %7.1lf\n",ann); 
	printf("ALS: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->plant).n_alloc_strg[h]);		ann+=(flux->plant).n_alloc_strg[h];	}		printf(" %7.1lf\n",ann); 
	printf("SAL: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->plant).n_salvage[h]);		ann+=(flux->plant).n_salvage[h];	}		printf(" %7.1lf\n",ann); 
	printf("RAL: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->plant).n_realloc[h]);		ann+=(flux->plant).n_realloc[h];	}		printf(" %7.1lf\n",ann); 
	printf("ABC: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->plant).n_abdn_cnpy[h]);		ann+=(flux->plant).n_abdn_cnpy[h];	}		printf(" %7.1lf\n",ann); 
	printf("ABS: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->plant).n_abdn_strg[h]);		ann+=(flux->plant).n_abdn_strg[h];	}		printf(" %7.1lf\n",ann); 
	printf("N4V: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->soil).n_nh3vlt[h]);		ann+=(flux->soil).n_nh3vlt[h];	}		printf(" %7.1lf\n",ann); 
	printf("N3L: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->soil).n_leach[h]);		ann+=(flux->soil).n_leach[h];	}		printf(" %7.1lf\n",ann); 
	printf("MLT: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->soil).n_minerlz_lttr[h]);		ann+=(flux->soil).n_minerlz_lttr[h];	}		printf(" %7.1lf\n",ann); 
	printf("MHM: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->soil).n_minerlz_hums[h]);		ann+=(flux->soil).n_minerlz_hums[h];	}		printf(" %7.1lf\n",ann); 
	printf("NTR: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->soil).n_nitrif[h]);		ann+=(flux->soil).n_nitrif[h];	}		printf(" %7.1lf\n",ann); 
	printf("MIM: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->soil).n_immbl[h]);		ann+=(flux->soil).n_immbl[h];	}		printf(" %7.1lf\n",ann); 
	printf("MAB: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->soil).n_mcrb_abdn[h]);		ann+=(flux->soil).n_mcrb_abdn[h];	}		printf(" %7.1lf\n",ann); */
	
/*	printf("N2ON: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->soil).d_n2o_ntr_ngas[h]);		ann+=(flux->soil).d_n2o_ntr_ngas[h];	}		printf(" %7.1lf\n",ann); 
	printf("N2OD: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->soil).d_n2o_dnt_ngas[h]);		ann+=(flux->soil).d_n2o_dnt_ngas[h];	}		printf(" %7.1lf\n",ann); 
	printf("N2N: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->soil).d_n2_ngas[h]);		ann+=(flux->soil).d_n2_ngas[h];	}		printf(" %7.1lf\n",ann); 
	printf("NOC: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->soil).d_no_casa[h]);		ann+=(flux->soil).d_no_casa[h];	}		printf(" %7.1lf\n",ann); 
	printf("N2C: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->soil).d_n2_casa[h]);		ann+=(flux->soil).d_n2_casa[h];	}		printf(" %7.1lf\n",ann); 
	printf("N2OC: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.0lf ",(flux->soil).d_n2o_casa[h]);		ann+=(flux->soil).d_n2o_casa[h];	}		printf(" %7.1lf\n",ann); */

/*	printf("MOR: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.1lf ",(flux->soil).ch4oxy_ridg[h]);		ann += (flux->soil).ch4oxy_ridg[h];		}		printf(" %.4lf\n",ann); 
	printf("MOC: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.1lf ",(flux->soil).ch4oxy_casa[h]);		ann += (flux->soil).ch4oxy_casa[h];		}		printf(" %.4lf\n",ann); 
	printf("MOD: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.1lf ",(flux->soil).ch4oxy_delgrosso[h]);		ann += (flux->soil).ch4oxy_delgrosso[h];	}	printf(" %.4lf\n",ann); */

/*	printf("DFR: ");		
	for(h=0;h<12;h++){		printf("%5.1lf ",flux->day_fire[h]);		}	printf("\n"); 
	printf("ABN: ");		
	for(h=0;h<12;h++){		printf("%5.1lf ",flux->a_burnt[h]*100.0);		}	printf("\n"); 
	printf("BCD: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.1lf ",flux->bb_co2_litter[h]+flux->bb_co2_leaf[h]+flux->bb_co2_wood[h]+flux->bb_co2_root[h]);		
							ann += flux->bb_co2_litter[h]+flux->bb_co2_leaf[h]+flux->bb_co2_wood[h]+flux->bb_co2_root[h];		}		
	printf(" %.2lf\n",ann); 
	printf("BCM: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.2lf ",flux->bb_co_litter[h]+flux->bb_co_leaf[h]+flux->bb_co_wood[h]+flux->bb_co_root[h]);		
							ann += flux->bb_co_litter[h]+flux->bb_co_leaf[h]+flux->bb_co_wood[h]+flux->bb_co_root[h];		}		
	printf(" %.2lf\n",ann); 
	printf("BMT: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.2lf ",flux->bb_ch4_litter[h]+flux->bb_ch4_leaf[h]+flux->bb_ch4_wood[h]+flux->bb_ch4_root[h]);		
							ann += flux->bb_ch4_litter[h]+flux->bb_ch4_leaf[h]+flux->bb_ch4_wood[h]+flux->bb_ch4_root[h];		}		
	printf(" %.2lf\n",ann); 
	printf("BNH: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.2lf ",flux->bb_nmhc_litter[h]+flux->bb_nmhc_leaf[h]+flux->bb_nmhc_wood[h]+flux->bb_nmhc_root[h]);		
							ann += flux->bb_nmhc_litter[h]+flux->bb_nmhc_leaf[h]+flux->bb_nmhc_wood[h]+flux->bb_nmhc_root[h];		}		
	printf(" %.2lf\n",ann); 
	printf("BOC: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.2lf ",flux->bb_oc_litter[h]+flux->bb_oc_leaf[h]+flux->bb_oc_wood[h]+flux->bb_oc_root[h]);		
							ann += flux->bb_oc_litter[h]+flux->bb_oc_leaf[h]+flux->bb_oc_wood[h]+flux->bb_oc_root[h];		}		
	printf(" %.2lf\n",ann); 
	printf("BBC: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.2lf ",flux->bb_bc_litter[h]+flux->bb_bc_leaf[h]+flux->bb_bc_wood[h]+flux->bb_bc_root[h]);		
							ann += flux->bb_bc_litter[h]+flux->bb_bc_leaf[h]+flux->bb_bc_wood[h]+flux->bb_bc_root[h];		}		
	printf(" %.2lf\n",ann); 
*/
/*	printf("MTE: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.1lf ", ((flux->soil).ch4emit_wetland_cao[h] +(flux->soil).ch4emit_paddy_cao[h]));		
							ann += ((flux->soil).ch4emit_wetland_cao[h] +(flux->soil).ch4emit_paddy_cao[h]);		}		
	printf(" %.1lf\n",ann); 

	printf("ISP: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.1lf ", flux->voc_isopr_g97[h]/1000000.0);		
							ann += flux->voc_isopr_g97[h]/1000000.0;		}		
	printf(" %.1lf\n",ann); 

	printf("NON: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.1lf ", (flux->soil).d_n2o_ngas[h]);		
							ann += (flux->soil).d_n2o_ngas[h];		}		
	printf(" %.1lf\n",ann); 

	printf("NOC: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.1lf ", (flux->soil).d_n2o_casa[h]);		
							ann += (flux->soil).d_n2o_casa[h];		}		
	printf(" %.1lf\n",ann); 

	printf("NDP: ");		ann=0.0;
	for(h=0;h<12;h++){		printf("%5.1lf ", loct->depo_no3[h]+loct->depo_nh4[h]);		
							ann += loct->depo_no3[h]+loct->depo_nh4[h];		}		
	printf(" %.1lf\n",ann); */

	fol = stm = rot = ltr = msl = 0.0;
	gpp = ar = npp = hr = nep = 0.0;
	for(h=0;h<12;h++){
		fol += (mass->plant).mfol[h]/12.0;
		stm += (mass->plant).mstm[h]/12.0;
		rot += (mass->plant).mrot[h]/12.0;
		ltr += (mass->soil).ltr_m[h]/12.0;
		msl += (mass->soil).msl_m[h]/12.0;
		
		gpp += (flux->plant).gpp[h];
		ar += (flux->plant).rp[h];
		npp += (flux->plant).npp[h];
		hr += (flux->soil).rS[h];
		nep += flux->nep[h];
	}
	printf("GPP: %.2lf AR: %.2lf NPP: %.2lf HR: %.2lf NEP: %.4lf\n", gpp,ar,npp,hr,nep); 
	printf("F: %.1lf C: %.1lf R: %.1lf L: %.1lf H: %.1lf\n",fol, stm,rot,ltr,msl); 
	printf("Erosion soil: %10.6lf org.mat:%10.6lf carbon:%10.6lf\n", flux->erod_soil, flux->erod_orgmat, flux->erod_carbon); 
	
	printf("Miami		%10.3lf\n", flux->npp_miami);
	printf("Montreal	%10.3lf\n", flux->npp_montreal);
	printf("Schuur		%10.3lf\n", flux->npp_schuur);
	printf("NCEAS		%10.3lf\n", flux->npp_nceas);
}
