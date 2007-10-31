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

#define TER_CON 0.001 /* criteria for determining the equilibrium, NEP value in Mg C ha-1 yr-1 */

extern double h_tmp[201], h_pre[201], h_dswr[201], h_aet[201], h_rof[201];
extern double h_gpp[201], h_npp[201], h_nep[201], h_plant[201], h_soil[201];
extern double h_sr[201], h_ersn_c[201];
extern double h_agrarea[201], h_luc[201];

extern double h_burnt_area[201];
extern double h_bioburn_co2[201], h_bioburn_ch4[201], h_bioburn_co[201];
extern double h_bioburn_nmhc[201], h_bioburn_oc[201], h_bioburn_bc[201];
extern double h_bioburn_nox[201], h_bioburn_so2[201], h_bioburn_pm25[201];
extern double h_bioburn_tpm[201], h_bioburn_tec[201];

extern double h_ch4ox1[201], h_ch4ox2[201], h_ch4ox3[201];
extern double h_ch4emit_cao_paddy[201], h_ch4emit_cao_wetland[201];
extern double h_n2o_emit_ngas[201], h_n2_emit_ngas[201];
extern double h_n2o_emit_casa[201], h_no_emit_casa[201], h_n2_emit_casa[201];
extern double h_nh3_emit[201];

extern double h_voc_isopr_g97[201], h_voc_monotrp_g97[201], h_voc_methanl_g97[201];
extern double h_voc_acetone_g97[201], h_voc_actaldhd_g97[201], h_voc_frmardhd_g97[201];
extern double h_voc_formacd_g97[201], h_voc_acetacd_g97[201], h_voc_co_g97[201];

extern double m_ch4ox1[12], m_ch4ox2[12], m_ch4ox3[12];
extern double m_bioburn_co2[12], m_bioburn_ch4[12], m_bioburn_co[12];
extern double m_bioburn_nmhc[12], m_bioburn_oc[12], m_bioburn_bc[12];

/**************** EQUILIBRIUM ******************/
void cal_stable(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux, 
	FILE *fp_carbon, 
	FILE *fp_nitrogen, 
	FILE *fp_ersn, 
	FILE *fp_ghg, 
	FILE *fp_bioburn, 
	FILE *fp_voc
){
	long f, nn, term_time;
	double acnep, plantmass;
	double pre_ann, ltr_ann, hrl_ann, msl_ann, hrm_ann;
	double ch4ox1_ann, ch4ox2_ann, ch4ox3_ann, sw1_ann, wfps_ann;
	
	/** maximum simulation times **/
	grid->phase = 0; /* spin-up */
	term_time = 4000;	
	
	(echar->soil).rl = (echar->soil).rl0;
	(echar->soil).rh = (echar->soil).rh0;
	
	f_cult_luc(grid);
	grid->f_crop_p = grid->fcrop[199];
	
	/********** roop to stable stage ***********/
	nn = 0; 
	acnep = 10.0;
	while(acnep>TER_CON){ /*** acnep>TER_CON nn<10 ***/
		grid->y = nn;
		
		plantmass = acnep = 0.0;
		pre_ann = ltr_ann = hrl_ann = msl_ann = hrm_ann = 0.0;
		ch4ox1_ann = ch4ox2_ann = ch4ox3_ann = sw1_ann = wfps_ann = 0.0;
		
		/***** empirical model NPP*****/
		if(grid->y==0){ /* for the first year */
			npp_empirical(grid, loct, flux);
		}
		
		for(f=0;f<12;f++){
			grid->m = f;
			
			/* initialize N fluxes ************/
			n_flux_zero(f, flux);

			cd_trend(grid);
			co2_in_canopy(grid, loct, mass, flux);
			
			/*** environmental condition ***/
			dynmcL(grid, loct, mass, echar);
			
			/***** vegetation processes *****/
			biome_processes(grid, loct, echar, mass, flux);
			plant_stand(grid, loct, mass, flux);
			
			(flux->soil).lL[f] = (flux->plant).lL[f];			
			(flux->soil).d13c_lL[f] = (flux->plant).d13c_lL[f];
			
			/***** soil processes *****/
			soil_processes(grid, loct, &(echar->soil), &(mass->soil), &(flux->soil));
			
			/* CH4 oxydation **************/
			f_ch4oxy_ridgewell(grid, loct, flux);
			f_ch4oxy_casa(grid, loct, flux);
			f_ch4oxy_delgrosso(grid, loct, flux);

			/* CH4 emission (wetlands) */
			f_ch4_emit_cao(grid, loct, flux);
			
			/* N2O emission */
			f_n2o_emit_ngas(grid, loct, mass, flux);
			f_n2o_emit_casa(grid, loct, mass, flux);
			
			/* VOC emission *****************/
			f_voc_emit_guenther97(grid, loct, mass, flux);

			/***** ecosystem mass balance *****/	
			/* net ecosystem production */
			flux->nep[f] = (flux->plant).npp[f] - (flux->soil).rS[f];
			/* total ecosystem carbon storage */
			mass->total[f] = (mass->c3).plant[f]*loct->C3ptn[f]+(mass->c4).plant[f]*loct->C4ptn[f]+(mass->soil).soil[f];
			/** net carbon balance taking crop harvest into account **/
			flux->ncb[f] = flux->nep[f] + (flux->plant).hvst[f];
			
			/* carbon isotope */
			d13c_efflux(grid, loct, flux);

			/* nitrogen budget */
			n_budget(grid, loct, mass, flux);
						
			/* annual average plant mass */
			plantmass += (mass->plant).plant[f]/12.0;
			
			if(grid->veg_olson!=29&&grid->veg_olson!=30&&grid->veg_olson!=31&&grid->veg_olson!=32){
				acnep += flux->nep[f]; /* annula NEP, flux->nep[f]; */
			}else{
				acnep += flux->ncb[f]; /* annula NCB, flux->nep[f]; */
			}
			
			/*** acclimation ***/
			(echar->soil).ft0_l[f] = (echar->soil).ft_l[f];
			(echar->soil).ft0_h[f] = (echar->soil).ft_h[f];
			(echar->soil).fm0_l[f] = (echar->soil).fm_l[f];
			(echar->soil).fm0_h[f] = (echar->soil).fm_h[f];
			flux->lL0[f] = (flux->plant).lL[f];
			
			pre_ann += grid->prate_sfc[f];
			ltr_ann += (mass->soil).ltr_m[f]*(double)(grid->mm[f])/365.0;
			msl_ann += (mass->soil).msl_m[f]*(double)(grid->mm[f])/365.0;
			hrl_ann += (flux->soil).rl[f];
			hrm_ann += (flux->soil).rh[f];
			
			ch4ox1_ann += (flux->soil).ch4oxy_ridg[grid->m];
			ch4ox2_ann += (flux->soil).ch4oxy_casa[grid->m];
			ch4ox3_ann += (flux->soil).ch4oxy_delgrosso[grid->m];
			sw1_ann += loct->sw30*(double)(grid->mm[f])/365.0;
			wfps_ann += loct->wfps[f]*(double)(grid->mm[f])/365.0;
		}
		
		/* biomass burning */
		f_biomassburning(grid, loct, mass, flux);
		
		/* erosion */
		f_erosion(grid, loct, echar, mass, flux);
		
		if(ERSN_CC==1){
			(mass->soil).ltr -= flux->erod_carbon;
			if((mass->soil).ltr<0.0){
				(mass->soil).ltr = 0.0;
			}
		}
		
		/****** terminal conditions ********/
		if(nn<200){	
			/* continued */
			acnep = 10.0; 
		}else if(nn>=200 && nn<term_time){	
			acnep = fabs(acnep); /**** 1. sufficiently stabilized ****/	
		}else if(nn>=term_time){
			break; /**** 3. stop by 2000 years ****/	
		}
		
		if(plantmass<0.0||plantmass>=500.0){
			vanish(mass, flux); /*** 2. excluding abnormal estimates ***/
			break;
		}
		grid->time = nn; /* simulation time of carbon budget */
		nn++;
	}
	/****** end of stsbilization roop *******/
	
	for(f=0;f<12;f++){
		h_tmp[0] += grid->tmp_2m[f]*(double)(grid->mm[f])/365.0 * grid->area;
		h_pre[0] += grid->prate_sfc[f]*0.1 * grid->area;
		h_dswr[0] += grid->gl_rad[f]*(double)(grid->mm[f])/365.0 * grid->area;
		h_aet[0] += (loct->evpr[f] + loct->trspr[f] + loct->incep[f])*0.1 * grid->area;
		h_rof[0] += loct->ro2[f]*0.1 * grid->area;
		h_gpp[0] += (flux->plant).gpp[f] * grid->area;
		h_npp[0] += (flux->plant).npp[f] * grid->area;
		h_nep[0] += flux->nep[f] * grid->area;
		h_plant[0] += ((mass->plant).mfol[f] + (mass->plant).mstm[f] + 
			(mass->plant).mrot[f])*(double)(grid->mm[f])/365.0 * grid->area;
		h_soil[0] += ((mass->soil).ltr_m[f] + 
			(mass->soil).msl_m[f])*(double)(grid->mm[f])/365.0 * grid->area;
		h_sr[0] += ((flux->plant).rrm[f] + (flux->plant).rrg[f] + (flux->soil).rS[f]) * grid->area;
		
		/* GHG */
		h_ch4ox1[0] += (flux->soil).ch4oxy_ridg[f] * grid->area *10000.0/1000.0;
		h_ch4ox2[0] += (flux->soil).ch4oxy_casa[f] * grid->area *10000.0/1000.0;
		h_ch4ox3[0] += (flux->soil).ch4oxy_delgrosso[f] * grid->area *10000.0/1000.0;
		h_ch4emit_cao_paddy[0] += (flux->soil).ch4emit_paddy_cao[f] * grid->area *10000.0/1000.0;
		h_ch4emit_cao_wetland[0] += (flux->soil).ch4emit_wetland_cao[f] * grid->area *10000.0/1000.0;
		h_n2o_emit_ngas[0] += (flux->soil).d_n2o_ngas[f] * grid->area;
		h_n2_emit_ngas[0] += (flux->soil).d_n2_ngas[f] * grid->area;
		h_n2o_emit_casa[0] += (flux->soil).d_n2o_casa[f] * grid->area;
		h_n2_emit_casa[0] += (flux->soil).d_n2_casa[f] * grid->area;
		h_no_emit_casa[0] += (flux->soil).d_no_casa[f] * grid->area;
		h_nh3_emit[0] += (flux->soil).n_nh3vlt[f] * grid->area;
		
		/* biomass burninig */
		h_burnt_area[0] += flux->a_burnt[f]*grid->area;
		h_bioburn_co2[0] += (flux->bb_co2_litter[f]+flux->bb_co2_leaf[f]+flux->bb_co2_wood[f]+flux->bb_co2_root[f]) * grid->area;
		h_bioburn_co[0] += (flux->bb_co_litter[f]+flux->bb_co_leaf[f]+flux->bb_co_wood[f]+flux->bb_co_root[f]) * grid->area;
		h_bioburn_ch4[0] += (flux->bb_ch4_litter[f]+flux->bb_ch4_leaf[f]+flux->bb_ch4_wood[f]+flux->bb_ch4_root[f]) * grid->area;
		h_bioburn_nmhc[0] += (flux->bb_nmhc_litter[f]+flux->bb_nmhc_leaf[f]+flux->bb_nmhc_wood[f]+flux->bb_nmhc_root[f]) * grid->area;
		h_bioburn_oc[0] += (flux->bb_oc_litter[f]+flux->bb_oc_leaf[f]+flux->bb_oc_wood[f]+flux->bb_oc_root[f]) * grid->area;
		h_bioburn_bc[0] += (flux->bb_bc_litter[f]+flux->bb_bc_leaf[f]+flux->bb_bc_wood[f]+flux->bb_bc_root[f]) * grid->area;
		h_bioburn_nox[0] += (flux->bb_bc_litter[f]+flux->bb_bc_leaf[f]+flux->bb_bc_wood[f]+flux->bb_bc_root[f]) * grid->area;
		h_bioburn_so2[0] += (flux->bb_so2_litter[f]+flux->bb_so2_leaf[f]+flux->bb_so2_wood[f]+flux->bb_so2_root[f]) * grid->area;
		h_bioburn_pm25[0] += (flux->bb_pm25_litter[f]+flux->bb_pm25_leaf[f]+flux->bb_pm25_wood[f]+flux->bb_pm25_root[f]) * grid->area;
		h_bioburn_tpm[0] += (flux->bb_tpm_litter[f]+flux->bb_tpm_leaf[f]+flux->bb_tpm_wood[f]+flux->bb_tpm_root[f]) * grid->area;
		h_bioburn_tec[0] += (flux->bb_tec_litter[f]+flux->bb_tec_leaf[f]+flux->bb_tec_wood[f]+flux->bb_tec_root[f]) * grid->area;

		/* VOC, g C month-1 */
		h_voc_isopr_g97[0] += flux->voc_isopr_g97[f] * grid->area *10000.0/1000000.0; 
		h_voc_monotrp_g97[0] += flux->voc_monotrp_g97[f] * grid->area *10000.0/1000000.0; 
		h_voc_methanl_g97[0] += flux->voc_methanl_g97[f] * grid->area *10000.0/1000000.0; 
		h_voc_acetone_g97[0] += flux->voc_acetone_g97[f] * grid->area *10000.0/1000000.0; 
		h_voc_actaldhd_g97[0] += flux->voc_actaldhd_g97[f] * grid->area *10000.0/1000000.0; 
		h_voc_frmardhd_g97[0] += flux->voc_frmardhd_g97[f] * grid->area *10000.0/1000000.0; 
		h_voc_formacd_g97[0] += flux->voc_formacd_g97[f] * grid->area *10000.0/1000000.0; 
		h_voc_acetacd_g97[0] += flux->voc_acetacd_g97[f] * grid->area *10000.0/1000000.0; 
		h_voc_co_g97[0] += flux->voc_co_g97[f] * grid->area *10000.0/1000000.0; 
	}
	
	/* erosion */
	h_ersn_c[0] += flux->erod_carbon * grid->area;

	/* land use change */
	f_luc_emit(grid, mass, flux);

	h_agrarea[0] += grid->area * grid->f_crop;
	h_luc[0] += (flux->lu_conv + flux->lu_ten + flux->lu_hund) * grid->area;
	
	/** output initial stable state **/
	publish_cbud(grid, loct, echar, mass, flux, fp_carbon);
	
	/* erosion */
	fprintf(fp_ersn,"1900 ");
	fprintf(fp_ersn,"%lf ", pre_ann);
	fprintf(fp_ersn,"%lf ", grid->f_crop);
	fprintf(fp_ersn,"%lf ", grid->f_erosion_r);
	fprintf(fp_ersn,"%lf ", grid->f_erosion_ls);
	fprintf(fp_ersn,"%lf ", grid->f_erosion_k);
	fprintf(fp_ersn,"%lf ", grid->f_erosion_c);
	fprintf(fp_ersn,"%lf ", grid->f_erosion_p);

	fprintf(fp_ersn,"%lf ", flux->erod_soil);
	fprintf(fp_ersn,"%lf ", flux->erod_orgmat);
	fprintf(fp_ersn,"%lf ", flux->erod_carbon);
	
	fprintf(fp_ersn,"%lf ", ltr_ann);
	fprintf(fp_ersn,"%lf ", hrl_ann);
	fprintf(fp_ersn,"%lf ", msl_ann);
	fprintf(fp_ersn,"%lf ", hrm_ann);
	fprintf(fp_ersn,"\n");
	
	/* nitrogen */
	fprintf(fp_nitrogen,"1900 ");
	for(f=0;f<12;f++){
		fprintf(fp_nitrogen,"%.2lf ", (mass->plant).n_cnpy_m[f]);
		fprintf(fp_nitrogen,"%.2lf ", (mass->plant).n_strg_m[f]);
		fprintf(fp_nitrogen,"%.2lf ", (mass->soil).n_no3_m[f]);
		fprintf(fp_nitrogen,"%.2lf ", (mass->soil).n_nh4_m[f]);
		fprintf(fp_nitrogen,"%.2lf ", (mass->soil).n_mcrb_m[f]);
		fprintf(fp_nitrogen,"%.2lf ", (mass->soil).n_lttr_m[f]);
		fprintf(fp_nitrogen,"%.2lf ", (mass->soil).n_hums_m[f]);
	
		fprintf(fp_nitrogen,"%.2lf ", (flux->plant).n_biofix[f]);
		fprintf(fp_nitrogen,"%.2lf ", (flux->plant).uptake_no3[f]);
		fprintf(fp_nitrogen,"%.2lf ", (flux->plant).uptake_nh4[f]);

		fprintf(fp_nitrogen,"%.2lf ", (flux->soil).n_leach[f]);
		fprintf(fp_nitrogen,"%.2lf ", (flux->soil).n_nh3vlt[f]);
	}
	fprintf(fp_nitrogen,"\n");
	
	/* GHG & trace gases */
	fprintf(fp_ghg,"1900 ");
	fprintf(fp_ghg,"%.3lf %.3lf ", sw1_ann, wfps_ann);
	for(f=0;f<12;f++){
		fprintf(fp_ghg,"%.3lf ", (flux->soil).ch4oxy_ridg[f]);
		fprintf(fp_ghg,"%.3lf ", (flux->soil).ch4oxy_casa[f]);
		fprintf(fp_ghg,"%.3lf ", (flux->soil).ch4oxy_delgrosso[f]);
		fprintf(fp_ghg,"%.3lf ", (flux->soil).ch4emit_paddy_cao[f]);
		fprintf(fp_ghg,"%.3lf ", (flux->soil).ch4emit_wetland_cao[f]);
		fprintf(fp_ghg,"%.3lf ", (flux->soil).d_n2o_ngas[f]);
		fprintf(fp_ghg,"%.3lf ", (flux->soil).d_n2_ngas[f]);
		fprintf(fp_ghg,"%.3lf ", (flux->soil).d_n2o_casa[f]);
		fprintf(fp_ghg,"%.3lf ", (flux->soil).d_n2_casa[f]);
		fprintf(fp_ghg,"%.3lf ", (flux->soil).d_no_casa[f]);
		fprintf(fp_ghg,"%.3lf ", (flux->soil).n_nh3vlt[f]);
	}
	fprintf(fp_ghg,"\n");

	/* biomass burning */
	fprintf(fp_bioburn,"1900 ");
	for(f=0;f<12;f++){
		fprintf(fp_bioburn,"%.3lf ", loct->msw30[f]);
		fprintf(fp_bioburn,"%.3lf ", flux->day_fire[f]);
		fprintf(fp_bioburn,"%.3lf ", flux->a_burnt[f]);
		fprintf(fp_bioburn,"%.3lf ", flux->bb_co2_litter[f]+flux->bb_co2_leaf[f]+flux->bb_co2_wood[f]+flux->bb_co2_root[f]);
		fprintf(fp_bioburn,"%.3lf ", flux->bb_co_litter[f]+flux->bb_co_leaf[f]+flux->bb_co_wood[f]+flux->bb_co_root[f]);
		fprintf(fp_bioburn,"%.3lf ", flux->bb_ch4_litter[f]+flux->bb_ch4_leaf[f]+flux->bb_ch4_wood[f]+flux->bb_ch4_root[f]);
		fprintf(fp_bioburn,"%.3lf ", flux->bb_nmhc_litter[f]+flux->bb_nmhc_leaf[f]+flux->bb_nmhc_wood[f]+flux->bb_nmhc_root[f]);
		fprintf(fp_bioburn,"%.3lf ", flux->bb_oc_litter[f]+flux->bb_oc_leaf[f]+flux->bb_oc_wood[f]+flux->bb_oc_root[f]);
		fprintf(fp_bioburn,"%.3lf ", flux->bb_bc_litter[f]+flux->bb_bc_leaf[f]+flux->bb_bc_wood[f]+flux->bb_bc_root[f]);
		fprintf(fp_bioburn,"%.3lf ", flux->bb_nox_litter[f]+flux->bb_nox_leaf[f]+flux->bb_nox_wood[f]+flux->bb_nox_root[f]);
		fprintf(fp_bioburn,"%.3lf ", flux->bb_so2_litter[f]+flux->bb_so2_leaf[f]+flux->bb_so2_wood[f]+flux->bb_so2_root[f]);
		fprintf(fp_bioburn,"%.3lf ", flux->bb_pm25_litter[f]+flux->bb_pm25_leaf[f]+flux->bb_pm25_wood[f]+flux->bb_pm25_root[f]);
		fprintf(fp_bioburn,"%.3lf ", flux->bb_tpm_litter[f]+flux->bb_tpm_leaf[f]+flux->bb_tpm_wood[f]+flux->bb_tpm_root[f]);
		fprintf(fp_bioburn,"%.3lf ", flux->bb_tec_litter[f]+flux->bb_tec_leaf[f]+flux->bb_tec_wood[f]+flux->bb_tec_root[f]);
	}
	fprintf(fp_bioburn,"\n");

	/* VOC */
	fprintf(fp_voc,"1900 ");
	for(f=0;f<12;f++){
		fprintf(fp_voc,"%.3lf ", flux->voc_isopr_g97[f]);
		fprintf(fp_voc,"%.3lf ", flux->voc_monotrp_g97[f]);
		fprintf(fp_voc,"%.3lf ", flux->voc_methanl_g97[f]);
		fprintf(fp_voc,"%.3lf ", flux->voc_acetone_g97[f]);
		fprintf(fp_voc,"%.3lf ", flux->voc_actaldhd_g97[f]);
		fprintf(fp_voc,"%.3lf ", flux->voc_frmardhd_g97[f]);
		fprintf(fp_voc,"%.3lf ", flux->voc_formacd_g97[f]);
		fprintf(fp_voc,"%.3lf ", flux->voc_acetacd_g97[f]);
		fprintf(fp_voc,"%.3lf ", flux->voc_co_g97[f]);
	}
	fprintf(fp_voc,"\n");
}

