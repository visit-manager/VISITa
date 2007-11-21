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
extern double h_ch4_emit_mass[201], h_ch4_emit_photo[201];

extern double h_voc_isopr_g97[201], h_voc_monotrp_g97[201], h_voc_methanl_g97[201];
extern double h_voc_acetone_g97[201], h_voc_actaldhd_g97[201], h_voc_frmardhd_g97[201];
extern double h_voc_formacd_g97[201], h_voc_acetacd_g97[201], h_voc_co_g97[201];

extern double m_ch4ox1[12], m_ch4ox2[12], m_ch4ox3[12];
extern double m_bioburn_co2[12], m_bioburn_ch4[12], m_bioburn_co[12];
extern double m_bioburn_nmhc[12], m_bioburn_oc[12], m_bioburn_bc[12];

extern double v_gpp[34], v_npp[34], v_nep[34];
extern double v_lai[34], v_fol[34], v_stm[34], v_rot[34], v_ltr[34], v_msl[34];

/* 1901-2000 ***********************************************************************/
void cal_cruclim(
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
	long f, g;
	double gpp_ann, ar_ann, npp_ann, hr_ann, nep_ann, ncb_ann, lf_ann, hrl_ann, hrm_ann;
	double fol_ann, stm_ann, rot_ann, ltr_ann, msl_ann, lai_max;
	double tmp_ann, pre_ann, dswr_ann, ipar_ann, apar_ann, aet_ann, rof_ann, rn_ann;
	double iparb_ann, ipard_ann, aparb_ann, apard_ann, ls_gpp, aparmono_ann;
	double ch4ox1_ann, ch4ox2_ann, ch4ox3_ann, sw1_ann, wfps_ann;
	double aaa, bbb;
	
	grid->phase = 1; /* history */
	
	(echar->soil).rl = (echar->soil).rl0;
	(echar->soil).rh = (echar->soil).rh0;
	
	/********** roop to dynamic stage****************************/
	for(g=0;g<CRU_PD;g++){
		/* AD1901 - 2000 */
		
		gpp_ann=ar_ann=npp_ann=hr_ann=nep_ann=ncb_ann = lf_ann = hrl_ann = hrm_ann = 0.0;
		fol_ann=stm_ann=rot_ann=ltr_ann=msl_ann=lai_max = 0.0;
		tmp_ann = pre_ann = dswr_ann = apar_ann = ipar_ann = aet_ann = rof_ann = rn_ann = 0.0;
		iparb_ann = ipard_ann = aparb_ann = apard_ann = ls_gpp = aparmono_ann = 0.0;
		ch4ox1_ann = ch4ox2_ann = ch4ox3_ann = sw1_ann = wfps_ann = 0.0;

		/** climate change ********************/
		grid->climy = PIVOT_CLIMY + g;
		if(grid->cru_exist==1){
			set_cru_clim(grid);
		}
		
		f_cult_luc(grid);
		
		/* CO2 change ********************/
		grid->CO2y = PIVOT_CO2Y + g; 

		for(f=0;f<12;f++){
			grid->m = f;
			
			/* initialize N fluxes ************/
			n_flux_zero(f, flux);

			cd_trend(grid);
			co2_in_canopy(grid, loct, mass, flux);
						
			/******************* environmental condition *******************/
			dynmcL(grid, loct, mass, echar);
			
			/***** vegetation processes *****/
			biome_processes(grid, loct, echar, mass, flux);
			plant_stand(grid, loct,mass,flux);
			
			if(BACC==3){
				(flux->plant).lL[f] = flux->lL0[f];
			}

			(flux->soil).lL[f] = (flux->plant).lL[f];			
			(flux->soil).d13c_lL[f] = (flux->plant).d13c_lL[f];
			
			/***** soil processes *****/
			soil_processes(grid, loct, &(echar->soil), &(mass->soil), &(flux->soil));

			/* CH4 oxydation (uplands) **************/
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
			/* Plant CH4 emission *****************/
			f_ch4_emit_veg(grid, loct, echar, mass, flux);

			/***** ecosystem mass balance *****/	
			/* net ecosystem production */
			flux->nep[f] = (flux->plant).npp[f]-(flux->soil).rS[f];
			/* total ecosystem carbon storage */
			mass->total[f] = (mass->c3).plant[f]*loct->C3ptn[f]+(mass->c4).plant[f]*loct->C4ptn[f]+(mass->soil).soil[f];
			/** net carbon balance taking crop harvest into account **/
			flux->ncb[f] = flux->nep[f]+(flux->plant).hvst[f];
			
			/* carbon isotope */
			d13c_efflux(grid, loct, flux);

			/* nitrogen budget */
			n_budget(grid, loct, mass, flux);
			
			tmp_ann += grid->tmp_2m[f]*(double)(grid->mm[f])/365.0;
			pre_ann += grid->prate_sfc[f];
			dswr_ann += grid->gl_rad[f]*(double)(grid->mm[f])/365.0;

			aaa = loct->C3ptn[f]*(echar->c3).albedo+loct->C4ptn[f]*(echar->c4).albedo;
			bbb = exp(-1.0*(loct->C3ptn[f]*(echar->c3).eK[f]*(mass->c3).lai[f] 
					+ loct->C4ptn[f]*(echar->c4).eK[f]*(mass->c4).lai[f])); /* *(1.0-0.12) */
			apar_ann += (1.0-aaa)*(1.0-bbb)*grid->par[f]*3600.0*grid->dlen[f]*0.5/1000000.0;
			ipar_ann += grid->par[f]*3600.0*grid->dlen[f]*0.5/1000000.0;
			
			iparb_ann += grid->par_bp[f]*3600.0*grid->dlen[f]*0.5/1000000.0;
			ipard_ann += grid->par_dp[f]*3600.0*grid->dlen[f]*0.5/1000000.0;
			aparb_ann += loct->apar_bp[f]*3600.0*grid->dlen[f]*0.5/1000000.0;
			apard_ann += loct->apar_dp[f]*3600.0*grid->dlen[f]*0.5/1000000.0;
			
			aparmono_ann += grid->par[f]*loct->fapar_mono[f]*3600.0*grid->dlen[f]*0.5/1000000.0;
			
			gpp_ann += (flux->plant).gpp[f];
			ar_ann += (flux->plant).rp[f];
			npp_ann += (flux->plant).npp[f];
			hr_ann += (flux->soil).rS[f];
			nep_ann += flux->nep[f];
			ncb_ann += flux->ncb[f];
			lf_ann += (flux->plant).lL[f];
			hrl_ann += (flux->soil).rl[f];
			hrm_ann += (flux->soil).rh[f];
			aet_ann += loct->evpr[f] + loct->trspr[f] + loct->incep[f];
			rn_ann += (loct->rad_net_p[f]+loct->rad_net_g[f])*(double)(grid->mm[f])/365.0;
			rof_ann += loct->ro2[grid->m];
			
			aaa = loct->C3ptn[f]*(mass->c3).lai[f]+loct->C4ptn[f]*(mass->c4).lai[f];
			lai_max = (aaa>lai_max)?aaa:lai_max;
			
			ls_gpp += lTs*grid->dlen[f]*(loct->C3ptn[f]*(echar->c3).psat[f]*(mass->c3).lai[f]
						+loct->C4ptn[f]*(echar->c4).psat[f]*(mass->c4).lai[f]);
			
			fol_ann += (mass->plant).mfol[f]*(double)(grid->mm[f])/365.0;
			stm_ann += (mass->plant).mstm[f]*(double)(grid->mm[f])/365.0;
			rot_ann += (mass->plant).mrot[f]*(double)(grid->mm[f])/365.0;
			ltr_ann += (mass->soil).ltr_m[f]*(double)(grid->mm[f])/365.0;
			msl_ann += (mass->soil).msl_m[f]*(double)(grid->mm[f])/365.0;
			
			ch4ox1_ann += (flux->soil).ch4oxy_ridg[grid->m];
			ch4ox2_ann += (flux->soil).ch4oxy_casa[grid->m];
			ch4ox3_ann += (flux->soil).ch4oxy_delgrosso[grid->m];

			h_tmp[grid->climy - PIVOT_CLIMY +1] += grid->tmp_2m[f]*(double)(grid->mm[f])/365.0 * grid->area;
			h_pre[grid->climy - PIVOT_CLIMY +1] += grid->prate_sfc[f]*0.1 * grid->area;
			h_dswr[grid->climy - PIVOT_CLIMY +1] += grid->gl_rad[f]*(double)(grid->mm[f])/365.0 * grid->area;
			h_aet[grid->climy - PIVOT_CLIMY +1] += (loct->evpr[f] + loct->trspr[f] + loct->incep[f])*0.1 * grid->area;
			h_rof[grid->climy - PIVOT_CLIMY +1] += loct->ro2[f]*0.1 * grid->area;
			h_gpp[grid->climy - PIVOT_CLIMY +1] += (flux->plant).gpp[f] * grid->area;
			h_npp[grid->climy - PIVOT_CLIMY +1] += (flux->plant).npp[f] * grid->area;
			h_nep[grid->climy - PIVOT_CLIMY +1] += flux->nep[f] * grid->area;
			h_plant[grid->climy - PIVOT_CLIMY +1] += ((mass->plant).mfol[f] + (mass->plant).mstm[f] + 
				(mass->plant).mrot[f])*(double)(grid->mm[f])/365.0 * grid->area;
			h_soil[grid->climy - PIVOT_CLIMY +1] += ((mass->soil).ltr_m[f] + 
				(mass->soil).msl_m[f])*(double)(grid->mm[f])/365.0 * grid->area;
			h_sr[grid->climy - PIVOT_CLIMY +1] += ((flux->plant).rrm[f] + (flux->plant).rrg[f] + (flux->soil).rS[f]) * grid->area;
			
			/* GHG */
			h_ch4ox1[grid->climy - PIVOT_CLIMY +1] += (flux->soil).ch4oxy_ridg[f] * grid->area *10000.0/1000.0;
			h_ch4ox2[grid->climy - PIVOT_CLIMY +1] += (flux->soil).ch4oxy_casa[f] * grid->area *10000.0/1000.0;
			h_ch4ox3[grid->climy - PIVOT_CLIMY +1] += (flux->soil).ch4oxy_delgrosso[f] * grid->area *10000.0/1000.0;			
			h_ch4emit_cao_paddy[grid->climy - PIVOT_CLIMY +1] += (flux->soil).ch4emit_paddy_cao[f] * grid->area *10000.0/1000.0;
			h_ch4emit_cao_wetland[grid->climy - PIVOT_CLIMY +1] += (flux->soil).ch4emit_wetland_cao[f] * grid->area *10000.0/1000.0;
			h_n2o_emit_ngas[grid->climy - PIVOT_CLIMY +1] += (flux->soil).d_n2o_ngas[f] * grid->area;
			h_n2_emit_ngas[grid->climy - PIVOT_CLIMY +1] += (flux->soil).d_n2_ngas[f] * grid->area;
			h_n2o_emit_casa[grid->climy - PIVOT_CLIMY +1] += (flux->soil).d_n2o_casa[f] * grid->area;
			h_n2_emit_casa[grid->climy - PIVOT_CLIMY +1] += (flux->soil).d_n2_casa[f] * grid->area;
			h_no_emit_casa[grid->climy - PIVOT_CLIMY +1] += (flux->soil).d_no_casa[f] * grid->area;
			h_nh3_emit[grid->climy - PIVOT_CLIMY +1] += (flux->soil).n_nh3vlt[f] * grid->area;
			h_ch4_emit_mass[grid->climy - PIVOT_CLIMY +1] += (flux->plant).emit_ch4_kirschbaum_mass[f]*10000.0*grid->area;
			h_ch4_emit_photo[grid->climy - PIVOT_CLIMY +1] += (flux->plant).emit_ch4_kirschbaum_photo[f]*10000.0*grid->area;

			/* VOC, g C month-1 */
			h_voc_isopr_g97[grid->climy - PIVOT_CLIMY +1] += flux->voc_isopr_g97[f] * grid->area *10000.0/1000000.0; 
			h_voc_monotrp_g97[grid->climy - PIVOT_CLIMY +1] += flux->voc_monotrp_g97[f] * grid->area *10000.0/1000000.0; 
			h_voc_methanl_g97[grid->climy - PIVOT_CLIMY +1] += flux->voc_methanl_g97[f] * grid->area *10000.0/1000000.0; 
			h_voc_acetone_g97[grid->climy - PIVOT_CLIMY +1] += flux->voc_acetone_g97[f] * grid->area *10000.0/1000000.0; 
			h_voc_actaldhd_g97[grid->climy - PIVOT_CLIMY +1] += flux->voc_actaldhd_g97[f] * grid->area *10000.0/1000000.0; 
			h_voc_frmardhd_g97[grid->climy - PIVOT_CLIMY +1] += flux->voc_frmardhd_g97[f] * grid->area *10000.0/1000000.0; 
			h_voc_formacd_g97[grid->climy - PIVOT_CLIMY +1] += flux->voc_formacd_g97[f] * grid->area *10000.0/1000000.0; 
			h_voc_acetacd_g97[grid->climy - PIVOT_CLIMY +1] += flux->voc_acetacd_g97[f] * grid->area *10000.0/1000000.0; 
			h_voc_co_g97[grid->climy - PIVOT_CLIMY +1] += flux->voc_co_g97[f] * grid->area *10000.0/1000000.0; 			
			
			sw1_ann += loct->sw30*(double)(grid->mm[f])/365.0;
			wfps_ann += loct->wfps[f]*(double)(grid->mm[f])/365.0;

			if(g>=90&&g<=99){
				m_ch4ox1[f] += (flux->soil).ch4oxy_ridg[f] * grid->area *10000.0/1000.0 / 10.0;
				m_ch4ox2[f] += (flux->soil).ch4oxy_casa[f] * grid->area *10000.0/1000.0 / 10.0;
				m_ch4ox3[f] += (flux->soil).ch4oxy_delgrosso[f] * grid->area *10000.0/1000.0 / 10.0;
				m_bioburn_co2[f] += (flux->bb_co2_litter[f]+flux->bb_co2_leaf[f]+flux->bb_co2_wood[f]+flux->bb_co2_root[f]) * grid->area / 10.0;
				m_bioburn_co[f] += (flux->bb_co_litter[f]+flux->bb_co_leaf[f]+flux->bb_co_wood[f]+flux->bb_co_root[f]) * grid->area / 10.0;
				m_bioburn_ch4[f] += (flux->bb_ch4_litter[f]+flux->bb_ch4_leaf[f]+flux->bb_ch4_wood[f]+flux->bb_ch4_root[f]) * grid->area / 10.0;
				m_bioburn_nmhc[f] += (flux->bb_nmhc_litter[f]+flux->bb_nmhc_leaf[f]+flux->bb_nmhc_wood[f]+flux->bb_nmhc_root[f]) * grid->area / 10.0;
				m_bioburn_oc[f] += (flux->bb_oc_litter[f]+flux->bb_oc_leaf[f]+flux->bb_oc_wood[f]+flux->bb_oc_root[f]) * grid->area / 10.0;
				m_bioburn_bc[f] += (flux->bb_bc_litter[f]+flux->bb_bc_leaf[f]+flux->bb_bc_wood[f]+flux->bb_bc_root[f]) * grid->area / 10.0;
				
				v_gpp[grid->veg_olson] += (flux->plant).gpp[f]/10.0 * grid->area;
				v_npp[grid->veg_olson] += (flux->plant).npp[f]/10.0 * grid->area;
				v_nep[grid->veg_olson] += flux->nep[f]/10.0 * grid->area;
				v_lai[grid->veg_olson] += (mass->plant).lai[f]*(double)(grid->mm[f])/365.0/10.0 * grid->area;
				v_fol[grid->veg_olson] += (mass->plant).mfol[f]*(double)(grid->mm[f])/365.0/10.0 * grid->area;
				v_stm[grid->veg_olson] += (mass->plant).mstm[f]*(double)(grid->mm[f])/365.0/10.0 * grid->area;
				v_rot[grid->veg_olson] += (mass->plant).mrot[f]*(double)(grid->mm[f])/365.0/10.0 * grid->area;
				v_ltr[grid->veg_olson] += (mass->soil).ltr_m[f]*(double)(grid->mm[f])/365.0/10.0 * grid->area;
				v_msl[grid->veg_olson] += (mass->soil).msl_m[f]*(double)(grid->mm[f])/365.0/10.0 * grid->area;
			}
		}
		/* empirical NPP models */
		npp_empirical(grid, loct, flux);
		
		/* biomass burning */
		f_biomassburning(grid, loct, mass, flux);
		for(f=0;f<12;f++){
			h_burnt_area[grid->climy - PIVOT_CLIMY +1] += flux->a_burnt[f]*grid->area;
			h_bioburn_co2[grid->climy - PIVOT_CLIMY +1] += (flux->bb_co2_litter[f]+flux->bb_co2_leaf[f]+flux->bb_co2_wood[f]+flux->bb_co2_root[f]) * grid->area;
			h_bioburn_co[grid->climy - PIVOT_CLIMY +1] += (flux->bb_co_litter[f]+flux->bb_co_leaf[f]+flux->bb_co_wood[f]+flux->bb_co_root[f]) * grid->area;
			h_bioburn_ch4[grid->climy - PIVOT_CLIMY +1] += (flux->bb_ch4_litter[f]+flux->bb_ch4_leaf[f]+flux->bb_ch4_wood[f]+flux->bb_ch4_root[f]) * grid->area;
			h_bioburn_nmhc[grid->climy - PIVOT_CLIMY +1] += (flux->bb_nmhc_litter[f]+flux->bb_nmhc_leaf[f]+flux->bb_nmhc_wood[f]+flux->bb_nmhc_root[f]) * grid->area;
			h_bioburn_oc[grid->climy - PIVOT_CLIMY +1] += (flux->bb_oc_litter[f]+flux->bb_oc_leaf[f]+flux->bb_oc_wood[f]+flux->bb_oc_root[f]) * grid->area;
			h_bioburn_bc[grid->climy - PIVOT_CLIMY +1] += (flux->bb_bc_litter[f]+flux->bb_bc_leaf[f]+flux->bb_bc_wood[f]+flux->bb_bc_root[f]) * grid->area;
			h_bioburn_nox[grid->climy - PIVOT_CLIMY +1] += (flux->bb_bc_litter[f]+flux->bb_bc_leaf[f]+flux->bb_bc_wood[f]+flux->bb_bc_root[f]) * grid->area;
			h_bioburn_so2[grid->climy - PIVOT_CLIMY +1] += (flux->bb_so2_litter[f]+flux->bb_so2_leaf[f]+flux->bb_so2_wood[f]+flux->bb_so2_root[f]) * grid->area;
			h_bioburn_pm25[grid->climy - PIVOT_CLIMY +1] += (flux->bb_pm25_litter[f]+flux->bb_pm25_leaf[f]+flux->bb_pm25_wood[f]+flux->bb_pm25_root[f]) * grid->area;
			h_bioburn_tpm[grid->climy - PIVOT_CLIMY +1] += (flux->bb_tpm_litter[f]+flux->bb_tpm_leaf[f]+flux->bb_tpm_wood[f]+flux->bb_tpm_root[f]) * grid->area;
			h_bioburn_tec[grid->climy - PIVOT_CLIMY +1] += (flux->bb_tec_litter[f]+flux->bb_tec_leaf[f]+flux->bb_tec_wood[f]+flux->bb_tec_root[f]) * grid->area;
		}

		/* erosion */
		f_erosion(grid, loct, echar, mass, flux);
		h_ersn_c[grid->climy - PIVOT_CLIMY +1] += flux->erod_carbon * grid->area;
		
		if(ERSN_CC==1){
			(mass->soil).ltr -= flux->erod_carbon;
			if((mass->soil).ltr<0.0){
				(mass->soil).ltr = 0.0;
			}
		}

		/* land use change */
		f_luc_emit(grid, mass, flux);
		
		h_agrarea[grid->climy - PIVOT_CLIMY +1] += grid->area * grid->f_crop;
		h_luc[grid->climy - PIVOT_CLIMY +1] += (flux->lu_conv + flux->lu_ten + flux->lu_hund) * grid->area;

		/* output */
		/* carbon budget */
		fprintf(fp_carbon,"%ld ", grid->climy);
		fprintf(fp_carbon,"%.2lf ", tmp_ann);
		fprintf(fp_carbon,"%.2lf ", pre_ann);
		fprintf(fp_carbon,"%.2lf ", dswr_ann);
		fprintf(fp_carbon,"%.2lf ", aet_ann);
		fprintf(fp_carbon,"%.2lf ", rof_ann);
		fprintf(fp_carbon,"%.2lf ", lai_max);
		fprintf(fp_carbon,"%.2lf ", fol_ann);
		fprintf(fp_carbon,"%.2lf ", stm_ann);
		fprintf(fp_carbon,"%.2lf ", rot_ann);
		fprintf(fp_carbon,"%.2lf ", ltr_ann);
		fprintf(fp_carbon,"%.2lf ", msl_ann);
		fprintf(fp_carbon,"%.2lf ", gpp_ann);
		fprintf(fp_carbon,"%.2lf ", ar_ann);
		fprintf(fp_carbon,"%.2lf ", hr_ann);
		fprintf(fp_carbon,"%.2lf ", ncb_ann);
		fprintf(fp_carbon,"\n"); /**/
		
		/* nitrogen */
		fprintf(fp_nitrogen,"%ld ", grid->climy);
		for(f=0;f<12;f++){
			fprintf(fp_nitrogen,"%.3lf ", (mass->plant).n_cnpy_m[f]);
			fprintf(fp_nitrogen,"%.3lf ", (mass->plant).n_strg_m[f]);
			fprintf(fp_nitrogen,"%.3lf ", (mass->soil).n_no3_m[f]);
			fprintf(fp_nitrogen,"%.3lf ", (mass->soil).n_nh4_m[f]);
			fprintf(fp_nitrogen,"%.3lf ", (mass->soil).n_mcrb_m[f]);
			fprintf(fp_nitrogen,"%.3lf ", (mass->soil).n_lttr_m[f]);
			fprintf(fp_nitrogen,"%.3lf ", (mass->soil).n_hums_m[f]);
			fprintf(fp_nitrogen,"%.3lf ", (flux->plant).n_biofix[f]);
			fprintf(fp_nitrogen,"%.3lf ", (flux->plant).uptake_no3[f]);
			fprintf(fp_nitrogen,"%.3lf ", (flux->plant).uptake_nh4[f]);
			fprintf(fp_nitrogen,"%.3lf ", (flux->soil).n_leach[f]);
			fprintf(fp_nitrogen,"%.3lf ", (flux->soil).n_nh3vlt[f]);
		}
		fprintf(fp_nitrogen,"\n");

		/* erosion */
		fprintf(fp_ersn,"%ld ", grid->climy);
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
		fprintf(fp_ersn,"%lf ", flux->erod_soil_crop);
		fprintf(fp_ersn,"%lf ", flux->erod_orgmat_crop);
		fprintf(fp_ersn,"%lf ", flux->erod_carbon_crop);
		fprintf(fp_ersn,"%lf ", ltr_ann);
		fprintf(fp_ersn,"%lf ", hrl_ann);
		fprintf(fp_ersn,"%lf ", msl_ann);
		fprintf(fp_ersn,"%lf ", hrm_ann);
		for(f=0;f<12;f++){
			fprintf(fp_ersn,"%lf ", loct->ro2[f]);
		}
		fprintf(fp_ersn,"\n"); /**/
		
		/* GHG & trace gases */
		fprintf(fp_ghg,"%ld ", grid->climy);
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
			fprintf(fp_ghg,"%.3lf ", (flux->plant).emit_ch4_kirschbaum_mass[f]);
			fprintf(fp_ghg,"%.3lf ", (flux->plant).emit_ch4_kirschbaum_photo[f]);
		}
		fprintf(fp_ghg,"\n");

		/* biomass burning */
		fprintf(fp_bioburn,"%ld ", grid->climy);
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
		fprintf(fp_voc,"%ld ", grid->climy);
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
		
		grid->f_crop_p = grid->f_crop;
	}
	fprintf(fp_carbon,"\n");
	fprintf(fp_ersn,"\n");
	fprintf(fp_ghg,"\n");
	fprintf(fp_bioburn,"\n");
	fprintf(fp_voc,"\n");
}
