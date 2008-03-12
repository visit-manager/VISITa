/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* History */
/* Revised November 24, 2007				*/

/* make a result-output file and show values in the console port */
#include<stdio.h>
#include<stdlib.h>
#include"structure.h"
#include"prototype.h"

extern double glandarea;
extern double h_tmp[201], h_pre[201], h_dswr[201], h_aet[201], h_rof[201];
extern double h_gpp[201], h_npp[201], h_nep[201], h_plant[201], h_soil[201];
extern double h_sr[201], h_ersn_c[201], h_agrersn_c[201], h_doc[201];
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
extern double h_nh3_emit[201], h_n2_biofix[201];
extern double h_ch4_emit_mass[201], h_ch4_emit_photo[201];

extern double h_voc_isopr_g97[201], h_voc_monotrp_g97[201], h_voc_methanl_g97[201];
extern double h_voc_acetone_g97[201], h_voc_actaldhd_g97[201], h_voc_frmardhd_g97[201];
extern double h_voc_formacd_g97[201], h_voc_acetacd_g97[201], h_voc_co_g97[201];

/* monthly results **********/
extern double m_ch4ox1[12], m_ch4ox2[12], m_ch4ox3[12];
extern double m_bioburn_co2[12], m_bioburn_ch4[12], m_bioburn_co[12];
extern double m_bioburn_nmhc[12], m_bioburn_oc[12], m_bioburn_bc[12];

/* vegetation (olson) results */
extern double v_area[34];
extern double v_gpp[34], v_npp[34], v_nep[34];
extern double v_lai[34], v_fol[34], v_stm[34], v_rot[34], v_ltr[34], v_msl[34];

/*********** make output file for stable state *************/
void publish_cbud(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux, 
	FILE *result
){
	long f;
	
	fprintf(result,"%.3lf %.3lf %ld %ld ",grid->lat, grid->lon, grid->time, grid->nnn);	
	/* fprintf(result,"%lf %lf %lf %lf ", grid->whc30, grid->whc, grid->sd, grid->hyd_cond); */
	fprintf(result,"\n");
	
	for(f=0;f<12;f++){
		fprintf(result,"%.2lf ", grid->bCO2[f]);	
		fprintf(result,"%lf %lf ", loct->C3ptn[f], loct->C4ptn[f]);
		
		fprintf(result,"%.2lf ", grid->tmp_2m[f]);
		fprintf(result,"%.2lf ", grid->prate_sfc[f]);
		fprintf(result,"%.2lf ", grid->gl_rad[f]);
		fprintf(result,"%.2lf ", loct->rad_net_short[f]);
		fprintf(result,"%.2lf ", loct->rad_net_long[f]);

		fprintf(result,"%.2lf ", loct->evpr[f]);
		fprintf(result,"%.2lf ", loct->incep[f]);
		fprintf(result,"%.2lf ", loct->trspr[f]);
		fprintf(result,"%.2lf ", loct->ro2[f]);
		fprintf(result,"%.2lf ", loct->msw30[f]);
		fprintf(result,"%.2lf ", loct->msww[f]);
		fprintf(result,"%.2lf ", loct->msnwa[f]);

		fprintf(result,"%ld %ld ", (echar->c3).season[f], (echar->c4).season[f]);

		fprintf(result,"%.2lf ", (mass->c3).lai[f]);
		fprintf(result,"%.2lf ", (mass->c3).mfol[f]);
		fprintf(result,"%.2lf ", (mass->c3).mstm[f]);
		fprintf(result,"%.2lf ", (mass->c3).mrot[f]);
		fprintf(result,"%.2lf ", (mass->c4).lai[f]);
		fprintf(result,"%.2lf ", (mass->c4).mfol[f]);
		fprintf(result,"%.2lf ", (mass->c4).mstm[f]);
		fprintf(result,"%.2lf ", (mass->c4).mrot[f]);
		fprintf(result,"%.2lf ", (mass->soil).ltr_m[f]);
		fprintf(result,"%.2lf ", (mass->soil).msl_m[f]);
		
		fprintf(result,"%.2lf ", (flux->c3).gpp[f]);
		fprintf(result,"%.2lf ", (flux->c3).spp[f]);
		fprintf(result,"%.2lf ", (flux->c3).epp[f]);
		fprintf(result,"%.2lf ", (flux->c3).npp[f]);
		fprintf(result,"%.2lf ", (flux->c4).gpp[f]);
		fprintf(result,"%.2lf ", (flux->c4).spp[f]);
		fprintf(result,"%.2lf ", (flux->c4).epp[f]);
		fprintf(result,"%.2lf ", (flux->c4).npp[f]);
		
		fprintf(result,"%.2lf ", (flux->c3).lf[f]);
		fprintf(result,"%.2lf ", (flux->c3).lc[f]);
		fprintf(result,"%.2lf ", (flux->c3).lr[f]);
		fprintf(result,"%.2lf ", (flux->c4).lf[f]);
		fprintf(result,"%.2lf ", (flux->c4).lc[f]);
		fprintf(result,"%.2lf ", (flux->c4).lr[f]);
		
		fprintf(result,"%.2lf ", (flux->soil).rl[f]);
		fprintf(result,"%.2lf ", (flux->soil).rh[f]);
		fprintf(result,"%.2lf ", (flux->soil).sf[f]);
		fprintf(result,"%.2lf ", (flux->soil).rS[f]);
		
		fprintf(result,"%.2lf ", flux->nep[f]);
		fprintf(result,"%.2lf ", flux->ncb[f]);
		
	/*	fprintf(result,"%lf %lf %lf %lf ",grid->tmp_sfc[f],grid->tmp_2m[f],grid->tcdc_clm[f],grid->wnd_10m[f]);
		fprintf(result,"%lf %lf %lf %lf %lf ",grid->prate_sfc[f],grid->spfh_2m[f],loct->vp[f],loct->vpd[f],loct->slope_vps[f]);
		fprintf(result,"%lf %lf %lf ",grid->dlen[f],grid->gl_rad[f],loct->gl_rad_g[f]);
		fprintf(result,"%lf %lf %lf ",loct->rad_net_g[f],loct->rad_net_p[f],loct->rad_net_long[f]);
		fprintf(result,"%.2lf %.2lf %.2lf %.2lf ", grid->tmp_2m[f], grid->prate_sfc[f], grid->par[f], loct->rad_net_g[f]+loct->rad_net_p[f]);
		fprintf(result,"%lf %lf %lf %lf %lf ",grid->par[f],grid->par_be[f],grid->par_de[f],grid->par_bp[f],grid->par_dp[f]);
		fprintf(result,"%lf %lf %lf %lf ",(echar->c3).albedo,(echar->c4).albedo,(echar->soil).albedo[f],loct->albedo_sfc[f]);	 
		fprintf(result,"%.2lf %.2lf %.2lf ",loct->canopy_con[f],(echar->c3).gs[f],(echar->c4).gs[f]);
		fprintf(result,"%.2lf %.2lf %.2lf ",loct->pm_evp[f],loct->pm_incep[f],loct->pm_trn[f]);
		fprintf(result,"%.2lf %.2lf %.2lf %.2lf ",(echar->c3).opt_lai[f],(echar->c3).psat[f],(echar->c3).eK[f],(echar->c3).lue[f]);
		fprintf(result,"%.2lf %.2lf %.2lf %.2lf %.2lf ",(echar->c3).ci[f],(echar->c3).gs[f],(echar->c3).gc[f],(echar->c3).qTf[f],(echar->c3).lf[f]);
		fprintf(result,"%.2lf %.2lf %.2lf %.2lf ",(echar->c4).opt_lai[f],(echar->c4).psat[f],(echar->c4).eK[f],(echar->c4).lue[f]);
		fprintf(result,"%.2lf %.2lf %.2lf %.2lf %.2lf ",(echar->c4).ci[f],(echar->c4).gs[f],(echar->c4).gc[f],(echar->c4).qTf[f],(echar->c4).lf[f]);
		fprintf(result,"%lf %lf ",(echar->c3).ft[f],(echar->c4).ft[f]);
		fprintf(result,"%lf %lf ",(echar->c3).fcd[f],(echar->c4).fcd[f]);
		fprintf(result,"%lf %lf ",(echar->c3).fsw[f],(echar->c4).fsw[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c3).rfm[f],(flux->c3).rcm[f],(flux->c3).rrm[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c3).rfg[f],(flux->c3).rcg[f],(flux->c3).rrg[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c3).rpm[f],(flux->c3).rpg[f],(flux->c3).rp[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c3).tpf[f],(flux->c3).tpc[f],(flux->c3).tpr[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c4).rfm[f],(flux->c4).rcm[f],(flux->c4).rrm[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c4).rfg[f],(flux->c4).rcg[f],(flux->c4).rrg[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c4).rpm[f],(flux->c4).rpg[f],(flux->c4).rp[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c4).tpf[f],(flux->c4).tpc[f],(flux->c4).tpr[f]); */

		fprintf(result,"\n");	
	}
	fprintf(result,"%lf ", flux->npp_miami); 
	fprintf(result,"%lf ", flux->npp_montreal); 
	fprintf(result,"%.2lf ", flux->erod_soil);
	fprintf(result,"%.2lf ", flux->erod_orgmat);
	fprintf(result,"%.2lf ", flux->erod_carbon);
	
/*	for(f=0;f<12;f++){
		fprintf(result,"%lf %lf ",grid->d13C_bCO2[f],loct->d13C_aCO2[f]);

		fprintf(result,"%lf %lf ",(echar->c3).photo_13c_frac[f], (flux->c3).d13c_gpp[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c3).d13c_rfm[f], (flux->c3).d13c_rcm[f], (flux->c3).d13c_rrm[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c3).d13c_rfg[f], (flux->c3).d13c_rcg[f], (flux->c3).d13c_rrg[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c3).d13c_tpf[f], (flux->c3).d13c_tpc[f], (flux->c3).d13c_tpr[f]);
		fprintf(result,"%lf %lf %lf ",(mass->c3).d13c_mfol[f], (mass->c3).d13c_mstm[f], (mass->c3).d13c_mrot[f]);

		fprintf(result,"%lf %lf ",(echar->c4).photo_13c_frac[f], (flux->c4).d13c_gpp[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c4).d13c_rfm[f], (flux->c4).d13c_rcm[f], (flux->c4).d13c_rrm[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c4).d13c_rfg[f], (flux->c4).d13c_rcg[f], (flux->c4).d13c_rrg[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c4).d13c_tpf[f], (flux->c4).d13c_tpc[f], (flux->c4).d13c_tpr[f]);
		fprintf(result,"%lf %lf %lf ",(mass->c4).d13c_mfol[f], (mass->c4).d13c_mstm[f], (mass->c4).d13c_mrot[f]);

		fprintf(result,"%lf %lf %lf ",(flux->plant).d13c_lL[f],(flux->soil).d13c_lL[f],(flux->soil).d13c_sf[f]);

		fprintf(result,"%lf %lf ",(flux->soil).d13c_rl[f], (flux->soil).d13c_rh[f]);
		fprintf(result,"%lf %lf ",(mass->soil).d13c_mltr[f], (mass->soil).d13c_mmsl[f]);

		fprintf(result,"%lf %lf ",(flux->c3).d13c_rp[f], (flux->c4).d13c_rp[f]);
		fprintf(result,"%lf %lf ",(flux->c3).d13c_lL[f], (flux->c4).d13c_lL[f]);
		fprintf(result,"%lf ",(flux->soil).d13c_rS[f]);

		fprintf(result," \n");
	} */
	
	fprintf(result," \n"); 	
}

/******************************************************************/
void f_output_ann(
	long year, 
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
	long f;
	double gpp_ann, ar_ann, npp_ann, hr_ann, nep_ann, ncb_ann, lf_ann, hrl_ann, hrm_ann;
	double fol_ann, stm_ann, rot_ann, ltr_ann, msl_ann, lai_max;
	double tmp_ann, pre_ann, dswr_ann, ipar_ann, apar_ann, aet_ann, rof_ann, rn_ann;
	double iparb_ann, ipard_ann, aparb_ann, apard_ann, ls_gpp, aparmono_ann;
	double ch4ox1_ann, ch4ox2_ann, ch4ox3_ann, sw1_ann, wfps_ann;
	
	gpp_ann = ar_ann = npp_ann = hr_ann = nep_ann = ncb_ann = lf_ann = hrl_ann = hrm_ann = 0.0;
	fol_ann = stm_ann = rot_ann = ltr_ann = msl_ann = lai_max = 0.0;
	tmp_ann = pre_ann = dswr_ann = ipar_ann = apar_ann = aet_ann = rof_ann = rn_ann = 0.0;
	iparb_ann = ipard_ann = aparb_ann = apard_ann = ls_gpp = aparmono_ann = 0.0;
	ch4ox1_ann = ch4ox2_ann = ch4ox3_ann = sw1_ann = wfps_ann = 0.0;

	/* carbon budget */
	fprintf(fp_carbon,"%ld ", year);
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
	fprintf(fp_nitrogen,"%ld ", year);
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
	fprintf(fp_ersn,"%ld ", year);
	fprintf(fp_ersn,"%lf ", pre_ann);
	fprintf(fp_ersn,"%lf ", grid->f_crop_con);
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
	fprintf(fp_ghg,"%ld ", year);
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
	fprintf(fp_bioburn,"%ld ", year);
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
	fprintf(fp_voc,"%ld ", year);
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

/*****************************************************************/
void f_set_history_data(
	long year, 
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Flux *flux
){
	long f;
	
	for(f=0;f<12;f++){
		h_tmp[year] += grid->tmp_2m[f]*(double)(grid->mm[f])/365.0 * grid->area;
		h_pre[year] += grid->prate_sfc[f]*0.1 * grid->area;
		h_dswr[year] += grid->gl_rad[f]*(double)(grid->mm[f])/365.0 * grid->area;
		h_aet[year] += (loct->evpr[f] + loct->trspr[f] + loct->incep[f])*0.1 * grid->area;
		h_rof[year] += loct->ro2[f]*0.1 * grid->area;
		h_gpp[year] += (flux->plant).gpp[f] * grid->area;
		h_npp[year] += (flux->plant).npp[f] * grid->area;
		h_nep[year] += flux->nep[f] * grid->area;
		h_plant[year] += ((mass->plant).mfol[f] + (mass->plant).mstm[f] + 
			(mass->plant).mrot[f])*(double)(grid->mm[f])/365.0 * grid->area;
		h_soil[year] += ((mass->soil).ltr_m[f] + 
			(mass->soil).msl_m[f])*(double)(grid->mm[f])/365.0 * grid->area;
		h_sr[year] += ((flux->plant).rrm[f] + (flux->plant).rrg[f] + (flux->soil).rS[f]) * grid->area;
		h_doc[year] += (flux->soil).doc_boyer[f] * grid->area;
		
		/* GHG */
		h_ch4ox1[year] += (flux->soil).ch4oxy_ridg[f] * grid->area *10000.0/1000.0;
		h_ch4ox2[year] += (flux->soil).ch4oxy_casa[f] * grid->area *10000.0/1000.0;
		h_ch4ox3[year] += (flux->soil).ch4oxy_delgrosso[f] * grid->area *10000.0/1000.0;
		h_ch4emit_cao_paddy[year] += (flux->soil).ch4emit_paddy_cao[f] * grid->area *10000.0/1000.0;
		h_ch4emit_cao_wetland[year] += (flux->soil).ch4emit_wetland_cao[f] * grid->area *10000.0/1000.0;
		h_n2o_emit_ngas[year] += (flux->soil).d_n2o_ngas[f] * grid->area;
		h_n2_emit_ngas[year] += (flux->soil).d_n2_ngas[f] * grid->area;
		h_n2o_emit_casa[year] += (flux->soil).d_n2o_casa[f] * grid->area;
		h_n2_emit_casa[year] += (flux->soil).d_n2_casa[f] * grid->area;
		h_no_emit_casa[year] += (flux->soil).d_no_casa[f] * grid->area;
		h_nh3_emit[year] += (flux->soil).n_nh3vlt[f] * grid->area;		
		h_ch4_emit_mass[year] += (flux->plant).emit_ch4_kirschbaum_mass[f]*10000.0*grid->area;
		h_ch4_emit_photo[year] += (flux->plant).emit_ch4_kirschbaum_photo[f]*10000.0*grid->area;
		h_n2_biofix[year] += (flux->plant).n_biofix[f]*grid->area;
		
		/* biomass burninig */
		h_burnt_area[year] += flux->a_burnt[f]*grid->area;
		h_bioburn_co2[year] += (flux->bb_co2_litter[f]+flux->bb_co2_leaf[f]+flux->bb_co2_wood[f]+flux->bb_co2_root[f]) * grid->area;
		h_bioburn_co[year] += (flux->bb_co_litter[f]+flux->bb_co_leaf[f]+flux->bb_co_wood[f]+flux->bb_co_root[f]) * grid->area;
		h_bioburn_ch4[year] += (flux->bb_ch4_litter[f]+flux->bb_ch4_leaf[f]+flux->bb_ch4_wood[f]+flux->bb_ch4_root[f]) * grid->area;
		h_bioburn_nmhc[year] += (flux->bb_nmhc_litter[f]+flux->bb_nmhc_leaf[f]+flux->bb_nmhc_wood[f]+flux->bb_nmhc_root[f]) * grid->area;
		h_bioburn_oc[year] += (flux->bb_oc_litter[f]+flux->bb_oc_leaf[f]+flux->bb_oc_wood[f]+flux->bb_oc_root[f]) * grid->area;
		h_bioburn_bc[year] += (flux->bb_bc_litter[f]+flux->bb_bc_leaf[f]+flux->bb_bc_wood[f]+flux->bb_bc_root[f]) * grid->area;
		h_bioburn_nox[year] += (flux->bb_bc_litter[f]+flux->bb_bc_leaf[f]+flux->bb_bc_wood[f]+flux->bb_bc_root[f]) * grid->area;
		h_bioburn_so2[year] += (flux->bb_so2_litter[f]+flux->bb_so2_leaf[f]+flux->bb_so2_wood[f]+flux->bb_so2_root[f]) * grid->area;
		h_bioburn_pm25[year] += (flux->bb_pm25_litter[f]+flux->bb_pm25_leaf[f]+flux->bb_pm25_wood[f]+flux->bb_pm25_root[f]) * grid->area;
		h_bioburn_tpm[year] += (flux->bb_tpm_litter[f]+flux->bb_tpm_leaf[f]+flux->bb_tpm_wood[f]+flux->bb_tpm_root[f]) * grid->area;
		h_bioburn_tec[year] += (flux->bb_tec_litter[f]+flux->bb_tec_leaf[f]+flux->bb_tec_wood[f]+flux->bb_tec_root[f]) * grid->area;

		/* VOC, g C month-1 */
		h_voc_isopr_g97[year] += flux->voc_isopr_g97[f] * grid->area *10000.0/1000000.0; 
		h_voc_monotrp_g97[year] += flux->voc_monotrp_g97[f] * grid->area *10000.0/1000000.0; 
		h_voc_methanl_g97[year] += flux->voc_methanl_g97[f] * grid->area *10000.0/1000000.0; 
		h_voc_acetone_g97[year] += flux->voc_acetone_g97[f] * grid->area *10000.0/1000000.0; 
		h_voc_actaldhd_g97[year] += flux->voc_actaldhd_g97[f] * grid->area *10000.0/1000000.0; 
		h_voc_frmardhd_g97[year] += flux->voc_frmardhd_g97[f] * grid->area *10000.0/1000000.0; 
		h_voc_formacd_g97[year] += flux->voc_formacd_g97[f] * grid->area *10000.0/1000000.0; 
		h_voc_acetacd_g97[year] += flux->voc_acetacd_g97[f] * grid->area *10000.0/1000000.0; 
		h_voc_co_g97[year] += flux->voc_co_g97[f] * grid->area *10000.0/1000000.0; 
	}
	
	/* erosion */
	h_ersn_c[year] += flux->erod_carbon * grid->area;
	h_agrersn_c[year] += flux->erod_carbon_crop * grid->area;

	h_agrarea[year] += grid->area * grid->f_crop_con;
	h_luc[year] += (flux->lu_conv + flux->lu_ten + flux->lu_hund) * grid->area;
}

/*********************************************************************************/
void f_output_result(
	long year, 
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux, 
	FILE *fp_o[7]
){
	long f;
	double pre_ann, hrl_ann, hrm_ann, ltr_ann, msl_ann;
	
	pre_ann = hrl_ann = hrm_ann = ltr_ann = msl_ann = 0.0;
	
	for(f=0;f<12;f++){
		pre_ann += grid->prate_sfc[f];
		
		hrl_ann += (flux->soil).rl[f];
		hrm_ann += (flux->soil).rh[f];

		ltr_ann += (mass->soil).ltr_m[f]*(double)(grid->mm[f])/365.0;
		msl_ann += (mass->soil).msl_m[f]*(double)(grid->mm[f])/365.0;
	}

	/* carbon budget */
	fprintf(fp_o[0],"%ld ", year);
	for(f=0;f<12;f++){
		fprintf(fp_o[0],"%.2lf ", grid->tmp_2m[f]);
		fprintf(fp_o[0],"%.2lf ", grid->prate_sfc[f]);
		fprintf(fp_o[0],"%.2lf ", grid->gl_rad[f]);
		
		fprintf(fp_o[0],"%.2lf ", loct->evpr[f] + loct->trspr[f] + loct->incep[f]);
		fprintf(fp_o[0],"%.2lf ", loct->ro2[grid->m]);
		
		fprintf(fp_o[0],"%.2lf ", loct->C3ptn[f]*(mass->c3).lai[f]+loct->C4ptn[f]*(mass->c4).lai[f]);
		fprintf(fp_o[0],"%.2lf ", (mass->plant).mfol[f]);
		fprintf(fp_o[0],"%.2lf ", (mass->plant).mstm[f]);
		fprintf(fp_o[0],"%.2lf ", (mass->plant).mrot[f]);
		fprintf(fp_o[0],"%.2lf ", (mass->soil).ltr_m[f]);
		fprintf(fp_o[0],"%.2lf ", (mass->soil).msl_m[f]);
		
		fprintf(fp_o[0],"%.2lf ", (flux->plant).gpp[f]);
		fprintf(fp_o[0],"%.2lf ", (flux->plant).rp[f]);
		fprintf(fp_o[0],"%.2lf ", (flux->soil).rS[f]);
		fprintf(fp_o[0],"%.2lf ",  flux->nep[f]);
		fprintf(fp_o[0],"%.2lf ",  flux->ncb[f]);
	}
	fprintf(fp_o[0],"%.2lf ", flux->lu_detr);
	fprintf(fp_o[0],"%.2lf ", flux->lu_conv);
	fprintf(fp_o[0],"%.2lf ", flux->lu_ten);
	fprintf(fp_o[0],"%.2lf ", flux->lu_hund);
	fprintf(fp_o[0],"\n"); /**/
	
	/* nitrogen */
	fprintf(fp_o[1],"%ld ", year);
	for(f=0;f<12;f++){
		fprintf(fp_o[1],"%.3lf ", (mass->plant).n_cnpy_m[f]);
		fprintf(fp_o[1],"%.3lf ", (mass->plant).n_strg_m[f]);
		fprintf(fp_o[1],"%.3lf ", (mass->soil).n_no3_m[f]);
		fprintf(fp_o[1],"%.3lf ", (mass->soil).n_nh4_m[f]);
		fprintf(fp_o[1],"%.3lf ", (mass->soil).n_mcrb_m[f]);
		fprintf(fp_o[1],"%.3lf ", (mass->soil).n_lttr_m[f]);
		fprintf(fp_o[1],"%.3lf ", (mass->soil).n_hums_m[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->plant).n_biofix[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->plant).uptake_no3[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->plant).uptake_nh4[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->soil).n_leach[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->soil).n_nh3vlt[f]);
	}
	fprintf(fp_o[1],"\n");

	/* erosion & DOC */
	fprintf(fp_o[2],"%ld ", year);
	fprintf(fp_o[2],"%lf ", pre_ann);
	fprintf(fp_o[2],"%lf ", grid->f_crop_con);
	fprintf(fp_o[2],"%lf ", grid->f_erosion_r);
	fprintf(fp_o[2],"%lf ", grid->f_erosion_ls);
	fprintf(fp_o[2],"%lf ", grid->f_erosion_k);
	fprintf(fp_o[2],"%lf ", grid->f_erosion_c);
	fprintf(fp_o[2],"%lf ", grid->f_erosion_p);
	fprintf(fp_o[2],"%lf ", flux->erod_soil);
	fprintf(fp_o[2],"%lf ", flux->erod_orgmat);
	fprintf(fp_o[2],"%lf ", flux->erod_carbon);		
	fprintf(fp_o[2],"%lf ", flux->erod_soil_crop);
	fprintf(fp_o[2],"%lf ", flux->erod_orgmat_crop);
	fprintf(fp_o[2],"%lf ", flux->erod_carbon_crop);
	fprintf(fp_o[2],"%.3lf ", ltr_ann);
	fprintf(fp_o[2],"%.3lf ", hrl_ann);
	fprintf(fp_o[2],"%.3lf ", msl_ann);
	fprintf(fp_o[2],"%.3lf ", hrm_ann);
	for(f=0;f<12;f++){
		fprintf(fp_o[2],"%lf ", loct->ro2[f]);
		fprintf(fp_o[2],"%lf ", (flux->soil).doc_boyer[f]);
	}
	fprintf(fp_o[2],"\n"); /**/
	
	/* GHG & trace gases */
	fprintf(fp_o[3],"%ld ", year);
	for(f=0;f<12;f++){
		fprintf(fp_o[3],"%.3lf ", loct->msww[f]);
		fprintf(fp_o[3],"%.3lf ", loct->wfps[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4oxy_ridg[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4oxy_casa[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4oxy_delgrosso[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4emit_paddy_cao[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4emit_wetland_cao[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).d_n2o_ngas[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).d_n2_ngas[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).d_n2o_casa[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).d_n2_casa[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).d_no_casa[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).n_nh3vlt[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->plant).emit_ch4_kirschbaum_mass[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->plant).emit_ch4_kirschbaum_photo[f]);
	}
	fprintf(fp_o[3],"\n");

	/* biomass burning */
	fprintf(fp_o[4],"%ld ", year);
	for(f=0;f<12;f++){
		fprintf(fp_o[4],"%.3lf ", loct->msw30[f]);
		fprintf(fp_o[4],"%.3lf ", flux->day_fire[f]);
		fprintf(fp_o[4],"%.3lf ", flux->a_burnt[f]);
		fprintf(fp_o[4],"%.3lf ", flux->bb_co2_litter[f]+flux->bb_co2_leaf[f]+flux->bb_co2_wood[f]+flux->bb_co2_root[f]);
		fprintf(fp_o[4],"%.3lf ", flux->bb_co_litter[f]+flux->bb_co_leaf[f]+flux->bb_co_wood[f]+flux->bb_co_root[f]);
		fprintf(fp_o[4],"%.3lf ", flux->bb_ch4_litter[f]+flux->bb_ch4_leaf[f]+flux->bb_ch4_wood[f]+flux->bb_ch4_root[f]);
		fprintf(fp_o[4],"%.3lf ", flux->bb_nmhc_litter[f]+flux->bb_nmhc_leaf[f]+flux->bb_nmhc_wood[f]+flux->bb_nmhc_root[f]);
		fprintf(fp_o[4],"%.3lf ", flux->bb_oc_litter[f]+flux->bb_oc_leaf[f]+flux->bb_oc_wood[f]+flux->bb_oc_root[f]);
		fprintf(fp_o[4],"%.3lf ", flux->bb_bc_litter[f]+flux->bb_bc_leaf[f]+flux->bb_bc_wood[f]+flux->bb_bc_root[f]);
		fprintf(fp_o[4],"%.3lf ", flux->bb_nox_litter[f]+flux->bb_nox_leaf[f]+flux->bb_nox_wood[f]+flux->bb_nox_root[f]);
		fprintf(fp_o[4],"%.3lf ", flux->bb_so2_litter[f]+flux->bb_so2_leaf[f]+flux->bb_so2_wood[f]+flux->bb_so2_root[f]);
		fprintf(fp_o[4],"%.3lf ", flux->bb_pm25_litter[f]+flux->bb_pm25_leaf[f]+flux->bb_pm25_wood[f]+flux->bb_pm25_root[f]);
		fprintf(fp_o[4],"%.3lf ", flux->bb_tpm_litter[f]+flux->bb_tpm_leaf[f]+flux->bb_tpm_wood[f]+flux->bb_tpm_root[f]);
		fprintf(fp_o[4],"%.3lf ", flux->bb_tec_litter[f]+flux->bb_tec_leaf[f]+flux->bb_tec_wood[f]+flux->bb_tec_root[f]);
	}
	fprintf(fp_o[4],"\n");

	/* VOC */
	fprintf(fp_o[5],"%ld ", year);
	for(f=0;f<12;f++){
		fprintf(fp_o[5],"%.3lf ", flux->voc_isopr_g97[f]);
		fprintf(fp_o[5],"%.3lf ", flux->voc_monotrp_g97[f]);
		fprintf(fp_o[5],"%.3lf ", flux->voc_methanl_g97[f]);
		fprintf(fp_o[5],"%.3lf ", flux->voc_acetone_g97[f]);
		fprintf(fp_o[5],"%.3lf ", flux->voc_actaldhd_g97[f]);
		fprintf(fp_o[5],"%.3lf ", flux->voc_frmardhd_g97[f]);
		fprintf(fp_o[5],"%.3lf ", flux->voc_formacd_g97[f]);
		fprintf(fp_o[5],"%.3lf ", flux->voc_acetacd_g97[f]);
		fprintf(fp_o[5],"%.3lf ", flux->voc_co_g97[f]);
	}
	fprintf(fp_o[5],"\n");
	
	/* hydrometeorology */
	fprintf(fp_o[6],"%ld ", year);
	for(f=0;f<12;f++){
		fprintf(fp_o[6],"%.3lf ", grid->tmp_2m[f]);
		fprintf(fp_o[6],"%.2lf ", grid->prate_sfc[f]);
		fprintf(fp_o[6],"%.2lf ", grid->gl_rad[f]);
		fprintf(fp_o[6],"%.2lf ", grid->par[f]);
		fprintf(fp_o[6],"%.2lf ", grid->par_bp[f]);
		fprintf(fp_o[6],"%.2lf ", grid->par_dp[f]);
		fprintf(fp_o[6],"%.2lf ", loct->vpd[f]);
		fprintf(fp_o[6],"%.2lf ", loct->msw30[f]);
		fprintf(fp_o[6],"%.2lf ", loct->msww[f]);
		fprintf(fp_o[6],"%.2lf ", loct->incep[f]);
		fprintf(fp_o[6],"%.2lf ", loct->evpr[f]);
		fprintf(fp_o[6],"%.2lf ", loct->trspr[f]);
		fprintf(fp_o[6],"%.2lf ", loct->ro2[f]);
	}
	fprintf(fp_o[6],"\n");
}

/*******************************************************************************/
void f_glosum_output(
	char s_date[25], 
	char s_case[25]
){
	short h;
	char filename[100];
	FILE *fp_glsum;

	/* global summary ******************************************/
	strcpy(filename, s_date);
	strcat(filename, s_case);
	strcat(filename, "glsum.dat");
	
	fp_glsum = fopen(filename,"wt");
	fprintf(fp_glsum,"%lf\n", glandarea);
	for(h=0;h<201;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		fprintf(fp_glsum,"%lf ", h_tmp[h]);
		fprintf(fp_glsum,"%lf ", h_pre[h]);
		fprintf(fp_glsum,"%lf ", h_dswr[h]);
		fprintf(fp_glsum,"%lf ", h_aet[h]);
		fprintf(fp_glsum,"%lf ", h_rof[h]);
		
		fprintf(fp_glsum,"%lf ", h_gpp[h]);
		fprintf(fp_glsum,"%lf ", h_npp[h]);
		fprintf(fp_glsum,"%lf ", h_nep[h]);
		fprintf(fp_glsum,"%lf ", h_sr[h]);
		fprintf(fp_glsum,"%lf ", h_plant[h]);
		fprintf(fp_glsum,"%lf ", h_soil[h]);
		
		fprintf(fp_glsum,"%lf ", h_ersn_c[h]);
		fprintf(fp_glsum,"%lf ", h_agrersn_c[h]);
		fprintf(fp_glsum,"%lf ", h_doc[h]);
		
		fprintf(fp_glsum,"%lf ", h_agrarea[h]);
		fprintf(fp_glsum,"%lf ", h_luc[h]);
		
		fprintf(fp_glsum,"%lf ", h_burnt_area[h]);
		fprintf(fp_glsum,"%lf ", h_bioburn_co2[h]);
		fprintf(fp_glsum,"%lf ", h_bioburn_co[h]);
		fprintf(fp_glsum,"%lf ", h_bioburn_ch4[h]);
		fprintf(fp_glsum,"%lf ", h_bioburn_nmhc[h]);
		fprintf(fp_glsum,"%lf ", h_bioburn_oc[h]);
		fprintf(fp_glsum,"%lf ", h_bioburn_bc[h]);
		fprintf(fp_glsum,"%lf ", h_bioburn_nox[h]);
		fprintf(fp_glsum,"%lf ", h_bioburn_so2[h]);
		fprintf(fp_glsum,"%lf ", h_bioburn_pm25[h]);
		fprintf(fp_glsum,"%lf ", h_bioburn_tpm[h]);
		fprintf(fp_glsum,"%lf ", h_bioburn_tec[h]);
		
		fprintf(fp_glsum,"%lf ", h_n2o_emit_ngas[h]);
		fprintf(fp_glsum,"%lf ", h_n2_emit_ngas[h]);
		fprintf(fp_glsum,"%lf ", h_n2o_emit_casa[h]);
		fprintf(fp_glsum,"%lf ", h_n2_emit_casa[h]);
		fprintf(fp_glsum,"%lf ", h_no_emit_casa[h]);
		fprintf(fp_glsum,"%lf ", h_nh3_emit[h]);
		fprintf(fp_glsum,"%lf ", h_n2_biofix[h]);
		fprintf(fp_glsum,"%lf ", h_ch4ox1[h]);
		fprintf(fp_glsum,"%lf ", h_ch4ox2[h]);
		fprintf(fp_glsum,"%lf ", h_ch4ox3[h]);
		fprintf(fp_glsum,"%lf ", h_ch4emit_cao_paddy[h]);
		fprintf(fp_glsum,"%lf ", h_ch4emit_cao_wetland[h]);
		fprintf(fp_glsum,"%lf ", h_ch4_emit_mass[h]);
		fprintf(fp_glsum,"%lf ", h_ch4_emit_photo[h]);
	
		fprintf(fp_glsum,"%lf ", h_voc_isopr_g97[h]);
		fprintf(fp_glsum,"%lf ", h_voc_monotrp_g97[h]);
		fprintf(fp_glsum,"%lf ", h_voc_methanl_g97[h]);
		fprintf(fp_glsum,"%lf ", h_voc_acetone_g97[h]);
		fprintf(fp_glsum,"%lf ", h_voc_actaldhd_g97[h]);
		fprintf(fp_glsum,"%lf ", h_voc_frmardhd_g97[h]);
		fprintf(fp_glsum,"%lf ", h_voc_formacd_g97[h]);
		fprintf(fp_glsum,"%lf ", h_voc_acetacd_g97[h]);
		fprintf(fp_glsum,"%lf ", h_voc_co_g97[h]);

		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<12;h++){
		fprintf(fp_glsum,"%ld ", h);
		fprintf(fp_glsum,"%lf ", m_ch4ox1[h]);
		fprintf(fp_glsum,"%lf ", m_ch4ox2[h]);
		fprintf(fp_glsum,"%lf ", m_ch4ox3[h]);
		
		fprintf(fp_glsum,"%lf ", m_bioburn_co2[h]);
		fprintf(fp_glsum,"%lf ", m_bioburn_co[h]);
		fprintf(fp_glsum,"%lf ", m_bioburn_ch4[h]);
		fprintf(fp_glsum,"%lf ", m_bioburn_nmhc[h]);
		fprintf(fp_glsum,"%lf ", m_bioburn_oc[h]);
		fprintf(fp_glsum,"%lf ", m_bioburn_bc[h]);
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<34;h++){
		fprintf(fp_glsum,"%lf ", v_area[h]);
		fprintf(fp_glsum,"%lf ", v_gpp[h]);
		fprintf(fp_glsum,"%lf ", v_npp[h]);
		fprintf(fp_glsum,"%lf ", v_nep[h]);
		fprintf(fp_glsum,"%lf ", v_lai[h]);
		fprintf(fp_glsum,"%lf ", v_fol[h]);
		fprintf(fp_glsum,"%lf ", v_stm[h]);
		fprintf(fp_glsum,"%lf ", v_rot[h]);
		fprintf(fp_glsum,"%lf ", v_ltr[h]);
		fprintf(fp_glsum,"%lf ", v_msl[h]);
		fprintf(fp_glsum,"\n");
	}
	fclose(fp_glsum);
}

/*******************************************************************************/
void f_output_file_open(
	short zone,
	char s_date[25], 
	char s_case[25], 
	char filename[100], 
	FILE *fp[6]
){
	char num[4];
	
	strcpy(filename, s_date);	
	strcat(filename, s_case);
	strcat(filename, "crbn_");
	snprintf(num, 4, "%02d", zone);
	strcat(filename, num);
	strcat(filename, ".dat");
	if( (fp[0] = fopen(filename,"wt"))==NULL ){  
		printf("No carbon output file\n");  
		exit(1); 
	}

	strcpy(filename, s_date);
	strcat(filename, s_case);
	strcat(filename, "nitr_");
	snprintf(num, 4, "%02d", zone);
	strcat(filename, num);
	strcat(filename, ".dat");
	if( (fp[1] = fopen(filename,"wt"))==NULL ){  
		printf("No nitrogen output file\n");  
		exit(1); 
	}

	strcpy(filename, s_date);
	strcat(filename, s_case);
	strcat(filename, "ersn_");
	snprintf(num, 4, "%02d", zone);
	strcat(filename, num);
	strcat(filename, ".dat");
	if( (fp[2] = fopen(filename,"wt"))==NULL ){  
		printf("No erosion output file\n");  
		exit(1); 
	}

	strcpy(filename, s_date);
	strcat(filename, s_case);
	strcat(filename, "ghg_");
	snprintf(num, 4, "%02d", zone);
	strcat(filename, num);
	strcat(filename, ".dat");
	if( (fp[3] = fopen(filename,"wt"))==NULL ){  
		printf("No GHG output file\n");  
		exit(1); 
	}

	strcpy(filename, s_date);
	strcat(filename, s_case);
	strcat(filename, "bb_");
	snprintf(num, 4, "%02d", zone);
	strcat(filename, num);
	strcat(filename, ".dat");
	if( (fp[4] = fopen(filename,"wt"))==NULL ){  
		printf("No biomass burning output file\n");  
		exit(1); 
	}

	strcpy(filename, s_date);
	strcat(filename, s_case);
	strcat(filename, "bvoc_");
	snprintf(num, 4, "%02d", zone);
	strcat(filename, num);
	strcat(filename, ".dat");
	if( (fp[5] = fopen(filename,"wt"))==NULL ){  
		printf("No BVOC output file\n");  
		exit(1); 
	}

	strcpy(filename, s_date);
	strcat(filename, s_case);
	strcat(filename, "hydmet_");
	snprintf(num, 4, "%02d", zone);
	strcat(filename, num);
	strcat(filename, ".dat");
	if( (fp[6] = fopen(filename,"wt"))==NULL ){  
		printf("No Hydrometorological output file\n");  
		exit(1); 
	}
}
