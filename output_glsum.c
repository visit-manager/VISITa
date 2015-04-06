/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* History */
/* Revised November 24, 2007				*/

/* make a result-output file and show values in the console port */
#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include"structure.h"
#include"prototype.h"

/*****************************************************************/
/* time-series of global total values */
void f_set_history_data(
	long year, 
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Flux *flux
){
	long f;
	extern double MDN[ASTEP];
	double fweight, fweight_nat, fweight_wet, fweight_pad;
	
	/* 2010/04/27 by A.Ito ***********/
	/* natural */
	if(loct->v_type == 1){
		if(REPL_OLSON_CROP == 1){
			fweight = 1.0 - grid->f_crop_con;
			
			if(NECB_LUC == 0 && (EX_CCPL != 3 && EX_CCPL != 8)){
				fweight = 1.0 - grid->fcrop_unh_hmnzed[200];
			}
			
			fweight_nat = 1.0;
			fweight_wet = 1.0;
			fweight_pad = 0.0;
		}else{
			fweight = 1.0;
			
			fweight_nat = 1.0;
			fweight_wet = 1.0;
			fweight_pad = 1.0;
		}
	}
	/* cropland */
	if(loct->v_type == 2){
		fweight = grid->f_crop_con;
		
		fweight_nat = 0.0;
		fweight_wet = 0.0;
		fweight_pad = 1.0;
		if(LANDUSE == 0){
			fweight_pad = 0.0;
		}
	}
	
	for(f=0;f<ASTEP;f++){
		/* climate */
		//h_tmp[year] += fweight * grid->tmp_2m[f]* MDN[f]/365.0 * grid->area;
		h_tmp[year] += fweight * grid->tmp_sfc[f]* MDN[f]/365.0 * grid->area;
		h_pre[year] += fweight * grid->prate_sfc[f] * grid->area;
		h_dswr[year] += fweight * grid->gl_rad[f]* MDN[f]/365.0 * grid->area;
		h_aet[year] += fweight * (loct->evpr[f] + loct->trspr[f] + loct->incep[f]) * grid->area;
		h_rof[year] += fweight * loct->ro2[f] * grid->area;
		
		h_trnsp[year] += fweight * loct->trspr[f] * grid->area;
		h_incepev[year] += fweight * loct->incep[f] * grid->area;
		h_ssurfev[year] += fweight * loct->evpr[f] * grid->area;
		
		h_sw1[year] += fweight * loct->msw30[f]* MDN[f]/365.0 * grid->area;
		h_sw2[year] += fweight * loct->msww[f]* MDN[f]/365.0 * grid->area;
        
        h_rns[year] += fweight * loct->rad_net_short[f]* MDN[f]/365.0 * grid->area;
        h_rnl[year] += fweight * loct->rad_net_long[f]* MDN[f]/365.0 * grid->area;
        h_rnsd[year] += fweight * loct->nsw_d[f]* MDN[f]/365.0 * grid->area;
        h_cld[year] += fweight * grid->tcdc_clm[f]* MDN[f]/365.0 * grid->area;
        
        h_apar[year] += fweight * loct->appfd_g[f]* MDN[f]/365.0 * grid->area;
        h_parb[year] += fweight * grid->par_bp[f]* MDN[f]/365.0 * grid->area;
        h_pard[year] += fweight * grid->par_dp[f]* MDN[f]/365.0 * grid->area;
        
		/* potential permafrost area: added by A.Ito (2010/03/27) */
		if(grid->tmp_sfc_am <= -2.0){
			h_pot_prmfrst[year] += fweight * grid->area;
		}
		
		if(DF97==1){
			h_gpp[year] += fweight * (flux->plant).gpp_df97[f] * grid->area;
			h_gpp_c4[year] += fweight * (flux->c4).gpp_df97[f] * grid->area*loct->c4ptn[f];
		}else{
			h_gpp[year] += fweight * (flux->plant).gpp[f] * grid->area;
			h_gpp_c4[year] += fweight * (flux->c4).gpp[f] * grid->area*loct->c4ptn[f];
		}
		h_npp[year] += fweight * (flux->plant).npp[f] * grid->area;
		h_nep[year] += fweight * flux->nep[f] * grid->area;
		h_nbp[year] += fweight * flux->nbp[f] * grid->area;
		h_hvst[year] += fweight * (flux->plant).hvst[f] * grid->area;
		
		h_plant[year] += fweight * ((mass->plant).mfol[f] + (mass->plant).mstm[f] + 
									(mass->plant).mrot[f]) * MDN[f] /365.0 * grid->area;
		h_soil[year] += fweight * ((mass->soil).ltr_m[f] + 
								   (mass->soil).msl_m[f])* MDN[f]/365.0 * grid->area;
		
        h_arm[year] += fweight * (flux->plant).arm[f] * grid->area;
		
		/* added by A.Ito (2011/12/16) */
		h_abgm[year] += fweight * ((mass->plant).mfol[f] + (mass->plant).mstm[f]) * MDN[f] /365.0 * grid->area;
		
		h_sr[year] += fweight * ((flux->plant).rrm[f] + (flux->plant).rrg[f] + (flux->soil).hr[f]) * grid->area;
		h_doc[year] += fweight * (flux->soil).doc_boyer[f] * grid->area;
		
		h_gpp_df97[year] += fweight * (flux->plant).gpp_df97[f] * grid->area;

		/* GHG */
		h_ch4ox1[year] += fweight * (flux->soil).ch4oxy_ridg[f] * grid->area *10000.0/1000.0;
		h_ch4ox2[year] += fweight * (flux->soil).ch4oxy_casa[f] * grid->area *10000.0/1000.0;
		h_ch4ox3[year] += fweight * (flux->soil).ch4oxy_delgrosso[f] * grid->area *10000.0/1000.0;
		h_ch4ox4[year] += fweight * (flux->soil).ch4oxy_curry[f] * grid->area *10000.0/1000.0;		
		
		h_ch4emit_cao_paddy[year] += fweight_pad * (flux->soil).ch4flux_paddy_cao[f] * grid->area *10000.0/1000.0;
		h_ch4emit_cao_wetland[year] += fweight_wet * (flux->soil).ch4flux_wetland_cao[f] * grid->area *10000.0/1000.0;

		h_ch4emit_paddy_wh_diff[year] += fweight_pad * (flux->soil).ch4_paddy_wh_diff[f] * grid->area *10000.0/1000.0;
		h_ch4emit_paddy_wh_plant[year] += fweight_pad * (flux->soil).ch4_paddy_wh_plant[f] * grid->area *10000.0/1000.0;
		h_ch4emit_paddy_wh_ebbl[year] += fweight_pad * (flux->soil).ch4_paddy_wh_ebull[f] * grid->area *10000.0/1000.0;
		h_ch4emit_paddy_wh_release[year] += fweight_pad * (flux->soil).ch4_paddy_wh_release[f] * grid->area *10000.0/1000.0;

		h_ch4emit_wetland_wh_diff[year] += fweight_wet * (flux->soil).ch4_wetland_wh_diff[f] * grid->area *10000.0/1000.0;
		h_ch4emit_wetland_wh_plant[year] += fweight_wet * (flux->soil).ch4_wetland_wh_plant[f] * grid->area *10000.0/1000.0;
		h_ch4emit_wetland_wh_ebbl[year] += fweight_wet * (flux->soil).ch4_wetland_wh_ebull[f] * grid->area *10000.0/1000.0;
		h_ch4emit_wetland_wh_release[year] += fweight_wet * (flux->soil).ch4_wetland_wh_release[f] * grid->area *10000.0/1000.0;
		
		hm_temp[year][f] += fweight * grid->tmp_2m[f] * grid->area;
		hm_prec[year][f] += fweight * grid->prate_sfc[f] * grid->area;
		hm_ch4_wh[year][f] += fweight_wet * ((flux->soil).ch4_wetland_wh_diff[f] + (flux->soil).ch4_wetland_wh_plant[f] + 
							(flux->soil).ch4_wetland_wh_ebull[f] + (flux->soil).ch4_wetland_wh_release[f]) * 
							grid->area *10000.0/1000.0;
		
		hm_inund[year][f] += fweight_wet * loct->f_inund_wet_wh[f] * grid->area;
		
		h_n2o_emit_ngas[year] += fweight * (flux->soil).d_n2o_ngas[f] * grid->area;
		h_n2_emit_ngas[year] += fweight * (flux->soil).d_n2_ngas[f] * grid->area;
		h_n2o_emit_casa[year] += fweight * (flux->soil).d_n2o_casa[f] * grid->area;
		h_n2_emit_casa[year] += fweight * (flux->soil).d_n2_casa[f] * grid->area;
		h_no_emit_casa[year] += fweight * (flux->soil).d_no_casa[f] * grid->area;
		h_nh3_emit[year] += fweight * (flux->soil).n_nh3vlt[f] * grid->area;		
		h_ch4_emit_mass[year] += fweight * (flux->plant).emit_ch4_kirschbaum_mass[f] * 1000.0 * grid->area;
		h_ch4_emit_photo[year] += fweight * (flux->plant).emit_ch4_kirschbaum_photo[f] * 1000.0 * grid->area;
		h_n2_biofix[year] += fweight * (flux->plant).n_biofix[f] * grid->area;
		h_n2o_n_emit_ngas[year] += fweight * (flux->soil).d_n2o_ntr_ngas[f] * grid->area;
		h_n2o_d_emit_ngas[year] += fweight * (flux->soil).d_n2o_dnt_ngas[f] * grid->area;
		h_no3_leach[year] += fweight * (flux->soil).n_leach[f] * grid->area;
		
		/* added by A.Ito (2010/05/02) */
		h_n_fertin[year] += fweight * (flux->soil).n_fertin[f] * grid->area;
		/* corrected by A.Ito (2013/11/07) */
		h_n_depoin[year] += fweight * (loct->depo_nh4[f] + loct->depo_no3[f]) * grid->area;
		
		if(loct->v_type == 1 && REPL_OLSON_CROP == 0){ /* added by A.Ito (2009/06/16) */
			if(grid->veg_olson==29 || grid->veg_olson==30 || grid->veg_olson==31 || grid->veg_olson==32){
				h_n2o_emit_ngas_agr[year] += (flux->soil).d_n2o_ngas[f] * grid->area;
				h_n2o_emit_casa_agr[year] += (flux->soil).d_n2o_casa[f] * grid->area;
				h_nh3_emit_agr[year] += (flux->soil).n_nh3vlt[f] * grid->area;
			}
		}
		if(loct->v_type == 2){
			h_n2o_emit_ngas_agr[year] += fweight * (flux->soil).d_n2o_ngas[f] * grid->area;
			h_n2o_emit_casa_agr[year] += fweight * (flux->soil).d_n2o_casa[f] * grid->area;
			h_nh3_emit_agr[year] += fweight * (flux->soil).n_nh3vlt[f] * grid->area;
		}
		
		/* biomass burninig *****************/
		h_burnt_area[year] += fweight * flux->a_burnt[f]*grid->area;
		h_bioburn_co2[year] += fweight * (flux->bb_co2_litter[f]+flux->bb_co2_leaf[f]+
										  flux->bb_co2_wood[f]+flux->bb_co2_root[f]) * grid->area;
		h_bioburn_co[year] += fweight * (flux->bb_co_litter[f]+flux->bb_co_leaf[f]+flux->bb_co_wood[f]+
										 flux->bb_co_root[f]) * grid->area;
		h_bioburn_ch4[year] += fweight * (flux->bb_ch4_litter[f]+flux->bb_ch4_leaf[f]+flux->bb_ch4_wood[f]+
										  flux->bb_ch4_root[f]) * grid->area;
		h_bioburn_nmhc[year] += fweight * (flux->bb_nmhc_litter[f]+flux->bb_nmhc_leaf[f]+
										   flux->bb_nmhc_wood[f]+flux->bb_nmhc_root[f]) * grid->area;
		h_bioburn_oc[year] += fweight * (flux->bb_oc_litter[f]+flux->bb_oc_leaf[f]+
										 flux->bb_oc_wood[f]+flux->bb_oc_root[f]) * grid->area;
		h_bioburn_bc[year] += fweight * (flux->bb_bc_litter[f]+flux->bb_bc_leaf[f]+
										 flux->bb_bc_wood[f]+flux->bb_bc_root[f]) * grid->area;
		h_bioburn_nox[year] += fweight * (flux->bb_nox_litter[f]+flux->bb_nox_leaf[f]+
										  flux->bb_nox_wood[f]+flux->bb_nox_root[f]) * grid->area;
		h_bioburn_so2[year] += fweight * (flux->bb_so2_litter[f]+flux->bb_so2_leaf[f]+
										  flux->bb_so2_wood[f]+flux->bb_so2_root[f]) * grid->area;
		h_bioburn_pm25[year] += fweight * (flux->bb_pm25_litter[f]+flux->bb_pm25_leaf[f]+
										   flux->bb_pm25_wood[f]+flux->bb_pm25_root[f]) * grid->area;
		h_bioburn_tpm[year] += fweight * (flux->bb_tpm_litter[f]+flux->bb_tpm_leaf[f]+
										  flux->bb_tpm_wood[f]+flux->bb_tpm_root[f]) * grid->area;
		h_bioburn_tec[year] += fweight * (flux->bb_tec_litter[f]+flux->bb_tec_leaf[f]+
										  flux->bb_tec_wood[f]+flux->bb_tec_root[f]) * grid->area;

		/* VOC, g C ha-1 month-1 **************/
		h_voc_isopr_g97[year] += fweight * flux->voc_isopr_g97[f] * grid->area *10000.0/1000000.0;
		h_voc_monotrp_g97[year] += fweight * flux->voc_monotrp_g97[f] * grid->area *10000.0/1000000.0;
		h_voc_methanl_g97[year] += fweight * flux->voc_methanl_g97[f] * grid->area *10000.0/1000000.0;
		h_voc_acetone_g97[year] += fweight * flux->voc_acetone_g97[f] * grid->area *10000.0/1000000.0;
		h_voc_actaldhd_g97[year] += fweight * flux->voc_actaldhd_g97[f] * grid->area *10000.0/1000000.0;
		h_voc_frmardhd_g97[year] += fweight * flux->voc_frmardhd_g97[f] * grid->area *10000.0/1000000.0;
		h_voc_formacd_g97[year] += fweight * flux->voc_formacd_g97[f] * grid->area *10000.0/1000000.0;
		h_voc_acetacd_g97[year] += fweight * flux->voc_acetacd_g97[f] * grid->area *10000.0/1000000.0;
		h_voc_co_g97[year] += fweight * flux->voc_co_g97[f] * grid->area *10000.0/1000000.0;
        /* added 2014/09/ */
		h_voc_afarnesene[year] += fweight * flux->voc_afarnesene[f] * grid->area *10000.0/1000000.0;
		h_voc_bcaryophyllene[year] += fweight * flux->voc_bcaryophyllene[f] * grid->area *10000.0/1000000.0;
		h_voc_othersesqui[year] += fweight * flux->voc_othersesqui[f] * grid->area *10000.0/1000000.0;
		
		/* d13c & d14c : added by A.Ito (2009/07/15) ***************/
		ci_aco2_d13c[year] = d13c_addition(loct->d13c_aco2[f],loct->aco2[f]*grid->area 
										   * MDN[f]/365.0, ci_aco2_d13c[year],ci_aco2[year]);
		ci_aco2_d14c[year] = (grid->d14c_bco2[f]*loct->aco2[f]*grid->area * MDN[f]/365.0 + ci_aco2_d14c[year]*ci_aco2[year]) / 
							(loct->aco2[f]*grid->area * MDN[f]/365.0 + ci_aco2[year]);
		ci_aco2[year] += loct->aco2[f]*grid->area * MDN[f]/365.0;
		
		/* GPP */
		if(DF97==1){
			if((flux->plant).gpp_df97[f] > 0.0){
				ci_gpp_d13c[year] = d13c_addition((flux->plant).d13c_gpp[f], (flux->plant).gpp_df97[f] 
											* grid->area, ci_gpp_d13c[year], ci_gpp[year]);
				ci_gpp_d14c[year] = ((flux->plant).d14c_gpp[f] * (flux->plant).gpp_df97[f] * grid->area 
									+ ci_gpp_d14c[year] * ci_gpp[year]) / ((flux->plant).gpp_df97[f] 
									* grid->area + ci_gpp[year]);
				ci_gpp[year] += (flux->plant).gpp_df97[f] * grid->area;
			}
		}else{
			if((flux->plant).gpp[f] > 0.0){
				ci_gpp_d13c[year] = d13c_addition((flux->plant).d13c_gpp[f], (flux->plant).gpp[f] 
												  * grid->area, ci_gpp_d13c[year], ci_gpp[year]);
				ci_gpp_d14c[year] = ((flux->plant).d14c_gpp[f] * (flux->plant).gpp[f] 
									 * grid->area + ci_gpp_d14c[year] * ci_gpp[year]) / 
									((flux->plant).gpp[f] * grid->area + ci_gpp[year]);
				ci_gpp[year] += (flux->plant).gpp[f] * grid->area;
			}
		}
		/* ER */
		if(flux->er[f] > 0.0){
			ci_er_d13c[year] = d13c_addition(flux->d13c_er[f], flux->er[f] * grid->area, ci_er_d13c[year], ci_er[year]);
			ci_er_d14c[year] = (flux->d14c_er[f] * flux->er[f] * grid->area + ci_er_d14c[year] * ci_er[year]) / 
							(flux->er[f] * grid->area + ci_er[year]);
			ci_er[year] += flux->er[f] * grid->area;
		}
		
		/* MASS */
		if((mass->plant).mfol[f] > 0.0){
			ci_f_d13c[year] = d13c_addition((mass->plant).d13c_mfol[f], (mass->plant).mfol[f] * 
											grid->area * MDN[f]/365.0, ci_f_d13c[year], ci_f[year]);
			ci_f_d14c[year] = ((mass->plant).d14c_mfol[f]*(mass->plant).mfol[f] * grid->area * 
							   MDN[f] /365.0 + ci_f_d14c[year]*ci_f[year]) / 
								((mass->plant).mfol[f] * grid->area * MDN[f]/365.0 + ci_f[year]);
			ci_f[year] += (mass->plant).mfol[f] * grid->area * MDN[f]/365.0;
		}
		if((mass->plant).mstm[f] > 0.0){
			ci_c_d13c[year] = d13c_addition((mass->plant).d13c_mstm[f], (mass->plant).mstm[f] * 
											grid->area * MDN[f]/365.0, ci_c_d13c[year], ci_c[year]);
			ci_c_d14c[year] = ((mass->plant).d14c_mstm[f]*(mass->plant).mstm[f] * grid->area * 
							   MDN[f] /365.0 + ci_c_d14c[year]*ci_c[year]) / 
								((mass->plant).mstm[f] * grid->area * MDN[f]/365.0 + ci_c[year]);
			ci_c[year] += (mass->plant).mstm[f] * grid->area * MDN[f]/365.0;
		}
		if((mass->plant).mrot[f] > 0.0){
			ci_r_d13c[year] = d13c_addition((mass->plant).d13c_mrot[f], (mass->plant).mrot[f] * 
											grid->area * MDN[f]/365.0, ci_r_d13c[year], ci_r[year]);
			ci_r_d14c[year] = ((mass->plant).d14c_mrot[f]*(mass->plant).mrot[f] * grid->area * 
							   MDN[f] /365.0 + ci_r_d14c[year]*ci_r[year]) / 
								((mass->plant).mrot[f] * grid->area * MDN[f]/365.0 + ci_r[year]);
			ci_r[year] += (mass->plant).mrot[f] * grid->area * MDN[f]/365.0;
		}
		
		if((mass->soil).ltr_m[f] > 0.0){
			ci_l_d13c[year] = d13c_addition((mass->soil).d13c_ltr_m[f], (mass->soil).ltr_m[f] * 
											grid->area * MDN[f]/365.0, ci_l_d13c[year], ci_l[year]);
			ci_l_d14c[year] = ((mass->soil).d14c_ltr_m[f]*(mass->soil).ltr_m[f] * 
							   grid->area * MDN[f] /365.0 + ci_l_d14c[year]*ci_l[year]) / 
								((mass->soil).ltr_m[f] * grid->area * MDN[f]/365.0 + ci_l[year]);
			ci_l[year] += (mass->soil).ltr_m[f] * grid->area * MDN[f]/365.0;
		}
		if((mass->soil).msl_m[f] > 0.0){
			ci_h_d13c[year] = d13c_addition((mass->soil).d13c_msl_m[f], (mass->soil).msl_m[f] * 
											grid->area * MDN[f]/365.0, ci_h_d13c[year], ci_h[year]);
			ci_h_d14c[year] = ((mass->soil).d14c_msl_m[f]*(mass->soil).msl_m[f] * 
							   grid->area * MDN[f] /365.0 + ci_h_d14c[year]*ci_h[year]) / 
								((mass->soil).msl_m[f] * grid->area * MDN[f]/365.0 + ci_h[year]);
			ci_h[year] += (mass->soil).msl_m[f] * grid->area * MDN[f]/365.0;
		}
		
		/* regional: added by A.Ito (2009/07/12) ****************/
		rh_temp[grid->reg_g][year] += grid->tmp_2m[f]* MDN[f]/365.0 * grid->area;
		rh_prec[grid->reg_g][year] += grid->prate_sfc[f] * grid->area;
		rh_dswrf[grid->reg_g][year] += grid->gl_rad[f]* MDN[f]/365.0 * grid->area;
		rh_rns[grid->reg_g][year] += loct->rad_net_long[f] * MDN[f]/365.0 * grid->area;
		rh_rnl[grid->reg_g][year] += loct->rad_net_short[f] * MDN[f]/365.0 * grid->area;
		rh_ipar[grid->reg_g][year] += grid->dlen[f] * 0.5 * grid->par[f] * MDN[f]/365.0 * grid->area;
		rh_apar[grid->reg_g][year] += loct->appfd_g[f] * MDN[f]/365.0 * grid->area;
		
		if(DF97==1){
			rh_gpp[grid->reg_g][year] += (flux->plant).gpp_df97[f] * grid->area;
		}else{
			rh_gpp[grid->reg_g][year] += (flux->plant).gpp[f] * grid->area;
		}
		rh_npp[grid->reg_g][year] += (flux->plant).npp[f] * grid->area;
		rh_nep[grid->reg_g][year] += flux->nep[f] * grid->area;
		rh_evpr[grid->reg_g][year] += loct->evpr[f] * grid->area;
		rh_trsp[grid->reg_g][year] += loct->trspr[f] * grid->area;
		rh_incp[grid->reg_g][year] += loct->incep[f] * grid->area;
		rh_rnof[grid->reg_g][year] += loct->ro2[f] * grid->area;
		
		/* added by A.Ito (2009/11/15) */
		rh_hvst[grid->reg_g][year] += (flux->plant).hvst[f] * grid->area;
		rh_ch4ox_curry[grid->reg_g][year] += (flux->soil).ch4oxy_curry[f] * grid->area *10000.0/1000.0;
		rh_ch4emit_wh_wet[grid->reg_g][year] += ((flux->soil).ch4_wetland_wh_diff[f] + (flux->soil).ch4_wetland_wh_plant[f] + 
												 (flux->soil).ch4_wetland_wh_ebull[f] + (flux->soil).ch4_wetland_wh_release[f]) 
												* grid->area *10000.0/1000.0;
		rh_ch4emit_wh_paddy[grid->reg_g][year] += ((flux->soil).ch4_paddy_wh_diff[f] + (flux->soil).ch4_paddy_wh_plant[f] + 
													 (flux->soil).ch4_paddy_wh_ebull[f] + (flux->soil).ch4_paddy_wh_release[f]) 
														* grid->area *10000.0/1000.0;
		rh_n2o_emit_ngas[grid->reg_g][year] += ((flux->soil).d_n2o_ntr_ngas[f] + (flux->soil).d_n2o_dnt_ngas[f]) * grid->area;
		if(CALC_OLSON == 1){
			if(grid->veg_olson==29 || grid->veg_olson==30 || grid->veg_olson==31 || grid->veg_olson==32){
				rh_n2o_emitagr_ngas[grid->reg_g][year] += ((flux->soil).d_n2o_ntr_ngas[f] + (flux->soil).d_n2o_dnt_ngas[f]) * grid->area;
			}
		}
		
		/* isotopes, added by A.Ito (2009/09/30) */
		/* GPP */
		if(DF97==1){
			if((flux->plant).gpp_df97[f] > 0.0){
				rh_ci_gpp_d13c[grid->reg_g][year] = d13c_addition((flux->plant).d13c_gpp[f], (flux->plant).gpp_df97[f] * grid->area, 
															rh_ci_gpp_d13c[grid->reg_g][year], rh_ci_gpp[grid->reg_g][year]);
				rh_ci_gpp_d14c[grid->reg_g][year] = ((flux->plant).d14c_gpp[f] * (flux->plant).gpp_df97[f] * grid->area 
													 + rh_ci_gpp_d14c[grid->reg_g][year] * rh_ci_gpp[grid->reg_g][year]) / 
													((flux->plant).gpp_df97[f] * grid->area + rh_ci_gpp[grid->reg_g][year]);
				rh_ci_gpp[grid->reg_g][year] += (flux->plant).gpp_df97[f] * grid->area;
			}
		}else{
			if((flux->plant).gpp[f] > 0.0){
				rh_ci_gpp_d13c[grid->reg_g][year] = d13c_addition((flux->plant).d13c_gpp[f], (flux->plant).gpp[f] * grid->area, 
															   rh_ci_gpp_d13c[grid->reg_g][year], rh_ci_gpp[grid->reg_g][year]);
				rh_ci_gpp_d14c[grid->reg_g][year] = ((flux->plant).d14c_gpp[f] * (flux->plant).gpp[f] * grid->area 
													 + rh_ci_gpp_d14c[grid->reg_g][year] * rh_ci_gpp[grid->reg_g][year]) / 
													((flux->plant).gpp[f] * grid->area + rh_ci_gpp[grid->reg_g][year]);
				rh_ci_gpp[grid->reg_g][year] += (flux->plant).gpp[f] * grid->area;
			}
		}
		/* ER */
		if(flux->er[f] > 0.0){
			rh_ci_er_d13c[grid->reg_g][year] = d13c_addition(flux->d13c_er[f], flux->er[f] * grid->area, 
															 rh_ci_er_d13c[grid->reg_g][year], rh_ci_er[grid->reg_g][year]);
			rh_ci_er_d14c[grid->reg_g][year] = (flux->d14c_er[f] * flux->er[f] * grid->area 
												+ rh_ci_er_d14c[grid->reg_g][year] * rh_ci_er[grid->reg_g][year]) / 
												(flux->er[f] * grid->area + rh_ci_er[grid->reg_g][year]);
			rh_ci_er[grid->reg_g][year] += flux->er[f] * grid->area;
		}
		
		/* MASS */
		if((mass->plant).mfol[f] > 0.0){
			rh_ci_f_d13c[grid->reg_g][year] = d13c_addition((mass->plant).d13c_mfol[f], (mass->plant).mfol[f] * 
								grid->area * MDN[f]/365.0, rh_ci_f_d13c[grid->reg_g][year], rh_ci_f[grid->reg_g][year]);
			rh_ci_f_d14c[grid->reg_g][year] = ((mass->plant).d14c_mfol[f]*(mass->plant).mfol[f] * grid->area * 
							   MDN[f] /365.0 + rh_ci_f_d14c[grid->reg_g][year]*rh_ci_f[grid->reg_g][year]) / 
								((mass->plant).mfol[f] * grid->area * MDN[f]/365.0 + rh_ci_f[grid->reg_g][year]);
			rh_ci_f[grid->reg_g][year] += (mass->plant).mfol[f] * grid->area * MDN[f]/365.0;
		}
		if((mass->plant).mstm[f] > 0.0){
			rh_ci_c_d13c[grid->reg_g][year] = d13c_addition((mass->plant).d13c_mstm[f], (mass->plant).mstm[f] * 
								grid->area * MDN[f]/365.0, rh_ci_c_d13c[grid->reg_g][year], rh_ci_c[grid->reg_g][year]);
			rh_ci_c_d14c[grid->reg_g][year] = ((mass->plant).d14c_mstm[f]*(mass->plant).mstm[f] * grid->area * 
							   MDN[f] /365.0 + rh_ci_c_d14c[grid->reg_g][year]*rh_ci_c[grid->reg_g][year]) / 
								((mass->plant).mstm[f] * grid->area * MDN[f]/365.0 + rh_ci_c[grid->reg_g][year]);
			rh_ci_c[grid->reg_g][year] += (mass->plant).mstm[f] * grid->area * MDN[f]/365.0;
		}
		if((mass->plant).mrot[f] > 0.0){
			rh_ci_r_d13c[grid->reg_g][year] = d13c_addition((mass->plant).d13c_mrot[f], (mass->plant).mrot[f] * 
                                grid->area * MDN[f]/365.0, rh_ci_r_d13c[grid->reg_g][year], rh_ci_r[grid->reg_g][year]);
			rh_ci_r_d14c[grid->reg_g][year] = ((mass->plant).d14c_mrot[f]*(mass->plant).mrot[f] * grid->area * 
							   MDN[f] /365.0 + rh_ci_r_d14c[grid->reg_g][year]*rh_ci_r[grid->reg_g][year]) / 
								((mass->plant).mrot[f] * grid->area * MDN[f]/365.0 + rh_ci_r[grid->reg_g][year]);
			rh_ci_r[grid->reg_g][year] += (mass->plant).mrot[f] * grid->area * MDN[f]/365.0;
		}
		
		if((mass->soil).ltr_m[f] > 0.0){
			rh_ci_l_d13c[grid->reg_g][year] = d13c_addition((mass->soil).d13c_ltr_m[f], (mass->soil).ltr_m[f] * 
								grid->area * MDN[f]/365.0, rh_ci_l_d13c[grid->reg_g][year], rh_ci_l[grid->reg_g][year]);
			rh_ci_l_d14c[grid->reg_g][year] = ((mass->soil).d14c_ltr_m[f]*(mass->soil).ltr_m[f] * 
							   grid->area * MDN[f] /365.0 + rh_ci_l_d14c[grid->reg_g][year]*rh_ci_l[grid->reg_g][year]) / 
								((mass->soil).ltr_m[f] * grid->area * MDN[f]/365.0 + rh_ci_l[grid->reg_g][year]);
			rh_ci_l[grid->reg_g][year] += (mass->soil).ltr_m[f] * grid->area * MDN[f]/365.0;
		}
		if((mass->soil).msl_m[f] > 0.0){
			rh_ci_h_d13c[grid->reg_g][year] = d13c_addition((mass->soil).d13c_msl_m[f], (mass->soil).msl_m[f] * 
								grid->area * MDN[f]/365.0, rh_ci_h_d13c[grid->reg_g][year], rh_ci_h[grid->reg_g][year]);
			rh_ci_h_d14c[grid->reg_g][year] = ((mass->soil).d14c_msl_m[f]*(mass->soil).msl_m[f] * 
							   grid->area * MDN[f] /365.0 + rh_ci_h_d14c[grid->reg_g][year]*rh_ci_h[grid->reg_g][year]) / 
								((mass->soil).msl_m[f] * grid->area * MDN[f]/365.0 + rh_ci_h[grid->reg_g][year]);
			rh_ci_h[grid->reg_g][year] += (mass->soil).msl_m[f] * grid->area * MDN[f]/365.0;
		}
	}
	
	/* erosion and land-use change: 2011/12/15 (A.Ito) */
	if(loct->v_type == 1){
		h_ersn_c[year] += fweight * flux->erod_carbon * grid->area;

		h_luc[year] += (flux->lu_conv + flux->lu_ten + flux->lu_hund) * grid->area;
		h_luc_1[year] += flux->lu_conv * grid->area;
		h_luc_2[year] += flux->lu_ten * grid->area;
		h_luc_3[year] += flux->lu_hund * grid->area;
		
		rh_luc[grid->reg_g][year] += (flux->lu_conv + flux->lu_ten + flux->lu_hund) * grid->area;
		
        /* 2014/12/10 by A.Ito */
		/* h_hvst_wood[year] += flux->hvst_wood * grid->area; */
		h_hvst_wood[year] += fweight * flux->hvst_wood * grid->area;
        
        h_wetarea[year] += grid->f_wetland * grid->area;

        h_deforest[year] += grid->area * grid->f_deforest;
	}
	if(loct->v_type == 2){
		h_agrersn_c[year] += fweight * grid->area * flux->erod_carbon;
		
		h_agrarea[year] += fweight * grid->area;
		h_paddyarea[year] += grid->f_paddy * grid->area;
	}
}

/*********************************************************************/
void f_glosum_output(
	char s_date[32], 
	char s_case[32]
){
	long h, i;
	char filename[128];
	FILE *fp_glsum;

	/* global summary ******************************************/
	strcpy(filename, s_date);
	strcat(filename, s_case);
	strcat(filename, "glsum.txt");
	
	fp_glsum = fopen(filename,"wt");
	fprintf(fp_glsum,"%.2lf %.2lf\n", go_landarea, gs_landarea);
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
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
		fprintf(fp_glsum,"%lf ", h_ch4ox4[h]);		/* 080704 */
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
		
		fprintf(fp_glsum,"%lf ", h_luc_1[h]);
		fprintf(fp_glsum,"%lf ", h_luc_2[h]);
		fprintf(fp_glsum,"%lf ", h_luc_3[h]);

		fprintf(fp_glsum,"%lf ", h_gpp_df97[h]);
		fprintf(fp_glsum,"%lf ", h_n2o_n_emit_ngas[h]);
		fprintf(fp_glsum,"%lf ", h_n2o_d_emit_ngas[h]);
		
		fprintf(fp_glsum,"%lf ", h_n2o_emit_ngas_agr[h]); /* added by A.Ito (2009/06/16) */
		fprintf(fp_glsum,"%lf ", h_n2o_emit_casa_agr[h]);
		fprintf(fp_glsum,"%lf ", h_nh3_emit_agr[h]);
		
		/* added by A.Ito (2009/08/31) */
		fprintf(fp_glsum,"%lf ", h_ch4emit_paddy_wh_diff[h]); 
		fprintf(fp_glsum,"%lf ", h_ch4emit_paddy_wh_plant[h]);
		fprintf(fp_glsum,"%lf ", h_ch4emit_paddy_wh_ebbl[h]);
		fprintf(fp_glsum,"%lf ", h_ch4emit_paddy_wh_release[h]);
		fprintf(fp_glsum,"%lf ", h_ch4emit_wetland_wh_diff[h]);
		fprintf(fp_glsum,"%lf ", h_ch4emit_wetland_wh_plant[h]);
		fprintf(fp_glsum,"%lf ", h_ch4emit_wetland_wh_ebbl[h]);
		fprintf(fp_glsum,"%lf ", h_ch4emit_wetland_wh_release[h]);

		fprintf(fp_glsum,"%lf ", h_gpp_c4[h]); /* added by A.Ito (2009/08/31) */
		fprintf(fp_glsum,"%lf ", h_pot_prmfrst[h]); /* added by A.Ito (2010/03/27) */
		fprintf(fp_glsum,"%lf ", h_no3_leach[h]); /* added by A.Ito (2010/03/29) */

		fprintf(fp_glsum,"%lf ", h_n_fertin[h]); /* added by A.Ito (2010/05/02) */
		fprintf(fp_glsum,"%lf ", h_n_depoin[h]); /* added by A.Ito (2010/05/02) */
		
		fprintf(fp_glsum,"%lf ", h_hvst_wood[h]); /* added by A.Ito (2010/11/09) */
		
		fprintf(fp_glsum,"%lf ", h_trnsp[h]);  /* added by A.Ito (2010/11/18) */
		fprintf(fp_glsum,"%lf ", h_incepev[h]); 
		fprintf(fp_glsum,"%lf ", h_ssurfev[h]); 

		fprintf(fp_glsum,"%lf ", h_nbp[h]);   /* added by A.Ito (2010/11/27) */
		fprintf(fp_glsum,"%lf ", h_hvst[h]);
		fprintf(fp_glsum,"%lf ", h_paddyarea[h]); /* added by A.Ito (2011/2/28) */
		fprintf(fp_glsum,"%lf ", h_abgm[h]); /* added by A.Ito (2011/12/16) */

		fprintf(fp_glsum,"%lf ", h_sw1[h]); /* added by A.Ito (2012/01/05) */
		fprintf(fp_glsum,"%lf ", h_sw2[h]); /* added by A.Ito (2011/01/05) */

		fprintf(fp_glsum,"%lf ", h_wetarea[h]); /* added by A.Ito (2012/10/26) */

		fprintf(fp_glsum,"%lf ", h_rns[h]); /* added by A.Ito (2013/01/02) */
		fprintf(fp_glsum,"%lf ", h_rnl[h]);
		fprintf(fp_glsum,"%lf ", h_rnsd[h]);
		fprintf(fp_glsum,"%lf ", h_cld[h]);

		fprintf(fp_glsum,"%lf ", h_apar[h]); /* added by A.Ito (2013/12/16) */
		fprintf(fp_glsum,"%lf ", h_parb[h]); /* added by A.Ito (2013/12/24) */
		fprintf(fp_glsum,"%lf ", h_pard[h]); /* added by A.Ito (2013/12/24) */

		fprintf(fp_glsum,"%lf ", h_arm[h]); /* added by A.Ito (2014/02/14) */

		fprintf(fp_glsum,"%lf ", h_voc_afarnesene[h]); /* added 2014/9/11 by A.Ito */
		fprintf(fp_glsum,"%lf ", h_voc_bcaryophyllene[h]);
		fprintf(fp_glsum,"%lf ", h_voc_othersesqui[h]);

		fprintf(fp_glsum,"%lf ", h_deforest[h]); /* added by A.Ito (2014/09/22) */

		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<ASTEP;h++){
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
		
		fprintf(fp_glsum,"%lf ", m_gpp[h]);
		fprintf(fp_glsum,"%lf ", m_npp[h]);
		fprintf(fp_glsum,"%lf ", m_nep[h]);

		fprintf(fp_glsum,"%lf ", m_ch4p_cao[h]);
		fprintf(fp_glsum,"%lf ", m_ch4p_wh[h]);
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<NVEG_OLSON;h++){
		fprintf(fp_glsum,"%lf ", vo_area[h]);
		fprintf(fp_glsum,"%lf ", vo_gpp[h]);
		fprintf(fp_glsum,"%lf ", vo_npp[h]);
		fprintf(fp_glsum,"%lf ", vo_nep[h]);
		fprintf(fp_glsum,"%lf ", vo_lai[h]);
		fprintf(fp_glsum,"%lf ", vo_fol[h]);
		fprintf(fp_glsum,"%lf ", vo_stm[h]);
		fprintf(fp_glsum,"%lf ", vo_rot[h]);
		fprintf(fp_glsum,"%lf ", vo_ltr[h]);
		fprintf(fp_glsum,"%lf ", vo_msl[h]);
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<NVEG_SAGE;h++){
		fprintf(fp_glsum,"%lf ", vs_area[h]);
		fprintf(fp_glsum,"%lf ", vs_gpp[h]);
		fprintf(fp_glsum,"%lf ", vs_npp[h]);
		fprintf(fp_glsum,"%lf ", vs_nep[h]);
		fprintf(fp_glsum,"%lf ", vs_lai[h]);
		fprintf(fp_glsum,"%lf ", vs_fol[h]);
		fprintf(fp_glsum,"%lf ", vs_stm[h]);
		fprintf(fp_glsum,"%lf ", vs_rot[h]);
		fprintf(fp_glsum,"%lf ", vs_ltr[h]);
		fprintf(fp_glsum,"%lf ", vs_msl[h]);
		fprintf(fp_glsum,"\n");
	}
	
	fprintf(fp_glsum,"YEAR ");
	for(i=0;i<N_REG;i++){
		fprintf(fp_glsum,"%lf ", rh_area[i]);
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_temp[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_prec[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_dswrf[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_rns[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_rnl[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ipar[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_apar[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_gpp[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_npp[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_nep[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_evpr[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_trsp[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_incp[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_rnof[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
		
	/* historical carbon isotopes: added by A.Ito (2009/07/15) */
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		fprintf(fp_glsum,"%lf ", ci_aco2[h]);
		fprintf(fp_glsum,"%lf ", ci_aco2_d13c[h]);
		fprintf(fp_glsum,"%lf ", ci_aco2_d14c[h]);

		fprintf(fp_glsum,"%lf ", ci_gpp[h]);
		fprintf(fp_glsum,"%lf ", ci_gpp_d13c[h]);
		fprintf(fp_glsum,"%lf ", ci_gpp_d14c[h]);
		fprintf(fp_glsum,"%lf ", ci_er[h]);
		fprintf(fp_glsum,"%lf ", ci_er_d13c[h]);
		fprintf(fp_glsum,"%lf ", ci_er_d14c[h]);
		
		fprintf(fp_glsum,"%lf ", ci_f[h]);
		fprintf(fp_glsum,"%lf ", ci_f_d13c[h]);
		fprintf(fp_glsum,"%lf ", ci_f_d14c[h]);
		fprintf(fp_glsum,"%lf ", ci_c[h]);
		fprintf(fp_glsum,"%lf ", ci_c_d13c[h]);
		fprintf(fp_glsum,"%lf ", ci_c_d14c[h]);
		fprintf(fp_glsum,"%lf ", ci_r[h]);
		fprintf(fp_glsum,"%lf ", ci_r_d13c[h]);
		fprintf(fp_glsum,"%lf ", ci_r_d14c[h]);
		fprintf(fp_glsum,"%lf ", ci_l[h]);
		fprintf(fp_glsum,"%lf ", ci_l_d13c[h]);
		fprintf(fp_glsum,"%lf ", ci_l_d14c[h]);
		fprintf(fp_glsum,"%lf ", ci_h[h]);
		fprintf(fp_glsum,"%lf ", ci_h_d13c[h]);
		fprintf(fp_glsum,"%lf ", ci_h_d14c[h]);
		
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");

	/* regional *************************************************/
	/* added by A.Ito (2009/09/30) */
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_gpp[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_gpp_d13c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_gpp_d14c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_er[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_er_d13c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_er_d14c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_f[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_f_d13c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_f_d14c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_c_d13c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_c_d14c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_r[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_r_d13c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_r_d14c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_l[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_l_d13c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_l_d14c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_h[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_h_d13c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_h_d14c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	
	/* added by A.Ito (2009/11/15) */
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_hvst[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_luc[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ch4ox_curry[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ch4emit_wh_wet[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_ch4emit_wh_paddy[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_n2o_emit_ngas[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");for(h=0;h<PD_SIM;h++){
		fprintf(fp_glsum,"%ld ", h+(PIVOT_CLIMY-1));
		for(i=0;i<N_REG;i++){
			fprintf(fp_glsum,"%lf ", rh_n2o_emitagr_ngas[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	
	for(h=0;h<PD_SIM;h++){
		for(i=0;i<ASTEP;i++){
			fprintf(fp_glsum,"%ld %ld %lf %lf %lf %lf\n", h+(PIVOT_CLIMY-1), i+1, 
				hm_temp[h][i], hm_prec[h][i], hm_ch4_wh[h][i], hm_inund[h][i]);
		}
	}
	fprintf(fp_glsum,"\n");
	
	/* perturbation factor */
	for(h=0;h<20;h++){
		fprintf(fp_glsum,"%ld %lf\n", h, f_pert[h]);
	}
    fprintf(fp_glsum,"\n");
    
    /* latitudinal */
    for(h=0;h<N_ROW;h++){
        fprintf(fp_glsum,"%lf ", glat_area[h]);
        
        for(i=0;i<ASTEP;i++){
            fprintf(fp_glsum,"%lf ", glat_gpp[i][h]);
        }
        for(i=0;i<ASTEP;i++){
            fprintf(fp_glsum,"%lf ", glat_npp[i][h]);
        }
        for(i=0;i<ASTEP;i++){
            fprintf(fp_glsum,"%lf ", glat_nep[i][h]);
        }
        for(i=0;i<ASTEP;i++){
            fprintf(fp_glsum,"%lf ", glat_ch4_cao[i][h]);
        }
        for(i=0;i<ASTEP;i++){
            fprintf(fp_glsum,"%lf ", glat_ch4_wh[i][h]);
        }
        
        fprintf(fp_glsum,"\n");
    }
	
	fclose(fp_glsum);
}
