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

/* clear all variables and parameters *******************************/
void f_clear(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux
){
	long f;
	
	grid->y = grid->m = 0;
	loct->time = loct->time_hyd = 0;
	
	vanish(mass, flux);
	
	(echar->c3).grw_pd = (echar->c4).grw_pd = 0.0;
	(echar->c3).gdd = (echar->c4).gdd = 0.0;
	
	grid->proj_prec_co = 0.0;	/* carry-over of negative precipitation change */
	for(f = 0;f<ASTEP;f++){
		grid->dlen[f] = grid->par[f] = 0.0;
		grid->par_bp[f] = grid->par_dp[f] = grid->par_be[f] = grid->par_de[f] = 0.0;
		grid->sl_dec[f] = grid->sl_hgt[f] = 0.0;
		grid->top_rad[f] = grid->gl_rad[f] = 0.0;
		
		loct->gl_rad_g[f] = 0.0;
		loct->rad_net_g[f] = loct->rad_net_p[f] = 0.0;
        loct->appfd_g[f] = loct->fappfd_g[f] = 0.0;
		loct->rad_net_long[f] = loct->albedo_sfc[f] = 0.0;
		loct->gd[f] = 0;
		loct->gdd[f] = 0.0;
		loct->prsr[f] = loct->dnsa[f] = loct->dnsa[f] = 0.0;
		loct->vp[f] = loct->vps[f] = loct->vpd[f] = 0.0;
		loct->slope_vps[f] = loct->r_aero[f] = 0.0;
		loct->pm_evp[f] = loct->pm_trn[f] = loct->pm_incep[f] = 0.0;
        loct->evpr[f] = loct->trspr[f] = loct->incep[f] = 0.0;
		loct->ro1[f] = loct->ro2[f] = loct->msnwa[f] = 0.0;
		loct->msww[f] = loct->msw30[f] = loct->snp[f] = 0.0;
		loct->snp[f] = loct->thaw[f] = 0.0;
		loct->vmc30[f] = loct->vmc[f] = 0.0;
		loct->npp_av[f] = 0.0;
		
		(echar->c3).mgdd[f] = (echar->c4).mgdd[f] = 0.0;
		(echar->c3).season[f] = (echar->c4).season[f] = 0;
		(echar->c3).psat[f] = (echar->c4).psat[f] = 0.0;
		(echar->c3).eK[f] = (echar->c4).eK[f] = 0.0;
		(echar->c3).lue[f] = (echar->c4).lue[f] = 0.0;
		(echar->c3).gs[f] = (echar->c4).gs[f] = 0.0;
		(echar->c3).gc[f] = (echar->c4).gc[f] = 0.0;
		(echar->c3).ci[f] = (echar->c4).ci[f] = 0.0;
		(echar->c3).opt_lai[f] = (echar->c4).opt_lai[f] = 0.0;
		(echar->c3).qTf[f] = (echar->c4).qTf[f] = 0.0;
		(echar->c3).qTc[f] = (echar->c4).qTc[f] = 0.0;
		(echar->c3).qTr[f] = (echar->c4).qTr[f] = 0.0;
		(echar->c3).lf[f] = (echar->c4).lf[f] = 0.0;
		(echar->c3).lc[f] = (echar->c4).lc[f] = 0.0;
		(echar->c3).lr[f] = (echar->c4).lr[f] = 0.0;
		(echar->c3).malloc_f[f] = (echar->c4).malloc_f[f] = 0.0;
		(echar->c3).malloc_c[f] = (echar->c4).malloc_c[f] = 0.0;
		(echar->c3).malloc_r[f] = (echar->c4).malloc_r[f] = 0.0;
		
		(echar->soil).albedo[f] = 0.0;
	}
	for(f=0;f<49;f++){
		(echar->c3).fleaf_age[f] = (echar->c4).fleaf_age[f] = 0.0;
	}
	grid->f_crop_con = grid->f_crop_p = grid->f_crop_trend = 0.0;
	grid->f_pasture_con = grid->f_pasture_p = grid->f_pasture_trend = 0.0;
	grid->f_deforest = 0.0;
	
	grid->f_erosion_r = 0.0;
	grid->f_erosion_ls = 0.0;
	grid->f_erosion_k = 0.0;
	grid->f_erosion_c = 0.0;
	grid->f_erosion_p = 0.0;
	
	for(f = 0;f<(N_SLAYER+2);f++){
		loct->prof_ch4[f] = 0.0;
	}
    
    grid->tmp_base_permaforst = 0.0;
    grid->tmp_soil_am = 0.0;
    
    loct->n_frtlz_in = loct->n_manure_in = 0.0;
    
    /* for(f=0;f<DL_BF;f++){
        grid->f_biofuel[f] = 0.0;
    } */
}

/* make plant fluxes vacant ***************************/
void plant_flux_zero(
	long month, 
	struct Pflx *flux
){
	/* mass */
	flux->gpp[month] = 0.0; 
	flux->epp[month] = 0.0; 
	flux->spp[month] = 0.0; 
	flux->npp[month] = 0.0; 
	
	flux->gpp_df97[month] = 0.0; 

	flux->rfm[month] = 0.0; 
	flux->rcm[month] = 0.0; 
	flux->rrm[month] = 0.0; 
	flux->arm[month] = 0.0; 
	flux->rfg[month] = 0.0; 
	flux->rcg[month] = 0.0; 
	flux->rrg[month] = 0.0; 	
	flux->arg[month] = 0.0; 	
	flux->ar[month] = 0.0; 
			
	flux->lf[month] = 0.0; 
	flux->lc[month] = 0.0; 
	flux->lr[month] = 0.0; 
	flux->lL[month] = 0.0; 
	
	flux->lf_c[month] = 0.0; 

	flux->tpp[month] = 0.0; 
	flux->tpf[month] = 0.0; 
	flux->tpc[month] = 0.0; 		
	flux->tpr[month] = 0.0; 	
	
	flux->net_crop[month] = 0.0; 
    flux->hvst_crop[month] = 0.0;

	flux->emit_ch4_kirschbaum_mass[month] = 0.0;
	flux->emit_ch4_kirschbaum_photo[month] = 0.0;
}

/* make plant fluxes vacant ********************************/
void ghg_flux_zero(
	long month, 
	struct Flux *flux
){
	flux->lu_detr = 0.0;
	flux->lu_conv = 0.0;
	flux->lu_ten = 0.0;
	flux->lu_hund = 0.0;
	
	flux->npp_miami = 0.0;
	flux->npp_montreal = 0.0;
	flux->npp_schuur = 0.0;
	flux->npp_nceas = 0.0;
	flux->npp_chikugo = 0.0;
	flux->npp_madison_gddswsi = 0.0;
	flux->npp_madison_parwsi = 0.0;
	flux->npp_madison_tp = 0.0;
	
	flux->erod_soil = 0.0;
	flux->erod_orgmat = 0.0;
	flux->erod_carbon = 0.0;
	flux->hvst_wood = 0.0; 

	flux->voc_acetacd_g97[month] = 0.0;
	flux->voc_acetone_g97[month] = 0.0;
	flux->voc_actaldhd_g97[month] = 0.0;
	flux->voc_afarnesene[month] = 0.0;
	flux->voc_bcaryophyllene[month] = 0.0;
	flux->voc_co_g97[month] = 0.0;
	flux->voc_formacd_g97[month] = 0.0;
	flux->voc_frmardhd_g97[month] = 0.0;
	flux->voc_isopr_g97[month] = 0.0;
	flux->voc_methanl_g97[month] = 0.0;
	flux->voc_monotrp_g97[month] = 0.0;
	flux->voc_othersesqui[month] = 0.0;

	(flux->c3).n_biofix[month] = 0.0;			
	(flux->c3).uptake_no3[month] = 0.0;			
	(flux->c3).uptake_nh4[month] = 0.0;	
	(flux->c3).n_alloc_cnpy[month] = 0.0;
	(flux->c3).n_alloc_strg[month] = 0.0;
	(flux->c3).n_salvage[month] = 0.0;
	(flux->c3).n_realloc[month] = 0.0;
	(flux->c3).n_abdn_cnpy[month] = 0.0;
	(flux->c3).n_abdn_strg[month] = 0.0;

	(flux->c4).n_biofix[month] = 0.0;			
	(flux->c4).uptake_no3[month] = 0.0;			
	(flux->c4).uptake_nh4[month] = 0.0;	
	(flux->c4).n_alloc_cnpy[month] = 0.0;
	(flux->c4).n_alloc_strg[month] = 0.0;
	(flux->c4).n_salvage[month] = 0.0;
	(flux->c4).n_realloc[month] = 0.0;
	(flux->c4).n_abdn_cnpy[month] = 0.0;
	(flux->c4).n_abdn_strg[month] = 0.0;

	(flux->plant).n_biofix[month] = 0.0;			
	(flux->plant).uptake_no3[month] = 0.0;			
	(flux->plant).uptake_nh4[month] = 0.0;	
	(flux->plant).n_alloc_cnpy[month] = 0.0;
	(flux->plant).n_alloc_strg[month] = 0.0;
	(flux->plant).n_salvage[month] = 0.0;
	(flux->plant).n_realloc[month] = 0.0;
	(flux->plant).n_abdn_cnpy[month] = 0.0;
	(flux->plant).n_abdn_strg[month] = 0.0;

	(flux->soil).n_abdn[month] = 0.0;
	(flux->soil).d_n2o_ntr_ngas[month] = 0.0;			
	(flux->soil).d_n2o_dnt_ngas[month] = 0.0;			
	(flux->soil).d_n2o_ngas[month] = 0.0;			
	(flux->soil).d_n2_ngas[month] = 0.0;			
	(flux->soil).d_n2o_ntr_ngas[month] = 0.0;
	(flux->soil).d_no_casa[month] = 0.0;
	(flux->soil).d_n2_casa[month] = 0.0;			
	(flux->soil).d_n2o_casa[month] = 0.0;			
	(flux->soil).n_nh3vlt[month] = 0.0;			
	(flux->soil).n_leach[month] = 0.0;			
	(flux->soil).n_minerlz_lttr[month] = 0.0;
	(flux->soil).n_minerlz_hums[month] = 0.0;
	(flux->soil).n_nitrif[month] = 0.0;
	(flux->soil).n_immbl_no3[month] = 0.0;
	(flux->soil).n_immbl_nh4[month] = 0.0;
	(flux->soil).n_mcrb_abdn[month] = 0.0;
	(flux->soil).doc_boyer[month] = 0.0;
	(flux->soil).n_fertin[month] = 0.0;
	(flux->soil).n_manurein[month] = 0.0;
	
	(flux->soil).ch4prod_wetland_cao[month] = 0.0;
	(flux->soil).ch4oxy_wetland_cao[month] = 0.0;
	(flux->soil).ch4flux_wetland_cao[month] = 0.0;
	(flux->soil).ch4prod_paddy_cao[month] = 0.0;
	(flux->soil).ch4oxy_paddy_cao[month] = 0.0;
	(flux->soil).ch4flux_paddy_cao[month] = 0.0;
	
	(flux->soil).ch4_wetland_wh_plant[month] = 0.0;
	(flux->soil).ch4_wetland_wh_ebull[month] = 0.0;
	(flux->soil).ch4_wetland_wh_diff[month] = 0.0;
	(flux->soil).ch4_wetland_wh_release[month] = 0.0;

	(flux->soil).ch4_paddy_wh_plant[month] = 0.0;
	(flux->soil).ch4_paddy_wh_ebull[month] = 0.0;
	(flux->soil).ch4_paddy_wh_diff[month] = 0.0;
	(flux->soil).ch4_paddy_wh_release[month] = 0.0;
}

/* make bare land without plant and soil *******************************/
void vanish(
	struct Mass *mass, 
	struct Flux *flux
){
	long k;
	
	/* tentative mass values */
	(mass->plant).fol = (mass->c3).fol = (mass->c4).fol = 0.0;
	(mass->plant).stm = (mass->c3).stm = (mass->c4).stm = 0.0;
	(mass->plant).rot = (mass->c3).rot = (mass->c4).rot = 0.0;
	(mass->soil).ltr = 0.0;
	(mass->soil).msl = 0.0;
	(mass->soil).doc = 0.0;
	
	for(k = 0;k<ASTEP;k++){
		/* monthly mass values */
		(mass->c3).mfol[k] = (mass->c4).mfol[k] = (mass->plant).mfol[k] = 0.0;
		(mass->c3).lai[k] = (mass->c4).lai[k] = (mass->plant).lai[k] = 0.0;
		(mass->c3).lai0[k] = (mass->c4).lai0[k] = (mass->plant).lai0[k] = 0.0;
		(mass->c3).mstm[k] = (mass->c4).mstm[k] = (mass->plant).mstm[k] = 0.0;
		(mass->c3).mrot[k] = (mass->c4).mrot[k] = (mass->plant).mrot[k] = 0.0;
		(mass->c3).plant[k] = (mass->c4).plant[k] = (mass->plant).plant[k] = 0.0;
		(mass->soil).ltr_m[k] = 0.0;
		(mass->soil).msl_m[k] = 0.0;
		(mass->soil).soil[k] = 0.0;
		(mass->soil).doc_m[k] = 0.0;
		
		/* fluxes */
		plant_flux_zero(k, &(flux->c3));
		plant_flux_zero(k, &(flux->c4));
		plant_flux_zero(k, &(flux->plant));		
		ghg_flux_zero(k, flux);
								
		(flux->soil).rl[k] = 0.0;
		(flux->soil).rh[k] = 0.0;
		(flux->soil).sf[k] = 0.0;
		(flux->soil).hr[k] = 0.0;
		
		flux->nep[k] = 0.0;
		flux->nbp[k] = 0.0;
		flux->ncb[k] = 0.0;
		flux->lL0[k] = 0.0;
	}
	mass->lai_p = 0.0;
	flux->efflux_p = 0.0;
	
	flux->lu_detr = 0.0;
	flux->lu_conv = 0.0;
	flux->lu_ten = 0.0;
	flux->lu_hund = 0.0;
	
	flux->npp_miami = 0.0;
	flux->npp_montreal = 0.0;
	flux->npp_schuur = 0.0;
	flux->npp_nceas = 0.0;
	
	flux->erod_soil = 0.0;
	flux->erod_orgmat = 0.0;
	flux->erod_carbon = 0.0;
	flux->hvst_wood = 0.0; 
}

/* make the biome type zero **********************************/
void vlzero(
	struct Grid *grid, 
	struct Pmas *mass, 
	struct Pflx *flux
){

	mass->fol = 0.0;
	mass->stm = 0.0;
	mass->rot = 0.0;
	
	mass->lai[grid->m] = mass->lai0[grid->m] = 0.0;

	mass->mfol[grid->m] = 0.0;
	mass->mstm[grid->m] = 0.0;
	mass->mrot[grid->m] = 0.0;
	mass->plant[grid->m] = 0.0;
	
	plant_flux_zero(grid->m, flux);
}

/* initialize carbon isotope variables ************************************/
void f_init_c_isotpes(
	struct Grid *grid, 
	struct Flux *flux, 
	struct Echar *echar, 
	struct Mass *mass
){
	long f;
	
	for(f = 0;f<ASTEP;f++){
		/* d13C **********************************************************************/
		/* mass */
		(mass->c3).d13c_fol = (mass->c4).d13c_fol = (mass->plant).d13c_fol = grid->d13c_bco2[f]; 	
		(mass->c3).d13c_stm = (mass->c4).d13c_stm = (mass->plant).d13c_stm = grid->d13c_bco2[f]; 	
		(mass->c3).d13c_rot = (mass->c4).d13c_rot = (mass->plant).d13c_rot = grid->d13c_bco2[f]; 	
		(mass->soil).d13c_ltr = grid->d13c_bco2[f]; 	
		(mass->soil).d13c_msl = grid->d13c_bco2[f]; 	
		
		flux->d13c_efflux_p = grid->d13c_bco2[f]; 
		
		/* ecophysiology */
		(echar->c3).photo_13c_frac[f] = 0.0;
		(echar->c4).photo_13c_frac[f] = 0.0;
		
		/* mass */
		(mass->c3).d13c_mfol[f] = (mass->c4).d13c_mfol[f] = (mass->plant).d13c_mfol[f] = grid->d13c_bco2[f]; 	
		(mass->c3).d13c_mstm[f] = (mass->c4).d13c_mstm[f] = (mass->plant).d13c_mstm[f] = grid->d13c_bco2[f]; 	
		(mass->c3).d13c_mrot[f] = (mass->c4).d13c_mrot[f] = (mass->plant).d13c_mrot[f] = grid->d13c_bco2[f]; 	
		(mass->c3).d13c_plant[f] = (mass->c4).d13c_plant[f] = (mass->plant).d13c_plant[f] = grid->d13c_bco2[f]; 	
		
		(mass->soil).d13c_ltr_m[f] = grid->d13c_bco2[f]; 	
		(mass->soil).d13c_msl_m[f] = grid->d13c_bco2[f]; 	
		(mass->soil).d13c_soil[f] = grid->d13c_bco2[f]; 	
		
		mass->d13c_total[f] = grid->d13c_bco2[f];
		
		/* flux */
		(flux->c3).d13c_gpp[f] = (flux->c4).d13c_gpp[f] = (flux->plant).d13c_gpp[f] = grid->d13c_bco2[f];
		(flux->c3).d13c_epp[f] = (flux->c4).d13c_epp[f] = (flux->plant).d13c_epp[f] = grid->d13c_bco2[f];
		(flux->c3).d13c_spp[f] = (flux->c4).d13c_spp[f] = (flux->plant).d13c_spp[f] = grid->d13c_bco2[f];
		(flux->c3).d13c_npp[f] = (flux->c4).d13c_npp[f] = (flux->plant).d13c_npp[f] = grid->d13c_bco2[f];
		
		(flux->c3).d13c_rfm[f] = (flux->c4).d13c_rfm[f] = (flux->plant).d13c_rfm[f] = grid->d13c_bco2[f];
		(flux->c3).d13c_rcm[f] = (flux->c4).d13c_rcm[f] = (flux->plant).d13c_rcm[f] = grid->d13c_bco2[f];
		(flux->c3).d13c_rrm[f] = (flux->c4).d13c_rrm[f] = (flux->plant).d13c_rrm[f] = grid->d13c_bco2[f];
		(flux->c3).d13c_arm[f] = (flux->c4).d13c_arm[f] = (flux->plant).d13c_arm[f] = grid->d13c_bco2[f];
		(flux->c3).d13c_rfg[f] = (flux->c4).d13c_rfg[f] = (flux->plant).d13c_rfg[f] = grid->d13c_bco2[f];
		(flux->c3).d13c_rcg[f] = (flux->c4).d13c_rcg[f] = (flux->plant).d13c_rcg[f] = grid->d13c_bco2[f];
		(flux->c3).d13c_rrg[f] = (flux->c4).d13c_rrg[f] = (flux->plant).d13c_rrg[f] = grid->d13c_bco2[f];
		(flux->c3).d13c_arg[f] = (flux->c4).d13c_arg[f] = (flux->plant).d13c_arg[f] = grid->d13c_bco2[f]; 
		(flux->c3).d13c_ar[f] = (flux->c4).d13c_ar[f] = (flux->plant).d13c_ar[f] = grid->d13c_bco2[f]; 
				
		(flux->c3).d13c_lf[f] = (flux->c4).d13c_lf[f] = (flux->plant).d13c_lf[f] = grid->d13c_bco2[f]; 
		(flux->c3).d13c_lc[f] = (flux->c4).d13c_lc[f] = (flux->plant).d13c_lc[f] = grid->d13c_bco2[f]; 
		(flux->c3).d13c_lr[f] = (flux->c4).d13c_lr[f] = (flux->plant).d13c_lr[f] = grid->d13c_bco2[f]; 
		(flux->c3).d13c_lL[f] = (flux->c4).d13c_lL[f] = (flux->plant).d13c_lL[f] = grid->d13c_bco2[f]; 
		
		(flux->c3).d13c_lf_c[f] = (flux->c4).d13c_lf_c[f] = (flux->plant).d13c_lf_c[f] = grid->d13c_bco2[f]; 

		(flux->c3).d13c_tpp[f] = (flux->c4).d13c_tpp[f] = (flux->plant).d13c_tpp[f] = grid->d13c_bco2[f]; 
		(flux->c3).d13c_tpf[f] = (flux->c4).d13c_tpf[f] = (flux->plant).d13c_tpf[f] = grid->d13c_bco2[f]; 
		(flux->c3).d13c_tpc[f] = (flux->c4).d13c_tpc[f] = (flux->plant).d13c_tpc[f] = grid->d13c_bco2[f]; 
		(flux->c3).d13c_tpr[f] = (flux->c4).d13c_tpr[f] = (flux->plant).d13c_tpr[f] = grid->d13c_bco2[f]; 
		
		(flux->c3).d13c_net_crop[f] = (flux->c4).d13c_net_crop[f] = (flux->plant).d13c_net_crop[f] = grid->d13c_bco2[f]; 	
		
		(flux->soil).d13c_lL[f] = grid->d13c_bco2[f]; 
		(flux->soil).d13c_rl[f] = grid->d13c_bco2[f]; 
		(flux->soil).d13c_rh[f] = grid->d13c_bco2[f]; 
		(flux->soil).d13c_hr[f] = grid->d13c_bco2[f]; 
		(flux->soil).d13c_sf[f] = grid->d13c_bco2[f]; 
		
		flux->d13c_nep[f] = grid->d13c_bco2[f]; 
		flux->d13c_ncb[f] = grid->d13c_bco2[f]; 
		
		/* D14C **********************************************************************/
		(mass->c3).d14c_fol = (mass->c4).d14c_fol = (mass->plant).d14c_fol = 0.0; 	
		(mass->c3).d14c_stm = (mass->c4).d14c_stm = (mass->plant).d14c_stm = 0.0; 	
		(mass->c3).d14c_rot = (mass->c4).d14c_rot = (mass->plant).d14c_rot = 0.0; 	
		(mass->soil).d14c_ltr = 0.0; 	
		(mass->soil).d14c_msl = 0.0; 	
		
		(mass->c3).d14c_mfol[f] = (mass->c4).d14c_mfol[f] = (mass->plant).d14c_mfol[f] = 0.0; 	
		(mass->c3).d14c_mstm[f] = (mass->c4).d14c_mstm[f] = (mass->plant).d14c_mstm[f] = 0.0; 	
		(mass->c3).d14c_mrot[f] = (mass->c4).d14c_mrot[f] = (mass->plant).d14c_mrot[f] = 0.0; 	
		(mass->soil).d14c_ltr_m[f] = 0.0; 
		(mass->soil).d14c_msl_m[f] = 0.0; 
		
		(flux->c3).d14c_gpp[f] = (flux->c3).d14c_lL[f] = 0.0;
		(flux->c4).d14c_gpp[f] = (flux->c4).d14c_lL[f] = 0.0;
		(flux->plant).d14c_gpp[f] = (flux->plant).d14c_lL[f] = 0.0;
		(flux->soil).d14c_lL[f] = 0.0;
		
		flux->d14c_sr[f] = 0.0; 
		flux->d14c_er[f] = 0.0; 
	}
}
