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

/*******************************************************************************/
void f_output_file_open(
						short vtype, 
						short zone,
						char s_date[25], 
						char s_case[25], 
						char filename[100], 
						FILE *fp[OFILES]
						){
	char num[4];
	char svtype[8]; /* corrected: A.Ito and E.Kato (2009/08/16) */
	
	switch(vtype){
		case 1:
			strcpy(svtype, "olson_");
			break;
		case 2:
			strcpy(svtype, "sage_");
			break;
		case 3:
			strcpy(svtype, "crop_");
			break;
		default:
			printf("Bad vegetation type !!!!!\n");
			exit (1);
	}
	
	strcpy(filename, s_date);	
	strcat(filename, s_case);
	strcat(filename, svtype);
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
	strcat(filename, svtype);
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
	strcat(filename, svtype);
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
	strcat(filename, svtype);
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
	strcat(filename, svtype);
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
	strcat(filename, svtype);
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
	strcat(filename, svtype);
	strcat(filename, "hydmet_");
	snprintf(num, 4, "%02d", zone);
	strcat(filename, num);
	strcat(filename, ".dat");
	if( (fp[6] = fopen(filename,"wt"))==NULL ){  
		printf("No Hydrometorological output file\n");  
		exit(1); 
	}
	
	strcpy(filename, s_date);
	strcat(filename, s_case);
	strcat(filename, svtype);
	strcat(filename, "cistp_");
	snprintf(num, 4, "%02d", zone);
	strcat(filename, num);
	strcat(filename, ".dat");
	if( (fp[7] = fopen(filename,"wt"))==NULL ){  
		printf("No Carbon isotope output file\n");  
		exit(1); 
	}
}

/* make output file for stable state *****************************************************/
void publish_cbud(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux, 
	FILE *result
){
	long f;
	
	/* corrected: A.Ito and E.Kato (2009/08/16) */
	fprintf(result,"%.3lf %.3lf %ld %ld %ld ",grid->lat, grid->lon, loct->time, grid->n_olson, grid->n_sage);	
	/* fprintf(result,"%lf %lf %lf %lf ", grid->whc30, grid->whc, grid->sd, grid->hyd_cond); */
	fprintf(result,"\n");
	
	for(f=0;f<ASTEP;f++){
		fprintf(result,"%.2lf ", grid->bco2[f]);	
		fprintf(result,"%lf %lf ", loct->c3ptn[f], loct->c4ptn[f]);
		
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
		fprintf(result,"%.2lf ", (flux->soil).hr[f]);
		
		fprintf(result,"%.2lf ", flux->nep[f]);
		fprintf(result,"%.2lf ", flux->nbp[f]);	/* added (2010/02/01 by A.Ito) */
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
	fprintf(result,"%lf ", flux->npp_schuur); /* added 2008/09/08 by A.Ito */
	fprintf(result,"%lf ", flux->npp_nceas);  /* added 2008/09/08 by A.Ito */
	
	fprintf(result,"%.2lf ", flux->erod_soil);
	fprintf(result,"%.2lf ", flux->erod_orgmat);
	fprintf(result,"%.2lf ", flux->erod_carbon);
	
/*	for(f=0;f<12;f++){
		fprintf(result,"%lf %lf ",grid->d13c_bco2[f],loct->d13c_aco2[f]);

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
	
	for(f=0;f<ASTEP;f++){
		h_tmp[year] += grid->tmp_2m[f]* MDN[f]/365.0 * grid->area;
		h_pre[year] += grid->prate_sfc[f] * grid->area;
		h_dswr[year] += grid->gl_rad[f]* MDN[f]/365.0 * grid->area;
		h_aet[year] += (loct->evpr[f] + loct->trspr[f] + loct->incep[f])*0.1 * grid->area;
		h_rof[year] += loct->ro2[f] * grid->area;
		if(DF97==1){
			h_gpp[year] += (flux->plant).gpp_df97[f] * grid->area;
			h_gpp_c4[year] += (flux->c4).gpp_df97[f] * grid->area*loct->c4ptn[f];
		}else{
			h_gpp[year] += (flux->plant).gpp[f] * grid->area;
			h_gpp_c4[year] += (flux->c4).gpp[f] * grid->area*loct->c4ptn[f];
		}
		h_npp[year] += (flux->plant).npp[f] * grid->area;
		h_nep[year] += flux->nep[f] * grid->area;
		h_plant[year] += ((mass->plant).mfol[f] + (mass->plant).mstm[f] + 
			(mass->plant).mrot[f]) * MDN[f] /365.0 * grid->area;
		h_soil[year] += ((mass->soil).ltr_m[f] + 
			(mass->soil).msl_m[f])* MDN[f]/365.0 * grid->area;
		h_sr[year] += ((flux->plant).rrm[f] + (flux->plant).rrg[f] + (flux->soil).hr[f]) * grid->area;
		h_doc[year] += (flux->soil).doc_boyer[f] * grid->area;
		
		h_gpp_df97[year] += (flux->plant).gpp_df97[f] * grid->area;

		/* GHG */
		h_ch4ox1[year] += (flux->soil).ch4oxy_ridg[f] * grid->area *10000.0/1000.0;
		h_ch4ox2[year] += (flux->soil).ch4oxy_casa[f] * grid->area *10000.0/1000.0;
		h_ch4ox3[year] += (flux->soil).ch4oxy_delgrosso[f] * grid->area *10000.0/1000.0;
		h_ch4ox4[year] += (flux->soil).ch4oxy_curry[f] * grid->area *10000.0/1000.0;		
		
		h_ch4emit_cao_paddy[year] += (flux->soil).ch4flux_paddy_cao[f] * grid->area *10000.0/1000.0;
		h_ch4emit_cao_wetland[year] += (flux->soil).ch4flux_wetland_cao[f] * grid->area *10000.0/1000.0;

		h_ch4emit_paddy_wh_diff[year] += (flux->soil).ch4_paddy_wh_diff[f] * grid->area *10000.0/1000.0;
		h_ch4emit_paddy_wh_plant[year] += (flux->soil).ch4_paddy_wh_plant[f] * grid->area *10000.0/1000.0;
		h_ch4emit_paddy_wh_ebbl[year] += (flux->soil).ch4_paddy_wh_ebull[f] * grid->area *10000.0/1000.0;
		h_ch4emit_paddy_wh_release[year] += (flux->soil).ch4_paddy_wh_release[f] * grid->area *10000.0/1000.0;

		h_ch4emit_wetland_wh_diff[year] += (flux->soil).ch4_wetland_wh_diff[f] * grid->area *10000.0/1000.0;
		h_ch4emit_wetland_wh_plant[year] += (flux->soil).ch4_wetland_wh_plant[f] * grid->area *10000.0/1000.0;
		h_ch4emit_wetland_wh_ebbl[year] += (flux->soil).ch4_wetland_wh_ebull[f] * grid->area *10000.0/1000.0;
		h_ch4emit_wetland_wh_release[year] += (flux->soil).ch4_wetland_wh_release[f] * grid->area *10000.0/1000.0;
		
		h_n2o_emit_ngas[year] += (flux->soil).d_n2o_ngas[f] * grid->area;
		h_n2_emit_ngas[year] += (flux->soil).d_n2_ngas[f] * grid->area;
		h_n2o_emit_casa[year] += (flux->soil).d_n2o_casa[f] * grid->area;
		h_n2_emit_casa[year] += (flux->soil).d_n2_casa[f] * grid->area;
		h_no_emit_casa[year] += (flux->soil).d_no_casa[f] * grid->area;
		h_nh3_emit[year] += (flux->soil).n_nh3vlt[f] * grid->area;		
		h_ch4_emit_mass[year] += (flux->plant).emit_ch4_kirschbaum_mass[f]*10000.0*grid->area;
		h_ch4_emit_photo[year] += (flux->plant).emit_ch4_kirschbaum_photo[f]*10000.0*grid->area;
		h_n2_biofix[year] += (flux->plant).n_biofix[f]*grid->area;
		h_n2o_n_emit_ngas[year] += (flux->soil).d_n2o_ntr_ngas[f] * grid->area;
		h_n2o_d_emit_ngas[year] += (flux->soil).d_n2o_dnt_ngas[f] * grid->area;
		
		if(CALC_OLSON == 1){ /* added by A.Ito (2009/06/16) */
			if(grid->veg_olson==29 || grid->veg_olson==30 || grid->veg_olson==31 || grid->veg_olson==32){
				h_n2o_emit_ngas_agr[year] += (flux->soil).d_n2o_ngas[f] * grid->area;
				h_n2o_emit_casa_agr[year] += (flux->soil).d_n2o_casa[f] * grid->area;
				h_nh3_emit_agr[year] += (flux->soil).n_nh3vlt[f] * grid->area;
			}
		}
		
		/* biomass burninig *****************/
		h_burnt_area[year] += flux->a_burnt[f]*grid->area;
		h_bioburn_co2[year] += (flux->bb_co2_litter[f]+flux->bb_co2_leaf[f]+flux->bb_co2_wood[f]+flux->bb_co2_root[f]) * grid->area;
		h_bioburn_co[year] += (flux->bb_co_litter[f]+flux->bb_co_leaf[f]+flux->bb_co_wood[f]+flux->bb_co_root[f]) * grid->area;
		h_bioburn_ch4[year] += (flux->bb_ch4_litter[f]+flux->bb_ch4_leaf[f]+flux->bb_ch4_wood[f]+flux->bb_ch4_root[f]) * grid->area;
		h_bioburn_nmhc[year] += (flux->bb_nmhc_litter[f]+flux->bb_nmhc_leaf[f]+flux->bb_nmhc_wood[f]+flux->bb_nmhc_root[f]) * grid->area;
		h_bioburn_oc[year] += (flux->bb_oc_litter[f]+flux->bb_oc_leaf[f]+flux->bb_oc_wood[f]+flux->bb_oc_root[f]) * grid->area;
		h_bioburn_bc[year] += (flux->bb_bc_litter[f]+flux->bb_bc_leaf[f]+flux->bb_bc_wood[f]+flux->bb_bc_root[f]) * grid->area;
		h_bioburn_nox[year] += (flux->bb_nox_litter[f]+flux->bb_nox_leaf[f]+flux->bb_nox_wood[f]+flux->bb_nox_root[f]) * grid->area;
		h_bioburn_so2[year] += (flux->bb_so2_litter[f]+flux->bb_so2_leaf[f]+flux->bb_so2_wood[f]+flux->bb_so2_root[f]) * grid->area;
		h_bioburn_pm25[year] += (flux->bb_pm25_litter[f]+flux->bb_pm25_leaf[f]+flux->bb_pm25_wood[f]+flux->bb_pm25_root[f]) * grid->area;
		h_bioburn_tpm[year] += (flux->bb_tpm_litter[f]+flux->bb_tpm_leaf[f]+flux->bb_tpm_wood[f]+flux->bb_tpm_root[f]) * grid->area;
		h_bioburn_tec[year] += (flux->bb_tec_litter[f]+flux->bb_tec_leaf[f]+flux->bb_tec_wood[f]+flux->bb_tec_root[f]) * grid->area;

		/* VOC, g C ha-1 month-1 **************/
		h_voc_isopr_g97[year] += flux->voc_isopr_g97[f] * grid->area *10000.0/1000000.0;
		h_voc_monotrp_g97[year] += flux->voc_monotrp_g97[f] * grid->area *10000.0/1000000.0;
		h_voc_methanl_g97[year] += flux->voc_methanl_g97[f] * grid->area *10000.0/1000000.0;
		h_voc_acetone_g97[year] += flux->voc_acetone_g97[f] * grid->area *10000.0/1000000.0;
		h_voc_actaldhd_g97[year] += flux->voc_actaldhd_g97[f] * grid->area *10000.0/1000000.0;
		h_voc_frmardhd_g97[year] += flux->voc_frmardhd_g97[f] * grid->area *10000.0/1000000.0;
		h_voc_formacd_g97[year] += flux->voc_formacd_g97[f] * grid->area *10000.0/1000000.0;
		h_voc_acetacd_g97[year] += flux->voc_acetacd_g97[f] * grid->area *10000.0/1000000.0;
		h_voc_co_g97[year] += flux->voc_co_g97[f] * grid->area *10000.0/1000000.0;
		
		/* d13c & d14c : added by A.Ito (2009/07/15) ***************/
		ci_aco2_d13c[year] = d13c_addition(loct->d13c_aco2[f],loct->aco2[f]*grid->area * MDN[f]/365.0, ci_aco2_d13c[year],ci_aco2[year]);
		ci_aco2_d14c[year] = (grid->d14c_bco2[f]*loct->aco2[f]*grid->area * MDN[f]/365.0 + ci_aco2_d14c[year]*ci_aco2[year]) / 
							(loct->aco2[f]*grid->area * MDN[f]/365.0 + ci_aco2[year]);
		ci_aco2[year] += loct->aco2[f]*grid->area * MDN[f]/365.0;
		
		/* GPP */
		if(DF97==1){
			if((flux->plant).gpp_df97[f] > 0.0){
				ci_gpp_d13c[year] = d13c_addition((flux->plant).d13c_gpp[f], (flux->plant).gpp_df97[f] * grid->area, ci_gpp_d13c[year], ci_gpp[year]);
				ci_gpp_d14c[year] = ((flux->plant).d14c_gpp[f] * (flux->plant).gpp_df97[f] * grid->area + ci_gpp_d14c[year] * ci_gpp[year]) / 
									((flux->plant).gpp_df97[f] * grid->area + ci_gpp[year]);
				ci_gpp[year] += (flux->plant).gpp_df97[f] * grid->area;
			}
		}else{
			if((flux->plant).gpp[f] > 0.0){
				ci_gpp_d13c[year] = d13c_addition((flux->plant).d13c_gpp[f], (flux->plant).gpp[f] * grid->area, ci_gpp_d13c[year], ci_gpp[year]);
				ci_gpp_d14c[year] = ((flux->plant).d14c_gpp[f] * (flux->plant).gpp[f] * grid->area + ci_gpp_d14c[year] * ci_gpp[year]) / 
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
		rh_apar[grid->reg_g][year] += grid->dlen[f] * 0.5 * loct->fapar_mono[f]* 
									grid->par[f] * MDN[f]/365.0 * grid->area;
		
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
												 (flux->soil).ch4_wetland_wh_ebull[f] + (flux->soil).ch4_wetland_wh_release[f]) * grid->area *10000.0/1000.0;
		rh_ch4emit_wh_paddy[grid->reg_g][year] += ((flux->soil).ch4_paddy_wh_diff[f] + (flux->soil).ch4_paddy_wh_plant[f] + 
													 (flux->soil).ch4_paddy_wh_ebull[f] + (flux->soil).ch4_paddy_wh_release[f]) * grid->area *10000.0/1000.0;
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
	
	/* erosion */
	h_ersn_c[year] += flux->erod_carbon * grid->area;
	h_agrersn_c[year] += flux->erod_carbon_crop * grid->area;

	h_agrarea[year] += grid->area * grid->f_crop_con;
	h_luc[year] += (flux->lu_conv + flux->lu_ten + flux->lu_hund) * grid->area;
	
	h_luc_1[year] += flux->lu_conv * grid->area;
	h_luc_2[year] += flux->lu_ten * grid->area;
	h_luc_3[year] += flux->lu_hund * grid->area;
	
	rh_luc[grid->reg_g][year] += (flux->lu_conv + flux->lu_ten + flux->lu_hund) * grid->area;
}

/*********************************************************************************/
void f_output_result(
	long year, 
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux, 
	FILE *fp_o[OFILES]
){
	long f;
	double pre_ann, hrl_ann, hrm_ann, ltr_ann, msl_ann;
	extern double MDN[ASTEP];
	
	pre_ann = hrl_ann = hrm_ann = ltr_ann = msl_ann = 0.0;
	
	for(f=0;f<ASTEP;f++){
		pre_ann += grid->prate_sfc[f];
		hrl_ann += (flux->soil).rl[f];
		hrm_ann += (flux->soil).rh[f];
		ltr_ann += (mass->soil).ltr_m[f] * MDN[f] /365.0;
		msl_ann += (mass->soil).msl_m[f] * MDN[f] /365.0;
	}

	/* carbon budget */
	fprintf(fp_o[0],"%ld ", year);
	for(f=0;f<ASTEP;f++){
		fprintf(fp_o[0],"%.4lf ", grid->tmp_2m[f]);
		fprintf(fp_o[0],"%.4lf ", grid->prate_sfc[f]);
		fprintf(fp_o[0],"%.4lf ", grid->gl_rad[f]);
		
		fprintf(fp_o[0],"%.4lf ", loct->evpr[f] + loct->trspr[f] + loct->incep[f]);
		fprintf(fp_o[0],"%.4lf ", loct->ro2[f]);
		
		fprintf(fp_o[0],"%.4lf ", loct->c3ptn[f]*(mass->c3).lai[f]+loct->c4ptn[f]*(mass->c4).lai[f]);
		fprintf(fp_o[0],"%.4lf ", (mass->plant).mfol[f]);
		fprintf(fp_o[0],"%.4lf ", (mass->plant).mstm[f]);
		fprintf(fp_o[0],"%.4lf ", (mass->plant).mrot[f]);
		fprintf(fp_o[0],"%.4lf ", (mass->soil).ltr_m[f]);
		fprintf(fp_o[0],"%.4lf ", (mass->soil).msl_m[f]);
		
		if(DF97==1){
			fprintf(fp_o[0],"%.4lf ", (flux->plant).gpp_df97[f]);
		}else{
			fprintf(fp_o[0],"%.4lf ", (flux->plant).gpp[f]);
		}
		fprintf(fp_o[0],"%.4lf ", (flux->plant).ar[f]);
		fprintf(fp_o[0],"%.4lf ", (flux->soil).hr[f]);
		fprintf(fp_o[0],"%.4lf ",  flux->nep[f]);
		fprintf(fp_o[0],"%.4lf ",  flux->nbp[f]);	/* added by A.Ito (2010/01/20) */
		fprintf(fp_o[0],"%.4lf ",  flux->ncb[f]);
	}
	fprintf(fp_o[0],"%.4lf ", flux->lu_detr);
	fprintf(fp_o[0],"%.4lf ", flux->lu_conv);
	fprintf(fp_o[0],"%.4lf ", flux->lu_ten);
	fprintf(fp_o[0],"%.4lf ", flux->lu_hund); /* */
	fprintf(fp_o[0],"\n");
	
	/* nitrogen */
	/* parameter added by A.Ito (2009/12/22) */
	fprintf(fp_o[1],"%ld ", year);
	for(f=0;f<ASTEP;f++){
		fprintf(fp_o[1],"%.3lf ", (mass->c3).n_cnpy_m[f]);
		fprintf(fp_o[1],"%.3lf ", (mass->c3).n_strg_m[f]);
		fprintf(fp_o[1],"%.3lf ", (mass->c4).n_cnpy_m[f]);
		fprintf(fp_o[1],"%.3lf ", (mass->c4).n_strg_m[f]);
		
		fprintf(fp_o[1],"%.3lf ", (mass->soil).n_no3_m[f]);
		fprintf(fp_o[1],"%.3lf ", (mass->soil).n_nh4_m[f]);
		fprintf(fp_o[1],"%.3lf ", (mass->soil).n_mcrb_m[f]);
		fprintf(fp_o[1],"%.3lf ", (mass->soil).n_lttr_m[f]);
		fprintf(fp_o[1],"%.3lf ", (mass->soil).n_hums_m[f]);
		
		fprintf(fp_o[1],"%.3lf ", (flux->plant).n_biofix[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->plant).uptake_no3[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->plant).uptake_nh4[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->plant).n_abdn_cnpy[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->plant).n_abdn_strg[f]);

		fprintf(fp_o[1],"%.3lf ", (flux->soil).n_leach[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->soil).n_nh3vlt[f]); /* */
		fprintf(fp_o[1],"%.3lf ", (flux->soil).d_n2o_ngas[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->soil).d_n2_ngas[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->soil).d_n2o_casa[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->soil).d_n2_casa[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->soil).d_no_casa[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->soil).d_n2o_ntr_ngas[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->soil).d_n2o_dnt_ngas[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->soil).n_nitrif[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->soil).n_immbl[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->soil).n_mcrb_abdn[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->soil).n_minerlz_lttr[f]);
		fprintf(fp_o[1],"%.3lf ", (flux->soil).n_minerlz_hums[f]);

		fprintf(fp_o[1],"%.3lf ", loct->n_frtlz_in*1000.0);
		fprintf(fp_o[1],"%.3lf ", loct->depo_no3[f]);
		fprintf(fp_o[1],"%.3lf ", loct->depo_nh4[f]);
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
	for(f=0;f<ASTEP;f++){
		fprintf(fp_o[2],"%lf ", loct->ro2[f]);
		fprintf(fp_o[2],"%lf ", (flux->soil).doc_boyer[f]);
	} /* */
	fprintf(fp_o[2],"\n");
	
	/* GHG & trace gases */
	fprintf(fp_o[3],"%ld ", year);
	for(f=0;f<ASTEP;f++){
		fprintf(fp_o[3],"%.3lf ", loct->msww[f]);
		fprintf(fp_o[3],"%.3lf ", loct->wfps[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4oxy_ridg[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4oxy_casa[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4oxy_delgrosso[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4oxy_curry[f]);
		
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4flux_paddy_cao[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4flux_wetland_cao[f]);
		
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4_wetland_wh_plant[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4_wetland_wh_ebull[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4_wetland_wh_diff[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4_wetland_wh_release[f]);
		
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4_paddy_wh_plant[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4_paddy_wh_ebull[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4_paddy_wh_diff[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).ch4_paddy_wh_release[f]);
		
		fprintf(fp_o[3],"%.3lf ", (flux->soil).d_n2o_ngas[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).d_n2_ngas[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).d_n2o_casa[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).d_n2_casa[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).d_no_casa[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->soil).n_nh3vlt[f]);
		
		fprintf(fp_o[3],"%.3lf ", (flux->plant).emit_ch4_kirschbaum_mass[f]);
		fprintf(fp_o[3],"%.3lf ", (flux->plant).emit_ch4_kirschbaum_photo[f]); /* */
	}
	fprintf(fp_o[3],"\n");

	/* biomass burning */
	fprintf(fp_o[4],"%ld ", year);
	for(f=0;f<ASTEP;f++){
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
		fprintf(fp_o[4],"%.3lf ", flux->bb_tec_litter[f]+flux->bb_tec_leaf[f]+flux->bb_tec_wood[f]+flux->bb_tec_root[f]);  /* */
	}
	fprintf(fp_o[4],"\n");

	/* VOC */
	fprintf(fp_o[5],"%ld ", year);
	for(f=0;f<ASTEP;f++){
		fprintf(fp_o[5],"%.3lf ", flux->voc_isopr_g97[f]);
		fprintf(fp_o[5],"%.3lf ", flux->voc_monotrp_g97[f]);
		fprintf(fp_o[5],"%.3lf ", flux->voc_methanl_g97[f]);
		fprintf(fp_o[5],"%.3lf ", flux->voc_acetone_g97[f]);
		fprintf(fp_o[5],"%.3lf ", flux->voc_actaldhd_g97[f]);
		fprintf(fp_o[5],"%.3lf ", flux->voc_frmardhd_g97[f]);
		fprintf(fp_o[5],"%.3lf ", flux->voc_formacd_g97[f]);
		fprintf(fp_o[5],"%.3lf ", flux->voc_acetacd_g97[f]);
		fprintf(fp_o[5],"%.3lf ", flux->voc_co_g97[f]); /* */
	}
	fprintf(fp_o[5],"\n");
	
	/* hydrometeorology */
	fprintf(fp_o[6],"%ld ", year);
	for(f=0;f<ASTEP;f++){
		fprintf(fp_o[6],"%.3lf ", grid->tmp_2m[f]);
		fprintf(fp_o[6],"%.2lf ", grid->prate_sfc[f]);
		fprintf(fp_o[6],"%.2lf ", grid->gl_rad[f]);
		fprintf(fp_o[6],"%.2lf ", grid->par[f]);
		fprintf(fp_o[6],"%.2lf ", grid->par_bp[f]);
		fprintf(fp_o[6],"%.2lf ", grid->par_dp[f]);
		fprintf(fp_o[6],"%.2lf ", loct->vpd[f]);
		fprintf(fp_o[6],"%.2lf ", loct->rad_net_long[f]);
		fprintf(fp_o[6],"%.2lf ", loct->rad_net_short[f]);
		fprintf(fp_o[6],"%.2lf ", loct->msw30[f]);
		fprintf(fp_o[6],"%.2lf ", loct->msww[f]);
		fprintf(fp_o[6],"%.2lf ", loct->incep[f]);
		fprintf(fp_o[6],"%.2lf ", loct->evpr[f]);
		fprintf(fp_o[6],"%.2lf ", loct->trspr[f]);
		fprintf(fp_o[6],"%.2lf ", loct->ro2[f]); /* */
	}
	fprintf(fp_o[6],"\n");
	
	/* carbon isotopes: added by A.Ito (2009/07/14) */
	fprintf(fp_o[7],"%ld ", year);
	for(f=0;f<ASTEP;f++){
		fprintf(fp_o[7],"%.4lf ", loct->c4ptn[f]);
		
		fprintf(fp_o[7],"%.3lf ", grid->d13c_bco2[f]);
		fprintf(fp_o[7],"%.3lf ", loct->d13c_aco2[f]);
		fprintf(fp_o[7],"%.3lf ", grid->d14c_bco2[f]);
		
		fprintf(fp_o[7],"%.3lf ", (echar->c3).photo_13c_frac[f]);
		fprintf(fp_o[7],"%.3lf ", (echar->c4).photo_13c_frac[f]);
		
		fprintf(fp_o[7],"%.3lf ", flux->d13c_sr[f]);
		fprintf(fp_o[7],"%.3lf ", flux->d13c_er[f]);
		fprintf(fp_o[7],"%.3lf ", flux->d14c_sr[f]);
		fprintf(fp_o[7],"%.3lf ", flux->d14c_er[f]);
		
		fprintf(fp_o[7],"%.4lf ", (flux->plant).d13c_gpp[f]*(flux->plant).gpp[f] - flux->d13c_er[f]*flux->er[f]);
		fprintf(fp_o[7],"%.4lf ", (flux->plant).d14c_gpp[f]*(flux->plant).gpp[f] - flux->d14c_er[f]*flux->er[f]);
		
		fprintf(fp_o[7],"%.3lf ", (mass->c3).d13c_mfol[f]);
		fprintf(fp_o[7],"%.3lf ", (mass->c3).d13c_mstm[f]);
		fprintf(fp_o[7],"%.3lf ", (mass->c3).d13c_mrot[f]);
		fprintf(fp_o[7],"%.3lf ", (mass->c4).d13c_mfol[f]);
		fprintf(fp_o[7],"%.3lf ", (mass->c4).d13c_mstm[f]);
		fprintf(fp_o[7],"%.3lf ", (mass->c4).d13c_mrot[f]);
		fprintf(fp_o[7],"%.3lf ", (mass->soil).d13c_ltr_m[f]);
		fprintf(fp_o[7],"%.3lf ", (mass->soil).d13c_msl_m[f]);
		
		fprintf(fp_o[7],"%.3lf ", (mass->c3).d14c_mfol[f]);
		fprintf(fp_o[7],"%.3lf ", (mass->c3).d14c_mstm[f]);
		fprintf(fp_o[7],"%.3lf ", (mass->c3).d14c_mrot[f]);
		fprintf(fp_o[7],"%.3lf ", (mass->c4).d14c_mfol[f]);
		fprintf(fp_o[7],"%.3lf ", (mass->c4).d14c_mstm[f]);
		fprintf(fp_o[7],"%.3lf ", (mass->c4).d14c_mrot[f]);
		fprintf(fp_o[7],"%.3lf ", (mass->soil).d14c_ltr_m[f]);
		fprintf(fp_o[7],"%.3lf ", (mass->soil).d14c_msl_m[f]); /* */
	}
	fprintf(fp_o[7],"\n");
}

/*********************************************************************/
void f_glosum_output(
	char s_date[25], 
	char s_case[25]
){
	long h, i;
	char filename[100];
	FILE *fp_glsum;

	/* global summary ******************************************/
	strcpy(filename, s_date);
	strcat(filename, s_case);
	strcat(filename, "glsum.dat");
	
	fp_glsum = fopen(filename,"wt");
	fprintf(fp_glsum,"%.2lf %.2lf\n", go_landarea, gs_landarea);
	for(h=0;h<HIST;h++){
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
	for(h=0;h<VEG_NUM_OLSON;h++){
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
	for(h=0;h<VEG_NUM_SAGE;h++){
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
	for(i=0;i<NREG;i++){
		fprintf(fp_glsum,"%lf ", rh_area[i]);
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_temp[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_prec[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_dswrf[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_rns[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_rnl[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ipar[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_apar[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_gpp[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_npp[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_nep[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_evpr[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_trsp[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_incp[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_rnof[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
		
	/* historical carbon isotopes: added by A.Ito (2009/07/15) */
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
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
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_gpp[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_gpp_d13c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_gpp_d14c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_er[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_er_d13c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_er_d14c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_f[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_f_d13c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_f_d14c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_c_d13c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_c_d14c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_r[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_r_d13c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_r_d14c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_l[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_l_d13c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_l_d14c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_h[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_h_d13c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ci_h_d14c[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	
	/* added by A.Ito (2009/11/15) */
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_hvst[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_luc[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ch4ox_curry[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ch4emit_wh_wet[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_ch4emit_wh_paddy[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_n2o_emit_ngas[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");for(h=0;h<HIST;h++){
		fprintf(fp_glsum,"%ld ", h+1900);
		for(i=0;i<NREG;i++){
			fprintf(fp_glsum,"%lf ", rh_n2o_emitagr_ngas[i][h]);
		}
		fprintf(fp_glsum,"\n");
	}
	fprintf(fp_glsum,"\n");
	
	fclose(fp_glsum);
}

/***********************************************************/
/* decadal average grid values */
void f_grid_av(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux
){
	short period;
	
	/* decadal periods **********/
	if(grid->climy>=1950 && grid->climy<1960){
		period = 0;
	}else if(grid->climy>=1990 && grid->climy<2000){
		period = 1;
	}else if(grid->climy>=2020 && grid->climy<2030){
		period = 2;
	}else if(grid->climy>=2050 && grid->climy<2060){
		period = 3;
	}else if(grid->climy>=2080 && grid->climy<2090){
		period = 4;
	}else{
		period = -1;
	}
	
	if(period>=0){
#if C13_GOUT==1				
		if(DF97==1){
			g_gpp13[period][grid->row][grid->col] = 
					d13c_addition((flux->plant).d13c_gpp[grid->m], (flux->plant).gpp_df97[grid->m] /10.0, 
					g_gpp13[period][grid->row][grid->col], g_gpp[period][grid->row][grid->col]);
		}else{
			g_gpp13[period][grid->row][grid->col] = 
					d13c_addition((flux->plant).d13c_gpp[grid->m], (flux->plant).gpp[grid->m] /10.0, 
					g_gpp13[period][grid->row][grid->col], g_gpp[period][grid->row][grid->col]);
		}
		g_er13[period][grid->row][grid->col] = d13c_addition(flux->d13c_er[grid->m], flux->er[grid->m] /10.0, 
												g_er13[period][grid->row][grid->col], g_er[period][grid->row][grid->col]);
#endif
		
#if C14_GOUT==1	
		if((flux->plant).gpp[grid->m]>0.0){
			g_gpp14[period][grid->row][grid->col] = (g_gpp14[period][grid->row][grid->col] * g_gpp[period][grid->row][grid->col] + 
												 (flux->plant).d14c_gpp[grid->m] * (flux->plant).gpp[grid->m] /10.0) / 
												(g_gpp[period][grid->row][grid->col] + (flux->plant).gpp[grid->m] /10.0);
		}
		if(flux->er[grid->m]>0.0){
			g_er14[period][grid->row][grid->col] = (g_er14[period][grid->row][grid->col] * g_er[period][grid->row][grid->col] + 
												flux->d14c_er[grid->m] * flux->er[grid->m] /10.0) / 
												(g_er[period][grid->row][grid->col] + flux->er[grid->m] /10.0);
		}
#endif
		
		g_tmp[period][grid->row][grid->col] += grid->tmp_2m[grid->m]* MDN[grid->m]/365.0 /10.0;
		g_prc[period][grid->row][grid->col] += grid->prate_sfc[grid->m] /10.0;
		g_swr[period][grid->row][grid->col] += grid->gl_rad[grid->m]* MDN[grid->m]/365.0 /10.0;
		if(DF97==1){
			g_gpp[period][grid->row][grid->col] += (flux->plant).gpp_df97[grid->m] /10.0;
		}else{
			g_gpp[period][grid->row][grid->col] += (flux->plant).gpp[grid->m] /10.0;
		}
		g_npp[period][grid->row][grid->col] += (flux->plant).npp[grid->m] /10.0;
		g_nep[period][grid->row][grid->col] += flux->nep[grid->m] /10.0;
		g_pmas[period][grid->row][grid->col] += ((mass->plant).mfol[grid->m] + (mass->plant).mstm[grid->m] + 
											(mass->plant).mrot[grid->m]) * MDN[grid->m] /365.0 /10.0;
		g_smas[period][grid->row][grid->col] += ((mass->soil).ltr_m[grid->m] + 
											(mass->soil).msl_m[grid->m])* MDN[grid->m]/365.0 /10.0;
		
		g_ch4e_cao[period][grid->row][grid->col] += ((flux->soil).ch4flux_paddy_cao[grid->m] + (flux->soil).ch4flux_wetland_cao[grid->m]) /10.0;
		g_ch4ep_cao[period][grid->row][grid->col] += (flux->soil).ch4flux_paddy_cao[grid->m] /10.0;
		g_ch4o_curry[period][grid->row][grid->col] += (flux->soil).ch4oxy_curry[grid->m] /10.0;
		g_n2oe[period][grid->row][grid->col] += (flux->soil).d_n2o_ngas[grid->m] /10.0;
		g_bbco2[period][grid->row][grid->col] += (flux->bb_co2_litter[grid->m]+flux->bb_co2_leaf[grid->m]+
												  flux->bb_co2_wood[grid->m]+flux->bb_co2_root[grid->m]) /10.0;
		g_isopr[period][grid->row][grid->col] += flux->voc_isopr_g97[grid->m] /10.0;
		g_sr[period][grid->row][grid->col] += ((flux->plant).rrm[grid->m] + (flux->plant).rrg[grid->m] + (flux->soil).hr[grid->m]) /10.0;
		g_er[period][grid->row][grid->col] += flux->er[grid->m] / 10.0;
		g_snh4[period][grid->row][grid->col] += (mass->soil).n_nh4*MDN[grid->m]/365.0 /10.0;
		g_sno3[period][grid->row][grid->col] += (mass->soil).n_no3*MDN[grid->m]/365.0 /10.0;

#if PHYS_GOUT==1
		g_lai[period][grid->row][grid->col] += (loct->c3ptn[grid->m]*(mass->c3).lai[grid->m]+
												loct->c4ptn[grid->m]*(mass->c4).lai[grid->m])*MDN[grid->m]/365.0 /10.0;
		g_parb[period][grid->row][grid->col] += grid->par_bp[grid->m]*MDN[grid->m]/365.0 /10.0;
		g_pard[period][grid->row][grid->col] += grid->par_dp[grid->m]*MDN[grid->m]/365.0 /10.0;
		g_apar[period][grid->row][grid->col] += loct->fapar_mono[grid->m] / 10.0;
		if(DF97==1){
			g_apar2[period][grid->row][grid->col] += loct->fapar_df[grid->m] / 10.0;
		}else{
			;
		}
		g_aet[period][grid->row][grid->col] += (loct->incep[grid->m]+loct->evpr[grid->m]+loct->trspr[grid->m]) /10.0;
		g_rof[period][grid->row][grid->col] += loct->ro2[grid->m] /10.0;
		g_rns[period][grid->row][grid->col] += loct->rad_net_short[grid->m]*MDN[grid->m]/365.0 /10.0;
		g_rnl[period][grid->row][grid->col] += loct->rad_net_long[grid->m]*MDN[grid->m]/365.0 /10.0;
		g_sw1[period][grid->row][grid->col] += loct->msw30[grid->m] /10.0;
		g_sw2[period][grid->row][grid->col] += loct->msww[grid->m] /10.0;
#endif
		
#if CH4_WH==1
		g_ch4ep_wh[period][grid->row][grid->col] += ((flux->soil).ch4_paddy_wh_plant[grid->m] + (flux->soil).ch4_paddy_wh_ebull[grid->m] + 
												(flux->soil).ch4_paddy_wh_diff[grid->m] + (flux->soil).ch4_paddy_wh_release[grid->m]) / 10.0;

		g_ch4ew_wh[period][grid->row][grid->col] += ((flux->soil).ch4_wetland_wh_plant[grid->m] + (flux->soil).ch4_wetland_wh_ebull[grid->m] + 
													 (flux->soil).ch4_wetland_wh_diff[grid->m] + (flux->soil).ch4_wetland_wh_release[grid->m]) / 10.0;
		if(period==1){
			/* monthly */
			gm_ch4ep_wh[grid->m][grid->row][grid->col] += ((flux->soil).ch4_paddy_wh_plant[grid->m] + (flux->soil).ch4_paddy_wh_ebull[grid->m] + 
													(flux->soil).ch4_paddy_wh_diff[grid->m] + (flux->soil).ch4_paddy_wh_release[grid->m]) / 10.0;
		}
#endif

#if C13_GOUT==1
		g_f13[period][grid->row][grid->col] += (mass->plant).d13c_mfol[grid->m]*MDN[grid->m]/365.0 /10.0; 
		g_c13[period][grid->row][grid->col] += (mass->plant).d13c_mstm[grid->m]*MDN[grid->m]/365.0 /10.0; 
		g_r13[period][grid->row][grid->col] += (mass->plant).d13c_mrot[grid->m]*MDN[grid->m]/365.0 /10.0; 
		g_l13[period][grid->row][grid->col] += (mass->soil).d13c_ltr_m[grid->m]*MDN[grid->m]/365.0 /10.0; 
		g_h13[period][grid->row][grid->col] += (mass->soil).d13c_msl_m[grid->m]*MDN[grid->m]/365.0 /10.0; 
#endif
		
#if C14_GOUT==1				
		g_f14[period][grid->row][grid->col] += (mass->plant).d14c_mfol[grid->m]*MDN[grid->m]/365.0 /10.0; 
		g_c14[period][grid->row][grid->col] += (mass->plant).d14c_mstm[grid->m]*MDN[grid->m]/365.0 /10.0; 
		g_r14[period][grid->row][grid->col] += (mass->plant).d14c_mrot[grid->m]*MDN[grid->m]/365.0 /10.0; 
		g_l14[period][grid->row][grid->col] += (mass->soil).d14c_ltr_m[grid->m]*MDN[grid->m]/365.0 /10.0; 
		g_h14[period][grid->row][grid->col] += (mass->soil).d14c_msl_m[grid->m]*MDN[grid->m]/365.0 /10.0; 
#endif				
	}
}
