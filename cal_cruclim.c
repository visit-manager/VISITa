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

/* HISTORICAL SIMULATION */
/* 1901-2000 (2002 / 2008 / 2009/ 2011) *******************************************/
void cal_cruclim(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux, 
	FILE *fp_o[OFILES]
){
	long f, g, dyr;
	double f_fert, total_hvst;
	extern double MDN[ASTEP];
	
    /* historical simulation */
	grid->phase = 1; 
	
	(echar->soil).rl = (echar->soil).rl0;
	(echar->soil).rh = (echar->soil).rh0;
	
	/* LOOP to dynamic stage *******************************************************/
	for(g=0; g<HIST_PD; g++){
		/* AD1901 - 2002 / 2008 / 2009 */
        /* ISIMIP: 1951-2099 */
		
		/* climate change ********************/
		grid->climy = PIVOT_CLIMY + g;
		if(grid->hist_exist == 1){
			set_hist_clim(grid);
		}
		
		/* land-use change ******************/
		f_cult_luc(grid);
		
		/* CO2 year ********************/
		grid->co2y = PIVOT_CO2Y + g; 
		/* sensitivity analysis: no CO2 rise */
		if(CC_CD==2){
			grid->co2y = PIVOT_CO2Y;
			/* PIVOT_CO2Y = 1901 (usual setting) */
		}
		
		/* historical change in fertilizer input: 2010/05/11 by A.Ito */
		if(grid->rank_nat==1){
			/* developing countries */
			f_fert = 2.0217112 / (1.0 + exp(0.049849599 * (2000.6575 - (double)grid->climy)))+0.0014929171;
		}else if(grid->rank_nat==2){
			/* developed countries */
			f_fert = 0.92939393 / (1.0 + exp(0.044112692 * (2000.0097 - (double)grid->climy)))+0.53533202;
		}
		
		/* seasonal (monthly) loop ********************************************/
		for(f=0;f<ASTEP;f++){
			grid->m = f;
			
			/* initialize N fluxes ************/
			ghg_flux_zero(f, flux);
			
			/* CO2 condition */
			f_co2_trend(grid);
			co2_in_canopy(grid, loct, mass, flux);
						
			/* environmental condition *******************/
			f_dyn_loct(grid, loct, mass, echar);
			
			/* vegetation processes *********************/
			f_biome_processes(grid, loct, echar, mass, flux);

			/* VOC emission *****************/
			f_voc_emit_guenther97(grid, loct, echar, mass, flux);
			/* Plant CH4 emission *****************/
			f_ch4_emit_veg(grid, loct, echar, mass, flux);

			f_plant_stand_budget(grid, loct,mass,flux);
			
			if(BACC==3){
				(flux->plant).lL[f] = flux->lL0[f];
			}

			(flux->soil).lL[f] = (flux->plant).lL[f];			
			(flux->soil).d13c_lL[f] = (flux->plant).d13c_lL[f];
			(flux->soil).d14c_lL[f] = (flux->plant).d14c_lL[f];
			
			/* soil processes *****************/
			soil_processes(grid, loct, &(echar->soil), &(mass->soil), &(flux->soil));
			flux->sr[f] = loct->c3ptn[f]*((flux->c3).rrm[f]+(flux->c3).rrg[f]) + 
							loct->c4ptn[f]*((flux->c4).rrm[f]+(flux->c4).rrg[f]) + 
							(flux->soil).hr[f];
			
			if(NECB_DOC==1){
				(mass->soil).msl -= (flux->soil).doc_boyer[f]/1000000.0;
				if((mass->soil).msl < 0.0){
					(mass->soil).msl = 0.0;
				}
			}

			/* fertilizaer input for croplands: revised by A.Ito (2009/06/04) */
			/* NH4:NO3 ratio is based on inventories */
			if((echar->soil).v_type == 1){
			   if(grid->veg_olson==29 || grid->veg_olson==30 || 
											 grid->veg_olson==31 || grid->veg_olson==32){
				   (flux->soil).n_fertin[f] = loct->n_frtlz_in * 1000.0 * f_fert;
				   (mass->soil).n_no3 += loct->n_frtlz_in * 0.2 * 1000.0 * f_fert;
				   (mass->soil).n_nh4 += loct->n_frtlz_in * 0.8 * 1000.0 * f_fert;
				}else{
				   (flux->soil).n_fertin[grid->m] = 0.0;
				}
			}
			
			/*
			if((echar->soil).v_type == 1){
				(mass->soil).n_no3 += grid->f_crop_con * loct->n_frtlz_in * 0.2 * 1000.0;
				(mass->soil).n_nh4 += grid->f_crop_con * loct->n_frtlz_in * 0.8 * 1000.0;
			} */
			
			if((echar->soil).v_type == 2){
				(flux->soil).n_fertin[f] = loct->n_frtlz_in * 1000.0 * f_fert;
				(mass->soil).n_no3 += loct->n_frtlz_in * 0.2 * 1000.0 * f_fert;
				(mass->soil).n_nh4 += loct->n_frtlz_in * 0.8 * 1000.0 * f_fert;
			}

			/* CH4 oxydation (uplands) ****************************/
			f_ch4oxy_ridgewell(grid, loct, flux);
			f_ch4oxy_casa(grid, loct, flux);
			f_ch4oxy_delgrosso(grid, loct, flux);
			f_ch4oxy_curry(grid, loct, flux);
			
			/* CH4 emission **************************************/
			/* Cao (paddy+wetlands) */
			f_ch4_emit_cao(grid, loct, flux);
			
			/* Walter & Heimann */
			/* wetlands */
			if(CH4_WH==1 && grid->f_wetland>0.0 && loct->v_type == 1){
				f_ch4_emit_walter(1, grid, loct, flux);
				f_ch4_emit_walter(2, grid, loct, flux);
			}else{
				loct->f_inund_wet_wh[f] = 0.0;
				(flux->soil).ch4_wetland_wh_plant[f] = 0.0;
				(flux->soil).ch4_wetland_wh_ebull[f] = 0.0;
				(flux->soil).ch4_wetland_wh_diff[f] = 0.0;
				(flux->soil).ch4_wetland_wh_release[f] = 0.0;
			}
			/* paddy fields */
			if(CH4_WH==1 && grid->f_paddy>0.0 && loct->v_type == 2){
				f_ch4_emit_walter(3, grid, loct, flux);
				f_ch4_emit_walter(4, grid, loct, flux);
			}else{
				loct->f_inund_pad_wh[f] = 0.0;
				(flux->soil).ch4_paddy_wh_plant[f] = 0.0;
				(flux->soil).ch4_paddy_wh_ebull[f] = 0.0;
				(flux->soil).ch4_paddy_wh_diff[f] = 0.0;
				(flux->soil).ch4_paddy_wh_release[f] = 0.0;
			}
			
			/* coupling carbon budget by CH4 */
			if(NECB_CH4 == 1){
					(mass->soil).msl += grid->f_upland * (flux->soil).ch4oxy_curry[f] * 0.00001
						- grid->f_paddy * ((flux->soil).ch4_paddy_wh_plant[f] + (flux->soil).ch4_paddy_wh_ebull[f] + 
									(flux->soil).ch4_paddy_wh_diff[f] + (flux->soil).ch4_paddy_wh_release[f]) * 0.00001
						- grid->f_wetland * ((flux->soil).ch4_wetland_wh_plant[f] + (flux->soil).ch4_wetland_wh_ebull[f] + 
									(flux->soil).ch4_wetland_wh_diff[f] + (flux->soil).ch4_wetland_wh_release[f]) * 0.00001;
				
				if((mass->soil).msl < 0.0){
					(mass->soil).msl = 0.0;
				}
			}
			
			/* N2O emission ***************************************/
			f_n2o_emit_ngas(grid, loct, mass, flux);
			f_n2o_emit_casa(grid, loct, mass, flux);
			
			/* ecosystem mass balance *****************************/	
			/* net ecosystem production */
			flux->nep[f] = (flux->plant).npp[f]-(flux->soil).hr[f];
			flux->er[f] = (flux->plant).ar[f] + (flux->soil).hr[f];
			/* total ecosystem carbon storage */
			mass->total[f] = (mass->c3).plant[f]*loct->c3ptn[f] + (mass->c4).plant[f]*loct->c4ptn[f] + (mass->soil).soil[f];
			/** net carbon balance taking crop harvest into account **/
			flux->ncb[f] = flux->nep[f] + (flux->plant).hvst[f];
			
			/* carbon isotope */
			f_cisotope_efflux(grid, loct, mass, flux);
			
			/* nitrogen budget */
			n_budget(grid, loct, mass, flux);
			
			/* average LAI: 2009/05/06 by A.Ito */
			if(grid->climy>=1990 && grid->climy<=1999){
				(mass->c3).lai0[f] += (mass->c3).lai[f]/10.0;
				(mass->c4).lai0[f] += (mass->c4).lai[f]/10.0;
				(mass->plant).lai0[f] += (mass->plant).lai[f]/10.0;
			}
			
			/* statistics *******************************************/
			if(g>=90 && g<=99){
				/* mean seasonal change *******/
				m_ch4ox1[f] += (flux->soil).ch4oxy_ridg[f] * grid->area *10000.0/1000.0 / 10.0;
				m_ch4ox2[f] += (flux->soil).ch4oxy_casa[f] * grid->area *10000.0/1000.0 / 10.0;
				m_ch4ox3[f] += (flux->soil).ch4oxy_delgrosso[f] * grid->area *10000.0/1000.0 / 10.0;
				
				m_gpp[f] += (flux->plant).gpp[f]/10.0 * grid->area;
				m_npp[f] += (flux->plant).npp[f]/10.0 * grid->area;
				m_nep[f] += flux->nep[f]/10.0 * grid->area;
				
				if(DF97==1){
					vo_gpp[grid->veg_olson] += (flux->plant).gpp_df97[f]/10.0 * grid->area;
				}else{
					vo_gpp[grid->veg_olson] += (flux->plant).gpp[f]/10.0 * grid->area;
				}
				
				/* mean biome budget *******/
				vo_npp[grid->veg_olson] += (flux->plant).npp[f]/10.0 * grid->area;
				vo_nep[grid->veg_olson] += flux->nep[f]/10.0 * grid->area;
				vo_lai[grid->veg_olson] += (mass->plant).lai[f]*MDN[f]/365.0/10.0 * grid->area;
				vo_fol[grid->veg_olson] += (mass->plant).mfol[f]*MDN[f]/365.0/10.0 * grid->area;
				vo_stm[grid->veg_olson] += (mass->plant).mstm[f]*MDN[f]/365.0/10.0 * grid->area;
				vo_rot[grid->veg_olson] += (mass->plant).mrot[f]*MDN[f]/365.0/10.0 * grid->area;
				vo_ltr[grid->veg_olson] += (mass->soil).ltr_m[f]*MDN[f]/365.0/10.0 * grid->area;
				vo_msl[grid->veg_olson] += (mass->soil).msl_m[f]*MDN[f]/365.0/10.0 * grid->area;
				
				if(DF97==1){
					vs_gpp[grid->veg_sage] += (flux->plant).gpp_df97[f]/10.0 * grid->area;
				}else{
					vs_gpp[grid->veg_sage] += (flux->plant).gpp[f]/10.0 * grid->area;
				}
				
				vs_npp[grid->veg_sage] += (flux->plant).npp[f]/10.0 * grid->area;
				vs_nep[grid->veg_sage] += flux->nep[f]/10.0 * grid->area;
				vs_lai[grid->veg_sage] += (mass->plant).lai[f]*MDN[f]/365.0/10.0 * grid->area;
				vs_fol[grid->veg_sage] += (mass->plant).mfol[f]*MDN[f]/365.0/10.0 * grid->area;
				vs_stm[grid->veg_sage] += (mass->plant).mstm[f]*MDN[f]/365.0/10.0 * grid->area;
				vs_rot[grid->veg_sage] += (mass->plant).mrot[f]*MDN[f]/365.0/10.0 * grid->area;
				vs_ltr[grid->veg_sage] += (mass->soil).ltr_m[f]*MDN[f]/365.0/10.0 * grid->area;
				vs_msl[grid->veg_sage] += (mass->soil).msl_m[f]*MDN[f]/365.0/10.0 * grid->area;
			}
			
			/******************/
			f_grid_av(grid, loct, echar, mass, flux);
		}
		
		/* empirical NPP models */
		npp_empirical(grid, loct, flux);
		
		/* biomass burning **************************/
		f_biomassburning(grid, loct, mass, flux);
		/* corrected: A.Ito and E.Kato (2009/08/16) */
		if(g>=90 && g<=99){
			for(f=0;f<ASTEP;f++){
				m_bioburn_co2[f] += (flux->bb_co2_litter[f]+flux->bb_co2_leaf[f]
									 +flux->bb_co2_wood[f]+flux->bb_co2_root[f]) * grid->area / 10.0;
				m_bioburn_co[f] += (flux->bb_co_litter[f]+flux->bb_co_leaf[f]
									+flux->bb_co_wood[f]+flux->bb_co_root[f]) * grid->area / 10.0;
				m_bioburn_ch4[f] += (flux->bb_ch4_litter[f]+flux->bb_ch4_leaf[f]
									 +flux->bb_ch4_wood[f]+flux->bb_ch4_root[f]) * grid->area / 10.0;
				m_bioburn_nmhc[f] += (flux->bb_nmhc_litter[f]+flux->bb_nmhc_leaf[f]
									  +flux->bb_nmhc_wood[f]+flux->bb_nmhc_root[f]) * grid->area / 10.0;
				m_bioburn_oc[f] += (flux->bb_oc_litter[f]+flux->bb_oc_leaf[f]
									+flux->bb_oc_wood[f]+flux->bb_oc_root[f]) * grid->area / 10.0;
				m_bioburn_bc[f] += (flux->bb_bc_litter[f]+flux->bb_bc_leaf[f]
									+flux->bb_bc_wood[f]+flux->bb_bc_root[f]) * grid->area / 10.0;
			}
		}
		if(g>=80 && g<=89){
			for(f=0;f<ASTEP;f++){
				m_ch4p_cao[f] += (flux->soil).ch4flux_paddy_cao[f] /10.0;
				m_ch4p_wh[f] += ((flux->soil).ch4_paddy_wh_plant[f] + (flux->soil).ch4_paddy_wh_ebull[f] + 
								 (flux->soil).ch4_paddy_wh_diff[f] + (flux->soil).ch4_paddy_wh_release[f]) / 10.0;
			}
		}

		/* erosion ****************************/
		f_erosion(grid, loct, echar, mass, flux);
		
		if(NECB_ERSN==1){
			(mass->soil).ltr -= flux->erod_carbon*0.20;
			if((mass->soil).ltr<0.0){
				(mass->soil).ltr = 0.0;
			}
		}

		/* land use change *********************/
		if(loct->v_type == 1){
			f_luc_emit(grid, mass, flux);
		}else{
			flux->lu_detr = 0.0;
			flux->lu_conv = 0.0;
			flux->lu_ten = 0.0;
			flux->lu_hund = 0.0;
		}
		
		/* wood harvest: 2010/10/15 by A.Ito ***************/
		total_hvst = 0.0;
		if((mass->c3).v_type == 1 && NECB_WHVST == 1){
			dyr = grid->climy - 1700;
			
			/* assumption for the period later than 2004: A.Ito (2010/11/11) */
			if(dyr>304){
				dyr = 304;	
			}
			
			total_hvst = grid->hvst_p1[dyr] + grid->hvst_p2[dyr] + grid->hvst_s1[dyr] 
						+ grid->hvst_s2[dyr] + grid->hvst_s3[dyr];
			
			total_hvst *= 1.0/1000.0 * 1.0/grid->area;
			
			if((mass->c3).stm > (total_hvst+1.0)){
				(mass->c3).stm -= total_hvst;
				flux->hvst_wood = total_hvst;
			}else{
				(mass->c3).stm = 1.0;
				flux->hvst_wood = 0.0; 
			}
			
			if((mass->c3).stm < 1.0){
				(mass->c3).stm = 1.0;
			}
		}else{
			flux->hvst_wood = 0.0;
		}
		
		/* net biome production (added by A.Ito: 2010/01/20) *************************/
		for(f=0;f<ASTEP;f++){
			flux->nbp[f] = flux->nep[f];
			
			if(NECB_LUC == 1){
				flux->nbp[f] -= (flux->lu_ten/12.0 + flux->lu_hund/12.0);
			}
			
			if(NECB_BB == 1){
				flux->nbp[f] -= (flux->bb_co2_litter[f] + flux->bb_co2_leaf[f] 
								 + flux->bb_co2_wood[f] + flux->bb_co2_root[f])/1000.0*12.0/44.0;
			}
		}
		
		/* history data */
		f_set_history_data(grid->climy - PIVOT_CLIMY +1, grid, loct, mass, flux);
		
		/* output */
		f_output_result(grid->climy, grid, loct, echar, mass, flux, fp_o); /* */
		
		/* setting previous land-use */
		grid->f_crop_p = grid->f_crop_con;
		grid->f_pasture_p = grid->f_pasture_con;
		
		if(grid->climy>=1950 && grid->climy<1960){
			g_ersn[0][grid->row][grid->col] += flux->erod_carbon /10.0;
			g_luc[0][grid->row][grid->col] += (flux->lu_conv + flux->lu_ten + flux->lu_hund) /10.0;
		}
		if(grid->climy>=1990 && grid->climy<2000){
			g_ersn[1][grid->row][grid->col] += flux->erod_carbon /10.0;
			g_luc[1][grid->row][grid->col] += (flux->lu_conv + flux->lu_ten + flux->lu_hund) /10.0;
		}
	}
	fprintf(fp_o[0],"\n");
	fprintf(fp_o[1],"\n");
	fprintf(fp_o[2],"\n");
	fprintf(fp_o[3],"\n");
	fprintf(fp_o[4],"\n");
	fprintf(fp_o[5],"\n");
	fprintf(fp_o[6],"\n");
}
