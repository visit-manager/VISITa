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

extern double m_ch4ox1[ASTEP], m_ch4ox2[ASTEP], m_ch4ox3[ASTEP];
extern double m_bioburn_co2[ASTEP], m_bioburn_ch4[ASTEP], m_bioburn_co[ASTEP];
extern double m_bioburn_nmhc[ASTEP], m_bioburn_oc[ASTEP], m_bioburn_bc[ASTEP];
extern double m_gpp[ASTEP], m_npp[ASTEP], m_nep[ASTEP];

extern double vo_area[34];
extern double vo_gpp[34], vo_npp[34], vo_nep[34];
extern double vo_lai[34], vo_fol[34], vo_stm[34], vo_rot[34], vo_ltr[34], vo_msl[34];
extern double vs_area[16];
extern double vs_gpp[16], vs_npp[16], vs_nep[16];
extern double vs_lai[16], vs_fol[16], vs_stm[16], vs_rot[16], vs_ltr[16], vs_msl[16];
extern short DF97;

extern float g_tmp[5][360][720];
extern float g_prc[5][360][720];
extern float g_swr[5][360][720];
extern float g_gpp[5][360][720];
extern float g_npp[5][360][720];
extern float g_nep[5][360][720];
extern float g_pmas[5][360][720];
extern float g_smas[5][360][720];
extern float g_ch4e[5][360][720];
extern float g_ch4o[5][360][720];
extern float g_n2oe[5][360][720];
extern float g_bbco2[5][360][720];
extern float g_ersn[5][360][720];
extern float g_isopr[5][360][720];
extern float g_sr[5][360][720];
extern float g_luc[5][360][720];

/* 1901-2000 (2002 / 2008) *****************************************************/
void cal_cruclim(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux, 
	FILE *fp_o[OFILES]
){
	long f, g;
	extern double MDN[ASTEP];
	
	grid->phase = 1; /* history */
	
	(echar->soil).rl = (echar->soil).rl0;
	(echar->soil).rh = (echar->soil).rh0;
	
	/* LOOP to dynamic stage ***********************************************/
	for(g=0; g<CRU_PD; g++){
		/* AD1901 - 2002 / 2008 */
		
		/* climate change ********************/
		grid->climy = PIVOT_CLIMY + g;
		if(grid->cru_exist == 1){
			set_cru_clim(grid);
		}
		
		/* land-use change ******************/
		f_cult_luc(grid);
		
		/* CO2 year ********************/
		grid->CO2y = PIVOT_CO2Y + g; 
		
		/* monthly roop ****************************************************/
		for(f=0;f<ASTEP;f++){
			grid->m = f;
			
			/* initialize N fluxes ************/
			n_flux_zero(f, flux);
			
			/* CO2 condition */
			cd_trend(grid);
			co2_in_canopy(grid, loct, mass, flux);
						
			/* environmental condition *******************/
			dynmcL(grid, loct, mass, echar);
			
			/***** vegetation processes *****/
			biome_processes(grid, loct, echar, mass, flux);

			/* VOC emission *****************/
			f_voc_emit_guenther97(grid, loct, echar, mass, flux);
			/* Plant CH4 emission *****************/
			f_ch4_emit_veg(grid, loct, echar, mass, flux);

			plant_stand(grid, loct,mass,flux);
			
			if(BACC==3){
				(flux->plant).lL[f] = flux->lL0[f];
			}

			(flux->soil).lL[f] = (flux->plant).lL[f];			
			(flux->soil).d13c_lL[f] = (flux->plant).d13c_lL[f];
			
			/***** soil processes *****/
			soil_processes(grid, loct, &(echar->soil), &(mass->soil), &(flux->soil));

			/* fertilizaer input */
			if((echar->soil).v_type == 3){
				(mass->soil).n_no3 += loct->n_frtlz_in * 0.5;
				(mass->soil).n_nh4 += loct->n_frtlz_in * 0.5;
			}

			/* CH4 oxydation (uplands) **************/
			f_ch4oxy_ridgewell(grid, loct, flux);
			f_ch4oxy_casa(grid, loct, flux);
			f_ch4oxy_delgrosso(grid, loct, flux);
			f_ch4oxy_curry(grid, loct, flux);
			
			/* CH4 emission (wetlands) */
			f_ch4_emit_cao(grid, loct, flux);
			
			/* N2O emission */
			f_n2o_emit_ngas(grid, loct, mass, flux);
			f_n2o_emit_casa(grid, loct, mass, flux);
						
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
			
			if(g>=90 && g<=99){
				/* mean seasonal change *******/
				m_ch4ox1[f] += (flux->soil).ch4oxy_ridg[f] * grid->area *10000.0/1000.0 / 10.0;
				m_ch4ox2[f] += (flux->soil).ch4oxy_casa[f] * grid->area *10000.0/1000.0 / 10.0;
				m_ch4ox3[f] += (flux->soil).ch4oxy_delgrosso[f] * grid->area *10000.0/1000.0 / 10.0;
				m_bioburn_co2[f] += (flux->bb_co2_litter[f]+flux->bb_co2_leaf[f]+flux->bb_co2_wood[f]+flux->bb_co2_root[f]) * grid->area / 10.0;
				m_bioburn_co[f] += (flux->bb_co_litter[f]+flux->bb_co_leaf[f]+flux->bb_co_wood[f]+flux->bb_co_root[f]) * grid->area / 10.0;
				m_bioburn_ch4[f] += (flux->bb_ch4_litter[f]+flux->bb_ch4_leaf[f]+flux->bb_ch4_wood[f]+flux->bb_ch4_root[f]) * grid->area / 10.0;
				m_bioburn_nmhc[f] += (flux->bb_nmhc_litter[f]+flux->bb_nmhc_leaf[f]+flux->bb_nmhc_wood[f]+flux->bb_nmhc_root[f]) * grid->area / 10.0;
				m_bioburn_oc[f] += (flux->bb_oc_litter[f]+flux->bb_oc_leaf[f]+flux->bb_oc_wood[f]+flux->bb_oc_root[f]) * grid->area / 10.0;
				m_bioburn_bc[f] += (flux->bb_bc_litter[f]+flux->bb_bc_leaf[f]+flux->bb_bc_wood[f]+flux->bb_bc_root[f]) * grid->area / 10.0;
				
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
			
			if(grid->climy>=1950 && grid->climy<1960){
				g_tmp[0][grid->row][grid->col] += grid->tmp_2m[grid->m]* MDN[grid->m]/365.0 /10.0;
				g_prc[0][grid->row][grid->col] += grid->prate_sfc[f] /10.0;
				g_swr[0][grid->row][grid->col] += grid->gl_rad[grid->m]* MDN[grid->m]/365.0 /10.0;
				if(DF97==1){
					g_gpp[0][grid->row][grid->col] += (flux->plant).gpp_df97[grid->m] /10.0;
				}else{
					g_gpp[0][grid->row][grid->col] += (flux->plant).gpp[grid->m] /10.0;
				}
				g_npp[0][grid->row][grid->col] += (flux->plant).npp[grid->m] /10.0;
				g_nep[0][grid->row][grid->col] += flux->nep[grid->m] /10.0;
				g_pmas[0][grid->row][grid->col] += ((mass->plant).mfol[grid->m] + (mass->plant).mstm[grid->m] + 
													(mass->plant).mrot[grid->m]) * MDN[grid->m] /365.0 /10.0;
				g_smas[0][grid->row][grid->col] += ((mass->soil).ltr_m[grid->m] + 
													(mass->soil).msl_m[grid->m])* MDN[grid->m]/365.0 /10.0;
				g_ch4e[0][grid->row][grid->col] += ((flux->soil).ch4flux_paddy_cao[grid->m] + (flux->soil).ch4flux_wetland_cao[grid->m]) /10.0;
				g_ch4o[0][grid->row][grid->col] += (flux->soil).ch4oxy_curry[grid->m] /10.0;
				g_n2oe[0][grid->row][grid->col] += (flux->soil).d_n2o_ngas[grid->m] /10.0;
				g_bbco2[0][grid->row][grid->col] += (flux->bb_co2_litter[grid->m]+flux->bb_co2_leaf[grid->m]+flux->bb_co2_wood[grid->m]+flux->bb_co2_root[grid->m]) /10.0;
				g_isopr[0][grid->row][grid->col] += flux->voc_isopr_g97[grid->m] /10.0;
				g_sr[0][grid->row][grid->col] += ((flux->plant).rrm[grid->m] + (flux->plant).rrg[grid->m] + (flux->soil).rS[grid->m]) /10.0;
			}
			if(grid->climy>=1990 && grid->climy<2000){
				g_tmp[1][grid->row][grid->col] += grid->tmp_2m[grid->m]* MDN[grid->m]/365.0 /10.0;
				g_prc[1][grid->row][grid->col] += grid->prate_sfc[f] /10.0;
				g_swr[1][grid->row][grid->col] += grid->gl_rad[grid->m]* MDN[grid->m]/365.0 /10.0;
				if(DF97==1){
					g_gpp[1][grid->row][grid->col] += (flux->plant).gpp_df97[grid->m] /10.0;
				}else{
					g_gpp[1][grid->row][grid->col] += (flux->plant).gpp[grid->m] /10.0;
				}
				g_npp[1][grid->row][grid->col] += (flux->plant).npp[grid->m] /10.0;
				g_nep[1][grid->row][grid->col] += flux->nep[grid->m] /10.0;
				g_pmas[1][grid->row][grid->col] += ((mass->plant).mfol[grid->m] + (mass->plant).mstm[grid->m] + 
													(mass->plant).mrot[grid->m]) * MDN[grid->m] /365.0 /10.0;
				g_smas[1][grid->row][grid->col] += ((mass->soil).ltr_m[grid->m] + 
													(mass->soil).msl_m[grid->m])* MDN[grid->m]/365.0 /10.0;
				g_ch4e[1][grid->row][grid->col] += ((flux->soil).ch4flux_paddy_cao[grid->m] + (flux->soil).ch4flux_wetland_cao[grid->m]) /10.0;
				g_ch4o[1][grid->row][grid->col] += (flux->soil).ch4oxy_curry[grid->m] /10.0;
				g_n2oe[1][grid->row][grid->col] += (flux->soil).d_n2o_ngas[grid->m] /10.0;
				g_bbco2[1][grid->row][grid->col] += (flux->bb_co2_litter[grid->m]+flux->bb_co2_leaf[grid->m]+flux->bb_co2_wood[grid->m]+flux->bb_co2_root[grid->m]) /10.0;
				g_isopr[1][grid->row][grid->col] += flux->voc_isopr_g97[grid->m] /10.0;
				g_sr[1][grid->row][grid->col] += ((flux->plant).rrm[grid->m] + (flux->plant).rrg[grid->m] + (flux->soil).rS[grid->m]) /10.0;
			}
		}
		/* empirical NPP models */
		npp_empirical(grid, loct, flux);
		
		/* biomass burning */
		f_biomassburning(grid, loct, mass, flux);

		/* erosion */
		f_erosion(grid, loct, echar, mass, flux);
		
		if(ERSN_CC==1){
			(mass->soil).ltr -= flux->erod_carbon*0.25;
			if((mass->soil).ltr<0.0){
				(mass->soil).ltr = 0.0;
			}
		}

		/* land use change */
		f_luc_emit(grid, mass, flux);
		
		/* history data */
		f_set_history_data(grid->climy - PIVOT_CLIMY +1, grid, loct, mass, flux);
		
		/* output */
		/* f_output_result(grid->climy, grid, loct, echar, mass, flux, fp_o); */
		
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
