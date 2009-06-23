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

extern long GCM, CO2S, GCM_R, GCM_C;
extern double MDN[12];
extern short DF97;
extern short TEMP_GC;

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

/* time-series with GCM climate scenarios : 2001-2099 ************************/
void cal_gcmclim2(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux, 
	FILE *fp_o[OFILES]
){
	long f, g, simyr;
	double rl_a;
	
	/* phase: prediction */
	grid->phase = 2; 
	
	/* set long-term average climate */
	initC(grid); 
	
	if(TEMP_GC==1 || TEMP_GC==2){
		simyr = 400;
	}else if(TEMP_GC==3 || TEMP_GC==4){
		simyr = 200;
	}else{
		simyr = GCM_PD;
	}
	
	/* LOOP to dynamic stage ***************************************/
	for(g=0;g<simyr;g++){ /*** AD 2001-2100 ***/
	
		/* climate change ********************/
		grid->climy = 2001+ g;
		if(GCM!=0){			
			set_gcm_clim(grid);
		}
		
		/* land-use change */
		if(TEMP_GC != 0){
			;
		}else{
			f_cult_luc(grid);
		}
		
		/* CO2 change ********************/
		if(CO2S==0){
			grid->CO2y = 2001; 
		}else if(CO2S==7){
			grid->CO2y = 2081;
		}else{
			grid->CO2y = 2001 + g; 
		}
		if(TEMP_GC != 0){
			grid->CO2y = 2001;
		}
				
		/* monthly loop *************************************/
		for(f=0;f<ASTEP;f++){
			grid->m = f;
									
			/* initialize N fluxes ************/
			n_flux_zero(f, flux);
			
			/* atmospheric CO2 */
			cd_trend(grid);
			
			/* simplified atm. CO2 change: added by A.Ito (2009/06/18) */
			if(TEMP_GC == 1){
				grid->bCO2[f] = grid->bCO2[f] * exp((-1.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 2){
				grid->bCO2[f] = grid->bCO2[f] * exp((1.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 3){
				grid->bCO2[f] = grid->bCO2[f] * exp((-2.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 4){
				grid->bCO2[f] = grid->bCO2[f] * exp((2.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 5){
				grid->bCO2[f] = grid->bCO2[f] * exp((-3.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 6){
				grid->bCO2[f] = grid->bCO2[f] * exp((3.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 7){
				grid->bCO2[f] = grid->bCO2[f] * exp((-4.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 8){
				grid->bCO2[f] = grid->bCO2[f] * exp((4.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 9){
				grid->bCO2[f] = grid->bCO2[f] * exp((-5.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 10){
				grid->bCO2[f] = grid->bCO2[f] * exp((5.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 11){
				grid->bCO2[f] = grid->bCO2[f] * exp((-6.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 12){
				grid->bCO2[f] = grid->bCO2[f] * exp((6.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}

			co2_in_canopy(grid, loct, mass, flux);
						
			/* environmental condition *******************/
			dynmcL(grid, loct, mass, echar);
			
			/* vegetation processes ***********************/
			biome_processes(grid, loct, echar, mass, flux);

			/* VOC emission *****************/
			f_voc_emit_guenther97(grid, loct, echar, mass, flux);
			/* Plant CH4 emission *****************/
			f_ch4_emit_veg(grid, loct, echar, mass, flux);
			
			/* aggregate plant mass and fluxes */
			plant_stand(grid, loct, mass, flux);
			
			if(BACC==3){
				(flux->plant).lL[f] = flux->lL0[f];
			}
			if(BACC==4){
				rl_a = (echar->soil).rl0*(1.0 - 0.001*(double)(g+1));
				if((mass->soil).ltr+(flux->plant).lL[f]){
					(echar->soil).rl = ((echar->soil).rl*(mass->soil).ltr + 
							rl_a*(flux->plant).lL[f])/((mass->soil).ltr+(flux->plant).lL[f]);
				}else{
					(echar->soil).rl = rl_a;
				}
			}

			(flux->soil).lL[f] = (flux->plant).lL[f];			
			(flux->soil).d13c_lL[f] = (flux->plant).d13c_lL[f];
			
			/* soil processes *****************************/
			soil_processes(grid, loct, &(echar->soil), &(mass->soil), &(flux->soil));
 
			/* fertilizaer input for croplands: revised by A.Ito (2009/06/04) */
			if((echar->soil).v_type == 1 && (grid->veg_olson!=29 || grid->veg_olson!=30 || 
					grid->veg_olson!=31 || grid->veg_olson!=32)){
				(mass->soil).n_no3 += loct->n_frtlz_in * 0.1 * 1000.0;
				(mass->soil).n_nh4 += loct->n_frtlz_in * 0.9 * 1000.0;
			}
			if((echar->soil).v_type == 3){
				(mass->soil).n_no3 += loct->n_frtlz_in * 0.1 * 1000.0;
				(mass->soil).n_nh4 += loct->n_frtlz_in * 0.9 * 1000.0;
			}

			/* CH4 oxydation **************/
			f_ch4oxy_ridgewell(grid, loct, flux);
			f_ch4oxy_casa(grid, loct, flux);
			f_ch4oxy_delgrosso(grid, loct, flux);
			f_ch4oxy_curry(grid, loct, flux);

			/* CH4 emission (wetlands) */
			f_ch4_emit_cao(grid, loct, flux);

			/* N2O emission */
			f_n2o_emit_ngas(grid, loct, mass, flux);
			f_n2o_emit_casa(grid, loct, mass, flux);
			
			/* ecosystem mass balance *************/	
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
			
			if(grid->climy>=2020 && grid->climy<2030){
				g_tmp[2][grid->row][grid->col] += grid->tmp_2m[grid->m]* MDN[grid->m]/365.0 /10.0;
				g_prc[2][grid->row][grid->col] += grid->prate_sfc[f] /10.0;
				g_swr[2][grid->row][grid->col] += grid->gl_rad[grid->m]* MDN[grid->m]/365.0 /10.0;
				if(DF97==1){
					g_gpp[2][grid->row][grid->col] += (flux->plant).gpp_df97[grid->m] /10.0;
				}else{
					g_gpp[2][grid->row][grid->col] += (flux->plant).gpp[grid->m] /10.0;
				}
				g_npp[2][grid->row][grid->col] += (flux->plant).npp[grid->m] /10.0;
				g_nep[2][grid->row][grid->col] += flux->nep[grid->m] /10.0;
				g_pmas[2][grid->row][grid->col] += ((mass->plant).mfol[grid->m] + (mass->plant).mstm[grid->m] + 
												   (mass->plant).mrot[grid->m]) * MDN[grid->m] /365.0 /10.0;
				g_smas[2][grid->row][grid->col] += ((mass->soil).ltr_m[grid->m] + 
												   (mass->soil).msl_m[grid->m])* MDN[grid->m]/365.0 /10.0;
				g_ch4e[2][grid->row][grid->col] += ((flux->soil).ch4flux_paddy_cao[grid->m] + (flux->soil).ch4flux_wetland_cao[grid->m]) /10.0;
				g_ch4o[2][grid->row][grid->col] += (flux->soil).ch4oxy_curry[grid->m] /10.0;
				g_n2oe[2][grid->row][grid->col] += (flux->soil).d_n2o_ngas[grid->m] /10.0;
				g_bbco2[2][grid->row][grid->col] += (flux->bb_co2_litter[grid->m]+flux->bb_co2_leaf[grid->m]+flux->bb_co2_wood[grid->m]+flux->bb_co2_root[grid->m]) /10.0;
				g_isopr[2][grid->row][grid->col] += flux->voc_isopr_g97[grid->m] /10.0;
				g_sr[2][grid->row][grid->col] += ((flux->plant).rrm[grid->m] + (flux->plant).rrg[grid->m] + (flux->soil).rS[grid->m]) /10.0;
			}
			if(grid->climy>=2050 && grid->climy<2060){
				g_tmp[3][grid->row][grid->col] += grid->tmp_2m[grid->m]* MDN[grid->m]/365.0 /10.0;
				g_prc[3][grid->row][grid->col] += grid->prate_sfc[f] /10.0;
				g_swr[3][grid->row][grid->col] += grid->gl_rad[grid->m]* MDN[grid->m]/365.0 /10.0;
				if(DF97==1){
					g_gpp[3][grid->row][grid->col] += (flux->plant).gpp_df97[grid->m] /10.0;
				}else{
					g_gpp[3][grid->row][grid->col] += (flux->plant).gpp[grid->m] /10.0;
				}
				g_npp[3][grid->row][grid->col] += (flux->plant).npp[grid->m] /10.0;
				g_nep[3][grid->row][grid->col] += flux->nep[grid->m] /10.0;
				g_pmas[3][grid->row][grid->col] += ((mass->plant).mfol[grid->m] + (mass->plant).mstm[grid->m] + 
													(mass->plant).mrot[grid->m]) * MDN[grid->m] /365.0 /10.0;
				g_smas[3][grid->row][grid->col] += ((mass->soil).ltr_m[grid->m] + 
													(mass->soil).msl_m[grid->m])* MDN[grid->m]/365.0 /10.0;
				g_ch4e[3][grid->row][grid->col] += ((flux->soil).ch4flux_paddy_cao[grid->m] + (flux->soil).ch4flux_wetland_cao[grid->m]) /10.0;
				g_ch4o[3][grid->row][grid->col] += (flux->soil).ch4oxy_curry[grid->m] /10.0;
				g_n2oe[3][grid->row][grid->col] += (flux->soil).d_n2o_ngas[grid->m] /10.0;
				g_bbco2[3][grid->row][grid->col] += (flux->bb_co2_litter[grid->m]+flux->bb_co2_leaf[grid->m]+flux->bb_co2_wood[grid->m]+flux->bb_co2_root[grid->m]) /10.0;
				g_isopr[3][grid->row][grid->col] += flux->voc_isopr_g97[grid->m] /10.0;
				g_sr[3][grid->row][grid->col] += ((flux->plant).rrm[grid->m] + (flux->plant).rrg[grid->m] + (flux->soil).rS[grid->m]) /10.0;
			}
			if(grid->climy>=2080 && grid->climy<2090){
				g_tmp[4][grid->row][grid->col] += grid->tmp_2m[grid->m]* MDN[grid->m]/365.0 /10.0;
				g_prc[4][grid->row][grid->col] += grid->prate_sfc[f] /10.0;
				g_swr[4][grid->row][grid->col] += grid->gl_rad[grid->m]* MDN[grid->m]/365.0 /10.0;
				if(DF97==1){
					g_gpp[4][grid->row][grid->col] += (flux->plant).gpp_df97[grid->m] /10.0;
				}else{
					g_gpp[4][grid->row][grid->col] += (flux->plant).gpp[grid->m] /10.0;
				}
				g_npp[4][grid->row][grid->col] += (flux->plant).npp[grid->m] /10.0;
				g_nep[4][grid->row][grid->col] += flux->nep[grid->m] /10.0;
				g_pmas[4][grid->row][grid->col] += ((mass->plant).mfol[grid->m] + (mass->plant).mstm[grid->m] + 
													(mass->plant).mrot[grid->m]) * MDN[grid->m] /365.0 /10.0;
				g_smas[4][grid->row][grid->col] += ((mass->soil).ltr_m[grid->m] + 
													(mass->soil).msl_m[grid->m])* MDN[grid->m]/365.0 /10.0;
				g_ch4e[4][grid->row][grid->col] += ((flux->soil).ch4flux_paddy_cao[grid->m] + (flux->soil).ch4flux_wetland_cao[grid->m]) /10.0;
				g_ch4o[4][grid->row][grid->col] += (flux->soil).ch4oxy_curry[grid->m] /10.0;
				g_n2oe[4][grid->row][grid->col] += (flux->soil).d_n2o_ngas[grid->m] /10.0;
				g_bbco2[4][grid->row][grid->col] += (flux->bb_co2_litter[grid->m]+flux->bb_co2_leaf[grid->m]+flux->bb_co2_wood[grid->m]+flux->bb_co2_root[grid->m]) /10.0;
				g_isopr[4][grid->row][grid->col] += flux->voc_isopr_g97[grid->m] /10.0;
				g_sr[4][grid->row][grid->col] += ((flux->plant).rrm[grid->m] + (flux->plant).rrg[grid->m] + (flux->soil).rS[grid->m]) /10.0;
			}
		}
		/* end of monthly loop ************************/
		
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
		f_output_result(grid->climy, grid, loct, echar, mass, flux, fp_o);  /*  */
		
		grid->f_crop_p = grid->f_crop_con;
		grid->f_pasture_p = grid->f_pasture_con;

		if(grid->climy>=2020 && grid->climy<2030){
			g_ersn[2][grid->row][grid->col] += flux->erod_carbon /10.0;
			g_luc[2][grid->row][grid->col] += (flux->lu_conv + flux->lu_ten + flux->lu_hund) /10.0;
		}
		if(grid->climy>=2050 && grid->climy<2060){
			g_ersn[3][grid->row][grid->col] += flux->erod_carbon /10.0;
			g_luc[3][grid->row][grid->col] += (flux->lu_conv + flux->lu_ten + flux->lu_hund) /10.0;
		}
		if(grid->climy>=2080 && grid->climy<2090){
			g_ersn[4][grid->row][grid->col] += flux->erod_carbon /10.0;
			g_luc[4][grid->row][grid->col] += (flux->lu_conv + flux->lu_ten + flux->lu_hund) /10.0;
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
