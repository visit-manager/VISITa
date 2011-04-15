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

extern short TEMP_GC;

/* time-series with GCM climate scenarios : 2001-2099 ************************/
void cal_gcmclim2(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux, 
	FILE *fp_o[OFILES]
){
	long f, g, simyr, dyr;
	double rl_a, f_fert, total_hvst;
	
	/* phase: prediction */
	grid->phase = 2; 
	
	/* set long-term average climate */
	f_init_clim(grid); 
	
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
			grid->co2y = 2001; 
		}else if(CO2S==7){
			grid->co2y = 2081;
		}else{
			grid->co2y = 2001 + g; 
		}
		if(TEMP_GC != 0){
			grid->co2y = 2001;
		}
		
		/* historical change in fertilizer input: 2010/05/11 by A.Ito */
		if(grid->rank_nat==1){
			/* developing countries */
			f_fert = 2.0217112 / (1.0 + exp(0.049849599 * (2000.6575 - (double)grid->climy)))+0.0014929171;
		}else if(grid->rank_nat==2){
			/* developed countries */
			f_fert = 0.92939393 / (1.0 + exp(0.044112692 * (2000.0097 - (double)grid->climy)))+0.53533202;
		}
				
		/* monthly loop *************************************/
		for(f=0;f<ASTEP;f++){
			grid->m = f;
									
			/* initialize N fluxes */
			ghg_flux_zero(f, flux);
			
			/* atmospheric CO2 */
			f_co2_trend(grid);
			
			/* simplified atm. CO2 change: added by A.Ito (2009/06/18) */
			if(TEMP_GC == 1){
				grid->bco2[f] = grid->bco2[f] * exp((-1.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 2){
				grid->bco2[f] = grid->bco2[f] * exp((1.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 3){
				grid->bco2[f] = grid->bco2[f] * exp((-2.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 4){
				grid->bco2[f] = grid->bco2[f] * exp((2.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 5){
				grid->bco2[f] = grid->bco2[f] * exp((-3.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 6){
				grid->bco2[f] = grid->bco2[f] * exp((3.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 7){
				grid->bco2[f] = grid->bco2[f] * exp((-4.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 8){
				grid->bco2[f] = grid->bco2[f] * exp((4.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 9){
				grid->bco2[f] = grid->bco2[f] * exp((-5.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 10){
				grid->bco2[f] = grid->bco2[f] * exp((5.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 11){
				grid->bco2[f] = grid->bco2[f] * exp((-6.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}else if(TEMP_GC == 12){
				grid->bco2[f] = grid->bco2[f] * exp((6.0/100.0 * (double)(grid->climy-2000)) / 6.0);
			}
			
			/* ambient CO2 in canopy */
			co2_in_canopy(grid, loct, mass, flux);
						
			/* environmental condition *******************/
			f_dyn_loct(grid, loct, mass, echar);
			
			/* vegetation processes ***********************/
			f_biome_processes(grid, loct, echar, mass, flux);

			/* VOC emission *****************/
			f_voc_emit_guenther97(grid, loct, echar, mass, flux);
			/* Plant CH4 emission *****************/
			f_ch4_emit_veg(grid, loct, echar, mass, flux);
			
			/* aggregate plant mass and fluxes */
			f_plant_stand_budget(grid, loct, mass, flux);
			
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
			(flux->soil).d14c_lL[f] = (flux->plant).d14c_lL[f];
			
			/* soil processes *****************************/
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
					(flux->soil).n_fertin[f] = 0.0;
				}
			}
			if((echar->soil).v_type == 2){
				(flux->soil).n_fertin[f] = loct->n_frtlz_in * 1000.0 * f_fert;
				(mass->soil).n_no3 += loct->n_frtlz_in * 0.2 * 1000.0 * f_fert;
				(mass->soil).n_nh4 += loct->n_frtlz_in * 0.8 * 1000.0 * f_fert;
			}

			/* CH4 oxydation **************/
			f_ch4oxy_ridgewell(grid, loct, flux);
			f_ch4oxy_casa(grid, loct, flux);
			f_ch4oxy_delgrosso(grid, loct, flux);
			f_ch4oxy_curry(grid, loct, flux);

			/* CH4 emission ***************/
			/* Cao (paddy+wetlands) */
			f_ch4_emit_cao(grid, loct, flux);
			
			/* Walter & Heimann (paddy) */
			if(CH4_WH==1 && grid->f_wetland>0.0 && loct->v_type == 1){
				f_ch4_emit_walter(1, grid, loct, flux);
				f_ch4_emit_walter(2, grid, loct, flux);
			}else{
				(flux->soil).ch4_wetland_wh_plant[f] = 0.0;
				(flux->soil).ch4_wetland_wh_ebull[f] = 0.0;
				(flux->soil).ch4_wetland_wh_diff[f] = 0.0;
				(flux->soil).ch4_wetland_wh_release[f] = 0.0;
			}
			if(CH4_WH==1 && grid->f_paddy>0.0 && loct->v_type == 2){
				f_ch4_emit_walter(3, grid, loct, flux);
				f_ch4_emit_walter(4, grid, loct, flux);
			}else{
				(flux->soil).ch4_paddy_wh_plant[f] = 0.0;
				(flux->soil).ch4_paddy_wh_ebull[f] = 0.0;
				(flux->soil).ch4_paddy_wh_diff[f] = 0.0;
				(flux->soil).ch4_paddy_wh_release[f] = 0.0;
			}
			if(NECB_CH4==1){
				(mass->soil).msl += grid->f_upland * (flux->soil).ch4oxy_curry[f] * 0.00001
					- grid->f_paddy * ((flux->soil).ch4_paddy_wh_plant[f] + (flux->soil).ch4_paddy_wh_ebull[f] + 
								   (flux->soil).ch4_paddy_wh_diff[f] + (flux->soil).ch4_paddy_wh_release[f]) * 0.00001
					- grid->f_wetland * ((flux->soil).ch4_wetland_wh_plant[f] + (flux->soil).ch4_wetland_wh_ebull[f] + 
									 (flux->soil).ch4_wetland_wh_diff[f] + (flux->soil).ch4_wetland_wh_release[f]) * 0.00001;
				
				if((mass->soil).msl < 0.0){
					(mass->soil).msl = 0.0;
				}
			}

			/* N2O emission */
			f_n2o_emit_ngas(grid, loct, mass, flux);
			f_n2o_emit_casa(grid, loct, mass, flux);
			
			/* ecosystem mass balance *************/	
			/* net ecosystem production */
			flux->nep[f] = (flux->plant).npp[f]-(flux->soil).hr[f];
			flux->er[f] = (flux->plant).ar[f] + (flux->soil).hr[f];
			/* total ecosystem carbon storage */
			mass->total[f] = (mass->c3).plant[f]*loct->c3ptn[f]+(mass->c4).plant[f]*loct->c4ptn[f]+(mass->soil).soil[f];
			/** net carbon balance taking crop harvest into account **/
			flux->ncb[f] = flux->nep[f]+(flux->plant).hvst[f];
			
			/* carbon isotope */
			f_cisotope_efflux(grid, loct, mass, flux);
			
			/* nitrogen budget */
			n_budget(grid, loct, mass, flux);
			
			/************************/
			f_grid_av(grid, loct, echar, mass, flux);
		}
		/* end of monthly loop ************************/
		
		/* empirical NPP models ********************/
		npp_empirical(grid, loct, flux);
				
		/* biomass burning */
		f_biomassburning(grid, loct, mass, flux);

		/* erosion */
		f_erosion(grid, loct, echar, mass, flux);
				
		if(NECB_ERSN==1){
			(mass->soil).ltr -= flux->erod_carbon*0.20;
			if((mass->soil).ltr < 0.0){
				(mass->soil).ltr = 0.0;
			}
		}

		/* land use change ************************/
		if(loct->v_type == 1){
			f_luc_emit(grid, mass, flux);
		}else{
			flux->lu_detr = 0.0;
			flux->lu_conv = 0.0;
			flux->lu_ten = 0.0;
			flux->lu_hund = 0.0;
		}
		
		/* wood harvest: 2010/10/15 by A.Ito ***********/
		total_hvst = 0.0;
		if((mass->c3).v_type == 1 && NECB_WHVST == 1){
			dyr = grid->climy - 1700;
			if(g>304){
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
		
		/* net biome production (added by A.Ito: 2010/01/20) */
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
