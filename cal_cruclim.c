/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed  in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* HISTORICAL SIMULATION */
/* 1901-2000 (2002 / 2008 / 2009/ 2011 / 2012) ***********************************/
void cal_historical(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux, 
	FILE *fp_o[OFILEN]
){
	long f, g, dyr, y_nin;
	double f_fert, fweight, total_hvst, f_nat, iweight, iweight3, avc3, prm_ensen;
    double mass_luc, icrop, aa, bb, base_nin;
	extern double MDN[ASTEP];
	
    /* phase: 1, historical simulation */
	grid->phase = 1; 
	
	(echar->soil).rl = (echar->soil).rl0;
	(echar->soil).rh = (echar->soil).rh0;
    
    /* all deforestation and regrowth: 2016/12/28  */
    if(EX_DEFOREST == 4 && (mass->soil).v_type == 1){
        
        mass_luc = 0.95* ((mass->c3).fol + (mass->c3).stm + (mass->c3).rot);
        
        (mass->soil).ltr += mass_luc;
        
        (mass->c3).fol *= 0.05;
        (mass->c3).stm *= 0.05;
        (mass->c3).rot *= 0.05;
    }
	
	/* LOOP to dynamic stage *******************************************************/
	for(g=0; g<PD_HIST; g++){
		/* AD1901 - 2002 / 2008 / 2009 */
        /* ISI-MIP1: 1950-2099 */
        /* GEOMIP: 1901-2005 */
        /* ISI-MIP2 (hist): 1901-2010 */
        /* ISI-MIP2b (1.5/2.0): 1661-2299 (2099) */
        /* NMIP: 1901-2012 => 1861–2015 */
        /* ISIMIP3a: 1850-2016 (FSY_HIST = 1850) */
        /* ISIMIP3b: 1601-2100 (FSY_HIST = 1601) */
        /* GCP2019: 1800-2021 (FSY_HIST = 1800) */

		/* simulation year ********************/
        /* updated: 2016/10/20 */
        grid->simy = FSY_HIST + g;
 
        /* climate change ********************/
		grid->climy = grid->simy;
        
        /* NMIP: 2015/11/19  **********/
        /* updated: 2016/10/20 */
        grid->niny = grid->simy;
        if(grid->niny < (FSY_HIST-1)){
            grid->niny = (FSY_HIST-1);
        }
        if(grid->niny > LSY_HIST){
            grid->niny = LSY_HIST;
        }
        
		/* CO2 year ********************/
		grid->co2y = grid->simy;
		/* sensitivity analysis: no CO2 rise */
		if(CC_CD == 2){
			grid->co2y = BGY_CO2Y;
			/* BGY_CO2Y = 1901 (usual setting) */
		}else if(CC_CD == 5){
            /* grid->co2y = 2000; */
            grid->co2y = 2006;
        }
        /* ISI-MIP no-co2-change run: 2012/07/31  */
        if((SCENARIO_ID==2005 ||SCENARIO_ID==2006 ||SCENARIO_ID==2007 ||SCENARIO_ID==2008 ||
            SCENARIO_ID==2015 ||SCENARIO_ID==2016 ||SCENARIO_ID==2017 ||SCENARIO_ID==2018 ||
            SCENARIO_ID==2025 ||SCENARIO_ID==2026 ||SCENARIO_ID==2027 ||SCENARIO_ID==2028 ||
            SCENARIO_ID==2035 ||SCENARIO_ID==2036 ||SCENARIO_ID==2037 ||SCENARIO_ID==2038 ||
            SCENARIO_ID==2045 ||SCENARIO_ID==2046 ||SCENARIO_ID==2047 ||SCENARIO_ID==2048) && grid->simy>=2000){
            /* fixed to AD2000 level */
            grid->co2y = 2000;
        }
        /* IIa: 2017/05/18  *****/
        if(ISIMIP2_FIXCD == 1 || (
            SCENARIO_ID==5011 || SCENARIO_ID==5021 || SCENARIO_ID==5031 || SCENARIO_ID==5041 ||
            SCENARIO_ID==5062 || SCENARIO_ID==5072 || SCENARIO_ID==5082 || SCENARIO_ID==5092 ||
            SCENARIO_ID==5064 || SCENARIO_ID==5074 || SCENARIO_ID==5084 || SCENARIO_ID==5094)){
            /* fix CO2 after 2006 */
            if(grid->simy >= 2006){
                grid->co2y = 2005;
            }else{
                ;
            }
        }
        if(GCP_FIXCD == 1){
            /* GCP-CH4: fix CO2 after 2006 */
            if(grid->simy >= 2007){
                grid->co2y = 2006;
            }else{
                ;
            }
        }
        
        if(ISIMIP_RUN == 6 &&
            (SCENARIO_ID==5126 || SCENARIO_ID==5130 || SCENARIO_ID==5134 || SCENARIO_ID==5135 || SCENARIO_ID==5136 ||
            SCENARIO_ID==5146 || SCENARIO_ID==5150 || SCENARIO_ID==5154 || SCENARIO_ID==5155 || SCENARIO_ID==5156 ||
            SCENARIO_ID==5166 || SCENARIO_ID==5170 || SCENARIO_ID==5174 || SCENARIO_ID==5175 || SCENARIO_ID==5176 ||
            SCENARIO_ID==5186 || SCENARIO_ID==5190 || SCENARIO_ID==5194 || SCENARIO_ID==5195 || SCENARIO_ID==5196 ||
            SCENARIO_ID==5206 || SCENARIO_ID==5210 || SCENARIO_ID==5214 || SCENARIO_ID==5215 || SCENARIO_ID==5216)
            ){
                if(grid->simy >= 2015){
                    /* fix CO2 after 2015 */
                    grid->co2y = 2015;
                }else{
                    ;
                }
        }
        if(grid->co2y < FDY_AGHG){
            grid->co2y = FDY_AGHG;
        }

        /* land-use year *****************/
        grid->lucy = grid->simy;
		
        /* NMIP **************/
        /* updated 2016/10/20  */
        if(NMIP_RUN == 1 || NMIP_RUN == 2 || NMIP_RUN == 3|| NMIP_RUN == 4){
            ;
        }else if(NMIP_RUN == 5){
            /* grid->climy = FSY_HIST;
            grid->co2y = FSY_HIST; */
            grid->niny = FSY_HIST;
            grid->lucy = FSY_HIST;
        }else if(NMIP_RUN == 6){
            /* grid->climy = FSY_HIST; */
            grid->co2y = FSY_HIST;
            grid->niny = FSY_HIST;
            grid->lucy = FSY_HIST;
        }else if(NMIP_RUN == 7){
            /* grid->climy = FSY_HIST; */ /* 1901 */
            grid->co2y = FSY_HIST;
            grid->niny = FSY_HIST;
            grid->lucy = FSY_HIST;
        }else if(NMIP_RUN == 20){
            grid->climy = 1901 + g%20; /* */ /* 1850 */
            grid->co2y = FSY_HIST;
            grid->niny = FSY_HIST;
            grid->lucy = FSY_HIST;
        }else if(NMIP_RUN == 21 || NMIP_RUN == 22 || NMIP_RUN == 23 || NMIP_RUN == 24 || NMIP_RUN == 25){
            if(g < 50){
                grid->climy = 1901 + g%20; /* */ /* 1850 */
            }
        }else if(NMIP_RUN == 26){
            if(g < 50){
                grid->climy = 1901 + g%20; /* */ /* 1850 */
            }
            grid->lucy = FSY_HIST;
        }else if(NMIP_RUN == 27){
            if(g < 50){
                grid->climy = 1901 + g%20; /* */ /* 1850 */
            }
            grid->co2y = FSY_HIST;
        }else if(NMIP_RUN == 28){
            grid->climy = 1901 + g%20; /* */ /* 1850 */
        }else if(NMIP_RUN == 29){
            grid->climy = 1901 + g%20; /* */ /* 1850 */
            grid->co2y = FSY_HIST;
            grid->lucy = FSY_HIST;
        }else if(NMIP_RUN == 30){
            if(g < 50){
                grid->climy = 1901 + g%20; /* */ /* 1850 */
            }
            grid->co2y = FSY_HIST;
            grid->niny = FSY_HIST;
            grid->lucy = FSY_HIST;
        }
        
        /* sensitivity run: 2017/06/24  */
        if(EX_NFERT == 101 && (EX_NFERT_SA == 4 || EX_NFERT_SA == 6)){
            if(grid->lucy >= 1950){
                grid->lucy = 1950;
            }
        }
        
        if(NMIP_RUN == 7){
            grid->climy = 1901; /* 1901 */
        }
        
        /* forced afforestation: 2020/11/09  */
        if((EX_FORCED_AFFOREST_2 >=1 && EX_FORCED_AFFOREST_2 <=12) && grid->lucy == EX_FORCED_AFFOREST_2_YR){
            
            (mass->c3).fol = (mass->c3).fol_p = (mass->c3).stm = (mass->c3).rot = INT_C;
            for(f=0;f<ASTEP;f++){
                (mass->c3).mfol[f] = (mass->c3).mstm[f] = (mass->c3).mrot[f] = INT_C;
            }

            (mass->c4).fol = (mass->c4).fol_p = (mass->c4).stm= (mass->c4).rot= INT_C;
            for(f=0;f<ASTEP;f++){
                (mass->c4).mfol[f] = (mass->c4).mstm[f] = (mass->c4).mrot[f] = INT_C;   
            }
            
            switch(EX_FORCED_AFFOREST_2){
                case 1: grid->veg_olson = 1; break;
                case 2: grid->veg_olson = 2; break;
                case 3: grid->veg_olson = 3; break;
                case 4: grid->veg_olson = 4; break;
                case 5: grid->veg_olson = 5; break;
                case 6: grid->veg_olson = 6; break;
                case 7: grid->veg_olson = 7; break;
                case 8: grid->veg_olson = 8; break;
                case 9: grid->veg_olson = 9; break;
                case 10: grid->veg_olson = 10; break;
                case 11: grid->veg_olson = 11; break;
                case 12: grid->veg_olson = 12; break;
                default: printf("WRONG AFFORESTATION\n"); break;
            }

            /* C3 */
            parameterC3(grid, &(echar->c3));
            /* C4 */
            parameterC4(grid, &(echar->c4));
        }
        
        /* experiments for trait modification: 2020/12/07 */
        if(EX_MOD_TRAIT_1 == 1){
            if(grid->simy >= 2020 && grid->simy <= 2029){
                (echar->c3).sla *= 1.01;
                (echar->c4).sla *= 1.01;
            }
        }
        if(EX_MOD_TRAIT_2 == 1){
            if(grid->simy >= 2020 && grid->simy <= 2029){
                (echar->soil).rl *= 0.99;
                (echar->soil).rh *= 1.01;
            }
        }

        /* ************/
        /* added : 2018/10/26 */
        if(EXTRA_CO2_FIX == 1){
            grid->co2y = 1901; /* */ /* 1901 */
        }
        if(EXTRA_CO2_FIX == 2){
            grid->lucy = 1950; /* */ /* 1950 */
        }
        if(EXTRA_CLIM_FIX == 1){
            grid->climy = 1901; /* 1901 */
        }
        
        if(EXTRA_LU_FIX == 1){
            grid->lucy = 1901; /* */ /* 1901 */
        }

        /* climate year modifications ************/
        if(grid->simy < BGY_CLIM){
            grid->climy = BGY_CLIM + g%20;
        }else if(grid->simy > (BGY_CLIM + DL_HCLIM - 1)){
            if(NCEP_RUN == 1 && grid->simy <(FDY_NCEP+DL_NCEP) ){
                ;
            }else{
                grid->climy = (BGY_CLIM + DL_HCLIM - 1);
            }
        }
        
        /* for considering leap years: 2014/09/29 */
        if(grid->simy%4 == 0){
            MDN[1] = 29.0;
            YDN = 366.0;
        }else{
            MDN[1] = 28.0;
            YDN = 365.0;
        }

 		if(grid->flag_histdata == 1){
            /* set climate variables  */
			set_hist_clim(grid);
		}
		
        /* **************************************************************/

		/* land-use change *************/
		f_cult_luc(grid);
		
		/*  Fertilizer input, historical change: 2010/05/11 */
		if(grid->rank_nat == 1){
			/* developing countries */
			f_fert = 2.0217112 / (1.0 + exp(0.049849599 * (2000.6575 - (double)grid->niny)))+0.0014929171;
		}else if(grid->rank_nat == 2){
			/* developed countries */
			f_fert = 0.92939393 / (1.0 + exp(0.044112692 * (2000.0097 - (double)grid->niny)))+0.53533202;
		}else{
            f_fert = 1.0;
        }
        
        if((echar->soil).v_type == 2){
            /* NMIP input: 2015/11/19 */
            if(NMIP_RUN >= 1 || EX_NFERT >= 1 || ISIMIP_RUN == 4 || ISIMIP_RUN == 5
                    || ISIMIP_RUN == 6 || EX_NFERT == 102){
                n_fertilizer_in(grid, loct);
                f_fert = 1.0; /* driven by data */
            }
        }else{
            loct->n_frtlz_in = 0.0;
            loct->n_manure_in = 0.0;
        }
		
		/* seasonal (monthly) loop ************************************************/
		for(f=0;f<ASTEP;f++){
			grid->m = f;
			
			/* initialize GHG fluxes ************/
			ghg_flux_zero(f, flux);
			
			/* CO2 condition */
			f_co2_trend(grid);
			co2_in_canopy(grid, loct, mass, flux);
						
			/* environmental condition *******************/
			f_dyn_loct(grid, loct, mass, echar);
            
            /* printf("%5.1lf ", grid->tmp_sfc[f]); ******/
			
			/* vegetation processes *********************/
			f_biome_processes(grid, loct, echar, mass, flux);

			/* VOC emission *****************/
			f_voc_emit_guenther97(grid, loct, echar, mass, flux);
			/* Plant CH4 emission *****************/
			f_ch4_emit_veg(grid, loct, echar, mass, flux);
   
            /* termite CH4 emission *****************/
            f_ch4_emit_termite(grid, loct, echar, mass, flux);

			f_plant_stand_budget(grid, loct, mass, flux);
			
			if(EX_ACCLM == 3){
				(flux->plant).lL[f] = flux->lL0[f];
			}

			(flux->soil).lL[f] = (flux->plant).lL[f];			
			(flux->soil).d13c_lL[f] = (flux->plant).d13c_lL[f];
			(flux->soil).d14c_lL[f] = (flux->plant).d14c_lL[f];
			
			/* soil processes *****************/
			soil_processes(grid, loct, &(echar->soil), &(mass->soil), &(flux->soil));
			flux->sr[f] = loct->c3ptn[f]*((flux->c3).rrm[f] + (flux->c3).rrg[f]) +
							loct->c4ptn[f]*((flux->c4).rrm[f] + (flux->c4).rrg[f]) + 
							(flux->soil).hr[f];
			
			if(NECB_DOC == 1 ){
				(mass->soil).msl -= (flux->soil).doc_boyer[f]/1000000.0;
				if((mass->soil).msl < 0.0){
					(mass->soil).msl = 0.0;
				}
			}

			/* fertilizaer input for croplands: revised  (2009/06/04) */
			/* NH4:NO3 ratio is based on inventories */
            /* this routine may not be activated when using REPLACE_OLSON_CROP option */
			if((echar->soil).v_type == 1){
			   if(grid->veg_olson == 29 || grid->veg_olson == 30 ||
                                    grid->veg_olson == 31 || grid->veg_olson == 32){
				   (flux->soil).n_fertin[f] = loct->n_frtlz_in * 1000.0 * f_fert;
                    if(NMIP_RUN >=20 && NMIP_RUN <=30){
                        (mass->soil).n_no3 += loct->n_frtlz_in_noy * 1000.0 * f_fert;
                        (mass->soil).n_nh4 += loct->n_frtlz_in_nh4 * 1000.0 * f_fert;
                    }else{
                        (mass->soil).n_no3 += loct->n_frtlz_in * 0.2 * 1000.0 * f_fert;
                        (mass->soil).n_nh4 += loct->n_frtlz_in * 0.8 * 1000.0 * f_fert;
                    }
                   (flux->soil).n_manurein[grid->m] = 0.0;
				}else{
				   (flux->soil).n_fertin[grid->m] = 0.0;
				   (flux->soil).n_manurein[grid->m] = 0.0;
				}
			}
			
			/*
			if((echar->soil).v_type == 1){
				(mass->soil).n_no3 += grid->f_crop_con * loct->n_frtlz_in * 0.2 * 1000.0;
				(mass->soil).n_nh4 += grid->f_crop_con * loct->n_frtlz_in * 0.8 * 1000.0;
			} */
            
            if(EX_NFERT == 101){
                if((echar->soil).v_type == 2){
                    /* 101: Nishina ESSD data: 2017/02/13  */
                    
                    if(grid->f_crop_con > 0.0){
                        icrop = 1.0 / grid->f_crop_con;
                    }else{
                        icrop = 0.0;
                    }
                    
                    /* sensitivity run: 2017/06/24  */
                    y_nin = grid->niny;
                    
                    if(EX_NFERT_SA == 1 || EX_NFERT_SA == 5 || EX_NFERT_SA == 6){
                        (flux->soil).n_fertin[grid->m] = icrop * (grid->nin_no3[0][grid->m]
                                                        + grid->nin_nh4[0][grid->m]) * 1000.0;
                        (mass->soil).n_no3 += icrop * grid->nin_no3[0][grid->m] * 1000.0;
                        (mass->soil).n_nh4 += icrop * grid->nin_nh4[0][grid->m] * 1000.0;
                    }else{
                        if(y_nin < 1960){
                            (flux->soil).n_fertin[grid->m] = icrop * (grid->nin_no3[0][grid->m]
                                                        + grid->nin_nh4[0][grid->m]) * 1000.0;
                            (mass->soil).n_no3 += icrop * grid->nin_no3[0][grid->m] * 1000.0;
                            (mass->soil).n_nh4 += icrop * grid->nin_nh4[0][grid->m] * 1000.0;
                        }else if(y_nin >= 1960 && y_nin <= 2009){
                            (flux->soil).n_fertin[grid->m] = icrop * (grid->nin_no3[y_nin - 1960][grid->m]
                                                        + grid->nin_nh4[y_nin - 1960][grid->m]) * 1000.0;
                            (mass->soil).n_no3 += icrop * grid->nin_no3[y_nin - 1960][grid->m] * 1000.0;
                            (mass->soil).n_nh4 += icrop * grid->nin_nh4[y_nin - 1960][grid->m] * 1000.0;
                        }else if(y_nin > 2009){
                            (flux->soil).n_fertin[grid->m] = icrop * (grid->nin_no3[49][grid->m]
                                                        + grid->nin_nh4[49][grid->m]) * 1000.0;
                            (mass->soil).n_no3 += icrop * grid->nin_no3[49][grid->m] * 1000.0;
                            (mass->soil).n_nh4 += icrop * grid->nin_nh4[49][grid->m] * 1000.0;
                        }
                    }
                    
                    base_nin = (flux->soil).n_fertin[grid->m];

                    /* sensitivity run: 2017/06/24  */
                    if(EX_NFERT_SA == 2 || EX_NFERT_SA == 5 || EX_NFERT_SA == 6){
                        /* fixed at 1960 */
                        base_nin = icrop * (grid->nin_no3[0][grid->m] + grid->nin_nh4[0][grid->m]) * 1000.0;
                    }else if(EX_NFERT_SA == 1){
                        /* change only for mature */
                        if(y_nin < 1960){
                            base_nin = icrop * (grid->nin_no3[0][grid->m]
                                            + grid->nin_nh4[0][grid->m]) * 1000.0;
                        }else if(y_nin >= 1960 && y_nin <= 2009){
                            base_nin = icrop * (grid->nin_no3[y_nin - 1960][grid->m]
                                            + grid->nin_nh4[y_nin - 1960][grid->m]) * 1000.0;
                        }else if(y_nin > 2009){
                            base_nin = icrop * (grid->nin_no3[49][grid->m]
                                            + grid->nin_nh4[49][grid->m]) * 1000.0;
                        }
                    }
                    
                    if(grid->nfert_potter > 0.0 || grid->nmanure_potter > 0.0){
                    
                        if(grid->nfert_potter > 0.0){
                            aa = grid->nmanure_potter / grid->nfert_potter;
                            bb = aa * (base_nin/1000.0);
                            
                            if(bb > 20.0*(grid->nmanure_potter * MDN[grid->m] / YDN) ){
                                bb = 20.0*(grid->nmanure_potter * MDN[grid->m] / YDN);
                            }
                        }else{
                            bb = icrop * grid->nmanure_potter  * MDN[grid->m] / YDN;
                        }
                        
                        if(bb < 0.0){
                            bb = 0.0;
                        }
                        
                        /* No manure: 2017/07/10  */
                        if(EX_NFERT_SA == 7){
                            bb = 0.0;
                        }
                        
                        loct->n_manure_in = bb * 1000.0;
                        (flux->soil).n_manurein[grid->m] = loct->n_manure_in;
                        (mass->soil).n_lttr += loct->n_manure_in;
                    }else{
                        loct->n_manure_in = (flux->soil).n_manurein[grid->m] = 0.0;
                    }
                }
            }else{
                if((echar->soil).v_type == 2){
                    (flux->soil).n_fertin[grid->m] = loct->n_frtlz_in * 1000.0 * f_fert;
                     if(NMIP_RUN >= 20 && NMIP_RUN <= 30){
                        (mass->soil).n_no3 += loct->n_frtlz_in_noy * 1000.0 * f_fert;
                        (mass->soil).n_nh4 += loct->n_frtlz_in_nh4 * 1000.0 * f_fert;
                    }else{
                        (mass->soil).n_no3 += loct->n_frtlz_in * 0.2 * 1000.0 * f_fert;
                        (mass->soil).n_nh4 += loct->n_frtlz_in * 0.8 * 1000.0 * f_fert;
                    }

                    /* 2016/10/20  */
                    (flux->soil).n_manurein[grid->m] = loct->n_manure_in * 1000.0 * f_fert;
                    (mass->soil).n_lttr += loct->n_manure_in * 1000.0 * f_fert;
                }
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
			if(CH4_WH == 1 && grid->f_wetland>0.0 && loct->v_type == 1){
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
			if(CH4_WH == 1 && grid->f_paddy>0.0 && loct->v_type == 2){
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
                /* revised (after comments by E.Kato): 2013/10/02  */
                (mass->soil).msl += 12.0/16.0 * (grid->f_upland * (flux->soil).ch4oxy_curry[f] * 0.00001
					- grid->f_paddy * ((flux->soil).ch4_paddy_wh_plant[f] + (flux->soil).ch4_paddy_wh_ebull[f] +
                            (flux->soil).ch4_paddy_wh_diff[f] + (flux->soil).ch4_paddy_wh_release[f]) * 0.00001
                    - grid->f_wetland * ((flux->soil).ch4_wetland_wh_plant[f] + (flux->soil).ch4_wetland_wh_ebull[f] +
                            (flux->soil).ch4_wetland_wh_diff[f] + (flux->soil).ch4_wetland_wh_release[f]) * 0.00001);
				
				if((mass->soil).msl < INT_C){
					(mass->soil).msl = INT_C;
				}
			}
			
			/* N2O emission ***************************************/
			f_n2o_emit_ngas(grid, loct, mass, flux);
			f_n2o_emit_casa(grid, loct, mass, flux);
			
			/* ecosystem mass balance *****************************/	
			/* net ecosystem production */
			flux->nep[f] = (flux->plant).npp[f] - (flux->soil).hr[f];
			flux->er[f] = (flux->plant).ar[f] + (flux->soil).hr[f];
			/* total ecosystem carbon storage */
			mass->total[f] = (mass->c3).plant[f]*loct->c3ptn[f]
                    + (mass->c4).plant[f]*loct->c4ptn[f] + (mass->soil).soil[f];
			/** net carbon balance taking crop harvest into account **/
			flux->ncb[f] = flux->nep[f] - (flux->plant).net_crop[f];
			
			/* carbon isotope */
			f_cisotope_efflux(grid, loct, mass, flux);
			
			/* nitrogen budget */
			n_budget(grid, loct, mass, flux);
			
			/* average LAI: 2009/05/06  */
			if(grid->niny >= 1990 && grid->niny <= 1999){
				(mass->c3).lai0[f] += (mass->c3).lai[f]/10.0;
				(mass->c4).lai0[f] += (mass->c4).lai[f]/10.0;
				(mass->plant).lai0[f] += (mass->plant).lai[f]/10.0;
			}
			
			/* statistics *******************************************/
			/* if(g>=90 && g<=99){ */
            if(grid->simy >= 2000 && grid->simy <= 2009){
				/* mean seasonal change *******/
				m_ch4ox1[f] += (flux->soil).ch4oxy_ridg[f] * grid->area *10000.0/1000.0 / 10.0;
				m_ch4ox2[f] += (flux->soil).ch4oxy_casa[f] * grid->area *10000.0/1000.0 / 10.0;
				m_ch4ox3[f] += (flux->soil).ch4oxy_delgrosso[f] * grid->area *10000.0/1000.0 / 10.0;
				
				m_gpp[f] += (flux->plant).gpp[f]/10.0 * grid->area;
				m_npp[f] += (flux->plant).npp[f]/10.0 * grid->area;
				m_nep[f] += flux->nep[f]/10.0 * grid->area;
				
				if(DF97 == 1){
					vo_gpp[grid->veg_olson] += (flux->plant).gpp_df97[f]/10.0 * grid->area;
				}else{
					vo_gpp[grid->veg_olson] += (flux->plant).gpp[f]/10.0 * grid->area;
				}
				
				/* mean biome budget *******/
				vo_npp[grid->veg_olson] += (flux->plant).npp[f]/10.0 * grid->area;
				vo_nep[grid->veg_olson] += flux->nep[f]/10.0 * grid->area;
				vo_lai[grid->veg_olson] += (mass->plant).lai[f] * MDN[f]/YDN/10.0 * grid->area;
				vo_fol[grid->veg_olson] += (mass->plant).mfol[f] * MDN[f]/YDN/10.0 * grid->area;
				vo_stm[grid->veg_olson] += (mass->plant).mstm[f] * MDN[f]/YDN/10.0 * grid->area;
				vo_rot[grid->veg_olson] += (mass->plant).mrot[f] * MDN[f]/YDN/10.0 * grid->area;
				vo_ltr[grid->veg_olson] += (mass->soil).ltr_m[f] * MDN[f]/YDN/10.0 * grid->area;
				vo_msl[grid->veg_olson] += (mass->soil).msl_m[f] * MDN[f]/YDN/10.0 * grid->area;
                
                vo_n_cnpy[grid->veg_olson] += (mass->plant).n_cnpy_m[f] * MDN[f]/YDN/10.0 * grid->area;
                vo_n_strg[grid->veg_olson] += (mass->plant).n_strg_m[f] * MDN[f]/YDN/10.0 * grid->area;
                vo_n_mcrb[grid->veg_olson] += (mass->soil).n_mcrb_m[f] * MDN[f]/YDN/10.0 * grid->area;
                vo_n_ltr[grid->veg_olson] += (mass->soil).n_lttr_m[f] * MDN[f]/YDN/10.0 * grid->area;
                vo_n_hms[grid->veg_olson] += (mass->soil).n_hums_m[f] * MDN[f]/YDN/10.0 * grid->area;
				
				if(DF97 == 1){
					vs_gpp[grid->veg_sage] += (flux->plant).gpp_df97[f]/10.0 * grid->area;
				}else{
					vs_gpp[grid->veg_sage] += (flux->plant).gpp[f]/10.0 * grid->area;
				}
				
				vs_npp[grid->veg_sage] += (flux->plant).npp[f]/10.0 * grid->area;
				vs_nep[grid->veg_sage] += flux->nep[f]/10.0 * grid->area;
				vs_lai[grid->veg_sage] += (mass->plant).lai[f]*MDN[f]/YDN/10.0 * grid->area;
				vs_fol[grid->veg_sage] += (mass->plant).mfol[f]*MDN[f]/YDN/10.0 * grid->area;
				vs_stm[grid->veg_sage] += (mass->plant).mstm[f]*MDN[f]/YDN/10.0 * grid->area;
				vs_rot[grid->veg_sage] += (mass->plant).mrot[f]*MDN[f]/YDN/10.0 * grid->area;
				vs_ltr[grid->veg_sage] += (mass->soil).ltr_m[f]*MDN[f]/YDN/10.0 * grid->area;
				vs_msl[grid->veg_sage] += (mass->soil).msl_m[f]*MDN[f]/YDN/10.0 * grid->area;
                
                if(loct->v_type == 1){
                    fweight = 1.0 - grid->f_crop_con;
                }else{
                    fweight = grid->f_crop_con;
                }
                glat_gpp[f][grid->row] += fweight * grid->area * (flux->plant).gpp[f]/10.0;
                glat_npp[f][grid->row] += fweight * grid->area * (flux->plant).npp[f]/10.0;
                glat_nep[f][grid->row] += fweight * grid->area * flux->nep[f]/10.0;

                glat_agb[grid->row] += fweight * grid->area * ((mass->plant).mfol[f] + (mass->plant).mstm[f])*MDN[f]/YDN/10.0;
                glat_soc[grid->row] += fweight * grid->area * ((mass->soil).ltr_m[f] + (mass->soil).msl_m[f])*MDN[f]/YDN/10.0;

                if(loct->v_type == 1){
                    glat_ch4_cao[f][grid->row] += grid->area * (flux->soil).ch4flux_wetland_cao[f] / 10.0;
                    glat_ch4_wh[f][grid->row] += grid->area * ((flux->soil).ch4_wetland_wh_plant[f]
                                + (flux->soil).ch4_wetland_wh_ebull[f] + (flux->soil).ch4_wetland_wh_diff[f]
								+ (flux->soil).ch4_wetland_wh_release[f]) / 10.0;
                }else{
                    glat_ch4_cao[f][grid->row] += grid->area * (flux->soil).ch4flux_paddy_cao[f] / 10.0;
                    glat_ch4_wh[f][grid->row] += grid->area * ((flux->soil).ch4_paddy_wh_plant[f]
                                + (flux->soil).ch4_paddy_wh_ebull[f] + (flux->soil).ch4_paddy_wh_diff[f]
								+ (flux->soil).ch4_paddy_wh_release[f]) / 10.0;
                }
			}
			
			/******************/
			f_grid_av(grid, loct, echar, mass, flux);
		}
        /* printf("\n"); */
		
		/* empirical NPP models */
		npp_empirical(grid, loct, flux);
		
		/* biomass burning **************************/
		f_biomassburning(grid, loct, mass, flux);
		/* corrected: A.Ito and E.Kato (2009/08/16) */
		if(g>=90 && g<=99){
			for(f=0;f<ASTEP;f++){
				m_bioburn_co2[f] += (flux->bb_co2_litter[f] + flux->bb_co2_leaf[f]
									 +flux->bb_co2_wood[f] + flux->bb_co2_root[f]) * grid->area / 10.0;
				m_bioburn_co[f] += (flux->bb_co_litter[f] + flux->bb_co_leaf[f]
									+flux->bb_co_wood[f] + flux->bb_co_root[f]) * grid->area / 10.0;
				m_bioburn_ch4[f] += (flux->bb_ch4_litter[f] + flux->bb_ch4_leaf[f]
									 +flux->bb_ch4_wood[f] + flux->bb_ch4_root[f]) * grid->area / 10.0;
				m_bioburn_nmhc[f] += (flux->bb_nmhc_litter[f] + flux->bb_nmhc_leaf[f]
									  +flux->bb_nmhc_wood[f] + flux->bb_nmhc_root[f]) * grid->area / 10.0;
				m_bioburn_oc[f] += (flux->bb_oc_litter[f] + flux->bb_oc_leaf[f]
									+flux->bb_oc_wood[f] + flux->bb_oc_root[f]) * grid->area / 10.0;
				m_bioburn_bc[f] += (flux->bb_bc_litter[f] + flux->bb_bc_leaf[f]
									+flux->bb_bc_wood[f] + flux->bb_bc_root[f]) * grid->area / 10.0;
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
		
        /* C-budget parameter ensemble: 2018/06/05  */
        if(PARAM_PTB == 20){
            prm_ensen = 1.0 + 0.3 * f_pert[7];
        }else{
            prm_ensen = 1.0;
        }
    
		if(NECB_ERSN == 1){
			(mass->soil).ltr -= flux->erod_carbon * (prm_ensen * 0.20);
			if((mass->soil).ltr < INT_C){
				(mass->soil).ltr = INT_C;
			}
		}

		/* land use change *********************/
		if(loct->v_type == 1){
            /* natural vegetation */
			f_luc_emit(grid, mass, flux);
		}else{
			flux->lu_detr = 0.0;
			flux->lu_conv = 0.0;
			flux->lu_ten = 0.0;
			flux->lu_hund = 0.0;
   
            flux->lu_fol = 0.0;
            flux->lu_stm = 0.0;
            flux->lu_rot = 0.0;
            flux->lu_ltr = 0.0;
            flux->lu_msl = 0.0;
		}
		
        /* fraction of natural vegetation area */
        f_nat = 1.0 - grid->f_crop_con;
        
        if(f_nat > 0.0){
            iweight = 1.0 / f_nat;
            /* factor from natural to whole grid */
        }else{
            iweight = 1.0;
        }
        avc3 = 0.0;
        for(f=0;f<ASTEP;f++){
            avc3 += (loct->c3ptn[f] * MDN[f]/YDN);
        }
        if(avc3 > 0.0){
            iweight3 = iweight * (1.0 / avc3);
            /* factor from C3-dominated natural to whole grid */
        }else{
            iweight3 = iweight;
        }
        
        if(NECB_LUC == 1 && EX_BECCS == 0){
            /* 2018/10/24  */
            (mass->c3).fol -= flux->lu_fol * iweight3;
            (mass->c3).stm -= flux->lu_stm * iweight3;
            (mass->c3).rot -= flux->lu_rot * iweight3 + flux->lu_ltr * iweight3;
            
            if((mass->c3).fol < INT_C){
                (mass->c3).fol = INT_C;
            }
            if((mass->c3).stm < INT_C){
                (mass->c3).stm = INT_C;
            }
            if((mass->c3).rot < INT_C){
                (mass->c3).rot = INT_C;
            }

            (mass->soil).ltr += flux->lu_ltr + flux->detr_ten[0] + flux->detr_hund[0];
            (mass->soil).ltr -= flux->lu_ten + flux->lu_hund;

            if((mass->soil).ltr < INT_C){
                (mass->soil).ltr = INT_C;
            }

            (mass->c3).mfol[ASTEP-1] = (mass->c3).fol;
            (mass->c3).mstm[ASTEP-1] = (mass->c3).stm;
            (mass->c3).mrot[ASTEP-1] = (mass->c3).rot;

            (mass->plant).fol = (mass->c3).fol*loct->c3ptn[ASTEP-1] + (mass->c4).fol*loct->c4ptn[ASTEP-1];
            (mass->plant).stm = (mass->c3).stm*loct->c3ptn[ASTEP-1] + (mass->c4).stm*loct->c4ptn[ASTEP-1];
            (mass->plant).rot = (mass->c3).rot*loct->c3ptn[ASTEP-1] + (mass->c4).rot*loct->c4ptn[ASTEP-1];
            (mass->plant).mfol[ASTEP-1] = (mass->plant).fol;
            (mass->plant).mstm[ASTEP-1] = (mass->plant).stm;
            (mass->plant).mrot[ASTEP-1] = (mass->plant).rot;

            (mass->c3).plant[ASTEP-1] = (mass->c3).fol + (mass->c3).stm + (mass->c3).rot;
            (mass->c4).plant[ASTEP-1] = (mass->c4).fol + (mass->c4).stm + (mass->c4).rot;
            (mass->plant).plant[ASTEP-1] = (mass->c3).plant[ASTEP-1]*loct->c3ptn[ASTEP-1]
                                        + (mass->c4).plant[ASTEP-1]*loct->c4ptn[ASTEP-1];

            (mass->c3).lai[ASTEP-1] = lai_mass(grid, &(mass->c3), &(echar->c3));
            (mass->c4).lai[ASTEP-1] = lai_mass(grid, &(mass->c4), &(echar->c4));
            (mass->plant).lai[ASTEP-1] = (mass->c3).lai[ASTEP-1] * loct->c3ptn[ASTEP-1]
                                        + (mass->c4).lai[ASTEP-1] * loct->c4ptn[ASTEP-1];
            loct->lai[ASTEP-1] = (mass->plant).lai[ASTEP-1];

        }else if(NECB_LUC == 1 && (EX_BECCS == 1 || EX_BECCS == 2 || EX_BECCS == 3)){
            /* used by BECCS experiment */
            if(grid->f_luc>0.0 && grid->f_luc<1.0){
                (mass->c3).fol *= (1.0 - grid->f_luc);
                (mass->c3).stm *= (1.0 - grid->f_luc);
                (mass->c3).rot *= (1.0 - grid->f_luc);
                (mass->c4).fol *= (1.0 - grid->f_luc);
                (mass->c4).stm *= (1.0 - grid->f_luc);
                (mass->c4).rot *= (1.0 - grid->f_luc);

                (mass->c3).mfol[ASTEP-1] = (mass->c3).fol;
                (mass->c3).mstm[ASTEP-1] = (mass->c3).stm;
                (mass->c3).mrot[ASTEP-1] = (mass->c3).rot;
                (mass->c4).mfol[ASTEP-1] = (mass->c4).fol;
                (mass->c4).mstm[ASTEP-1] = (mass->c4).stm;
                (mass->c4).mrot[ASTEP-1] = (mass->c4).rot;

                (mass->plant).fol = (mass->c3).fol*loct->c3ptn[ASTEP-1] + (mass->c4).fol*loct->c4ptn[ASTEP-1];
                (mass->plant).stm = (mass->c3).stm*loct->c3ptn[ASTEP-1] + (mass->c4).stm*loct->c4ptn[ASTEP-1];
                (mass->plant).rot = (mass->c3).rot*loct->c3ptn[ASTEP-1] + (mass->c4).rot*loct->c4ptn[ASTEP-1];
                (mass->plant).mfol[ASTEP-1] = (mass->plant).fol;
                (mass->plant).mstm[ASTEP-1] = (mass->plant).stm;
                (mass->plant).mrot[ASTEP-1] = (mass->plant).rot;

                (mass->c3).plant[ASTEP-1] = (mass->c3).fol + (mass->c3).stm + (mass->c3).rot;
                (mass->c4).plant[ASTEP-1] = (mass->c4).fol + (mass->c4).stm + (mass->c4).rot;
                (mass->plant).plant[ASTEP-1] = (mass->c3).plant[ASTEP-1]*loct->c3ptn[ASTEP-1]
                                            + (mass->c4).plant[ASTEP-1]*loct->c4ptn[ASTEP-1];

                (mass->c3).lai[ASTEP-1] = lai_mass(grid, &(mass->c3), &(echar->c3));
                (mass->c4).lai[ASTEP-1] = lai_mass(grid, &(mass->c4), &(echar->c4));
                (mass->plant).lai[ASTEP-1] = (mass->c3).lai[ASTEP-1] * loct->c3ptn[ASTEP-1]
                                            + (mass->c4).lai[ASTEP-1] * loct->c4ptn[ASTEP-1];
                loct->lai[ASTEP-1] = (mass->plant).lai[ASTEP-1];
            }else{
                ;
            }
        }
        
		/* wood harvest: 2010/10/15  ***************/
		total_hvst = 0.0;
		if((mass->c3).v_type == 1){
			dyr = grid->lucy - FDY_LUC;
			
			/* assumption for the period later than 2004: A.Ito (2010/11/11) */
			if( (LANDUSE != 10 && LANDUSE != 11 && LANDUSE != 12 && LANDUSE != 13) &&
                    grid->lucy > (FDY_LUC+DL_LUC-1)){
				dyr = (FDY_LUC + DL_LUC - 1);
			}

            /* assumption for the period later than 2016: A.Ito (2019/02/11) */
            if( (LANDUSE == 26 || LANDUSE == 27 || LANDUSE == 28) &&
                    grid->lucy > (FDY_LUC+DL_LUC-1)){
                dyr = (DL_LUC - 1);
            }

            /* from total grid */
			/* total_hvst = grid->hvst_p1[dyr] + grid->hvst_p2[dyr] + grid->hvst_s1[dyr]
						+ grid->hvst_s2[dyr] + grid->hvst_s3[dyr]; */

            total_hvst = grid->hvst_p1[dyr] + grid->hvst_s1[dyr];

            /* parameter ensemble: 2014/11/19  */
            prm_ensen = 1.0;
            if(PARAM_PTB == 9){
                if(PARAM_ENS==1){
                    prm_ensen *= 0.7;
                }
                if(PARAM_ENS==2){
                    prm_ensen *= 0.8;
                }
                if(PARAM_ENS==3){
                    prm_ensen *= 0.9;
                }
                if(PARAM_ENS==4){
                    prm_ensen *= 1.1;
                }
                if(PARAM_ENS==5){
                    prm_ensen *= 1.2;
                }
                if(PARAM_ENS==6){
                    prm_ensen *= 1.3;
                }
            }
			
            /* C-budget parameter ensemble: 2018/06/05  */
            if(PARAM_PTB == 20){
                prm_ensen = 1.0 + 0.3 * f_pert[1];
            }
            
            /* Forest management: 2019/10/17  */
            if(EX_FORMAN == 2){
                prm_ensen = 0.5;
            }
            if(EX_FORMAN == 3){
                prm_ensen = 2.0;
            }

            /* into MgC/ha */
			total_hvst *= 1.0/1000.0 * 1.0/grid->area * prm_ensen;
   
            /* Note: only C3 tree species */
            if(NECB_WHVST == 1){
                if((mass->c3).stm > (total_hvst + INT_C)){
            
                    if((mass->c3).stm > (total_hvst * iweight3 + INT_C)){
                        (mass->c3).stm -= total_hvst * iweight3;
                        flux->hvst_wood = total_hvst * iweight3;
                    }else{
                        (mass->c3).stm -= total_hvst;
                        flux->hvst_wood = total_hvst;
                    }
                }else{
                    flux->hvst_wood = total_hvst - INT_C;
                    if(flux->hvst_wood < 0.0){
                        flux->hvst_wood = 0.0;
                    }
                    (mass->c3).stm = INT_C;
                }
                
                if((mass->c3).stm < INT_C){
                    (mass->c3).stm = INT_C;
                }
            }
		}else{
			flux->hvst_wood = 0.0;
		}
        
		/* net biome production (added : 2010/01/20) *************************/
        flux->gpp_ann = 0.0;      
		for(f=0;f<ASTEP;f++){
            grid->m = f;
            
            /* base */
			flux->nbp[f] = flux->nep[f];
   
            flux->gpp_ann += (flux->plant).gpp[f];
            
            /* altered: 2018/10/16  */
			if((mass->c3).v_type == 1 && NECB_LUC == 1){
				flux->nbp[f] -= iweight * (flux->lu_conv/(double)ASTEP + flux->lu_ten/(double)ASTEP + flux->lu_hund/(double)ASTEP);
			}
			
            if(NECB_WHVST == 1){
                flux->nbp[f] -= flux->hvst_wood /(double)ASTEP;
            }

			if(NECB_BB == 1){
                /* revised (after comments by E.Kato): 2013/10/02  */
            
				flux->nbp[f] -= (flux->bb_co2_litter[f] + flux->bb_co2_leaf[f] 
								 + flux->bb_co2_wood[f] + flux->bb_co2_root[f])/1000.0 * 12.0/44.0;

                flux->nbp[f] -= (flux->bb_co_litter[f] + flux->bb_co_leaf[f] 
                                 + flux->bb_co_wood[f] + flux->bb_co_root[f])/1000.0 * 12.0/28.0;

                flux->nbp[f] -= (flux->bb_ch4_litter[f] + flux->bb_ch4_leaf[f]
                                 + flux->bb_ch4_wood[f] + flux->bb_ch4_root[f])/1000.0 * 12.0/16.0;

                flux->nbp[f] -= (flux->bb_bc_litter[f] + flux->bb_bc_leaf[f]
                                 + flux->bb_bc_wood[f] + flux->bb_bc_root[f])/1000.0;
			}
            
            if(NECB_DOC == 1){
                /* revised (after comments by E.Kato): 2013/10/02  */
                flux->nbp[f] -= (flux->soil).doc_boyer[f]/1000000.0;
            }
            
            if(NECB_CH4 == 1){
                /* revised (after comments by E.Kato): 2013/10/02  */
                flux->nbp[f] += 12.0/16.0 * (grid->f_upland * (flux->soil).ch4oxy_curry[f] * 0.00001
						- grid->f_paddy * ((flux->soil).ch4_paddy_wh_plant[f] + (flux->soil).ch4_paddy_wh_ebull[f] + 
									(flux->soil).ch4_paddy_wh_diff[f] + (flux->soil).ch4_paddy_wh_release[f]) * 0.00001
						- grid->f_wetland * ((flux->soil).ch4_wetland_wh_plant[f] + (flux->soil).ch4_wetland_wh_ebull[f] + 
									(flux->soil).ch4_wetland_wh_diff[f] + (flux->soil).ch4_wetland_wh_release[f]) * 0.00001);
            }
            
            if(NECB_ERSN == 1){
                /* C-budget parameter ensemble: 2018/06/14  */
                if(PARAM_PTB == 20){
                    prm_ensen = 1.0 + 0.3 * f_pert[7];
                }else{
                    prm_ensen = 1.0;
                }
    
                /* revised (after comments by E.Kato): 2013/10/02  */
                flux->nbp[f] -= flux->erod_carbon * (prm_ensen * 0.20) /(double)ASTEP;
            }
            
            if(NECB_BVOC == 1){
                /* revised (after comments by E.Kato): 2013/10/02  */
                flux->nbp[f] -= (flux->voc_isopr_g97[f] + flux->voc_monotrp_g97[f] + flux->voc_methanl_g97[f] +
                    flux->voc_acetone_g97[f] + flux->voc_actaldhd_g97[f] + flux->voc_frmardhd_g97[f] +
                    flux->voc_formacd_g97[f] + flux->voc_acetacd_g97[f] + flux->voc_co_g97[f] +
                    flux->voc_afarnesene[f] + flux->voc_bcaryophyllene[f] + flux->voc_othersesqui[f]) * 10000.0/1000000.0/1000000.0;
            }
            
            if(NECB_CROP == 1){
                /* revised (after comments by E.Kato): 2013/10/02  */
                flux->nbp[f] -= 1.0 * (flux->plant).net_crop[f]; /* ! hvst is positive */
            }
		}
        
		/* history data */
		f_set_history_data(grid->simy - FSY_HIST +1, grid, loct, mass, flux);
		
		/* output */
		f_output_result(grid->simy, grid, loct, echar, mass, flux, fp_o); /* */
		
		/* setting previous land-use */
		grid->f_crop_p = grid->f_crop_con;
		grid->f_pasture_p = grid->f_pasture_con;
		
	}
    
    for(f=0;f<OFILEN;f++){
        fprintf(fp_o[f],"\n");
    }
}
