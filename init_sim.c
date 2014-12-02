/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* Revised August 15, 2007 by A.Ito */
/* Revised August 19, 2007 by A.Ito */
/* Revised July 1, 2008 by A.Ito */
/* Separated February 17, 2009 by A.Ito */

/* initialize environmental characteristics of each grid */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* initial carbon pool size */
#define INT_C 0.01

/* initialize general parameters used in the simulation ******************/
void f_init_sim(
	struct Grid *grid
){
	long f, g, h, year;
	double data;
    float rdat[36*72],albvar;
	FILE *fpi;
		
	/**********************************************/
	set_rowcol_gcm();	/* -> vegetdeal.c */

	/* grid ID, sequential number */
	grid->n_olson = 0;
	grid->n_sage = 0;
	grid->n_crop = 0;
    
	/**********************************************/
    /* albedo perturbation */
    if(EX_ALBEDO >= 1){
        /* CMIP5: by A.Ito */
        fpi = fopen("./data/albedo_cmip_5deg_2.flt","rb");
        
        for(f=0;f<12;f++){
            fread(rdat,sizeof(float),36*72, fpi);
            for(g=0;g<36;g++){
                for(h=0;h<72;h++){
                    grid->albedo_av[f][g][h] = rdat[g*72+h];
                }
            }
            
            fread(rdat,sizeof(float),36*72, fpi);
            for(g=0;g<36;g++){
                for(h=0;h<72;h++){
                    if(rdat[g*72+h]>0.0 && rdat[g*72+h]<1.0){
                        albvar = sqrt(rdat[g*72+h]);
                    }else{
                        albvar = 0.0;
                    }
                    grid->albedo_sd[f][g][h] = albvar;
                }
            }
            
            fread(rdat,sizeof(float),36*72, fpi);
            for(g=0;g<36;g++){
                for(h=0;h<72;h++){
                    grid->albedo_max[f][g][h] = rdat[g*72+h];
                }
            }
            
            fread(rdat,sizeof(float),36*72, fpi);
            for(g=0;g<36;g++){
                for(h=0;h<72;h++){
                    grid->albedo_min[f][g][h] = rdat[g*72+h];
                }
            }
        }
        fclose(fpi);
        
    }else{
        grid->albedo_av[0][0][0] = 0.2;
        grid->albedo_sd[0][0][0] = 0.1;
    }
		
	/* atm. CO2 scenario ****************************/
	printf("reading CO2 data...");
    if(ISIMIP_RUN == 0){
        if(CO2S==1){
            if((fpi = fopen("./data/SRES_A1.dat","rt"))==NULL){
                printf("No SRES_A1.dat\n");
                exit(1);
            }
        }else if(CO2S==2){
            if((fpi = fopen("./data/SRES_A1FI.dat","rt"))==NULL){
                printf("No SRES_A1FI.dat\n");
                exit(1);
            }
        }else if(CO2S==3){
            if((fpi = fopen("./data/SRES_A1T.dat","rt"))==NULL){
                printf("No SRES_A1T.dat\n");
                exit(1);
            }
        }else if(CO2S==4){
            if((fpi = fopen("./data/SRES_A2.dat","rt"))==NULL){
                printf("No SRES_A2.dat\n");
                exit(1);
            }
        }else if(CO2S==5){
            if((fpi = fopen("./data/SRES_B1.dat","rt"))==NULL){
                printf("No SRES_B1.dat\n");
                exit(1);
            }
        }else if(CO2S==6){
            if((fpi = fopen("./data/SRES_B2.dat","rt"))==NULL){
                printf("No SRES_B2.dat\n");
                exit(1);
            }
        }else{
            if((fpi = fopen("./data/SRES_A1.dat","rt"))==NULL){
                printf("No SRES_A1.dat\n");
                exit(1);
            }
        }
        
        /* SRES Scenario CO2 */
        for(f=0;f<111;f++){
            fscanf(fpi,"%ld %lf", &year, &data);
            sres_co2[f] = data;
        }
        
        fclose(fpi);
        
        /* source: http://crga.atmos.uiuc.edu/research/post-sres.html
            M.E.Schlesinger and S.Malyshev			*/
        if((fpi = fopen("./data/AtmGHG_timeseries.dat","rt"))==NULL){
            printf("No AtmGHG_timeseries.dat\n");
            exit(1);
        }
        for(f=0;f<DL_AGHG;f++){
            fscanf(fpi,"%ld", &year);
            /* CO2, ppmv */
            fscanf(fpi,"%lf", &aco2_a1[f]);
            fscanf(fpi,"%lf", &aco2_b1[f]);
            fscanf(fpi,"%lf", &aco2_a2[f]);
            fscanf(fpi,"%lf", &aco2_b2[f]);
            /* CH4, ppbv */
            fscanf(fpi,"%lf", &ach4_a1[f]);
            fscanf(fpi,"%lf", &ach4_b1[f]);
            fscanf(fpi,"%lf", &ach4_a2[f]);
            fscanf(fpi,"%lf", &ach4_b2[f]);
            /* N2O, ppbv */
            fscanf(fpi,"%lf", &an2o_a1[f]);
            fscanf(fpi,"%lf", &an2o_b1[f]);
            fscanf(fpi,"%lf", &an2o_a2[f]);
            fscanf(fpi,"%lf", &an2o_b2[f]);
        }
        fclose(fpi);
    }else if(ISIMIP_RUN==1 || ISIMIP_RUN==2 || ISIMIP_RUN==3 || GEOMIP_RUN==1){
        if((fpi = fopen("./data/rcp_ghg.txt","rt"))==NULL){
            printf("No rcp_co2.txt\n");
            exit(1);
        }
        for(f=0;f<DL_AGHG;f++){
            fscanf(fpi,"%ld", &year); /* 1765-2500 */
            /* CO2, ppmv */
            /* CH4, ppbv */
            /* N2O, ppbv */
            
            /* RCP2.6 */
            fscanf(fpi,"%lf", &aco2_b1[f]);
            fscanf(fpi,"%lf", &ach4_b1[f]);
            fscanf(fpi,"%lf", &an2o_b1[f]);
            /* RCP4.5 */
            fscanf(fpi,"%lf", &aco2_b2[f]);
            fscanf(fpi,"%lf", &ach4_b2[f]);
            fscanf(fpi,"%lf", &an2o_b2[f]);
            /* RCP6.0 */
            fscanf(fpi,"%lf", &aco2_a1[f]);
            fscanf(fpi,"%lf", &ach4_a1[f]);
            fscanf(fpi,"%lf", &an2o_a1[f]);
            /* RCP8.5 */
            fscanf(fpi,"%lf", &aco2_a2[f]);
            fscanf(fpi,"%lf", &ach4_a2[f]);
            fscanf(fpi,"%lf", &an2o_a2[f]);
        }
        fclose(fpi);
    }
	
	/* global analysis initialization ********************************/
	go_landarea = gs_landarea = 0.0;
	for(f=0;f<PD_SIM;f++){
		h_tmp[f] = h_pre[f] = h_dswr[f] = h_aet[f] = h_rof[f] = 0.0;
		h_gpp[f] = h_npp[f] = h_nep[f] = h_plant[f] = h_soil[f] = 0.0;
		h_sr[f] = h_ersn_c[f] = h_agrersn_c[f] = h_doc[f] = 0.0;
		h_gpp_df97[f] = h_gpp_c4[f] = 0.0;
		h_pot_prmfrst[f] = 0.0;
		h_trnsp[f] = h_incepev[f] = h_ssurfev[f] = 0.0;
		h_nbp[f] = h_hvst[f] = h_abgm[f] = 0.0;
		h_sw1[f] = h_sw2[f] = 0.0;
        h_rns[f] = h_rnl[f] = 0.0; /* added by A.Ito (2013/01/02) */
		h_rnsd[f] = h_cld[f] = h_apar[f] = 0.0;
        h_parb[f] = h_pard[f] = 0.0;
        h_arm[f] = 0.0;
        
		h_agrarea[f] = h_paddyarea[f] = h_luc[f] = 0.0;
		h_luc_1[f] = h_luc_2[f] = h_luc_3[f] = 0.0;
		h_burnt_area[f] = 0.0;
		h_bioburn_co2[f] = h_bioburn_co[f] = h_bioburn_ch4[f] = 0.0;
		h_bioburn_nmhc[f] = h_bioburn_oc[f] = h_bioburn_bc[f] = 0.0;
		h_bioburn_nox[f] = h_bioburn_so2[f] = h_bioburn_pm25[f] = 0.0;
		h_bioburn_tpm[f] = h_bioburn_tec[f] = 0.0;
		
		h_ch4ox1[f] = h_ch4ox2[f] = h_ch4ox3[f] = h_ch4ox4[f] = 0.0;
		h_ch4emit_cao_paddy[f] = h_ch4emit_cao_wetland[f] = 0.0;
		h_ch4emit_paddy_wh_diff[f] = h_ch4emit_paddy_wh_plant[f] = 0.0;
		
		h_ch4emit_paddy_wh_ebbl[f] = h_ch4emit_paddy_wh_release[f] = 0.0;
		h_ch4emit_wetland_wh_diff[f] = h_ch4emit_wetland_wh_plant[f] = 0.0;
		h_ch4emit_wetland_wh_ebbl[f] = h_ch4emit_wetland_wh_release[f] = 0.0;

		h_n2o_emit_ngas[f] = h_n2_emit_ngas[f] = 0.0;
		h_n2o_emit_casa[f] = h_no_emit_casa[f] = h_n2_emit_casa[f] = 0.0;
		h_nh3_emit[f] = h_n2_biofix[f] = 0.0;
		h_ch4_emit_mass[f] = h_ch4_emit_photo[f] = 0.0;
		h_n2o_d_emit_ngas[f] = h_n2o_n_emit_ngas[f] = 0.0;
		h_n2o_emit_ngas_agr[f] = h_n2o_emit_casa_agr[f] = 0.0;
		h_nh3_emit_agr[f] = 0.0;
		h_no3_leach[f] = 0.0;
		h_n_fertin[f] = h_n_depoin[f] = 0.0;
		h_hvst_wood[f] = h_wetarea[f] = h_deforest[f] = 0.0;

		h_voc_isopr_g97[f] = h_voc_monotrp_g97[f] = h_voc_methanl_g97[f] = 0.0;
		h_voc_acetone_g97[f] = h_voc_actaldhd_g97[f] = h_voc_frmardhd_g97[f] = 0.0;
		h_voc_formacd_g97[f] = h_voc_acetacd_g97[f] = h_voc_co_g97[f] = 0.0;
		h_voc_afarnesene[f] = h_voc_bcaryophyllene[f] = h_voc_othersesqui[f] = 0.0;
		
		ci_aco2[f] = ci_aco2_d13c[f] = ci_aco2_d14c[f] = 0.0;
		ci_gpp[f] = ci_gpp_d13c[f] = ci_gpp_d14c[f] = 0.0;
		ci_er[f] = ci_er_d13c[f] = ci_er_d14c[f] = 0.0;
		ci_f[f] = ci_f_d13c[f] = ci_f_d14c[f] = 0.0;
		ci_c[f] = ci_c_d13c[f] = ci_c_d14c[f] = 0.0;
		ci_r[f] = ci_r_d13c[f] = ci_r_d14c[f] = 0.0;
		ci_l[f] = ci_l_d13c[f] = ci_l_d14c[f] = 0.0;
		ci_h[f] = ci_h_d13c[f] = ci_h_d14c[f] = 0.0;
		
		for(g=0;g<ASTEP;g++){
			hm_temp[f][g] = 0.0;
			hm_prec[f][g] = 0.0;
			hm_ch4_wh[f][g] = 0.0;
			hm_inund[f][g] = 0.0;
		}
	}
	for(f=0;f<ASTEP;f++){
		m_ch4ox1[f] = m_ch4ox2[f] = m_ch4ox3[f] = 0.0;
		m_bioburn_co2[f] = m_bioburn_ch4[f] = m_bioburn_co[f] = 0.0;
		m_bioburn_nmhc[f] = m_bioburn_oc[f] = m_bioburn_bc[f] = 0.0;
		m_gpp[f] = m_npp[f] = m_nep[f] = 0.0;
		m_ch4p_cao[f] = m_ch4p_wh[f] = 0.0;
	}
	for(f=0;f<NVEG_OLSON;f++){
		vo_area[f] = 0.0;
		vo_gpp[f] = vo_npp[f] = vo_nep[f] = 0.0;
		vo_lai[f] = vo_fol[f] = vo_stm[f] = vo_rot[f] = vo_ltr[f] = vo_msl[f] = 0.0;
	}
	for(f=0;f<NVEG_SAGE;f++){
		vs_area[f] = 0.0;
		vs_gpp[f] = vs_npp[f] = vs_nep[f] = 0.0;
		vs_lai[f] = vs_fol[f] = vs_stm[f] = vs_rot[f] = vs_ltr[f] = vs_msl[f] = 0.0;
	}

	for(g=0;g<N_ROW;g++){
		for(h=0;h<N_COL;h++){
			for(f=0;f<5;f++){
				g_tmp[f][g][h] = 0.0;
				g_prc[f][g][h] = 0.0;
				g_swr[f][g][h] = 0.0;
				g_gpp[f][g][h] = 0.0;
				g_npp[f][g][h] = 0.0;
				g_nep[f][g][h] = 0.0;
				g_pmas[f][g][h] = 0.0;
				g_smas[f][g][h] = 0.0;
				g_ch4e_cao[f][g][h] = 0.0;
				g_ch4ep_cao[f][g][h] = 0.0;
				g_ch4o_curry[f][g][h] = 0.0;
				g_n2oe[f][g][h] = 0.0;
				g_bbco2[f][g][h] = 0.0;
				g_ersn[f][g][h] = 0.0;
				g_isopr[f][g][h] = 0.0;
				g_sr[f][g][h] = 0.0;
				g_luc[f][g][h] = 0.0;
				g_er[f][g][h] = 0.0; 
				g_snh4[f][g][h] = 0.0; 
				g_sno3[f][g][h] = 0.0; 
				
#if C13_GOUT==1				
				g_f13[f][g][h] = 0.0; 
				g_c13[f][g][h] = 0.0; 
				g_r13[f][g][h] = 0.0; 
				g_l13[f][g][h] = 0.0; 
				g_h13[f][g][h] = 0.0; 
				g_gpp13[f][g][h] = 0.0; 
				g_er13[f][g][h] = 0.0; 
#endif
				
#if C14_GOUT==1				
				g_f14[f][g][h] = 0.0; 
				g_c14[f][g][h] = 0.0; 
				g_r14[f][g][h] = 0.0; 
				g_l14[f][g][h] = 0.0; 
				g_h14[f][g][h] = 0.0; 
				g_gpp14[f][g][h] = 0.0; 
				g_er14[f][g][h] = 0.0; 
#endif				
				
#if PHYS_GOUT==1
				g_lai[f][g][h] = 0.0;
				g_parb[f][g][h] = 0.0;
				g_pard[f][g][h] = 0.0;
				g_apar[f][g][h] = 0.0;
				g_apar2[f][g][h] = 0.0;
				g_aet[f][g][h] = 0.0;
				g_rof[f][g][h] = 0.0;
				g_rns[f][g][h] = 0.0;
				g_rnl[f][g][h] = 0.0;
				g_sw1[f][g][h] = 0.0;
				g_sw2[f][g][h] = 0.0;
				g_rnsd[f][g][h] = 0.0;
#endif
			}
		}
	}  /* */
	
#if CH4_WH==1
	for(g=0;g<N_ROW;g++){
		for(h=0;h<N_COL;h++){
			for(f=0;f<5;f++){
				g_ch4ep_wh[f][g][h] = 0.0;
				g_ch4ew_wh[f][g][h] = 0.0;
			}
			for(f=0;f<ASTEP;f++){
				gm_ch4ep_wh[f][g][h] = 0.0;
			}
		}
	}
#endif	
	
	/* regional historical */
	for(f=0;f<N_REG;f++){
		rh_area[f] = 0.0;
		for(g=0;g<PD_SIM;g++){
			rh_temp[f][g] = rh_prec[f][g] = rh_dswrf[f][g] = 0.0;
			rh_rns[f][g] = rh_rnl[f][g] = 0.0;
			rh_ipar[f][g] = rh_apar[f][g] = 0.0;
			rh_gpp[f][g] = rh_npp[f][g] = rh_nep[f][g] = 0.0;
			rh_evpr[f][g] = rh_trsp[f][g] = rh_incp[f][g] = rh_rnof[f][g] = 0.0;
			rh_hvst[f][g] = rh_luc[f][g] = 0.0;
			rh_ch4ox_curry[f][g] = rh_ch4emit_wh_wet[f][g] = rh_ch4emit_wh_paddy[f][g] = 0.0;
			rh_n2o_emit_ngas[f][g] = rh_n2o_emitagr_ngas[f][g] = 0.0;
			
			rh_ci_gpp[f][g] = rh_ci_gpp_d13c[f][g] = rh_ci_gpp_d14c[f][g] = 0.0;
			rh_ci_er[f][g] = rh_ci_er_d13c[f][g] = rh_ci_er_d14c[f][g] = 0.0;
			rh_ci_f[f][g] = rh_ci_f_d13c[f][g] = rh_ci_f_d14c[f][g] = 0.0;
			rh_ci_c[f][g] = rh_ci_c_d13c[f][g] = rh_ci_c_d14c[f][g] = 0.0;
			rh_ci_r[f][g] = rh_ci_r_d13c[f][g] = rh_ci_r_d14c[f][g] = 0.0;
			rh_ci_l[f][g] = rh_ci_l_d13c[f][g] = rh_ci_l_d14c[f][g] = 0.0;
			rh_ci_h[f][g] = rh_ci_h_d13c[f][g] = rh_ci_h_d14c[f][g] = 0.0;
		}
	}
    
    for(h=0;h<N_COL;h++){
        glat_area[h] = 0.0;
    }
	for(g=0;g<ASTEP;g++){
		for(h=0;h<N_COL;h++){
            glat_gpp[g][h] = glat_npp[g][h] = glat_nep[g][h] = 0.0;
            glat_ch4_cao[g][h] = glat_ch4_wh[g][h] = 0.0;
        }
    }
}
