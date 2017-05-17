/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* initialize environmental locations */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

#define TERM_HYD 0.1

/* initialization of climatic conditions (Primary data) ******************************/
void f_init_clim(
	struct Grid *grid
){
	short h, impex_p, impex_t, imp_case;
	double aaa, bbb, tmp_var, pre_var;
   /* IMPRESSIONS phase2: 2017/04/20 by A.Ito */
    double tw_scotland[10][ASTEP] = {
        { -3.1 , -3.22 , -3.34 , -3.34 , -3.22 , -2.99 , -2.83 , -2.73 , -2.71 , -2.75 , -2.85 , -2.98 },
        { -2.07 , -2.15 , -2.23 , -2.23 , -2.15 , -1.99 , -1.89 , -1.82 , -1.8 , -1.83 , -1.89 , -1.99 },
        { -1.03 , -1.07 , -1.11 , -1.11 , -1.07 , -1 , -0.94 , -0.91 , -0.9 , -0.91 , -0.95 , -0.99 },
        { 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 },
        { 0.48 , 0.46 , 0.44 , 0.44 , 0.47 , 0.51 , 0.53 , 0.55 , 0.55 , 0.54 , 0.52 , 0.5 },
        { 0.97 , 0.93 , 0.89 , 0.89 , 0.93 , 1 , 1.06 , 1.09 , 1.1 , 1.09 , 1.05 , 1.01 },
        { 1.93 , 1.85 , 1.77 , 1.77 , 1.85 , 2.01 , 2.11 , 2.18 , 2.2 , 2.17 , 2.11 , 2.01 },
        { 2.9 , 2.78 , 2.66 , 2.66 , 2.78 , 3.01 , 3.17 , 3.27 , 3.29 , 3.25 , 3.15 , 3.02 },
        { 3.87 , 3.71 , 3.55 , 3.55 , 3.7 , 4.01 , 4.23 , 4.36 , 4.39 , 4.34 , 4.2 , 4.04 },
        { 4.84 , 4.63 , 4.44 , 4.44 , 4.63 , 5.02 , 5.29 , 5.45 , 5.49 , 5.42 , 5.26 , 5.05 }};
    double pw_scotland[10][ASTEP] = {
        { -1.4 , -2.59 , -5.85 , -10.2 , -15.72 , -22.3 , -23.88 , -20.39 , -11.85 , -6 , -2.84 , -2.3 },
        { 1.6 , 0.41 , -2.85 , -7.2 , -12.72 , -19.3 , -20.88 , -17.39 , -8.85 , -3 , 0.16 , 0.7 },
        { 4.6 , 3.41 , 0.15 , -4.2 , -9.72 , -16.3 , -17.88 , -14.39 , -5.85 , 0 , 3.16 , 3.7 },
        { 10.6 , 9.41 , 6.15 , 1.8 , -3.72 , -10.3 , -11.88 , -8.39 , 0.15 , 6 , 9.16 , 9.7 },
        { 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 },
        { 13.6 , 12.41 , 9.15 , 4.8 , -0.72 , -7.3 , -8.88 , -5.39 , 3.15 , 9 , 12.16 , 12.7 },
        { 16.6 , 15.41 , 12.15 , 7.8 , 2.28 , -4.3 , -5.88 , -2.39 , 6.15 , 12 , 15.16 , 15.7 },
        { 19.6 , 18.41 , 15.15 , 10.8 , 5.28 , -1.3 , -2.88 , 0.61 , 9.15 , 15 , 18.16 , 18.7 },
        { 22.6 , 21.41 , 18.15 , 13.8 , 8.28 , 1.7 , 0.12 , 3.61 , 12.15 , 18 , 21.16 , 21.7 },
        { 25.6 , 24.41 , 21.15 , 16.8 , 11.28 , 4.7 , 3.12 , 6.61 , 15.15 , 21 , 24.16 , 24.7 }};
    double tw_iberia[10][ASTEP] = {
        { -1.28 , -1.29 , -1.23 , -1.12 , -0.96 , -0.76 , -0.66 , -0.68 , -0.8 , -0.92 , -1.05 , -1.18 },
        { 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 },
        { 0.37 , 0.36 , 0.39 , 0.44 , 0.52 , 0.62 , 0.67 , 0.66 , 0.6 , 0.54 , 0.48 , 0.41 },
        { 0.72 , 0.71 , 0.77 , 0.88 , 1.04 , 1.24 , 1.34 , 1.32 , 1.2 , 1.08 , 0.95 , 0.82 },
        { 1.46 , 1.43 , 1.55 , 1.77 , 2.08 , 2.49 , 2.68 , 2.66 , 2.4 , 2.16 , 1.91 , 1.64 },
        { 2.18 , 2.14 , 2.33 , 2.66 , 3.13 , 3.74 , 4.02 , 3.98 , 3.6 , 3.24 , 2.86 , 2.45 },
        { 2.92 , 2.86 , 3.11 , 3.55 , 4.17 , 4.99 , 5.36 , 5.31 , 4.8 , 4.32 , 3.82 , 3.28 },
        { 3.64 , 3.57 , 3.89 , 4.43 , 5.21 , 6.23 , 6.7 , 6.64 , 6 , 5.42 , 4.78 , 4.1 },
        { 4.37 , 4.28 , 4.66 , 5.32 , 6.26 , 7.48 , 8.05 , 7.97 , 7.2 , 6.5 , 5.73 , 4.91 },
        { 5.1 , 5 , 5.44 , 6.21 , 7.3 , 8.72 , 9.38 , 9.29 , 8.4 , 7.58 , 6.69 , 5.74 }};
    double pw_iberia[11][ASTEP] = {
        { -34.8 , -36.0 , -40.3 , -44.1 , -47.37 , -50.2 , -51.1 , -50.2 , -47.29 , -44.3 , -41.0 , -36.61 },
        { -28.8 , -30.0 , -34.3 , -38.1 , -41.37 , -44.2 , -45.1 , -44.2 , -41.29 , -38.3 , -35.0 , -30.61 },
        { -22.8 , -24.0 , -28.3 , -32.1 , -35.37 , -38.2 , -39.1 , -38.2 , -35.29 , -32.3 , -29.0 , -24.61 },
        { -16.8 , -18.0 , -22.3 , -26.1 , -29.37 , -32.2 , -33.1 , -32.2 , -29.29 , -26.3 , -23.0 , -18.61 },
        { -10.8 , -12.0 , -16.3 , -20.1 , -23.37 , -26.2 , -27.1 , -26.2 , -23.29 , -20.3 , -17.0 , -12.61 },
        { -4.8 , -6 , -10.3 , -14.1 , -17.37 , -20.2 , -21.1 , -20.2 , -17.29 , -14.3 , -11 , -6.61 },
        { 1.2 , 0.0 , -4.3 , -8.1 , -11.37 , -14.2 , -15.1 , -14.2 , -11.29 , -8.3 , -5 , -0.61 },
        { 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 },
        { 13.2 , 12.0 , 7.7 , 3.9 , 0.63 , -2.2 , -3.1 , -2.2 , 0.71 , 3.7 , 7.0 , 11.39 },
        { 19.2 , 18.0 , 13.7 , 9.9 , 6.63 , 3.8 , 2.9 , 3.8 , 6.71 , 9.7 , 13.0 , 17.39 },
        { 25.2 , 24.0 , 19.7 , 15.9 , 12.63 , 9.8 , 8.9 , 9.8 , 12.71 , 15.7 , 19.0 , 23.39 }};
    double tw_hungary[11][ASTEP] = {
        { -1.05 , -1.12 , -1.22 , -1.22 , -1.12 , -0.91 , -0.81 , -0.82 , -0.94 , -1.02 , -1.06 , -1.01 },
        { 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 },
        { 0.48 , 0.44 , 0.39 , 0.39 , 0.45 , 0.55 , 0.6 , 0.59 , 0.52 , 0.48 , 0.46 , 0.49 },
        { 0.95 , 0.88 , 0.78 , 0.78 , 0.88 , 1.09 , 1.19 , 1.18 , 1.06 , 0.98 , 0.94 , 0.99 },
        { 1.9 , 1.76 , 1.55 , 1.55 , 1.76 , 2.18 , 2.37 , 2.35 , 2.1 , 1.94 , 1.87 , 1.97 },
        { 2.87 , 2.65 , 2.33 , 2.32 , 2.64 , 3.27 , 3.57 , 3.53 , 3.16 , 2.92 , 2.8 , 2.97 },
        { 3.82 , 3.54 , 3.1 , 3.1 , 3.52 , 4.36 , 4.76 , 4.71 , 4.21 , 3.88 , 3.73 , 3.96 },
        { 4.77 , 4.42 , 3.88 , 3.87 , 4.4 , 5.45 , 5.95 , 5.88 , 5.26 , 4.86 , 4.67 , 4.94 },
        { 5.72 , 5.3 , 4.66 , 4.65 , 5.28 , 6.54 , 7.13 , 7.06 , 6.31 , 5.82 , 5.6 , 5.93 },
        { 6.67 , 6.18 , 5.43 , 5.42 , 6.16 , 7.63 , 8.32 , 8.23 , 7.37 , 6.8 , 6.53 , 6.92 },
        { 7.64 , 7.07 , 6.21 , 6.2 , 7.03 , 8.72 , 9.51 , 9.41 , 8.41 , 7.76 , 7.46 , 7.91 }};
    double pw_hungary[9][ASTEP] = {
        { -21.7 , -22 , -24.18 , -28.7 , -35.59 , -44.8 , -47.8 , -44.6 , -35.2 , -28.5 , -24.49 , -23.5 },
        { -16.7 , -17 , -19.18 , -23.7 , -30.59 , -39.8 , -42.8 , -39.6 , -30.2 , -23.5 , -19.49 , -18.5 },
        { -11.7 , -12 , -14.18 , -18.7 , -25.59 , -34.8 , -37.8 , -34.6 , -25.2 , -18.5 , -14.49 , -13.5 },
        { -6.7 , -7.0 , -9.18 , -13.7 , -20.59 , -29.8 , -32.8 , -29.6 , -20.2 , -13.5 , -9.49 , -8.5 },
        { -1.7 , -2.0 , -4.18 , -8.7 , -15.59 , -24.8 , -27.8 , -24.6 , -15.2 , -8.5 , -4.49 , -3.5 },
        { 3.3 , 3.0 , 0.82 , -3.7 , -10.59 , -19.8 , -22.8 , -19.6 , -10.2 , -3.5 , 0.51 , 1.5 },
        { 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 },
        { 13.3 , 13 , 10.82 , 6.3 , -0.59 , -9.8 , -12.8 , -9.6 , -0.2 , 6.5 , 10.51 , 11.5 },
        { 18.3 , 18 , 15.82 , 11.3 , 4.41 , -4.8 , -7.8 , -4.6 , 4.8 , 11.5 , 15.51 , 16.5 }};
	
	/* in 1950 :311 ppmv*/
	/* in 1990 : 352.7 ppmv*/
	grid->co2y = BGY_CO2Y; 
	if(CO2S == 7){
		grid->co2y = 2081; /* in 2081 : 700 ppmv*/
	}
    if(CC_CD == 5){
        grid->co2y = 2000;
    }
	grid->climy = BGY_CLIM;
	grid->lucy = BGY_CLIM;

	for(h=0;h<ASTEP;h++){
		grid->m = h;
		
		/* ambient CO2 condition *******/
		f_co2_trend(grid);

		/* climatic conditions */
		/* NCEP/NCAR mean climate data ***************************************/
		grid->tmp_sfc[h] = grid->tmp_sfc_a[h]; /* temperature, degree C */
		grid->tmp_2m[h] = grid->tmp_2m_a[h];
		grid->tmp10_soil[h] = grid->tmp10_soil_a[h];
		grid->tmp200_soil[h] = grid->tmp200_soil_a[h];
		grid->tcdc_clm[h] = grid->tcdc_clm_a[h]; /* cloudiness, 0-1 fraction */
		grid->prate_sfc[h] = grid->prate_sfc_a[h]; /* precipitation, mm / month */
		
		/*** UEA/CRU data ***/
		if(grid->flag_histdata == 1){ /* for grids data are available */
			/* temporary */
			grid->tmp_sfc[h] = grid->hist_tmp_b[h] + (grid->tmp_sfc_a[h] - grid->tmp_2m_a[h]);
			grid->tmp10_soil[h] = grid->hist_tmp_b[h] + (grid->tmp10_soil_a[h] - grid->tmp_2m_a[h]);
			grid->tmp200_soil[h] = grid->hist_tmp_b[h] + (grid->tmp200_soil_a[h] - grid->tmp_2m_a[h]);
			grid->tmp_2m[h] = grid->hist_tmp_b[h];
			grid->tcdc_clm[h] = grid->hist_cld_b[h];
			grid->prate_sfc[h] = grid->hist_pre_b[h];
			
			/* long-term average */
			grid->tmp_sfc_a[h] = grid->hist_tmp_b[h] + (grid->tmp_sfc_a[h] - grid->tmp_2m_a[h]);
			grid->tmp10_soil_a[h] = grid->hist_tmp_b[h] + (grid->tmp10_soil_a[h] - grid->tmp_2m_a[h]);
			grid->tmp200_soil_a[h] = grid->hist_tmp_b[h] + (grid->tmp200_soil_a[h] - grid->tmp_2m_a[h]);
			grid->tmp_2m_a[h] = grid->hist_tmp_b[h];
			grid->tcdc_clm_a[h] = grid->hist_cld_b[h];
			grid->prate_sfc_a[h] = grid->hist_pre_b[h];
		}
		
		/* alternative precipitation data *************************/
		if(grid->prec_sub_a[h] >= 0.0){
			grid->prate_sfc[h] = grid->prec_sub_a[h];
		}else if(grid->prec_sub_a[h] < 0.0){
			grid->prate_sfc[h] = grid->prate_sfc_a[h];
		}
		
		grid->spfh_2m[h] = grid->spfh_2m_a[h];
		
		aaa = grid->ugrd_10m_a[h] * grid->ugrd_10m_a[h];
		bbb = grid->vgrd_10m_a[h] * grid->vgrd_10m_a[h];
		grid->wnd_10m[h] = sqrt(aaa + bbb);
        
        /*******************************************/
        /* IMPRESSIONS IRS: 2015/07/17 by A.Ito */
        if(IMPRESSIONS_RUN == 1){
        
            if(SCENARIO_ID<6001 || SCENARIO_ID>7000){
                printf("BAD experimental ID\n");
                exit(1);
            }
            
            /*******/
            impex_t = (short)(SCENARIO_ID%50);
            
            switch(impex_t){
                case 1: tmp_var = -3.0; break;
                case 2: tmp_var = -2.0; break;
                case 3: tmp_var = -1.0; break;
                case 4: tmp_var = 0.0; break;
                case 5: tmp_var = 1.0; break;
                case 6: tmp_var = 2.0; break;
                case 7: tmp_var = 3.0; break;
                case 8: tmp_var = 4.0; break;
                case 9: tmp_var = 5.0; break;
                case 10: tmp_var = 6.0; break;
                case 11: tmp_var = 7.0; break;
                case 12: tmp_var = 8.0; break;
                case 13: tmp_var = 9.0; break;
                case 14: tmp_var = 10.0; break;
                case 15: tmp_var = 11.0; break;
                default: tmp_var = 0.0; break;
            }
            
            grid->tmp_sfc[h] += tmp_var;
            grid->tmp_2m[h] += tmp_var;
            grid->tmp10_soil[h] += tmp_var;
            grid->tmp200_soil[h] += tmp_var;
            
            /*******/
            impex_p = (short)((SCENARIO_ID - 6000)/50);
            
            switch(impex_p){
                case 0: pre_var = 0.4; break;
                case 1: pre_var = 0.5; break;
                case 2: pre_var = 0.6; break;
                case 3: pre_var = 0.7; break;
                case 4: pre_var = 0.8; break;
                case 5: pre_var = 0.9; break;
                case 6: pre_var = 1.0; break;
                case 7: pre_var = 1.1; break;
                case 8: pre_var = 1.2; break;
                case 9: pre_var = 1.3; break;
                case 10: pre_var = 1.4; break;
                case 11: pre_var = 1.5; break;
                case 12: pre_var = 1.6; break;
                default: pre_var = 1.0; break;
            }
            
            grid->prate_sfc[h] *= pre_var;
        }else if(IMPRESSIONS_RUN == 2 || IMPRESSIONS_RUN == 3){
            /* IMPRESSIONS 2: 2017/05/02 by A.Ito */
            
            if(SCENARIO_ID<60000 || SCENARIO_ID>70000){
                printf("BAD experimental ID\n");
                exit(1);
            }
            
            impex_t = (short)(SCENARIO_ID%100);
            
            tmp_var = 0.0;
            switch(impex_t){
                case 0: tmp_var = -3.0; break;
                case 1: tmp_var = -2.0; break;
                case 2: tmp_var = -1.0; break;
                case 3: tmp_var = 0.0; break;
                case 4: tmp_var = 0.5; break;
                case 5: tmp_var = 1.0; break;
                case 6: tmp_var = 2.0; break;
                case 7: tmp_var = 3.0; break;
                case 8: tmp_var = 4.0; break;
                case 9: tmp_var = 5.0; break;
                case 10: tmp_var = 6.0; break;
                case 11: tmp_var = 7.0; break;
                case 12: tmp_var = 8.0; break;
            }
            
            /* Hungary */
            if(grid->impressions_mask == 1){
                if(impex_t == 2){
                    imp_case = 0;
                }else if(impex_t == 3){
                    imp_case = 1;
                }else if(impex_t == 4){
                    imp_case = 2;
                }else if(impex_t == 5){
                    imp_case = 3;
                }else if(impex_t == 6){
                    imp_case = 4;
                }else if(impex_t == 7){
                    imp_case = 5;
                }else if(impex_t == 8){
                    imp_case = 6;
                }else if(impex_t == 9){
                    imp_case = 7;
                }else if(impex_t == 10){
                    imp_case = 8;
                }else if(impex_t == 11){
                    imp_case = 9;
                }else if(impex_t == 12){
                    imp_case = 10;
                }else{
                    imp_case = -999;
                }
                
                for(h=0;h<ASTEP;h++){
                    if(imp_case>=0){ tmp_var = tw_hungary[imp_case][h]; }
                
                    grid->tmp_sfc[h] += tmp_var;
                    grid->tmp_2m[h] += tmp_var;
                    grid->tmp10_soil[h] += tmp_var;
                    grid->tmp200_soil[h] += tmp_var;
                }
            }
            
            /* Iberia */
            if(grid->impressions_mask == 2){
                if(impex_t == 2){
                    imp_case = 0;
                }else if(impex_t == 3){
                    imp_case = 1;
                }else if(impex_t == 4){
                    imp_case = 2;
                }else if(impex_t == 5){
                    imp_case = 3;
                }else if(impex_t == 6){
                    imp_case = 4;
                }else if(impex_t == 7){
                    imp_case = 5;
                }else if(impex_t == 8){
                    imp_case = 6;
                }else if(impex_t == 9){
                    imp_case = 7;
                }else if(impex_t == 10){
                    imp_case = 8;
                }else if(impex_t == 11){
                    imp_case = 9;
                }else{
                    imp_case = -999;
                }
                
                for(h=0;h<ASTEP;h++){
                    if(imp_case>=0){ tmp_var = tw_iberia[imp_case][h]; }
                
                    grid->tmp_sfc[h] += tmp_var;
                    grid->tmp_2m[h] += tmp_var;
                    grid->tmp10_soil[h] += tmp_var;
                    grid->tmp200_soil[h] += tmp_var;
                }
            }
            
            /* Scotland */
            if(grid->impressions_mask == 3){
                if(impex_t == 0){
                    imp_case = 0;
                }else if(impex_t == 1){
                    imp_case = 1;
                }else if(impex_t == 2){
                    imp_case = 2;
                }else if(impex_t == 3){
                    imp_case = 3;
                }else if(impex_t == 4){
                    imp_case = 4;
                }else if(impex_t == 5){
                    imp_case = 5;
                }else if(impex_t == 6){
                    imp_case = 6;
                }else if(impex_t == 7){
                    imp_case = 7;
                }else if(impex_t == 8){
                    imp_case = 8;
                }else if(impex_t == 9){
                    imp_case = 9;
                }else{
                    imp_case = -999;
                }
                
                for(h=0;h<ASTEP;h++){
                    if(imp_case>=0){ tmp_var = tw_scotland[imp_case][h]; }
                
                    grid->tmp_sfc[h] += tmp_var;
                    grid->tmp_2m[h] += tmp_var;
                    grid->tmp10_soil[h] += tmp_var;
                    grid->tmp200_soil[h] += tmp_var;
                }
            }
            
            impex_p = (short)((SCENARIO_ID - 60000)/100);
            
            pre_var = 1.0;
            switch(impex_p){
                case 0: pre_var = 0.58; break;
                case 1: pre_var = 0.64; break;
                case 2: pre_var = 0.70; break;
                case 3: pre_var = 0.75; break;
                case 4: pre_var = 0.76; break;
                case 5: pre_var = 0.80; break;
                case 6: pre_var = 0.82; break;
                case 7: pre_var = 0.85; break;
                case 8: pre_var = 0.88; break;
                case 9: pre_var = 0.90; break;
                case 10: pre_var = 0.91; break;
                case 11: pre_var = 0.94; break;
                case 12: pre_var = 0.95; break;
                case 13: pre_var = 0.97; break;
                case 14: pre_var = 1.00; break;
                case 15: pre_var = 1.03; break;
                case 16: pre_var = 1.05; break;
                case 17: pre_var = 1.06; break;
                case 18: pre_var = 1.09; break;
                case 19: pre_var = 1.10; break;
                case 20: pre_var = 1.12; break;
                case 21: pre_var = 1.15; break;
                case 22: pre_var = 1.18; break;
            }
            
            /* Hungary */
            if(grid->impressions_mask == 1){
                
                if(impex_t == 2){
                    imp_case = 0;
                }else if(impex_t == 3){
                    imp_case = 1;
                }else if(impex_t == 5){
                    imp_case = 2;
                }else if(impex_t == 7){
                    imp_case = 3;
                }else if(impex_t == 9){
                    imp_case = 4;
                }else if(impex_t == 12){
                    imp_case = 5;
                }else if(impex_t == 14){
                    imp_case = 6;
                }else if(impex_t == 16){
                    imp_case = 7;
                }else if(impex_t == 19){
                    imp_case = 8;
                }else{
                    imp_case = -999;
                }
                
                for(h=0;h<ASTEP;h++){
                    if(imp_case>=0){ pre_var =  (100.0 + pw_hungary[imp_case][h])/100.0; }
                
                    grid->prate_sfc[h] *= pre_var;
                }
            }

            /* Iberia */
            if(grid->impressions_mask == 2){
                
                if(impex_t == 0){
                    imp_case = 0;
                }else if(impex_t == 1){
                    imp_case = 1;
                }else if(impex_t == 2){
                    imp_case = 2;
                }else if(impex_t == 4){
                    imp_case = 3;
                }else if(impex_t == 6){
                    imp_case = 4;
                }else if(impex_t == 8){
                    imp_case = 5;
                }else if(impex_t == 11){
                    imp_case = 6;
                }else if(impex_t == 14){
                    imp_case = 7;
                }else if(impex_t == 17){
                    imp_case = 8;
                }else if(impex_t == 20){
                    imp_case = 9;
                }else if(impex_t == 22){
                    imp_case = 10;
                }else{
                    imp_case = -999;
                }
                
                for(h=0;h<ASTEP;h++){
                    if(imp_case>=0){ pre_var =  (100.0 + pw_iberia[imp_case][h])/100.0; }
                
                    grid->prate_sfc[h] *= pre_var;
                }
            }

            /* Scotland */
            if(grid->impressions_mask == 3){
                
                if(impex_t == 10){
                    imp_case = 0;
                }else if(impex_t == 11){
                    imp_case = 1;
                }else if(impex_t == 13){
                    imp_case = 2;
                }else if(impex_t == 14){
                    imp_case = 3;
                }else if(impex_t == 15){
                    imp_case = 4;
                }else if(impex_t == 17){
                    imp_case = 5;
                }else if(impex_t == 18){
                    imp_case = 6;
                }else if(impex_t == 20){
                    imp_case = 7;
                }else if(impex_t == 21){
                    imp_case = 8;
                }else if(impex_t == 22){
                    imp_case = 9;
                }else{
                    imp_case = -999;
                }
                
                for(h=0;h<ASTEP;h++){
                    if(imp_case>=0){ pre_var =  (100.0 + pw_scotland[imp_case][h])/100.0; }
                
                    grid->prate_sfc[h] *= pre_var;
                }
            }
        }

        grid->tmp_soil_am += grid->tmp200_soil_a[h] * MDN[h] / 365.0;

		/* solar decrination and solar height */
		grid->sl_dec[h] = f_solar_decl(grid);
		grid->sl_hgt[h] = f_solar_hgt(grid);
		
		/* radiation fluxes and day-length */
		grid->dlen[h] = f_day_length(grid);		
		grid->top_rad[h] = f_top_rad(grid, 0); 	

		grid->gl_rad[h] = f_gl_rad(grid); 	
		grid->par[h] = f_par(grid); 
		
		grid->rad_a[h] = grid->gl_rad[h];
		grid->par_a[h] = grid->par[h];
        
		/* sensitivity analysis *************************/
		if(TM == 1){
			grid->tmp_sfc[h] += 1.0;
			grid->tmp_2m[h] += 1.0;
			grid->tmp10_soil[h] += 1.0;
			grid->tmp200_soil[h] += 1.0;
		}else if(TM == 5){
			grid->tmp_sfc[h] -= 1.0;
			grid->tmp_2m[h] -= 1.0;
			grid->tmp10_soil[h] -= 1.0;
			grid->tmp200_soil[h] -= 1.0;
		}
		if(PR == 1){
			grid->prate_sfc[h] *= 1.1;
		}else if(PR == 5){
			grid->prate_sfc[h] *= 0.9;
		}
	}
}

/* location conditions derived from the primary data (Secondary data1) *******************/
void f_init_loct(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Flux *flux, 
	struct Echar *echar
){
	long h, nn;
	double tem_grow, d_smc_a;
	double ftmp10b, ftmp200b, ftmp10, ftmp200;
	extern double MDN[12];
	
	/* long-term average annual climate condition */
	tem_grow = 0.0;
	grid->tmp_sfc_am = 0.0;
	grid->tmp_sfc_mx = -100.0;
	grid->tmp_sfc_mn = 100.0;
	grid->gp_atem = 0.0;
	grid->prate_sfc_ann = 0.0;
	nn = 0;
	for(h=0;h<ASTEP;h++){
		/* annual mean temperature */
		grid->tmp_sfc_am += grid->tmp_sfc[h]/12.0; 
		/* annual maximum */
		grid->tmp_sfc_mx = (grid->tmp_sfc[h]>grid->tmp_sfc_mx)?grid->tmp_sfc[h]:grid->tmp_sfc_mx; 
		/* annual minimum */
		grid->tmp_sfc_mn = (grid->tmp_sfc[h]<grid->tmp_sfc_mn)?grid->tmp_sfc[h]:grid->tmp_sfc_mn; 
		/* annual total precipitation */
		grid->prate_sfc_ann += grid->prate_sfc[h]; 
		
		if(grid->tmp_sfc[h]>5.0){
			nn += MDN[h];
			grid->gp_atem += (grid->tmp_sfc[h]-5.0)*MDN[h];
			tem_grow += grid->tmp_sfc[h]*MDN[h];
			grid->gp_pre += grid->prate_sfc[h]; 
		}
	}
	grid->gp_tem = (nn>0)?tem_grow/(double)nn:0.0;

	/* portion of C3 and C4 plant ***************/
	c34composition((echar->c3).v_type, grid, loct);
	
	/* sensitivity analysis */
	if(T_D==1 || T_D==2 || T_D==3 || T_D==4){
		ftmp10b = ftmp200b = ftmp10 = ftmp200 = 0.0;
		for(h=0;h<ASTEP;h++){
			if(grid->tmp10_soil[h] > -20.0){
				ftmp10b += 0.05 + 0.95*exp(308.56*(1.0/56.02-1.0/(grid->tmp10_soil[h]+46.02))); 
				if(T_D==0){
					ftmp10 += 0.05 + 0.95*exp(308.56*(1.0/56.02-1.0/(grid->tmp10_soil[h]+46.02)));
				}else if(T_D==1){
					ftmp10 += 0.05 + 0.95*exp(308.56*1.3*(1.0/56.02-1.0/(grid->tmp10_soil[h]+46.02)));
				}else if(T_D==2){
					ftmp10 += 0.05 + 0.95*exp(308.56*0.7*(1.0/56.02-1.0/(grid->tmp10_soil[h]+46.02)));
				}else if(T_D==3){
					ftmp10 += 0.05 + 0.95*exp(308.56*(1.0/56.02-1.0/(grid->tmp10_soil[h]+46.02*1.3)));
				}else if(T_D==4){
					ftmp10 += 0.05 + 0.95*exp(308.56*(1.0/56.02-1.0/(grid->tmp10_soil[h]+46.02*0.7)));
				}
			}else{
				ftmp10b += 0.05;
				ftmp10 += 0.05;
			}
			if(grid->tmp200_soil[h] > -20.0){
				ftmp200b += 0.05 + 0.95*exp(308.56*(1.0/56.02-1.0/(grid->tmp200_soil[h]+46.02)));
				if(T_D==0){
					ftmp200 += 0.05 + 0.95*exp(308.56*(1.0/56.02-1.0/(grid->tmp200_soil[h]+46.02)));
				}else if(T_D==1){
					ftmp200 += 0.05 + 0.95*exp(308.56*1.3*(1.0/56.02-1.0/(grid->tmp200_soil[h]+46.02)));
				}else if(T_D==2){
					ftmp200 += 0.05 + 0.95*exp(308.56*0.7*(1.0/56.02-1.0/(grid->tmp200_soil[h]+46.02)));
				}else if(T_D==3){
					ftmp200 += 0.05 + 0.95*exp(308.56*(1.0/56.02-1.0/(grid->tmp200_soil[h]+46.02*1.3)));
				}else if(T_D==4){
					ftmp200 += 0.05 + 0.95*exp(308.56*(1.0/56.02-1.0/(grid->tmp200_soil[h]+46.02*0.7)));
				}
			}else{
				ftmp200b += 0.05;
				ftmp200 += 0.05;
			}
		}
		(echar->soil).rl0 *= ftmp10b/ftmp10;
		(echar->soil).rh0 *= ftmp200b/ftmp200;
	}
	
	/* water condition - Sim-HYDRO *****************/
	loct->sw30 = 0.5*grid->field_cap1;
	loct->sww = 0.5*grid->field_cap2; 
	loct->snwa = 0.0;
	loct->m_m_pre = 0.24;
	nn = 0; 
	d_smc_a = 10.0;
	while(d_smc_a > TERM_HYD){
		d_smc_a = loct->sww;
		for(h=0;h<ASTEP;h++){
			grid->m = h;
		/*	(mass->c3).lai[grid->m]=loct->c3ptn[grid->m];
			(mass->c4).lai[grid->m]=loct->c4ptn[grid->m]; */
			
			co2_in_canopy(grid, loct, mass, flux);
			
			f_dyn_loct(grid, loct, mass, echar);
		}
		d_smc_a = fabs(loct->sww - d_smc_a);
		loct->time_hyd = nn; /* simulation time of carbon budget */
		nn++;
		
		if(nn < 6){
			d_smc_a = 10.0; /* at least 5 years */
		}
		if(nn > 50){
			break;
		}
	}
	
	/* average fertilizer-N input for each county, kg N ha-1 yr-1 */
	n_fertilizer_in(grid, loct);
    
    /* empirical NPP and meta-analysis-based optimal LAI */
    npp_empirical(grid, loct, flux);
}

/* dynamic estimation of environmnetal conditions (Secondary data2) *********************/
void f_dyn_loct(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Echar *echar
){
	long h, offset;
	double alt, k_c, vpres_var, tmp_ann, tmp_var, wet_var;
	
    /* initialization */
    if(grid->m == 0){
        for(h=0;h<ASTEP;h++){
            loct->xx1[h] = 0.0;
            loct->xx2[h] = 0.0;
            loct->xx3[h] = 0.0;
            loct->xx4[h] = 0.0;
            loct->xx5[h] = 0.0;
            loct->xx6[h] = 0.0;
            loct->xx7[h] = 0.0;
            loct->xx8[h] = 0.0;
            loct->xx9[h] = 0.0;
        }
        loct->est_maxlai = 0.0;
    }
    for(h=0;h<DSTEP;h++){
        loct->ppfd_h[h] = 0.0;
        loct->ppfdb_h[h] = 0.0;
        loct->ppfdd_h[h] = 0.0;
    }

	/* solar constant sensitivity */
	if(SC==3 || SC==4){
		grid->top_rad[grid->m] = f_top_rad(grid, 0); 	
		grid->gl_rad[grid->m] = f_gl_rad(grid); 	
		grid->par[grid->m] = f_par(grid); 
	}else if(SC == 5){
		grid->par[grid->m] = f_par(grid); 
		grid->par[grid->m] += 10.0;
	}
    
	/* radiatin for cal_historical: 1901-2000 */
	if(grid->flag_histdata == 1){
    
        /* daily-mean and hourly: added: 2013/01/10 by A.Ito */
        loct->grad_d[grid->m] = 0.0;
        for(h=0;h<DSTEP;h++){
            grid->h = h;
        
            /* atmosphere-top, hourly */
            grid->top_rad[grid->m] = f_top_rad(grid, -180 + h*15);
            /* surface, hourly */
            grid->gl_rad[grid->m] = f_gl_rad(grid);
            /* daily mean */
            loct->grad_d[grid->m] += f_gl_rad(grid)/(double)DSTEP;
            
            /* hourly PAR, micro mol photon m-2 s-1 */
            f_par_h(grid, loct);
        }
        
        /* midday */
        grid->top_rad[grid->m] = f_top_rad(grid, 0);
		grid->gl_rad[grid->m] = f_gl_rad(grid); 
		grid->par[grid->m] = f_par(grid);
	}
	
	if(grid->m == 0){
		grid->tmp_sfc_am = 0.0;
		grid->tmp_sfc_mx = -100.0;
		grid->tmp_sfc_mn = 100.0;
		grid->prate_sfc_ann = 0.0;
		grid->tmp_soil_mean = 0.0;
		for(h=0;h<ASTEP;h++){
			/* annual mean temperature */
			grid->tmp_sfc_am += grid->tmp_sfc[h] * MDN[h]/365.0;
			/* annual maximum */
			grid->tmp_sfc_mx = (grid->tmp_sfc[h]>grid->tmp_sfc_mx)?grid->tmp_sfc[h]:grid->tmp_sfc_mx; 
			/* annual minimum */
			grid->tmp_sfc_mn = (grid->tmp_sfc[h]<grid->tmp_sfc_mn)?grid->tmp_sfc[h]:grid->tmp_sfc_mn; 
			/* annual total precipitation */
			grid->prate_sfc_ann += grid->prate_sfc[h]; 
			/* annual mean soil temperature */
			grid->tmp_soil_mean += grid->tmp10_soil[h] * MDN[h]/365.0;
		}
		loct->cum_dprec = 0.0;
        
        /* change in wetland area due to permafrost melting: 2012/10/26 by A.Ito ****/
        if(VAR_PFMWET == 1){
            
            tmp_ann = 0.0;
            for(h=0;h<ASTEP;h++){
                tmp_ann = grid->tmp_2m[h]*MDN[h]/365.0;
            }
            
            tmp_var = tmp_ann - grid->tmp_base_permaforst;
            
            wet_var = 0.0;
            if(tmp_ann > -2.0 && tmp_var > 0.0 ){
                switch(grid->type_permaforst){
                    case 1: case 5: case 9: case 13: case 17:
                        /* continuous permafrost */
                        wet_var = 0.05 * tmp_var;
                        break;
                    case 2: case 6: case 10: case 14: case 18:
                        /* discontinuous permafrost */
                        wet_var = 0.03 * tmp_var;
                        break;
                    case 3: case 7: case 11: case 15: case 19:
                        /* sporadic permafrost */
                        wet_var = 0.01 * tmp_var;
                        break;
                    case 4: case 8: case 12: case 16: case 20:
                        /* isolated permafrost */
                        wet_var = 0.003 * tmp_var;
                        break;
                    default:
                        /* no permafrost */
                        wet_var = 0.0;
               }
            }
            if(wet_var < 0.0){
                wet_var = 0.0;
            }
            
            grid->f_wetland += wet_var;
            
            if(grid->f_wetland > 1.0){
                grid->f_wetland = 1.0;
            }
        }
	}
	/* cumulative precipitation anomaly: 2010/07/23 by A.Ito */
	/* if(grid->climy<2005){ */
		loct->cum_dprec += (grid->prate_sfc[grid->m] - grid->hist_pre_b[grid->m]);
	/* } */
	
	/** air conditions **/
	/* altitude */
	alt = (grid->topo>=0.0)?grid->topo:0.0; 
	/* air pressure */
	loct->prsr[grid->m] = 1013.25 * exp(-1.0 * (28.964 * 0.001) * GAC * alt / (UGC * (grid->tmp_2m[grid->m] + ZAT)));
	
	/* saturated vapour pressure, hPa */
	loct->vps[grid->m] = vap_pre_sat(grid); 
	/* slope of saturated vapour pressure as a function of temperature */
	loct->slope_vps[grid->m] = slope_vps(grid);	
	 /* aerodynamic resistance */ 
	loct->r_aero[grid->m] = r_aero(grid);	 

	/* initial soil CH4 concentration */
	for(h=0;h<=(N_SLAYER+1);h++){
		loct->prof_ch4[h] = ach4_a1[grid->co2y - 1750]/1000.0 
			* loct->prsr[grid->m] / (UGC * (grid->tmp10_soil[grid->m] + ZAT));
	}
    
    offset = 0;
    if(ISIMIP_RUN == 3){
        if(grid->phase == 0){
            offset = 0;
        }else if(grid->phase == 1){
            /* skip spin-up data */
            offset = 30;
        }
    }

	if(grid->flag_histdata == 1){
		/* vapour pressure, hPa */
		if(grid->phase == 0){
			/* spin-up */
			loct->vp[grid->m] = grid->hist_vap_b[grid->m];
		}else if(grid->phase == 1){
			if(grid->climy < (BGY_CLIM + DL_HCLIM)){
				/* based on UEA/CRU or ISI-MIP data */
				loct->vp[grid->m] = grid->hist_vap[grid->climy - BGY_CLIM + offset][grid->m];
			}else{
				/* based on NCEP/NCAR */
				vpres_var = grid->ncep_vpres[grid->climy - FDY_NCEP][grid->m][grid->ncep_lat][grid->ncep_lon] 
								- grid->ncep_vpres_b[grid->m][grid->ncep_lat][grid->ncep_lon];
				
				loct->vp[grid->m] = grid->hist_vap_b[grid->m] + vpres_var;
			}
		}else if(grid->phase == 2){
			/* prediction using AOGCM */
			/* loct->vp[grid->m] = loct->prsr[grid->m]*grid->spfh_2m[grid->m]/(0.622 + 0.378*grid->spfh_2m[grid->m]);  */
			
			/* revided by A.Ito (2009/08/17) */
			vpres_var = grid->proj_hum[grid->climy - FDY_GCM-1][grid->m][grid->gcm_row][grid->gcm_col] -
							grid->proj_hum_b[grid->m][grid->gcm_row][grid->gcm_col];
			
			loct->vp[grid->m] = grid->hist_vap_b[grid->m] + vpres_var;
            
            if(ISIMIP_RUN == 2){
                loct->vp[grid->m] = grid->proj_hum[grid->climy - FDY_GCM][grid->m][0][0];
            }
		}
		if(loct->vp[grid->m] < 0.0){
			loct->vp[grid->m] = 0.0;
		}
		
		if(loct->vp[grid->m]<=loct->vps[grid->m]){
			loct->vpd[grid->m] = loct->vps[grid->m] - loct->vp[grid->m];
		}else{
			loct->vpd[grid->m] = 0.0;
		}
	}else{
		/* vapour pressure, hPa */
		loct->vp[grid->m] = loct->prsr[grid->m] * grid->spfh_2m[grid->m]/(0.622 + 0.378*grid->spfh_2m[grid->m]);
		/* vapour pressure deficit */
		loct->vpd[grid->m] = (loct->vps[grid->m] >= loct->vp[grid->m])?loct->vps[grid->m] - loct->vp[grid->m]:0.0;
	}

	/* air density */
	loct->dnsa[grid->m] = air_density(grid, loct); 
	
	/** ecophysiology : ecophysiology.c **/
	/* C3 plants */
	f_ecophysiology(grid, loct, &(echar->c3), &(mass->c3));
	/* C4 plants */
	f_ecophysiology(grid, loct, &(echar->c4), &(mass->c4));
	
	/** net radiation **/
	f_net_rad(grid, loct, mass, echar);
    	
	(mass->plant).lai[grid->m] = (mass->c3).lai[grid->m] * loct->c3ptn[grid->m]
					+ (mass->c4).lai[grid->m] * loct->c4ptn[grid->m];
	loct->lai[grid->m] = (mass->plant).lai[grid->m];
    
	/** hydrological water budget **/
	f_waterbudget(grid, loct, echar);
	
	/* fractional vegetation cover */
	/* light attenuation by non-photosynthetic organs */
	/*
	Kurachi, N., A. Hagihara, and K. Hozumi. 1989. Effect of light interception by 
	non-photosynthetic organs on canopy photosynthetic production. 
	Ecological Research 4:187-197.
	*/
	if(grid->veg_sage >= 1 && grid->veg_sage <= 8){
		k_c = 0.221;	/* Kurachi and Hagihara */
	}else{
		k_c = 0.001;
	}

	loct->f_vegcov[grid->m] = 1.0 - loct->c3ptn[grid->m] * exp(-((echar->c3).eK0 + k_c)*(mass->c3).lai[grid->m]) 
								- loct->c4ptn[grid->m] * exp(-((echar->c4).eK0 + 0.001)*(mass->c4).lai[grid->m]);
	if(loct->f_vegcov[grid->m]<0.0){
		loct->f_vegcov[grid->m] = 0.0;
	}
	if(loct->f_vegcov[grid->m]>1.0){
		loct->f_vegcov[grid->m] = 1.0;
	}
	
	/* CASA water scheme **************/
	loct->m_sw[grid->m] = loct->sw30;
	loct->m_vmc[grid->m] = loct->sw30/300.0;
	loct->m_pet[grid->m] = loct->pm_incep[grid->m] + loct->pm_evp[grid->m] + loct->pm_trn[grid->m];
	
	loct->wfps[grid->m] = ((loct->m_vmc[grid->m]*100.0) / ((1.0 - grid->bulkdens/2.65)*100.0));
	
	if(loct->wfps[grid->m] > 1.5){
		loct->wfps[grid->m] = 1.5;
	}
	if(loct->wfps[grid->m] < 0.05){
		loct->wfps[grid->m] = 0.05;
	}
	/* soil moisture index */
	f_casa_moisture(grid, loct);
	
	/* nitrogen deposition ***********/
	f_n_deposit(grid, loct);
	
	/* decay of 14C: added by A.Ito (2009/06/24) ************************/
	(mass->c3).d14c_fol = f_decay_14c((mass->c3).d14c_fol);
	(mass->c3).d14c_stm = f_decay_14c((mass->c3).d14c_stm);
	(mass->c3).d14c_rot = f_decay_14c((mass->c3).d14c_rot);
	(mass->c4).d14c_fol = f_decay_14c((mass->c4).d14c_fol);
	(mass->c4).d14c_stm = f_decay_14c((mass->c4).d14c_stm);
	(mass->c4).d14c_rot = f_decay_14c((mass->c4).d14c_rot);
	(mass->soil).d14c_ltr = f_decay_14c((mass->soil).d14c_ltr);
	(mass->soil).d14c_msl = f_decay_14c((mass->soil).d14c_msl);
}
