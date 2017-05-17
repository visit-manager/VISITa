/*	VISIT: Vegetation Integrative SImulation Tool						*/
/*  Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/*  Carbon cycle, erosion, biomass burning, land-use change,			*/
/*  CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

extern short CC_R;
extern short TEMP_GC;

/**************************************************************************/
void set_hist_clim(
	struct Grid *grid
){
    short offset = 0, impex_t, impex_p, imp_case;
	long h, cru_te;
	double tmp_var, pre_var, tcdc_var;
    
    /* IMPRESSIONS phase2: 2017/04/20 by A.Ito */
    double tw_scotland[10][ASTEP] = {
        { -3.1 , -3.22 , -3.34 , -3.34 , -3.22 , -2.99 , -2.83 , -2.73 , -2.71 , -2.75 , -2.85 , -2.98 },
        { -2.07 , -2.15 , -2.23 , -2.23 , -2.15 , -1.99 , -1.89 , -1.82 , -1.8 , -1.83 , -1.89 , -1.99 },
        { -1.03 , -1.07 , -1.11 , -1.11 , -1.07 , -1 , -0.94 , -0.91 , -0.9 , -0.91 , -0.95 , -0.99 },
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
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
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
        { 13.6 , 12.41 , 9.15 , 4.8 , -0.72 , -7.3 , -8.88 , -5.39 , 3.15 , 9 , 12.16 , 12.7 },
        { 16.6 , 15.41 , 12.15 , 7.8 , 2.28 , -4.3 , -5.88 , -2.39 , 6.15 , 12 , 15.16 , 15.7 },
        { 19.6 , 18.41 , 15.15 , 10.8 , 5.28 , -1.3 , -2.88 , 0.61 , 9.15 , 15 , 18.16 , 18.7 },
        { 22.6 , 21.41 , 18.15 , 13.8 , 8.28 , 1.7 , 0.12 , 3.61 , 12.15 , 18 , 21.16 , 21.7 },
        { 25.6 , 24.41 , 21.15 , 16.8 , 11.28 , 4.7 , 3.12 , 6.61 , 15.15 , 21 , 24.16 , 24.7 }};
    double tw_iberia[10][ASTEP] = {
        { -1.28 , -1.29 , -1.23 , -1.12 , -0.96 , -0.76 , -0.66 , -0.68 , -0.8 , -0.92 , -1.05 , -1.18 },
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
        { 0.37 , 0.36 , 0.39 , 0.44 , 0.52 , 0.62 , 0.67 , 0.66 , 0.6 , 0.54 , 0.48 , 0.41 },
        { 0.72 , 0.71 , 0.77 , 0.88 , 1.04 , 1.24 , 1.34 , 1.32 , 1.2 , 1.08 , 0.95 , 0.82 },
        { 1.46 , 1.43 , 1.55 , 1.77 , 2.08 , 2.49 , 2.68 , 2.66 , 2.4 , 2.16 , 1.91 , 1.64 },
        { 2.18 , 2.14 , 2.33 , 2.66 , 3.13 , 3.74 , 4.02 , 3.98 , 3.6 , 3.24 , 2.86 , 2.45 },
        { 2.92 , 2.86 , 3.11 , 3.55 , 4.17 , 4.99 , 5.36 , 5.31 , 4.8 , 4.32 , 3.82 , 3.28 },
        { 3.64 , 3.57 , 3.89 , 4.43 , 5.21 , 6.23 , 6.7 , 6.64 , 6 , 5.42 , 4.78 , 4.1 },
        { 4.37 , 4.28 , 4.66 , 5.32 , 6.26 , 7.48 , 8.05 , 7.97 , 7.2 , 6.5 , 5.73 , 4.91 },
        { 5.1 , 5 , 5.44 , 6.21 , 7.3 , 8.72 , 9.38 , 9.29 , 8.4 , 7.58 , 6.69 , 5.74 }};
    double pw_iberia[11][ASTEP] = {
        { -34.8 , -36 , -40.3 , -44.1 , -47.37 , -50.2 , -51.1 , -50.2 , -47.29 , -44.3 , -41 , -36.61 },
        { -28.8 , -30 , -34.3 , -38.1 , -41.37 , -44.2 , -45.1 , -44.2 , -41.29 , -38.3 , -35 , -30.61 },
        { -22.8 , -24 , -28.3 , -32.1 , -35.37 , -38.2 , -39.1 , -38.2 , -35.29 , -32.3 , -29 , -24.61 },
        { -16.8 , -18 , -22.3 , -26.1 , -29.37 , -32.2 , -33.1 , -32.2 , -29.29 , -26.3 , -23 , -18.61 },
        { -10.8 , -12 , -16.3 , -20.1 , -23.37 , -26.2 , -27.1 , -26.2 , -23.29 , -20.3 , -17 , -12.61 },
        { -4.8 , -6 , -10.3 , -14.1 , -17.37 , -20.2 , -21.1 , -20.2 , -17.29 , -14.3 , -11 , -6.61 },
        { 1.2 , 0 , -4.3 , -8.1 , -11.37 , -14.2 , -15.1 , -14.2 , -11.29 , -8.3 , -5 , -0.61 },
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
        { 13.2 , 12 , 7.7 , 3.9 , 0.63 , -2.2 , -3.1 , -2.2 , 0.71 , 3.7 , 7 , 11.39 },
        { 19.2 , 18 , 13.7 , 9.9 , 6.63 , 3.8 , 2.9 , 3.8 , 6.71 , 9.7 , 13 , 17.39 },
        { 25.2 , 24 , 19.7 , 15.9 , 12.63 , 9.8 , 8.9 , 9.8 , 12.71 , 15.7 , 19 , 23.39 }};
    double tw_hungary[11][ASTEP] = {
        { -1.05 , -1.12 , -1.22 , -1.22 , -1.12 , -0.91 , -0.81 , -0.82 , -0.94 , -1.02 , -1.06 , -1.01 },
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
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
        { -6.7 , -7 , -9.18 , -13.7 , -20.59 , -29.8 , -32.8 , -29.6 , -20.2 , -13.5 , -9.49 , -8.5 },
        { -1.7 , -2 , -4.18 , -8.7 , -15.59 , -24.8 , -27.8 , -24.6 , -15.2 , -8.5 , -4.49 , -3.5 },
        { 3.3 , 3 , 0.82 , -3.7 , -10.59 , -19.8 , -22.8 , -19.6 , -10.2 , -3.5 , 0.51 , 1.5 },
        { 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 , 0 },
        { 13.3 , 13 , 10.82 , 6.3 , -0.59 , -9.8 , -12.8 , -9.6 , -0.2 , 6.5 , 10.51 , 11.5 },
        { 18.3 , 18 , 15.82 , 11.3 , 4.41 , -4.8 , -7.8 , -4.6 , 4.8 , 11.5 , 15.51 , 16.5 }};
	
	/* last year of CRU-data calculation */
	cru_te = DL_HCLIM + (BGY_CLIM-1);
    
    impex_t = impex_p = -999;
    
    if(grid->climy < BGY_CLIM){
        printf("Bad climate-data year: %ld\n", grid->climy);
        exit(1);
    }
	
    if(ISIMIP_RUN == 0){
        if(grid->climy <= cru_te){
            /* 1901-2000:CRU TS2.1 (20th century) */
            /* 1901-2002:CRU TS2.1 */
            /* 1901-2005:CRU TS3.0 */
            /* 1901-2009:CRU TS3.1 */
            /* 1901-2011:CRU TS3.2 */
            /* 1901-2012:CRU TS3.21 */
            /* 1901-2013:CRU TS3.22 */
            /* 1901-2014:CRU TS3.23 */
            /* 1901-2015:CRU TS3.24 */
            for(h=0;h<ASTEP;h++){
                grid->tmp_sfc[h] = grid->hist_tmp[grid->climy - BGY_CLIM][h] 
                                + (grid->tmp_sfc_a[h] - grid->tmp_2m_a[h]);
                grid->tmp_2m[h] = grid->hist_tmp[grid->climy - BGY_CLIM][h];
                grid->tmp10_soil[h] = grid->hist_tmp[grid->climy - BGY_CLIM][h] 
                                + (grid->tmp10_soil_a[h] - grid->tmp_2m_a[h]);
                grid->tmp200_soil[h] = grid->hist_tmp[grid->climy - BGY_CLIM][h] 
                                + (grid->tmp200_soil_a[h] - grid->tmp_2m_a[h]);
                grid->tcdc_clm[h] = grid->hist_cld[grid->climy - BGY_CLIM][h];
                grid->prate_sfc[h] = grid->hist_pre[grid->climy - BGY_CLIM][h];
                
           }
        }else{
            /* extention by NCEP/NCAR reanalysis data */
            /* 2003-2008: extrapolation using NCEP/NCAR data: 2009/01/05 by A.Ito */
            /* 2006-2009: extrapolation using NCEP/NCAR data: 2010/01/04 by A.Ito */
            /* 2006-2010: extrapolation using NCEP/NCAR data: 2011/03/XX by A.Ito */
            /* 2012-2012: extrapolation using NCEP/NCAR data: 2012/04/14 by A.Ito */
            /* 2013-2013: extrapolation using NCEP/NCAR data: 2013/01/10 by A.Ito */
            /* 2014-2014: extrapolation using NCEP/NCAR data: 2015/01/05 by A.Ito */
            /* 2015-2015: extrapolation using NCEP/NCAR data: 2016/01/04 by A.Ito */
            for(h=0;h<ASTEP;h++){
                /* temperature */
                tmp_var = grid->ncep_tmp2m[grid->climy - FDY_NCEP][h][grid->ncep_lat][grid->ncep_lon] 
                        - grid->ncep_tmp2m_b[h][grid->ncep_lat][grid->ncep_lon];
                grid->tmp_sfc[h] = grid->tmp_sfc_a[h] + tmp_var;
                grid->tmp_2m[h] = grid->tmp_2m_a[h] + tmp_var;
                grid->tmp10_soil[h] = grid->tmp10_soil_a[h] + tmp_var*0.3;
                grid->tmp200_soil[h] = grid->tmp200_soil_a[h] + tmp_var*0.1;
                
                /* precipitation */
                pre_var = grid->ncep_prate[grid->climy - FDY_NCEP][h][grid->ncep_lat][grid->ncep_lon] 
                            - grid->ncep_prate_b[h][grid->ncep_lat][grid->ncep_lon];
                grid->prate_sfc[h] = grid->prate_sfc_a[h] + pre_var;
                if(grid->prate_sfc[h]<0.0){
                    grid->prate_sfc[h] = 0.0;
                }
                
                /* cloudiness */
                tcdc_var = grid->ncep_tcdc[grid->climy - FDY_NCEP][h][grid->ncep_lat][grid->ncep_lon] 
                        - grid->ncep_tcdc_b[h][grid->ncep_lat][grid->ncep_lon];
                grid->tcdc_clm[h] = grid->tcdc_clm_a[h] + tcdc_var;
                if(grid->tcdc_clm[h] < 0.0){
                    grid->tcdc_clm[h] = 0.0;
                }
                if(grid->tcdc_clm[h] > 1.0){
                    grid->tcdc_clm[h] = 1.0;
                }
            }
        }
        
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
            
            for(h=0;h<ASTEP;h++){
                grid->tmp_sfc[h] += tmp_var;
                grid->tmp_2m[h] += tmp_var;
                grid->tmp10_soil[h] += tmp_var;
                grid->tmp200_soil[h] += tmp_var;
            }
            
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
            
            for(h=0;h<ASTEP;h++){
                grid->prate_sfc[h] *= pre_var;
            }
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
        
    }else if(ISIMIP_RUN == 1){
        
        /* ISI-MIP climate data: 2012/06/28 by A.Ito */
        offset = 0;
        if(grid->phase == 0){
            offset = 0;
        }else if(grid->phase == 1 || grid->phase == 2){
            /* skip spin-up data */
            offset = 30;
        }
        
        for(h=0;h<ASTEP;h++){
            grid->tmp_sfc[h] = grid->hist_tmp[grid->climy - BGY_CLIM + offset][h]
                            + (grid->tmp_sfc_a[h] - grid->tmp_2m_a[h]);
            grid->tmp_2m[h] = grid->hist_tmp[grid->climy - BGY_CLIM + offset][h];
            grid->tmp10_soil[h] = grid->hist_tmp[grid->climy - BGY_CLIM + offset][h] 
                            + (grid->tmp10_soil_a[h] - grid->tmp_2m_a[h]);
            grid->tmp200_soil[h] = grid->hist_tmp[grid->climy - BGY_CLIM + offset][h] 
                            + (grid->tmp200_soil_a[h] - grid->tmp_2m_a[h]);
            grid->tcdc_clm[h] = grid->hist_cld[grid->climy - BGY_CLIM + offset][h];
            grid->prate_sfc[h] = grid->hist_pre[grid->climy - BGY_CLIM + offset][h];  
        }
    }else if(ISIMIP_RUN == 2){

        /* ISI-MIP climate data: 2014/07/31 by A.Ito */
        if(grid->phase == 0){
            offset = 0;
        }else if(grid->phase == 1 || grid->phase == 2){
            offset = 30;
        }
        
        for(h=0;h<ASTEP;h++){
            grid->tmp_sfc[h] = grid->hist_tmp[grid->climy - BGY_CLIM + offset][h] 
                            + (grid->tmp_sfc_a[h] - grid->tmp_2m_a[h]);
            grid->tmp_2m[h] = grid->hist_tmp[grid->climy - BGY_CLIM + offset][h];
            grid->tmp10_soil[h] = grid->hist_tmp[grid->climy - BGY_CLIM + offset][h] 
                            + (grid->tmp10_soil_a[h] - grid->tmp_2m_a[h]);
            grid->tmp200_soil[h] = grid->hist_tmp[grid->climy - BGY_CLIM + offset][h] 
                            + (grid->tmp200_soil_a[h] - grid->tmp_2m_a[h]);
            grid->tcdc_clm[h] = grid->hist_cld[grid->climy - BGY_CLIM + offset][h];
            grid->prate_sfc[h] = grid->hist_pre[grid->climy - BGY_CLIM + offset][h];  
        }
    }else if(ISIMIP_RUN == 3){
        
        /* ISI-MIP2 climate data: 2014/11/30 by A.Ito */
        if(grid->phase == 0){
            offset = 0;
        }else if(grid->phase == 1 || grid->phase == 2){
            /* skip spin-up data */
            offset = 30;
        }
        
        for(h=0;h<ASTEP;h++){
            grid->tmp_sfc[h] = grid->hist_tmp[grid->climy - BGY_CLIM + offset][h] 
                            + (grid->tmp_sfc_a[h] - grid->tmp_2m_a[h]);
            grid->tmp_2m[h] = grid->hist_tmp[grid->climy - BGY_CLIM + offset][h];
            grid->tmp10_soil[h] = grid->hist_tmp[grid->climy - BGY_CLIM + offset][h] 
                            + (grid->tmp10_soil_a[h] - grid->tmp_2m_a[h]);
            grid->tmp200_soil[h] = grid->hist_tmp[grid->climy - BGY_CLIM + offset][h] 
                            + (grid->tmp200_soil_a[h] - grid->tmp_2m_a[h]);
            grid->tcdc_clm[h] = grid->hist_cld[grid->climy - BGY_CLIM + offset][h];
            grid->prate_sfc[h] = grid->hist_pre[grid->climy - BGY_CLIM + offset][h];  
        }
    }else if(ISIMIP_RUN == 4){
        
        /* ISI-MIP2b climate data: 2016/12/24 by A.Ito */
        if(grid->phase == 0){
            offset = 0;
        }else if(grid->phase == 1 || grid->phase == 2){
            /* skip spin-up data */
            offset = 0;
        }
        
        for(h=0;h<ASTEP;h++){
            grid->tmp_sfc[h] = grid->hist_tmp[grid->climy - BGY_CLIM + offset][h] 
                            + (grid->tmp_sfc_a[h] - grid->tmp_2m_a[h]);
            grid->tmp_2m[h] = grid->hist_tmp[grid->climy - BGY_CLIM + offset][h];
            grid->tmp10_soil[h] = grid->hist_tmp[grid->climy - BGY_CLIM + offset][h] 
                            + (grid->tmp10_soil_a[h] - grid->tmp_2m_a[h]);
            grid->tmp200_soil[h] = grid->hist_tmp[grid->climy - BGY_CLIM + offset][h] 
                            + (grid->tmp200_soil_a[h] - grid->tmp_2m_a[h]);
            grid->tcdc_clm[h] = grid->hist_cld[grid->climy - BGY_CLIM + offset][h];
            grid->prate_sfc[h] = grid->hist_pre[grid->climy - BGY_CLIM + offset][h];  
        }
    }
	
	/* perturbation for uncertainty analysis: 2010/05/17 (A.Ito) ***************/
	if(PRT_CLIM == 1){
		if(PARAM_PTB == 1){
			for(h=0;h<ASTEP;h++){
				grid->tmp_sfc[h] += 0.2*f_pert[11];
				grid->tmp_2m[h] += 0.2*f_pert[11];
				grid->tmp10_soil[h] += 0.2*f_pert[11];
				grid->tmp200_soil[h] += 0.2*f_pert[11];
				
				grid->prate_sfc[h] *= 1.0 + 0.1*f_pert[12];
				grid->tcdc_clm[h] *= 1.0 + 0.1*f_pert[13];
				if(grid->tcdc_clm[h]<0.0){
					grid->tcdc_clm[h]=0.0;
				}
				if(grid->tcdc_clm[h]>1.0){
					grid->tcdc_clm[h]=1.0;
				}
			}
		}
	}
}

/* prepare GCM-derived climate condition ******************************/
void set_gcm_clim(
	struct Grid *grid
){
	long h;
	double tmp_var, pre_var, rad_var;
	double alt, apres, rvap, est_cld;
    
    /***************************************************/
    for(h=0;h<ASTEP;h++){
        grid->m = h;
        
        /* temperature ***************************************************/
        /* gradual temperature change: added by A.Ito (2009/06/15) */
        if(TEMP_GC == 1){
            tmp_var = -1.0/100.0 * (double)(grid->climy-2000);
        }else if(TEMP_GC == 2){
            tmp_var = 1.0/100.0 * (double)(grid->climy-2000);
        }else if(TEMP_GC == 3){
            tmp_var = -2.0/100.0 * (double)(grid->climy-2000);
        }else if(TEMP_GC == 4){
            tmp_var = 2.0/100.0 * (double)(grid->climy-2000);
        }else if(TEMP_GC == 5){
            tmp_var = -3.0/100.0 * (double)(grid->climy-2000);
        }else if(TEMP_GC == 6){
            tmp_var = 3.0/100.0 * (double)(grid->climy-2000);
        }else if(TEMP_GC == 7){
            tmp_var = -4.0/100.0 * (double)(grid->climy-2000);
        }else if(TEMP_GC == 8){
            tmp_var = 4.0/100.0 * (double)(grid->climy-2000);
        }else if(TEMP_GC == 9){
            tmp_var = -5.0/100.0 * (double)(grid->climy-2000);
        }else if(TEMP_GC == 10){
            tmp_var = 5.0/100.0 * (double)(grid->climy-2000);
        }else if(TEMP_GC == 11){
            tmp_var = -6.0/100.0 * (double)(grid->climy-2000);
        }else if(TEMP_GC == 12){
            tmp_var = 6.0/100.0 * (double)(grid->climy-2000);
        }else if(TEMP_GC ==0){
            if(CC_T == 1){
                tmp_var = grid->proj_tmp2m[grid->climy - FDY_GCM][h][grid->gcm_row][grid->gcm_col] -
                            grid->proj_tmp2m_b[h][grid->gcm_row][grid->gcm_col];
            }else{
                tmp_var = 0.0;
            }
        }else{
            tmp_var = 0.0;
        }
        
        if(CC_T_A == 1){
            grid->tmp_sfc[h] = grid->tmp_sfc_a[h] + tmp_var;
            grid->tmp_2m[h] = grid->tmp_2m_a[h] + tmp_var;
        }else{
            grid->tmp_sfc[h] = grid->tmp_sfc_a[h];
            grid->tmp_2m[h] = grid->tmp_2m_a[h];
        }
        if(CC_T_B == 1){
            grid->tmp10_soil[h] = grid->tmp10_soil_a[h] + tmp_var;
            grid->tmp200_soil[h] = grid->tmp200_soil_a[h] + tmp_var;
        }else{
            grid->tmp10_soil[h] = grid->tmp10_soil_a[h];
            grid->tmp200_soil[h] = grid->tmp200_soil_a[h];
        }
        
        if(ISIMIP_RUN == 2){
            /* PLUME: 2014/07/31 by A.Ito */
            grid->tmp_2m[h] = grid->proj_tmp2m[grid->climy - FDY_GCM][h][grid->gcm_row][grid->gcm_col];
            grid->tmp_sfc[h] = grid->proj_tmp2m[grid->climy - FDY_GCM][h][grid->gcm_row][grid->gcm_col]
                                + (grid->tmp_sfc_a[h] - grid->tmp_2m_a[h]);
            grid->tmp10_soil[h] = grid->proj_tmp2m[grid->climy - FDY_GCM][h][grid->gcm_row][grid->gcm_col]
                                + (grid->tmp10_soil_a[h] - grid->tmp_2m_a[h]);
            grid->tmp200_soil[h] = grid->proj_tmp2m[grid->climy - FDY_GCM][h][grid->gcm_row][grid->gcm_col]
                                + (grid->tmp200_soil_a[h] - grid->tmp_2m_a[h]);
       }
        
        /* precipitation ***************************************************/
        if(TEMP_GC != 0){
            pre_var = 0.0;
        }else{
            if(CC_P == 1){
                pre_var = grid->proj_prec[grid->climy - FDY_GCM][h][grid->gcm_row][grid->gcm_col] -
                            grid->proj_prec_b[h][grid->gcm_row][grid->gcm_col];
            }else{
                pre_var = 0.0;
            }
        }
        grid->prate_sfc[h] = grid->prate_sfc_a[h] + pre_var;
       
        if(ISIMIP_RUN == 2){
            /* PLUME: 2014/07/31 by A.Ito */
            grid->prate_sfc[h] = grid->proj_prec[grid->climy-FDY_GCM][h][grid->gcm_row][grid->gcm_col];
        }
        
        if(grid->prate_sfc[h] <= 0.0){
            /* carry over of negative precipitation */
            grid->proj_prec_co += grid->prate_sfc[h];
            grid->prate_sfc[h] = 0.0;
        }else if(grid->prate_sfc[h] > 0.0){
            grid->prate_sfc[h] += grid->proj_prec_co;
            grid->proj_prec_co = 0.0;
            if(grid->prate_sfc[h] <= 0.0){
                /* carry over of negative precipitation */
                grid->proj_prec_co += grid->prate_sfc[h];
                grid->prate_sfc[h] = 0.0;
            }
        }						

        /* air humidity ******************************************************************/
        /* revised by A.Ito (2009/08/17) */
        if(TEMP_GC != 0){
            grid->spfh_2m[h] = grid->spfh_2m_a[h];
        }else{
            alt = (grid->topo>=0.0)?grid->topo:0.0; 
            
            /* atmospheric pressure, hPa */
            apres = 1013.25*exp(-1.0 * (28.964*0.001) * GAC * alt/(UGC * (grid->tmp_2m[h] + ZAT)));
            
            if(CC_H == 1){
                rvap = grid->proj_hum[grid->climy - FDY_GCM][h][grid->gcm_row][grid->gcm_col] / apres;
                grid->spfh_2m[h] = 0.622 * rvap / (1.0 - 0.378 * rvap);
                
                /* shm_var = grid->proj_hum[grid->climy-FDY_GCM-1][h][grid->gcm_row][grid->gcm_col] - 
                            grid->proj_hum_b[h][grid->gcm_row][grid->gcm_col]; */
            }else{
                rvap = grid->proj_hum_b[h][grid->gcm_row][grid->gcm_col] / apres;
                /* shm_var = 0.0; */
            }
            grid->spfh_2m[h] = 0.622 * rvap / (1.0 - 0.378 * rvap);
        }
        
        /* grid->spfh_2m[h] = grid->spfh_2m_a[h] + shm_var; */
        if(grid->spfh_2m[h] < 0.0){
            grid->spfh_2m[h] = 0.0;
        }

        /* surface shortwave radiation ***************************************************/
        /*  0: no change  */
        /*  1: SW and PAR change  */
        /*  2: PAR change but mean SW  */
        /*  3: SW change but mean PAR  */
        if(TEMP_GC != 0){
            rad_var = 0.0;
        }else{
            if(CC_R == 1 || CC_R == 2){
                rad_var = grid->proj_rad[grid->climy-FDY_GCM][h][grid->gcm_row][grid->gcm_col] - 
                            grid->proj_rad_b[h][grid->gcm_row][grid->gcm_col];
            }else if(CC_R == 0 || CC_R == 3){
                rad_var = 0.0;	/* mean SW & PAR */
            }else{
                rad_var = 0.0;	
            }
        }

        /* experiment for SRM by reflector */
        /* added: 2014/07/06 by A.Ito     */
        if(SCENARIO_ID == 3313 || SCENARIO_ID == 3913){
            grid->top_rad[h] = f_top_rad(grid, 0);
        }
        
        grid->gl_rad[h] = grid->rad_a[h] + rad_var;
        
        if(ISIMIP_RUN == 2){
            /* PLUME: 2014/07/31 by A.Ito */
            grid->gl_rad[h] = grid->proj_rad[grid->climy-FDY_GCM][h][grid->gcm_row][grid->gcm_col];
        }
        
        if(grid->gl_rad[h]<0.0){
            grid->gl_rad[h] = 0.0;
        }
        
        /* estimation of cloudiness: 2013/12/03 by A.Ito */
        if(grid->top_rad[h] > 0.0){
            est_cld = grid->gl_rad[h] / grid->top_rad[h];
            if(est_cld < 0.01){
                est_cld = 0.01;
            }else if(est_cld >= 0.8964){
                est_cld = 0.8964;
            }
            
            est_cld = (0.8964 - est_cld)/0.5392;
            if(est_cld < 0.01){
                est_cld = 0.01;
            }else if(est_cld >= 0.99){
                est_cld = 0.99;
            }
        }else{
            est_cld = 0.5;
        }
        grid->tcdc_clm[h] = est_cld;
        
        /* PAR */
        grid->par[h] = f_par(grid);
        
        if(CC_R == 2){
            grid->gl_rad[h] = grid->rad_a[h];	/* mean SW */
        }
        if(CC_R == 3){
            rad_var = grid->proj_rad[grid->climy-FDY_GCM-1][h][grid->gcm_row][grid->gcm_col] - 
                    grid->proj_rad_b[h][grid->gcm_row][grid->gcm_col];
                    
            grid->gl_rad[h] = grid->rad_a[h] + rad_var;
            if(grid->gl_rad[h]<0.0){
                grid->gl_rad[h] = 0.0;	/* SW change */
            }
        }
    }
}

