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
    short offset, impex_t, impex_p;
	long h, cru_te;
	double tmp_var, pre_var, tcdc_var;
	
	/* last year of CRU-data calculation */
	cru_te = DL_CRU + (PIVOT_CLIMY-1);
    
    impex_t = impex_p = -999;
	
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
            for(h=0;h<ASTEP;h++){
                grid->tmp_sfc[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY][h] 
                                + (grid->tmp_sfc_a[h] - grid->tmp_2m_a[h]);
                grid->tmp_2m[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY][h];
                grid->tmp10_soil[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY][h] 
                                + (grid->tmp10_soil_a[h] - grid->tmp_2m_a[h]);
                grid->tmp200_soil[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY][h] 
                                + (grid->tmp200_soil_a[h] - grid->tmp_2m_a[h]);
                grid->tcdc_clm[h] = grid->hist_cld[grid->climy - PIVOT_CLIMY][h];
                grid->prate_sfc[h] = grid->hist_pre[grid->climy - PIVOT_CLIMY][h];  
            }
        }else{
            /* extention by NCEP/NCAR data */
            /* 2003-2008: extrapolation using NCEP/NCAR data: 2009/01/05 by A.Ito */
            /* 2006-2009: extrapolation using NCEP/NCAR data: 2010/01/04 by A.Ito */
            /* 2006-2010: extrapolation using NCEP/NCAR data: 2011/03/XX by A.Ito */
            /* 2012-2012: extrapolation using NCEP/NCAR data: 2012/04/14 by A.Ito */
            /* 2013-2013: extrapolation using NCEP/NCAR data: 2013/01/10 by A.Ito */
            /* 2014-2014: extrapolation using NCEP/NCAR data: 2015/01/05 by A.Ito */
            for(h=0;h<ASTEP;h++){
                /* temperature */
                tmp_var = grid->ncep_tmp2m[grid->climy - PIVOT_NCEP][h][grid->ncep_lat][grid->ncep_lon] 
                        - grid->ncep_tmp2m_b[h][grid->ncep_lat][grid->ncep_lon];
                grid->tmp_sfc[h] = grid->tmp_sfc_a[h] + tmp_var;
                grid->tmp_2m[h] = grid->tmp_2m_a[h] + tmp_var;
                grid->tmp10_soil[h] = grid->tmp10_soil_a[h] + tmp_var*0.3;
                grid->tmp200_soil[h] = grid->tmp200_soil_a[h] + tmp_var*0.1;
                
                /* precipitation */
                pre_var = grid->ncep_prate[grid->climy - PIVOT_NCEP][h][grid->ncep_lat][grid->ncep_lon] 
                            - grid->ncep_prate_b[h][grid->ncep_lat][grid->ncep_lon];
                grid->prate_sfc[h] = grid->prate_sfc_a[h] + pre_var;
                if(grid->prate_sfc[h]<0.0){
                    grid->prate_sfc[h] = 0.0;
                }
                
                /* cloudiness */
                tcdc_var = grid->ncep_tcdc[grid->climy - PIVOT_NCEP][h][grid->ncep_lat][grid->ncep_lon] 
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
        
            if(GCM_ID<6001 || GCM_ID>7000){
                printf("BAD experimental ID\n");
                exit(1);
            }
            
            /*******/
            impex_t = (short)(GCM_ID%50);
            
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
            impex_p = (short)((GCM_ID - 6000)/50);
            
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
        }
        
    }else if(ISIMIP_RUN == 1){
        
        /* ISI-MIP climate data: 2012/06/28 by A.Ito */
        if(grid->phase == 0){
            offset = 0;
        }else if(grid->phase == 1 || grid->phase == 2){
            /* skip spin-up data */
            offset = 30;
        }
        
        for(h=0;h<ASTEP;h++){
            grid->tmp_sfc[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY + offset][h] 
                            + (grid->tmp_sfc_a[h] - grid->tmp_2m_a[h]);
            grid->tmp_2m[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY + offset][h];
            grid->tmp10_soil[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY + offset][h] 
                            + (grid->tmp10_soil_a[h] - grid->tmp_2m_a[h]);
            grid->tmp200_soil[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY + offset][h] 
                            + (grid->tmp200_soil_a[h] - grid->tmp_2m_a[h]);
            grid->tcdc_clm[h] = grid->hist_cld[grid->climy - PIVOT_CLIMY + offset][h];
            grid->prate_sfc[h] = grid->hist_pre[grid->climy - PIVOT_CLIMY + offset][h];  
        }
    }else if(ISIMIP_RUN == 2){

        /* ISI-MIP climate data: 2014/07/31 by A.Ito */
        if(grid->phase == 0){
            offset = 0;
        }else if(grid->phase == 1 || grid->phase == 2){
            offset = 30;
        }
        
        for(h=0;h<ASTEP;h++){
            grid->tmp_sfc[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY + offset][h] 
                            + (grid->tmp_sfc_a[h] - grid->tmp_2m_a[h]);
            grid->tmp_2m[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY + offset][h];
            grid->tmp10_soil[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY + offset][h] 
                            + (grid->tmp10_soil_a[h] - grid->tmp_2m_a[h]);
            grid->tmp200_soil[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY + offset][h] 
                            + (grid->tmp200_soil_a[h] - grid->tmp_2m_a[h]);
            grid->tcdc_clm[h] = grid->hist_cld[grid->climy - PIVOT_CLIMY + offset][h];
            grid->prate_sfc[h] = grid->hist_pre[grid->climy - PIVOT_CLIMY + offset][h];  
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
            grid->tmp_sfc[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY + offset][h] 
                            + (grid->tmp_sfc_a[h] - grid->tmp_2m_a[h]);
            grid->tmp_2m[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY + offset][h];
            grid->tmp10_soil[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY + offset][h] 
                            + (grid->tmp10_soil_a[h] - grid->tmp_2m_a[h]);
            grid->tmp200_soil[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY + offset][h] 
                            + (grid->tmp200_soil_a[h] - grid->tmp_2m_a[h]);
            grid->tcdc_clm[h] = grid->hist_cld[grid->climy - PIVOT_CLIMY + offset][h];
            grid->prate_sfc[h] = grid->hist_pre[grid->climy - PIVOT_CLIMY + offset][h];  
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
                tmp_var = grid->proj_tmp2m[grid->climy - PIVOT_GCMY][h][grid->gcm_row][grid->gcm_col] -
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
            grid->tmp_2m[h] = grid->proj_tmp2m[grid->climy - PIVOT_GCMY][h][grid->gcm_row][grid->gcm_col];
            grid->tmp_sfc[h] = grid->proj_tmp2m[grid->climy - PIVOT_GCMY][h][grid->gcm_row][grid->gcm_col]
                                + (grid->tmp_sfc_a[h] - grid->tmp_2m_a[h]);
            grid->tmp10_soil[h] = grid->proj_tmp2m[grid->climy - PIVOT_GCMY][h][grid->gcm_row][grid->gcm_col]
                                + (grid->tmp10_soil_a[h] - grid->tmp_2m_a[h]);
            grid->tmp200_soil[h] = grid->proj_tmp2m[grid->climy - PIVOT_GCMY][h][grid->gcm_row][grid->gcm_col]
                                + (grid->tmp200_soil_a[h] - grid->tmp_2m_a[h]);
       }
        
        /* precipitation ***************************************************/
        if(TEMP_GC != 0){
            pre_var = 0.0;
        }else{
            if(CC_P == 1){
                pre_var = grid->proj_prec[grid->climy-PIVOT_GCMY][h][grid->gcm_row][grid->gcm_col] - 
                            grid->proj_prec_b[h][grid->gcm_row][grid->gcm_col];
            }else{
                pre_var = 0.0;
            }
        }
        grid->prate_sfc[h] = grid->prate_sfc_a[h] + pre_var;
       
        if(ISIMIP_RUN == 2){
            /* PLUME: 2014/07/31 by A.Ito */
            grid->prate_sfc[h] = grid->proj_prec[grid->climy-PIVOT_GCMY][h][grid->gcm_row][grid->gcm_col];
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
            apres = 1013.25*exp(-1.0 * (28.964*0.001) * GAC * alt/(UGC * (grid->tmp_2m[h]+  ZAT)));
            
            if(CC_H == 1){
                rvap = grid->proj_hum[grid->climy - PIVOT_GCMY][h][grid->gcm_row][grid->gcm_col] / apres;
                grid->spfh_2m[h] = 0.622 * rvap / (1.0 - 0.378 * rvap);
                
                /* shm_var = grid->proj_hum[grid->climy-PIVOT_GCMY-1][h][grid->gcm_row][grid->gcm_col] - 
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
                rad_var = grid->proj_rad[grid->climy-PIVOT_GCMY][h][grid->gcm_row][grid->gcm_col] - 
                            grid->proj_rad_b[h][grid->gcm_row][grid->gcm_col];
            }else if(CC_R == 0 || CC_R == 3){
                rad_var = 0.0;	/* mean SW & PAR */
            }else{
                rad_var = 0.0;	
            }
        }

        /* experiment for SRM by reflector */
        /* added: 2014/07/06 by A.Ito     */
        if(GCM_ID == 3313 || GCM_ID == 3913){
            grid->top_rad[h] = f_top_rad(grid, 0);
        }
        
        grid->gl_rad[h] = grid->rad_a[h] + rad_var;
        
        if(ISIMIP_RUN == 2){
            /* PLUME: 2014/07/31 by A.Ito */
            grid->gl_rad[h] = grid->proj_rad[grid->climy-PIVOT_GCMY][h][grid->gcm_row][grid->gcm_col];
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
            rad_var = grid->proj_rad[grid->climy-PIVOT_GCMY-1][h][grid->gcm_row][grid->gcm_col] - 
                    grid->proj_rad_b[h][grid->gcm_row][grid->gcm_col];
                    
            grid->gl_rad[h] = grid->rad_a[h] + rad_var;
            if(grid->gl_rad[h]<0.0){
                grid->gl_rad[h] = 0.0;	/* SW change */
            }
        }
    }
}

