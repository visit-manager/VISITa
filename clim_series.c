/*	VISIT: Vegetation Integrative SImulation Tool						*/
/*  Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/*  Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/*  Carbon cycle, erosion, biomass burning, land-use change,			*/
/*  CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

extern short CC_R;

/*******************************/
void set_cru_clim(
	struct Grid *grid
){
	long h;
	double tmp_var, pre_var, tcdc_var;
	
	if(grid->climy<=2002){
		for(h=0;h<ASTEP;h++){
			grid->tmp_sfc[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY][h] + (grid->tmp_sfc_a[h] - grid->tmp_2m_a[h]);
			grid->tmp_2m[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY][h];
			grid->tmp10_soil[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY][h] + (grid->tmp10_soil_a[h] - grid->tmp_2m_a[h]);
			grid->tmp200_soil[h] = grid->hist_tmp[grid->climy - PIVOT_CLIMY][h] + (grid->tmp200_soil_a[h] - grid->tmp_2m_a[h]);
			grid->tcdc_clm[h] = grid->hist_cld[grid->climy - PIVOT_CLIMY][h];
			grid->prate_sfc[h] = grid->hist_pre[grid->climy - PIVOT_CLIMY][h];  
		}
	}else{
		/* 2003-2008: extrapolation using NCEP/NCAR data: 2009/01/05 by A.Ito */
		for(h=0;h<ASTEP;h++){
			/* temperature */
			tmp_var = grid->ncep_tmp2m[grid->climy - 1948][h][grid->ncep_lat][grid->ncep_lon] - grid->ncep_tmp2m_b[h][grid->ncep_lat][grid->ncep_lon];
			grid->tmp_sfc[h] = grid->tmp_sfc_a[h] + tmp_var;
			grid->tmp_2m[h] = grid->tmp_2m_a[h] + tmp_var;
			grid->tmp10_soil[h]=grid->tmp10_soil_a[h] + tmp_var;
			grid->tmp200_soil[h]=grid->tmp200_soil_a[h] + tmp_var;
			
			/* precipitation */
			pre_var = grid->ncep_prate[grid->climy - 1948][h][grid->ncep_lat][grid->ncep_lon] - grid->ncep_prate_b[h][grid->ncep_lat][grid->ncep_lon];
			grid->prate_sfc[h] = grid->prate_sfc_a[h] + pre_var;
			if(grid->prate_sfc[h]<0.0){
				grid->prate_sfc[h] = 0.0;
			}
			
			/* cloudiness */
			tcdc_var = grid->ncep_tcdc[grid->climy - 1948][h][grid->ncep_lat][grid->ncep_lon] - grid->ncep_tcdc_b[h][grid->ncep_lat][grid->ncep_lon];
			grid->tcdc_clm[h] = grid->tcdc_clm_a[h] + tcdc_var;
			if(grid->tcdc_clm[h]<0.0){
				grid->tcdc_clm[h]=0.0;
			}
			if(grid->tcdc_clm[h]>1.0){
				grid->tcdc_clm[h]=1.0;
			}
		}
	}
}

/*******************************/
void set_gcm_clim(
	struct Grid *grid
){
	long h;
	double tmp_var, pre_var, shm_var, rad_var;
	
	/***************************************************/
	for(h=0;h<ASTEP;h++){
		grid->m = h;
		/****** temperature ******/
		if(CC_T == 1){
			tmp_var = grid->proj_tmp2m[grid->climy-PIVOT_GCMY-1][h][grid->gcm_row][grid->gcm_col] - 
					grid->proj_tmp2m_b[h][grid->gcm_row][grid->gcm_col];
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
		
		/****** precipitation ******/
		if(CC_P==1){
			pre_var = grid->proj_prec[grid->climy-PIVOT_GCMY-1][h][grid->gcm_row][grid->gcm_col] - 
						grid->proj_prec_b[h][grid->gcm_row][grid->gcm_col];
		}else{
			pre_var = 0.0;
		}
		grid->prate_sfc[h] = grid->prate_sfc_a[h] + pre_var;
		if(grid->prate_sfc[h]<=0.0){
			/* carry over of negative precipitation */
			grid->proj_prec_co += grid->prate_sfc[h];
			grid->prate_sfc[h] = 0.0;
		}else if(grid->prate_sfc[h]>0.0){
			grid->prate_sfc[h] += grid->proj_prec_co;
			grid->proj_prec_co = 0.0;
			if(grid->prate_sfc[h] <= 0.0){
				/* carry over of negative precipitation */
				grid->proj_prec_co += grid->prate_sfc[h];
				grid->prate_sfc[h] = 0.0;
			}
		}						

		/****** air humidity ******/
		if(CC_H==1){
			shm_var = grid->proj_shum[grid->climy-PIVOT_GCMY-1][h][grid->gcm_row][grid->gcm_col] - 
						grid->proj_shum_b[h][grid->gcm_row][grid->gcm_col];
		}else{
			shm_var = 0.0;
		}
		grid->spfh_2m[h] = grid->spfh_2m_a[h] + shm_var;
		if(grid->spfh_2m[h]<0.0){
			grid->spfh_2m[h] = 0.0;
		}
		
		/****** surface shortwave radiation ******/
		/*  0: no change  */
		/*  1: SW and PAR change  */
		/*  2: PAR change but mean SW  */
		/*  3: SW change but mean PAR  */
		if(CC_R==1||CC_R==2){
			rad_var = grid->proj_rad[grid->climy-PIVOT_GCMY-1][h][grid->gcm_row][grid->gcm_col] - 
						grid->proj_rad_b[h][grid->gcm_row][grid->gcm_col];
		}else if(CC_R==0||CC_R==3){	
			rad_var = 0.0;	/* mean SW & PAR */
		}
		grid->gl_rad[h] = grid->rad_a[h] + rad_var;
		if(grid->gl_rad[h]<0.0){
			grid->gl_rad[h] = 0.0;
		}
		grid->par[h] = par(grid);
		
		if(CC_R==2){
			grid->gl_rad[h] = grid->rad_a[h];	/* mean SW */
		}
		if(CC_R==3){
			rad_var = grid->proj_rad[grid->climy-PIVOT_GCMY-1][h][grid->gcm_row][grid->gcm_col] - 
					grid->proj_rad_b[h][grid->gcm_row][grid->gcm_col];
					
			grid->gl_rad[h] = grid->rad_a[h] + rad_var;
			if(grid->gl_rad[h]<0.0){
				grid->gl_rad[h] = 0.0;	/* SW change */
			}
		}
	}
}

