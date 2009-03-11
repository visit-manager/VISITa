/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
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

extern long GCM, CO2S, GCM_R, GCM_C;
extern double aco2_a1[553], aco2_a2[553], aco2_b1[553], aco2_b2[553];
extern double ach4_a1[553], ach4_a2[553], ach4_b1[553], ach4_b2[553];
extern double an2o_a1[553], an2o_a2[553], an2o_b1[553], an2o_b2[553];

extern double glandarea;
extern double h_tmp[201], h_pre[201], h_dswr[201], h_aet[201], h_rof[201];
extern double h_gpp[201], h_npp[201], h_nep[201], h_plant[201], h_soil[201];
extern double h_sr[201], h_ersn_c[201], h_agrersn_c[201], h_doc[201];
extern double h_agrarea[201], h_luc[201];
extern double h_luc_1[201], h_luc_2[201], h_luc_3[201];
extern double h_gpp_df97[201];

extern double h_burnt_area[201];
extern double h_bioburn_co2[201], h_bioburn_ch4[201], h_bioburn_co[201];
extern double h_bioburn_nmhc[201], h_bioburn_oc[201], h_bioburn_bc[201];
extern double h_bioburn_nox[201], h_bioburn_so2[201], h_bioburn_pm25[201];
extern double h_bioburn_tpm[201], h_bioburn_tec[201];

extern double h_ch4ox1[201], h_ch4ox2[201], h_ch4ox3[201], h_ch4ox4[201];
extern double h_ch4emit_cao_paddy[201], h_ch4emit_cao_wetland[201];
extern double h_n2o_emit_ngas[201], h_n2_emit_ngas[201];
extern double h_n2o_emit_casa[201], h_no_emit_casa[201], h_n2_emit_casa[201];
extern double h_nh3_emit[201], h_n2_biofix[201];
extern double h_ch4_emit_mass[201], h_ch4_emit_photo[201];

extern double h_voc_isopr_g97[201], h_voc_monotrp_g97[201], h_voc_methanl_g97[201];
extern double h_voc_acetone_g97[201], h_voc_actaldhd_g97[201], h_voc_frmardhd_g97[201];
extern double h_voc_formacd_g97[201], h_voc_acetacd_g97[201], h_voc_co_g97[201];

/* monthly results **********/
extern double m_ch4ox1[12], m_ch4ox2[12], m_ch4ox3[12];
extern double m_bioburn_co2[12], m_bioburn_ch4[12], m_bioburn_co[12];
extern double m_bioburn_nmhc[12], m_bioburn_oc[12], m_bioburn_bc[12];
extern double m_gpp[ASTEP], m_npp[ASTEP], m_nep[ASTEP];

/* vegetation (olson) results */
extern double go_landarea, gs_landarea;
extern double vo_area[34];
extern double vo_gpp[34], vo_npp[34], vo_nep[34];
extern double vo_lai[34], vo_fol[34], vo_stm[34], vo_rot[34], vo_ltr[34], vo_msl[34];
extern double vs_area[16];
extern double vs_gpp[16], vs_npp[16], vs_nep[16];
extern double vs_lai[16], vs_fol[16], vs_stm[16], vs_rot[16], vs_ltr[16], vs_msl[16];

/* initialize general parameters used in the simulation ******************/
void initSim(
	struct Grid *grid
){
	long f, year;
	double data;
	FILE *fp_co2;
	
	/* grid ID, sequential number */
	grid->n_olson = 0;
	grid->n_sage = 0;
	grid->n_crop = 0;
		
	/* atm. CO2 scenario ****************************/
	printf("reading CO2 data...");
	if(CO2S==1){
		if((fp_co2 = fopen("./data/SRES_A1.dat","rt"))==NULL){
			printf("No SRES_A1.dat\n");
			exit(1);
		}
	}else if(CO2S==2){
		if((fp_co2 = fopen("./data/SRES_A1FI.dat","rt"))==NULL){
			printf("No SRES_A1FI.dat\n");
			exit(1);
		}
	}else if(CO2S==3){
		if((fp_co2 = fopen("./data/SRES_A1T.dat","rt"))==NULL){
			printf("No SRES_A1T.dat\n");
			exit(1);
		}
	}else if(CO2S==4){
		if((fp_co2 = fopen("./data/SRES_A2.dat","rt"))==NULL){
			printf("No SRES_A2.dat\n");
			exit(1);
		}
	}else if(CO2S==5){
		if((fp_co2 = fopen("./data/SRES_B1.dat","rt"))==NULL){
			printf("No SRES_B1.dat\n");
			exit(1);
		}
	}else if(CO2S==6){
		if((fp_co2 = fopen("./data/SRES_B2.dat","rt"))==NULL){
			printf("No SRES_B2.dat\n");
			exit(1);
		}
	}else{
		if((fp_co2 = fopen("./data/SRES_A1.dat","rt"))==NULL){
			printf("No SRES_A1.dat\n");
			exit(1);
		}
	}
	
	/* SRES Scenario CO2 */
	for(f=0;f<111;f++){
		fscanf(fp_co2,"%ld %lf", &year, &data);
		sres_co2[f] = data;
	}
	
	fclose(fp_co2);
	
	/* source: http://crga.atmos.uiuc.edu/research/post-sres.html
		M.E.Schlesinger and S.Malyshev			*/
	if((fp_co2 = fopen("./data/AtmGHG_timeseries.dat","rt"))==NULL){
		printf("No AtmGHG_timeseries.dat\n");
		exit(1);
	}
	for(f=0;f<553;f++){
		fscanf(fp_co2,"%ld", &year);
		/* CO2, ppmv */
		fscanf(fp_co2,"%lf", &aco2_a1[f]);
		fscanf(fp_co2,"%lf", &aco2_b1[f]);
		fscanf(fp_co2,"%lf", &aco2_a2[f]);
		fscanf(fp_co2,"%lf", &aco2_b2[f]);
		/* CH4, ppbv */
		fscanf(fp_co2,"%lf", &ach4_a1[f]);
		fscanf(fp_co2,"%lf", &ach4_b1[f]);
		fscanf(fp_co2,"%lf", &ach4_a2[f]);
		fscanf(fp_co2,"%lf", &ach4_b2[f]);
		/* N2O, ppbv */
		fscanf(fp_co2,"%lf", &an2o_a1[f]);
		fscanf(fp_co2,"%lf", &an2o_b1[f]);
		fscanf(fp_co2,"%lf", &an2o_a2[f]);
		fscanf(fp_co2,"%lf", &an2o_b2[f]);
	}
	fclose(fp_co2);
	
	/* global analysis initialization ********************************/
	go_landarea = gs_landarea = 0.0;
	for(f=0;f<201;f++){
		h_tmp[f] = h_pre[f] = h_dswr[f] = h_aet[f] = h_rof[f] = 0.0;
		h_gpp[f] = h_npp[f] = h_nep[f] = h_plant[f] = h_soil[f] = 0.0;
		h_sr[f] = h_ersn_c[f] = h_agrersn_c[f] = h_doc[f] = 0.0;
		h_gpp_df97[f] = 0.0;
		
		h_agrarea[f] = h_luc[f] = 0.0;
		h_luc_1[f] = h_luc_2[f] = h_luc_3[f] = 0.0;
		h_burnt_area[f] = 0.0;
		h_bioburn_co2[f] = h_bioburn_co[f] = h_bioburn_ch4[f] = 0.0;
		h_bioburn_nmhc[f] = h_bioburn_oc[f] = h_bioburn_bc[f] = 0.0;
		h_bioburn_nox[f] = h_bioburn_so2[f] = h_bioburn_pm25[f] = 0.0;
		h_bioburn_tpm[f] = h_bioburn_tec[f] = 0.0;
		
		h_ch4ox1[f] = h_ch4ox2[f] = h_ch4ox3[f] = h_ch4ox4[f] = 0.0;
		h_ch4emit_cao_paddy[f] = h_ch4emit_cao_wetland[f] = 0.0;
		h_n2o_emit_ngas[f] = h_n2_emit_ngas[f] = 0.0;
		h_n2o_emit_casa[f] = h_no_emit_casa[f] = h_n2_emit_casa[f] = 0.0;
		h_nh3_emit[f] = h_n2_biofix[f] = 0.0;
		h_ch4_emit_mass[f] = h_ch4_emit_photo[f] = 0.0;

		h_voc_isopr_g97[f] = h_voc_monotrp_g97[f] = h_voc_methanl_g97[f] = 0.0;
		h_voc_acetone_g97[f] = h_voc_actaldhd_g97[f] = h_voc_frmardhd_g97[f] = 0.0;
		h_voc_formacd_g97[f] = h_voc_acetacd_g97[f] = h_voc_co_g97[f] = 0.0;
	}
	for(f=0;f<ASTEP;f++){
		m_ch4ox1[f] = m_ch4ox2[f] = m_ch4ox3[f] = 0.0;
		m_bioburn_co2[f] = m_bioburn_ch4[f] = m_bioburn_co[f] = 0.0;
		m_bioburn_nmhc[f] = m_bioburn_oc[f] = m_bioburn_bc[f] = 0.0;
		m_gpp[f] = m_npp[f] = m_nep[f] = 0.0;
	}
	for(f=0;f<VEG_NUM_OLSON;f++){
		vo_area[f] = 0.0;
		vo_gpp[f] = vo_npp[f] = vo_nep[f] = 0.0;
		vo_lai[f] = vo_fol[f] = vo_stm[f] = vo_rot[f] = vo_ltr[f] = vo_msl[f] = 0.0;
	}
	for(f=0;f<VEG_NUM_SAGE;f++){
		vs_area[f] = 0.0;
		vs_gpp[f] = vs_npp[f] = vs_nep[f] = 0.0;
		vs_lai[f] = vs_fol[f] = vs_stm[f] = vs_rot[f] = vs_ltr[f] = vs_msl[f] = 0.0;
	}
}

