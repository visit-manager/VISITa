/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* History */
/* Revised August 15, 2007				*/
/* Revised August 19, 2007				*/
/* Revised 2007 / 12 / 25 by A.Ito				*/
/* Revised 2008 / 02 / 15 by A.Ito				*/
/* Revised 2008 / 03 / 10 by A.Ito				*/

/* VISIT a: global model	*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"structure.h"
#include"prototype.h"

#define FROWS 10

/* structure ************************************************/
struct Grid grid;
struct Loct loct;
struct Echar echar;
struct Mass mass;
struct Flux flux;

struct Loct loct_agr;
struct Echar echar_agr;
struct Mass mass_agr;
struct Flux flux_agr;

struct Loct loct_nat;
struct Echar echar_nat;
struct Mass mass_nat;
struct Flux flux_nat;

/* global variables ********************************************/
/* GCM variables ************/
long GCM, CO2S, GCM_R, GCM_C; /* */

/* atm. GHG scenario */
double aco2_a1[553], aco2_a2[553], aco2_b1[553], aco2_b2[553];
double ach4_a1[553], ach4_a2[553], ach4_b1[553], ach4_b2[553];
double an2o_a1[553], an2o_a2[553], an2o_b1[553], an2o_b2[553];

double glandarea;
/* historical results *************/
double h_tmp[201], h_pre[201], h_dswr[201], h_aet[201], h_rof[201];
double h_gpp[201], h_npp[201], h_nep[201], h_plant[201], h_soil[201];
double h_sr[201], h_ersn_c[201], h_agrersn_c[201], h_doc[201];
double h_agrarea[201], h_luc[201];

double h_burnt_area[201];
double h_bioburn_co2[201], h_bioburn_ch4[201], h_bioburn_co[201];
double h_bioburn_nmhc[201], h_bioburn_oc[201], h_bioburn_bc[201];
double h_bioburn_nox[201], h_bioburn_so2[201], h_bioburn_pm25[201];
double h_bioburn_tpm[201], h_bioburn_tec[201];

double h_ch4ox1[201], h_ch4ox2[201], h_ch4ox3[201];
double h_ch4emit_cao_paddy[201], h_ch4emit_cao_wetland[201];
double h_n2o_emit_ngas[201], h_n2_emit_ngas[201];
double h_n2o_emit_casa[201], h_no_emit_casa[201], h_n2_emit_casa[201];
double h_nh3_emit[201], h_n2_biofix[201];
double h_ch4_emit_mass[201], h_ch4_emit_photo[201];

double h_voc_isopr_g97[201], h_voc_monotrp_g97[201], h_voc_methanl_g97[201];
double h_voc_acetone_g97[201], h_voc_actaldhd_g97[201], h_voc_frmardhd_g97[201];
double h_voc_formacd_g97[201], h_voc_acetacd_g97[201], h_voc_co_g97[201];

/* monthly results **********/
double m_ch4ox1[12], m_ch4ox2[12], m_ch4ox3[12];
double m_bioburn_co2[12], m_bioburn_ch4[12], m_bioburn_co[12];
double m_bioburn_nmhc[12], m_bioburn_oc[12], m_bioburn_bc[12];

/* vegetation (olson) results */
double v_area[34];
double v_gpp[34], v_npp[34], v_nep[34];
double v_lai[34], v_fol[34], v_stm[34], v_rot[34], v_ltr[34], v_msl[34];

/*** main simulation roop ******************************************/
int main(
	int argc, 
	char *argv[]
){
	short zone;
	long f, g, h;
	/* global land area */
	/* file name strings */
	char filename[100];
	char s_date[25];
	char s_case[25];
	/* file pointer */
	FILE *fp_s[IFILEN];
	FILE *fp_c[4];
	FILE *fp_o[7];
	
	/* setting by arguments ***************************/
	if(argc<3){
		printf("No scenario specified !!!\n");
		exit(1);
	}
	
	/* climate scenario used: see setting.h */
	GCM = atol(argv[1]);
	/* GCM = 1003; */
	if(GCM>=0 && GCM<=2000){
		;
	}else{
		printf("Bad scenario specified !!!\n");
		exit(1);
	}
	
	/***************************************************/
	set_rowcol_gcm();	/* -> vegetdeal.c */
	
	printf("Open files...");
	/** open source files ***************************************/
	open_input(fp_s, fp_c);	/* -> open_input.c */
	
	/** open result file ******************************************/
	set_gcm_index(s_case);	/* -> vegetdeal.c */
		
	/***********/
	strcpy(s_date, argv[2]);
	strcat(s_date, "_");
				
	printf("done\n");
	
	/*******************************************************************/
	printf("Initialize simulation...");
	
	/** initialize configuration **/
	initSim(&grid);
	printf("done\n");

	echar.v_type = 1; /* olson */
	echar_nat.v_type = 2; /* sage */
	
	/*** read GCM climate scenario ***/
	printf("Reading GCM climate projection...");
	read_gcm_clim(&grid); /* */
	printf("done\n");
	
	/******************************************************************/
	/* latitude roop **************/
	printf("Start simulation...\n");
	for(f=0;f<360;f++){
		/* zone number: e.g., 5-deg latitude */
		zone = (short)f/FROWS + 1;
		
		if(f%FROWS == 0){
			/* open output files */	/* -> vegetdeal.c */
			f_output_file_open(zone, s_date, s_case, filename, fp_o);
		}
		
		/* grid area */
		grid.area = grid_area(90.0-0.5*(double)f, 89.5-0.5*(double)f, 0.0, 0.5);

		/* longitude roop **************/
		for(g=0;g<720;g++){
			/** roop of longitudinal rows, from W to E at the 0.5deg resolution **/	
			grid.row = f; 
			grid.col = g;	
			
			/* initialize grid condition */
			initG(fp_s, &grid);
			
			/* read CRU climate data */
			read_cru_clim(fp_c, &grid);
			
			printf("%3ld %3ld: %7.2lf %7.2lf: %2ld: %1ld :area %.1lf\n", 
				grid.row, grid.col, grid.lat, grid.lon, grid.veg_olson, grid.cru_exist, grid.area); /* */
			
			fprintf(fp_o[0],"%ld %ld %ld\n", grid.row, grid.col, grid.veg_olson); 
			fprintf(fp_o[1],"%ld %ld %ld\n", grid.row, grid.col, grid.veg_olson); 
			fprintf(fp_o[2],"%ld %ld %ld %ld %ld %lf\n", 
				grid.row, grid.col, grid.veg_olson, grid.veg_sage, grid.rvbasin, grid.f_crop_trend); 
			fprintf(fp_o[3],"%ld %ld %ld\n", grid.row, grid.col, grid.veg_olson); 
			fprintf(fp_o[4],"%ld %ld %ld\n", grid.row, grid.col, grid.veg_olson); 
			fprintf(fp_o[5],"%ld %ld %ld\n", grid.row, grid.col, grid.veg_olson); 
			fprintf(fp_o[6],"%ld %ld %ld\n", grid.row, grid.col, grid.veg_olson); 
			
			if(grid.veg_olson!=0 && grid.veg_olson!=33 && (g+0)%1==0){ /* */ /* && (g+0)%10==0 */
			/* if( (grid.veg_olson==10 || grid.veg_olson==12) && (g+0)%1==0){  */
				grid.nnn++;
				glandarea += grid.area;
				v_area[grid.veg_olson] += grid.area;
				
				/*************************************/
				/** clear all parameters **/
				clear(&grid, &loct, &echar, &mass, &flux);
				clear(&grid, &loct_nat, &echar_nat, &mass_nat, &flux_nat); 
				
				/** initialize vegatation and soil conditions **/ 
				/* input parameter */
				initVS(&grid, &loct, &mass, &flux, &echar);				
				initVS(&grid, &loct_nat, &mass_nat, &flux_nat, &echar_nat);

				/** initialize climate conditions **/
				initC(&grid);
				
				/** initialize location conditions **/
				initL(&grid, &loct, &mass, &flux, &echar);
				initL(&grid, &loct_nat, &mass_nat, &flux_nat, &echar_nat);
				
				/** initialize stable carbon isotope **/
				init_d13c(&grid, &flux, &echar, &mass);
				init_d13c(&grid, &flux_nat, &echar_nat, &mass_nat);
				
				/** stabilization roop **/
				cal_stable(&grid, &loct, &echar, &mass, &flux, fp_o);
				
				/** outputs **/
				screenshow(&grid, &loct, &mass, &flux, &echar);

				/* printf("Lan use change: %.3lf<-%.3lf Emit: %.3lf\n", grid.f_crop_p, 
					grid.f_crop, (flux.lu_conv+flux.lu_ten+flux.lu_hund)); */
			
				/** dynamic roop **/
				/* 1901-2000 */
				cal_cruclim(&grid, &loct, &echar, &mass, &flux, fp_o);  /* */
				/* printf("Lan use change: %.3lf<-%.3lf Emit: %.3lf\n", grid.f_crop_p, 
					grid.f_crop, (flux.lu_conv+flux.lu_ten+flux.lu_hund)); */
			
				/* 2001-2100 */
				/* cal_gcmclim2(&grid, &loct, &echar, &mass, &flux, fp_o); */
				
				/* printf("Lan use change: %.3lf<-%.3lf Emit: %.3lf\n", grid.f_crop_p, 
					grid.f_crop, (flux.lu_conv+flux.lu_ten+flux.lu_hund)); */
				
				/* printf("CUM LUC: %.1lf %.1lf %.1lf\n\n", h_luc[0], h_luc[100], h_luc[200]); */
				
				printf("\n");
			}
		}
		
		/* output of global total values */
		f_glosum_output(s_date, s_case);
		
		if(f%FROWS == (FROWS-1)){
			/* close output files */
			fclose(fp_o[0]); 
			fclose(fp_o[1]); 
			fclose(fp_o[2]);
			fclose(fp_o[3]);
			fclose(fp_o[4]);
			fclose(fp_o[5]);
			fclose(fp_o[6]);
		}
	}
	
	/* close files */
	for(h=0;h<IFILEN;h++){
		fclose(fp_s[h]); 
	}
	for(h=0;h<4;h++){
		fclose(fp_c[h]); 
	}
	
	printf("Simulation ended\n");
	
	return 0;
}

