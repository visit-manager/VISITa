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
/* Revised 2008 / 07 / 01 by A.Ito				*/
/* Revised 2008 / 09 / 24 by A.Ito	(based on E.Kato's comments)	*/
/* Revised 2008 / 12 / 05 by A.Ito	radiation sensitivity analysis	*/

/* a previous version, Sim-CYCLE was described in
Ito, A. and Oikawa, T., 2002. A simulation model of the carbon cycle in land 
ecosystems (Sim-CYCLE): A description based on dry-matter production theory 
and plot-scale validation. Ecological Modelling, 151:147-179.
*/

/* VISIT a: global model	*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"structure.h"
#include"prototype.h"

#define FROWS 10

/* structure *****************************************************/
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

/* historical results: 201 years, spinup+1901-2100 *************/
double h_tmp[201], h_pre[201], h_dswr[201], h_aet[201], h_rof[201];
double h_gpp[201], h_npp[201], h_nep[201], h_plant[201], h_soil[201];
double h_sr[201], h_ersn_c[201], h_agrersn_c[201], h_doc[201];
double h_agrarea[201], h_luc[201];
double h_luc_1[201], h_luc_2[201], h_luc_3[201];
double h_gpp_df97[201];

double h_burnt_area[201];
double h_bioburn_co2[201], h_bioburn_ch4[201], h_bioburn_co[201];
double h_bioburn_nmhc[201], h_bioburn_oc[201], h_bioburn_bc[201];
double h_bioburn_nox[201], h_bioburn_so2[201], h_bioburn_pm25[201];
double h_bioburn_tpm[201], h_bioburn_tec[201];

double h_ch4ox1[201], h_ch4ox2[201], h_ch4ox3[201], h_ch4ox4[201];
double h_ch4emit_cao_paddy[201], h_ch4emit_cao_wetland[201];
double h_n2o_emit_ngas[201], h_n2_emit_ngas[201];
double h_n2o_emit_casa[201], h_no_emit_casa[201], h_n2_emit_casa[201];
double h_nh3_emit[201], h_n2_biofix[201];
double h_ch4_emit_mass[201], h_ch4_emit_photo[201];

double h_voc_isopr_g97[201], h_voc_monotrp_g97[201], h_voc_methanl_g97[201];
double h_voc_acetone_g97[201], h_voc_actaldhd_g97[201], h_voc_frmardhd_g97[201];
double h_voc_formacd_g97[201], h_voc_acetacd_g97[201], h_voc_co_g97[201];

/* monthly mean results **********/
double m_ch4ox1[ASTEP], m_ch4ox2[ASTEP], m_ch4ox3[ASTEP];
double m_bioburn_co2[ASTEP], m_bioburn_ch4[ASTEP], m_bioburn_co[ASTEP];
double m_bioburn_nmhc[ASTEP], m_bioburn_oc[ASTEP], m_bioburn_bc[ASTEP];
double m_gpp[ASTEP], m_npp[ASTEP], m_nep[ASTEP];

/* vegetation (olson) mean results */
double go_landarea, gs_landarea;
double vo_area[34];
double vo_gpp[34], vo_npp[34], vo_nep[34];
double vo_lai[34], vo_fol[34], vo_stm[34], vo_rot[34], vo_ltr[34], vo_msl[34];
double vs_area[16];
double vs_gpp[16], vs_npp[16], vs_nep[16];
double vs_lai[16], vs_fol[16], vs_stm[16], vs_rot[16], vs_ltr[16], vs_msl[16];

double MDN[ASTEP] = {31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0};

/* 0: 1950s */
/* 1: 1990s */
/* 2: 2020s */
/* 3: 2050s */
/* 4: 2080s */
float g_tmp[5][360][720];
float g_prc[5][360][720];
float g_swr[5][360][720];
float g_gpp[5][360][720];
float g_npp[5][360][720];
float g_nep[5][360][720];
float g_pmas[5][360][720];
float g_smas[5][360][720];
float g_ch4e[5][360][720];
float g_ch4o[5][360][720];
float g_n2oe[5][360][720];
float g_bbco2[5][360][720];
float g_ersn[5][360][720];
float g_isopr[5][360][720];
float g_sr[5][360][720];
float g_luc[5][360][720];  /* */

short RAD_SENS;
/* 0: control */
/* 1: +10 diffuse PAR fraction */
/* 2: -10 diffuse PAR fraction */
/* 3: +10 soil albedo */
/* 4: -10 soil albedo */
/* 5: +10 canopy scattering coefficient  */
/* 6: -10 canopy scattering coefficient */
/* 7: +10 diffuse attenuation coefficient */
/* 8: -10 diffuse attenuation coefficient */
/* 9: +10 direct attenuation coefficient */
/* 10: -10 direct attenuation coefficient */

short DF97;
/* 0: off - Monsi-Saeki */
/* 1: on De Pury-Farquhar */

/* future solar radiation change */
short CC_R;
/* 0: off */
/* 1: on  */

/* main simulation roop *************************************************/
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
	FILE *fp_o1[OFILES], *fp_o2[OFILES], *fp_o3[OFILES];
	FILE *fp_binout;
	
	/* setting by arguments ******************/
	if(argc<6){
		printf("No scenario specified !!!\n");
		exit(1);
	}
	/* argv[]: ARGUMEMTS ***********
	  0: execution filename itself (e.g., visita)
	  1: experiomental scenario ID number (see setting.h)
	  2: file identifier (arbitrary phrase such as date, your name, etc.)
	  3: code for sensitivity analysis
	  4: photosynthesis model 0(Monsi-Saeki) or 1(DePury-Farquhar)
	  5: future solar radiation change
	*/
	
	/* climate scenario used: see setting.h */
	GCM = atol(argv[1]);
	if(GCM>=0 && GCM<=2000){
		;
	}else{
		printf("Bad scenario ID specified !!!\n");
		exit(1);
	}
	
	RAD_SENS = (short)atol(argv[3]);
	
	DF97 = (short)atol(argv[4]);
	
	CC_R = (short)atol(argv[5]);
	
	/************************************************************/
	set_rowcol_gcm();	/* -> vegetdeal.c */
	
	printf("Open files...");
	/* open source files *****************************************/
	open_input(fp_s, fp_c);	/* -> open_input.c */
	
	/* open result file ******************************************/
	set_gcm_index(s_case);	/* -> vegetdeal.c */
		
	/***********/
	strcpy(s_date, argv[2]);
	strcat(s_date, "_");
				
	printf("done\n");
	
	strcpy(filename, argv[2]);
	strcat(filename, "_");
	strcat(filename, s_case);
	strcat(filename,"-output.bin");
	fp_binout = fopen(filename, "wb");
	
	for(f=0;f<5;f++){
		for(g=0;g<360;g++){
			for(h=0;h<720;h++){
				g_tmp[f][g][h] = 0.0;
				g_prc[f][g][h] = 0.0;
				g_swr[f][g][h] = 0.0;
				g_gpp[f][g][h] = 0.0;
				g_npp[f][g][h] = 0.0;
				g_nep[f][g][h] = 0.0;
				g_pmas[f][g][h] = 0.0;
				g_smas[f][g][h] = 0.0;
				g_ch4e[f][g][h] = 0.0;
				g_ch4o[f][g][h] = 0.0;
				g_n2oe[f][g][h] = 0.0;
				g_bbco2[f][g][h] = 0.0;
				g_ersn[f][g][h] = 0.0;
				g_isopr[f][g][h] = 0.0;
				g_sr[f][g][h] = 0.0;
				g_luc[f][g][h] = 0.0;
			}
		}
	}
	
	/*******************************************************************/
	printf("Initialize simulation...");
	
	/** initialize configuration **/
	initSim(&grid);
	printf("done\n");
	
	/* vegetation type: olson */
	(echar.c3).v_type = 1; 
	(echar.c4).v_type = 1;
	(echar.soil).v_type = 1;
	/* vegetation type: SAGE - U.Wisconsin */
	(echar_nat.c3).v_type = 2; 
	(echar_nat.c4).v_type = 2;
	(echar_nat.soil).v_type = 2;
	/* vegetation type: cropland */
	(echar_agr.c3).v_type = 3; 
	(echar_agr.c4).v_type = 3;
	(echar_agr.soil).v_type = 3;
		
	/*** read GCM climate scenario ***/
	printf("Reading GCM climate projection...");
	read_gcm_clim(&grid); /* */
	/* read_ncep_clim(&grid); */
	printf("done\n");
	
	/************************************************************************/
	/* latitude loop: north to south **************/
	printf("Start simulation...\n");
	for(f=0;f<360;f++){
		/* zone number: e.g., 5-deg latitude */
		zone = (short)f/FROWS + 1;
		
		if(f%FROWS == 0){
			/* open output files */	/* -> output_file.c */
			if(CALC_OLSON == 1){
				f_output_file_open(1, zone, s_date, s_case, filename, fp_o1);
			}
			if(CALC_SAGE == 1){
				f_output_file_open(2, zone, s_date, s_case, filename, fp_o2);
			}
			if(CALC_CROP == 1){
				f_output_file_open(3, zone, s_date, s_case, filename, fp_o3);
			}
		}
		
		/* grid area */	/* -> vegetdeal.c */
		grid.area = grid_area(90.0-0.5*(double)f, 89.5-0.5*(double)f, 0.0, 0.5);

		/* longitude loop: west to east ****************************************/
		for(g=0;g<720;g++){
			/* grid position */
			grid.row = f; 
			grid.col = g;	
			
			/* initialize grid condition: data setting */
			init_grid(fp_s, &grid);
			
			/* read CRU TS2.1 climate data */
			read_cru_clim(fp_c, &grid);
			
			printf("%3ld %3ld: %7.2lf %7.2lf: %2ld %2ld %2ld: %1ld\n", 
				grid.row, grid.col, grid.lat, grid.lon, grid.veg_olson, grid.veg_sage, 
				grid.veg_crop, grid.cru_exist); /* */
			
			/* head record */
			for(h=0;h<OFILES;h++){
				if(CALC_OLSON == 1){
					fprintf(fp_o1[h],"%ld %ld %ld %ld\n", grid.row, grid.col, grid.veg_olson, grid.veg_sage); 
				}
				if(CALC_SAGE == 1){
					fprintf(fp_o2[h],"%ld %ld %ld %ld\n", grid.row, grid.col, grid.veg_olson, grid.veg_sage); 
				}
				if(CALC_CROP == 1){
					fprintf(fp_o3[h],"%ld %ld %ld %ld\n", grid.row, grid.col, grid.veg_olson, grid.veg_sage); 
				}
			}
			
			/* calculation for lands *******************************************/
			/* Olson map */
			if(CALC_OLSON == 1){
				if(grid.veg_olson!=0 && grid.veg_olson!=33 && (g+0)%1==0){
					/* sequential number */
					grid.n_olson++;
					/* total area */
					go_landarea += grid.area;
					vo_area[grid.veg_olson] += grid.area;

					/* clear all parameters ******/
					clear(&grid, &loct, &echar, &mass, &flux);

					/* initialize vegatation and soil conditions ****/ 
					initVS(&grid, &loct, &mass, &flux, &echar);				
				
					/* initialize climate conditions ****/
					initC(&grid);
				
					/* initialize location conditions ****/
					initL(&grid, &loct, &mass, &flux, &echar);
				
					/* initialize stable carbon isotope ****/
					init_d13c(&grid, &flux, &echar, &mass);

					/* spin-up: stabilization roop ***************************/
					cal_stable(&grid, &loct, &echar, &mass, &flux, fp_o1);

					/* snap shot for checking *****/
					screenshow(&grid, &loct, &mass, &flux, &echar); 
					
					/* experiment *******************************************/
					/* past: 1901-2000 */
					cal_cruclim(&grid, &loct, &echar, &mass, &flux, fp_o1);	

					/* future: 2001-2100 */
					if(CRU_PD==102 || CRU_PD==108){
						;
					}else{
						cal_gcmclim2(&grid, &loct, &echar, &mass, &flux, fp_o1);
					}
					
					printf("\n");
				}
			}
			
			/* SAGE map */
			if(CALC_SAGE == 1){
				if(grid.veg_sage!=0 && grid.veg_sage!=15 && (g+0)%1==0){
					/* sequential number */
					grid.n_sage++;
					/* total area */
					gs_landarea += grid.area;
					vs_area[grid.veg_sage] += grid.area;

					/* clear all parameters ******/
					clear(&grid, &loct_nat, &echar_nat, &mass_nat, &flux_nat); 

					/* initialize vegatation and soil conditions ****/ 
					initVS(&grid, &loct_nat, &mass_nat, &flux_nat, &echar_nat);
				
					/* initialize climate conditions ****/
					initC(&grid);
				
					/* initialize location conditions ****/
					initL(&grid, &loct_nat, &mass_nat, &flux_nat, &echar_nat);
				
					/* initialize stable carbon isotope ****/
					init_d13c(&grid, &flux_nat, &echar_nat, &mass_nat);

					/* spin-up: stabilization roop ***************************/
					cal_stable(&grid, &loct_nat, &echar_nat, &mass_nat, &flux_nat, fp_o2);

					/* snap shot for checking *****/
					screenshow(&grid, &loct_nat, &mass_nat, &flux_nat, &echar_nat);
					
					/* experiment **************************************/
					/* past: 1901-2000 */
					cal_cruclim(&grid, &loct_nat, &echar_nat, &mass_nat, &flux_nat, fp_o2);

					/* future: 2001-2100 */
					cal_gcmclim2(&grid, &loct_nat, &echar_nat, &mass_nat, &flux_nat, fp_o2);
					
					printf("\n");
				}
			}
			
			/* croplands */
			if(CALC_CROP == 1){
				if((g+0)%1==0 && (grid.veg_sage!=0 && grid.veg_sage!=15)){
					/* sequential number */
					grid.n_crop++;
					/* total area */

					/* clear all parameters ******/
					clear(&grid, &loct_agr, &echar_agr, &mass_agr, &flux_agr); 

					/* initialize vegatation and soil conditions ****/ 
					initVS(&grid, &loct_agr, &mass_agr, &flux_agr, &echar_agr);
				
					/* initialize climate conditions ****/
					initC(&grid);
				
					/* initialize location conditions ****/
					initL(&grid, &loct_agr, &mass_agr, &flux_agr, &echar_agr);
				
					/* initialize stable carbon isotope ****/
					init_d13c(&grid, &flux_agr, &echar_agr, &mass_agr);

					/* spin-up: stabilization roop ***************************/
					cal_stable(&grid, &loct_agr, &echar_agr, &mass_agr, &flux_agr, fp_o3);

					/* snap shot for checking *****/
					screenshow(&grid, &loct_agr, &mass_agr, &flux_agr, &echar_agr);
					
					/* experiment **************************************/
					/* past: 1901-2000 */
				/*	cal_cruclim(&grid, &loct_agr, &echar_agr, &mass_agr, &flux_agr, fp_o3);  */

					/* future: 2001-2100 */
				/*	cal_gcmclim2(&grid, &loct_agr, &echar_agr, &mass_agr, &flux_agr, fp_o3);  */
					
					printf("\n");
				}
			}
		}
		/* end of longitudinal loop *************************/
		
		/* output of global total values ******/
		f_glosum_output(s_date, s_case);
		
		/* close output files */
		if(f%FROWS == (FROWS-1)){
			for(h=0;h<OFILES;h++){
				if(CALC_OLSON == 1){
					fclose(fp_o1[h]);
				}
				if(CALC_SAGE == 1){
					fclose(fp_o2[h]);
				}
				if(CALC_CROP == 1){
					fclose(fp_o3[h]);
				}
			}
		}
	}
	/* end of latitudinal loop *************************************/
	
	fwrite(g_tmp, sizeof(float), 5*360*720, fp_binout);  // 0-4
	fwrite(g_prc, sizeof(float), 5*360*720, fp_binout);  // 5-9
	fwrite(g_swr, sizeof(float), 5*360*720, fp_binout);  // 10-14
	fwrite(g_gpp, sizeof(float), 5*360*720, fp_binout);  // 15-19
	fwrite(g_npp, sizeof(float), 5*360*720, fp_binout);  // 20-24
	fwrite(g_nep, sizeof(float), 5*360*720, fp_binout);  // 25-29
	fwrite(g_pmas, sizeof(float), 5*360*720, fp_binout);  // 30-34
	fwrite(g_smas, sizeof(float), 5*360*720, fp_binout);  // 35-39
	fwrite(g_ch4e, sizeof(float), 5*360*720, fp_binout);  // 40-44
	fwrite(g_ch4o, sizeof(float), 5*360*720, fp_binout);  // 45-49
	fwrite(g_n2oe, sizeof(float), 5*360*720, fp_binout);  // 50-54
	fwrite(g_bbco2, sizeof(float), 5*360*720, fp_binout);  // 55-59
	fwrite(g_ersn, sizeof(float), 5*360*720, fp_binout);  // 60-64
	fwrite(g_isopr, sizeof(float), 5*360*720, fp_binout);  // 65-69
	fwrite(g_sr, sizeof(float), 5*360*720, fp_binout);  // 70-74
	fwrite(g_luc, sizeof(float), 5*360*720, fp_binout);  // 75-79
		
	/* close files */
	for(h=0;h<IFILEN;h++){
		fclose(fp_s[h]); 
	}
	for(h=0;h<4;h++){
		fclose(fp_c[h]); 
	}
	
	fclose(fp_binout);
	
	printf("Simulation ended\n");
	
	return 0;
}

