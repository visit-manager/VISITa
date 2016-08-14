/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
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
/* Revised 2008 / 08 / 17 by A.Ito	(based on E.Kato's comments)	*/
/* Revised 2010 / 03 / 21 by A.Ito	*/

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
#include<time.h>

/* global variables */
#include"global_var.h"
/* structure */
#include"structure.h"
/* prototypes */
#include"prototype.h"

#define FROWS 10

/*  how to run the model?
 
 <prompt> ./visita
 
*/

/* main simulation roop ******************************************************/
int main(
	int argc, 
	char *argv[]
){
	short zone;
	long f, g, h, l_config, rpert, flag_calc;
	double area_t, area_b, area_r, area_l;
	/* global land area */
	/* file name strings */
	char filename[128];
	char s_date[32];
	char s_case[32];
	char s_config[16];
	char num[4];
	/* file pointer */
	FILE *fp_s[IFILEN];
	FILE *fp_c[4];
	FILE *fp_c2[4];
	FILE *fp_o1[OFILEN], *fp_o2[OFILEN];
	FILE *fp_binout;
	FILE *fp_setting;
	
	/* read configure (instead of arguments) by A.Ito (2009/09/01) ************/
	if((fp_setting = fopen("setting.txt","rt")) == NULL){
	   printf("No configuration file !!!!!!!!!\n");
	   exit(1);
	}
	
	/* config: 1 experiomental scenario ID number (see setting.h) */
	fscanf(fp_setting,"%s %ld", s_config, &l_config);
	printf("config  1: %s %ld\n", s_config, l_config);
	GCM_ID = l_config;
	   if(GCM_ID>=0 && GCM_ID<=9999){
	   ;
	}else{
	   printf("Bad scenario ID specified !!!\n");
	   exit(1);
	}   
	set_gcm_index(s_case);	/* -> vegetdeal.c */
	
	/* config: 2 file identifier (arbitrary phrase such as date, your name, etc.) */
	fscanf(fp_setting,"%s %s", s_config, s_date);
	strcat(s_date, "_");
	printf("config  2: %s %s\n", s_config, s_date);
	   
	/* config: 3 code for radiation sensitivity analysis */
	fscanf(fp_setting,"%s %ld", s_config, &l_config);   
	printf("config  3: %s %ld\n", s_config, l_config);
	SENS_RAD = l_config;
	
	/* config: 4 photosynthesis model 0(Monsi-Saeki) or 1(DePury-Farquhar) */
	fscanf(fp_setting,"%s %ld", s_config, &l_config);
	printf("config  4: %s %ld\n", s_config, l_config);
	DF97 = l_config;
	/* 0: Monsi-Saeki */
	/* 1: de Pury-Farquhar */
	   
	/* config: 5 future solar radiation change */
	fscanf(fp_setting,"%s %ld", s_config, &l_config);
	printf("config  5: %s %ld\n", s_config, l_config);
	CC_R = l_config;
	/* 0: no radiation change */
	/* 1: with radiation change */

	/* config: 6 simple temperature change scenario */
	fscanf(fp_setting,"%s %ld", s_config, &l_config);
	printf("config  6: %s %ld\n", s_config, l_config);
	TEMP_GC = l_config;
		
	/* config: 7 parameter perturbation */
    /* note: no perturbation for PARAM_PTB<=0 */
    if(ENSEMBLE_RUN == 1){
        fscanf(fp_setting,"%s %ld %ld", s_config, &l_config, &rpert);
        printf("config  7: %s %ld %ld\n", s_config, l_config, rpert);
        PARAM_PTB = (long)atol(argv[1]);
        PARAM_ENS = (long)atol(argv[2]);
        EX_CCPL = (long)atol(argv[3]);
    }else{
        fscanf(fp_setting,"%s %ld", s_config, &l_config);
        printf("config  7: %s %ld\n", s_config, l_config);
        //PARAM_PTB = 0;
        PARAM_PTB = l_config;
        PARAM_ENS = 0;
        EX_CCPL = 0;
        rpert = 0;
    }
    
    if(EX_ALBEDO == 1){
        srand((long)atol(argv[1]) + clock()%1000);
    
        snprintf(num, 4, "%03d", (short)atol(argv[1]));
		strcat(s_date, "E");
		strcat(s_date, num);
		strcat(s_date, "_");
        for(f=0;f<N_PARA_ENS;f++){
            f_pert[f] = 0.0;
        }
    }else{
        if(PARAM_PTB >= 1){
            srand((unsigned int)(rpert + clock()%1000));
            rand();
            for(f=0;f<N_PARA_ENS;f++){
                f_pert[f] = 0.0;
                for(g=0;g<12;g++){
                    f_pert[f] += (double)rand() / (double)RAND_MAX;
                }
                f_pert[f] -= 6.0;
                
                if(f_pert[f] > 3.0){
                    f_pert[f] = 3.0;
                }
                if(f_pert[f] < -3.0){
                    f_pert[f] = -3.0;
                }
            }
            strcat(s_date, "E");
            switch(PARAM_PTB){
                case 2: strcat(s_date, "02"); break;
                case 3: strcat(s_date, "03"); break;
                case 4: strcat(s_date, "04"); break;
                case 5: strcat(s_date, "05"); break;
                case 6: strcat(s_date, "06"); break;
                case 7: strcat(s_date, "07"); break;
                case 8: strcat(s_date, "08"); break;
                case 9: strcat(s_date, "09"); break;
                case 10: strcat(s_date, "10"); break;
                case 11: strcat(s_date, "11"); break;
                default: break;
            }
            strcat(s_date, "_");
            snprintf(num, 4, "%03d", (short)PARAM_ENS);
            strcat(s_date, num);
            strcat(s_date, "_");
        }else{
            for(f=0;f<N_PARA_ENS;f++){
                f_pert[f] = 0.0;
            }
        }
    }
    
    /* small carbon flow coupling/decoupling runs */
    switch(EX_CCPL){
        case 1: strcat(s_date, "UC1_"); break;
        case 2: strcat(s_date, "UC2_"); break;
        case 3: strcat(s_date, "UC3_"); break;
        case 4: strcat(s_date, "UC4_"); break;
        case 5: strcat(s_date, "UC5_"); break;
        case 6: strcat(s_date, "UC6_"); break;
        case 7: strcat(s_date, "UC7_"); break;
        case 8: strcat(s_date, "UC8_"); break;
        default: break;
    }
	
	/* config: 8 CH4 experiment */
	fscanf(fp_setting,"%s %ld %ld %ld", s_config, &EX_CH4_1, &EX_CH4_2, &EX_CH4_3);
	printf("config  8: %s %ld %ld %ld\n", s_config, EX_CH4_1, EX_CH4_2, EX_CH4_3);
	
	/* config: 9 simulation area */
	fscanf(fp_setting,"%s %lf %lf %lf %lf", s_config, &area_t, &area_b, &area_l, &area_r);
	printf("config  9: %s %lf %lf %lf %lf\n", s_config, area_t, area_b, area_l, area_r);
    
	/* config: 10 SRM experiment */
	fscanf(fp_setting,"%s %ld", s_config, &EX_SRM);
	printf("config  10: %s %ld\n", s_config, EX_SRM);
	
    /* close setting.txt */
    fclose(fp_setting);
	
	/* open source files *************************************************/
	printf("Open input files...");
	open_input(fp_s, fp_c);	/* -> open_input.c */
	printf("done\n");
	
	/* open result file **************************************************/
	printf("Create output binary file...");
	strcpy(filename, s_date);
	strcat(filename, s_case);
	strcat(filename,"output.bin");
	fp_binout = fopen(filename, "wb");
	printf("done\n");
	
	/**********************************************************************/
	printf("Initialize simulation...");
	/* initialize simulation configuration **/
	/* vegetation type: olson */
	loct.v_type = 1;
	(echar.c3).v_type = 1; 
	(echar.c4).v_type = 1;
	(echar.soil).v_type = 1;
	(mass.c3).v_type = 1; 
	(mass.c4).v_type = 1;
	(mass.soil).v_type = 1;
	(flux.c3).v_type = 1; 
	(flux.c4).v_type = 1;
	(flux.soil).v_type = 1;
	/* vegetation type: cropland */
	loct_agr.v_type = 2;
	(echar_agr.c3).v_type = 2; 
	(echar_agr.c4).v_type = 2;
	(echar_agr.soil).v_type = 2;
	(mass_agr.c3).v_type = 2; 
	(mass_agr.c4).v_type = 2;
	(mass_agr.soil).v_type = 2;
	(flux_agr.c3).v_type = 2; 
	(flux_agr.c4).v_type = 2;
	(flux_agr.soil).v_type = 2;
	
    /* initilization of simulation */
	f_init_sim(&grid);
	printf("done\n");
	
	/* read climate scenario 2010/01/04 (A.Ito) ***********/
	if(NCEP_RUN == 1){
		printf("Reading NCEP climate data...");
		read_ncep_clim(&grid);
	}
	if(GCM_RUN == 1){
		printf("Reading GCM climate projection...");
		read_gcm_clim(fp_c2, &grid);
	}
    
    /* N deposition by CHASER: 2014/11/19 */
    f_read_chaser_ndepo(fp_s, &grid);
    
	printf("done\n");
	
	/*************************************************************************/
	/* latitude loop: north to south **************/
	printf("Start simulation...\n");
	for(f=0;f<N_ROW;f++){
		/* zone number: e.g., 5-deg latitude */
		zone = (short)f/FROWS + 1;
		
		if(f%FROWS == 0){
			/* open output files */	/* -> output_file.c */
			if(CALC_OLSON == 1){
				f_output_file_open(1, zone, s_date, s_case, filename, fp_o1);
			}
			if(CALC_CROP == 1){
				f_output_file_open(2, zone, s_date, s_case, filename, fp_o2);
			}
		}
		
		/* grid area */	/* -> vegetdeal.c */
		grid.area = grid_area(90.0-0.5*(double)f, 89.5-0.5*(double)f, 0.0, 0.5);

		/* longitude loop: west to east ****************************************/
		for(g=0;g<N_COL;g++){
			/* grid position */
			grid.row = f; 
			grid.col = g;	
			
			/* initialize grid condition: data setting */
			f_init_grid(fp_s, &grid);
			
			/* read CRU TS2.X/TS3.X climate data */
			read_cru_clim(fp_c, fp_c2, &grid);
			
			printf("%3ld %3ld: %7.2lf %7.2lf: %2ld %2ld %2ld: %1ld\n", 
				grid.row, grid.col, grid.lat, grid.lon, grid.veg_olson, grid.veg_sage, 
				grid.veg_crop, grid.flag_histdata); /* */
            
            /****/
            /* printf("*************%lf %lf %lf\n", grid.f_biofuel[0], grid.f_biofuel[10], grid.f_biofuel[50]); */
			
			/* header information of output files */
			for(h=0;h<OFILEN;h++){
				if(CALC_OLSON == 1){
					fprintf(fp_o1[h],"%ld %ld %ld %ld %ld\n",
							grid.row, grid.col, grid.veg_olson, grid.veg_sage, grid.flag_histdata);
					
					fprintf(fp_o1[h],"%lf %lf %lf\n", 
							grid.field_cap1, grid.field_cap2, grid.bulkdens);
				}
				if(CALC_CROP == 1){
					/* modified: 2011/02/04 (A.Ito) */
					fprintf(fp_o2[h],"%ld %ld %ld %ld %ld %ld\n",
							grid.row, grid.col, grid.veg_olson, grid.veg_sage, grid.veg_crop, grid.flag_histdata);
					
					fprintf(fp_o2[h],"%lf %lf %lf\n", 
							grid.field_cap1, grid.field_cap2, grid.bulkdens);
				}
			}
			
			/* selection of calculation grids ***************************************/
			flag_calc = 0;
			if(grid.veg_olson!=0 && grid.veg_olson!=33 && grid.flag_histdata==1
					/* && (g+5)%10==0  */
					&& (g+0)%1==0 /* */
					&& grid.lat<=area_t&&grid.lat>=area_b && grid.lon>=area_l&&grid.lon<=area_r
					/* && grid.lat<90.0&&grid.lat>-90.0 && grid.lon>-180.0&&grid.lon<180.0 */
					/* && grid.lat<65.0&&grid.lat>55.0 && grid.lon>65.0&&grid.lon<85.0 */
					/* && grid.lat<=60.75&&grid.lat>=58.25 && grid.lon>=67.25&&grid.lon<=73.75 */   /*DEM*/
					/* && grid.lat<=59.25&&grid.lat>=57.25 && grid.lon>=79.25&&grid.lon<=85.25 */   /*KRS*/
			){
				flag_calc = 1;
			}
			
			/* calculation for lands *******************************************/
            if(flag_calc == 1){
                glat_area[f] += grid.area;  /* latitudinal */
                
                /* Olson map */
                if(CALC_OLSON == 1){
                    /* sequential number */
                    grid.n_olson++;
                    /* total area */
                    go_landarea += grid.area; /* total land */
                    vo_area[grid.veg_olson] += grid.area; /* vegetation */
                    rh_area[grid.reg_g] += grid.area; /* regional */

                    /* clear all parameters ******/
                    f_clear(&grid, &loct, &echar, &mass, &flux);

                    /* initialize vegatation and soil conditions ****/ 
                    initVS(&grid, &loct, &mass, &flux, &echar);		
                
                    /* initialize climate & CO2 conditions ****/
                    f_init_clim(&grid);
                
                    /* initialize location conditions ****/
                    f_init_loct(&grid, &loct, &mass, &flux, &echar);
                    
                    /* initialize stable carbon isotope ****/
                    f_init_c_isotpes(&grid, &flux, &echar, &mass);
                    
                    grid.phase = 0;  /* confirmation */
                    /* spin-up: stabilization roop ***************************/
                    cal_spinup(&grid, &loct, &echar, &mass, &flux, fp_o1);

                    /* snap shot for checking *****/
                    screenshow(&grid, &loct, &mass, &flux, &echar); 
                    
                    /* experiment *******************************************/
                    /* historical: 1901-2000/2009 */
                    /* ISI-MIP: 1950-2099 */
                    cal_historical(&grid, &loct, &echar, &mass, &flux, fp_o1);	

                    /* future: 2001-2100 */
                    if(GCM_RUN==1){
                        cal_projection(&grid, &loct, &echar, &mass, &flux, fp_o1);
                    }
                    
                    printf("\n");
                }
                            
                /* croplands */
                if(CALC_CROP == 1){
                    /* generic cropland */
                    
                    grid.veg_olson = 31;
                    vo_area[grid.veg_olson] += grid.area;
                    
                    /* sequential number */
                    grid.n_crop++;
                    /* total area */

                    /* clear all parameters ******/
                    f_clear(&grid, &loct_agr, &echar_agr, &mass_agr, &flux_agr); 

                    /* initialize vegatation and soil conditions ****/ 
                    initVS(&grid, &loct_agr, &mass_agr, &flux_agr, &echar_agr);
                
                    /* initialize climate conditions ****/
                    f_init_clim(&grid);
                
                    /* initialize location conditions ****/
                    f_init_loct(&grid, &loct_agr, &mass_agr, &flux_agr, &echar_agr);
                
                    /* initialize stable carbon isotope ****/
                    f_init_c_isotpes(&grid, &flux_agr, &echar_agr, &mass_agr);

                    /* spin-up: stabilization roop ***************************/
                    grid.phase = 0;  /* confirmation */
                    cal_spinup(&grid, &loct_agr, &echar_agr, &mass_agr, &flux_agr, fp_o2);

                    /* snap shot for checking *****/
                    screenshow(&grid, &loct_agr, &mass_agr, &flux_agr, &echar_agr);
                    
                    /* experiment **************************************/
                    /* historical: 1901-2000 */
                    /* ISI-MIP: 1950-2099 */
                    cal_historical(&grid, &loct_agr, &echar_agr, &mass_agr, &flux_agr, fp_o2); 

                    /* future: 2001-2100 */
                    if(GCM_RUN){
                        cal_projection(&grid, &loct_agr, &echar_agr, &mass_agr, &flux_agr, fp_o2);
                    }
                    
                    printf("\n");
                }
            }
		}
		/* end of longitudinal loop *************************/
		
		/* output of global total values ******/
		f_glosum_output(s_date, s_case);
		
		/* close output files */
		if(f%FROWS == (FROWS-1)){
			for(h=0;h<OFILEN;h++){
				if(CALC_OLSON == 1){
					fclose(fp_o1[h]);
				}
				if(CALC_CROP == 1){
					fclose(fp_o2[h]);
				}
			}
		}
	}
	/* end of latitudinal loop *************************************/
	
	/* decadal average grid values */
    /* used for IMPRESSIONS analysis: 2015/08/11 by A.Ito */
	fwrite(g_tmp, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 0-4
	fwrite(g_prc, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 5-9
	fwrite(g_swr, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 10-14
	fwrite(g_gpp, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 15-19
	fwrite(g_npp, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 20-24
	fwrite(g_nep, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 25-29
	fwrite(g_pmas, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 30-34
	fwrite(g_smas, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 35-39
	fwrite(g_ch4e_cao, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 40-44
	fwrite(g_ch4o_curry, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 45-49
	fwrite(g_n2oe, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 50-54
	fwrite(g_bbco2, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 55-59
	fwrite(g_ersn, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 60-64
	fwrite(g_isopr, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 65-69
	fwrite(g_sr, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 70-74
	fwrite(g_luc, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 75-79
	fwrite(g_er, sizeof(float), 5*N_ROW*N_COL, fp_binout);		// 80-84
	fwrite(g_ch4ep_cao, sizeof(float), 5*N_ROW*N_COL, fp_binout);	// 85-89

#if C13_GOUT==1	
	fwrite(g_gpp13, sizeof(float), 5*N_ROW*N_COL, fp_binout);	// 90
	fwrite(g_er13, sizeof(float), 5*N_ROW*N_COL, fp_binout);	// 95
	fwrite(g_f13, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 100
	fwrite(g_c13, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 105
	fwrite(g_r13, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 110
	fwrite(g_l13, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 115
	fwrite(g_h13, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 120
#endif
	
#if C14_GOUT==1	
	fwrite(g_gpp14, sizeof(float), 5*N_ROW*N_COL, fp_binout);	// 125
	fwrite(g_er14, sizeof(float), 5*N_ROW*N_COL, fp_binout);	// 130
	fwrite(g_f14, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 135
	fwrite(g_c14, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 140
	fwrite(g_r14, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 145
	fwrite(g_l14, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 150
	fwrite(g_h14, sizeof(float), 5*N_ROW*N_COL, fp_binout);  // 155
#endif

#if PHYS_GOUT==1
	fwrite(g_lai, sizeof(float), 5*N_ROW*N_COL, fp_binout);     // 160 90
	fwrite(g_parb, sizeof(float), 5*N_ROW*N_COL, fp_binout);	// 165 95
	fwrite(g_pard, sizeof(float), 5*N_ROW*N_COL, fp_binout);	// 170 100
	fwrite(g_apar, sizeof(float), 5*N_ROW*N_COL, fp_binout);	// 175 105
	fwrite(g_apar2, sizeof(float), 5*N_ROW*N_COL, fp_binout);	// 180 110
	fwrite(g_aet, sizeof(float), 5*N_ROW*N_COL, fp_binout);     // 185
	fwrite(g_rof, sizeof(float), 5*N_ROW*N_COL, fp_binout);     // 190 120
	fwrite(g_rns, sizeof(float), 5*N_ROW*N_COL, fp_binout);     // 195
	fwrite(g_rnl, sizeof(float), 5*N_ROW*N_COL, fp_binout);     // 200 130
	fwrite(g_sw1, sizeof(float), 5*N_ROW*N_COL, fp_binout);     // 205
	fwrite(g_sw2, sizeof(float), 5*N_ROW*N_COL, fp_binout);     // 210 140
	fwrite(g_snh4, sizeof(float), 5*N_ROW*N_COL, fp_binout);	// 215
	fwrite(g_sno3, sizeof(float), 5*N_ROW*N_COL, fp_binout);	// 220 150
	fwrite(g_rnsd, sizeof(float), 5*N_ROW*N_COL, fp_binout);	// 225 added: 2013/01/10 by A.Ito
#endif
	
#if CH4_WH==1	
	fwrite(g_ch4ep_wh, sizeof(float), 5*N_ROW*N_COL, fp_binout);	// 225
	fwrite(g_ch4ew_wh, sizeof(float), 5*N_ROW*N_COL, fp_binout);	// 230
	fwrite(gm_ch4ep_wh, sizeof(float), 12*N_ROW*N_COL, fp_binout);	// 235
#endif
	
	fclose(fp_binout);
	
	/* close files *************/
	for(h=0;h<4;h++){
		fclose(fp_c[h]);
        
        /* revised 2015/8/12 by A.Ito */
        if(GCM_RUN == 1 && GCM_ID >= 1){
            fclose(fp_c2[h]);
        }
	}
	for(h=0;h<IFILEN;h++){
		fclose(fp_s[h]); 
	}
	
	printf("Simulation ended\n");
	
	return 0;
}
