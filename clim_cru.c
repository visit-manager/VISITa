/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* UEA/CRU data  **********************************************/
/*
New, M., Lister, D., Hulme, M. and Makin, I., 2002. A high-resolution data set of 
 surface climate over global land areas. Climate Research, 21:1-25.
Mitchell, T. D., and P. D. Jones (2005), An improved method of constructing a database 
 of monthly climate observations and associated high-resolution grids, 
 International Journal of Climatology, 25, 693-712.
Harris, I., P. D. Jones, T. J. Osborn, and D. H. Lister, 2014: Updated
 high-resolution grids of monthly climatic observations – the CRU TS3.10 Dataset.
 International Journal of Climatology, 34, 623–642.
*/
void read_cru_clim(
	FILE *fp_c[4], 
	FILE *fp_c2[4],
	struct Grid *grid
){
	long kk[4], flag;
	long f, g, h, i;
	double data, alt, aa, bb, vps;
    double apres, shum, avtas, avpr, drad;
    float r_isimip_data[ASTEP*DL_ISIMIP]; /* corrected: 2012/08/05 */
    float r_gcm_data[ASTEP*DL_ISIMIP]; /* added: 2014/07/31 */

	/*  printf("reading CRU data...");  */
	/* read CRU TS3.0 Vapor-pressure data: 2010/01/04 (A.Ito) */
	
    for(g=0;g<ASTEP;g++){
        grid->hist_cld_b[g] = grid->hist_pre_b[g] = 0.0;
        grid->hist_tmp_b[g] = grid->hist_vap_b[g] = 0.0;
    }

    if(ISIMIP_RUN == 0){
        flag = 0;
        
        /* read CRU TS Cloud data */
        fscanf(fp_c[0],"%ld", &kk[0]);
        if(kk[0]!=0){
            flag++;
            for(h=0;h<DL_CRU;h++){
                for(g=0;g<ASTEP;g++){
                    fscanf(fp_c[0],"%lf", &data);
                    
                    grid->hist_cld[h][g] = data/100.0;
                }
            }
        }
        
        /* read CRU TS precipitation data */
        fscanf(fp_c[1],"%ld", &kk[1]);
        if(kk[1]!=0){
            flag++;
            for(h=0;h<DL_CRU;h++){
                for(g=0;g<ASTEP;g++){
                    fscanf(fp_c[1],"%lf", &data);
                    
                    grid->hist_pre[h][g] = data;
                }
            }
        }
        
        /* read CRU TS temperature data */
        fscanf(fp_c[2],"%ld", &kk[2]);
        if(kk[2]!=0){
            flag++;
            for(h=0;h<DL_CRU;h++){
                for(g=0;g<ASTEP;g++){
                    fscanf(fp_c[2],"%lf", &data);
                    
                    grid->hist_tmp[h][g] = data;
                }
            }
        }
        
        /* read CRU TS vapor-pressure data */
        fscanf(fp_c[3],"%ld", &kk[3]);
        if(kk[3]!=0){
            flag++;
            for(h=0;h<DL_CRU;h++){
                for(g=0;g<ASTEP;g++){
                    fscanf(fp_c[3],"%lf", &data);
                    
                    grid->hist_vap[h][g] = data;
                }
            }
        }
        
        /*  printf("********* %ld\n",cru_flag);  */
        
        /* if valid CRU climate data are all available **/
        if(flag == 4){
            /* data available */
            grid->flag_histdata = 1;
            alt = (grid->topo>=0.0)?grid->topo:0.0; 

            /* base climate (average 1971 - 2000) ******************/
            for(f=0;f<30;f++){
                for(g=0;g<ASTEP;g++){
                    grid->hist_cld_b[g] += grid->hist_cld[70+f][g]/30.0;
                    grid->hist_pre_b[g] += grid->hist_pre[70+f][g]/30.0;
                    grid->hist_vap_b[g] += grid->hist_vap[70+f][g]/30.0;
                    grid->hist_tmp_b[g] += grid->hist_tmp[70+f][g]/30.0;
                }
            }
        }else{
            /* unavailable CRU TS data, for example on ocean */
            grid->flag_histdata = 0;
        }
    }else if(ISIMIP_RUN == 1 ||ISIMIP_RUN == 2 ||ISIMIP_RUN == 3){
        
        /* ISI-MIP: 2012/06/27 by A.Ito ****************/
        /* also for ICARUS */
        /* 1950-1979-detrended: spi-up */
        /* 1950-2005:           historical */
        /* 2006-2099:           future projection */
        
        /* PLUME (ISI-MIP2): 2014/07/31 by A.Ito ****************/
        /* 1901-1930-detrended: spi-up */
        /* 1901-2005:           historical */
        
        /* ISI-MIP2: 2014/11/30 by A.Ito ****************/
        /* 1901-1930-detrended: spi-up */
        /* 1901-2010:           historical */

        /* ait tempetaure, deg-C */
        fread(r_isimip_data,sizeof(float),ASTEP*DL_ISIMIP, fp_c[0]);
        avtas = 0.0;
        for(h=0;h<DL_ISIMIP;h++){
            for(g=0;g<ASTEP;g++){
                /* K => degC */
                grid->hist_tmp[h][g] = (double)r_isimip_data[h*ASTEP + g] - ZAT;
                avtas += grid->hist_tmp[h][g] / (double)DL_ISIMIP / (double)ASTEP;
            }
        }
        
        /* precipitation, mm month-1 */
        fread(r_isimip_data,sizeof(float),ASTEP*DL_ISIMIP, fp_c[1]);
        avpr = 0.0;
        for(h=0;h<DL_ISIMIP;h++){
            for(g=0;g<ASTEP;g++){
                grid->hist_pre[h][g] = (double)r_isimip_data[h*ASTEP+g] * (double)MDN[g] *24.0*3600.0;
                avpr += grid->hist_pre[h][g] / (double)DL_ISIMIP;
                grid->hist_pre[h][g] = (grid->hist_pre[h][g]>0.0)?grid->hist_pre[h][g]:0.0;
            }
        }
        
        /* effective data availability */
        if(avtas > -50.0 && avtas < 50.0 && avpr >= 0.0 && avpr < 10000.0){
            grid->flag_histdata = 1;
        }else{
            grid->flag_histdata = 0;
        }
        
        fread(r_isimip_data, sizeof(float), ASTEP*DL_ISIMIP, fp_c[2]);
        if(grid->flag_histdata == 1){
        
            /* relative humidity (%) => vapor pressure (hPa) */
            for(h=0;h<DL_ISIMIP;h++){
                for(g=0;g<ASTEP;g++){
                
                    /* relative humidity to vapor pressure */
                    /* revided by A.Ito (2012/06/28) */
                    
                    if(ISIMIP_RUN == 1 ||ISIMIP_RUN == 2){
                        /* saturated water vapor pressure */
                        if(grid->hist_tmp[h][g] > 0.0){ /* at water surface */
                            vps = 6.1078*pow(10.0, (7.5 * grid->hist_tmp[h][g])/(237.3 + grid->hist_tmp[h][g]));
                        }else{ /* at ice surface */  /*  if(grid->tmp_2m[grid->m]<=0.0) */
                            vps = 6.1078*pow(10.0, (9.5 * grid->hist_tmp[h][g])/(265.3 + grid->hist_tmp[h][g]));
                        }
                        vps = (vps>=0.0)?vps:0.0;
                        
                        grid->hist_vap[h][g] = vps * (double)r_isimip_data[h*ASTEP+g] / 100.0;
                    }else if(ISIMIP_RUN == 3){
                        /* altitude */
                        alt = (grid->topo>=0.0)?grid->topo:0.0;
                        
                        /* air pressure, hPa */
                        apres = 1013.25*exp(-1.0*(28.964*0.001) * GAC * alt / (UGC*(grid->hist_tmp[h][g] + ZAT)));

                        /* vapour pressure, hPa */
                        shum = (double)r_isimip_data[h*ASTEP+g];
                        if(shum < 0.0){
                            shum = 0.0;
                        }
                        grid->hist_vap[h][g] = apres * shum /(0.622 + 0.378 * shum);
                    }
                    
                    grid->hist_vap[h][g] = (grid->hist_vap[h][g]>0.0)?grid->hist_vap[h][g]:0.0;
                }
            }
        }else{
            for(h=0;h<DL_ISIMIP;h++){
                for(g=0;g<ASTEP;g++){
                    grid->hist_vap[h][g] = 0.0;
                }
            }
        }
        
        fread(r_isimip_data, sizeof(float), ASTEP*DL_ISIMIP, fp_c[3]);
        if(grid->flag_histdata == 1){
            
            /* radiation => cloudiness, fraction */
            for(h=0;h<DL_ISIMIP;h++){
                for(g=0;g<ASTEP;g++){
                    /* average downward-shortwave radiation */
                    /* 2012/06/29 by A.Ito */
                    drad = 0.0;
                    grid->m = g;
                    for (f=0;f<24;f++) {
                        drad += f_top_rad(grid, -180+15*f) / 24.0; 
                    }
                    
                    /* inverse estimation of cloudiness */
                    if(drad > 0.0){
                        aa = (double)r_isimip_data[h*ASTEP + g] / drad;
                    }else{
                        aa = 0.0;
                    }
                    bb = (0.8964 - aa) / 0.5392;
                    if(bb < 0.0){
                        bb = 0.0;
                    }
                    if(bb > 1.0){
                        bb = 1.0;
                    }
                
                    grid->hist_cld[h][g] = bb;
                    grid->hist_cld[h][g] = (grid->hist_cld[h][g]>0.0)?grid->hist_cld[h][g]:0.0;
                }
            }
        }else{
            for(h=0;h<DL_ISIMIP;h++){
                for(g=0;g<ASTEP;g++){
                    grid->hist_cld[h][g] = 0.0;
                }
            }
        }
        
        /* climatology *************************************/
        for(f=0;f<30;f++){
            for(g=0;g<ASTEP;g++){
                /* ISI-MIP1: average of 1951–1980 (historial detrended) data */
                /* ISI-MIP2: average of 1901–1930 (historial detrended) data */
                grid->hist_cld_b[g] += grid->hist_cld[f][g] / 30.0;
                grid->hist_pre_b[g] += grid->hist_pre[f][g] / 30.0;
                grid->hist_vap_b[g] += grid->hist_vap[f][g] / 30.0;
                grid->hist_tmp_b[g] += grid->hist_tmp[f][g] / 30.0;
            }
        }
        
        /* base temperature for permafrost dynamics */
        for(f=0;f<20;f++){
            for(g=0;g<ASTEP;g++){
                grid->tmp_base_permaforst += grid->hist_tmp[f + 10][g] / 240.0;
            }
        }
    }
    
    /* projection data ************************************************/
    if(ISIMIP_RUN == 2){
        /* PLUME (ISI-MIP2): 2014/07/31 by A.Ito ****************/
        /* 2006-2091:   future projection */
        
        /* ait tempetaure, deg-C */
        fread(r_gcm_data,sizeof(float), ASTEP*DL_GCM, fp_c2[0]);
        
        for(h=0;h<DL_GCM;h++){
            for(g=0;g<ASTEP;g++){
                /* K => degC */
                grid->proj_tmp2m[h][g][0][0] = (double)r_gcm_data[h*ASTEP + g] - ZAT;
            }
        }
        
        /* precipitation, mm month-1 */
        fread(r_gcm_data,sizeof(float), ASTEP*DL_GCM, fp_c2[1]);
        avpr = 0.0;
        for(h=0;h<DL_GCM;h++){
            for(g=0;g<ASTEP;g++){
                grid->proj_prec[h][g][0][0] = (double)r_gcm_data[h*ASTEP+g] * (double)MDN[g] *24.0*3600.0;
                grid->proj_prec[h][g][0][0] = (grid->proj_prec[h][g][0][0]>0.0)?grid->proj_prec[h][g][0][0]:0.0;
            }
        }
        
        /* relative humidity (%) => vapor pressure (hPa) */
        fread(r_gcm_data,sizeof(float), ASTEP*DL_GCM, fp_c2[2]);
        if(grid->flag_histdata == 1){
            for(h=0;h<DL_GCM;h++){
                for(g=0;g<ASTEP;g++){
                
                    /* relative humidity to vapor pressure */
                    /* revided by A.Ito (2012/06/28) */
                    
                    /* saturated water vapor pressure */
                    if(grid->proj_tmp2m[h][g][0][0] > 0.0){ /* at water surface */
                        vps = 6.1078*pow(10.0, (7.5*grid->proj_tmp2m[h][g][0][0])/(237.3 + grid->proj_tmp2m[h][g][0][0]));
                    }else{ /* at ice surface */  /*  if(grid->tmp_2m[grid->m]<=0.0) */
                        vps = 6.1078*pow(10.0, (9.5*grid->proj_tmp2m[h][g][0][0])/(265.3 + grid->proj_tmp2m[h][g][0][0]));
                    }
                    vps = (vps>=0.0)?vps:0.0;
        
                    grid->proj_hum[h][g][0][0] = vps * (double)r_gcm_data[h*ASTEP+g]/100.0;
                    grid->proj_hum[h][g][0][0] = (grid->proj_hum[h][g][0][0]>0.0)?grid->proj_hum[h][g][0][0]:0.0;
                }
            }
        }else{
            for(h=0;h<DL_GCM;h++){
                for(g=0;g<ASTEP;g++){
                    grid->proj_hum[h][g][0][0] = 0.0;
                }
            }
        }
        
        /* radiation => cloudiness, fraction */
        fread(r_gcm_data, sizeof(float), ASTEP*DL_GCM, fp_c2[3]);
        if(grid->flag_histdata == 1){
            for(h=0;h<DL_GCM;h++){
                for(g=0;g<ASTEP;g++){
                    grid->proj_rad[h][g][0][0] = (double)r_gcm_data[h*ASTEP+g];
                    grid->proj_rad[h][g][0][0] = (grid->proj_rad[h][g][0][0]>0.0)?grid->proj_rad[h][g][0][0]:0.0;
                }
            }
        }else{
            for(h=0;h<DL_GCM;h++){
                for(g=0;g<ASTEP;g++){
                    grid->proj_rad[h][g][0][0] = 0.0;
                }
            }
        }
        
        /**********************************************************/
        for(f=0;f<30;f++){ /* 2006-2035 */
            for(g=0;g<ASTEP;g++){
                for(h=0;h<GCM_R;h++){
                    for(i=0;i<GCM_C;i++){
                        if(f == 0){
                            grid->proj_tmp2m_b[g][h][i] = 0.0;
                            grid->proj_prec_b[g][h][i] = 0.0;
                            grid->proj_hum_b[g][h][i] = 0.0;
                            grid->proj_rad_b[g][h][i] = 0.0;
                        }
                    
                        grid->proj_tmp2m_b[g][h][i] += grid->proj_tmp2m[f][g][h][i]/30.0;
                        grid->proj_prec_b[g][h][i] += grid->proj_prec[f][g][h][i]/30.0;
                        grid->proj_hum_b[g][h][i] += grid->proj_hum[f][g][h][i]/30.0;
                        grid->proj_rad_b[g][h][i] += grid->proj_rad[f][g][h][i]/30.0;
                    }
                }
            }
        }
        
    }
	
	/*  printf("done\n");  */
}
