/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed  in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/**************************************************************************************/
void read_ncep_clim(
	struct Grid  *grid
){
	long f, g, h, i;
	long yr, mon;
	FILE *fp_t, *fp_p, *fp_h, *fp_r;
	
	/* if( (fp_t=fopen("./data/air.2m.gauss.mon4809.dat","rt"))==NULL ){  
		printf("No air.2m.gauss.mon4809.dat\n");  exit(1); }
	if( (fp_p=fopen("./data/prate.sfc.gauss.mon4809.dat","rt"))==NULL ){  
		printf("No prate.sfc.gauss.mon4809.dat\n");  exit(1); }
	if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4809.dat","rt"))==NULL ){  
		printf("No vpres.sfc.gauss.mon4809.dat\n");  exit(1); }
	if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4809.dat","rt"))==NULL ){  
		printf("No tcdc.eatm.gauss.mon4809.dat\n");  exit(1); } */
	
	/* if( (fp_t=fopen("./data/air.2m.gauss.mon4810.dat","rt"))==NULL ){  
		printf("No air.2m.gauss.mon4810.dat\n");  exit(1); }
	if( (fp_p=fopen("./data/prate.sfc.gauss.mon4810.dat","rt"))==NULL ){  
		printf("No prate.sfc.gauss.mon4810.dat\n");  exit(1); }
	if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4810.dat","rt"))==NULL ){  
		printf("No vpres.sfc.gauss.mon4810.dat\n");  exit(1); }
	if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4810.dat","rt"))==NULL ){  
		printf("No tcdc.eatm.gauss.mon4810.dat\n");  exit(1); } */
	
    /* extension to 2011: 2012/01/26  */
	/* if( (fp_t=fopen("./data/air.2m.gauss.mon4811.dat","rt"))==NULL ){
		printf("No air.2m.gauss.mon4811.dat\n");  exit(1); }
	if( (fp_p=fopen("./data/prate.sfc.gauss.mon4811.dat","rt"))==NULL ){  
		printf("No prate.sfc.gauss.mon4811.dat\n");  exit(1); }
	if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4811.dat","rt"))==NULL ){  
		printf("No vpres.sfc.gauss.mon4811.dat\n");  exit(1); }
	if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4811.dat","rt"))==NULL ){  
		printf("No tcdc.eatm.gauss.mon4811.dat\n");  exit(1); } */
	
	/* if( (fp_t=fopen("./data/air.2m.gauss.mon4812.dat","rt"))==NULL ){
		printf("No air.2m.gauss.mon4812.dat\n");  exit(1); }
	if( (fp_p=fopen("./data/prate.sfc.gauss.mon4812.dat","rt"))==NULL ){  
		printf("No prate.sfc.gauss.mon4812.dat\n");  exit(1); }
	if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4812.dat","rt"))==NULL ){  
		printf("No vpres.sfc.gauss.mon4812.dat\n");  exit(1); }
	if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4812.dat","rt"))==NULL ){  
		printf("No tcdc.eatm.gauss.mon4812.dat\n");  exit(1); } */
	
	/* if( (fp_t=fopen("./data/air.2m.gauss.mon4813.dat","rt"))==NULL ){
		printf("No air.2m.gauss.mon4813.dat\n");  exit(1); }
	if( (fp_p=fopen("./data/prate.sfc.gauss.mon4813.dat","rt"))==NULL ){  
		printf("No prate.sfc.gauss.mon4813.dat\n");  exit(1); }
	if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4813.dat","rt"))==NULL ){  
		printf("No vpres.sfc.gauss.mon4813.dat\n");  exit(1); }
	if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4813.dat","rt"))==NULL ){  
		printf("No tcdc.eatm.gauss.mon4813.dat\n");  exit(1); } */
	
    /* added 2014: 2015/01/02  */
	/* if( (fp_t=fopen("./data/air.2m.gauss.mon4814.dat","rt"))==NULL ){
		printf("No air.2m.gauss.mon4814.dat\n");  exit(1); }
	if( (fp_p=fopen("./data/prate.sfc.gauss.mon4814.dat","rt"))==NULL ){
		printf("No prate.sfc.gauss.mon4814.dat\n");  exit(1); }
	if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4814.dat","rt"))==NULL ){
		printf("No vpres.sfc.gauss.mon4814.dat\n");  exit(1); }
	if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4814.dat","rt"))==NULL ){
		printf("No tcdc.eatm.gauss.mon4814.dat\n");  exit(1); } */
	
    /* added 2015: 2016/01/04  */
	/* if( (fp_t=fopen("./data/air.2m.gauss.mon4815.dat","rt"))==NULL ){
		printf("No air.2m.gauss.mon4814.dat\n");  exit(1); }
	if( (fp_p=fopen("./data/prate.sfc.gauss.mon4815.dat","rt"))==NULL ){
		printf("No prate.sfc.gauss.mon4814.dat\n");  exit(1); }
	if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4815.dat","rt"))==NULL ){
		printf("No vpres.sfc.gauss.mon4814.dat\n");  exit(1); }
	if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4815.dat","rt"))==NULL ){
		printf("No tcdc.eatm.gauss.mon4814.dat\n");  exit(1); } */
	
    /* added 2016: 2017/01/03  */
	/* if( (fp_t=fopen("./data/air.2m.gauss.mon4816.dat","rt"))==NULL ){
		printf("No air.2m.gauss.mon4816.dat\n");  exit(1); }
	if( (fp_p=fopen("./data/prate.sfc.gauss.mon4816.dat","rt"))==NULL ){
		printf("No prate.sfc.gauss.mon4816.dat\n");  exit(1); }
	if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4816.dat","rt"))==NULL ){
		printf("No vpres.sfc.gauss.mon4816.dat\n");  exit(1); }
	if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4816.dat","rt"))==NULL ){
		printf("No tcdc.eatm.gauss.mon4816.dat\n");  exit(1); } */
    
    /* added 2017: 2018/01/03  */
    /* if( (fp_t=fopen("./data/air.2m.gauss.mon4817.dat","rt"))==NULL ){
        printf("No air.2m.gauss.mon4817.dat\n");  exit(1); }
    if( (fp_p=fopen("./data/prate.sfc.gauss.mon4817.dat","rt"))==NULL ){
        printf("No prate.sfc.gauss.mon4817.dat\n");  exit(1); }
    if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4817.dat","rt"))==NULL ){
        printf("No vpres.sfc.gauss.mon4817.dat\n");  exit(1); }
    if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4817.dat","rt"))==NULL ){
        printf("No tcdc.eatm.gauss.mon4817.dat\n");  exit(1); } */

    /* added 2018: 2019/01/03  */
    /* if( (fp_t=fopen("./data/air.2m.gauss.mon4818.dat","rt"))==NULL ){
        printf("No air.2m.gauss.mon4818.dat\n");  exit(1); }
    if( (fp_p=fopen("./data/prate.sfc.gauss.mon4818.dat","rt"))==NULL ){
        printf("No prate.sfc.gauss.mon4818.dat\n");  exit(1); }
    if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4818.dat","rt"))==NULL ){
        printf("No vpres.sfc.gauss.mon4818.dat\n");  exit(1); }
    if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4818.dat","rt"))==NULL ){
        printf("No tcdc.eatm.gauss.mon4818.dat\n");  exit(1); } */

    /* added 2019: 2020/01/02  */
    /* if( (fp_t=fopen("./data/air.2m.gauss.mon4819.dat","rt"))==NULL ){
        printf("No air.2m.gauss.mon4819.dat\n");  exit(1); }
    if( (fp_p=fopen("./data/prate.sfc.gauss.mon4819.dat","rt"))==NULL ){
        printf("No prate.sfc.gauss.mon4819.dat\n");  exit(1); }
    if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4819.dat","rt"))==NULL ){
        printf("No vpres.sfc.gauss.mon4819.dat\n");  exit(1); }
    if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4819.dat","rt"))==NULL ){
        printf("No tcdc.eatm.gauss.mon4819.dat\n");  exit(1); } */

    /* added 2020: 2021/01/XX  */
    /* if( (fp_t=fopen("./data/air.2m.gauss.mon4820.txt","rt"))==NULL ){
        printf("No air.2m.gauss.mon4820.dat\n");  exit(1); }
    if( (fp_p=fopen("./data/prate.sfc.gauss.mon4820.txt","rt"))==NULL ){
        printf("No prate.sfc.gauss.mon4820.dat\n");  exit(1); }
    if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4820.txt","rt"))==NULL ){
        printf("No vpres.sfc.gauss.mon4820.dat\n");  exit(1); }
    if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4820.txt","rt"))==NULL ){
        printf("No tcdc.eatm.gauss.mon4820.dat\n");  exit(1); } */

    /* added 2021: 2022/01/04  */
    /* if( (fp_t=fopen("./data/air.2m.gauss.mon4821.txt","rt"))==NULL ){
        printf("No air.2m.gauss.mon4821.dat\n");  exit(1); }
    if( (fp_p=fopen("./data/prate.sfc.gauss.mon4821.txt","rt"))==NULL ){
        printf("No prate.sfc.gauss.mon4821.dat\n");  exit(1); }
    if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4821.txt","rt"))==NULL ){
        printf("No vpres.sfc.gauss.mon4821.dat\n");  exit(1); }
    if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4821.txt","rt"))==NULL ){
        printf("No tcdc.eatm.gauss.mon4821.dat\n");  exit(1); } */

    /* added 2022: 2023/02/08  */
    /* if( (fp_t=fopen("./data/air.2m.gauss.mon4822.txt","rt"))==NULL ){
        printf("No air.2m.gauss.mon4822.dat\n");  exit(1); }
    if( (fp_p=fopen("./data/prate.sfc.gauss.mon4822.txt","rt"))==NULL ){
        printf("No prate.sfc.gauss.mon4822.dat\n");  exit(1); }
    if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4822.txt","rt"))==NULL ){
        printf("No vpres.sfc.gauss.mon4822.dat\n");  exit(1); }
    if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4822.txt","rt"))==NULL ){
        printf("No tcdc.eatm.gauss.mon4822.dat\n");  exit(1); } */

    /* added 2023: 2024/01/05  */
    if( (fp_t=fopen("./data/air.2m.gauss.mon4823.txt","rt"))==NULL ){
        printf("No air.2m.gauss.mon4822.dat\n");  exit(1); }
    if( (fp_p=fopen("./data/prate.sfc.gauss.mon4823.txt","rt"))==NULL ){
        printf("No prate.sfc.gauss.mon4822.dat\n");  exit(1); }
    if( (fp_h=fopen("./data/vpres.sfc.gauss.mon4823.txt","rt"))==NULL ){
        printf("No vpres.sfc.gauss.mon4822.dat\n");  exit(1); }
    if( (fp_r=fopen("./data/tcdc.eatm.gauss.mon4823.txt","rt"))==NULL ){
        printf("No tcdc.eatm.gauss.mon4822.dat\n");  exit(1); }

	for(f=0;f<DL_NCEP;f++){   /*  1948-20XX */
		for(g=0;g<ASTEP;g++){
			fscanf(fp_t,"%ld %ld", &yr, &mon);
			fscanf(fp_p,"%ld %ld", &yr, &mon);
			fscanf(fp_h,"%ld %ld", &yr, &mon);
			fscanf(fp_r,"%ld %ld", &yr, &mon);
			for(h=0;h<94;h++){
				for(i=0;i<192;i++){
					fscanf(fp_t,"%f", &(grid->ncep_tmp2m[f][g][h][i]));
					fscanf(fp_p,"%f", &(grid->ncep_prate[f][g][h][i]));
					fscanf(fp_h,"%f", &(grid->ncep_vpres[f][g][h][i]));
					fscanf(fp_r,"%f", &(grid->ncep_tcdc[f][g][h][i]));
					
					grid->ncep_tmp2m[f][g][h][i] -= ZAT;
					grid->ncep_tcdc[f][g][h][i] *= 0.01;
				}
			}
		}
	}
	
	/* average climatology 1970-1999 */
	/* average climatology 1980-2009 */
	for(f=0;f<30;f++){
		for(g=0;g<ASTEP;g++){
			for(h=0;h<94;h++){
				for(i=0;i<192;i++){
					if(f==0){
						grid->ncep_tmp2m_b[g][h][i] = 0.0;
						grid->ncep_prate_b[g][h][i] = 0.0;
						grid->ncep_vpres_b[g][h][i] = 0.0;
						grid->ncep_tcdc_b[g][h][i] = 0.0;
					}
					
					grid->ncep_tmp2m_b[g][h][i] += grid->ncep_tmp2m[f+32][g][h][i]/30.0;
					grid->ncep_prate_b[g][h][i] += grid->ncep_prate[f+32][g][h][i]/30.0;
					grid->ncep_vpres_b[g][h][i] += grid->ncep_vpres[f+32][g][h][i]/30.0;
					grid->ncep_tcdc_b[g][h][i] += grid->ncep_tcdc[f+32][g][h][i]/30.0;
				}
			}
		}
	}
	
	fclose(fp_t);
	fclose(fp_p);
	fclose(fp_h);
	fclose(fp_r);
}
