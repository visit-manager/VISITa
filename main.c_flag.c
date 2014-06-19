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

/* main simulation roop ******************************************************/
int main(
	int argc, 
	char *argv[]
){
    long f,g,h;
	FILE *fp_s[IFILEN];
	FILE *fp_c[4], *fpo;
	
	/* open source files *************************************************/
	printf("Open input files...");
	open_input(fp_s, fp_c);	/* -> open_input.c */
	printf("done\n");
    
    fpo = fopen("sim_flag.txt","wt");
	
	/*************************************************************************/
	for(f=0;f<N_ROW;f++){
		for(g=0;g<N_COL;g++){
			/* grid position */
			grid.row = f; 
			grid.col = g;	
			
			/* read CRU TS2.X/TS3.X climate data */
			read_cru_clim(fp_c, &grid);
			
            fprintf(fpo,"%ld %ld %ld\n",f,g,grid.flag_histdata);
        }
	}
	
	/* close files *************/
	for(h=0;h<IFILEN;h++){
		fclose(fp_s[h]); 
	}
	for(h=0;h<4;h++){
		fclose(fp_c[h]); 
	}
    
    fclose(fpo);
	
	return 0;
}
