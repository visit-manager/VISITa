/*	VISIT: Vegetation Integrative SImulation Tool						*/
/*  Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/*  Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/*  Carbon cycle, erosion, biomass burning, land-use change,			*/
/*  CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/*  main roop for 0.5-degree version of Sim-CYCLE				*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* land use change ******************************************/
void f_cult_luc(
	struct Grid *grid
){
	if(grid->climy<=1990){
		grid->f_crop = grid->fcrop[grid->climy - 1700];
	}else{
		
		if(LANDUSE==1){
			/* fixed 1900 */
			grid->f_crop = grid->fcrop[200]; /*  */
		}else if(LANDUSE==2){
			/* fixed 1990 */
			grid->f_crop = grid->fcrop[290]; /*  */
		}else if(LANDUSE==3){
			/* conventional scenario : ex b, d */
			grid->f_crop = grid->fcrop[290] + grid->f_crop_trend*(double)(grid->climy - 1990);
		}else if(LANDUSE==4){
			/* high scenario : ex e */
			if(grid->f_crop_trend>0.0){
				grid->f_crop = grid->fcrop[290] + grid->f_crop_trend*(1.0+0.01*((double)(grid->climy - 1990))) *(double)(grid->climy - 1990);
			}else{
				grid->f_crop = grid->fcrop[290];
			}
		}else if(LANDUSE==5){
			grid->f_crop = grid->fcrop[290] + ((grid->fcrop3_image[grid->climy - 1990]+grid->fcrop4_image[grid->climy - 1990]) 
				- (grid->fcrop3_image[0]+grid->fcrop4_image[0]))/100.0;
		}
	}
	
	/*********************************/
	if(LANDUSE==0){
		/* Natural condition */
		grid->f_crop = 0.0;
	}
	
	if(grid->f_crop < 0.0){
		grid->f_crop = 0.0;
	}
	if(grid->f_crop > 0.975){
		grid->f_crop = 0.975;
	}

}

/* emission from land use change *********************************/
void f_luc_emit(
	struct Grid *grid, 
	struct Mass *mass, 
	struct Flux *flux
){
	long f;
	double f_luc;
	double fe_conv;		/* fraction of conversion flux (1year) */
	double fe_ten;		/* fraction of 10-year pool flux */
	double fe_hund;		/* fraction of 100-year pool flux */
	double fe_detr;		/* fraction of detritus flux */
	
	switch(grid->veg_sage){
		case 1: case 2:
			fe_detr = 0.33;
			fe_conv = 0.40;
			fe_ten = 0.27;
			fe_hund = 0.0;
			break;
		case 3: case 4: case 5: case 6: case 7: case 8:
			fe_detr = 0.33;
			fe_conv = 0.40;
			fe_ten = 0.20;
			fe_hund = 0.07;
			break;
		case 10: case 13: case 14: case 15:
			fe_detr = 0.50;
			fe_conv = 0.50;
			fe_ten = 0.0;
			fe_hund = 0.0;
			break;
		case 9: case 11: case 12:
			fe_detr = 0.50;
			fe_conv = 0.40;
			fe_ten = 0.10;
			fe_hund = 0.0;
			break;
		default:
			fe_detr = 0.50;
			fe_conv = 0.40;
			fe_ten = 0.10;
			fe_hund = 0.0;
			break;
	}
	
	/* spin-up: fluxes for 1801-1900 **/
	if(grid->phase==0){
		for(f=1891;f<=1900;f++){
			f_luc = grid->fcrop[f-1700]-grid->fcrop[f-1700-1];
			
			if(f_luc > 0.0){
				fe_ten = f_luc * ((mass->plant).fol + (mass->plant).stm) * fe_ten/(fe_conv+fe_ten+fe_hund);
				
				flux->detr_ten[1900-f] = fe_ten;
			}else{
				flux->detr_ten[1900-f] = 0.0;
			}
		}
		
		for(f=1801;f<=1900;f++){
			f_luc = grid->fcrop[f-1700]-grid->fcrop[f-1700-1];

			if(f_luc > 0.0){
				fe_hund = f_luc * ((mass->plant).fol + (mass->plant).stm) * fe_hund/(fe_conv+fe_ten+fe_hund);
				
				flux->detr_hund[1900-f] = fe_hund;
			}else{
				flux->detr_hund[1900-f] = 0.0;
			}
		}
	}
	
	flux->lu_ten = 0.0;
	for(f=1;f<10;f++){
		flux->lu_ten += 0.1 * flux->detr_ten[f-1];
		flux->detr_ten[f] = flux->detr_ten[f-1];
	}
	
	flux->lu_hund = 0.0;
	for(f=1;f<100;f++){
		flux->lu_hund += 0.01 * flux->detr_hund[f-1];
		flux->detr_hund[f] = flux->detr_hund[f-1];
	}
	
	/* annual land use change */
	f_luc = grid->f_crop - grid->f_crop_p;
	
	if(f_luc > 0.0){ /* deforested */
		fe_detr = f_luc * 0.2*(mass->plant).rot;
		fe_conv = f_luc * ((mass->plant).fol + (mass->plant).stm + 0.8*(mass->plant).rot) * fe_conv/(fe_conv+fe_ten+fe_hund);
		fe_ten = f_luc * ((mass->plant).fol + (mass->plant).stm + 0.8*(mass->plant).rot) * fe_ten/(fe_conv+fe_ten+fe_hund);
		fe_hund = f_luc * ((mass->plant).fol + (mass->plant).stm + 0.8*(mass->plant).rot) * fe_hund/(fe_conv+fe_ten+fe_hund);
		
		flux->lu_conv = fe_conv;
		
		flux->detr_ten[0] = fe_ten;
		flux->lu_ten += 0.1 * fe_ten;

		flux->detr_hund[0] = fe_hund;
		flux->lu_hund += 0.01 * fe_hund;
	}else{
		flux->lu_conv = 0.0;
		flux->detr_ten[0] = 0.0;
		flux->detr_hund[0] = 0.0;
	}
	
	/* maximum */
	/*  if(f_luc > 0.0){
		flux->lu_conv = f_luc * ((mass->plant).fol + (mass->plant).stm + (mass->plant).rot + (mass->soil).ltr + (mass->soil).msl);
		flux->lu_ten = 0.0;
		flux->lu_hund =0.0;
	}else{
		flux->lu_conv = 0.0;
		flux->lu_ten = 0.0;
		flux->lu_hund =0.0;
	}  */
}

