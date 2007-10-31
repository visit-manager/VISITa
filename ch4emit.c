/*	VISIT: Vegetation Integrative SImulation Tool						*/
/*  Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/*  Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/*  Carbon cycle, erosion, biomass burning, land-use change,			*/
/*  CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

#define NOTICE 0

/* CH4 emission model by Cao **********************************/
/*
Cao, M., Marshall, S. and Gregson, K., 1996. Global carbon exchange and 
methane emissions from natural wetlands: Application of a process-based model. 
Journal of Geophysical Research 101, 14399-14414.
*/
void f_ch4_emit_cao(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Flux *flux
){
	double wtable;
	double f_temp, f_wtable;
	double hr_decomp;
	
	/* soil decomposition rate. Mg C ha-1 month-1 */
	hr_decomp = (flux->soil).rS[grid->m];
	if(hr_decomp<0.0){
		hr_decomp = 0.0;
	}
	
	/* temperature (deg C) coefficient */
	f_temp = exp(grid->tmp10_soil[grid->m] * 0.0693) / 7.996;
	if(f_temp<0.0){
		f_temp = 0.0;
	}
	
	/* CH4 emission, mg CH4 m-2 month-1 */

	/* water table (cm relative to surface) coefficient */
	/* wetland */
	wtable = -10.0;
	f_wtable = 0.383 * exp(0.096 * wtable);
	if(f_wtable<0.0){
		f_wtable = 0.0;
	}
	(flux->soil).ch4emit_wetland_cao[grid->m] = hr_decomp * f_temp * f_wtable;
	(flux->soil).ch4emit_wetland_cao[grid->m] *= 16.0/12.0 * 1000000000.0 / 10000.0;
	(flux->soil).ch4emit_wetland_cao[grid->m] *= grid->f_wetland;
	
	/* paddy field */
	if(grid->tmp_sfc[grid->m] > 15.0){
		wtable = 0.0;	
	}else{
		wtable = -20.0;
	}
	f_wtable = 0.383 * exp(0.096 * wtable);
	if(f_wtable<0.0){
		f_wtable = 0.0;
	}
	(flux->soil).ch4emit_paddy_cao[grid->m] = hr_decomp * f_temp * f_wtable;
	(flux->soil).ch4emit_paddy_cao[grid->m] *= 16.0/12.0 * 1000000000.0 / 10000.0;
	(flux->soil).ch4emit_paddy_cao[grid->m] *= grid->f_paddy;
}
