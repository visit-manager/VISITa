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

/* CH4 emission model by Cao **********************************************/
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
	double wtable;	/* water table, cm */
	double f_temp, f_wtable, f_wtable_lake;
	double hr_decomp, gpp_factor;
	
	/* soil decomposition rate. Mg C ha-1 month-1 */
	hr_decomp = (flux->soil).hr[grid->m];
	if(hr_decomp<0.0){
		hr_decomp = 0.0;
	}
	
	/* temperature (deg C) coefficient */
	/* eq.7 */
	f_temp = exp(grid->tmp10_soil[grid->m] * 0.0693) / 7.996;
	if(f_temp<0.0){
		f_temp = 0.0;
	}
	
	/* CH4 emission, mg CH4 m-2 month-1 */

	/* water table (cm relative to surface) coefficient */
	/* wetland *********************************/
	/* eq.6 */
	/* wtable = 5.0; */
	/* f_wtable = 0.383 * (0.5*exp(0.096 * wtable) + 0.5*exp(0.096 * -5.0)); OLD */
	/* revised by A.Ito (2009/07/13) **/
	f_wtable = 0.383 * (grid->inundation_ssmi[grid->m]*exp(0.096 * 0.0) 
						+ (1.0-grid->inundation_ssmi[grid->m])*exp(0.096 * -25.0));
	if(f_wtable<0.0){
		f_wtable = 0.0;
	}
	/* lake: added by A.Ito (2009/07/14) */
	f_wtable_lake = 0.383 * exp(0.096 * 5.0);
	
	/* Mg C ha-1 month-1 */
	(flux->soil).ch4prod_wetland_cao[grid->m] = hr_decomp * f_temp * 
					(f_wtable*grid->f_wetland + f_wtable_lake*0.2*grid->f_lake); /* 0.2: 090717 */
	/* mg CH4 m-2 month-1 */
	(flux->soil).ch4prod_wetland_cao[grid->m] *= 16.0/12.0 * 1000000000.0 / 10000.0;
	
	/* CH4 oxidation */
	if((flux->soil).ch4prod_wetland_cao[grid->m]>0.0 && loct->gppmax>0.0){
		gpp_factor = (flux->plant).gpp[grid->m]/loct->gppmax;
		gpp_factor = (gpp_factor>0.0)?gpp_factor:0.0;
		gpp_factor = (gpp_factor<1.0)?gpp_factor:1.0;
		
		(flux->soil).ch4oxy_wetland_cao[grid->m] = (flux->soil).ch4prod_wetland_cao[grid->m] 
			* (0.60 + 0.30*gpp_factor);
	}else{
		(flux->soil).ch4oxy_wetland_cao[grid->m] = 0.0;
	}
	
	/* net flux */
	(flux->soil).ch4flux_wetland_cao[grid->m] = (flux->soil).ch4prod_wetland_cao[grid->m] 
			- (flux->soil).ch4oxy_wetland_cao[grid->m];
	
	/* paddy field *********************************/
	/* eq.6 */
	/* if(grid->tmp_2m[grid->m] > 15.0 && grid->prate_sfc[grid->m] > 50.0){
		wtable = 0.0;
	}else{
		wtable = -20.0;
	}
	f_wtable = 0.383 * exp(0.096 * wtable); */
	/* revised by A.Ito (2009/07/13) */
	f_wtable = 0.383 * (grid->inundation_ssmi[grid->m]*exp(0.096 * 3.0) 
						+ (1.0-grid->inundation_ssmi[grid->m])*exp(0.096 * -50.0));	
	if(f_wtable<0.0){
		f_wtable = 0.0;
	}
	/* Mg C ha-1 month-1 */
	/* assuming low decomposition at paddy field: 0.4, 2008/06/11 */
	(flux->soil).ch4prod_paddy_cao[grid->m] = 0.4*hr_decomp * f_temp * f_wtable;  
	/* mg CH4 m-2 month-1 */
	(flux->soil).ch4prod_paddy_cao[grid->m] *= 16.0/12.0 * 1000000000.0 / 10000.0;
	(flux->soil).ch4prod_paddy_cao[grid->m] *= grid->f_paddy;
	
	/* CH4 oxidation */
	if((flux->soil).ch4prod_paddy_cao[grid->m]>0.0 && loct->gppmax>0.0){
		gpp_factor = (flux->plant).gpp[grid->m]/loct->gppmax;
		gpp_factor = (gpp_factor>0.0)?gpp_factor:0.0;
		gpp_factor = (gpp_factor<1.0)?gpp_factor:1.0;

		(flux->soil).ch4oxy_paddy_cao[grid->m] = (flux->soil).ch4prod_paddy_cao[grid->m] 
			* (0.60 + 0.30*gpp_factor);
	}else{
		(flux->soil).ch4oxy_paddy_cao[grid->m] = 0.0;
	}
	
	/* net flux */
	(flux->soil).ch4flux_paddy_cao[grid->m] = (flux->soil).ch4prod_paddy_cao[grid->m] 
			- (flux->soil).ch4oxy_paddy_cao[grid->m];
}

/* aerobic CH4 emission ************************************************/
/*
Keppler, F., Hamilton, J.T.G., Bra, M. and Rkmann, T., 2006. Methane emissions from 
terrestrial plants under aerobic conditions. Nature 439, 187-191.
*/
/*
Kirschbaum, M.U.F., Bruhn, D., Etheridge, D.M., Evans, J.R., Farquhar, G.D., Gifford, R.M., 
Paul, K.I. and Winters, A.J., 2006. A comment on the quantitative significance of aerobic 
methane releasse by plants. Functional Plant Biology 33, 521-530.
*/
void f_ch4_emit_veg(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux
){
	double femit_sun, femit_shade;
	double sunshine;
	extern double MDN[ASTEP];
	
	/* base emission rate */
	femit_sun = 374.0;		/* ng gdw-1 h-1 */
	femit_shade = 119.0;	/* ng gdw-1 h-1 */
	
	sunshine = grid->dlen[grid->m] * (1.0 - grid->tcdc_clm[grid->m]);
	
	/* MASS-based scaling up **/
	/* C3, g m-2 month-1 */
	if((echar->c3).season[grid->m]!=0){
		(flux->c3).emit_ch4_kirschbaum_mass[grid->m] = ((mass->c3).mfol[grid->m]*dmTc*100.0) * 
			(sunshine*femit_sun + (24.0 - sunshine)*femit_shade) * pow(10.0, -9.0) * MDN[grid->m];
	}else{
		(flux->c3).emit_ch4_kirschbaum_mass[grid->m] = 0.0;
	}

	/* C4, g m-2 month-1 */
	if((echar->c4).season[grid->m]!=0){
		(flux->c4).emit_ch4_kirschbaum_mass[grid->m] = ((mass->c4).mfol[grid->m]*dmTc*100.0) * 
			(sunshine*femit_sun + (24.0 - sunshine)*femit_shade) * pow(10.0, -9.0) * MDN[grid->m];
	}else{
		(flux->c4).emit_ch4_kirschbaum_mass[grid->m] = 0.0;
	}
	
	/* PHOTO-based scaling up **/
	if(sunshine > 0.0){
		(flux->c3).emit_ch4_kirschbaum_photo[grid->m] = 2.0 * (16.0/12.0) * 
			((flux->c3).npp[grid->m]*1000.0) / 30000.0 * (1.0 + (24.0 - sunshine)/sunshine * 
			femit_shade / femit_sun);
		(flux->c4).emit_ch4_kirschbaum_photo[grid->m] = 2.0 * (16.0/12.0) * 
			((flux->c4).npp[grid->m]*1000.0) / 30000.0 * (1.0 + (24.0 - sunshine)/sunshine * 
			femit_shade / femit_sun);
	}else{
		(flux->c3).emit_ch4_kirschbaum_photo[grid->m] = 0.0;
		(flux->c4).emit_ch4_kirschbaum_photo[grid->m] = 0.0;
	}
}

/* wetland CH4 emission by Walter & Heimann */
/*
	Walter, B. P., M. Heimann, et al. (2001). "Modeling modern methane emissions 
	from natural wetlands 1. Model description and results." 
	Journal of Geophysical Research 106(D24): 34189-34206.
*/
void f_ch4_emit_walter(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Flux *flux
){
	
}
