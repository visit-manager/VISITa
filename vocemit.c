/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/*                  Revised August 14, 2007					*/

#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* VOC emission model by Guenther **********************************/
/*
Guenther, A., 1997. Seasonal and spatial variations in natural volatile 
organic compound emissions. Ecological Applications 7, 34-45.
*/
void f_voc_emit_guenther97(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Flux *flux
){
	/* biogenic emission factors, micro g C/gdm/hour   */
	/* Latherie et al. (2006) */
	double emit_potent_isopr[18] = {0.0,
		24.0, 24.0, 16.0, 8.0, 45.0, 8.0, 8.0, 8.0, 
		24.0, 16.0, 16.0, 24.0, 16.0, 24.0, 16.0, 
		5.0, 5.0};
	double emit_potent_monotrp[18] = {0.0,
		0.8, 0.8, 1.2, 2.4, 0.8, 2.4, 2.4, 2.4, 
		1.2, 0.8, 0.8, 1.2, 0.8, 1.2, 0.8, 
		0.2, 0.2};
	double emit_potent_methanl[18] = {0.0,
		0.6, 0.6, 0.9, 1.8, 0.6, 1.8, 1.8, 1.8, 
		0.9, 0.6, 0.6, 0.9, 0.6, 0.9, 0.6, 
		2.0, 2.0};
	double emit_potent_acetone[18] = {0.0,
		0.29, 0.29, 0.43, 0.87, 0.29, 0.87, 0.87, 0.87, 
		0.43, 0.29, 0.29, 0.43, 0.29, 0.43, 0.29, 
		0.07, 0.07};
	double emit_potent_actaldhd[18] = {0.0,
		0.1, 0.1, 0.15, 0.3, 0.1, 0.3, 0.3, 0.3, 
		0.15, 0.1, 0.1, 0.15, 0.1, 0.15, 0.1, 
		0.025, 0.025};
	double emit_potent_frmardhd[18] = {0.0,
		0.07, 0.07, 0.1, 0.2, 0.07, 0.2, 0.2, 0.2, 
		0.1, 0.07, 0.07, 0.1, 0.07, 0.1, 0.07, 
		0.017, 0.017};
	double emit_potent_formacd[18] = {0.0,
		0.01, 0.01, 0.015, 0.03, 0.01, 0.03, 0.03, 0.03, 
		0.015, 0.01, 0.01, 0.015, 0.01, 0.015, 0.01, 
		0.0025, 0.0025};
	double emit_potent_acetacd[18] = {0.0,
		0.002, 0.002, 0.003, 0.006, 0.002, 0.006, 0.006, 0.006, 
		0.003, 0.002, 0.002, 0.003, 0.002, 0.003, 0.002, 
		0.0005, 0.0005};
	/* Tao and Jain (2005) **/
	double emit_potent_co[18] = {0.0,
		0.3, 0.3, 0.3, 0.3, 0.3, 0.36, 0.3, 0.3, 
		0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 
		0.3, 0.3};
	
	double foliar_dens, leaf_temp;
	double f_ppfd, f_temp_isopr, f_temp_monotrp, f_phenology;
	double aa, bb, cc;
	
	/* veg_sage:
	0	Ocean
	1	Tropical Evergreen Forest/Woodland
	2	Tropical Deciduous Forest/Woodland
	3	Temperate Broadleaf Evergreen Forest/Woodland
	4	Temperate Needleleaf Evergreen Forest/Woodland
	5	Temperate Deciduous Forest/Woodland
	6	Boreal Evergreen Forest/Woodland
	7	Boreal Deciduous Forest/Woodland
	8	Evergreen/Deciduous Mixed Forest/Woodland
	9	Savanna
	10	Grassland/Steppe
	11	Dense Shrubland
	12	Open Shrubland
	13	Tundra
	14	Desert
	15	Polar Desert/Rock/Ice
	*/

	/* doliar density, g d.m. C / m2   */
	foliar_dens = (mass->plant).fol * 100.0 * dmTc;
	
	/* light factor */
	f_ppfd = 0.0027 * 1.066 * grid->par[grid->m] / 
		sqrt(1.0 + 0.0027 * 0.0027 * grid->par[grid->m] * grid->par[grid->m]);
	
	/* temperature factor */
	leaf_temp = grid->tmp_sfc[grid->m] + ZAT;
	aa = exp(95000.0*(leaf_temp - 293.15)/(8.314 * leaf_temp * 293.15));
	bb = 0.961 + exp(230000.0*(leaf_temp - 314.0)/(8.314 * leaf_temp * 314.0));
	f_temp_isopr = aa/bb;
	
	f_temp_monotrp = exp(0.09*(leaf_temp - 293.15)); 
	
	/* other factor */
	f_phenology = 1.0;
	
	cc = foliar_dens * (double)grid->mm[grid->m] * grid->dlen[grid->m];
	
	/* VOC emission, micro g C m-2 month-1  */
	flux->voc_isopr_g97[grid->m] = emit_potent_isopr[grid->veg_sage] * cc * f_ppfd * f_temp_isopr;
	flux->voc_monotrp_g97[grid->m] = emit_potent_monotrp[grid->veg_sage] * cc * f_temp_monotrp;
	flux->voc_methanl_g97[grid->m] = emit_potent_methanl[grid->veg_sage] * cc * f_temp_monotrp;
	flux->voc_acetone_g97[grid->m] = emit_potent_acetone[grid->veg_sage] * cc * f_temp_monotrp;
	flux->voc_actaldhd_g97[grid->m] = emit_potent_actaldhd[grid->veg_sage] * cc * f_temp_monotrp;
	flux->voc_frmardhd_g97[grid->m] = emit_potent_frmardhd[grid->veg_sage] * cc * f_temp_monotrp;
	flux->voc_formacd_g97[grid->m] = emit_potent_formacd[grid->veg_sage] * cc * f_temp_monotrp;
	flux->voc_acetacd_g97[grid->m] = emit_potent_acetacd[grid->veg_sage] * cc * f_temp_monotrp;
	flux->voc_co_g97[grid->m] = emit_potent_co[grid->veg_sage] * cc * f_temp_monotrp;
}
