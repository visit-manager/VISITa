/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed  in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* file created: 01/07/2008  */

/*** ecophysiological and morphological parameters ***/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* 0: others */
/* 1: rice */
/* 2: wheat */
/* 3: maize */
/* 4: biofuel (added 2015/08/21 ) */
/* 5: agroforestry trees */
/* 6: agroforestry crops */

/* Setting characteristics for crops ***************************/
void parameterCrop(
	struct Grid *grid, 
	struct Pchar *crop
){
	short phototypev[NVEG_CROP] = {3, 3, 4, 4, 3, 3};
	double albcv[NVEG_CROP] = {0.15, 0.14, 0.14, 0.14, 0.14, 0.14};
	double alloc_assv[NVEG_CROP] = {0.21, 0.21, 0.21, 0.21, 0.21, 0.21};
	double alloc_abgv[NVEG_CROP] = {0.52, 0.51, 0.52, 0.52, 0.52, 0.52};
	double slav[NVEG_CROP] = {170.0, 160.0, 150.0, 150.0, 140.0, 150.0};
	double eK0v[NVEG_CROP] = {0.50, 0.50, 0.50, 0.50, 0.50, 0.50};
	double luev[NVEG_CROP] = {0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
	double pmaxv[NVEG_CROP] = {10.8, 10.6, 10.6, 10.6};
	double topt0v[NVEG_CROP] = {25.0, 25.0, 25.0, 25.0};
	double tminv[NVEG_CROP] = {11.0, 9.0, 11.0, 11.0};
	double tmaxv[NVEG_CROP] = {45.0, 45.0, 45.0, 45.0};
	double gs_b0v[NVEG_CROP] = {10.0, 10.0, 10.0, 10.0};
	double gs_b1v[NVEG_CROP] = {140000.0, 140000.0, 140000.0, 140000.0};
	double gs_b2v[NVEG_CROP] = {4.5, 4.5, 4.8, 4.8};
	double km_nstlv[NVEG_CROP] = {0.34, 0.35, 0.34, 0.34};
	double kmciv[NVEG_CROP] = {33.0, 33.0, 33.0, 33.0};
	double cmpcdv[NVEG_CROP] = {50.0, 50.0, 50.0, 50.0};
	double rgfv[NVEG_CROP] = { 0.50, 0.50, 0.50, 0.50};
	double rgcv[NVEG_CROP] = {0.14, 0.14, 0.14, 0.14};
	double rgrv[NVEG_CROP] = {0.22, 0.23, 0.23, 0.23};
	double rmf0v[NVEG_CROP] = {1.37, 1.38, 1.37, 1.37};
	double rmc_sv[NVEG_CROP] = {0.054, 0.054, 0.054, 0.054};
	double rmr_sv[NVEG_CROP] = {0.217, 0.227, 0.232, 0.232};
	double rmc_hv[NVEG_CROP] = {0.0026, 0.0028, 0.0029, 0.0029};
	double rmr_hv[NVEG_CROP] = {0.0161, 0.0166, 0.0166, 0.0166};
	double qTfv[NVEG_CROP] = {2.0, 2.0, 2.0, 2.0};
	double qTcv[NVEG_CROP] = {2.0, 2.0, 2.0, 2.0};
	double qTrv[NVEG_CROP] = {2.0, 2.0, 2.0, 2.0};
	double lfv[NVEG_CROP] = {0.00238, 0.00249, 0.00241, 0.00241};
	double lcv[NVEG_CROP] = {0.000060, 0.000060, 0.000060, 0.000060};
	double lrv[NVEG_CROP] = {0.00036, 0.00036, 0.00035, 0.00035};
	double dcdv[NVEG_CROP] = {0.00, 0.00, 0.50, 0.50};
	double n_salvagev[NVEG_CROP] = {0.50, 0.50, 0.50, 0.50};
	double cn0_strgv[NVEG_CROP] = {20.0, 20.0, 20.0, 20.0};
	double kn_nphotov[NVEG_CROP] = {100.0, 100.0, 100.0, 100.0};
	double amax_nphotov[NVEG_CROP] = {30.0, 30.0, 30.0, 30.0};
	double rd_nv[NVEG_CROP] = {0.05, 0.05, 0.05, 0.05};
	double crit_tempv[NVEG_CROP] = {5.0, 5.0, 5.0, 5.0};
	double crit_gddv[NVEG_CROP] = {300.0, 300.0, 300.0, 300.0};
	double root_dist_av[NVEG_CROP] = {5.558, 5.558, 5.558, 5.558};
	double root_dist_bv[NVEG_CROP] = {2.614, 2.614, 2.614, 2.614};
	double root_depthv[NVEG_CROP] = {1.5, 1.5, 1.5, 1.5};
	double n_leafv[NVEG_CROP] = {120.0, 120.0, 120.0, 120.0};
	
	crop->phototype = phototypev[grid->type_crop -1];
	crop->albedo = albcv[grid->type_crop -1];
	crop->alloc_ass = alloc_assv[grid->type_crop -1];
	crop->alloc_abg = alloc_abgv[grid->type_crop -1];
	
	crop->sla = slav[grid->type_crop -1];
	crop->eK0 = eK0v[grid->type_crop -1];
	crop->lue0 = luev[grid->type_crop -1];
	crop->pmax = pmaxv[grid->type_crop -1];
	crop->topt0 = topt0v[grid->type_crop -1];
	crop->tmin = tminv[grid->type_crop -1];
	crop->tmax = tmaxv[grid->type_crop -1];
	crop->gs_b0 = gs_b0v[grid->type_crop -1];
	crop->gs_b1 = gs_b1v[grid->type_crop -1];
	crop->gs_b2 = gs_b2v[grid->type_crop -1];
	crop->km_nstl = km_nstlv[grid->type_crop -1];
	crop->kmci = kmciv[grid->type_crop -1];
	crop->cmpcd0 = cmpcdv[grid->type_crop -1];
	
	crop->rgf = rgfv[grid->type_crop -1];
	crop->rgc = rgcv[grid->type_crop -1];
	crop->rgr = rgrv[grid->type_crop -1];
	crop->rmf0 = rmf0v[grid->type_crop -1];
	crop->rmc_s = rmc_sv[grid->type_crop -1];
	crop->rmr_s = rmr_sv[grid->type_crop -1] ;
	crop->rmc_h = rmc_hv[grid->type_crop -1];
	crop->rmr_h = rmr_hv[grid->type_crop -1] ;
	crop->qTf0 = qTfv[grid->type_crop -1];
	crop->qTc0 = qTcv[grid->type_crop -1];
	crop->qTr0 = qTrv[grid->type_crop -1];
	
	crop->lf0 = lfv[grid->type_crop -1];
	crop->lc0 = lcv[grid->type_crop -1];
	crop->lr0 = lrv[grid->type_crop -1];
	crop->dcd = dcdv[grid->type_crop -1];
	crop->n_salvage = n_salvagev[grid->type_crop -1];
	crop->cn0_strg = cn0_strgv[grid->type_crop -1];
	crop->kn_nphoto = kn_nphotov[grid->type_crop -1];
	crop->amax_nphoto = amax_nphotov[grid->type_crop -1];
	crop->rd_n = rd_nv[grid->type_crop -1];
	
	crop->crit_temp = crit_tempv[grid->type_crop -1];
	crop->crit_gdd = crit_gddv[grid->type_crop -1];

	crop->root_dist_a = root_dist_av[grid->type_crop -1];
	crop->root_dist_b = root_dist_bv[grid->type_crop -1];
	crop->root_depth = root_depthv[grid->type_crop -1];

	crop->n_leaf_df97 = n_leafv[grid->type_crop -1];
}

/* setting soil parameters ***********************/
void parameterSoil_crop(
	struct Grid *grid, 
	struct Schar *soil
){
	double albgv[NVEG_CROP]={0.13, 0.13, 0.13, 0.13};
	double rlv[NVEG_CROP]={1.40, 1.40, 1.30, 1.30};
	double rhv[NVEG_CROP]={0.105, 0.110, 0.104, 0.104};
	double qTlv[NVEG_CROP]={2.0, 2.0, 2.0, 2.0};
	double qThv[NVEG_CROP]={2.0, 2.0, 2.0, 2.0};
	double kmlv[NVEG_CROP]={0.21, 0.20, 0.21, 0.21};
	double kmslv[NVEG_CROP]={0.11, 0.11, 0.11, 0.11};
	double kmhv[NVEG_CROP]={0.23, 0.20, 0.23, 0.23};
	double kmshv[NVEG_CROP]={0.10, 0.10, 0.10, 0.10};
	double mev[NVEG_CROP]={1.38, 1.38, 1.38, 1.38};
	
	soil->albedo0 = albgv[grid->type_crop -1];
	soil->rl0=soil->rl = rlv[grid->type_crop -1];
	soil->rh0=soil->rh = rhv[grid->type_crop -1]*0.9;
	soil->qTl = qTlv[grid->type_crop -1];
	soil->qTh = qThv[grid->type_crop -1];
	soil->kml = kmlv[grid->type_crop -1];
	soil->kmh = kmhv[grid->type_crop -1];
	soil->kmsl = kmslv[grid->type_crop -1];
	soil->kmsh = kmshv[grid->type_crop -1];
	soil->me = mev[grid->type_crop -1];
}
