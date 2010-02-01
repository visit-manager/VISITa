/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* file created: 01/07/2008 by A.Ito */

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

/* Setting characteristics for crops ***************************/
void parameterCrop(
	struct Grid *grid, 
	struct Pchar *crop
){
	short phototypev[VEG_NUM_CROP] = {3, 3, 4};
	double albcv[VEG_NUM_CROP] = {0.15, 0.14, 0.14};
	double alloc_assv[VEG_NUM_CROP] = {0.21, 0.21, 0.21};
	double alloc_abgv[VEG_NUM_CROP] = {0.52, 0.51, 0.52}; 
	double slav[VEG_NUM_CROP] = {170.0, 160.0, 150.0};
	double eK0v[VEG_NUM_CROP] = {0.50, 0.50, 0.50};
	double luev[VEG_NUM_CROP] = {0.05, 0.05, 0.05};
	double pmaxv[VEG_NUM_CROP] = {10.8, 10.6, 10.6};
	double topt0v[VEG_NUM_CROP] = {25.0, 25.0, 25.0};
	double tminv[VEG_NUM_CROP] = {11.0, 9.0, 11.0}; 
	double tmaxv[VEG_NUM_CROP] = {45.0, 45.0, 45.0};
	double gs_b0v[VEG_NUM_CROP] = {10.0, 10.0, 10.0};
	double gs_b1v[VEG_NUM_CROP] = {140000.0, 140000.0, 140000.0};
	double gs_b2v[VEG_NUM_CROP] = {4.5, 4.5, 4.8};
	double km_nstlv[VEG_NUM_CROP] = {0.34, 0.35, 0.34};
	double kmciv[VEG_NUM_CROP] = {33.0, 33.0, 33.0};
	double cmpcdv[VEG_NUM_CROP] = {50.0, 50.0, 50.0};
	double rgfv[VEG_NUM_CROP] = { 0.50, 0.50, 0.50};
	double rgcv[VEG_NUM_CROP] = {0.14, 0.14, 0.14};
	double rgrv[VEG_NUM_CROP] = {0.22, 0.23, 0.23};
	double rmf0v[VEG_NUM_CROP] = {1.37, 1.38, 1.37}; 
	double rmc_sv[VEG_NUM_CROP] = {0.054, 0.054, 0.054};
	double rmr_sv[VEG_NUM_CROP] = {0.217, 0.227, 0.232};
	double rmc_hv[VEG_NUM_CROP] = {0.0026, 0.0028, 0.0029};
	double rmr_hv[VEG_NUM_CROP] = {0.0161, 0.0166, 0.0166};
	double qTfv[VEG_NUM_CROP] = {2.0, 2.0, 2.0};
	double qTcv[VEG_NUM_CROP] = {2.0, 2.0, 2.0};
	double qTrv[VEG_NUM_CROP] = {2.0, 2.0, 2.0};
	double lfv[VEG_NUM_CROP] = {0.00238, 0.00249, 0.00241};
	double lcv[VEG_NUM_CROP] = {0.000060, 0.000060, 0.000060};
	double lrv[VEG_NUM_CROP] = {0.00036, 0.00036, 0.00035};
	double dcdv[VEG_NUM_CROP] = {0.00, 0.00, 0.50};
	double n_salvagev[VEG_NUM_CROP] = {0.50, 0.50, 0.50};
	double cn0_strgv[VEG_NUM_CROP] = {20.0, 20.0, 20.0};
	double kn_nphotov[VEG_NUM_CROP] = {100.0, 100.0, 100.0};
	double amax_nphotov[VEG_NUM_CROP] = {30.0, 30.0, 30.0};
	double rd_nv[VEG_NUM_CROP] = {0.05, 0.05, 0.05};
	double crit_tempv[VEG_NUM_CROP] = {5.0, 5.0, 5.0};
	double crit_gddv[VEG_NUM_CROP] = {300.0, 300.0, 300.0};
	double root_dist_av[VEG_NUM_CROP] = {5.558, 5.558, 5.558};
	double root_dist_bv[VEG_NUM_CROP] = {2.614, 2.614, 2.614};
	double root_depthv[VEG_NUM_CROP] = {1.5, 1.5, 1.5};
	double n_leafv[VEG_NUM_CROP] = {120.0, 120.0, 120.0};
	
	crop->phototype = phototypev[grid->veg_crop];
	crop->albedo = albcv[grid->veg_crop];
	crop->alloc_ass = alloc_assv[grid->veg_crop]; 
	crop->alloc_abg = alloc_abgv[grid->veg_crop]; 
	
	crop->sla = slav[grid->veg_crop];
	crop->eK0 = eK0v[grid->veg_crop];
	crop->lue0 = luev[grid->veg_crop]; 
	crop->pmax = pmaxv[grid->veg_crop]; 
	crop->topt0 = topt0v[grid->veg_crop]; 
	crop->tmin = tminv[grid->veg_crop]; 
	crop->tmax = tmaxv[grid->veg_crop];
	crop->gs_b0 = gs_b0v[grid->veg_crop]; 
	crop->gs_b1 = gs_b1v[grid->veg_crop]; 
	crop->gs_b2 = gs_b2v[grid->veg_crop]; 
	crop->km_nstl = km_nstlv[grid->veg_crop]; 
	crop->kmci = kmciv[grid->veg_crop]; 
	crop->cmpcd0 = cmpcdv[grid->veg_crop]; 
	
	crop->rgf = rgfv[grid->veg_crop]; 
	crop->rgc = rgcv[grid->veg_crop]; 
	crop->rgr = rgrv[grid->veg_crop];
	crop->rmf0 = rmf0v[grid->veg_crop]; 
	crop->rmc_s = rmc_sv[grid->veg_crop]; 
	crop->rmr_s = rmr_sv[grid->veg_crop] ;
	crop->rmc_h = rmc_hv[grid->veg_crop]; 
	crop->rmr_h = rmr_hv[grid->veg_crop] ;
	crop->qTf0 = qTfv[grid->veg_crop]; 
	crop->qTc0 = qTcv[grid->veg_crop]; 
	crop->qTr0 = qTrv[grid->veg_crop];
	
	crop->lf0 = lfv[grid->veg_crop]; 
	crop->lc0 = lcv[grid->veg_crop];
	crop->lr0 = lrv[grid->veg_crop];
	crop->dcd = dcdv[grid->veg_crop];
	crop->n_salvage = n_salvagev[grid->veg_crop];
	crop->cn0_strg = cn0_strgv[grid->veg_crop];
	crop->kn_nphoto = kn_nphotov[grid->veg_crop];
	crop->amax_nphoto = amax_nphotov[grid->veg_crop];
	crop->rd_n = rd_nv[grid->veg_crop];
	
	crop->crit_temp = crit_tempv[grid->veg_crop];
	crop->crit_gdd = crit_gddv[grid->veg_crop];

	crop->root_dist_a = root_dist_av[grid->veg_crop];
	crop->root_dist_b = root_dist_bv[grid->veg_crop];
	crop->root_depth = root_depthv[grid->veg_crop];

	crop->n_leaf_df97 = n_leafv[grid->veg_crop];
}

/* setting soil parameters ***********************/
void parameterSoil_crop(
	struct Grid *grid, 
	struct Schar *soil
){
	double albgv[VEG_NUM_CROP]={0.13, 0.13, 0.13};
	double rlv[VEG_NUM_CROP]={1.40, 1.40, 1.30};
	double rhv[VEG_NUM_CROP]={0.105, 0.110, 0.104};
	double qTlv[VEG_NUM_CROP]={2.0, 2.0, 2.0};
	double qThv[VEG_NUM_CROP]={2.0, 2.0, 2.0};
	double kmlv[VEG_NUM_CROP]={0.21, 0.20, 0.21};
	double kmslv[VEG_NUM_CROP]={0.11, 0.11, 0.11};
	double kmhv[VEG_NUM_CROP]={0.23, 0.20, 0.23};
	double kmshv[VEG_NUM_CROP]={0.10, 0.10, 0.10};
	double mev[VEG_NUM_CROP]={1.38, 1.38, 1.38};
	
	soil->albedo0 = albgv[grid->veg_crop]; 	 				 
	soil->rl0=soil->rl = rlv[grid->veg_crop]; 
	soil->rh0=soil->rh = rhv[grid->veg_crop]*0.9;
	soil->qTl = qTlv[grid->veg_crop]; 
	soil->qTh = qThv[grid->veg_crop];
	soil->kml = kmlv[grid->veg_crop]; 
	soil->kmh = kmhv[grid->veg_crop];
	soil->kmsl = kmslv[grid->veg_crop]; 
	soil->kmsh = kmshv[grid->veg_crop];
	soil->me = mev[grid->veg_crop];	
}
