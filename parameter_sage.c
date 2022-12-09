/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed  in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/***  ecophysiological  and  morphological  parameters  ***/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* Setting characteristics for C3 plants ***************************/
void  parameterC3_sage(
	struct Grid *grid, 
	struct Pchar *c3
){
	double  albcv[NVEG_SAGE]={ 
					0.15,  0.14,  0.14,  0.14,  0.15,  0.15,  0.14,  0.14, 
					0.14,  0.14,  0.14,  0.14,  0.14,  0.20,  0.20,  0.19};
	double  alloc_assv[NVEG_SAGE]={
					0.21,  0.21,  0.21,  0.20,  0.20,  0.22,  0.21,  0.19,  
					0.19,  0.19,  0.19,  0.19,  0.23,  0.23,  0.25,  0.25};
	double  alloc_abgv[NVEG_SAGE]={
					0.52,  0.51,  0.52,  0.51,  0.52,  0.50,  0.52,  0.52,  
					0.50,  0.49,  0.50,  0.49,  0.45,  0.30,  0.13,  0.13}; 
	double slav[NVEG_SAGE]={
					170.0,  160.0,  150.0,  160.0,  160.0,  140.0,  130.0,  125.0, 
					115.0,  115.0,  110.0,  110.0,  130.0,  135.0,  135.0,  135.0};
	double eK0v[NVEG_SAGE]={
					0.50,  0.50,  0.50,  0.52,  0.52,  0.52,  0.53,  0.54, 
					0.54,  0.54,  0.54,  0.54,  0.50,  0.42,  0.43,  0.42};
	double luev[NVEG_SAGE]={
					0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
					0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
	double pmaxv[NVEG_SAGE]={
					10.8, 10.6, 10.6,  9.7,  9.7,  10.5,   9.7,   8.8, 
					  8.8,  10.2,   8.8,  10.2,   8.8,   8.8,   9.7,   9.7};
	double topt0v[NVEG_SAGE]={
					25.0, 25.0, 25.0, 20.0, 20.0, 24.0, 21.0, 18.0, 
					18.0, 18.0, 18.0, 18.0, 22.0, 22.0, 23.0, 23.0};
	double tminv[NVEG_SAGE]={
					11.0,  9.0,  11.0,  7.0,  7.0,  7.0,  1.0,  2.0, 
					0.0,  -2.0,  0.5,  -2.0,  5.0,  5.0,  5.0,  8.0};  
	double tmaxv[NVEG_SAGE]={
					45.0, 45.0, 45.0, 42.0, 42.0, 44.0, 42.0, 40.0, 
					40.0, 40.0, 40.0, 40.0, 40.0, 40.0, 43.0, 45.0};
	double gs_b0v[NVEG_SAGE]={
					10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 
					10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
	double gs_b1v[NVEG_SAGE]={
					140000.0, 140000.0, 140000.0, 140000.0, 140000.0, 135000.0, 121000.0, 110000.0, 
					110000.0, 110000.0, 110000.0, 110000.0, 130000.0, 130000.0, 120000.0, 130000.0};
	double gs_b2v[NVEG_SAGE]={
					4.5, 4.5, 4.8, 4.7, 4.7, 5.3, 4.5, 4.7, 
					4.7, 4.9, 4.7, 4.9, 4.7, 4.2, 4.2, 4.2};
	double km_nstlv[NVEG_SAGE]={
					0.34, 0.35, 0.34, 0.34, 0.34, 0.31, 0.31, 0.34, 
					0.34, 0.32, 0.34, 0.32, 0.32, 0.32, 0.32, 0.32};
	double kmciv[NVEG_SAGE]={
					33.0, 33.0, 33.0, 30.0, 30.0, 30.0, 30.0, 30.0, 
					30.0, 30.0, 30.0, 30.0, 30.0, 40.0, 35.0, 40.0};
	double cmpcdv[NVEG_SAGE]={
					50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 
					50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0, 50.0};
	double rgfv[NVEG_SAGE]={ 
					0.50, 0.50, 0.50, 0.50, 0.50, 0.49, 0.495, 0.49, 
					0.49, 0.51, 0.49, 0.51, 0.52, 0.52, 0.51, 0.52};
	double rgcv[NVEG_SAGE]={
					0.14, 0.14, 0.14, 0.14, 0.14, 0.14, 0.142, 0.14, 
					0.14, 0.14, 0.14, 0.14, 0.21, 0.23, 0.21, 0.23};
	double rgrv[NVEG_SAGE]={
					0.22, 0.23, 0.23, 0.22, 0.22, 0.22, 0.232, 0.23, 
					0.23, 0.23, 0.23, 0.23, 0.31, 0.32, 0.29, 0.32};
	double rmf0v[NVEG_SAGE]={
					1.37, 1.38, 1.37, 1.37, 1.37, 1.27, 1.26, 1.25, 
					1.25, 1.27, 1.25, 1.27, 1.48, 1.54, 1.49, 1.51};   
	double rmc_sv[NVEG_SAGE]={
					0.054, 0.054, 0.054, 0.056, 0.056, 0.055, 0.055, 0.056, 
					0.056, 0.056, 0.056, 0.056, 0.070, 0.080, 0.076, 0.080};
	double rmr_sv[NVEG_SAGE]={
					0.217, 0.227, 0.232, 0.226, 0.226, 0.219, 0.227, 0.232, 
					0.232, 0.232, 0.232, 0.232, 0.318, 0.478, 0.394, 0.408};
	double rmc_hv[NVEG_SAGE]={
					0.0026, 0.0028, 0.0029, 0.0029, 0.0029, 0.0030, 0.0029, 0.0035, 
					0.0035, 0.0035, 0.0037, 0.0037, 0.0064, 0.0074, 0.0069, 0.0084};
	double rmr_hv[NVEG_SAGE]={
					0.0161, 0.0166, 0.0166, 0.0166, 0.0166, 0.0162, 0.0168, 0.0166, 
					0.0176, 0.0179, 0.0181, 0.0181, 0.0291, 0.0431, 0.0351, 0.0381};
	double qTfv[NVEG_SAGE]={
					2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 
					2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
	double qTcv[NVEG_SAGE]={
					2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 
					2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
	double qTrv[NVEG_SAGE]={
					2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 
					2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
	double lfv[NVEG_SAGE]={
					0.00238, 0.00249, 0.00241, 0.00208, 0.00193, 0.00130, 0.00137, 0.00133, 
					0.00133, 0.00120, 0.00133, 0.00120, 0.00231, 0.00211, 0.00211, 0.00251};
	double lcv[NVEG_SAGE]={
					0.000060, 0.000060, 0.000060, 0.000062, 0.000063, 0.000062, 0.000062, 0.000064, 
					0.000064, 0.000060, 0.000064, 0.000062, 0.000085, 0.000102, 0.000124, 0.000101};
	double lrv[NVEG_SAGE]={
					0.00036, 0.00036, 0.00035, 0.00036, 0.00037, 0.00037, 0.00037, 0.00038, 
					0.00038, 0.00038, 0.00038, 0.00038, 0.00080, 0.00100, 0.00107, 0.00085};
	double dcdv[NVEG_SAGE]={
					0.00, 0.00, 0.50, 0.40, 0.00, 0.00, 0.00, 0.00, 
					0.00, 0.55, 0.00, 0.55, 0.00, 0.60, 0.60, 0.50};
	double n_salvagev[NVEG_SAGE]={
					0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 
					0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50};
	double cn0_strgv[NVEG_SAGE]={
					20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 
					20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0};
	double kn_nphotov[NVEG_SAGE]={
					100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 
					100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0};
	double amax_nphotov[NVEG_SAGE]={
					30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 
					30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0};
	double rd_nv[NVEG_SAGE]={
					0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
					0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
	double crit_tempv[NVEG_SAGE]={
					5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 
					5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
	double crit_gddv[NVEG_SAGE]={
					300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 
					300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0};
	double root_dist_av[NVEG_SAGE]={
					0.95, 7.344, 5.990, 7.344, 7.344, 5.990, 6.706, 7.066, 
					4.453, 8.235, 10.74, 6.326, 6.326, 8.992, 4.372, 4.372};
	double root_dist_bv[NVEG_SAGE]={
					0.95, 1.303, 1.955, 1.303, 1.303, 1.955, 2.175, 1.953, 
					1.631, 1.627, 2.608, 1.567, 1.567, 8.992, 0.978, 0.978};
	double root_depthv[NVEG_SAGE]={
					0.95, 3.0, 2.0, 3.0, 3.0, 2.0, 1.8, 2.0, 
					2.4, 2.4, 1.5, 2.5, 2.5, 0.5, 4.0, 4.0};
	double n_leafv[NVEG_SAGE]={
					0.0, 90.0, 90.0, 82.9, 44.1, 44.1, 44.1, 44.1, 
					60.0, 116.5, 116.5, 116.5, 116.5, 116.5, 264.7, 4.0};
	
	c3->phototype = 3;
	c3->albedo = albcv[grid->veg_sage];
	c3->alloc_ass = alloc_assv[grid->veg_sage];  
	c3->alloc_abg = alloc_abgv[grid->veg_sage];  
	
	c3->sla = slav[grid->veg_sage];
	c3->eK0 = eK0v[grid->veg_sage];
	c3->lue0 = luev[grid->veg_sage]; 
	c3->pmax = pmaxv[grid->veg_sage]; 
	c3->topt0 = topt0v[grid->veg_sage];  
	c3->tmin = tminv[grid->veg_sage];  
	c3->tmax = tmaxv[grid->veg_sage];
	c3->gs_b0 = gs_b0v[grid->veg_sage];   
	c3->gs_b1 = gs_b1v[grid->veg_sage];   
	c3->gs_b2 = gs_b2v[grid->veg_sage];   
	c3->km_nstl = km_nstlv[grid->veg_sage];   
	c3->kmci = kmciv[grid->veg_sage];  
	c3->cmpcd0 = cmpcdv[grid->veg_sage];  
	
	c3->rgf = rgfv[grid->veg_sage]; 
	c3->rgc = rgcv[grid->veg_sage]; 
	c3->rgr = rgrv[grid->veg_sage];
	c3->rmf0 = rmf0v[grid->veg_sage]; 
	c3->rmc_s = rmc_sv[grid->veg_sage]; 
	c3->rmr_s = rmr_sv[grid->veg_sage] ;
	c3->rmc_h = rmc_hv[grid->veg_sage]; 
	c3->rmr_h = rmr_hv[grid->veg_sage] ;
	c3->qTf0 = qTfv[grid->veg_sage];  
	c3->qTc0 = qTcv[grid->veg_sage];  
	c3->qTr0 = qTrv[grid->veg_sage];
	
	c3->lf0 = lfv[grid->veg_sage];  
	c3->lc0 = lcv[grid->veg_sage];
	c3->lr0 = lrv[grid->veg_sage];
	c3->dcd = dcdv[grid->veg_sage];
	c3->n_salvage = n_salvagev[grid->veg_sage];
	c3->cn0_strg = cn0_strgv[grid->veg_sage];
	c3->kn_nphoto = kn_nphotov[grid->veg_sage];
	c3->amax_nphoto = amax_nphotov[grid->veg_sage];
	c3->rd_n = rd_nv[grid->veg_sage];

	c3->crit_temp = crit_tempv[grid->veg_sage];
	c3->crit_gdd = crit_gddv[grid->veg_sage];

	c3->root_dist_a = root_dist_av[grid->veg_sage];
	c3->root_dist_b = root_dist_bv[grid->veg_sage];
	c3->root_depth = root_depthv[grid->veg_sage];

	c3->n_leaf_df97 = n_leafv[grid->veg_sage];
}

/* Setting characteristics for C4 plants ***************************/
void parameterC4_sage(
	struct Grid *grid,
	struct Pchar *c4
){
	double  albcv[NVEG_SAGE]={
					0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  
					0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  0.20,  0.20};
	double  alloc_assv[NVEG_SAGE]={
					0.20,  0.20,  0.20,  0.20,  0.20,  0.20,   0.20,   0.20,  
					0.20,  0.20,  0.20,  0.20,  0.20,  0.30,   0.30,   0.30};
	double  alloc_abgv[NVEG_SAGE]={
					0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  0.10,  
					0.10,  0.10,  0.10,  0.10,  0.10,  0.09,  0.07,  0.07}; 
	double slav[NVEG_SAGE]={
					150.0,  150.0,  150.0,  150.0,  150.0,  150.0,  150.0,  150.0,
					150.0,  150.0,  150.0,  150.0,  150.0,  120.0,  120.0,  120.0};
	double eK0v[NVEG_SAGE]={
					0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,  0.50,
					0.50,  0.50,  0.50,  0.50,  0.50,  0.40,  0.40,  0.40};
	double luev[NVEG_SAGE]={
					0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05,  
					0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
	double pmaxv[NVEG_SAGE]={
					13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 13.0, 
					13.0, 13.0, 13.0, 13.0, 13.0,   9.3,   8.9,   8.9};
	double topt0v[NVEG_SAGE]={
					20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 
					20.0, 20.0, 20.0, 20.0, 20.0, 30.0, 30.0, 30.0};
	double tminv[NVEG_SAGE]={
					11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 11.0,
					11.0, 11.0, 11.0, 11.0, 11.0, 11.0, 10.0, 11.0};  
	double tmaxv[NVEG_SAGE]={
					45.0, 45.0, 45.0, 45.0, 45.0, 45.0, 45.0, 45.0, 
					45.0, 45.0, 45.0, 45.0, 45.0, 45.0, 45.0, 45.0};
	double gs_b0v[NVEG_SAGE]={
					10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 
					10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
	double gs_b1v[NVEG_SAGE]={
					130000.0, 130000.0, 130000.0, 130000.0, 130000.0, 130000.0, 130000.0, 130000.0, 
					130000.0, 130000.0, 130000.0, 130000.0, 130000.0, 130000.0, 130000.0, 130000.0};
	double gs_b2v[NVEG_SAGE]={
					5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 
					5.0, 5.0, 5.0, 5.0, 5.0, 6.9, 6.9, 6.9};
	double km_nstlv[NVEG_SAGE]={
					0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 0.30, 
					0.30, 0.30, 0.30, 0.30, 0.30, 0.25, 0.25, 0.25};
	double kmciv[NVEG_SAGE]={
					10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 
					10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0, 10.0};
	double cmpcdv[NVEG_SAGE]={
					5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 
					5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
	double rgfv[NVEG_SAGE]={
					0.53, 0.53, 0.53, 0.53, 0.53, 0.53, 0.53, 0.53,
					0.53, 0.53, 0.53, 0.53, 0.53, 0.53, 0.53, 0.53};
	double rgcv[NVEG_SAGE]={
					0.32, 0.32, 0.33, 0.32, 0.32, 0.32, 0.32, 0.32, 
					0.32, 0.32, 0.33, 0.32, 0.32, 0.31, 0.31, 0.31};
	double rgrv[NVEG_SAGE]={
					0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 
					0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41, 0.41};
	double rmf0v[NVEG_SAGE]={
					1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 1.57, 
					1.57, 1.57, 1.57, 1.57, 1.57, 1.65, 1.65, 1.65};   
	double rmc_sv[NVEG_SAGE]={
					0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 
					0.08, 0.08, 0.08, 0.08, 0.08, 0.30, 0.30, 0.30};
	double rmr_sv[NVEG_SAGE]={
					0.48, 0.48, 0.48, 0.48, 0.48, 0.48, 0.48, 0.48, 
					0.48, 0.48, 0.48, 0.48, 0.48, 0.68, 0.68, 0.68};
	double rmc_hv[NVEG_SAGE]={
					0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 0.010, 
					0.010, 0.010, 0.010, 0.010, 0.010, 0.070, 0.070, 0.070};
	double rmr_hv[NVEG_SAGE]={
					0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 0.050, 
					0.050, 0.050, 0.050, 0.050, 0.050, 0.170, 0.170, 0.170};
	double qTfv[NVEG_SAGE]={
					2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 
					2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
	double qTcv[NVEG_SAGE]={
					2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 
					2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
	double qTrv[NVEG_SAGE]={
					2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 
					2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
	double lfv[NVEG_SAGE]={
					0.0015, 0.0015, 0.0015, 0.0015, 0.0015, 0.0015, 0.0015, 0.0015, 
					0.0015, 0.0015, 0.0015, 0.0015, 0.0015, 0.0040, 0.0040, 0.0040};
	double lcv[NVEG_SAGE]={
					0.00010, 0.00010, 0.00010, 0.00010, 0.00010, 0.00010, 0.00010, 0.00010, 
					0.00010, 0.00010, 0.00010, 0.00010, 0.00010, 0.00055, 0.00055, 0.00055};
	double lrv[NVEG_SAGE]={
					0.0010, 0.0010, 0.0010, 0.0010, 0.0010, 0.0010, 0.0010, 0.0010, 
					0.0010, 0.0010, 0.0010, 0.0010, 0.0010, 0.0025, 0.0025, 0.0026};
	double dcdv[NVEG_SAGE]={
					0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
					0.0, 0.0, 0.0, 0.0, 0.0, 0.85, 0.85, 0.85};
	double n_salvagev[NVEG_SAGE]={
					0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 
					0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50, 0.50};
	double cn0_strgv[NVEG_SAGE]={
					20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 
					20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0, 20.0};
	double kn_nphotov[NVEG_SAGE]={
					100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 
					100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0, 100.0};
	double amax_nphotov[NVEG_SAGE]={
					30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 
					30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0, 30.0};
	double rd_nv[NVEG_SAGE]={
					0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 
					0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05, 0.05};
	double crit_tempv[NVEG_SAGE]={
					5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 
					5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0, 5.0};
	double crit_gddv[NVEG_SAGE]={
					300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 
					300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0, 300.0};
	double root_dist_av[NVEG_SAGE]={
					0.95, 7.344, 5.990, 7.344, 7.344, 5.990, 6.706, 7.066, 
					4.453, 8.235, 10.74, 6.326, 6.326, 8.992, 4.372, 4.372};
	double root_dist_bv[NVEG_SAGE]={
					0.95, 1.303, 1.955, 1.303, 1.303, 1.955, 2.175, 1.953, 
					1.631, 1.627, 2.608, 1.567, 1.567, 8.992, 0.978, 0.978};
	double root_depthv[NVEG_SAGE]={
					0.95, 3.0, 2.0, 3.0, 3.0, 2.0, 1.8, 2.0, 
					2.4, 2.4, 1.5, 2.5, 2.5, 0.5, 4.0, 4.0};
	double n_leafv[NVEG_SAGE]={
					0.0, 90.0, 90.0, 82.9, 44.1, 44.1, 44.1, 44.1, 
					60.0, 116.5, 116.5, 116.5, 116.5, 116.5, 264.7, 4.0};

	c4->phototype = 4;
	c4->albedo = albcv[grid->veg_sage];
	c4->alloc_ass = alloc_assv[grid->veg_sage];  
	c4->alloc_abg = alloc_abgv[grid->veg_sage];  
	
	c4->sla = slav[grid->veg_sage];
	c4->eK0 = eK0v[grid->veg_sage];
	c4->lue0 = luev[grid->veg_sage]; 
	c4->pmax = pmaxv[grid->veg_sage];
	c4->topt0 = topt0v[grid->veg_sage];  
	c4->tmin = tminv[grid->veg_sage];  
	c4->tmax = tmaxv[grid->veg_sage];
	c4->gs_b0 = gs_b0v[grid->veg_sage];   
	c4->gs_b1 = gs_b1v[grid->veg_sage];   
	c4->gs_b2 = gs_b2v[grid->veg_sage];   
	c4->km_nstl = km_nstlv[grid->veg_sage];   
	c4->kmci = kmciv[grid->veg_sage];  
	c4->cmpcd0 = cmpcdv[grid->veg_sage];  
	
	c4->rgf = rgfv[grid->veg_sage]; 
	c4->rgc = rgcv[grid->veg_sage]; 
	c4->rgr = rgrv[grid->veg_sage];
	c4->rmf0 = rmf0v[grid->veg_sage]; 
	c4->rmc_s = rmc_sv[grid->veg_sage]; 
	c4->rmr_s = rmr_sv[grid->veg_sage] ;
	c4->rmc_h = rmc_hv[grid->veg_sage]; 
	c4->rmr_h = rmr_hv[grid->veg_sage] ;
	c4->qTf0 = qTfv[grid->veg_sage];  
	c4->qTc0 = qTcv[grid->veg_sage];  
	c4->qTr0 = qTrv[grid->veg_sage];
	
	c4->lf0 = lfv[grid->veg_sage];  
	c4->lc0 = lcv[grid->veg_sage];  
	c4->lr0 = lrv[grid->veg_sage];
	c4->dcd = dcdv[grid->veg_sage];
	c4->n_salvage = n_salvagev[grid->veg_sage];
	c4->cn0_strg = cn0_strgv[grid->veg_sage];
	c4->kn_nphoto = kn_nphotov[grid->veg_sage];
	c4->amax_nphoto = amax_nphotov[grid->veg_sage];
	c4->rd_n = rd_nv[grid->veg_sage];

	c4->crit_temp = crit_tempv[grid->veg_sage];
	c4->crit_gdd = crit_gddv[grid->veg_sage];

	c4->root_dist_a = root_dist_av[grid->veg_sage];
	c4->root_dist_b = root_dist_bv[grid->veg_sage];
	c4->root_depth = root_depthv[grid->veg_sage];

	c4->n_leaf_df97 = n_leafv[grid->veg_sage];
}

/* setting soil parameters ******************************/
void  parameterSoil_sage(
	struct Grid *grid, 
	struct Schar *soil
){
	double albgv[NVEG_SAGE]={
					0.13, 0.13, 0.13, 0.13, 0.13, 0.13, 0.13, 0.13, 
					0.13, 0.13, 0.13, 0.13, 0.18, 0.18, 0.18, 0.18};
	double rlv[NVEG_SAGE]={
					1.40, 1.40, 1.30, 1.40, 1.40, 1.35, 1.40, 1.40, 
					1.40, 1.50, 1.40, 1.50, 1.32, 1.32, 1.27, 1.27};
	double rhv[NVEG_SAGE]={
					0.105, 0.110, 0.104, 0.120, 0.130, 0.100, 0.145, 0.170, 
					0.155, 0.165, 0.155, 0.165, 0.094, 0.090, 0.080, 0.080};
	double qTlv[NVEG_SAGE]={
					2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
					2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
	double qThv[NVEG_SAGE]={
					2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0,
					2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0, 2.0};
	double kmlv[NVEG_SAGE]={
					0.21, 0.20, 0.21, 0.20, 0.20, 0.20, 0.20, 0.20, 
					0.22, 0.21, 0.21, 0.21, 0.22, 0.20, 0.23, 0.23};
	double kmslv[NVEG_SAGE]={
					0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.11, 0.10, 
					0.10, 0.10, 0.10, 0.10, 0.11, 0.11, 0.11, 0.11};
	double kmhv[NVEG_SAGE]={
					0.23, 0.20, 0.23, 0.22, 0.22, 0.20, 0.22, 0.22, 
					0.21, 0.15, 0.21, 0.15, 0.21, 0.20, 0.20, 0.20};
	double kmshv[NVEG_SAGE]={
					0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.10, 0.08, 
					0.08, 0.07, 0.08, 0.07, 0.10, 0.10, 0.10, 0.10};
	double mev[NVEG_SAGE]={
					1.38, 1.38, 1.38, 1.38, 1.38, 1.38, 1.38, 1.38, 
					1.28, 1.28, 1.28, 1.23, 1.47, 1.47, 1.47, 1.47};
	
	soil->albedo0 = albgv[grid->veg_sage]; 	    				     
	soil->rl0 = soil->rl = rlv[grid->veg_sage]; 
	soil->rh0 = soil->rh = rhv[grid->veg_sage];
	soil->qTl = qTlv[grid->veg_sage];  
	soil->qTh = qThv[grid->veg_sage];
	soil->kml = kmlv[grid->veg_sage];  
	soil->kmh = kmhv[grid->veg_sage];
	soil->kmsl = kmslv[grid->veg_sage];  
	soil->kmsh = kmshv[grid->veg_sage];
	soil->me = mev[grid->veg_sage];	
}
