/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed  in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* [history]
	create: 2008 / 02 / 27 by Akihiko Ito
*/

/* header files */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<string.h>
#include"structure.h"
#include"prototype.h"

/* dissolved organic carbon (DOC) *************************************/
/*
Boyer, E.W., Hornberger, G.M., Bencala, K.E. and McKnight, D., 1996. 
Overview of a simple model describing variation of dissolved organic 
carbon in an upland catchment. Ecological Modelling, 86:183-188.
*/
void f_doc_boyer(
	struct Grid *grid,  
	struct Loct *loct,  
	struct Smas *mass, 
	struct Sflx *flux
){
	double aa, kk;
	double stemp, doc_est, prm_ensen;
	
	aa = 0.11;
	kk = 0.002;
	
	stemp = (grid->tmp10_soil[grid->m] + grid->tmp200_soil[grid->m])/2.0;
	
	if(loct->sww > 0.0){
		doc_est = mass->doc + aa*pow(10.0, 0.04*stemp)
			 - (1.0 - exp(-kk*stemp)) * mass->doc
			 - loct->ro2[grid->m]/loct->sww * mass->doc;
		if(doc_est <= 0.0){
			doc_est = 0.0;
		}
		if(doc_est >= 100.0){
			doc_est = 100.0;
		}
	}else{
		doc_est = 0.0;
		flux->doc_boyer[grid->m] = 0.0;
	}
    
    /* parameter ensemble: 2014/11/19  */
    prm_ensen = 1.0;
    if(PARAM_PTB == 11){
        if(PARAM_ENS == 1){
            prm_ensen *= 0.7;
        }
        if(PARAM_ENS == 2){
            prm_ensen *= 0.8;
        }
        if(PARAM_ENS == 3){
            prm_ensen *= 0.9;
        }
        if(PARAM_ENS == 4){
            prm_ensen *= 1.1;
        }
        if(PARAM_ENS == 5){
            prm_ensen *= 1.2;
        }
        if(PARAM_ENS == 6){
            prm_ensen *= 1.3;
        }
    }
	
    /* C-budget parameter ensemble: 2018/06/05  */
    if(PARAM_PTB == 20){
        /* prm_ensen = 1.0 + 0.3 * f_pert[4]; */
        prm_ensen = 1.0 + 0.15 * f_pert[4];
    }
    
	/* DOC mass, mg/L */
	mass->doc = doc_est * prm_ensen;
	mass->doc_m[grid->m] = mass->doc;
 
	/* runoff DOC, mg/L -> gC/ha/month */
	flux->doc_boyer[grid->m] = mass->doc * loct->ro2[grid->m] * 10000.0 / 1000.0;
	
	/*  flux->doc_boyer[grid->m] *= (double)(grid->mm[grid->m]);  */
}
