/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* annual NPP estimated with empirical models *************************/
void npp_empirical(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Flux *flux
){
	long f, n;
	double lhvp, aet_ann, pet_ann, rn_ann, npp_tem, npp_pre;
	double gdd, wsi, par, rdi;
    double wi, mat, lll;
	extern double MDN[ASTEP];
	
	/* annual climatology *********************************/
	aet_ann = pet_ann = rn_ann = loct->pet_prty_ann = 0.0;
	gdd = wsi = par = 0.0;
	n = 0;
	for(f=0;f<ASTEP;f++){
		/* annual AET or PET */
		aet_ann += loct->evpr[f]+loct->incep[f]+loct->trspr[f]; 
		pet_ann += loct->pm_evp[f]+loct->pm_incep[f]+loct->pm_trn[f];
		
		/* annual mean net radiation, W m-2 */
		loct->rad_net[f] = loct->rad_net_short[f] * 0.636619 * grid->dlen[f] / 24.0 - loct->rad_net_long[f];
		rn_ann += loct->rad_net[f] * MDN[grid->m]/365.0;
		
		/* PRIESTRIE-TAYLOR PET model, mm ***/
		lhvp = 1000000.0* (2.501 - 0.012/5.0*grid->tmp_2m[f]);
		loct->pet_prty[f] = 1.26 * 0.667/(loct->slope_vps[f] + 0.667) * 
				loct->rad_net[f]/lhvp*24.0*3600.0 * MDN[grid->m];
		loct->pet_prty_ann += loct->pet_prty[f];
		
		if(grid->tmp_2m[f] >= 0.0){
			n++;
			par += grid->par[f] * 0.636619 * grid->dlen[f] / 24.0;
			gdd += (grid->tmp_2m[f] - 5.0)* MDN[grid->m];
		}
	}
	if(n>0){
		par *= 1.0/(double)n;
	}
	
	if(pet_ann > 0.0){
		wsi = aet_ann / pet_ann;
	}else{
		wsi = 0.0;
	}
    
    if(grid->prate_sfc_ann > 0.0){
        rdi = (rn_ann * 24.0*3600.0*365.0) / 2500000.0 / grid->prate_sfc_ann;
    }else{
        rdi = 10.0;
    }
    
    /* Chikugo **************/
    flux->npp_chikugo = cTdm * 0.29 * (exp(-0.216*rdi*rdi)) * (rn_ann*24.0*3600.0*365.0 / pow(10.0, 7.0) / 4.1868);
	
	/* Lieth, H., 1975. Modeling the primary productivity of the world. 
	In: H. Lieth and R.H. Whittaker (Editor), Primary productivity of the biosphere. 
	Springer-Verlag, pp. 237-263.
	*/
	/* MIAMI model *************/
	npp_tem = cTdm*30.0/(1.0 + exp(1.315 - 0.119 * grid->tmp_sfc_am));
	npp_pre = cTdm*30.0*(1.0 - exp(-0.000664 * grid->prate_sfc_ann));
	flux->npp_miami = (npp_tem<npp_pre)?npp_tem:npp_pre;
	/*** MONTREAL model ***/
	flux->npp_montreal = cTdm*30.0*(1.0 - exp(-0.0009695*(aet_ann - 20.0)));
	
	/* SCHUUR NPP model ***********/
	/* Schuur, E.A.G., 2003. Productivity and global climate revisited; 
	the sensitivity of tropical forest growth to precipitation. 
	Ecology, 84:1165-1170.
	*/
	npp_tem = 17.6243/(1.0 + exp(1.3496 - grid->tmp_sfc_am * 0.071514));
	npp_pre = 0.005212*pow(grid->prate_sfc_ann, 1.12363)/exp(0.000459532 * grid->prate_sfc_ann);
	flux->npp_schuur = (npp_tem<npp_pre)?npp_tem:npp_pre;
	
	/* ROSENZWEIG model ***********/
	/* Rosenzweig, M., 1968. Net primary productivity of terrestrial environments: 
	predictions from climatological data. American Naturalist, 102:67-74.
	*/
	flux->npp_rosenzweig = cTdm*0.219*pow(aet_ann, 1.66);
	
	/* Madison model */
	/* Zaks, D. P. M., et al. (2007), From Miami to Madison: Investigating the 
	relationship between climate and terrestrial net primary production, 
	Global Biogeochemical Cycles, 21(GB3004), 10.1029/2006GB002705. */
	
	flux->npp_madison_parwsi = (0.5*par + 0.6*wsi - 0.5) * 0.01;
	if(flux->npp_madison_parwsi < 0.0){
		flux->npp_madison_parwsi = 0.0;
	}
	
    if(gdd>0.0 && gdd<10.0){
        flux->npp_madison_gddswsi = (3.96 / (1.0 + exp(6.33 - 1.5*gdd))) * (39.58 * wsi - 14.52) * 0.01;
    }else{
        flux->npp_madison_gddswsi = 3.96 * (39.58 * wsi - 14.52) * 0.01;
    }
	if(flux->npp_madison_gddswsi < 0.0){
		flux->npp_madison_gddswsi = 0.0;
	}
	
	flux->npp_madison_tp = 20.13/(1.0 + exp(9.5 - 2.25*grid->tmp_sfc_am)) 
						* 45.83 * (1.0 - exp(-3.5 * grid->prate_sfc_ann/365.0)) * 0.01;
	
	/* NCEAS model ****/
	/* Del Grosso, S., W. Parton, T. Stohlgren, D. Zhang, D. Bachelet, S. Prince, 
	K. Hibbard, and R. Olson. 2008. Global potential net primary production predicted 
	from vegetation class, precipitation, and temperature. Ecology 89:2117-2126.
	*/
	switch(grid->veg_sage){
		case 9: case 10: case 11: case 12: case 13: case 14: case 15: 
			if(grid->prate_sfc_ann>=0.0){
				flux->npp_nceas = 61.160 * (1.0 - exp(-6.05*0.00001*grid->prate_sfc_ann));
			}else{
				flux->npp_nceas = 0.0;
			}
			break;
		case 1: case 2: case 3: case 4: case 5: case 6: case 7: case 8:
			npp_tem = 25.4 / (1.0 + exp(1.584 - 0.0622*grid->tmp_sfc_am));
			npp_pre = 0.551 * pow(grid->prate_sfc_ann, 1.055) / exp(0.000306*grid->prate_sfc_ann)/100.0;
			flux->npp_nceas = (npp_tem<npp_pre)?npp_tem:npp_pre;
			break;
	}
    
    /* estimated max LAI: 2014/05/20 by A.Ito */
    /* Iio, A., K. Hikosaka, N. P. R. Anten, Y. Nakagawa, and A. Ito. 2014. 
    Global dependence of field-observed leaf area index on climate in woody 
    species: Systematic review. Global Ecology and Biogeography 3:274–285. */
    
    if(pet_ann > 0.0){
        wi = grid->prate_sfc_ann / pet_ann;
    }else{
        wi = 1.0;
    }
    if(wi < 0.01){
        wi = 0.01;
    }
    wi = log10(wi);
    
    mat = grid->tmp_sfc_am;
    
    /* model4 */
    lll = 0.571 + 0.637*wi - 0.498*wi*wi*wi - 0.254*wi*wi + 0.003*mat;
    if(lll < -2.0){
        lll = -2.0;
    }
    if(lll > 2.0){
        lll = 2.0;
    }
    
    /* model5 */
    lll = 0.532191 + 0.593486*wi - 0.441269*wi*wi*wi - 0.245266*wi*wi + 0.0164922*mat - 0.000620463*mat*mat;
    if(lll < -2.0){
        lll = -2.0;
    }
    if(lll > 2.0){
        lll = 2.0;
    }
    
    loct->est_maxlai = pow(10.0, lll);
}
