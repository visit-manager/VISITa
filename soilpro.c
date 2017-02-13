/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* functions of soil process */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* yearly prcesses in belowground ***********************************/
void soil_processes(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Schar *schar, 
	struct Smas *mass, 
	struct Sflx *flux
){
	double nn;
	extern double MDN[ASTEP];
	
	nn = MDN[grid->m];
	
	/** clear **/
	schar->ft_l[grid->m] = 0.0;
	schar->ft_h[grid->m] = 0.0;
	schar->fm_l[grid->m] = 0.0;
	schar->fm_h[grid->m] = 0.0;

	if(grid->veg_olson==0 || grid->veg_olson==33){
		/* no soil carbon in bare lands */
		flux->rl[grid->m] = 0.0;
		flux->rh[grid->m] = 0.0;
		flux->sf[grid->m] = 0.0;
		/* stable carbon isotope */
		flux->d13c_rl[grid->m] = mass->d13c_ltr;
		flux->d13c_rh[grid->m] = mass->d13c_msl;
		flux->d13c_sf[grid->m] = mass->d13c_ltr;
		
		mass->ltr = 0.0;
		mass->msl = 0.0;
        
        /* DOC leaching */
		flux->doc_boyer[grid->m] = 0.0;
		mass->doc = 0.0;
	}else{
		/* soil respiration of litter layer */
		flux->rl[grid->m] = nn * frl(grid, loct, schar, mass);
		/* soil respiration of mineral soil and humus */
		flux->rh[grid->m] = nn * frh(grid, loct, schar, mass);
		/* soil decomposition from upper litter to lower mineral soil */
		flux->sf[grid->m] = fsf(grid, schar, flux);
		/* stable carbon isotope */
		flux->d13c_rl[grid->m] = mass->d13c_ltr;
		flux->d13c_rh[grid->m] = mass->d13c_msl;
		flux->d13c_sf[grid->m] = mass->d13c_ltr;
		
		/* stable carbon isotope */
		mass->d13c_ltr = d13c_addition(mass->d13c_ltr, mass->ltr, flux->d13c_lL[grid->m], flux->lL[grid->m]);
		mass->d13c_msl = d13c_addition(mass->d13c_msl, mass->msl, flux->d13c_sf[grid->m], flux->sf[grid->m]);	
		
		/* mass balance */
		mass->ltr += flux->lL[grid->m]-flux->rl[grid->m]-flux->sf[grid->m]; 
		mass->msl += flux->sf[grid->m]-flux->rh[grid->m];
		
		mass->ltr = (mass->ltr>=0.0)?mass->ltr:0.0;
		mass->msl = (mass->msl>=0.0)?mass->msl:0.0;
		
        /* DOC leaching */
		f_doc_boyer(grid, loct, mass, flux);
	}
	/* total soil respiration */
	flux->hr[grid->m] = flux->rl[grid->m] + flux->rh[grid->m];
	flux->d13c_hr[grid->m] = d13c_addition(flux->d13c_rl[grid->m], flux->rl[grid->m], 
										   flux->d13c_rh[grid->m], flux->rh[grid->m]);
	
	/* monthly values */
	mass->ltr_m[grid->m] = mass->ltr;
	mass->msl_m[grid->m] = mass->msl;
	mass->soil[grid->m] = mass->ltr + mass->msl;
	
	/* stable carbon isotope */
	mass->d13c_ltr_m[grid->m] = mass->d13c_ltr;
	mass->d13c_msl_m[grid->m] = mass->d13c_msl;
	mass->d13c_soil[grid->m] = d13c_addition(mass->ltr, mass->d13c_ltr, mass->msl, mass->d13c_msl);
	
	/* d14c: added by A.Ito (2009/07/12) *********/
	/* d14c: revised by A.Ito (2009/11/17) *********/
	if((flux->sf[grid->m] + mass->msl) > 0.0){
		mass->d14c_msl = (mass->d14c_ltr*flux->sf[grid->m] + mass->d14c_msl*mass->msl) / 
						(flux->sf[grid->m] + mass->msl);	
	}else{
		mass->d14c_msl = grid->d14c_bco2[grid->m];
	}
	
	if((flux->lL[grid->m] + mass->ltr) > 0.0){
		mass->d14c_ltr = (flux->d14c_lL[grid->m]*flux->lL[grid->m] + mass->d14c_ltr*mass->ltr) / 
						(flux->lL[grid->m] + mass->ltr);		
	}else{
		mass->d14c_ltr = grid->d14c_bco2[grid->m];
	}
	
	mass->d14c_msl_m[grid->m] = mass->d14c_msl;
	mass->d14c_ltr_m[grid->m] = mass->d14c_ltr;
	
	/* nitrogen flows *********************************************/
	/* mineralization of organic N */
	f_n_mineralz(grid, loct, mass, flux);
	
	/* NH3 volatilization */
	f_nh3_volatilization(grid, loct, schar, mass, flux);
	
	/* NO3- leaching */
	f_n_leaching(grid, loct, mass, flux);
	
	/* microbial immobilization */
	f_n_immoblz(grid, loct, schar, mass, flux);
	
	/* microbial abandonment */
	f_n_mcrb_abdn(grid, loct, schar, mass, flux);
}

/* soil respiration of litter layer *******************************/
double frl(
	struct Grid *grid, 
	struct Loct *loct,
	struct Schar *schar, 
	struct Smas *mass
){
	double rlto, to, rl, ftl, fwl, fal, fsm, p_scale;
	
	rlto = schar->rl/1000.0; 
	to = 15.0; /* specific rate at 15 deg C */
	/* temperature effect, exponential */
	/* ftl=exp(log(soil->qTl)/10.0*(grid->tmp10_soil[grid->m]-to)); */
	if(grid->tmp10_soil[grid->m]>-20.0){
		if(T_D == 0){
			ftl = 0.05 + 0.95 * exp(308.56*(1.0/56.02-1.0/(grid->tmp10_soil[grid->m] + 46.02))); /* control */
		}else if(T_D == 1){
			ftl = 0.05 + 0.95 * exp(308.56*1.3*(1.0/56.02-1.0/(grid->tmp10_soil[grid->m] + 46.02)));
		}else if(T_D == 2){
			ftl = 0.05 + 0.95 * exp(308.56*0.7*(1.0/56.02-1.0/(grid->tmp10_soil[grid->m] + 46.02)));
		}else if(T_D == 3){
			ftl = 0.05 + 0.95 * exp(308.56*(1.0/56.02-1.0/(grid->tmp10_soil[grid->m] + 46.02*1.3)));
		}else if(T_D == 4){
			ftl = 0.05 + 0.95 * exp(308.56*(1.0/56.02-1.0/(grid->tmp10_soil[grid->m] + 46.02*0.7)));
		}
	}else{
		ftl=0.05;
	}
    
    if(EX_SDTD == 1){
        if(grid->tmp10_soil[grid->m]>-20.0){
            ftl = 0.05 + 0.95 * exp(log(2.0)/10.0 * (grid->tmp10_soil[grid->m] - 10.0));
        }else{
            ftl = 0.05;
        }
    }else if(EX_SDTD == 2){
        if(grid->tmp10_soil[grid->m]>-20.0){
            ftl = 0.05 + 0.95 * exp(log(1.5)/10.0 * (grid->tmp10_soil[grid->m] - 10.0));
        }else{
            ftl = 0.05;
        }
    }
    
	/* acclimation */
	if(EX_ACCLM==1 && grid->phase>=1){
		ftl = schar->ft0_l[grid->m];
	}
	
    if(PARAM_PTB == 4){
        p_scale = 0.10;
        ftl *= 1.0 + p_scale*f_pert[13];
    }
	
	schar->ft_l[grid->m] = ftl;

	/* soil moisture effect, saturating */
	fwl = 0.8 * loct->sw30/(schar->kml*grid->field_cap1 + loct->sw30) + 0.2;
	/* soil apparence effect */
	/* 2009/02/22 by A.Ito */
	/* fal = 0.4*loct->soil_appr30*(1.0*schar->kmsl)/(schar->kmsl+loct->soil_appr30)+0.6; */
	fal = 0.4 * loct->soil_appr30/(schar->kmsl + loct->soil_appr30) + 0.6;
	/* take minimum */
	fsm = (fwl>fal)?fal:fwl;
	
	/* acclimation */
	if(EX_ACCLM==2 && grid->phase>=1){
		fsm = schar->fm0_l[grid->m];
	}
	
	schar->fm_l[grid->m] = fsm;

	/* printf("%lf %lf %lf\n",ftl,fwl,fsm); */

	rl = mass->ltr * rlto * ftl * fsm; 
	
	/* in case of too much emission, in order to avoid negative mass value */ 
	if((mass->ltr - rl*(1.0 + schar->me))<0.0){
		rl = mass->ltr/(1.0 + schar->me);
	}
	return(rl);	
}

/* soil respiration of mineral soil layer *****************************/
double frh(
	struct Grid *grid, 
	struct Loct *loct,
	struct Schar *schar, 
	struct Smas *mass
){
	double rhto, to, rh, fth, fwh, fah, fsm, p_scale;
	
	rhto = schar->rh/1000.0; 
	to = 15.0; /* specific rate at 15 deg C*/
	/* temperature effect, exponential */
	/* fth=exp(log(soil->qTh)/10.0*(grid->tmp200_soil[grid->m]-to)); */
	if(grid->tmp200_soil[grid->m] > -20.0){
		if(T_D == 0){
			fth = 0.05 + 0.95*exp(308.56*(1.0/56.02 - 1.0/(grid->tmp200_soil[grid->m] + 46.02))); /* control */
		}else if(T_D == 1){
			fth = 0.05 + 0.95*exp(308.56*1.3*(1.0/56.02 - 1.0/(grid->tmp200_soil[grid->m] + 46.02)));
		}else if(T_D == 2){
			fth = 0.05 + 0.95*exp(308.56*0.7*(1.0/56.02 - 1.0/(grid->tmp200_soil[grid->m] + 46.02)));
		}else if(T_D == 3){
			fth = 0.05 + 0.95*exp(308.56*(1.0/56.02 - 1.0/(grid->tmp200_soil[grid->m] + 46.02*1.3)));
		}else if(T_D == 4){
			fth = 0.05 + 0.95*exp(308.56*(1.0/56.02 - 1.0/(grid->tmp200_soil[grid->m] + 46.02*0.7)));
		}
	}else{
		fth = 0.05;
	}
    
    if(EX_SDTD == 1){
        if(grid->tmp200_soil[grid->m]>-20.0){
            fth = 0.05 + 0.95 * exp(log(2.0)/10.0 * (grid->tmp200_soil[grid->m] - 10.0));
        }else{
            fth = 0.05;
        }
    }else if(EX_SDTD == 2){
        if(grid->tmp200_soil[grid->m]>-20.0){
            fth = 0.05 + 0.95 * exp(log(2.5)/10.0 * (grid->tmp200_soil[grid->m] - 10.0));
        }else{
            fth = 0.05;
        }
    }

	/* acclimation */
	if(EX_ACCLM==1 && grid->phase>=1){
		fth = schar->ft0_h[grid->m];
	}
	
    if(PARAM_PTB == 4){
        p_scale = 0.10;
        fth *= 1.0 + p_scale*f_pert[13];
    }

	schar->ft_h[grid->m] = fth;

	/* soil moisture effect, saturating */
	fwh = 0.8 * loct->sww/(schar->kmh * grid->field_cap2 + loct->sww) + 0.2;
	/* soil apparence effect */
	/* 2009/02/22 by A.Ito */
	/* fah = 0.4*loct->soil_apprw*(1.0*schar->kmsh)/(schar->kmsh + loct->soil_apprw) + 0.6;  */
	fah = 0.4 * loct->soil_apprw/(schar->kmsh + loct->soil_apprw) + 0.6;
    /* take minimum */
	fsm = (fwh>fah)?fah:fwh;
	
	/* acclimation */
	if(EX_ACCLM==2 && grid->phase>=1){
		fsm = schar->fm0_h[grid->m];
	}
	
	schar->fm_h[grid->m] = fsm;

	rh = mass->msl * rhto * fth * fsm; 
	
	/* in case of too much emission, in order to avoid negative mass value */ 
	if((mass->msl - rh) < 0.0){
		rh = mass->msl;
	}
	
	return(rh);
}

/* humus formation from litter to mineral soil ****************************/
double fsf(
	struct Grid *grid, 
	struct Schar *schar, 
	struct Sflx *flux
){
	double sf;
	
	/* humus formation is assumed to proceed in paralell with litter respiration, proportionally */
	sf = schar->me * flux->rl[grid->m];
	
	return(sf);	
}
