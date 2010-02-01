/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/*  initialize  environmental  characteristics  of  each  grid*/

#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* Atmospheric CO2 *********************************************/
void co2_trend(
	struct Grid *grid
){
	double base, inc, lgrd, season;
	double amplitude;
	double time, aa0, aa1, aa2, aa3, aa4, aa5;
	
	/** time of CO2 level  **/
	time = (double)(grid->co2y);
	
	/** BASE **/
	if(time<1990.0){
		/* fitting curve from observational data */
		aa0 = 1904299.0;
		aa1 = -3322.4242*pow(time, 1.0);
		aa2 = 1.6541596*pow(time, 2.0);
		aa3 = 1.3362655*pow(time, 3.0)/10000.0;
		aa4 = -3.0828809*pow(time, 4.0)/10000000.0;
		aa5 = 6.2121261*pow(time, 5.0)/100000000000.0;
		base = aa0+aa1+aa2+aa3+aa4+aa5;

		/** INC **/
		inc = (-52.11 + 0.026984*time)*((double)(grid->m)-5.5)/12.0;
		inc = (inc>=0.0)?inc:0.0;
	}else if(time>=1990.0 && time<=2100.0){
		/* IPCC SRES scenarios */
		if(CC_CD==0 && time>=2001.0){
			/* constant CO2 level */
			base = sres_co2[11];
		}else{
			base = sres_co2[grid->co2y-1990];
		}
		
		inc = 0.0;
	}

	/** latitudinal gradient **/
	lgrd = 1.6*(grid->lat/85.0);

	/** seasonal change **/
	amplitude = exp(0.04*grid->lat);
	if(grid->lat>=0.0){
		season = amplitude/2.0*sin(((double)(grid->m)-0.0)/12.0*2.0*PI);
	}else{   /*  if(grid->lat<0.0) */
		season = amplitude/2.0*sin(((double)(grid->m)+6.0)/12.0*2.0*PI);
	}
	
	/*  grid->bco2[grid->m]=base+lgrd+season;   */
	grid->bco2[grid->m] = base + lgrd + season;  /*   + 350.0  */
		
	/*  grid->d13c_bco2[grid->m]=-7.0;  */
	grid->d13c_bco2[grid->m] = -6.0 + (-0.02 * (base-280.0)) + (0.05*season);
	
	/* radiocarbon (d14C, D14C) of atmospheric CO2 *********************/
	/* added 2009/06/23 by A.Ito */
	/* fitting curves derived from data by U.S. ORNL CDIAC */
	if(time < 1964.0){
		grid->d14c_bco2[grid->m] = exp(0.6655532*(time - 1953.686));
	}else{
		grid->d14c_bco2[grid->m] = 1462.196 * exp(-0.06327057 * (time - 1954.282));
	}
}

/* Intra-canopy CO2-d13C *************************************/
void co2_in_canopy(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Flux *flux
){
	double closure_factor, source_factor;
	double ambient_d13c;

	/* CO2 concentration */
	loct->aco2[grid->m] = grid->bco2[grid->m];
	/* d13C - CO2 */
	closure_factor = 1.5*0.75*mass->lai_p/(1.5 + 0.75*mass->lai_p);
	source_factor = 2.0*0.5*flux->efflux_p/(2.0 + 0.5*flux->efflux_p);
	
	loct->cnpy_co2_recyc = 0.2*closure_factor*source_factor;  
	/*  0.3  011017  */  /*  0.1  011020  */  /*  0.5  011022  */
	
	if(loct->cnpy_co2_recyc>=0.0 && loct->cnpy_co2_recyc<=0.3){
		ambient_d13c = d13c_addition(grid->d13c_bco2[grid->m], 1.0-loct->cnpy_co2_recyc,  
							flux->d13c_efflux_p, loct->cnpy_co2_recyc);
		loct->d13c_aco2[grid->m] = ambient_d13c;
	}else{
		loct->d13c_aco2[grid->m] = grid->d13c_bco2[grid->m];
	}
}

