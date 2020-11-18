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
void f_co2_trend(
	struct Grid *grid
){
	double base, inc, lgrd, season;
	double amplitude;
	double time, aa0, aa1, aa2, aa3, aa4, aa5;
	
	/** time of CO2 level  **/
	time = (double)(grid->co2y);
	
    base = 350.0;
    if(ISIMIP_RUN == 0 && GEOMIP_RUN == 0 && IMPRESSIONS_RUN <= 1){
        /** BASE **/
        if(time < 1990.0){
            /* fitting curve from observational data */
            aa0 = 1904299.0;
            aa1 = -3322.4242 * pow(time, 1.0);
            aa2 = 1.6541596 * pow(time, 2.0);
            aa3 = 1.3362655 * pow(time, 3.0)/10000.0;
            aa4 = -3.0828809 * pow(time, 4.0)/10000000.0;
            aa5 = 6.2121261 * pow(time, 5.0)/100000000000.0;
            base = aa0 + aa1 + aa2 + aa3 + aa4 + aa5;

            /** INC **/
            inc = (-52.11 + 0.026984 * time)*((double)(grid->m) - 5.5)/12.0;
            inc = (inc>=0.0)?inc:0.0;
        }else if(time >= 1990.0 && time <= 2100.0){
            /* IPCC SRES scenarios */
            if(CC_CD == 0 && time >= 2001.0){
                /* constant CO2 level */
                base = sres_co2[11];
            }else{
                base = sres_co2[grid->co2y - 1990];
            }
            
            inc = 0.0;
        }
    }else if(ISIMIP_RUN == 1 || ISIMIP_RUN == 2 || ISIMIP_RUN == 3 || ISIMIP_RUN == 4
            || GEOMIP_RUN == 1){
        /** ISI-MIP or GeoMIP or ICARUS runs **/
        switch(CO2S){
            case 1: /* rcp2.6 */
                base = aco2_3[grid->co2y - FDY_AGHG];
                break;
            case 2: /* rcp4.5 */
                base = aco2_4[grid->co2y - FDY_AGHG]; /* corrected: 2012/08/06 */
                break;
            case 3: /* rcp6.0 */
                base = aco2_1[grid->co2y - FDY_AGHG]; /* corrected: 2012/08/06 */
                break;
            case 4: /* rcp8.5 */
                base = aco2_2[grid->co2y - FDY_AGHG]; /* corrected: 2012/08/06 */
                break;
        }
    }else if(ISIMIP_RUN == 4){
        /** ISI-MIP 2b runs **/
        switch(CO2S){
            case 1: /* picontrol */
                base = aco2_3[grid->co2y - FDY_AGHG];
                break;
            case 2: /* pi-hist-rcp2.6 */
                base = aco2_4[grid->co2y - FDY_AGHG];
                break;
            case 3: /* pi-hist-rcp6.0 */
                base = aco2_1[grid->co2y - FDY_AGHG];
                break;
            case 4: /* pi-hist-rcp6.0 */
                base = aco2_2[grid->co2y - FDY_AGHG];
                break;
        }
    }else if(ISIMIP_RUN == 5){
        /** ISI-MIP 3a runs **/
        switch(SCENARIO_ID){
            case 5100: case 5101: case 5102:   /* obsclim */
            case 5110: case 5111: case 5112:
                base = aco2_1[grid->co2y - FDY_AGHG];
                break;
            case 5103: case 5104: case 5105: case 5106: case 5107: /* counterclim */
            case 5113: case 5114: case 5115: case 5116: case 5117:
                base = aco2_2[grid->co2y - FDY_AGHG];
                break;
        }
    }else if(ISIMIP_RUN == 6){
        /** ISI-MIP 3b runs **/
        switch(SCENARIO_ID){
            case 5120: case 5121: case 5122:
            case 5140: case 5141: case 5142:
            case 5160: case 5161: case 5162:
            case 5180: case 5181: case 5182:
            case 5200: case 5201: case 5202:
                base = aco2_1[grid->co2y - FDY_AGHG]; /* Pi control */
                break;
            case 5123: case 5124: case 5125: case 5126:
            case 5143: case 5144: case 5145: case 5146:
            case 5163: case 5164: case 5165: case 5166:
            case 5183: case 5184: case 5185: case 5186:
            case 5203: case 5204: case 5205: case 5206:
                base = aco2_2[grid->co2y - FDY_AGHG]; /* RCP2.6 */
                break;
            case 5127: case 5128: case 5129: case 5130:
            case 5147: case 5148: case 5149: case 5150:
            case 5167: case 5168: case 5169: case 5170:
            case 5187: case 5188: case 5189: case 5190:
            case 5207: case 5208: case 5209: case 5210:
                base = aco2_3[grid->co2y - FDY_AGHG]; /* RCP7.0 */
                break;
            case 5131: case 5132: case 5133: case 5134: case 5135: case 5136:
            case 5151: case 5152: case 5153: case 5154: case 5155: case 5156:
            case 5171: case 5172: case 5173: case 5174: case 5175: case 5176:
            case 5191: case 5192: case 5193: case 5194: case 5195: case 5196:
            case 5211: case 5212: case 5213: case 5214: case 5215: case 5216:
                base = aco2_4[grid->co2y - FDY_AGHG]; /* RCP8.5 */
                break;
        }
    }else if(IMPRESSIONS_RUN == 2 || IMPRESSIONS_RUN == 3){
        /** IMPRESSIONS time-dependent runs **/
        switch(CO2S){
            case 1:
                base = aco2_3[grid->co2y - FDY_AGHG];
                break;
            case 2:
                base = aco2_4[grid->co2y - FDY_AGHG];
                break;
            case 3:
                base = aco2_1[grid->co2y - FDY_AGHG];
                break;
            case 4:
                base = aco2_2[grid->co2y - FDY_AGHG];
                break;
        }
        
        switch(IMP_CD){
            case 1: base = 360.0; break;
            
            case 2: base = 423.0; break;
            case 3: base = 448.0; break;
            case 4: base = 474.0; break;
            case 5: base = 498.0; break;
            case 6: base = 517.0; break;
            case 7: base = 528.0; break;
            case 8: base = 532.0; break;

            case 9: base = 432.0; break;
            case 10: base = 469.0; break;
            case 11: base = 515.0; break;
            case 12: base = 572.0; break;
            case 13: base = 640.0; break;
            case 14: base = 718.0; break;
            case 15: base = 802.0; break;
            default: break;
        }
        
    }else{
        printf("BAD experimental setting!!\n");
    }
    
    /* assuming SRM + CDR : 2014/06/18 by A.Ito */
    if(GEOMIP_RUN == 1 && CC_CD == 3){
        if(grid->co2y > 2020){
        
            switch(CO2S){
                case 1:
                    base = aco2_3[grid->co2y - FDY_AGHG - (long)((grid->co2y-2020)/5)];
                    break;
                case 2:
                    base = aco2_4[grid->co2y - FDY_AGHG - (long)((grid->co2y-2020)/5)];
                    break;
                case 3:
                    base = aco2_1[grid->co2y - FDY_AGHG - (long)((grid->co2y-2020)/5)];
                    break;
                case 4:
                    base = aco2_2[grid->co2y - FDY_AGHG - (long)((grid->co2y-2020)/5)];
                    break;
            }
        }
    }
    
    if(GEOMIP_RUN == 1 && CC_CD == 4){
        /* cap at 450ppmv CO2: guardrail for ocean acidification */
        if(base > 450.0){
            base = 450.0;
        }
    }

	/** latitudinal gradient **/
	lgrd = 1.6*(grid->lat / 85.0);

	/** seasonal change **/
	amplitude = exp(0.04 * grid->lat);
	if(grid->lat >= 0.0){
		season = amplitude/2.0 * sin(((double)(grid->m) - 0.0)/12.0 * 2.0 * PI);
	}else{   /*  if(grid->lat<0.0) */
		season = amplitude/2.0 * sin(((double)(grid->m) + 6.0)/12.0 * 2.0 * PI);
	}
	
	/*  grid->bco2[grid->m]=base+lgrd+season;   */
	grid->bco2[grid->m] = base + lgrd + season;  /*   + 350.0  */
		
	/*  grid->d13c_bco2[grid->m]=-7.0;  */
    if(grid->co2y>=1850 && grid->co2y<=2015){
        grid->d13c_bco2[grid->m] = d13c_graven[grid->co2y - 1850];
    }else{
	    grid->d13c_bco2[grid->m] = -6.0 + (-0.02 * (base - 280.0)) + (0.05 * season);
    }
	
	/* radiocarbon (d14C, D14C) of atmospheric CO2 *********************/
	/* added 2009/06/23 by A.Ito */
	/* fitting curves derived from data by U.S. ORNL CDIAC */
	/* if(time < 1964.0){
		grid->d14c_bco2[grid->m] = exp(0.6655532*(time - 1953.686));
	}else{
		grid->d14c_bco2[grid->m] = 1462.196 * exp(-0.06327057 * (time - 1954.282));
	} */
    if(grid->co2y >= 1850 && grid->co2y <= 2015){
        if(grid->lat >= 30.0){
            grid->d14c_bco2[grid->m] = d14c1_graven[grid->co2y - 1850];
        }else if(grid->lat < 30.0 && grid->lat > -30.0){
            grid->d14c_bco2[grid->m] = d14c2_graven[grid->co2y - 1850];
        }else if(grid->lat <= -30.0){
            grid->d14c_bco2[grid->m] = d14c3_graven[grid->co2y - 1850];
        }
    }else if(grid->co2y < 1850){
        if(grid->lat >= 30.0){
            grid->d14c_bco2[grid->m] = d14c1_graven[0];
        }else if(grid->lat < 30.0 && grid->lat > -30.0){
            grid->d14c_bco2[grid->m] = d14c2_graven[0];
        }else if(grid->lat <= -30.0){
            grid->d14c_bco2[grid->m] = d14c3_graven[0];
        }
    }else if(grid->co2y > 2015){
        if(grid->lat >= 30.0){
            grid->d14c_bco2[grid->m] = d14c1_graven[165];
        }else if(grid->lat < 30.0 && grid->lat > -30.0){
            grid->d14c_bco2[grid->m] = d14c2_graven[165];
        }else if(grid->lat <= -30.0){
            grid->d14c_bco2[grid->m] = d14c3_graven[165];
        }
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
	closure_factor = 1.5 * 0.75 * mass->lai_p/(1.5 + 0.75 * mass->lai_p);
	source_factor = 2.0 * 0.5 * flux->efflux_p/(2.0 + 0.5 * flux->efflux_p);
	
	loct->cnpy_co2_recyc = 0.2 * closure_factor * source_factor;
	/*  0.3  011017  */  /*  0.1  011020  */  /*  0.5  011022  */
	
	if(loct->cnpy_co2_recyc>=0.0 && loct->cnpy_co2_recyc<=0.3){
		ambient_d13c = d13c_addition(grid->d13c_bco2[grid->m], 1.0-loct->cnpy_co2_recyc,  
							flux->d13c_efflux_p, loct->cnpy_co2_recyc);
		loct->d13c_aco2[grid->m] = ambient_d13c;
	}else{
		loct->d13c_aco2[grid->m] = grid->d13c_bco2[grid->m];
	}
}

