/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/*                  Revised August 14, 2007					*/
/* Revised 2008 / 09 / 02 by A.Ito				*/

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
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux
){
	/* biogenic emission factors, micro g C/gdm/hour   */
	/* 
	Lathiére, J., D. A. Hauglustaine, A. D. Friend, N. de Noblet-Ducoudré, 
	N. Viovy, and G. A. Folberth, 2006: Impact of climate variability and 
	land use changes on global biogenic volatile organic compound emissions. 
	Atmospheric Chemistry and Physics, 6, 2129-2146.
	*/
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
	/* CO: Tao, Z. and A. K. Jain, 2005: Modeling of global biogenic emissions of 
	 key indirect greenhouse gases and their response to atmospheric CO2 increases
	 and changes in land cover and climate. 
	 Journal of Geophysical Research, 110, 10.1029/2005JD005874.
	*/
	double emit_potent_co[18] = {0.0,
		0.3, 0.3, 0.3, 0.3, 0.3, 0.36, 0.3, 0.3, 
		0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 0.3, 
		0.3, 0.3};
    /* Sesquiterpenes: micro g/m2/h */
    /* Guenther AB, Jiang X, Heald CL, Sakulyanontvittaya T, Duhl T, 
    Emmons LK, Wang X (2012) The Model of Emissions of Gases and Aerosols from 
    Nature version 2.1 (MEGAN2.1): an extended and updated framework for modeling 
    biogenic emissions. Geoscientific Model Development 5: 1471–1492. 
    doi: 10.5194/gmd-5-1471-2012    */
    /* added 2014/9/11 by A.Ito */
    double emit_potent_afarnesene[18] = {0.0,
        60.0, 60.0, 40.0, 40.0, 40.0, 40.0, 40.0, 40.0,
        3.0, 3.0, 40.0, 40.0, 3.0, 3.0, 3.0,
        4.0, 40.0};
    double emit_potent_bcaryophyllene[18] = {0.0,
        60.0, 60.0, 40.0, 80.0, 40.0, 80.0, 80.0, 50.0,
        1.0, 1.0, 50.0, 50.0, 1.0, 1.0, 1.0,
        4.0, 50.0};
    double emit_potent_othersesqui[18] = {0.0,
        120.0, 20.0, 100.0, 120.0, 100.0, 120.0, 120.0, 100.0,
        2.0, 2.0, 100.0, 100.0, 2.0, 2.0, 2.0,
        2.0, 100.0};
	long f, idveg;
	double foliar_dens, leaf_temp, parday, prm_ensen;
	double f_ppfd, f_temp_isopr, f_temp_monotrp, f_phenology;
	double aa, bb, cc, cc2, dd, laiage[49], t_lai, total_closs;
	extern double MDN[ASTEP];
	
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
    /* revised by A.Ito: 2014/03/27 */
    if(loct->v_type == 1){
        idveg = grid->veg_sage;
    }else if(loct->v_type == 2){
        idveg = 16;
    }else{
        idveg = grid->veg_sage;
    }

	/* foliar density, g d.m. C / m2   *********************/
	foliar_dens = ((mass->c3).fol * loct->c3ptn[grid->m] 
				+ (mass->c4).fol * loct->c4ptn[grid->m]) * 100.0 * dmTc;
	
	/* light factor *********************/
	parday = grid->par[grid->m];
	f_ppfd = 0.0027 * 1.066 * parday / 
		sqrt(1.0 + 0.0027 * 0.0027 * parday * parday) * 0.5;
	
	/* temperature factor ***********/
    /* leaf temperature in K */
	leaf_temp = grid->tmp_sfc[grid->m] + ZAT;
    if(leaf_temp > 150.0 && leaf_temp < 350.0){
        /* isoprene */
        dd = 95000.0*(leaf_temp - 303.15)/(8.314 * leaf_temp * 303.15);
        aa = exp(dd);
        dd = 230000.0*(leaf_temp - 314.0)/(8.314 * leaf_temp * 303.15);
        bb = 0.961 + exp(dd);
        f_temp_isopr = aa/bb;
        /* monoterpene */
        f_temp_monotrp = exp(0.09*(leaf_temp - 303.15));
    }else{
        /* exceptional leaf temperature */
        f_temp_isopr = 0.0;
        f_temp_monotrp = 0.0;
    }
	
	/* leaf aging factor ***********/
	/* f_phenology = 0.75; */ /* conventional value */
	
	t_lai = 0.0;
	for(f=0;f<=48;f++){
		laiage[f] = (echar->c3).fleaf_age[f]*loct->c3ptn[grid->m] 
					+ (echar->c4).fleaf_age[f]*loct->c4ptn[grid->m];
		t_lai += laiage[f];
	}
	for(f=0;f<=48;f++){
		if(t_lai > 0.0){
			laiage[f] *= 1.0/t_lai;
		}else{
			laiage[f] = 0.0;
		}
	}
	
	switch(idveg){
		/* evergreen */
		case 1: case 3: case 4: case 6: case 8: case 9: case 10: 
		case 11: case 12: case 13: case 14: case 15: 
			
			f_phenology = 0.05 * laiage[0];
			f_phenology += 0.5 * (laiage[1] + laiage[2]);
			for(f=3;f<=24;f++){  /* changed f<=36 to f<=12 to f<=24: 2008/09/17 */
				/* f_phenology += 1.1 * laiage[f]; */
				f_phenology += 1.2 * laiage[f]; /* 2014/04/04 by A.Ito */
			}
			for(f=25;f<=48;f++){  /* changed f=37 to f<=13 to f<=25: 2008/09/17 */
				/* f_phenology += 0.4 * laiage[f]; */ 
				f_phenology += 0.5 * laiage[f]; /* 2014/04/04 by A.Ito */
			}
			
			break;
		/* deciduous */
		case 2: case 5: case 7: case 16:

			f_phenology = 0.05 * laiage[0];
			f_phenology += 0.5 * laiage[1];
			for(f=2;f<=10;f++){  /* changed f<=8 to f<=6: 080613 */
				/* f_phenology += 1.1 * laiage[f]; */
				f_phenology += 1.2 * laiage[f]; /* 2014/04/04 by A.Ito */
			}
			for(f=11;f<=18;f++){  /* changed f=9 to f=7: 080613 */
				/* f_phenology += 0.4 * laiage[f]; */
				f_phenology += 0.5 * laiage[f];
			}
			for(f=19;f<=48;f++){  /* changed f=9 to f=7: 080613 */
				f_phenology += 0.1 * laiage[f];
			}

			break;
		default:
			f_phenology = 0.75;
			break;
	}
	
	/***********************************************************************/
	cc = foliar_dens * MDN[grid->m] * grid->dlen[grid->m];
    cc2 = loct->lai[grid->m] * MDN[grid->m] * grid->dlen[grid->m];
    
    /* parameter ensemble: 2014/11/19 by A.Ito */
    prm_ensen = 1.0;
    if(PARAM_PTB == 8){
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
	
	/* VOC emission, micro g C m-2 month-1  */
	flux->voc_isopr_g97[grid->m] = emit_potent_isopr[idveg] * cc * f_ppfd * f_temp_isopr * f_phenology * prm_ensen;
	flux->voc_monotrp_g97[grid->m] = emit_potent_monotrp[idveg] * cc * f_temp_monotrp * f_phenology * prm_ensen;
	flux->voc_methanl_g97[grid->m] = emit_potent_methanl[idveg] * cc * f_temp_monotrp * f_phenology * prm_ensen;
	flux->voc_acetone_g97[grid->m] = emit_potent_acetone[idveg] * cc * f_temp_monotrp * f_phenology * prm_ensen;
	flux->voc_actaldhd_g97[grid->m] = emit_potent_actaldhd[idveg] * cc * f_temp_monotrp * f_phenology * prm_ensen;
	flux->voc_frmardhd_g97[grid->m] = emit_potent_frmardhd[idveg] * cc *f_temp_monotrp * f_phenology * prm_ensen;
	flux->voc_formacd_g97[grid->m] = emit_potent_formacd[idveg] * cc * f_temp_monotrp * f_phenology * prm_ensen;
	flux->voc_acetacd_g97[grid->m] = emit_potent_acetacd[idveg] * cc * f_temp_monotrp * f_phenology * prm_ensen;
	flux->voc_co_g97[grid->m] = emit_potent_co[idveg] * cc * f_temp_monotrp * f_phenology * prm_ensen;
    /* added 2014/9/11 by A.Ito */
	flux->voc_afarnesene[grid->m] = emit_potent_afarnesene[idveg] * cc2 * f_temp_monotrp * f_phenology * prm_ensen;
	flux->voc_bcaryophyllene[grid->m] = emit_potent_bcaryophyllene[idveg] * cc2 * f_temp_monotrp * f_phenology * prm_ensen;
	flux->voc_othersesqui[grid->m] = emit_potent_othersesqui[idveg] * cc2 * f_temp_monotrp * f_phenology * prm_ensen;
	
	/* carbon loss by BVOC emission: 2008/10/09 */
	if(NECB_BVOC == 1 && (EX_CCPL != 4 && EX_CCPL != 8)){
        /* revised (after comments by E.Kato): 2013/10/02 by A.Ito */
         
		total_closs = flux->voc_isopr_g97[grid->m] + flux->voc_monotrp_g97[grid->m] + flux->voc_methanl_g97[grid->m] + 
			flux->voc_acetone_g97[grid->m] + flux->voc_actaldhd_g97[grid->m] + flux->voc_frmardhd_g97[grid->m] + 
			flux->voc_formacd_g97[grid->m] + flux->voc_acetacd_g97[grid->m] + flux->voc_co_g97[grid->m] +
            flux->voc_afarnesene[grid->m] + flux->voc_bcaryophyllene[grid->m] + flux->voc_othersesqui[grid->m];
		
		(mass->c3).fol -= loct->c3ptn[grid->m] * total_closs/100000000.0;
        (mass->c3).mfol[grid->m] = (mass->c3).fol;
		if((mass->c3).fol < INT_C){
			(mass->c3).fol = INT_C;
		}
		
		(mass->c4).fol -= loct->c4ptn[grid->m] * total_closs/100000000.0;
        (mass->c4).mfol[grid->m] = (mass->c4).fol;
		if((mass->c4).fol < INT_C){
			(mass->c4).fol = INT_C;
		}
	}
}
