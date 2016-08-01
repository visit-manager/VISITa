/*	VISIT: Vegetation Integrative SImulation Tool						*/
/*  Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/*  Carbon cycle, erosion, biomass burning, land-use change,			*/
/*  CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/*	N2O code by M.Inatomi					*/

#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

#define NOTICE 0

/* formula of DAILY N2O emission  ****************************************/
/* NGAS by Parton et al. (1996)  */
/*
Parton WJ, Mosier AR, Ojima DS, Valentine DW, Schimel DS, 
Weier K, Kulmala AE (1996) Generalized model for N2 and N2O 
production from nitrification and denitrification. 
Global Biogeochemical Cycles 10:401-412
*/
void f_n2o_emit_ngas(
	struct Grid *grid, 
	struct Loct *loct,  
	struct Mass *mass, 
	struct Flux *flux
){
	double aa, bb, cc, dd, ee;
	double n_h2o, n_t, n_ph, n_nh4, day_n_n2o;
	double fd_wfps, fd_no3, fd_co2, dt;
	double fr_wfps, fr_no3, fr_co2;
	double day_d_n2, day_d_n2o;
	double nh4_soil;		/* micro g g-1*/
	double no3_soil;		
	double kmax;		/* g N ha-1 day-1 */ /* assumed */
	double nmax;							/* assumed */
	extern double MDN[ASTEP];
	double vmw_b, wfps_b;
    
    /* 2016/07/08 by A.Ito */
    /* 2016/08/01 by A.Ito */
    ee = 1.0 - exp(-4.0 * 0.3);
	
    /* 2016/06/12 by A.Ito */
    if(EX_NITROGEN == 3){
    
        if( grid->field_cap1 > 0.0 && grid->field_cap2 > 0.0){
            /* vmw_b = (loct->sw30 + loct->sww) / (grid->field_cap1 + grid->field_cap2); */
            vmw_b = ee * loct->sw30 / grid->field_cap1 + (1.0-ee) * loct->sww / grid->field_cap2;
        }else{
            /* vmw_b = (loct->sw30 + loct->sww) / 1500.0; */
            vmw_b = loct->sw30 / 300.0;
        }
        
    }else{
        /* vmw_b = (loct->sw30 + loct->sww) / 1500.0; */
        vmw_b = (loct->sw30 + loct->sww) / 1000.0;
    }
    
    wfps_b = vmw_b;
	/* wfps_b = vmw_b / (1.0 - grid->bulkdens/2.65); */
    
	if(wfps_b > 0.95){
		wfps_b = 0.95;
	}
	if(wfps_b < 0.05){
		wfps_b = 0.05;
	}
	
	/* added by A.Ito (2009/06/16) */
	/* revised 2014/11/27 by A.Ito */
	/* if(CALC_OLSON == 1 && (grid->veg_olson==29 || grid->veg_olson==30 || grid->veg_olson==31 || grid->veg_olson==32)){ */
	if(CALC_OLSON == 1 && (loct->v_type == 2)){
        /* cropland */
		/* kmax = 28.6; */
		/* kmax = 18.0; */ /* 22.5=>20.0=>18.0 2014/12/02 by A.Ito */
		/* nmax = 23.0; */ /* 30.0=>25.0=>23.0 2014/11/30 by A.Ito */
		/* kmax = 18.0; */ /* 22.0: 2016/05/30 by A.Ito */
		/* nmax = 30.0; */ /* 30.0: 2016/05/30 by A.Ito */
		kmax = 12.0; /* 22.0: 2016/05/30 by A.Ito */
		nmax = 20.0; /* 30.0: 2016/05/30 by A.Ito */
		/* 2009/06/15 by A.Ito */
        /* micro g g-1*/
		nh4_soil = (mass->soil).n_nh4*1000000.0/10000.0 /(grid->bulkdens*1000.0*1000.0);	
		no3_soil = (mass->soil).n_no3*1000000.0/10000.0 /(grid->bulkdens*1000.0*1000.0); /* */ /* low */
		/* nh4_soil = (mass->soil).n_no3*1000000.0/10000.0 /(grid->bulkdens*300.0*1000.0);	
		no3_soil = (mass->soil).n_nh4*1000000.0/10000.0 /(grid->bulkdens*300.0*1000.0); */	 /* high */
		
	}else{
		/* natural */
		/* kmax = 3.8; */ 
		/* kmax = 3.3; */   /* 3.8=>3.4=>3.3 2014/11/30 by A.Ito */
		/* nmax = 23.0; */  /* 30.0=>25.0=>23.0 2014/11/30 by A.Ito */
		/* kmax = 4.0; */   /* 3.8: 2016/05/30 by A.Ito */
		/* nmax = 30.0; */  /* 30.0: 2016/05/30 by A.Ito */
		kmax = 2.0;   /* 3.8: 2016/05/30 by A.Ito */
		nmax = 10.0;  /* 30.0: 2016/05/30 by A.Ito */
		/* 2009/06/15 by A.Ito */
        /* micro g g-1*/
		nh4_soil = (mass->soil).n_nh4*1000000.0/10000.0 /(grid->bulkdens*1000000.0);
		no3_soil = (mass->soil).n_no3*1000000.0/10000.0 /(grid->bulkdens*1000000.0); /* */ /* low */
		/* nh4_soil = (mass->soil).n_no3*1000000.0/10000.0 /(grid->bulkdens*600.0*1000.0);	
		no3_soil = (mass->soil).n_nh4*1000000.0/10000.0 /(grid->bulkdens*600.0*1000.0); */	 /* high */
	}
	
	/*
	nh4_soil = grid->total_n_1m * 0.3 * 0.01*1000000.0 /(grid->bulkdens*300.0*1000.0)*0.5;
	no3_soil = grid->total_n_1m * 0.3 * 0.01*1000000.0 /(grid->bulkdens*300.0*1000.0)*0.5;  */
    
    /* 2016/07/13 */
    /* nh4_soil *= 0.1;
    no3_soil *= 0.1; */
	
	/* NITRIFICATION **************************************************/
	/* water factor */
	/* Fig.(2a) in Parton et al. (1996) */
	if(grid->soiltexture == 1 || grid->soiltexture == 2){
		/* sandy */
		aa = (wfps_b - 1.70) / (0.55 - 1.70);
		bb = 3.22 * ((1.70 - 0.55) / (0.55 + 0.007));
		cc = (wfps_b + 0.007) / (0.55 + 0.007);
		n_h2o = pow(aa, bb) * pow(cc, 3.22);
	}else{
		/* medium and fine */
		aa = (wfps_b - 1.27) / (0.60 - 1.27);
		bb = 2.84 * ((1.27 - 0.60) / (0.60 - 0.0012));
		cc = (wfps_b - 0.0012) / (0.60 - 0.0012);
		n_h2o = pow(aa, bb) * pow(cc, 2.84);
	}
	
	/* temperature factor */
	/* Fig.(2b) in Parton et al. (1996) */
	n_t = -0.06 + 0.13 * exp(0.07 * grid->tmp10_soil[grid->m]);
	if(n_t < 0.0){
		n_t = 0.0;
	}
	
	/* pH factor */
	/* Fig.(2c) in Parton et al. (1996) */
	n_ph = 0.56 + atan(PI * 0.45 * (-5.0 + grid->soil_ph))/PI;
	if(n_ph < 0.0){
		n_ph = 0.0;
	}
	
	/* ammonium factor */
	/* Fig.(2d) in Parton et al. (1996) */
	n_nh4 = 1.0 - exp(-0.0105 * nh4_soil);
	if(n_nh4 < 0.0){
		n_nh4 = 0.0;
	}
	
	/* N2O emission through nitrification	*/
	/* Eq.(1) in Parton et al. (1996)		*/
	day_n_n2o = n_h2o * n_ph * n_t * (kmax + nmax * n_nh4);
	
	/* DENITRIFICATION ************************************************/
	/* Fig.(3a) in Parton et al. (1996) */
	if(grid->soiltexture == 1 || grid->soiltexture == 2){
		/* sandy */
		fd_wfps = 1.56 / pow(12.0, (16.0 / pow(12.0, 2.01*wfps_b)));
	}else if(grid->soiltexture == 4 || grid->soiltexture == 5){
		/* fine */
		fd_wfps = 60.0 / pow(18.0, (22.0 / pow(18.0, 1.06*wfps_b)));
	}else{
		/* medium */
		fd_wfps = 4.82 / pow(14.0, (16.0 / pow(14.0, 1.39*wfps_b)));
	}
	if(fd_wfps < 0.0){
		fd_wfps = 0.0;
	}
	
	/* Fig.(3b) in Parton et al. (1996) */
	fd_no3 = 11000.0 + (40000.0 * atan(PI * 0.002 * (no3_soil - 180.0)))/PI;
	if(fd_no3 < 0.0){
		fd_no3 = 0.0;
	}
	/* Fig.(3c) in Parton et al. (1996) */
	fd_co2 = (24000.0 / (1.0 + 200.0 / exp(0.35 * (flux->soil).hr[grid->m]*1000.0/MDN[grid->m]))) - 100.0;
	if(fd_co2 < 0.0){
		fd_co2 = 0.0;
	}
	
	/* Gas flux, Eq.(2) in Parton et al. (1996) */
	/* g N ha-1 day-1 */
	dt = fd_wfps * ((fd_no3>fd_co2)?fd_co2:fd_no3);
	if(dt < 0.0){
		dt = 0.0;
	}
	
	/* Fig.(5a) in Parton et al. (1996) */
	fr_wfps = 1.4 / pow(13.0, 17.0/pow(13.0, 2.2*wfps_b));
	if(fr_wfps < 0.0){
		fr_wfps = 0.0;
	}
	/* Fig.(5b) in Parton et al. (1996) */
	fr_no3 = (1.0 - (0.5 + (1.0 * atan(PI * 0.01 * (no3_soil - 190.0))) / PI)) * 25.0;
	if(fr_no3 < 0.0){
		fr_no3 = 0.0;
	}
	/* Fig.(5c) in Parton et al. (1996) */
	fr_co2 = 13.0 + (30.78 * atan(PI * 0.07 * ((flux->soil).hr[grid->m]*1000.0/MDN[grid->m] - 13.0))) / PI;
	if(fr_co2 < 0.0){
		fr_co2 = 0.0;
	}
	
	loct->xx1[grid->m] = n_nh4;
	loct->xx2[grid->m] = fd_wfps;
	loct->xx3[grid->m] = fd_no3;
	loct->xx4[grid->m] = fd_co2;
	loct->xx5[grid->m] = dt;
	loct->xx6[grid->m] = fr_wfps;
	loct->xx7[grid->m] = fr_no3;
	loct->xx8[grid->m] = fr_co2;
	loct->xx9[grid->m] = wfps_b; /* */
		
	/* N2O emission */
	/* Eqs.(3+4) in Parton et al. (1996) */
	day_d_n2o = dt / (1.0 + fr_wfps * ((fr_no3>fr_co2)?fr_co2:fr_no3));
	/* Eqs.(3+5) in Parton et al. (1996) */
    dd = (fr_no3>fr_co2)?fr_co2:fr_no3;
    if((fr_wfps * dd) > 0.0){
        ee = 1.0 + 1.0/(fr_wfps * dd);
        day_d_n2 = dt / ee;
    }else{
        day_d_n2 = 0.0;
    }
	
	/* g N20 ha-1 month-1 */
	(flux->soil).d_n2o_ntr_ngas[grid->m] = day_n_n2o * 44.0/28.0 * MDN[grid->m];
	(flux->soil).d_n2o_dnt_ngas[grid->m] = day_d_n2o * 44.0/28.0 * MDN[grid->m];
	
	/* g N2 ha-1 month-1 */
	(flux->soil).d_n2_ngas[grid->m] = day_d_n2 * MDN[grid->m];
	
	/* nitrification */
	/* (flux->soil).n_nitrif[grid->m] = day_n_n2o / 0.01 * MDN[grid->m]; */ /* revised by A.Ito (2009/07/18) */
	/* (flux->soil).n_nitrif[grid->m] = day_n_n2o / 0.02 * MDN[grid->m]; */ /* 2009/07/23 */
	/* (flux->soil).n_nitrif[grid->m] = day_n_n2o / 0.012 * MDN[grid->m]; */ /* 2010/03/30 */
	/* (flux->soil).n_nitrif[grid->m] = day_n_n2o / 0.005 * MDN[grid->m]; */ /* 2016/05/30 */
	/* (flux->soil).n_nitrif[grid->m] = day_n_n2o / 0.00291 * MDN[grid->m]; */ /* 2016/07/06 */
	(flux->soil).n_nitrif[grid->m] = day_n_n2o / 0.005 * MDN[grid->m]; /* 2016/07/06 */
    
    aa = (flux->soil).d_n2o_ntr_ngas[grid->m] + (flux->soil).n_nitrif[grid->m];
    if(aa > 0.5*(mass->soil).n_nh4){
        bb = 0.5*(mass->soil).n_nh4 / aa;
    
        (flux->soil).d_n2o_ntr_ngas[grid->m] *= bb;
        (flux->soil).n_nitrif[grid->m] *= bb;
    }
    
    aa = (flux->soil).d_n2o_dnt_ngas[grid->m] + (flux->soil).d_n2_ngas[grid->m];
    if(aa > 0.5*(mass->soil).n_no3){
        bb = 0.5*(mass->soil).n_no3 / aa;
    
        (flux->soil).d_n2o_dnt_ngas[grid->m] *= bb;
        (flux->soil).d_n2_ngas[grid->m] *= bb;
    }
    
	/* total ***************************************************/
	(flux->soil).d_n2o_ngas[grid->m] = (flux->soil).d_n2o_ntr_ngas[grid->m] + (flux->soil).d_n2o_dnt_ngas[grid->m];
}

/* Daily step CASA nitrogen trace gas emission from soil ****************************/
/* Potter et al. (1996) */
/*
Potter CS, Matson PA, Vitousek PM, Davidson EA (1996) 
Process modeling of controls on nitrogen trace gas emissions from soils 
worldwide. Journal of Geophysical Research 101:1361-1377
*/
void f_n2o_emit_casa(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Flux *flux
){
	double d_n_min;		/* nitrogen gross mineralization, g N ha-1 day-1 */
	double d_no, d_n2o, d_n2;		/* soil trace gas emission, g N ha-1 day-1 */
	double f_emit;
		
	/* fraction of gas emission per mineralization */
	/* 2009/06/15 by A.Ito */
	/* f_emit = 0.01; */ /* low */
	f_emit = 0.02; /* control */
	
	d_no = d_n2o = 0.0;
		
	/* N mineralization rate */
	/* net */
	/* d_n_min = (flux->soil).n_minerlz_lttr[grid->m] + (flux->soil).n_minerlz_hums[grid->m]
			- (flux->soil).n_immbl[grid->m]; */
	/* gross: 2010/04/06 by A.Ito */
	d_n_min = (flux->soil).n_minerlz_lttr[grid->m] + (flux->soil).n_minerlz_hums[grid->m];
	if(d_n_min < 0.0){
		d_n_min = 0.0;
	}
							
	/* d_n_min = 0.5 * ((flux->soil).rl[grid->m] / (mass->soil).ltr * ((mass->soil).ltr/40.0)
					+ (flux->soil).rh[grid->m] / (mass->soil).msl * ((mass->soil).msl/20.0))*1000000.0;  */
	
	/* Fig.3 of Potter et al. (1996) */
	if(loct->i_w[grid->m]<10.0){
		d_no = 0.1 * loct->i_w[grid->m] * d_n_min * f_emit;
		d_n2o = 0.0;
		d_n2 = 0.0;
	}else if(10.0 <= loct->i_w[grid->m] && loct->i_w[grid->m]<30.0){
		d_no = (-0.005 * loct->i_w[grid->m] + 1.05) * d_n_min * f_emit;
		d_n2o = (0.005 * loct->i_w[grid->m] - 0.05) * d_n_min * f_emit;
		d_n2 = 0.0;
	}else if(30.0 <= loct->i_w[grid->m] && loct->i_w[grid->m]<60.0){
		d_no = (-0.01 * loct->i_w[grid->m] + 1.2) * d_n_min * f_emit;
		d_n2o = (0.01 * loct->i_w[grid->m] - 0.2) * d_n_min * f_emit;
		d_n2 = 0.0;
	}else if(60.0 <= loct->i_w[grid->m] && loct->i_w[grid->m]<70.0){
		d_no = (-0.04 * loct->i_w[grid->m] + 3.0) * d_n_min * f_emit;
		d_n2o = (0.04 * loct->i_w[grid->m] - 2.0) * d_n_min * f_emit;
		d_n2 = 0.0;
	}else if(70.0 <= loct->i_w[grid->m] && loct->i_w[grid->m]<80.0){
		d_no = (-0.02 * loct->i_w[grid->m] + 1.6) * d_n_min * f_emit;
		d_n2o = (0.02 * loct->i_w[grid->m] - 0.6) * d_n_min * f_emit;
		d_n2 = 0.0;
	}else if(80.0 <= loct->i_w[grid->m] && loct->i_w[grid->m]<90.0){
		d_no = 0.0;
		d_n2o = (-0.02 * loct->i_w[grid->m] + 2.6) * d_n_min * f_emit;
		d_n2 = (0.02 * loct->i_w[grid->m] - 1.6) * d_n_min * f_emit;
	}else{  /*  if(90.0 <= loct->i_w[grid->m]) */
		d_no = 0.0;
		d_n2o = (-0.08 * loct->i_w[grid->m] + 8.0) * d_n_min * f_emit;
		d_n2 = (0.08 * loct->i_w[grid->m] - 7.0) * d_n_min * f_emit;
	}
	
	/* g NO ha-1 month-1 */
	(flux->soil).d_no_casa[grid->m] = d_no*30.0/14.0;
	/* g N2 ha-1 month-1 */
	(flux->soil).d_n2_casa[grid->m] = d_n2;

	/* g N20 ha-1 month-1 */
	(flux->soil).d_n2o_casa[grid->m] = d_n2o*44.0/28.0;
}
