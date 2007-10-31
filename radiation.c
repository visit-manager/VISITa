/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/*********** functions of light environment ************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/*** solar declination at the middle day (15th day) of month ***/
double sl_dec(
	struct Grid *grid
){
	double sl_dec;
	/* day of the year of the 15th day in each month*/
	double doy[12] = {15.0, 46.0, 74.0, 105.0, 135.0, 166.0, 196.0, 227.0, 258.0, 288.0, 319.0, 349.0}; 

	sl_dec = 23.45*sin((doy[grid->m]-80.0)*360.0/370.0*dTr); /** **/
	
	return (sl_dec);	
}

/*********** solar hight at midday, in degree ************/
double sl_hgt(
	struct Grid *grid
){
	double aaa, bbb;
	
	aaa = sin(grid->lat*dTr)*sin(grid->sl_dec[grid->m]*dTr)+cos(grid->lat*dTr)*cos(grid->sl_dec[grid->m]*dTr);
	aaa = (aaa<=1.0)?aaa:1.0; 
	aaa = (aaa>=-1.0)?aaa:-1.0;
	bbb = asin(aaa)*rTd;
	
	return (bbb);
}

/************ day length, or photoperiod of vegetation production *************/
double dlen(
	struct Grid *grid
){
	double ho, sr, ss, dl, ha;
	
	/** hour-angle when the sun-angle equals to zero **/
	ho = -sin(grid->lat*dTr)*sin(grid->sl_dec[grid->m]*dTr)/cos(grid->lat*dTr)/cos(grid->sl_dec[grid->m]*dTr); 
	ho = (ho<=1.0)?ho:1.0; ho=(ho>=-1.0)?ho:-1.0;
	ha = acos(ho);
	
	sr = 12.0-ha*rTd/15.0; /** time of sunrise **/
	ss = 12.0+ha*rTd/15.0; /** time of sunset **/
	dl = ss-sr; /** day length **/
	
	return (dl);
}

/*********** shortwave radiation at the atmosphere-top *************/
double top_rad(
	struct Grid *grid
){
	double doy[12] = {15.0, 46.0, 74.0, 105.0, 135.0, 166.0, 
		196.0, 227.0, 258.0, 288.0, 319.0, 349.0};
	double ge, dlt, dtc, ho;
	double aa, bb, cc, dd, ee, gg, hh, ii;
	
	ge = 2.0*PI/365.0*doy[grid->m];
	
	aa=-0.399912*cos(ge)+0.070257*sin(ge);
	bb=-0.006758*cos(2.0*ge)+0.000907*sin(2.0*ge);
	cc=-0.002697*cos(3.0*ge)+0.00148*sin(3.0*ge);
	dlt=(0.006918+aa+bb+cc);
	
	/** relative distance between the sun and the earth **/
	dd=0.034221*cos(ge)+0.00128*sin(ge);
	ee=0.000719*cos(2.0*ge)+0.000077*sin(2.0*ge);
	dtc=1.00011+dd+ee;
		 
	/* solar constant = 4.921 MJ/m2 =1367 W/m2=1.96 cal/cm2/min
		 = 6151.5 micro mol photons/ m2 / s */	
	gg = 1367.0; 
	
	if(SC==1){
		gg *= 1.01;
	}else if(SC==2){
		gg *= 0.99;
	}else if(SC==3){
		if(grid->climy>=1990){
			gg *= 1.03;
		}
	}else if(SC==4){
		if(grid->climy>=1990){
			gg *= 0.97;
		}
	}
	
	/* holizontally incident radiation at the top of the atmosphere */
	hh=sin(dlt)*sin(grid->lat*dTr); 
	ii=cos(dlt)*cos(grid->lat*dTr); 
	ho=gg*dtc*(hh+ii); 
	ho=(ho>=0.0)?ho:0.0;
			
	return(ho);
}

/********** global radiation at the ground surface ************/
double gl_rad(
	struct Grid *grid
){
	double cloudiness, jj, hh;
		
	/* radiation at the ground, including cloud effect based on the empirical Equbal's equation */
	cloudiness = grid->tcdc_clm[grid->m];
	/* jj=0.803-0.34*cloudiness-0.458*cloudiness*cloudiness; */ /* Black's */
	jj = 0.8964-0.5392*cloudiness; /* new regression based on NCEP/NCAR data*/
	jj = (jj<=1.0)?jj:1.0; 
	jj = (jj>=0.0)?jj:0.0;
	
	hh = grid->top_rad[grid->m]*jj; 
	
	return(hh);
}

/********* photosynthetically active radiation ***********/
double par(
	struct Grid *grid
){
	double kt, hd, dd, par;
	double e2p_b, e2p_d;
	
	/* constant, after McCree (1971) */
	e2p_b = 4.6; /* W/m2 to micro-mol photon /m2/s for diffused radiation*/
	e2p_d = 4.2; /* W/m2 to micro-mol photon /m2/s for beam radiation*/
				
	/* photosynthetically active radiation (par) in the global radiation, 
	based on the empirical Tooming's equation */	
	if(grid->top_rad[grid->m]!=0.0){
		kt = grid->gl_rad[grid->m]/grid->top_rad[grid->m];
		dd = 0.958-0.982*kt;
		dd = (dd>0.01)?dd:0.01;
		dd = (dd<=1.0)?dd:1.0;
		hd = grid->gl_rad[grid->m]*dd; /** fraction of diffused radiation**/
		
		/* variable conversion factor after Dye et al. (2003) */
		if(D_PAR==1){
			e2p_d = 4.5886*dd/(0.010773+dd);
			e2p_d = (e2p_d>4.2)?e2p_d:4.2;
			
			e2p_b = 4.576 - 0.033144*dd;		/* 050409 */
		}
		
		/* fraction of PAR */
		grid->par_be[grid->m] = 0.43*(grid->gl_rad[grid->m]-hd); 
		grid->par_de[grid->m] = 0.57*hd; 
		
		/* conversion from W/m2 to micro-mol photon /m2/s */
		grid->par_bp[grid->m] = grid->par_be[grid->m]*e2p_b; 
		grid->par_dp[grid->m] = grid->par_de[grid->m]*e2p_d; 
		
		par = 0.43*(grid->gl_rad[grid->m]-hd)*e2p_b + 0.57*hd*e2p_d;
	}else{
		par = 0.0;
	}
	
	return (par);
}

/************* net radiation of canopy and ground surface *******************/
void net_rad(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Echar *echar
){
	double aaa, bbb, ccc, ddd1, ddd2, eee, ee_c3, ee_c4, fff;
	double net_long, rad_net_p, rad_net_g, c3_canopy, c4_canopy, kmono_c3, kmono_c4;
	double transmittance, ground;
	
	/* transmittance=0.1; */ /*2003-06-27*/
	transmittance=0.12;

	/** longwave budget : modified 2002/12/25, based on Budyko (1971) **/
	aaa=pow((grid->tmp_2m[grid->m]+ZAT), 4.0) * 5.6703 / 100000000.0;
	if(loct->vp[grid->m]>0.1&&loct->vp[grid->m]<40.0){
		bbb = 0.39 - 0.058*sqrt( loct->vp[grid->m]*760.0/1013.0 );
	}else if(loct->vp[grid->m]<=0.1){
		bbb = 0.39 - 0.058*sqrt( 0.1*760.0/1013.0 );
	}else if(loct->vp[grid->m]>=40.0){
		bbb = 0.39 - 0.058*sqrt( 40.0 );
	}
	ccc = 1.0-0.65*grid->tcdc_clm[grid->m];
	net_long = aaa*bbb*ccc;
	loct->rad_net_long[grid->m]=net_long;
	
	/** soil surface albedo **/
	(echar->soil).albedo[grid->m] = albedo_soil(loct, &(echar->soil));
	
	ee_c3 = loct->C3ptn[grid->m]*(echar->c3).eK[grid->m]*(mass->c3).lai[grid->m];
	ee_c4 = loct->C4ptn[grid->m]*(echar->c4).eK[grid->m]*(mass->c4).lai[grid->m];
	
	eee = ee_c3+ee_c4;
	ground = exp(-1.0*eee);
	c3_canopy = (1.0-ground)*loct->C3ptn[grid->m];
	c4_canopy = (1.0-ground)*loct->C4ptn[grid->m];
	loct->albedo_sfc[grid->m] = (echar->soil).albedo[grid->m]*ground + 
					(echar->c3).albedo*c3_canopy + (echar->c4).albedo*c4_canopy;
	
	ddd1 = exp(-1.0*eee*(1.0-transmittance)); /*2003-06-27*/
	ddd2 = exp(-1.0*eee); /*2003-06-27*/
	
	loct->rad_net_short[grid->m] = (1.0 - loct->albedo_sfc[grid->m])*grid->gl_rad[grid->m];
	kmono_c3 = irr_attn(grid, loct, &(echar->c3));
	kmono_c4 = irr_attn(grid, loct, &(echar->c4));
	loct->fapar_mono[grid->m] = loct->C3ptn[grid->m]*(1.0-(echar->c3).albedo)*(1.0-exp(kmono_c3*(mass->c3).lai[grid->m])) 
						+ loct->C4ptn[grid->m]*(1.0-(echar->c4).albedo)*(1.0-exp(kmono_c4*(mass->c4).lai[grid->m]));
	
	/** global radiation under the canopy or at the soil surface **/
	loct->gl_rad_g[grid->m] = grid->gl_rad[grid->m]*ddd1;
	
	/** net radiation of plant canopy **/
	fff = loct->C3ptn[grid->m]*(echar->c3).albedo+loct->C4ptn[grid->m]*(echar->c4).albedo;
	rad_net_p = (1.0-fff)*(1.0-ddd1)*grid->gl_rad[grid->m]-net_long*(1.0-ddd2);
	rad_net_p = (rad_net_p>=0.0)?rad_net_p:0.0;
	loct->rad_net_p[grid->m] = rad_net_p;
	
	/** net radiation of soil surface **/
	rad_net_g = (1.0-(echar->soil).albedo[grid->m])*ddd1*grid->gl_rad[grid->m] - net_long*ddd2;
	rad_net_g = (rad_net_g>=0.0)?rad_net_g:0.0;
	loct->rad_net_g[grid->m] = rad_net_g;
}

/*********** albedo of soil surface **************/
double albedo_soil(
	struct Loct *loct, 
	struct Schar *schar
){
	double albedo;
	
	/* a function of snow accumulation */
	albedo = schar->albedo0 + (0.7-schar->albedo0)/(1+exp(-0.05*(loct->snwa-70.0)));
	albedo = (albedo>0.05)?albedo:0.05; 
	albedo = (albedo<0.75)?albedo:0.75;
	
	return(albedo);
}

