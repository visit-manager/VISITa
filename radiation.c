/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* functions of light environment ******************************************/
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

extern short SENS_RAD;

/* solar declination at the middle day (15th day) of month *********/
double f_solar_decl(
	struct Grid *grid
){
	double sl_dec;
	/* day of the year of the 15th day in each month*/
	double doy[ASTEP] = {15.0, 46.0, 74.0, 105.0, 135.0, 166.0, 
				196.0, 227.0, 258.0, 288.0, 319.0, 349.0}; 

	sl_dec = 23.45 * sin((doy[grid->m] - 80.0)*360.0/370.0*dTr); /** **/
	
	return (sl_dec);	
}

/* solar hight at midday, in degree ********************************/
double f_solar_hgt(
	struct Grid *grid
){
	double aaa, bbb;
	
	aaa = sin(grid->lat*dTr) * sin(grid->sl_dec[grid->m]*dTr)
			+ cos(grid->lat*dTr) * cos(grid->sl_dec[grid->m]*dTr);
	aaa = (aaa<=1.0)?aaa:1.0; 
	aaa = (aaa>=-1.0)?aaa:-1.0;
	bbb = asin(aaa)*rTd;
	
	return (bbb);
}

/* day length, or photoperiod of vegetation production *************/
double f_day_length(
	struct Grid *grid
){
	double ho, sr, ss, dl, ha;
	
	/** hour-angle when the sun-angle equals to zero **/
	ho = -sin(grid->lat * dTr) * sin(grid->sl_dec[grid->m] * dTr)/
			cos(grid->lat * dTr) / cos(grid->sl_dec[grid->m] * dTr);
	ho = (ho<=1.0)?ho:1.0; ho=(ho>=-1.0)?ho:-1.0;
	ha = acos(ho);
	
	sr = 12.0 - ha*rTd/15.0; /** time of sunrise **/
	ss = 12.0 + ha*rTd/15.0; /** time of sunset **/
	dl = ss-sr; /** day length **/
	
	return (dl);
}

/* shortwave radiation at the atmosphere-top ***********************/
double f_top_rad(
	struct Grid *grid, 
    short ha
){
	double doy[ASTEP] = {15.0, 46.0, 74.0, 105.0, 135.0, 166.0, 
		196.0, 227.0, 258.0, 288.0, 319.0, 349.0};
	double ge, dlt, dtc, ho;
	double aa, bb, cc, dd, ee, gg, hh, ii;
	
	ge = 2.0 * PI / 365.0 * doy[grid->m];
	
	aa = -0.399912 * cos(ge) + 0.070257 * sin(ge);
	bb = -0.006758 * cos(2.0 * ge) + 0.000907 * sin(2.0 * ge);
	cc = -0.002697 * cos(3.0 * ge) + 0.00148 * sin(3.0 * ge);
	dlt = (0.006918 + aa + bb + cc);
	
	/** relative distance between the sun and the earth **/
	dd = 0.034221 * cos(ge) + 0.00128 * sin(ge);
	ee = 0.000719 * cos(2.0*ge) + 0.000077 * sin(2.0*ge);
	dtc = 1.00011 + dd + ee;
		 
	/* solar constant = 4.921 MJ/m2 =1367 W/m2=1.96 cal/cm2/min
		 = 6151.5 micro mol photons/ m2 / s */	
	gg = SLC;
    
    /* EX SRM: 2013/06/04 by A.Ito *******************/
	if(EX_SRM == 1 && grid->climy >= 2000){
        gg -= 2.6;
    }
    if(EX_SRM == 2 && grid->climy >= 2000){
        gg -= 4.5;
    }
    if(EX_SRM == 3 && grid->climy >= 2000){
        gg -= 6.0;
    }
    if(EX_SRM == 4 && grid->climy >= 2000){
        gg -= 8.5;
    }
    
    /* experiment for SRM by reflector */
    /* added: 2014/07/06 by A.Ito     */
    /* if(SCENARIO_ID == 3012 || SCENARIO_ID == 3033 || SCENARIO_ID == 3091){
        if(grid->climy >= 2020){
            gg -= ((double)grid->climy-2020.0) * 4.5/50.0;
        }
    } */
    
	if(SC == 1){
		gg *= 1.01;
	}else if(SC == 2){
		gg *= 0.99;
	}else if(SC == 3){
		if(grid->climy>=1990){
			gg *= 1.03;
		}
	}else if(SC == 4){
		if(grid->climy>=1990){
			gg *= 0.97;
		}
	}
    
	/* holizontally incident radiation at the top of the atmosphere */
	hh = sin(dlt) * sin(grid->lat * dTr); 
	ii = cos(dlt) * cos(grid->lat * dTr) * cos((double)ha * dTr); 
	ho = gg * dtc * (hh + ii); 
	ho = (ho >= 0.0)?ho:0.0;
			
	return(ho);
}

/* global radiation at the ground surface *********************/
double f_gl_rad(
	struct Grid *grid
){
	double cloudiness, jj, hh;
		
	/* radiation at the ground, including cloud effect based on
	 the empirical Equbal's equation */
	cloudiness = grid->tcdc_clm[grid->m];
	/* jj=0.803-0.34*cloudiness-0.458*cloudiness*cloudiness; */ /* Black's */
	jj = 0.8964 - 0.5392 * cloudiness; /* new regression based on NCEP/NCAR data*/
	jj = (jj<=1.0)?jj:1.0; 
	jj = (jj>=0.0)?jj:0.0;
	
	hh = grid->top_rad[grid->m] * jj; 

    /* EX SRM: 2013/06/04 by A.Ito *******************/
    if(EX_SRM == 11 && grid->climy >= 2000){
        hh -= 2.6;
    }
    if(EX_SRM == 12 && grid->climy >= 2000){
        hh -= 4.5;
    }
    if(EX_SRM == 13 && grid->climy >= 2000){
        hh -= 6.0;
    }
    if(EX_SRM == 14 && grid->climy >= 2000){
        hh -= 8.5;
    }
	
	hh = (hh<=SLC)?hh:SLC;
	hh = (hh>=0.0)?hh:0.0;
	
	return(hh);
}

/* photosynthetically active radiation ***************************/
double f_par(
	struct Grid *grid
){
	double kt, hd, dd, par;
	double e2p_b, e2p_d;
	
	/* constant, after McCree (1981) */
	/* McCree, K. J. 1981. Photosynthetically active radiation. 
	Pages 41-55 in O. L. Lange, P. S. Nobel, and C. B. Osmond, editors. 
	Encyclopedia of Plant Physiology. Springer, Berlin. */
	e2p_b = 4.6; /* W/m2 to micro-mol photon /m2/s for diffused radiation*/
	e2p_d = 4.2; /* W/m2 to micro-mol photon /m2/s for beam radiation*/
	
	/* photosynthetically active radiation (par) in the global radiation, 
	based on the empirical Tooming's equation */	
	/* ref. Iqbal, M. 1983. An introduction to solar radiation. 
	Academic Press, Toronto. */
	if(grid->top_rad[grid->m] > 0.0){
		/* surface / top ratio */
		kt = grid->gl_rad[grid->m]/grid->top_rad[grid->m];
		
		/* new estimation of diffuse radiation: 2008/09/08 by A.Ito */
		if(DIF_SRB == 1){
			if((grid->srb_dif_rr * grid->srb_dif_rr) > 0.25){
				dd = grid->srb_dif_aa + grid->srb_dif_bb * kt;
			}else{
				/* global average */
				dd = 1.306833 - 1.250070 * kt;
			}
		}else{
			dd = 0.958 - 0.982 * kt;
		}
		
		if(SENS_RAD == 1){
			dd *= 1.1;
		}
		if(SENS_RAD == 2){
			dd *= 0.9;
		}
		
		dd = (dd>0.01)?dd:0.01;
		dd = (dd<=1.0)?dd:1.0;
		
		/** diffused radiation **/
		hd = grid->gl_rad[grid->m] * dd;
		
		/* variable conversion factor after Dye et al. (2004) */
		/*
		Dye, D.G., 2004. Spectral composition and quantum-to-energy 
		ratio of diffuse photosynthetically active radiation under 
		diverse cloud conditions. Journal of Geophysical Research 
		109, 10.1029/2003JD004251.
		*/
		if(D_PAR == 1){
			e2p_d = 4.5886 * dd / (0.010773 + dd);
			e2p_d = (e2p_d>4.2)?e2p_d:4.2;
			
			e2p_b = 4.576 - 0.033144 * dd;		/* 050409 */
		}
		
		/* fraction of PAR ******/
		/* beam */
		grid->par_be[grid->m] = 0.43 * (grid->gl_rad[grid->m] - hd);
		/* diffuse */
		grid->par_de[grid->m] = 0.57 * hd;
		
		/* conversion from W/m2 to micro-mol photon /m2/s */
		/* beam */
		grid->par_bp[grid->m] = grid->par_be[grid->m] * e2p_b;
		/* diffuse */
		grid->par_dp[grid->m] = grid->par_de[grid->m] * e2p_d; 
		
		/* total: micro-mol photon /m2/s */
		par = grid->par_bp[grid->m] + grid->par_dp[grid->m];
	}else{
		par = 0.0;
	}
	
	return (par);
}

/* photosynthetically active radiation ***************************/
void f_par_h(
	struct Grid *grid,
    struct Loct *loct
){
	double kt, hd, dd;
	double e2p_b, e2p_d, par_be, par_de;
	
	/* constant, after McCree (1981) */
	/* McCree, K. J. 1981. Photosynthetically active radiation. 
	Pages 41-55 in O. L. Lange, P. S. Nobel, and C. B. Osmond, editors. 
	Encyclopedia of Plant Physiology. Springer, Berlin. */
	e2p_b = 4.6; /* W/m2 to micro-mol photon /m2/s for diffused radiation*/
	e2p_d = 4.2; /* W/m2 to micro-mol photon /m2/s for beam radiation*/
	
	/* photosynthetically active radiation (par) in the global radiation, 
	based on the empirical Tooming's equation */	
	/* ref. Iqbal, M. 1983. An introduction to solar radiation. 
	Academic Press, Toronto. */
	if(grid->top_rad[grid->m] > 0.0){
		/* surface / top ratio */
		kt = grid->gl_rad[grid->m]/grid->top_rad[grid->m];
		
		/* new estimation of diffuse radiation: 2008/09/08 by A.Ito */
		if(DIF_SRB == 1){
			if((grid->srb_dif_rr * grid->srb_dif_rr) > 0.25){
				dd = grid->srb_dif_aa + grid->srb_dif_bb * kt;
			}else{
				/* global average */
				dd = 1.306833 - 1.250070 * kt;
			}
		}else{
			dd = 0.958 - 0.982 * kt;
		}
		
		if(SENS_RAD == 1){
			dd *= 1.1;
		}
		if(SENS_RAD == 2){
			dd *= 0.9;
		}
		
		dd = (dd>0.01)?dd:0.01;
		dd = (dd<=1.0)?dd:1.0;
		
		/** diffused radiation **/
		hd = grid->gl_rad[grid->m] * dd;
		
		/* variable conversion factor after Dye et al. (2004) */
		/*
		Dye, D.G., 2004. Spectral composition and quantum-to-energy 
		ratio of diffuse photosynthetically active radiation under 
		diverse cloud conditions. Journal of Geophysical Research 
		109, 10.1029/2003JD004251.
		*/
		if(D_PAR == 1){
			e2p_d = 4.5886 * dd / (0.010773 + dd);
			e2p_d = (e2p_d>4.2)?e2p_d:4.2;
			
			e2p_b = 4.576 - 0.033144 * dd;		/* 050409 */
		}
		
		/* fraction of PAR ******/
		/* beam */
		par_be = 0.43 * (grid->gl_rad[grid->m] - hd);
		/* diffuse */
		par_de = 0.57 * hd;
		
		/* conversion from W/m2 to micro-mol photon /m2/s */
		/* beam */
		loct->ppfdb_h[grid->h] = par_be * e2p_b;
		/* diffuse */
		loct->ppfdd_h[grid->h] = par_de * e2p_d;
		
		/* total: micro-mol photon /m2/s */
		loct->ppfd_h[grid->h] = loct->ppfdb_h[grid->h] + loct->ppfdd_h[grid->h];
	}else{
		loct->ppfdb_h[grid->h] = 0.0;
		loct->ppfdd_h[grid->h] = 0.0;
		loct->ppfd_h[grid->h] = 0.0;
	}
}

/* net radiation of canopy and ground surface *******************************/
void f_net_rad(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Mass *mass, 
	struct Echar *echar
){
    short nn;
	double aaa, bbb, ccc, ddd1, ddd2, eee, ee_c3, ee_c4, fff;
	double net_long, rad_net_p, rad_net_g, c3_canopy, c4_canopy;
	double transmittance, ground;
    double albedo_base, albedo_var, rn_short_base, rn_short_var;
    double tsfc, tsfc_base, tsfc_var,crit, tt1, tt2, latheat, snsheat, dtsfc;
	
	/* transmittance=0.1; */ /*2003-06-27*/
	transmittance = 0.12;
    
	/** longwave budget : modified 2002/12/25, based on Budyko (1971) **/
	aaa = pow((grid->tmp_2m[grid->m] + ZAT), 4.0) * SBC;
    bbb = 0.39;
	if(loct->vp[grid->m]>0.1 && loct->vp[grid->m]<40.0){
		bbb = 0.39 - 0.058 * sqrt(loct->vp[grid->m]*  760.0/1013.0 );
	}else if(loct->vp[grid->m] <= 0.1){
		bbb = 0.39 - 0.058 * sqrt( 0.1*760.0/1013.0 );
	}else if(loct->vp[grid->m] >= 40.0){
		bbb = 0.39 - 0.058 * sqrt( 40.0*760.0/1013.0 );
	}
	ccc = 1.0 - 0.65 * grid->tcdc_clm[grid->m];
	net_long = aaa * bbb * ccc;
	loct->rad_net_long[grid->m] = net_long;
	
	/** soil surface albedo **/
	(echar->soil).albedo[grid->m] = albedo_soil(loct, &(echar->soil));
	
	ee_c3 = loct->c3ptn[grid->m] * (echar->c3).eK[grid->m] * (mass->c3).lai[grid->m];
	ee_c4 = loct->c4ptn[grid->m] * (echar->c4).eK[grid->m] * (mass->c4).lai[grid->m];
	
	eee = ee_c3 + ee_c4;
	ground = exp(-1.0*eee);
	c3_canopy = (1.0 - ground) * loct->c3ptn[grid->m];
	c4_canopy = (1.0 - ground) * loct->c4ptn[grid->m];
	loct->albedo_sfc[grid->m] = (echar->soil).albedo[grid->m] * ground +
					(echar->c3).albedo * c3_canopy + (echar->c4).albedo * c4_canopy;
    albedo_base = loct->albedo_sfc[grid->m];
    /* loct->xx8[grid->m] = albedo_base; */
	
    /* albedo perturbation: 2012/12/29 by A.Ito */
    if(EX_ALBEDO == 1){
        loct->albedo_sfc[grid->m] += grid->albedo_pert[grid->m];
    }else if(EX_ALBEDO == 2){
        loct->albedo_sfc[grid->m] = grid->albedo_max[grid->m][grid->row/10][grid->col/10];
    }else if(EX_ALBEDO == 3){
        loct->albedo_sfc[grid->m] = grid->albedo_min[grid->m][grid->row/10][grid->col/10];
    }else if(EX_ALBEDO == 4){
        loct->albedo_sfc[grid->m] = grid->albedo_av[grid->m][grid->row/10][grid->col/10];
    }else if(EX_ALBEDO == 5){
        loct->albedo_sfc[grid->m] = grid->glbalbedo[grid->m];
    }
    
    if(loct->albedo_sfc[grid->m] > 0.99){
        loct->albedo_sfc[grid->m] = 0.99;
    }
    if(loct->albedo_sfc[grid->m] < 0.1){
        loct->albedo_sfc[grid->m] = 0.1;
    }
    albedo_var = loct->albedo_sfc[grid->m];
    /* loct->xx5[grid->m] = albedo_var; */
    
    /* 2015/08/10 by A.Ito *****/
    /* if(grid->phase == 0){
        grid->tmp_sfc[grid->m] = grid->tmp_sfc_a[grid->m];
        grid->tmp10_soil[grid->m] = grid->tmp10_soil_a[grid->m];
        grid->tmp200_soil[grid->m] = grid->tmp200_soil_a[grid->m];
    } */
    
    /* temperature change due to albedo: 2014/5/19 by A.Ito */
    if(EX_TVAR == 1 && EX_ALBEDO>=1){
 
        if(grid->phase == 0){
            grid->tmp_sfc[grid->m] = grid->tmp_sfc_a[grid->m];
            grid->tmp10_soil[grid->m] = grid->tmp10_soil_a[grid->m];
            grid->tmp200_soil[grid->m] = grid->tmp200_soil_a[grid->m];
        }
        
        /* latent heat, W m-2, approximated by the previous year's value */
        latheat = (loct->incep[grid->m] + loct->trspr[grid->m] + loct->evpr[grid->m])
                * (2.5*1000000.0) / DHN / HSN / MDN[grid->m];
        
        /* net shortwave radiation, W m-2 */
        //rn_short_base = (1.0 - albedo_base) * grid->gl_rad[grid->m];
        //rn_short_var = (1.0 - albedo_var) * grid->gl_rad[grid->m];

        rn_short_base = (1.0 - albedo_base) * loct->grad_d[grid->m];
        rn_short_var = (1.0 - albedo_var) * loct->grad_d[grid->m];
        
        /* base temperature: default albedo */
        nn = 0; crit = 10.0;
        tsfc = grid->tmp_sfc[grid->m];
        if(grid->phase == 0){
            tsfc = grid->tmp_sfc_a[grid->m];
        }
        while(nn<12 && crit>0.2){
        
            aaa = pow((tsfc + ZAT), 4.0) * SBC;
            snsheat = 1210.0 * (tsfc - grid->tmp_2m[grid->m]) / loct->r_aero[grid->m];
            tt1 = (rn_short_base - latheat) - aaa*bbb*ccc - snsheat;
            
            aaa = pow(((tsfc + 0.5) + ZAT), 4.0) * SBC;
            snsheat = 1210.0 * ((tsfc+0.5) - grid->tmp_2m[grid->m]) / loct->r_aero[grid->m];
            tt2 = (rn_short_base - latheat) - aaa*bbb*ccc - snsheat;
            
            if( fabs(tt2 - tt1) >= 0.5){
                tsfc = tsfc - tt1*0.5 / (tt2 - tt1);
                crit = fabs( tt1*0.5 / (tt2 - tt1) );
            }else{
                break;
            }
            nn++;
        }
        tsfc_base = tsfc;
        /* loct->xx1[grid->m] = tsfc;
        loct->xx6[grid->m] = snsheat; */
        
        /* varied temperature: default albedo */
        nn = 0; crit = 10.0;
        tsfc = grid->tmp_sfc[grid->m];
        if(grid->phase == 0){
            tsfc = grid->tmp_sfc_a[grid->m];
        }
        while(nn<12 && crit>0.2){
        
            aaa = pow((tsfc + ZAT), 4.0) * SBC;
            snsheat = 1210.0 * (tsfc - grid->tmp_2m[grid->m]) / loct->r_aero[grid->m];
            tt1 = (rn_short_var - latheat) - aaa*bbb*ccc - snsheat;
            
            aaa = pow(((tsfc + 0.5) + ZAT), 4.0) * SBC;
            snsheat = 1210.0 * ((tsfc+0.5) - grid->tmp_2m[grid->m]) / loct->r_aero[grid->m];
            tt2 = (rn_short_var - latheat) - aaa*bbb*ccc - snsheat;
            
            if( fabs(tt2 - tt1) >= 0.5){
                tsfc = tsfc - tt1*0.5 / (tt2 - tt1);
                crit = fabs( tt1*0.5 / (tt2 - tt1) );
            }else{
                break;
            }
            nn++;
        }
        tsfc_var = tsfc;
        /* loct->xx2[grid->m] = tsfc;
        loct->xx7[grid->m] = snsheat; */
        
        /* loct->xx8[grid->m] = loct->r_aero[grid->m];  */
        
        /* temperature change */
        dtsfc = tsfc_var - tsfc_base;
        
        if(dtsfc>=-10.0 && dtsfc<=10.0){
            ;
        }else if(dtsfc < -10.0){
            dtsfc = -10.0;
        }else if(dtsfc > 10.0){
            dtsfc = 10.0;
        }else{
            dtsfc = 0.0;
            /*  printf("*********************bad dtsfc %lf\n", dtsfc);  */
        }        
        /* loct->xx3[grid->m] = dtsfc; */
       
        grid->tmp_sfc[grid->m] += dtsfc;
        grid->tmp10_soil[grid->m] += dtsfc;
        grid->tmp200_soil[grid->m] += dtsfc; /* */
        
        /* loct->xx4[grid->m] = grid->tmp_sfc[grid->m];
        loct->xx5[grid->m] = grid->tmp10_soil[grid->m];

        loct->xx9[grid->m] = latheat; */

       /* Assumption: this surface/sub-surface temperature change does not
                       affect air temperature and humidity */
        
        aaa = pow(((grid->tmp_sfc[grid->m] + 0.5) + ZAT), 4.0) * SBC;
        loct->rad_net_long[grid->m] = net_long = aaa * bbb * ccc;
        
        loct->rad_net_short[grid->m] = rn_short_var;
    }else{
        loct->rad_net_short[grid->m] = (1.0 - loct->albedo_sfc[grid->m]) * grid->gl_rad[grid->m];
    }
    
    /****/
	ddd1 = exp(-1.0 * eee * (1.0 - transmittance)); /*2003-06-27*/
	ddd2 = exp(-1.0 * eee); /*2003-06-27*/
   
    /* absorbed PAR */
    loct->appfd_g[grid->m] = loct->c3ptn[grid->m] * (echar->c3).appfd_db[grid->m]
                           + loct->c4ptn[grid->m] * (echar->c4).appfd_db[grid->m];
 
    loct->ippfd_g[grid->m] = loct->c3ptn[grid->m] * (echar->c3).ppfd_db[grid->m]
           + loct->c4ptn[grid->m] * (echar->c4).ppfd_db[grid->m];

    /* fapar */
    if(loct->ippfd_g[grid->m] > 0.0){
        loct->fappfd_g[grid->m] = loct->appfd_g[grid->m] / loct->ippfd_g[grid->m];
    }else{
        loct->fappfd_g[grid->m] = 0.0;
    }
 
    /* added: 2013/01/10 by A.Ito */
    loct->nsw_d[grid->m] = (1.0 - loct->albedo_sfc[grid->m]) * loct->grad_d[grid->m];
	
	/** global radiation under the canopy or at the soil surface **/
	loct->gl_rad_g[grid->m] = grid->gl_rad[grid->m]*ddd1;
	
	/** net radiation of plant canopy, W m-2 **/
	fff = loct->c3ptn[grid->m]*(echar->c3).albedo + loct->c4ptn[grid->m]*(echar->c4).albedo;
    
    if(EX_ALBEDO == 1){
        fff += grid->albedo_pert[grid->m];
    }else if(EX_ALBEDO == 2){
        fff = grid->albedo_max[grid->m][grid->row/10][grid->col/10];
    }else if(EX_ALBEDO == 3){
        fff = grid->albedo_min[grid->m][grid->row/10][grid->col/10];
    }else if(EX_ALBEDO == 4){
        fff = grid->albedo_av[grid->m][grid->row/10][grid->col/10];
    }else if(EX_ALBEDO == 5){
        fff = grid->glbalbedo[grid->m];
    }
    fff = (fff<0.99)?fff:0.99;
    fff = (fff>0.01)?fff:0.01;
	rad_net_p = (1.0 - fff)*(1.0 - ddd1) * grid->gl_rad[grid->m] - net_long * (1.0 - ddd2);
	rad_net_p = (rad_net_p>=0.0)?rad_net_p:0.0;
	loct->rad_net_p[grid->m] = rad_net_p;
	
	/** net radiation of soil surface, W m-2 **/
    fff = (echar->soil).albedo[grid->m];
    if(EX_ALBEDO == 1){
        fff += grid->albedo_pert[grid->m];
    }else if(EX_ALBEDO == 2){
        fff = grid->albedo_max[grid->m][grid->row/10][grid->col/10];
    }else if(EX_ALBEDO == 3){
        fff = grid->albedo_min[grid->m][grid->row/10][grid->col/10];
    }else if(EX_ALBEDO == 4){
        fff = grid->albedo_av[grid->m][grid->row/10][grid->col/10];
    }else if(EX_ALBEDO == 5){
        fff = grid->glbalbedo[grid->m];
    }
    fff = (fff<0.99)?fff:0.99;
    fff = (fff>0.01)?fff:0.01;
	rad_net_g = (1.0 - fff) * ddd1 * grid->gl_rad[grid->m] - net_long*ddd2;
	rad_net_g = (rad_net_g >= 0.0)?rad_net_g:0.0;
	loct->rad_net_g[grid->m] = rad_net_g;
}

/* albedo of soil surface ********************************************/
double albedo_soil(
	struct Loct *loct, 
	struct Schar *schar
){
	double albedo;
	
	/* a function of snow accumulation */
	/* albedo = schar->albedo0 + (0.7 - schar->albedo0)/(1.0 + exp(-0.05*(loct->snwa - 70.0))); */
    /* revised: 2012/12/29 by A.Ito */
	albedo = schar->albedo0 + (0.95 - schar->albedo0)/(1.0 + exp(-0.05 * (loct->snwa - 75.0)));
	
	if(SENS_RAD == 3){
		albedo *= 1.1;
	}
	if(SENS_RAD == 4){
		albedo *= 0.9;
	}
    
	/* albedo = (albedo>0.05)?albedo:0.05;
	albedo = (albedo<0.75)?albedo:0.75; */
    
    /* revised: 2012/12/29 by A.Ito */
    albedo = (albedo>0.01)?albedo:0.01;
	albedo = (albedo<0.95)?albedo:0.95;
	
	return(albedo);
}

