/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* initialize hydrology */
#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

/* water condition *******************************************************/
void f_waterbudget(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar
){
	double snow, rain, thaw;
	double gain1, gain2; 
	double dry_inx, retran;
	double aa, bb, cc;
	double cap_incep, baseflow;
	extern double MDN[12];
	
	/** snow/rain portions **/
	loct->snp[grid->m] = 1.0/(1.0+exp(0.75*(grid->tmp_2m[grid->m]-2.0)));
	snow = loct->snp[grid->m]*grid->prate_sfc[grid->m];
	rain = (1.0-loct->snp[grid->m])*grid->prate_sfc[grid->m];
	
	/** thaw water **/
	thaw = 1.0/(1.0 + exp(-0.375 * grid->tmp_sfc[grid->m]));	
	loct->thaw[grid->m] = loct->snwa*thaw;
	if(loct->thaw[grid->m] < 0.0){
		loct->thaw[grid->m] = 0.0;
	}
	
	/** canopy interception **/
	loct->pm_incep[grid->m] = pm_interception(grid, loct); /* potential */
	/* cap_incep=loct->lai[grid->m]*(double)(grid->mm[grid->m])*0.2; */
	/* cap_incep=loct->lai[grid->m]*(double)(grid->mm[grid->m])*0.25; */
	cap_incep = loct->lai[grid->m] * MDN[grid->m]*0.125; 
	loct->pm_incep[grid->m] = (loct->pm_incep[grid->m]<cap_incep)?loct->pm_incep[grid->m]:cap_incep;
	/* aa=0.75; */ /*2003-06-27*/
	/* aa=0.80; */ /*2003-06-27*/
	aa = 0.75; 
	bb = (rain)+loct->pm_incep[grid->m];
	cc = (rain)*loct->pm_incep[grid->m];
	loct->incep[grid->m] = (bb-sqrt(bb*bb-4.0*aa*cc))/(2.0*aa); /* actual */
	if(loct->incep[grid->m] < 0.0){
		loct->incep[grid->m] = 0.0;
	}
	
	rain -= loct->incep[grid->m];
	
	/** water balance 1 **/
	loct->snwa += snow - loct->thaw[grid->m];
	
	/* base flow from lower soil water *****/
	if(grid->tmp200_soil[grid->m] > 0.0){
		baseflow = 0.001 * loct->sww * MDN[grid->m];   /* after Sellers et al. (1996) */
	}else{
		baseflow = 0.0;
	}
	loct->sww -= baseflow;
	
	/** runoff 1, estimated by the bucket model **/
	dry_inx = grid->field_cap1 - loct->sw30;
	gain1 = rain+loct->thaw[grid->m];
	aa = gain1*gain1*gain1;
	bb = dry_inx*dry_inx*dry_inx;
	cc = aa+bb;
	cc = (cc>=0.0)?cc:0.0;
	loct->ro1[grid->m] = pow(cc,0.33333) - dry_inx;
	loct->ro1[grid->m] = (loct->ro1[grid->m]>=0.0)?loct->ro1[grid->m]:0.0;
	
	/** water balance 2 **/
	loct->sw30 += (rain + loct->thaw[grid->m]) - loct->ro1[grid->m];
	
	/** runoff 2, estimated by the bucket model **/
	dry_inx = grid->field_cap2 - loct->sww;
	gain2 = loct->ro1[grid->m]; 
	aa = gain2*gain2*gain2;
	bb = dry_inx*dry_inx*dry_inx;
	cc = aa+bb;
	cc = (cc>=0.0)?cc:0.0;
	loct->ro2[grid->m] = pow(cc, 0.33333) - dry_inx + loct->sww*0.0012*MDN[grid->m];
	loct->ro2[grid->m] = (loct->ro2[grid->m]>=0.0)?loct->ro2[grid->m]:0.0;
	
	/** water balance 3 **/
	loct->sww += loct->ro1[grid->m] - loct->ro2[grid->m]; 

	loct->ro2[grid->m] += baseflow;
	
	/** evaporation and transpiration **/
	loct->pm_evp[grid->m] = pm_evaporation(grid, loct); /* potential */
	/* aa=0.90; */ /*2003-06-27*/
	/* aa = 0.95; */ /*2003-11-25*/
	aa = 0.83; /* 2007-09-05*/
	bb = (loct->sw30) + loct->pm_evp[grid->m];
	cc = (loct->sw30) * loct->pm_evp[grid->m];
	loct->evpr[grid->m] = (bb - sqrt(bb*bb-4.0*aa*cc))/(2.0*aa); /* actual */
	loct->evpr[grid->m] = (loct->evpr[grid->m] > 0.0)?loct->evpr[grid->m]:0.0;
	
	loct->canopy_con[grid->m] = loct->c3ptn[grid->m]*(echar->c3).gc[grid->m] 
								+ loct->c4ptn[grid->m]*(echar->c4).gc[grid->m];
	loct->pm_trn[grid->m] = pm_transpiration(grid, loct); /* potential */
	/* aa=0.85; */ /*2003-06-27*/
	/* aa = 0.87; */ /*2003-06-27*/
	aa = 0.83;
	bb = (loct->sww)+loct->pm_trn[grid->m];
	cc = (loct->sww)*loct->pm_trn[grid->m];
	loct->trspr[grid->m] = (bb-sqrt(bb*bb-4.0*aa*cc))/(2.0*aa); /* actual */
	loct->trspr[grid->m] = (loct->trspr[grid->m]>0.0)?loct->trspr[grid->m]:0.0;
		
	/** water balance 4 **/
	loct->sw30 -= loct->evpr[grid->m];
	loct->sww -= loct->trspr[grid->m]; 

	/* re translocation ************/
	retran = (loct->sww*grid->field_cap1/grid->field_cap2 - loct->sw30)
				/(1.0+grid->field_cap1/grid->field_cap2);
	if(retran > 0.0){
		/* retran *= 0.4; */
		
		retran *= 0.35;  /* 080215 */
	}
	
	/** water balance 5 **/
	loct->sw30 += retran;
	loct->sww -= retran; 
	
	if(loct->sw30 < 0.0){
		loct->sw30 = 0.0;
	}
	if(loct->sww < 0.0){
		loct->sww = 0.0;
	}
	
	/** soil apperture **/
	loct->soil_appr30 = (grid->field_cap1 - loct->sw30)/grid->field_cap1;
	loct->soil_appr30 = (loct->soil_appr30<1.0)?loct->soil_appr30:1.0;
	loct->soil_appr30 = (loct->soil_appr30>0.0)?loct->soil_appr30:0.0;
	loct->soil_apprw = (grid->field_cap2 - loct->sww)/grid->field_cap2;
	loct->soil_apprw = (loct->soil_apprw<1.0)?loct->soil_apprw:1.0;
	loct->soil_apprw = (loct->soil_apprw>0.0)?loct->soil_apprw:0.0;
	
	/** put monthly variables **/
	loct->msnwa[grid->m] = loct->snwa;
	loct->msww[grid->m] = loct->sww;
	loct->msw30[grid->m] = loct->sw30;
}

/* air density *********************************/
double air_density(
	struct Grid *grid, 
	struct Loct *loct
){
	double aaa, bbb, ccc, air_density;
	
	aaa = ZAT/(grid->tmp_2m[grid->m]+ZAT);
	bbb = loct->prsr[grid->m]/1013.25;
	ccc = 1.0-0.378*loct->vp[grid->m]/loct->prsr[grid->m];
	air_density = 1.293*aaa*bbb*ccc;
	
	return(air_density);
}

/* saturation vapour pressure (in hPa) ****************/
double vap_pre_sat(
	struct Grid *grid
){
	double vps;
	
	/* printf("%ld %ld %lf\n", grid->climy, grid->m, grid->tmp_2m[grid->m]); */
	if(grid->tmp_2m[grid->m]>0.0){ /* at water surface */
		vps = 6.1078*pow(10.0, (7.5*grid->tmp_2m[grid->m])/(237.3+grid->tmp_2m[grid->m]));
	}else{ /* at ice surface */  /*  if(grid->tmp_2m[grid->m]<=0.0) */
		vps = 6.1078*pow(10.0, (9.5*grid->tmp_2m[grid->m])/(265.3+grid->tmp_2m[grid->m]));
	}
	vps = (vps>=0.0)?vps:0.0;
	
	return(vps);
}

/* slope of the saturation vapour pressure curve ***************/
double slope_vps(
	struct Grid *grid
){
	double slope, aaa, bbb, ccc;
	
	if(grid->tmp_2m[grid->m] > 0.0){ /* at water surface */
		aaa = 6.1078*(2500.0 - 2.4*grid->tmp_2m[grid->m]);
		bbb = 0.4615*(ZAT + grid->tmp_2m[grid->m])*(ZAT + grid->tmp_2m[grid->m]);
		ccc = pow(10.0, (7.5*grid->tmp_2m[grid->m])/(237.3 + grid->tmp_2m[grid->m]));
	}else{ /* at ice surface */  /* if(grid->tmp_2m[grid->m]<=0.0) */
		aaa = 6.1078*2834.0;
		bbb = 0.4615*(ZAT+grid->tmp_2m[grid->m])*(ZAT+grid->tmp_2m[grid->m]);
		ccc = pow(10.0, (9.5*grid->tmp_2m[grid->m])/(265.3+grid->tmp_2m[grid->m]));
	}
	slope = (aaa/bbb)*ccc;
	
	return(slope);
}

/* Aerodynamic Resistance ********************************/
double r_aero(
	struct Grid *grid
){
	double k_con;
	double r_aero;
	
	k_con = 0.41; /* von Karman's constant */
	
	grid->wnd_10m[grid->m] = (grid->wnd_10m[grid->m]>=0.1)?grid->wnd_10m[grid->m]:0.1;
	
	r_aero = (log(10.0)*log(10.0))/(k_con*k_con*grid->wnd_10m[grid->m]);
	r_aero = (r_aero>=0.0)?r_aero:0.0;
	
	return(r_aero);
}

/* Penman-Monteith Evaporation ***************************/
double pm_evaporation(
	struct Grid *grid, 
	struct Loct *loct
){
	double cp, psycon, spwt, lht, eta;
	double ggc, rc_g, evaporation;
	double aaa,bbb;
	extern double MDN[12];
	
	lht = 695.0; /** latent heat of water, in W h kg-1 **/
	spwt = loct->dnsa[grid->m]; /** density of air, in kg m-3 **/
	cp = 0.2813; /** specific heat of air, in W h kg-1 K-1 **/
	psycon = 0.667; /** psychlometer constant, in hPa K-1 **/
	eta = 0.0224*1.0/1000.0; /** unit conversion of conductance from mmol H2O m-2 s-1 to m s-1 **/
	
	/* ground resistance */
	/* ggc=1000.0*(1.0-(grid->whc-loct->sww)/grid->whc)+10.0; */ /*2003-06-27*/
	/* ggc=1000.0*(1.0-((grid->whc30+grid->whc)-(loct->sw30+loct->sww))/(grid->whc30+grid->whc))+10.0; */
	/* ggc=1000.0*(1.0-((grid->whc30+grid->whc)-(loct->sw30+loct->sww))/(grid->whc30+grid->whc))+100.0; */ /*2003-06-27*/
	ggc = 500.0*(1.0-((grid->field_cap1+grid->field_cap2)-(loct->sw30+loct->sww))/
				 (grid->field_cap1+grid->field_cap2))+10.0; /*2003-06-27*/
	rc_g = 1.0/(ggc*eta);
	
	aaa = (loct->slope_vps[grid->m]*loct->rad_net_g[grid->m])+(cp*spwt*loct->vpd[grid->m]/loct->r_aero[grid->m]);
	bbb = loct->slope_vps[grid->m]+psycon*(1.0+rc_g/loct->r_aero[grid->m]);	
	
	evaporation = MDN[grid->m]*grid->dlen[grid->m]*aaa/bbb/lht;
	evaporation = (evaporation>=0.0)?evaporation:0.0;
	
	return(evaporation);
}

/* Penman-Monteith Transpiration *******************************/
double pm_transpiration(
	struct Grid *grid, 
	struct Loct *loct
){
	double cp,psycon,spwt,lht,eta;
	double rc_p, transpiration;
	double aaa,bbb;
	double rn_transp;
	extern double MDN[12];
	
	lht = 695.0; /** latent heat of water, in W h kg-1 **/
	spwt = loct->dnsa[grid->m]; /** density of air, in kg m-3 **/
	cp = 0.2813; /** specific heat of air, in W h kg-1 K-1 **/
	psycon = 0.667; /** psychlometer constant, in hPa K-1 **/
	eta = 0.0224*1.0/1000.0; /** unit conversion of conductance from mmol H2O m-2 s-1 to m s-1 **/
	
	rn_transp = loct->rad_net_p[grid->m]-loct->incep[grid->m]/MDN[grid->m]/24.0*lht;

	if(loct->canopy_con[grid->m] > 0.0 && rn_transp > 0.0){
		/** canopy resistance **/
		rc_p = 1.0/(loct->canopy_con[grid->m]*eta);
		
		aaa = (loct->slope_vps[grid->m]*rn_transp) + (cp*spwt*loct->vpd[grid->m]/loct->r_aero[grid->m]);
		bbb = loct->slope_vps[grid->m]+psycon*(1.0 + rc_p/loct->r_aero[grid->m]);	
		transpiration = MDN[grid->m]*grid->dlen[grid->m]*aaa/bbb/lht;
	}else{
		transpiration = 0.0;
	}
	transpiration = (transpiration>=0.0)?transpiration:0.0;
	
	return(transpiration);
}

/* Penman-Monteith Transpiration *******************************/
double pm_interception(
	struct Grid *grid, 
	struct Loct *loct
){
	double cp,psycon,spwt,lht,eta;
	double rc_p, interception;
	double aaa,bbb;
	extern double MDN[12];
	
	lht = 695.0; /** latent heat of water, in W h kg-1 **/
	spwt = loct->dnsa[grid->m]; /** density of air, in kg m-3 **/
	cp = 0.2813; /** specific heat of air, in W h kg-1 K-1 **/
	psycon = 0.667; /** psychlometer constant, in hPa K-1 **/
	eta = 0.0224*1.0/1000.0; /** unit conversion of conductance from mmol H2O m-2 s-1 to m s-1 **/
	
	/** no vegetation resistance **/
	rc_p = 0.0;
	
	aaa = (loct->slope_vps[grid->m]*loct->rad_net_p[grid->m])
		+ (cp*spwt*loct->vpd[grid->m]/loct->r_aero[grid->m]);
	bbb = loct->slope_vps[grid->m] + psycon*(1.0 + rc_p/loct->r_aero[grid->m]);	
	interception = MDN[grid->m]*grid->dlen[grid->m]*aaa/bbb/lht;
	interception = (interception>=0.0)?interception:0.0;
	
	return(interception);
}

