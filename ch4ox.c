/*	VISIT: Vegetation Integrative SImulation Tool						*/
/*  Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/*  Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/*  Carbon cycle, erosion, biomass burning, land-use change,			*/
/*  CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

#include<stdio.h>
#include<math.h>
#include"structure.h"
#include"prototype.h"

#define NOTICE 0

extern double ach4_a1[553], ach4_a2[553], ach4_b1[553], ach4_b2[553];

/* Daily CASA soil moisture submodel **********************************/
void f_casa_mositure(
	struct Grid *grid, 
	struct Loct *loct
){
	double aa, bb;
	double fc;		/* field capacity */
	double pc;		/* pore space capacity */
		
	/* from SSSAJ86_Saxton */
	/*  aa = 0.002576;
	bb = -8.990332;  */
	
	switch(grid->soiltexture){
		case 1:	/* coarse */
			aa = 0.002;	bb = -5.48;		
			break;
		case 2:	/* medium coarse */
			aa = 0.002;	bb = -6.54;		
			break;
		case 3:	/* medium */
			aa = 0.013;	bb = -6.57;		
			break;
		case 4:	/* fine medium */
			aa = 0.006;	bb = -9.47;		
			break;
		case 5:	/* fine */
			aa = 0.004;	bb = -13.75;		
			break;
		default:	/* others: ocean, ice, organic */
			aa = 0.013;	bb = -6.57;		
			break;
	}
	
	/*  fc = 0.349224;  */
	fc = grid->field_cap1/300.0;  /* mm to fraction */
	pc = grid->pore_cap1/300.0;
	if(pc<fc){
		pc  = fc + 0.05;
	}
	
	loct->m_rdr[grid->m] = (1.0 + grid->a_sw)/(1.0 + grid->a_sw * pow(loct->m_vmc[grid->m], grid->b_sw));
	
	/* m_m in meter */
	if(grid->prate_sfc[grid->m] > loct->m_pet[grid->m]){
		loct->m_m[grid->m] =  loct->m_m_pre + ((grid->prate_sfc[grid->m] - loct->m_pet[grid->m])/1000.0);
	}else if(grid->prate_sfc[grid->m] < loct->m_pet[grid->m]){
		loct->m_m[grid->m] = loct->m_m_pre - ((loct->m_pet[grid->m] - grid->prate_sfc[grid->m])/1000.0) * (loct->m_rdr[grid->m]);
	}
	
	/* never negative: 070529 */
	if(loct->m_m[grid->m] < 0.0){
		loct->m_m[grid->m] = 0.0;
	}
	
	if(loct->m_m[grid->m] >= fc){
		loct->m_m[grid->m] = fc;
	}
	
	loct->m_e[grid->m] = (((grid->prate_sfc[grid->m] - loct->m_pet[grid->m])/1000.0) 
			- (fc - loct->m_m_pre)/1000.0)/(double)(grid->mm[grid->m]);

	/* i_w, moinsture scalar */
	loct->i_w[grid->m] = 0.0;
	if(loct->m_e[grid->m] > 0.0){
		/* Eq.3a in Potter et al. (1996) */
		loct->i_w[grid->m] = ((loct->m_e[grid->m] + fc) / pc) * 100.0;
	}else if(loct->m_e==0){
		/* Eq.3b in Potter et al. (1996) */
		loct->i_w[grid->m] = (loct->m_m[grid->m] / pc) * 100.0;
	}else{
		/* Eq.3b in Potter et al. (1996) */
		/*i_w = 0.0;*/
		loct->i_w[grid->m] = (loct->m_m[grid->m] / pc) * 100.0;
	}
	
	loct->i_w[grid->m] = (loct->i_w[grid->m]<0)?0:loct->i_w[grid->m];
	loct->i_w[grid->m] = (loct->i_w[grid->m]>=100.0)?100.0:loct->i_w[grid->m];
	
	loct->m_m_pre = loct->m_m[grid->m];
}

/* CH4 oxydation by Ridgewell et al. (1999) **********************************/
/*
Ridgwell AJ, Marshall SJ, Gregson K (1999) 
Consumption of atmospheric methane by soils: A process-based model. 
Global Biogeochemical Cycles 13:59-70
*/
void f_ch4oxy_ridgewell(
	struct Grid *grid, 
	struct Loct *loct,  
	struct Flux *flux
){
	double r_sm;
	double r_t;
	double r_n;
	double k_d;
	double co_ch4 = 1.72;	/* ppmv */
	double zd = 30.0;		/* depth, cm */
	double f = 616.9;		/* mg CH4 ppmv-1 cm-1 */
	double pd = 2.65;		/* particle density, mg cm-3 */
	double phi;				/* total pore volume, cm3 cm-3 */
	double e;				/* air-filled porosity, cm3 cm-3 */
	double g_soil;
	double g_t;
	double d0_ch4 = 0.196;
	double d_ch4;			/* diffusivity of CH4 in free air, STP of 0.196 cm2 s-1 */
	double k0 = 0.00087;	/* base oxydation rate constant for uncultivated moist soil at 0 deg C */
	double lcult;			/* fractional intensity of cultivation by SAGE */
	double d_j_ch4;			/*  */
	
	/*Atmospheric CH4 */
	co_ch4 = ach4_a1[grid->CO2y-1750]/1000.0;
	
	phi = 1.0 - (grid->bulkdens/pd);		/* total pore volume, cm3 cm-3 */
	if(NOTICE==1 && phi<=0.0){
		printf("!! ERROR: negative phi (ch4_oxy.c): %lf\n", phi);
	}
	
	/* moisture factor */
	if(((grid->prate_sfc[grid->m] + loct->m_sw[grid->m])/loct->m_pet[grid->m])>1.0){
		/* Eq. (11a) in Ridgewell et al. (1999) */
		r_sm = 1.0;
	}else{
		/* Eq. (11b) in Ridgewell et al. (1999) */
		r_sm = (grid->prate_sfc[grid->m] + loct->m_sw[grid->m])/loct->m_pet[grid->m];
	}
	/* temperature factor */
	if(grid->tmp10_soil[grid->m] < 0.0){
		/* Eq. (9a) in Ridgewell et al. (1999) */
		r_t = 0.0;
	}else{
		/* Eq. (9b) in Ridgewell et al. (1999) */
		r_t = exp(0.0693 * grid->tmp10_soil[grid->m] - 8.56 * pow(10.0, -7.0) * pow(grid->tmp10_soil[grid->m], 4.0));
	}
	/* nitrogen factor */
	/* Eq. (10) in Ridgewell et al. (1999) */
	lcult = grid->f_crop_con;
	r_n = 1.0 - (0.75 * lcult);
	
	/* Eq. (8) in Ridgewell et al. (1999) */
	k_d = k0 * r_sm * r_t * r_n;

	e = (phi - ((loct->m_sw[grid->m]/10.0)/30.0));
	/* Eq. (6) in Ridgewell et al. (1999) */
	if(e>0.0){
		g_soil = pow(phi, 4.0/3.0) * pow(e/phi, (1.5+3.0/grid->b_sw));
	}else{
		g_soil = 0.0;
	}
	/* Eq. (7) in Ridgewell et al. (1999) */
	g_t = 1.0 + 0.0055 * grid->tmp10_soil[grid->m];
	
	/* Eq. (5) in Ridgewell et al. (1999) */
	d_ch4 = g_soil * g_t * d0_ch4;
	/* Eq. (4) in Ridgewell et al. (1999) */
	if((d_ch4 + k_d*zd)>0.0){
		d_j_ch4 = ((co_ch4 * d_ch4)/ zd) * (1.0 - (d_ch4/(d_ch4 + k_d*zd))) * f;
	}else{
		d_j_ch4 = ((co_ch4 * d_ch4)/ zd) * f;
	}
	
	/* mg CH4 m-2 day-1 */
	(flux->soil).ch4oxy_ridg[grid->m] = d_j_ch4;
	/* mg CH4 m-2 day-1 */
	(flux->soil).ch4oxy_ridg[grid->m] *= (double)grid->mm[grid->m];
	
	/* removal of wetlands and paddy fields */
	(flux->soil).ch4oxy_ridg[grid->m] *= grid->f_upland;
}

/* CH4 oxydation by Potter et al. (1996) **********************************/
/*  
Potter CS, Davidson EA, Verchot LV (1996) Estimation of global 
biogeochemical controls and seasonality in soil methane consumption. 
Chemosphere 32:2219-2246
*/
void f_ch4oxy_casa(
	struct Grid *grid, 
	struct Loct *loct,  
	struct Flux *flux
){
	double j_ch4;		/* methane consumption flux, mg CH4 m-2 day-1 */
	double d_ch4;		/* diffusion coefficient of CH4 in the topsoil(0-0.3m) */
	double c_ch4;		/* concentration gradient of CH4 in the topsoil profile, ppmv cm-1 */
	double d0_ch4;		/* diffusion coefficient of CH4 in air */
	double d_d0;		/* normalized diffusivity in aggregated media */
	double ss;			/* solid phase of the aggregated medium, m3 m-3 bed volume */
	double aa;			/* intra-aggregate pore space, m3 m-3 bed volume */
	double pp;			/* inter-aggregate pore space, m3 m-3 bed volume */
	double th_a;		/* volume of water per unit bed volume contained in intra-arregate pore space, m3 m-3 */
	double th_p;		/* volume of water per unit bed volume contained in inter-arregate pore space, m3 m-3 */
	double s_wa, s_wp;		/* fractional liquid saturation of aa and pp components of the total void volume */
	double xx, yy, zz;		/* exponent terms */
	double fc;		/* field capacity */
	double pc;		/* pore space */
	double aaa, bbb, ccc, ddd, eee, fff, hhh, iii, jjj;
	
	/* from CHEM96_Potter */
	c_ch4 = 0.04;		
	/*邵ｲﾂfc = aa = grid->whc30/300.0;邵ｲﾂ*/
	fc = aa = grid->field_cap1/300.0;  
	pc = grid->pore_cap1/300.0;
	if(pc<fc){
		pc  = fc + 0.05;
	}
	pp = pc - aa;	/* inter-aggregate pore space */
	ss = 1.0 - pc;
			
	th_a = loct->m_m[grid->m];
	if(loct->m_m[grid->m] <= aa){
		th_p = 0.0;
	}else if(loct->m_m[grid->m] > aa && pp > loct->m_e[grid->m]){
		th_p = loct->m_e[grid->m];
	}else if(loct->m_m[grid->m] > aa && pp <= loct->m_e[grid->m]){
		th_p = pp;
	}
	
	s_wa = th_a / aa;
	s_wp = th_p / pp;
	
	aaa = pp - th_p;
	bbb = (aa - th_a) / (aa + ss);
	
	/* Eq. 4a in Potter (1996) */
	xx = 0.477 * pow(pp, 3.0) - 0.596 * pow(pp, 2.0) + 0.437 * pp + 0.564;
	/* Eq. 4b in Potter (1996) */
	yy = 0.477 * pow(aaa, 3.0) - 0.596 * pow(aaa, 2.0) + 0.437 * aaa + 0.564;
	/* Eq. 4c in Potter (1996) */
	zz = 0.477 * pow(bbb, 3.0) - 0.596 * pow(bbb, 2.0) + 0.437 * bbb + 0.564; /* 070427 */
	
	/* Eq. 3 in Potter (1996) */
	fff = pow((1.0 - s_wa), 2.0);
	jjj = (aa - th_a / aa + ss);
	hhh = (1.0 - pow(pp, 2.0 * xx));
	iii = (aaa - pow(aaa, 2.0 * yy));
	
	ccc = fff * pow(jjj, 2.0 * zz) * hhh * iii; /*************/
	ddd = fff * pow(jjj, 2.0) * hhh + iii;
	eee = pow((1.0 - s_wp), 2.0) * pow(aaa, 2.0 * yy);
	
	if(jjj>0){
		d_d0 = ccc / ddd + eee;
	}else if(jjj<=0){
		d_d0 = eee;
	}
	
	/* Eq. 2 in Potter (1996) */
	d0_ch4 = (0.194 * (0.9734 + 0.0055 * grid->tmp10_soil[grid->m]));
	
	d_ch4 = d0_ch4 * d_d0; 
	
	if(grid->tmp10_soil[grid->m] > -10.0){
		/* Eq. 1 in Potter (1996) */
		j_ch4 = d_ch4 * c_ch4 * (2.7 * pow(10.0, -4.0))* 30.0 * 24.0 * 3600.0;  /*  */
	}else{
		j_ch4 = 0.0;
	}

	/* mg CH4 m-2 day-1 */
	(flux->soil).ch4oxy_casa[grid->m] = j_ch4;
	/* mg CH4 m-2 month-1 */
	(flux->soil).ch4oxy_casa[grid->m] *= (double)grid->mm[grid->m];
	
	(flux->soil).ch4oxy_casa[grid->m] *= grid->f_upland;
	
	/*  if((flux->soil).ch4oxy_casa[grid->m]>1000.0 || (flux->soil).ch4oxy_casa[grid->m]<0.0){
		(flux->soil).ch4oxy_casa[grid->m] = 0.0;
	}  */
}

/* CH4 oxydation by Del Grosso et al. (2000) **********************************/
/*
Del Grosso, S. J., et al. (2000), General CH4 oxidation model and 
comparisons of CH4 oxidation in natural and managed systems, 
Global Biogeochemical Cycles, 14, 999-1019.
*/
void f_ch4oxy_delgrosso(
	struct Grid *grid, 
	struct Loct *loct,  
	struct Flux *flux
){
	double max_ch4oxy;
	double f_swater, f_tsoil, f_dopt, f_wfps;
	double fc, pc, vac, dopt;
	double wmin, wmax, wopt;
	double aa, bb, cc, dd, ee, ff;
	double tp1, tp2, tp3, tp4, tp5, tp6, tp7, tp8;
	double pfc, vfrac, theta_V, theta_P, theta_A, sw_p, s_wat;
	double ch4ox_nat, ch4ox_agr, ch4ox_com;
	
	/* grassland, tropical and coniferous forests ********************/
	
	pc = grid->pore_cap1/300.0;  /* porosity */
	fc = grid->field_cap1/300.0;  
	if(pc<fc){
		pc  = fc + 0.05;
	}
	
	/* inter-aggregate pore (air) space */
	vac = pc - fc;	
	/* percent water fraction, % */
	pfc = loct->wfps[grid->m]*100.0 / (fc / pc); 
	/* volumetric fraction */
	vfrac = (pfc > 100.0)?((loct->wfps[grid->m]*pc - fc)/vac):0.0;
	
	aa = (vfrac>1.0)?1.0:vfrac;
	theta_V = (pfc < 100.0)? (fc + aa * vac) : (pfc/100.0 * fc);
	
	theta_P = (theta_V < fc)? 0.0 : (theta_V - fc);
	
	theta_A = (theta_V > fc)? fc : theta_V;
	
	s_wat = ((theta_V / fc)<1.0)?(theta_V / fc):1.0;
	
	sw_p =  ((theta_P / vac)<1.0)?(theta_P / vac):1.0;
	
	if((1.0 - s_wat)>0.0){
		tp1 = pow((1.0 - s_wat), 2.0);
	}else{
		tp1 = 0.0;
	}
	
	tp2 = (fc - theta_A) / (fc + (1.0 - pc));
	
	if(tp2>0.0){
		tp3 = pow(tp2, (0.5 * tp2 + 1.16));
	}else{
		tp3 = 0.0;
	}
	
	tp4 = 1.0 - pow(vac, (0.5 * vac + 1.16));
	
	tp5 = vac - theta_P;
	
	if(tp5>0.0){
		tp6 = pow(tp5, (0.5 * tp5 + 1.16));
	}else{
		tp6 = 0.0;
	}
	
	tp7 = pow((1.0 - sw_p), 2.0);
	
	tp8 = tp1 * tp3 * tp4 * (tp5 - tp6) / (0.000001 + (tp1 * tp3 * tp4) + tp5 - tp6) * 10000000.0;
	if(tp8<0.0){
		tp8 = 0.0;
	}
	
	dopt = tp8/10000000.0 + tp7 * tp6;
	dopt = (dopt>0.0)?dopt:0.0;
	dopt = (dopt<1.0)?dopt:1.0;
	
	/* in Fig. 5a */
	max_ch4oxy = 53.8 * dopt + 0.58;
	
	/* Fig.5a */
	wmin = 3.0 * fc - 0.28;   /* Wmin */
	wopt = 6.3 * fc - 0.58;   /* Wopt */
	wmax = 10.6 * fc + 1.9;   /* Wmax */
	
	if(loct->m_vmc[grid->m]*10.0<wmin || loct->m_vmc[grid->m]*10.0>wmax){
		f_swater = 0.1;
	}else{
		
		dd = (10.0 * loct->m_vmc[grid->m] - wmin)/(wopt - wmin);
		ee = (10.0 * loct->m_vmc[grid->m] - wmax)/(wopt - wmax);
		ff = 0.4 * (wmax - wopt)/(wopt - wmin);
		f_swater = pow(dd, 0.4) * pow(ee, ff);
	}
	
	/* Fig.5b */
	aa = (dopt>0.11)?dopt:0.11;
	f_tsoil = grid->tmp10_soil[grid->m] * aa * 0.095 + 0.9;
	
	/* natural ecosystems */
	ch4ox_nat = max_ch4oxy * 1.0  * f_tsoil * f_swater;

	/* croplands */
	/* Fig. 5c */
	if(dopt<0.1){
		f_dopt = 0.9;
	}else if(dopt>0.28){
		f_dopt = 0.28;
	}else{
		f_dopt = -4.6 * dopt + 1.6;
	}
	ch4ox_agr = max_ch4oxy * f_dopt  * f_tsoil * f_swater;
	
	/* deciduous forests ****************************/
	if(grid->veg_sage==2||grid->veg_sage==5||grid->veg_sage==7){
		
		/*  bulk_dens = 1.3;  */
		max_ch4oxy = 40.0 - 18.3 * grid->bulkdens;
		
		aa = (0.1 * loct->wfps[grid->m]*100.0 - 0.5)/(1.84 - 0.5);
		bb = (0.1 * loct->wfps[grid->m]*100.0 - 55.0)/(1.84 - 55.0);
		cc = 0.13*(55.0 - 1.84)/(1.84 - 0.5);
		
		f_wfps = pow(aa, 0.13) * pow(bb, cc);
		if(f_wfps < 0.1){
			f_wfps = 0.1;
		}
		
		f_tsoil = 0.0209 * grid->tmp10_soil[grid->m] + 0.845;
		if(f_tsoil < 0.0){
			f_tsoil = 0.0;
		}
		
		/* g C ha-1 day-1 */
		ch4ox_nat = max_ch4oxy * f_wfps * f_tsoil;
	}

	ch4ox_com = (1.0 - grid->f_crop_con)*ch4ox_nat + grid->f_crop_con*ch4ox_agr;

	/* mg CH4 m-2 day-1 */
	(flux->soil).ch4oxy_delgrosso[grid->m] = ch4ox_com * 1000.0/10000.0 * 16.0/12.0;
	/* mg CH4 m-2 month-1 */
	(flux->soil).ch4oxy_delgrosso[grid->m] *= (double)grid->mm[grid->m];
	
	(flux->soil).ch4oxy_delgrosso[grid->m] *= grid->f_upland;

	/*  if((flux->soil).ch4oxy_delgrosso[grid->m]>1000.0 || (flux->soil).ch4oxy_delgrosso[grid->m]<0.0){
		(flux->soil).ch4oxy_delgrosso[grid->m] = 0.0;
	}  */
}

/* upland CH4 oxidation by Curry (2007) */
/*
	Curry, C. L. (2007). "Modeling the soil consumption of atmospheric methane 
	at the global scale." Global Biogeochemical Cycles 21(GB4012): 10.1029/2006GB002818.
*/
void f_ch4oxy_curry(
	struct Grid *grid, 
	struct Loct *loct,  
	struct Flux *flux
){
	(flux->soil).ch4oxy_curry[grid->m] = 0.0;
}
