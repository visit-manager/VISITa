/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/* make a result-output file and show values in the console port */
#include<stdio.h>
#include<stdlib.h>
#include"structure.h"
#include"prototype.h"

/*********** make output file for stable state *************/
void publish_cbud(
	struct Grid *grid, 
	struct Loct *loct, 
	struct Echar *echar, 
	struct Mass *mass, 
	struct Flux *flux, 
	FILE *result
){
	long f;
	
	fprintf(result,"%.3lf %.3lf %ld %ld ",grid->lat, grid->lon, grid->time, grid->nnn);	
	/* fprintf(result,"%lf %lf %lf %lf ", grid->whc30, grid->whc, grid->sd, grid->hyd_cond); */
	fprintf(result,"\n");
	
	for(f=0;f<12;f++){
		fprintf(result,"%.2lf ", grid->bCO2[f]);	
		fprintf(result,"%lf %lf ", loct->C3ptn[f],loct->C4ptn[f]);
		
		fprintf(result,"%.2lf %.2lf ", grid->tmp_2m[f], grid->prate_sfc[f]);
		fprintf(result,"%.2lf %.2lf %.2lf %.2lf ", loct->evpr[f],loct->incep[f],loct->trspr[f],loct->ro2[f]);
		fprintf(result,"%.2lf %.2lf %.2lf ", loct->msw30[f],loct->msww[f],loct->msnwa[f]);
		
		fprintf(result,"%ld %ld ",(echar->c3).season[f],(echar->c4).season[f]);
		fprintf(result,"%.2lf %.2lf %.2lf %.2lf ",(mass->c3).lai[f],(mass->c3).mfol[f],(mass->c3).mstm[f],(mass->c3).mrot[f]);
		fprintf(result,"%.2lf %.2lf %.2lf %.2lf ",(mass->c4).lai[f],(mass->c4).mfol[f],(mass->c4).mstm[f],(mass->c4).mrot[f]);
		fprintf(result,"%.2lf %.2lf ",(mass->soil).ltr_m[f],(mass->soil).msl_m[f]);
		
		fprintf(result,"%.2lf %.2lf %.2lf %.2lf ",(flux->c3).gpp[f],(flux->c3).spp[f],(flux->c3).epp[f],(flux->c3).npp[f]);
		fprintf(result,"%.2lf %.2lf %.2lf %.2lf ",(flux->c4).gpp[f],(flux->c4).spp[f],(flux->c4).epp[f],(flux->c4).npp[f]);
		
		fprintf(result,"%.2lf %.2lf %.2lf ",(flux->c3).lf[f],(flux->c3).lc[f],(flux->c3).lr[f]);
		fprintf(result,"%.2lf %.2lf %.2lf ",(flux->c4).lf[f],(flux->c4).lc[f],(flux->c4).lr[f]);

		fprintf(result,"%.2lf %.2lf %.2lf %.2lf ",(flux->soil).rl[f],(flux->soil).rh[f],(flux->soil).sf[f],(flux->soil).rS[f]);
		fprintf(result,"%.2lf %.2lf ",flux->nep[f],flux->ncb[f]);
		
		
	/*	fprintf(result,"%lf %lf %lf %lf ",grid->tmp_sfc[f],grid->tmp_2m[f],grid->tcdc_clm[f],grid->wnd_10m[f]);
		fprintf(result,"%lf %lf %lf %lf %lf ",grid->prate_sfc[f],grid->spfh_2m[f],loct->vp[f],loct->vpd[f],loct->slope_vps[f]);
		fprintf(result,"%lf %lf %lf ",grid->dlen[f],grid->gl_rad[f],loct->gl_rad_g[f]);
		fprintf(result,"%lf %lf %lf ",loct->rad_net_g[f],loct->rad_net_p[f],loct->rad_net_long[f]);
		fprintf(result,"%.2lf %.2lf %.2lf %.2lf ", grid->tmp_2m[f], grid->prate_sfc[f], grid->par[f], loct->rad_net_g[f]+loct->rad_net_p[f]);
		fprintf(result,"%lf %lf %lf %lf %lf ",grid->par[f],grid->par_be[f],grid->par_de[f],grid->par_bp[f],grid->par_dp[f]);
		fprintf(result,"%lf %lf %lf %lf ",(echar->c3).albedo,(echar->c4).albedo,(echar->soil).albedo[f],loct->albedo_sfc[f]);	 
		fprintf(result,"%.2lf %.2lf %.2lf ",loct->canopy_con[f],(echar->c3).gs[f],(echar->c4).gs[f]);
		fprintf(result,"%.2lf %.2lf %.2lf ",loct->pm_evp[f],loct->pm_incep[f],loct->pm_trn[f]);
		fprintf(result,"%.2lf %.2lf %.2lf %.2lf ",(echar->c3).opt_lai[f],(echar->c3).psat[f],(echar->c3).eK[f],(echar->c3).lue[f]);
		fprintf(result,"%.2lf %.2lf %.2lf %.2lf %.2lf ",(echar->c3).ci[f],(echar->c3).gs[f],(echar->c3).gc[f],(echar->c3).qTf[f],(echar->c3).lf[f]);
		fprintf(result,"%.2lf %.2lf %.2lf %.2lf ",(echar->c4).opt_lai[f],(echar->c4).psat[f],(echar->c4).eK[f],(echar->c4).lue[f]);
		fprintf(result,"%.2lf %.2lf %.2lf %.2lf %.2lf ",(echar->c4).ci[f],(echar->c4).gs[f],(echar->c4).gc[f],(echar->c4).qTf[f],(echar->c4).lf[f]);
		fprintf(result,"%lf %lf ",(echar->c3).ft[f],(echar->c4).ft[f]);
		fprintf(result,"%lf %lf ",(echar->c3).fcd[f],(echar->c4).fcd[f]);
		fprintf(result,"%lf %lf ",(echar->c3).fsw[f],(echar->c4).fsw[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c3).rfm[f],(flux->c3).rcm[f],(flux->c3).rrm[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c3).rfg[f],(flux->c3).rcg[f],(flux->c3).rrg[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c3).rpm[f],(flux->c3).rpg[f],(flux->c3).rp[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c3).tpf[f],(flux->c3).tpc[f],(flux->c3).tpr[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c4).rfm[f],(flux->c4).rcm[f],(flux->c4).rrm[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c4).rfg[f],(flux->c4).rcg[f],(flux->c4).rrg[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c4).rpm[f],(flux->c4).rpg[f],(flux->c4).rp[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c4).tpf[f],(flux->c4).tpc[f],(flux->c4).tpr[f]); */

		fprintf(result,"\n");	
	}
	fprintf(result,"%lf %lf\n", flux->npp_miami, flux->npp_montreal); 
	fprintf(result,"%.2lf %.2lf %.2lf ", flux->erod_soil, flux->erod_orgmat, flux->erod_carbon);
	
/*	for(f=0;f<12;f++){
		fprintf(result,"%lf %lf ",grid->d13C_bCO2[f],loct->d13C_aCO2[f]);

		fprintf(result,"%lf %lf ",(echar->c3).photo_13c_frac[f], (flux->c3).d13c_gpp[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c3).d13c_rfm[f], (flux->c3).d13c_rcm[f], (flux->c3).d13c_rrm[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c3).d13c_rfg[f], (flux->c3).d13c_rcg[f], (flux->c3).d13c_rrg[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c3).d13c_tpf[f], (flux->c3).d13c_tpc[f], (flux->c3).d13c_tpr[f]);
		fprintf(result,"%lf %lf %lf ",(mass->c3).d13c_mfol[f], (mass->c3).d13c_mstm[f], (mass->c3).d13c_mrot[f]);

		fprintf(result,"%lf %lf ",(echar->c4).photo_13c_frac[f], (flux->c4).d13c_gpp[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c4).d13c_rfm[f], (flux->c4).d13c_rcm[f], (flux->c4).d13c_rrm[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c4).d13c_rfg[f], (flux->c4).d13c_rcg[f], (flux->c4).d13c_rrg[f]);
		fprintf(result,"%lf %lf %lf ",(flux->c4).d13c_tpf[f], (flux->c4).d13c_tpc[f], (flux->c4).d13c_tpr[f]);
		fprintf(result,"%lf %lf %lf ",(mass->c4).d13c_mfol[f], (mass->c4).d13c_mstm[f], (mass->c4).d13c_mrot[f]);

		fprintf(result,"%lf %lf %lf ",(flux->plant).d13c_lL[f],(flux->soil).d13c_lL[f],(flux->soil).d13c_sf[f]);

		fprintf(result,"%lf %lf ",(flux->soil).d13c_rl[f], (flux->soil).d13c_rh[f]);
		fprintf(result,"%lf %lf ",(mass->soil).d13c_mltr[f], (mass->soil).d13c_mmsl[f]);

		fprintf(result,"%lf %lf ",(flux->c3).d13c_rp[f], (flux->c4).d13c_rp[f]);
		fprintf(result,"%lf %lf ",(flux->c3).d13c_lL[f], (flux->c4).d13c_lL[f]);
		fprintf(result,"%lf ",(flux->soil).d13c_rS[f]);

		fprintf(result," \n");
	} */
	
	fprintf(result," \n"); 	
}
