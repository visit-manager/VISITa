/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in June 16, 2009							*/

/* definition of global variables ***************************************/

#include"setting.h"

/* structure *****************************************************/
struct Grid grid;

struct Loct loct;
struct Echar echar;
struct Mass mass;
struct Flux flux;

struct Loct loct_agr;
struct Echar echar_agr;
struct Mass mass_agr;
struct Flux flux_agr;

/* global variables ********************************************/

double MDN[ASTEP] = {31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0}, YDN=365.0;

/* experimental variables ************/
long SCENARIO_ID, CO2S, GCM_R, GCM_C; /* */
long PARAM_PTB, PARAM_ENS; /* added by A.Ito (2010/05/10) */
long EX_CH4_1, EX_CH4_2, EX_CH4_3; /* added by A.Ito (2010/07/02) */
long EX_SRM;
long IMP_DT, IMP_DP, IMP_CD;
/* experiment: geoengeneering, solar radiation management *********/
/* EX SRM: 2013/06/04 by A.Ito *******************/
/* 0: off */
/* 1: solar constant -2.6 W m-2 */
/* 2: solar constant -4.5 W m-2 */
/* 3: solar constant -6.0 W m-2 */
/* 4: solar constant -8.5 W m-2 */

/* 11: surface radiation -2.6 W m-2 */
/* 12: surface radiation -4.5 W m-2 */
/* 13: surface radiation -6.0 W m-2 */
/* 14: surface radiation -8.5 W m-2 */

/* perturbation term */
double f_pert[N_PARA_ENS];

/* atm. GHG scenario */
double aco2_1[DL_AGHG], aco2_2[DL_AGHG], aco2_3[DL_AGHG], aco2_4[DL_AGHG];
double ach4_a1[DL_AGHG], ach4_a2[DL_AGHG], ach4_b1[DL_AGHG], ach4_b2[DL_AGHG];
double an2o_a1[DL_AGHG], an2o_a2[DL_AGHG], an2o_b1[DL_AGHG], an2o_b2[DL_AGHG];
/* for ISI-MIP: 2012/06/27 by A.Ito */
/* a1 => RCP6.0 */
/* a2 => RCP8.5 */
/* b1 => RCP2.6 */
/* b2 => RCP4.5 */
/* Atmospheric d13C, D14C by Graven: 2019/1/17 by A.Ito */
double d13c_graven[166], d14c1_graven[166], d14c2_graven[166], d14c3_graven[166];

/* historical results: 201 years, spinup+1901-2100 *************/
double h_tmp[PD_SIM], h_pre[PD_SIM], h_dswr[PD_SIM], h_aet[PD_SIM], h_rof[PD_SIM];
double h_gpp[PD_SIM], h_npp[PD_SIM], h_nep[PD_SIM], h_plant[PD_SIM], h_soil[PD_SIM];
double h_sr[PD_SIM], h_ersn_c[PD_SIM], h_agrersn_c[PD_SIM], h_doc[PD_SIM];
double h_agrarea[PD_SIM], h_paddyarea[PD_SIM], h_luc[PD_SIM];
double h_luc_0[PD_SIM], h_luc_1[PD_SIM], h_luc_2[PD_SIM], h_luc_3[PD_SIM]; /* added by A.Ito (2018/10/23) */
double h_gpp_df97[PD_SIM];
double h_gpp_c4[PD_SIM];	/* added by A.Ito (2009/08/31) */
double h_pot_prmfrst[PD_SIM];
double h_trnsp[PD_SIM], h_incepev[PD_SIM], h_ssurfev[PD_SIM];
double h_nbp[PD_SIM], h_net_crop[PD_SIM], h_hvst_crop[PD_SIM], h_abgm[PD_SIM];
double h_sw1[PD_SIM], h_sw2[PD_SIM];
double h_rns[PD_SIM], h_rnl[PD_SIM];	/* added by A.Ito (2013/01/02) */
double h_rnsd[PD_SIM], h_cld[PD_SIM], h_apar[PD_SIM], h_ipar[PD_SIM];
double h_parb[PD_SIM], h_pard[PD_SIM];
double h_arm[PD_SIM], h_bco2[PD_SIM]; /* added by A.Ito: 2018/10/22 */
double h_lL[PD_SIM];

/* seasonal-cycle amplitude: 2019/03/02 by A.Ito */
double hm_sca_gpp_nh[PD_SIM][12], hm_sca_re_nh[PD_SIM][12], hm_sca_nep_nh[PD_SIM][12];

double h_burnt_area[PD_SIM];
double h_burnt_area_wood[PD_SIM];
double h_bioburn_co2[PD_SIM], h_bioburn_ch4[PD_SIM], h_bioburn_co[PD_SIM];
double h_bioburn_nmhc[PD_SIM], h_bioburn_oc[PD_SIM], h_bioburn_bc[PD_SIM];
double h_bioburn_nox[PD_SIM], h_bioburn_so2[PD_SIM], h_bioburn_pm25[PD_SIM];
double h_bioburn_tpm[PD_SIM], h_bioburn_tec[PD_SIM], h_bioburn_n2o[PD_SIM];

double h_ch4ox1[PD_SIM], h_ch4ox2[PD_SIM], h_ch4ox3[PD_SIM], h_ch4ox4[PD_SIM];
double h_ch4emit_cao_paddy[PD_SIM], h_ch4emit_cao_wetland[PD_SIM];
double h_ch4emit_paddy_wh_diff[PD_SIM], h_ch4emit_paddy_wh_plant[PD_SIM];
double h_ch4emit_paddy_wh_ebbl[PD_SIM], h_ch4emit_paddy_wh_release[PD_SIM];
double h_ch4emit_wetland_wh_diff[PD_SIM], h_ch4emit_wetland_wh_plant[PD_SIM];
double h_ch4emit_wetland_wh_ebbl[PD_SIM], h_ch4emit_wetland_wh_release[PD_SIM];

double h_n2o_emit_ngas[PD_SIM], h_n2_emit_ngas[PD_SIM];
double h_n2o_emit_casa[PD_SIM], h_no_emit_casa[PD_SIM], h_n2_emit_casa[PD_SIM];
double h_nh3_emit[PD_SIM], h_n2_biofix[PD_SIM];
double h_ch4_emit_mass[PD_SIM], h_ch4_emit_photo[PD_SIM];
double h_n2o_d_emit_ngas[PD_SIM], h_n2o_n_emit_ngas[PD_SIM];
double h_n2o_emit_ngas_agr[PD_SIM], h_n2o_emit_casa_agr[PD_SIM];
double h_nh3_emit_agr[PD_SIM];
double h_no3_leach[PD_SIM];
double h_n_fertin[PD_SIM], h_n_manurein[PD_SIM], h_n_depoin[PD_SIM]; /* added by A.Ito (2010/05/02) */
double h_n_mcrb[PD_SIM], h_n_no3[PD_SIM],h_n_nh4[PD_SIM]; /* 2016/06/23 by A.Ito */
double h_n_cnpy[PD_SIM], h_n_strg[PD_SIM], h_n_lttr[PD_SIM], h_n_hums[PD_SIM];
double h_n_immbl[PD_SIM], h_n_lmnrl[PD_SIM],h_n_hmnrl[PD_SIM],h_n_cabdn[PD_SIM];
double h_n_sabdn[PD_SIM], h_n_uptk[PD_SIM];

double h_voc_isopr_g97[PD_SIM], h_voc_monotrp_g97[PD_SIM], h_voc_methanl_g97[PD_SIM];
double h_voc_acetone_g97[PD_SIM], h_voc_actaldhd_g97[PD_SIM], h_voc_frmardhd_g97[PD_SIM];
double h_voc_formacd_g97[PD_SIM], h_voc_acetacd_g97[PD_SIM], h_voc_co_g97[PD_SIM];
double h_voc_afarnesene[PD_SIM], h_voc_bcaryophyllene[PD_SIM], h_voc_othersesqui[PD_SIM];

double h_hvst_wood[PD_SIM], h_wetarea[PD_SIM], h_deforest[PD_SIM];

/* Tropical-Extratropical (Schimel et al. 2015): 2019/03/01 by A.Ito */
double h_gpp_trp[PD_SIM], h_npp_trp[PD_SIM], h_nep_trp[PD_SIM], h_nbp_trp[PD_SIM];
double h_luc_trp[PD_SIM], h_bb_trp[PD_SIM];

double ci_aco2[PD_SIM], ci_aco2_d13c[PD_SIM], ci_aco2_d14c[PD_SIM];
double ci_gpp[PD_SIM], ci_gpp_d13c[PD_SIM], ci_gpp_d14c[PD_SIM];
double ci_er[PD_SIM], ci_er_d13c[PD_SIM], ci_er_d14c[PD_SIM];
double ci_f[PD_SIM], ci_f_d13c[PD_SIM], ci_f_d14c[PD_SIM];
double ci_c[PD_SIM], ci_c_d13c[PD_SIM], ci_c_d14c[PD_SIM];
double ci_r[PD_SIM], ci_r_d13c[PD_SIM], ci_r_d14c[PD_SIM];
double ci_l[PD_SIM], ci_l_d13c[PD_SIM], ci_l_d14c[PD_SIM];
double ci_h[PD_SIM], ci_h_d13c[PD_SIM], ci_h_d14c[PD_SIM];

double hm_temp[PD_SIM][ASTEP];
double hm_prec[PD_SIM][ASTEP];
double hm_ch4_wh[PD_SIM][ASTEP];
double hm_inund[PD_SIM][ASTEP];

/* regional historical */
double rh_area[N_REG];
double rh_temp[N_REG][PD_SIM], rh_prec[N_REG][PD_SIM], rh_dswrf[N_REG][PD_SIM];
double rh_rns[N_REG][PD_SIM], rh_rnl[N_REG][PD_SIM];
double rh_ipar[N_REG][PD_SIM], rh_apar[N_REG][PD_SIM];
double rh_gpp[N_REG][PD_SIM], rh_npp[N_REG][PD_SIM], rh_nep[N_REG][PD_SIM];
double rh_evpr[N_REG][PD_SIM], rh_trsp[N_REG][PD_SIM], rh_incp[N_REG][PD_SIM], rh_rnof[N_REG][PD_SIM];

/* added by A.Ito (2009/11/15) */
double rh_net_crop[N_REG][PD_SIM], rh_luc[N_REG][PD_SIM];
double rh_ch4ox_curry[N_REG][PD_SIM], rh_ch4emit_wh_wet[N_REG][PD_SIM], rh_ch4emit_wh_paddy[N_REG][PD_SIM];
double rh_n2o_emit_ngas[N_REG][PD_SIM], rh_n2o_emitagr_ngas[N_REG][PD_SIM];

/* added by A.Ito (2009/09/30) */
double rh_ci_gpp[N_REG][PD_SIM], rh_ci_gpp_d13c[N_REG][PD_SIM], rh_ci_gpp_d14c[N_REG][PD_SIM];
double rh_ci_er[N_REG][PD_SIM], rh_ci_er_d13c[N_REG][PD_SIM], rh_ci_er_d14c[N_REG][PD_SIM];
double rh_ci_f[N_REG][PD_SIM], rh_ci_f_d13c[N_REG][PD_SIM], rh_ci_f_d14c[N_REG][PD_SIM];
double rh_ci_c[N_REG][PD_SIM], rh_ci_c_d13c[N_REG][PD_SIM], rh_ci_c_d14c[N_REG][PD_SIM];
double rh_ci_r[N_REG][PD_SIM], rh_ci_r_d13c[N_REG][PD_SIM], rh_ci_r_d14c[N_REG][PD_SIM];
double rh_ci_l[N_REG][PD_SIM], rh_ci_l_d13c[N_REG][PD_SIM], rh_ci_l_d14c[N_REG][PD_SIM];
double rh_ci_h[N_REG][PD_SIM], rh_ci_h_d13c[N_REG][PD_SIM], rh_ci_h_d14c[N_REG][PD_SIM];

/* monthly mean results **********/
double m_ch4ox1[ASTEP], m_ch4ox2[ASTEP], m_ch4ox3[ASTEP];
double m_bioburn_co2[ASTEP], m_bioburn_ch4[ASTEP], m_bioburn_co[ASTEP];
double m_bioburn_nmhc[ASTEP], m_bioburn_oc[ASTEP], m_bioburn_bc[ASTEP];
double m_gpp[ASTEP], m_npp[ASTEP], m_nep[ASTEP];
double m_ch4p_cao[ASTEP], m_ch4p_wh[ASTEP];

/* vegetation (olson) mean results */
double go_landarea, gs_landarea;
double vo_area[NVEG_OLSON];
double vo_gpp[NVEG_OLSON], vo_npp[NVEG_OLSON], vo_nep[NVEG_OLSON];
double vo_lai[NVEG_OLSON], vo_fol[NVEG_OLSON], vo_stm[NVEG_OLSON], vo_rot[NVEG_OLSON], vo_ltr[NVEG_OLSON], vo_msl[NVEG_OLSON];
/* vegetation (SAGE) mean results */
double vs_area[NVEG_SAGE];
double vs_gpp[NVEG_SAGE], vs_npp[NVEG_SAGE], vs_nep[NVEG_SAGE];
double vs_lai[NVEG_SAGE], vs_fol[NVEG_SAGE], vs_stm[NVEG_SAGE], vs_rot[NVEG_SAGE], vs_ltr[NVEG_SAGE], vs_msl[NVEG_SAGE];

double vo_n_cnpy[NVEG_OLSON], vo_n_strg[NVEG_OLSON], vo_n_mcrb[NVEG_OLSON], vo_n_ltr[NVEG_OLSON], vo_n_hms[NVEG_OLSON];

double glat_area[N_ROW];
double glat_gpp[ASTEP][N_ROW],glat_npp[ASTEP][N_ROW],glat_nep[ASTEP][N_ROW];
double glat_ch4_cao[ASTEP][N_ROW], glat_ch4_wh[ASTEP][N_ROW];
double glat_agb[N_ROW],glat_soc[N_ROW];

/* float gs_gpp_1[12][N_ROW][N_COL], gs_gpp_2[12][N_ROW][N_COL], gs_gpp_3[12][N_ROW][N_COL];
float gs_ipar_1[12][N_ROW][N_COL], gs_ipar_2[12][N_ROW][N_COL], gs_ipar_3[12][N_ROW][N_COL];
float gs_ipard_1[12][N_ROW][N_COL], gs_ipard_2[12][N_ROW][N_COL], gs_ipard_3[12][N_ROW][N_COL]; */

/* CHASER4.0 monthly, by A.Ito (2014/11/19) */
double	ndepo_chaser4_nhx_h[ASTEP][64][128];		/* NHx */
double	ndepo_chaser4_noy_h[ASTEP][64][128];		/* NOy */
double	ndepo_chaser4_ont_h[ASTEP][64][128];		/* Org NOx */
double	ndepo_chaser4_nhx_p[ASTEP][64][128];		/* NHx */
double	ndepo_chaser4_noy_p[ASTEP][64][128];		/* NOy */
double	ndepo_chaser4_ont_p[ASTEP][64][128];		/* Org NOx */

short SENS_RAD;
/* 0: control */
/* 1: +10 diffuse PAR fraction */
/* 2: -10 diffuse PAR fraction */
/* 3: +10 soil albedo */
/* 4: -10 soil albedo */
/* 5: +10 canopy scattering coefficient  */
/* 6: -10 canopy scattering coefficient */
/* 7: +10 diffuse attenuation coefficient */
/* 8: -10 diffuse attenuation coefficient */
/* 9: +10 direct attenuation coefficient */
/* 10: -10 direct attenuation coefficient */

short DF97;
/* 0: off - Monsi-Saeki */
/* 1: on De Pury-Farquhar */

/* future solar radiation change */
short CC_R;
/* 0: off */
/* 1: on  */

/* gradual temperature change */
short TEMP_GC;
/* 0: control */
/* 1: -1 deg C per century */
/* 2: +1 deg C per century */
/* 3: -2 deg C per century */
/* 4: +2 deg C per century */
/* 5: -3 deg C per century */
/* 6: +3 deg C per century */
/* 7: -4 deg C per century */
/* 8: +4 deg C per century */
/* 9: -5 deg C per century */
/* 10: +5 deg C per century */
/* 11: -6 deg C per century */
/* 12: +6 deg C per century */

/* CHASER 2001 monthly, by A.Ito (2010/05/21) */
double	ndepo_chaser_dnhx[ASTEP][64][128];		/* NHx, dry */
double	ndepo_chaser_dnoy[ASTEP][64][128];		/* NOy, dry */
double	ndepo_chaser_wnhx[ASTEP][64][128];		/* NHx, wet */
double	ndepo_chaser_wnoy[ASTEP][64][128];		/* NOy, wet */

/* Burnt area by GFED4s: 2018/05/19 by A.Ito */
double bf_gfed4s[20][N_REG+1];

