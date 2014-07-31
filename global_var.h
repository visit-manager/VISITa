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

double MDN[ASTEP] = {31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0};

/* experimental variables ************/
long GCM, CO2S, GCM_R, GCM_C; /* */
long PARAM_PTB; /* added by A.Ito (2010/05/10) */
long EX_CH4_1, EX_CH4_2, EX_CH4_3; /* added by A.Ito (2010/07/02) */
long EX_SRM;
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
double f_pert[20];

/* atm. GHG scenario */
double aco2_a1[DL_AGHG], aco2_a2[DL_AGHG], aco2_b1[DL_AGHG], aco2_b2[DL_AGHG];
double ach4_a1[DL_AGHG], ach4_a2[DL_AGHG], ach4_b1[DL_AGHG], ach4_b2[DL_AGHG];
double an2o_a1[DL_AGHG], an2o_a2[DL_AGHG], an2o_b1[DL_AGHG], an2o_b2[DL_AGHG];
/* for ISI-MIP: 2012/06/27 by A.Ito */
/* a1 => RCP6.0 */
/* a2 => RCP8.5 */
/* b1 => RCP2.6 */
/* b2 => RCP4.5 */

/* historical results: 201 years, spinup+1901-2100 *************/
double h_tmp[PD_SIM], h_pre[PD_SIM], h_dswr[PD_SIM], h_aet[PD_SIM], h_rof[PD_SIM];
double h_gpp[PD_SIM], h_npp[PD_SIM], h_nep[PD_SIM], h_plant[PD_SIM], h_soil[PD_SIM];
double h_sr[PD_SIM], h_ersn_c[PD_SIM], h_agrersn_c[PD_SIM], h_doc[PD_SIM];
double h_agrarea[PD_SIM], h_paddyarea[PD_SIM], h_luc[PD_SIM];
double h_luc_1[PD_SIM], h_luc_2[PD_SIM], h_luc_3[PD_SIM];
double h_gpp_df97[PD_SIM];
double h_gpp_c4[PD_SIM];	/* added by A.Ito (2009/08/31) */
double h_pot_prmfrst[PD_SIM];
double h_trnsp[PD_SIM], h_incepev[PD_SIM], h_ssurfev[PD_SIM];
double h_nbp[PD_SIM], h_hvst[PD_SIM], h_abgm[PD_SIM];
double h_sw1[PD_SIM], h_sw2[PD_SIM];
double h_rns[PD_SIM], h_rnl[PD_SIM];	/* added by A.Ito (2013/01/02) */
double h_rnsd[PD_SIM], h_cld[PD_SIM], h_apar[PD_SIM];
double h_parb[PD_SIM], h_pard[PD_SIM];
double h_arm[PD_SIM];

double h_burnt_area[PD_SIM];
double h_bioburn_co2[PD_SIM], h_bioburn_ch4[PD_SIM], h_bioburn_co[PD_SIM];
double h_bioburn_nmhc[PD_SIM], h_bioburn_oc[PD_SIM], h_bioburn_bc[PD_SIM];
double h_bioburn_nox[PD_SIM], h_bioburn_so2[PD_SIM], h_bioburn_pm25[PD_SIM];
double h_bioburn_tpm[PD_SIM], h_bioburn_tec[PD_SIM];

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
double h_n_fertin[PD_SIM], h_n_depoin[PD_SIM]; /* added by A.Ito (2010/05/02) */

double h_voc_isopr_g97[PD_SIM], h_voc_monotrp_g97[PD_SIM], h_voc_methanl_g97[PD_SIM];
double h_voc_acetone_g97[PD_SIM], h_voc_actaldhd_g97[PD_SIM], h_voc_frmardhd_g97[PD_SIM];
double h_voc_formacd_g97[PD_SIM], h_voc_acetacd_g97[PD_SIM], h_voc_co_g97[PD_SIM];

double h_hvst_wood[PD_SIM], h_wetarea[PD_SIM];

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
double rh_hvst[N_REG][PD_SIM], rh_luc[N_REG][PD_SIM];
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

/* 0: 1950s */
/* 1: 1990s */
/* 2: 2020s */
/* 3: 2050s */
/* 4: 2080s */
float g_tmp[5][N_ROW][N_COL];
float g_prc[5][N_ROW][N_COL];
float g_swr[5][N_ROW][N_COL];
float g_gpp[5][N_ROW][N_COL];
float g_npp[5][N_ROW][N_COL];
float g_nep[5][N_ROW][N_COL];
float g_pmas[5][N_ROW][N_COL];
float g_smas[5][N_ROW][N_COL];
float g_ch4e_cao[5][N_ROW][N_COL];
float g_ch4ep_cao[5][N_ROW][N_COL];
float g_ch4o_curry[5][N_ROW][N_COL];
float g_n2oe[5][N_ROW][N_COL];
float g_bbco2[5][N_ROW][N_COL];
float g_ersn[5][N_ROW][N_COL];
float g_isopr[5][N_ROW][N_COL];
float g_sr[5][N_ROW][N_COL];
float g_luc[5][N_ROW][N_COL]; 
float g_er[5][N_ROW][N_COL]; 
float g_snh4[5][N_ROW][N_COL]; 
float g_sno3[5][N_ROW][N_COL]; 

#if C13_GOUT==1
float g_f13[5][N_ROW][N_COL]; 
float g_c13[5][N_ROW][N_COL]; 
float g_r13[5][N_ROW][N_COL]; 
float g_l13[5][N_ROW][N_COL]; 
float g_h13[5][N_ROW][N_COL]; 
float g_gpp13[5][N_ROW][N_COL]; 
float g_er13[5][N_ROW][N_COL]; 
#endif

#if C14_GOUT==1
float g_f14[5][N_ROW][N_COL]; 
float g_c14[5][N_ROW][N_COL]; 
float g_r14[5][N_ROW][N_COL]; 
float g_l14[5][N_ROW][N_COL]; 
float g_h14[5][N_ROW][N_COL]; 
float g_gpp14[5][N_ROW][N_COL]; 
float g_er14[5][N_ROW][N_COL]; 
#endif

#if PHYS_GOUT==1
float g_lai[5][N_ROW][N_COL]; 
float g_parb[5][N_ROW][N_COL]; 
float g_pard[5][N_ROW][N_COL]; 
float g_apar[5][N_ROW][N_COL]; 
float g_apar2[5][N_ROW][N_COL]; 
float g_aet[5][N_ROW][N_COL]; 
float g_rof[5][N_ROW][N_COL]; 
float g_rns[5][N_ROW][N_COL];
float g_rnl[5][N_ROW][N_COL]; 
float g_sw1[5][N_ROW][N_COL]; 
float g_sw2[5][N_ROW][N_COL]; 
float g_rnsd[5][N_ROW][N_COL];
#endif

#if CH4_WH==1
float g_ch4ep_wh[5][N_ROW][N_COL];
float g_ch4ew_wh[5][N_ROW][N_COL];
float gm_ch4ep_wh[12][N_ROW][N_COL];
#endif

double glat_area[N_ROW];
double glat_gpp[ASTEP][N_ROW],glat_npp[ASTEP][N_ROW],glat_nep[ASTEP][N_ROW];
double glat_ch4_cao[ASTEP][N_ROW], glat_ch4_wh[ASTEP][N_ROW];

/* float gs_gpp_1[12][N_ROW][N_COL], gs_gpp_2[12][N_ROW][N_COL], gs_gpp_3[12][N_ROW][N_COL];
float gs_ipar_1[12][N_ROW][N_COL], gs_ipar_2[12][N_ROW][N_COL], gs_ipar_3[12][N_ROW][N_COL];
float gs_ipard_1[12][N_ROW][N_COL], gs_ipard_2[12][N_ROW][N_COL], gs_ipard_3[12][N_ROW][N_COL]; */

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
