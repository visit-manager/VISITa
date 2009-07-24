/*	VISIT: Vegetation Integrative SImulator for Trace gases				*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
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

struct Loct loct_nat;
struct Echar echar_nat;
struct Mass mass_nat;
struct Flux flux_nat;

/* global variables ********************************************/
/* GCM variables ************/
long GCM, CO2S, GCM_R, GCM_C; /* */

/* atm. GHG scenario */
double aco2_a1[553], aco2_a2[553], aco2_b1[553], aco2_b2[553];
double ach4_a1[553], ach4_a2[553], ach4_b1[553], ach4_b2[553];
double an2o_a1[553], an2o_a2[553], an2o_b1[553], an2o_b2[553];

/* historical results: 201 years, spinup+1901-2100 *************/
double h_tmp[HIST], h_pre[HIST], h_dswr[HIST], h_aet[HIST], h_rof[HIST];
double h_gpp[HIST], h_npp[HIST], h_nep[HIST], h_plant[HIST], h_soil[HIST];
double h_sr[HIST], h_ersn_c[HIST], h_agrersn_c[HIST], h_doc[HIST];
double h_agrarea[HIST], h_luc[HIST];
double h_luc_1[HIST], h_luc_2[HIST], h_luc_3[HIST];
double h_gpp_df97[HIST];

double h_burnt_area[HIST];
double h_bioburn_co2[HIST], h_bioburn_ch4[HIST], h_bioburn_co[HIST];
double h_bioburn_nmhc[HIST], h_bioburn_oc[HIST], h_bioburn_bc[HIST];
double h_bioburn_nox[HIST], h_bioburn_so2[HIST], h_bioburn_pm25[HIST];
double h_bioburn_tpm[HIST], h_bioburn_tec[HIST];

double h_ch4ox1[HIST], h_ch4ox2[HIST], h_ch4ox3[HIST], h_ch4ox4[HIST];
double h_ch4emit_cao_paddy[HIST], h_ch4emit_cao_wetland[HIST];
double h_n2o_emit_ngas[HIST], h_n2_emit_ngas[HIST];
double h_n2o_emit_casa[HIST], h_no_emit_casa[HIST], h_n2_emit_casa[HIST];
double h_nh3_emit[HIST], h_n2_biofix[HIST];
double h_ch4_emit_mass[HIST], h_ch4_emit_photo[HIST];
double h_n2o_d_emit_ngas[HIST], h_n2o_n_emit_ngas[HIST];
double h_n2o_emit_ngas_agr[HIST], h_n2o_emit_casa_agr[HIST];
double h_nh3_emit_agr[HIST];

double h_voc_isopr_g97[HIST], h_voc_monotrp_g97[HIST], h_voc_methanl_g97[HIST];
double h_voc_acetone_g97[HIST], h_voc_actaldhd_g97[HIST], h_voc_frmardhd_g97[HIST];
double h_voc_formacd_g97[HIST], h_voc_acetacd_g97[HIST], h_voc_co_g97[HIST];

double ci_aco2[HIST], ci_aco2_d13c[HIST], ci_aco2_d14c[HIST];
double ci_gpp[HIST], ci_gpp_d13c[HIST], ci_gpp_d14c[HIST];
double ci_er[HIST], ci_er_d13c[HIST], ci_er_d14c[HIST];
double ci_f[HIST], ci_f_d13c[HIST], ci_f_d14c[HIST];
double ci_c[HIST], ci_c_d13c[HIST], ci_c_d14c[HIST];
double ci_r[HIST], ci_r_d13c[HIST], ci_r_d14c[HIST];
double ci_l[HIST], ci_l_d13c[HIST], ci_l_d14c[HIST];
double ci_h[HIST], ci_h_d13c[HIST], ci_h_d14c[HIST];

/* regional historical */
double rh_area[NREG];
double rh_temp[NREG][HIST], rh_prec[NREG][HIST], rh_dswrf[NREG][HIST];
double rh_rns[NREG][HIST], rh_rnl[NREG][HIST];
double rh_ipar[NREG][HIST], rh_apar[NREG][HIST];
double rh_gpp[NREG][HIST], rh_npp[NREG][HIST], rh_nep[NREG][HIST];
double rh_evpr[NREG][HIST], rh_trsp[NREG][HIST], rh_incp[NREG][HIST], rh_rnof[NREG][HIST];

/* monthly mean results **********/
double m_ch4ox1[ASTEP], m_ch4ox2[ASTEP], m_ch4ox3[ASTEP];
double m_bioburn_co2[ASTEP], m_bioburn_ch4[ASTEP], m_bioburn_co[ASTEP];
double m_bioburn_nmhc[ASTEP], m_bioburn_oc[ASTEP], m_bioburn_bc[ASTEP];
double m_gpp[ASTEP], m_npp[ASTEP], m_nep[ASTEP];

/* vegetation (olson) mean results */
double go_landarea, gs_landarea;
double vo_area[34];
double vo_gpp[34], vo_npp[34], vo_nep[34];
double vo_lai[34], vo_fol[34], vo_stm[34], vo_rot[34], vo_ltr[34], vo_msl[34];
double vs_area[16];
double vs_gpp[16], vs_npp[16], vs_nep[16];
double vs_lai[16], vs_fol[16], vs_stm[16], vs_rot[16], vs_ltr[16], vs_msl[16];

double MDN[ASTEP] = {31.0, 28.0, 31.0, 30.0, 31.0, 30.0, 31.0, 31.0, 30.0, 31.0, 30.0, 31.0};

/* 0: 1950s */
/* 1: 1990s */
/* 2: 2020s */
/* 3: 2050s */
/* 4: 2080s */
float g_tmp[5][360][720];
float g_prc[5][360][720];
float g_swr[5][360][720];
float g_gpp[5][360][720];
float g_npp[5][360][720];
float g_nep[5][360][720];
float g_pmas[5][360][720];
float g_smas[5][360][720];
float g_ch4e[5][360][720];
float g_ch4o[5][360][720];
float g_n2oe[5][360][720];
float g_bbco2[5][360][720];
float g_ersn[5][360][720];
float g_isopr[5][360][720];
float g_sr[5][360][720];
float g_luc[5][360][720]; 

float g_f13[5][360][720]; 
float g_c13[5][360][720]; 
float g_r13[5][360][720]; 
float g_l13[5][360][720]; 
float g_h13[5][360][720]; 

float g_f14[5][360][720]; 
float g_c14[5][360][720]; 
float g_r14[5][360][720]; 
float g_l14[5][360][720]; 
float g_h14[5][360][720]; 

float g_er[5][360][720]; 
float g_gpp13[5][360][720]; 
float g_er13[5][360][720]; 
float g_gpp14[5][360][720]; 
float g_er14[5][360][720]; 

float g_snh4[5][360][720]; 
float g_sno3[5][360][720]; 

/* float gs_gpp_1[12][360][720], gs_gpp_2[12][360][720], gs_gpp_3[12][360][720];
float gs_ipar_1[12][360][720], gs_ipar_2[12][360][720], gs_ipar_3[12][360][720];
float gs_ipard_1[12][360][720], gs_ipard_2[12][360][720], gs_ipard_3[12][360][720]; */

short RAD_SENS;
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
