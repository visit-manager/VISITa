/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & RIGC/JAMSTEC						*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/*  Updated August 7, 2007								*/
/*  Updated November 29, 2007							*/
/*  Updated by A.Ito (2010/05/02)						*/

#include"setting.h"

/* #define IFILEN 59 */  /* normal case */
#define IFILEN 89 /* */  /* normal case */
#define OFILEN 9

extern short DF97;
extern double MDN[ASTEP];
extern long GCM_ID, CO2S, GCM_R, GCM_C;
extern long PARAM_PTB, PARAM_ENS;   /* added by A.Ito (2010/05/10) */
extern long EX_CH4_1, EX_CH4_2, EX_CH4_3;   /* added by A.Ito (2010/07/02) */
extern long EX_SRM, EX_CCPL;

extern double f_pert[20];
extern double aco2_a1[DL_AGHG], aco2_a2[DL_AGHG], aco2_b1[DL_AGHG], aco2_b2[DL_AGHG];
extern double ach4_a1[DL_AGHG], ach4_a2[DL_AGHG], ach4_b1[DL_AGHG], ach4_b2[DL_AGHG];
extern double an2o_a1[DL_AGHG], an2o_a2[DL_AGHG], an2o_b1[DL_AGHG], an2o_b2[DL_AGHG];

extern double glandarea;
extern double h_tmp[PD_SIM], h_pre[PD_SIM], h_dswr[PD_SIM], h_aet[PD_SIM], h_rof[PD_SIM];
extern double h_gpp[PD_SIM], h_npp[PD_SIM], h_nep[PD_SIM], h_plant[PD_SIM], h_soil[PD_SIM];
extern double h_sr[PD_SIM], h_ersn_c[PD_SIM], h_agrersn_c[PD_SIM], h_doc[PD_SIM];
extern double h_agrarea[PD_SIM], h_paddyarea[PD_SIM], h_luc[PD_SIM];
extern double h_luc_1[PD_SIM], h_luc_2[PD_SIM], h_luc_3[PD_SIM];
extern double h_gpp_df97[PD_SIM], h_gpp_c4[PD_SIM]; /* added by A.Ito (2009/08/31) */
extern double h_pot_prmfrst[PD_SIM];
extern double h_trnsp[PD_SIM], h_incepev[PD_SIM], h_ssurfev[PD_SIM];
extern double h_nbp[PD_SIM], h_hvst[PD_SIM], h_abgm[PD_SIM];
extern double h_sw1[PD_SIM], h_sw2[PD_SIM];
extern double h_rns[PD_SIM], h_rnl[PD_SIM];	/* added by A.Ito (2013/01/02) */
extern double h_rnsd[PD_SIM], h_cld[PD_SIM], h_apar[PD_SIM], h_ipar[PD_SIM];
extern double h_parb[PD_SIM], h_pard[PD_SIM];
extern double h_arm[PD_SIM];

extern double h_burnt_area[PD_SIM];
extern double h_bioburn_co2[PD_SIM], h_bioburn_ch4[PD_SIM], h_bioburn_co[PD_SIM];
extern double h_bioburn_nmhc[PD_SIM], h_bioburn_oc[PD_SIM], h_bioburn_bc[PD_SIM];
extern double h_bioburn_nox[PD_SIM], h_bioburn_so2[PD_SIM], h_bioburn_pm25[PD_SIM];
extern double h_bioburn_tpm[PD_SIM], h_bioburn_tec[PD_SIM];

extern double h_ch4ox1[PD_SIM], h_ch4ox2[PD_SIM], h_ch4ox3[PD_SIM], h_ch4ox4[PD_SIM];
extern double h_ch4emit_cao_paddy[PD_SIM], h_ch4emit_cao_wetland[PD_SIM];
extern double h_ch4emit_paddy_wh_diff[PD_SIM],  h_ch4emit_paddy_wh_plant[PD_SIM];
extern double h_ch4emit_paddy_wh_ebbl[PD_SIM],  h_ch4emit_paddy_wh_release[PD_SIM];
extern double h_ch4emit_wetland_wh_diff[PD_SIM],  h_ch4emit_wetland_wh_plant[PD_SIM];
extern double h_ch4emit_wetland_wh_ebbl[PD_SIM],  h_ch4emit_wetland_wh_release[PD_SIM];

extern double h_n2o_emit_ngas[PD_SIM], h_n2_emit_ngas[PD_SIM];
extern double h_n2o_emit_casa[PD_SIM], h_no_emit_casa[PD_SIM], h_n2_emit_casa[PD_SIM];
extern double h_nh3_emit[PD_SIM], h_n2_biofix[PD_SIM];
extern double h_ch4_emit_mass[PD_SIM], h_ch4_emit_photo[PD_SIM];
extern double h_n2o_d_emit_ngas[PD_SIM], h_n2o_n_emit_ngas[PD_SIM];
extern double h_n2o_emit_ngas_agr[PD_SIM], h_n2o_emit_casa_agr[PD_SIM];
extern double h_nh3_emit_agr[PD_SIM];
extern double h_no3_leach[PD_SIM];
extern double h_n_fertin[PD_SIM], h_n_depoin[PD_SIM]; /* added by A.Ito (2010/05/02) */
extern double h_n_mcrb[PD_SIM],h_n_no3[PD_SIM],h_n_nh4[PD_SIM]; /* 2016/06/23 by A.Ito */
extern double h_n_cnpy[PD_SIM], h_n_strg[PD_SIM], h_n_lttr[PD_SIM], h_n_hums[PD_SIM];
extern double h_n_immbl[PD_SIM],h_n_lmnrl[PD_SIM],h_n_hmnrl[PD_SIM],h_n_cabdn[PD_SIM],h_n_sabdn[PD_SIM],h_n_uptk[PD_SIM];

extern double h_voc_isopr_g97[PD_SIM], h_voc_monotrp_g97[PD_SIM], h_voc_methanl_g97[PD_SIM];
extern double h_voc_acetone_g97[PD_SIM], h_voc_actaldhd_g97[PD_SIM], h_voc_frmardhd_g97[PD_SIM];
extern double h_voc_formacd_g97[PD_SIM], h_voc_acetacd_g97[PD_SIM], h_voc_co_g97[PD_SIM];
extern double h_voc_afarnesene[PD_SIM], h_voc_bcaryophyllene[PD_SIM], h_voc_othersesqui[PD_SIM];

extern double h_hvst_wood[PD_SIM], h_wetarea[PD_SIM], h_deforest[PD_SIM];

extern double ci_aco2[PD_SIM], ci_aco2_d13c[PD_SIM], ci_aco2_d14c[PD_SIM];
extern double ci_gpp[PD_SIM], ci_gpp_d13c[PD_SIM], ci_gpp_d14c[PD_SIM];
extern double ci_er[PD_SIM], ci_er_d13c[PD_SIM], ci_er_d14c[PD_SIM];
extern double ci_f[PD_SIM], ci_f_d13c[PD_SIM], ci_f_d14c[PD_SIM];
extern double ci_c[PD_SIM], ci_c_d13c[PD_SIM], ci_c_d14c[PD_SIM];
extern double ci_r[PD_SIM], ci_r_d13c[PD_SIM], ci_r_d14c[PD_SIM];
extern double ci_l[PD_SIM], ci_l_d13c[PD_SIM], ci_l_d14c[PD_SIM];
extern double ci_h[PD_SIM], ci_h_d13c[PD_SIM], ci_h_d14c[PD_SIM];

extern double hm_temp[PD_SIM][ASTEP], hm_prec[PD_SIM][ASTEP], hm_ch4_wh[PD_SIM][ASTEP], hm_inund[PD_SIM][ASTEP];

/* monthly results **********/
extern double m_ch4ox1[ASTEP], m_ch4ox2[12], m_ch4ox3[12];
extern double m_bioburn_co2[12], m_bioburn_ch4[12], m_bioburn_co[12];
extern double m_bioburn_nmhc[12], m_bioburn_oc[12], m_bioburn_bc[12];
extern double m_gpp[ASTEP], m_npp[ASTEP], m_nep[ASTEP];
extern double m_ch4p_cao[ASTEP], m_ch4p_wh[ASTEP];

/* vegetation (olson) results */
extern double go_landarea, gs_landarea;
extern double vo_area[NVEG_OLSON];
extern double vo_gpp[NVEG_OLSON], vo_npp[NVEG_OLSON], vo_nep[NVEG_OLSON];
extern double vo_lai[NVEG_OLSON], vo_fol[NVEG_OLSON], vo_stm[NVEG_OLSON], vo_rot[NVEG_OLSON], vo_ltr[NVEG_OLSON], vo_msl[NVEG_OLSON];
extern double vs_area[16];
extern double vs_gpp[16], vs_npp[16], vs_nep[16];
extern double vs_lai[16], vs_fol[16], vs_stm[16], vs_rot[16], vs_ltr[16], vs_msl[16];

extern double vo_n_cnpy[NVEG_OLSON], vo_n_strg[NVEG_OLSON], vo_n_mcrb[NVEG_OLSON], vo_n_ltr[NVEG_OLSON], vo_n_hms[NVEG_OLSON];

extern float g_tmp[5][N_ROW][N_COL];
extern float g_prc[5][N_ROW][N_COL];
extern float g_swr[5][N_ROW][N_COL];
extern float g_gpp[5][N_ROW][N_COL];
extern float g_npp[5][N_ROW][N_COL];
extern float g_nep[5][N_ROW][N_COL];
extern float g_pmas[5][N_ROW][N_COL];
extern float g_smas[5][N_ROW][N_COL];
extern float g_ch4e_cao[5][N_ROW][N_COL];
extern float g_ch4o_curry[5][N_ROW][N_COL];
extern float g_n2oe[5][N_ROW][N_COL];
extern float g_bbco2[5][N_ROW][N_COL];
extern float g_ersn[5][N_ROW][N_COL];
extern float g_isopr[5][N_ROW][N_COL];
extern float g_sr[5][N_ROW][N_COL];
extern float g_luc[5][N_ROW][N_COL];

#if C13_GOUT==1
extern float g_f13[5][N_ROW][N_COL]; 
extern float g_c13[5][N_ROW][N_COL]; 
extern float g_r13[5][N_ROW][N_COL]; 
extern float g_l13[5][N_ROW][N_COL]; 
extern float g_h13[5][N_ROW][N_COL]; 
extern float g_gpp13[5][N_ROW][N_COL]; 
extern float g_er13[5][N_ROW][N_COL]; 
#endif

#if C14_GOUT==1
extern float g_f14[5][N_ROW][N_COL]; 
extern float g_c14[5][N_ROW][N_COL]; 
extern float g_r14[5][N_ROW][N_COL]; 
extern float g_l14[5][N_ROW][N_COL]; 
extern float g_h14[5][N_ROW][N_COL]; 
extern float g_gpp14[5][N_ROW][N_COL]; 
extern float g_er14[5][N_ROW][N_COL]; 
#endif

extern float g_er[5][N_ROW][N_COL]; 
extern float g_snh4[5][N_ROW][N_COL]; 
extern float g_sno3[5][N_ROW][N_COL]; 

#if PHYS_GOUT==1
extern float g_lai[5][N_ROW][N_COL]; 
extern float g_parb[5][N_ROW][N_COL]; 
extern float g_pard[5][N_ROW][N_COL]; 
extern float g_apar[5][N_ROW][N_COL]; 
extern float g_apar2[5][N_ROW][N_COL];
extern float g_aet[5][N_ROW][N_COL]; 
extern float g_rof[5][N_ROW][N_COL]; 
extern float g_rns[5][N_ROW][N_COL]; 
extern float g_rnl[5][N_ROW][N_COL]; 
extern float g_sw1[5][N_ROW][N_COL]; 
extern float g_sw2[5][N_ROW][N_COL]; 
extern float g_rnsd[5][N_ROW][N_COL];
#endif

#if CH4_WH==1
extern float g_ch4ep_wh[5][N_ROW][N_COL]; 
extern float g_ch4ew_wh[5][N_ROW][N_COL]; 
extern float gm_ch4ep_wh[12][N_ROW][N_COL];
#endif
extern float g_ch4ep_cao[5][N_ROW][N_COL];

extern double glat_area[N_ROW];
extern double glat_gpp[ASTEP][N_ROW],glat_npp[ASTEP][N_ROW],glat_nep[ASTEP][N_ROW];
extern double glat_ch4_cao[ASTEP][N_ROW], glat_ch4_wh[ASTEP][N_ROW];

/* regional historical */
extern double rh_area[N_REG];
extern double rh_temp[N_REG][PD_SIM], rh_prec[N_REG][PD_SIM], rh_dswrf[N_REG][PD_SIM];
extern double rh_rns[N_REG][PD_SIM], rh_rnl[N_REG][PD_SIM];
extern double rh_ipar[N_REG][PD_SIM], rh_apar[N_REG][PD_SIM];
extern double rh_gpp[N_REG][PD_SIM], rh_npp[N_REG][PD_SIM], rh_nep[N_REG][PD_SIM];
extern double rh_evpr[N_REG][PD_SIM], rh_trsp[N_REG][PD_SIM], rh_incp[N_REG][PD_SIM], rh_rnof[N_REG][PD_SIM];
extern double rh_ci_gpp[N_REG][PD_SIM], rh_ci_gpp_d13c[N_REG][PD_SIM], rh_ci_gpp_d14c[N_REG][PD_SIM];
extern double rh_ci_er[N_REG][PD_SIM], rh_ci_er_d13c[N_REG][PD_SIM], rh_ci_er_d14c[N_REG][PD_SIM];
extern double rh_ci_f[N_REG][PD_SIM], rh_ci_f_d13c[N_REG][PD_SIM], rh_ci_f_d14c[N_REG][PD_SIM];
extern double rh_ci_c[N_REG][PD_SIM], rh_ci_c_d13c[N_REG][PD_SIM], rh_ci_c_d14c[N_REG][PD_SIM];
extern double rh_ci_r[N_REG][PD_SIM], rh_ci_r_d13c[N_REG][PD_SIM], rh_ci_r_d14c[N_REG][PD_SIM];
extern double rh_ci_l[N_REG][PD_SIM], rh_ci_l_d13c[N_REG][PD_SIM], rh_ci_l_d14c[N_REG][PD_SIM];
extern double rh_ci_h[N_REG][PD_SIM], rh_ci_h_d13c[N_REG][PD_SIM], rh_ci_h_d14c[N_REG][PD_SIM];

extern double rh_hvst[N_REG][PD_SIM], rh_luc[N_REG][PD_SIM];
extern double rh_ch4ox_curry[N_REG][PD_SIM], rh_ch4emit_wh_wet[N_REG][PD_SIM], rh_ch4emit_wh_paddy[N_REG][PD_SIM];
extern double rh_n2o_emit_ngas[N_REG][PD_SIM], rh_n2o_emitagr_ngas[N_REG][PD_SIM];

/* CHASER 2001 monthly, by A.Ito (2010/05/21) */
extern double	ndepo_chaser_dnhx[ASTEP][64][128];		/* NHx, dry */
extern double	ndepo_chaser_dnoy[ASTEP][64][128];		/* NOy, dry */
extern double	ndepo_chaser_wnhx[ASTEP][64][128];		/* NHx, wet */
extern double	ndepo_chaser_wnoy[ASTEP][64][128];		/* NOy, wet */

/* CHASER4.0 monthly, by A.Ito (2014/11/19) */
extern double	ndepo_chaser4_nhx_h[ASTEP][64][128];		/* NHx */
extern double	ndepo_chaser4_noy_h[ASTEP][64][128];		/* NOy */
extern double	ndepo_chaser4_ont_h[ASTEP][64][128];		/* Org NOx */
extern double	ndepo_chaser4_nhx_p[ASTEP][64][128];		/* NHx */
extern double	ndepo_chaser4_noy_p[ASTEP][64][128];		/* NOy */
extern double	ndepo_chaser4_ont_p[ASTEP][64][128];		/* Org NOx */

/* CLEARANCE *****************************************************/
void f_clear(struct Grid *grid, struct Loct *loct, struct Echar *echar, 
		   struct Mass *mass, struct Flux *flux);
void plant_flux_zero(long month, struct Pflx *flux);
void vanish(struct Mass *mass, struct Flux *flux);
void vlzero(struct Grid *grid, struct Pmas *mass, struct Pflx *flux);

/* INITIALIZATION *********************************************/
void open_input(FILE *fp_s[IFILEN], FILE *fp_c[4]);
void f_output_file_open(short vtype, short zone, char s_date[32], char s_case[32], 
	char filename[128], FILE *fp[OFILEN]);
void f_init_sim(struct Grid *grid);
void f_init_grid(FILE *fp_r[IFILEN], struct Grid *grid); 
void f_init_clim(struct Grid *grid);
void f_init_loct(struct Grid *grid, struct Loct *loct, struct Mass *mass, 
	struct Flux *flux, struct Echar *echar);
void f_dyn_loct(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Echar *echar);

void f_co2_trend(struct Grid *grid);
void read_gcm_clim(FILE *fp_c2[4], struct Grid *grid);
void read_ncep_clim(struct Grid *grid);
void set_gcm_clim(struct Grid *grid);
void set_hist_clim(struct Grid *grid);
void read_cru_clim(FILE *fp_c[4], FILE *fp_c2[4], struct Grid *grid);
void f_cult_luc(struct Grid *grid);
long f_basin_id_trip(long original);
long region_giorgi(double lat, double lon);
void f_read_chaser_ndepo(FILE *fp_s[IFILEN],struct Grid  *grid);

void f_parameter_perturbation(long iseed,struct Grid *grid,struct Loct *loct,struct Echar *echar,
	double f_prtrb[20]);

/* MASS & PARAMETERS INITIALIZATION *****************************/
void initVS(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux, 
	struct Echar *echar);
void parameterC3(struct Grid *grid, struct Pchar *C3);
void parameterC4(struct Grid *grid, struct Pchar *C4);
void parameterSoil(struct Grid *grid, struct Schar *Soil);
void parameterC3_sage(struct Grid *grid, struct Pchar *C3);
void parameterC4_sage(struct Grid *grid, struct Pchar *C4);
void parameterSoil_sage(struct Grid *grid, struct Schar *Soil);
void parameterCrop(struct Grid *grid, struct Pchar *C3);
void parameterSoil_crop(struct Grid *grid, struct Schar *Soil);

/* EXPERIMENTAL STEPS *********************************************/
void cal_spinup(struct Grid *grid, struct Loct *loct, 
	struct Echar *echar, struct Mass *mass, struct Flux *flux, FILE *fp[OFILEN]);
void cal_historical(struct Grid *grid, struct Loct *loct, 
	struct Echar *echar, struct Mass *mass, struct Flux *flux, FILE *fp[OFILEN]);
void cal_projection(struct Grid *grid, struct Loct *loct, 
	struct Echar *echar, struct Mass *mass, struct Flux *flux, FILE *fp[OFILEN]);

/* RADIATION *****************************************/
double f_solar_decl(struct Grid *grid);
double f_solar_hgt(struct Grid *grid);
double f_day_length(struct Grid *grid);
double f_top_rad(struct Grid *grid, short ha);
double f_gl_rad(struct Grid *grid);
double f_par(struct Grid *grid);
void f_par_h(struct Grid *grid, struct Loct *loct);
void f_net_rad(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Echar *echar);
double albedo_soil(struct Loct *loct, struct Schar *schar);

/* PHYTOGEOGRAPHY ****************************************/
void c34composition(short v_type, struct Grid *grid, struct Loct *loct);
void alternate_c34(struct Grid *grid, struct Loct *loct, struct Echar *echar, 
	struct Mass *mass, struct Flux *flux);
void phenol_tran_c34(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
	struct Pmas *mass, struct Pflx *flux);

/* HYDROLOGY *********************************************/
void f_waterbudget(struct Grid *grid, struct Loct *loct, struct Echar *echar);
double air_density(struct Grid *grid, struct Loct *loct);
double vap_pre_sat(struct Grid *grid);
double slope_vps(struct Grid *grid);
double r_aero(struct Grid *grid);
double pm_evaporation(struct Grid *grid, struct Loct *loct);
double pm_transpiration(struct Grid *grid, struct Loct *loct);
double pm_interception(struct Grid *grid, struct Loct *loct);

/* BIOLOGICAL SCHEMES *******************************************/
void f_biome_processes(struct Grid *grid, struct Loct *loct, 
		struct Echar *echar, struct Mass *mass, struct Flux *flux);
void forest_process(struct Grid *grid, struct Loct* loct, struct Pflx *pflx,
		struct Pchar *pchar, struct Pmas *mass);
void grass_process(struct Grid *grid, struct Loct *loct, struct Pflx *pflx, 
		struct Pchar *c34char, struct Pmas *mass);
void agri_process(struct Grid *grid, struct Loct *loct, struct Pflx *pflx, 
		struct Pchar *c3char, struct Pmas *mass);
void soil_processes(struct Grid *grid, struct Loct *loct, struct Schar *schar, 
		struct Smas *mass, struct Sflx *flux);

/* SUB-SCHEMES *************************/
void set_rowcol_gcm(void);
void set_gcm_index(char gcmindex[]);
void f_before_deal(struct Grid *grid, struct Pflx *flux);
void f_after_deal(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *pmas, struct Pflx *flux);
void f_plant_stand_budget(struct Grid *grid,struct Loct *loct,struct Mass *mass, struct Flux *flux);
double grid_area(double lat1, double lat2, double lon1, double lon2);

/* ECOPHYSIOLOGY *****************************/
void f_ecophysiology(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mass);
void incel_cdc(struct Grid *grid, struct Loct *loct, struct Pchar *plant);
void quantum_yield(struct Grid *grid, struct Pchar *plant);
void opt_lai(struct Grid *grid,struct Loct *loct, struct Pchar *plant);
double lai_mass(struct Grid *grid, struct Pmas *mass, struct Pchar *pchar);
double irr_attn(struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void f_qten_ar(struct Grid *grid, struct Pchar *plant);
void spcfc_res_mass(struct Pchar *plant, struct Pmas *mass);
void mortality(struct Grid *grid, struct Pchar *plant);
void stom_cond(struct Grid *grid,struct Loct *loct, struct Pchar *plant);
double canopy_cond(struct Grid *grid,struct Loct *loct, struct Pchar *plant, struct Pmas *mass);
void pc_sat(struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void f_leaf_age(short update, struct Pchar *pchar, struct Pmas *mass, double addshed);

/* PHENOLOGICAL CYCLES *****************************************/
void growthperiod(struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void phenology(struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void phenology_evergreen(struct Grid *grid, struct Pchar *pchar);
void phenology_colddeciduous(struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void phenology_ariddeciduous(struct Grid *grid, struct Pchar *pchar);
void phenology_bareland(struct Grid *grid, struct Pchar *pchar);
void phenology_grass(struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void phenology_agriculture(struct Grid *grid, struct Loct *loct, struct Pchar *pchar);

/* C FLUXES *********************************************************************/
/* photosynthesis */
double fgpp(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mas);
double f_df97_gpp(short mode, struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
		struct Pmas *pmas);
void leaf_ansolv(double vpd, double aco2, double lai, double ft_rd, double vcmax, 
	double k_effc, double comp_co2, double j_et, double *assim, double *rdark,
	double *stocon, double *inco2, double *debug1, double *debug2, double *debug3);
/* respiration */
double frfm(struct Grid *grid, struct Pchar *pchar, struct Pmas *pmas);
double frfg(struct Grid *grid, struct Pchar *pchar, struct Pflx *flux);
double frcm(struct Grid *grid, struct Pchar *pchar, struct Pmas *pmas);
double frcg(struct Grid *grid, struct Pchar *pchar, struct Pflx *flux);
double frrm(struct Grid *grid, struct Pchar *pchar, struct Pmas *pmas);
double frrg(struct Grid *grid, struct Pchar *pchar, struct Pflx *flux);
/* litterfall */
double flf(struct Grid *grid, struct Pchar *pchar, struct Pmas *pmas);
double flc(struct Grid *grid, struct Pchar *pchar, struct Pmas *pmas);
double flr(struct Grid *grid, struct Pchar *pchar, struct Pmas *pmas);
/* allocation */
void allocation(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mass, struct Pflx *flux);
void recluit(struct Grid *grid, struct Pchar *pchar, struct Pmas *mass);
void reallocation_survival(struct Grid *grid, struct Pchar *pchar, struct Pmas *mass);
/* soil */
double frl(struct Grid *grid, struct Loct *loct,struct Schar *schar, struct Smas *smas);
double frh(struct Grid *grid, struct Loct *loct,struct Schar *schar, struct Smas *smas);
double fsf(struct Grid *grid, struct Schar *schar, struct Sflx *sflx);

void npp_empirical(struct Grid *grid, struct Loct *loct, struct Flux *flux);

/* CROP CYCLES *****************************************************/
void planting(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
	struct Pmas *mass, struct Pflx *flux);
void harvesting(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
	struct Pmas *mass, struct Pflx *flux);
void interval(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
	struct Pmas *mass, struct Pflx *flux);

/* CARBON DYNAMICS *****************************************************/
void greenperiod(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
	struct Pmas *mass, struct Pflx *flux);					
void leaffall(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
	struct Pmas *mass, struct Pflx *flux);
void noleafperiod(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
	struct Pmas *mass, struct Pflx *flux);
void leafemergence(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
	struct Pmas *mass, struct Pflx *flux);

/* STABLE CARBON ISOTOPE *************************************************/
void f_init_c_isotpes(struct Grid *grid, struct Flux *flux, struct Echar *echar, struct Mass *mass);
void photo_13c_frac(struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
double deltaTratio(double delta);
double ratioTdelta(double ratio);
double d13c_addition(double d13c_a, double mass_a, double d13c_b, double mass_b);
double d13c_addition3(double d13c_a, double mass_a, double d13c_b, double mass_b, double d13c_c, double mass_c);
void co2_in_canopy(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);	
void f_cisotope_efflux(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
double f_decay_14c(double init_val);

/* MINOR PROCESSES ****************************************************/
void f_erosion(struct Grid *grid, struct Loct *loct, struct Echar *echar, struct Mass *mass, struct Flux *flux);
void f_biomassburning(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
void f_luc_emit(struct Grid *grid, struct Mass *mass, struct Flux *flux);
void f_voc_emit_guenther97(struct Grid *grid, struct Loct *loct, struct Echar *echar, 
	struct Mass *mass, struct Flux *flux);
void f_doc_boyer(struct Grid *grid, struct Loct *loct, struct Smas *mass, struct Sflx *sflx);

/* GHG PROCESSES ****************************************************/
void f_casa_moisture(struct Grid *grid, struct Loct *loct);
void f_ch4oxy_ridgewell(struct Grid *grid, struct Loct *loct, struct Flux *flux);
void f_ch4oxy_casa(struct Grid *grid, struct Loct *loct, struct Flux *flux);
void f_ch4oxy_delgrosso(struct Grid *grid, struct Loct *loct, struct Flux *flux);
void f_ch4oxy_curry(struct Grid *grid, struct Loct *loct, struct Flux *flux);
void f_ch4_emit_cao(struct Grid *grid, struct Loct *loct, struct Flux *flux);
void f_ch4_emit_walter(short mode, struct Grid *grid, struct Loct *loct, struct Flux *flux);
void f_n2o_emit_ngas(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
void f_n2o_emit_casa(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
void f_ch4_emit_veg(struct Grid *grid, struct Loct *loct, struct Echar *echar, struct Mass *mass, struct Flux *flux);

void ghg_flux_zero(long month, struct Flux *flux);
void n_budget(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
void f_n_deposit(struct Grid *grid, struct Loct *loct);
void n_fertilizer_in(struct Grid *grid, struct Loct *loct);
void f_biolfix(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
void f_n_uptake(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
void f_n_abandon_salvage(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
		struct Pmas *mass, struct Pflx *flux);
void f_n_alloc(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mass, struct Pflx *flux);
void f_n_realloc(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mass, struct Pflx *flux);
void f_nh3_volatilization(struct Grid *grid, struct Loct *loct, struct Schar *schar, struct Smas *mass, struct Sflx *sflx);
void f_n_leaching(struct Grid *grid, struct Loct *loct, struct Smas *mass, struct Sflx *sflx);
void f_n_mineralz(struct Grid *grid, struct Loct *loct, struct Smas *mass, struct Sflx *sflx);
void f_n_leaf_conc(struct Grid *grid, struct Pchar *pchar, struct Pmas *mass);
void f_n_mcrb_abdn(struct Grid *grid, struct Loct *loct, struct Schar *schar, struct Smas *smas, struct Sflx *flux);
void f_n_immoblz(struct Grid *grid, struct Loct *loct, struct Schar *schar, struct Smas *smas, struct Sflx *flux);

/* OUTPUT *************************************/
void f_set_history_data(long year, struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
void f_output_result(long year, struct Grid *grid, struct Loct *loct, struct Echar *echar, struct Mass *mass, 
	struct Flux *flux, FILE *fp_o[OFILEN]);
void screenshow(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux, struct Echar *echar);
void publish_cbud(struct Grid*grid, struct Loct *loct, struct Echar *echar, 
	struct Mass *mass, struct Flux *flux, FILE *result);
void f_glosum_output(char sdate[32], char scase[32]);
void f_grid_av(struct Grid *grid, struct Loct *loct, struct Echar *echar, struct Mass *mass, struct Flux *flux);
