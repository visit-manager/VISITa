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

#define IFILEN 52
#define OFILES 8

extern short DF97;
extern double MDN[ASTEP];
extern long GCM, CO2S, GCM_R, GCM_C;
extern double aco2_a1[553], aco2_a2[553], aco2_b1[553], aco2_b2[553];
extern double ach4_a1[553], ach4_a2[553], ach4_b1[553], ach4_b2[553];
extern double an2o_a1[553], an2o_a2[553], an2o_b1[553], an2o_b2[553];

extern double glandarea;
extern double h_tmp[HIST], h_pre[HIST], h_dswr[HIST], h_aet[HIST], h_rof[HIST];
extern double h_gpp[HIST], h_npp[HIST], h_nep[HIST], h_plant[HIST], h_soil[HIST];
extern double h_sr[HIST], h_ersn_c[HIST], h_agrersn_c[HIST], h_doc[HIST];
extern double h_agrarea[HIST], h_luc[HIST];
extern double h_luc_1[HIST], h_luc_2[HIST], h_luc_3[HIST];
extern double h_gpp_df97[HIST], h_gpp_c4[HIST]; /* added by A.Ito (2009/08/31) */
extern double h_pot_prmfrst[HIST];

extern double h_burnt_area[HIST];
extern double h_bioburn_co2[HIST], h_bioburn_ch4[HIST], h_bioburn_co[HIST];
extern double h_bioburn_nmhc[HIST], h_bioburn_oc[HIST], h_bioburn_bc[HIST];
extern double h_bioburn_nox[HIST], h_bioburn_so2[HIST], h_bioburn_pm25[HIST];
extern double h_bioburn_tpm[HIST], h_bioburn_tec[HIST];

extern double h_ch4ox1[HIST], h_ch4ox2[HIST], h_ch4ox3[HIST], h_ch4ox4[HIST];
extern double h_ch4emit_cao_paddy[HIST], h_ch4emit_cao_wetland[HIST];
extern double h_ch4emit_paddy_wh_diff[HIST],  h_ch4emit_paddy_wh_plant[HIST];
extern double h_ch4emit_paddy_wh_ebbl[HIST],  h_ch4emit_paddy_wh_release[HIST];
extern double h_ch4emit_wetland_wh_diff[HIST],  h_ch4emit_wetland_wh_plant[HIST];
extern double h_ch4emit_wetland_wh_ebbl[HIST],  h_ch4emit_wetland_wh_release[HIST];

extern double h_n2o_emit_ngas[HIST], h_n2_emit_ngas[HIST];
extern double h_n2o_emit_casa[HIST], h_no_emit_casa[HIST], h_n2_emit_casa[HIST];
extern double h_nh3_emit[HIST], h_n2_biofix[HIST];
extern double h_ch4_emit_mass[HIST], h_ch4_emit_photo[HIST];
extern double h_n2o_d_emit_ngas[HIST], h_n2o_n_emit_ngas[HIST];
extern double h_n2o_emit_ngas_agr[HIST], h_n2o_emit_casa_agr[HIST];
extern double h_nh3_emit_agr[HIST];
extern double h_no3_leach[HIST];
extern double h_n_fertin[HIST], h_n_depoin[HIST]; /* added by A.Ito (2010/05/02) */

extern double h_voc_isopr_g97[HIST], h_voc_monotrp_g97[HIST], h_voc_methanl_g97[HIST];
extern double h_voc_acetone_g97[HIST], h_voc_actaldhd_g97[HIST], h_voc_frmardhd_g97[HIST];
extern double h_voc_formacd_g97[HIST], h_voc_acetacd_g97[HIST], h_voc_co_g97[HIST];

extern double ci_aco2[HIST], ci_aco2_d13c[HIST], ci_aco2_d14c[HIST];
extern double ci_gpp[HIST], ci_gpp_d13c[HIST], ci_gpp_d14c[HIST];
extern double ci_er[HIST], ci_er_d13c[HIST], ci_er_d14c[HIST];
extern double ci_f[HIST], ci_f_d13c[HIST], ci_f_d14c[HIST];
extern double ci_c[HIST], ci_c_d13c[HIST], ci_c_d14c[HIST];
extern double ci_r[HIST], ci_r_d13c[HIST], ci_r_d14c[HIST];
extern double ci_l[HIST], ci_l_d13c[HIST], ci_l_d14c[HIST];
extern double ci_h[HIST], ci_h_d13c[HIST], ci_h_d14c[HIST];

/* monthly results **********/
extern double m_ch4ox1[12], m_ch4ox2[12], m_ch4ox3[12];
extern double m_bioburn_co2[12], m_bioburn_ch4[12], m_bioburn_co[12];
extern double m_bioburn_nmhc[12], m_bioburn_oc[12], m_bioburn_bc[12];
extern double m_gpp[ASTEP], m_npp[ASTEP], m_nep[ASTEP];
extern double m_ch4p_cao[ASTEP], m_ch4p_wh[ASTEP];

/* vegetation (olson) results */
extern double go_landarea, gs_landarea;
extern double vo_area[34];
extern double vo_gpp[34], vo_npp[34], vo_nep[34];
extern double vo_lai[34], vo_fol[34], vo_stm[34], vo_rot[34], vo_ltr[34], vo_msl[34];
extern double vs_area[16];
extern double vs_gpp[16], vs_npp[16], vs_nep[16];
extern double vs_lai[16], vs_fol[16], vs_stm[16], vs_rot[16], vs_ltr[16], vs_msl[16];

extern float g_tmp[5][360][720];
extern float g_prc[5][360][720];
extern float g_swr[5][360][720];
extern float g_gpp[5][360][720];
extern float g_npp[5][360][720];
extern float g_nep[5][360][720];
extern float g_pmas[5][360][720];
extern float g_smas[5][360][720];
extern float g_ch4e_cao[5][360][720];
extern float g_ch4o_curry[5][360][720];
extern float g_n2oe[5][360][720];
extern float g_bbco2[5][360][720];
extern float g_ersn[5][360][720];
extern float g_isopr[5][360][720];
extern float g_sr[5][360][720];
extern float g_luc[5][360][720];

#if C13_GOUT==1
extern float g_f13[5][360][720]; 
extern float g_c13[5][360][720]; 
extern float g_r13[5][360][720]; 
extern float g_l13[5][360][720]; 
extern float g_h13[5][360][720]; 
extern float g_gpp13[5][360][720]; 
extern float g_er13[5][360][720]; 
#endif

#if C14_GOUT==1
extern float g_f14[5][360][720]; 
extern float g_c14[5][360][720]; 
extern float g_r14[5][360][720]; 
extern float g_l14[5][360][720]; 
extern float g_h14[5][360][720]; 
extern float g_gpp14[5][360][720]; 
extern float g_er14[5][360][720]; 
#endif

extern float g_er[5][360][720]; 
extern float g_snh4[5][360][720]; 
extern float g_sno3[5][360][720]; 

#if PHYS_GOUT==1
extern float g_lai[5][360][720]; 
extern float g_parb[5][360][720]; 
extern float g_pard[5][360][720]; 
extern float g_apar[5][360][720]; 
extern float g_apar2[5][360][720];
extern float g_aet[5][360][720]; 
extern float g_rof[5][360][720]; 
extern float g_rns[5][360][720]; 
extern float g_rnl[5][360][720]; 
extern float g_sw1[5][360][720]; 
extern float g_sw2[5][360][720]; 
#endif

#if CH4_WH==1
extern float g_ch4ep_wh[5][360][720]; 
extern float g_ch4ew_wh[5][360][720]; 
extern float gm_ch4ep_wh[12][360][720];
#endif
extern float g_ch4ep_cao[5][360][720]; 

/* regional historical */
extern double rh_area[NREG];
extern double rh_temp[NREG][HIST], rh_prec[NREG][HIST], rh_dswrf[NREG][HIST];
extern double rh_rns[NREG][HIST], rh_rnl[NREG][HIST];
extern double rh_ipar[NREG][HIST], rh_apar[NREG][HIST];
extern double rh_gpp[NREG][HIST], rh_npp[NREG][HIST], rh_nep[NREG][HIST];
extern double rh_evpr[NREG][HIST], rh_trsp[NREG][HIST], rh_incp[NREG][HIST], rh_rnof[NREG][HIST];
extern double rh_ci_gpp[NREG][HIST], rh_ci_gpp_d13c[NREG][HIST], rh_ci_gpp_d14c[NREG][HIST];
extern double rh_ci_er[NREG][HIST], rh_ci_er_d13c[NREG][HIST], rh_ci_er_d14c[NREG][HIST];
extern double rh_ci_f[NREG][HIST], rh_ci_f_d13c[NREG][HIST], rh_ci_f_d14c[NREG][HIST];
extern double rh_ci_c[NREG][HIST], rh_ci_c_d13c[NREG][HIST], rh_ci_c_d14c[NREG][HIST];
extern double rh_ci_r[NREG][HIST], rh_ci_r_d13c[NREG][HIST], rh_ci_r_d14c[NREG][HIST];
extern double rh_ci_l[NREG][HIST], rh_ci_l_d13c[NREG][HIST], rh_ci_l_d14c[NREG][HIST];
extern double rh_ci_h[NREG][HIST], rh_ci_h_d13c[NREG][HIST], rh_ci_h_d14c[NREG][HIST];

extern double rh_hvst[NREG][HIST], rh_luc[NREG][HIST];
extern double rh_ch4ox_curry[NREG][HIST], rh_ch4emit_wh_wet[NREG][HIST], rh_ch4emit_wh_paddy[NREG][HIST];
extern double rh_n2o_emit_ngas[NREG][HIST], rh_n2o_emitagr_ngas[NREG][HIST];

/* CLEARANCE *****************************************************/
void clear(struct Grid *grid, struct Loct *loct, struct Echar *echar, 
		   struct Mass *mass, struct Flux *flux);
void plant_flux_zero(long month, struct Pflx *flux);
void vanish(struct Mass *mass, struct Flux *flux);
void vlzero(struct Grid *grid, struct Pmas *mass, struct Pflx *flux);

/* INITIALIZATION *********************************************/
void open_input(FILE *fp_s[IFILEN], FILE *fp_c[4]);
void f_output_file_open(short vtype, short zone, char s_date[25], char s_case[25], 
	char filename[100], FILE *fp[OFILES]);
void f_init_sim(struct Grid *grid);
void f_init_grid(FILE *fp_r[IFILEN], struct Grid *grid); 
void initC(struct Grid *grid);
void initL(struct Grid *grid, struct Loct *loct, struct Mass *mass, 
	struct Flux *flux, struct Echar *echar);
void dynmcL(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Echar *echar);

void co2_trend(struct Grid *grid);
void read_gcm_clim(struct Grid *grid);
void read_ncep_clim(struct Grid *grid);
void set_gcm_clim(struct Grid *grid);
void set_cru_clim(struct Grid *grid);
void read_cru_clim(FILE *fp_c[4], struct Grid *grid);
void f_cult_luc(struct Grid *grid);
long basin_id_trip(long original);
long region_giorgi(double lat, double lon);

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
void cal_stable(struct Grid *grid, struct Loct *loct, 
	struct Echar *echar, struct Mass *mass, struct Flux *flux, FILE *fp[OFILES]);
void cal_cruclim(struct Grid *grid, struct Loct *loct, 
	struct Echar *echar, struct Mass *mass, struct Flux *flux, FILE *fp[OFILES]);
void cal_gcmclim2(struct Grid *grid, struct Loct *loct, 
	struct Echar *echar, struct Mass *mass, struct Flux *flux, FILE *fp[OFILES]);

/* RADIATION *****************************************/
double sl_dec(struct Grid *grid);
double sl_hgt(struct Grid *grid);
double dlen(struct Grid *grid);
double top_rad(struct Grid *grid);
double gl_rad(struct Grid *grid);
double par(struct Grid *grid);
void net_rad(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Echar *echar);
double albedo_soil(struct Loct *loct, struct Schar *schar);

/* PHYTOGEOGRAPHY ****************************************/
void c34composition(short v_type, struct Grid *grid, struct Loct *loct);
void alternate_c34(struct Grid *grid, struct Loct *loct, struct Echar *echar, 
	struct Mass *mass, struct Flux *flux);
void phenol_tran_c34(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
	struct Pmas *mass, struct Pflx *flux);

/* HYDROLOGY *********************************************/
void waterbudget(struct Grid *grid, struct Loct *loct, struct Echar *echar);
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
void f_after_deal(struct Grid *grid, struct Pchar *pchar, struct Pmas *pmas, struct Pflx *flux);
void f_plant_stand_budget(struct Grid *grid,struct Loct *loct,struct Mass *mass, struct Flux *flux);
double grid_area(double lat1, double lat2, double lon1, double lon2);

/* ECOPHYSIOLOGY *****************************/
void ecophysiology(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mass);
void incel_cdc(struct Grid *grid, struct Loct *loct, struct Pchar *plant);
void quantum_yield(struct Grid *grid, struct Pchar *plant);
void opt_lai(struct Grid *grid,struct Loct *loct, struct Pchar *plant);
double lai_mass(struct Grid *grid, struct Pmas *mass, struct Pchar *pchar);
double irr_attn(struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void qten_ar(struct Grid *grid, struct Pchar *plant);
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
void allocation(struct Grid *grid, struct Pchar *pchar, struct Pmas *mass, struct Pflx *flux);
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
void f_casa_mositure(struct Grid *grid, struct Loct *loct);
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
void f_biolfix(struct Grid *grid, struct Loct *loct, struct Flux *flux);
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
	struct Flux *flux, FILE *fp_o[OFILES]);
void screenshow(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux, struct Echar *echar);
void publish_cbud(struct Grid*grid, struct Loct *loct, struct Echar *echar, 
	struct Mass *mass, struct Flux *flux, FILE *result);
void f_glosum_output(char sdate[25], char scase[25]);
void f_grid_av(struct Grid *grid, struct Loct *loct, struct Echar *echar, struct Mass *mass, struct Flux *flux);
