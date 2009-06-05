/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/*  Updated August 7, 2007								*/
/*  Updated November 29, 2007							*/

#define IFILEN 48
#define OFILES 7

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
void initSim(struct Grid *grid);
void init_grid(FILE *fp_r[IFILEN], struct Grid *grid); 
void initC(struct Grid *grid);
void initL(struct Grid *grid, struct Loct *loct, struct Mass *mass, 
	struct Flux *flux, struct Echar *echar);
void dynmcL(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Echar *echar);

void cd_trend(struct Grid *grid);
void read_gcm_clim(struct Grid *grid);
void read_ncep_clim(struct Grid *grid);
void set_gcm_clim(struct Grid *grid);
void set_cru_clim(struct Grid *grid);
void read_cru_clim(FILE *fp_c[4], struct Grid *grid);
void f_cult_luc(struct Grid *grid);
long basin_id_trip(long original);

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
void biome_processes(struct Grid *grid, struct Loct *loct, 
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
void beforedeal(struct Grid *grid, struct Pflx *flux);
void afterdeal(struct Grid *grid, struct Pchar *pchar, struct Pmas *pmas, struct Pflx *flux);
void plant_stand(struct Grid *grid,struct Loct *loct,struct Mass *mass, struct Flux *flux);
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
void init_d13c(struct Grid *grid, struct Flux *flux, struct Echar *echar, struct Mass *mass);
void photo_13c_frac(struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
double deltaTratio(double delta);
double ratioTdelta(double ratio);
double d13c_addition(double d13c_a, double mass_a, double d13c_b, double mass_b);
double d13c_addition3(double d13c_a, double mass_a, double d13c_b, double mass_b, double d13c_c, double mass_c);
void co2_in_canopy(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);	
void d13c_efflux(struct Grid *grid, struct Loct *loct, struct Flux *flux);

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
void f_ch4_emit_walter(struct Grid *grid, struct Loct *loct, struct Flux *flux);
void f_n2o_emit_ngas(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
void f_n2o_emit_casa(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
void f_ch4_emit_veg(struct Grid *grid, struct Loct *loct, struct Echar *echar, struct Mass *mass, struct Flux *flux);

void n_flux_zero(long month, struct Flux *flux);
void n_budget(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
void f_n_deposit(struct Grid *grid, struct Loct *loct);
void n_fertilizer_in(struct Grid *grid, struct Loct *loct);
void f_biolfix(struct Grid *grid, struct Loct *loct, struct Flux *flux);
void f_n_uptake(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
void f_n_abandon_salvage(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
		struct Pmas *mass, struct Pflx *flux);
void f_n_alloc(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mass, struct Pflx *flux);
void f_n_realloc(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mass, struct Pflx *flux);
void f_nh3_volatilization(struct Grid *grid, struct Loct *loct, struct Smas *mass, struct Sflx *sflx);
void f_n_leaching(struct Grid *grid, struct Loct *loct, struct Smas *mass, struct Sflx *sflx);
void f_n_mineralz(struct Grid *grid, struct Loct *loct, struct Smas *mass, struct Sflx *sflx);
void f_n_leaf_conc(struct Grid *grid, struct Pchar *pchar, struct Pmas *mass);
void f_n_mcrb_abdn(struct Grid *grid, struct Loct *loct, struct Schar *schar, struct Smas *smas, struct Sflx *flux);
void f_n_immoblz(struct Grid *grid, struct Loct *loct, struct Schar *schar, struct Smas *smas, struct Sflx *flux);

/* OUTPUT *************************************/
void f_set_history_data(long year, struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
void f_output_result(long year, struct Grid *grid, struct Loct *loct, struct Echar *echar, struct Mass *mass, 
	struct Flux *flux, FILE *fp_o[6]);
void screenshow(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux, struct Echar *echar);
void publish_cbud(struct Grid*grid, struct Loct *loct, struct Echar *echar, 
	struct Mass *mass, struct Flux *flux, FILE *result);
void f_glosum_output(char sdate[25], char scase[25]);
