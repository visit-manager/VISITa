/*	VISIT: Vegetation Integrative SImulation Tool						*/
/* Old name: Simulation model of Carbon cYCle in Land Ecosystems		*/
/* Developed by A.Ito in CGER/NIES & EAIMG/ECRP/FRSGC					*/
/* Carbon cycle, erosion, biomass burning, land-use change,				*/
/* CH4 emission and oxidation, N2O emission,,,,,						*/
/*	version 1.0.0	cerated in August 14, 2007							*/

/*   Last updated August 7, 2007								*/

/* prototype declaration */
#define dTr 0.0174533 /* angle conversion, from degree to radian */
#define rTd 57.29577951 /* angle conversion, from radian to degree */
#define PI 3.141592653 /** pai **/
#define cdTc 0.272727 /* from CO2-base to Crabon-base */
#define dmTc 2.2 /* from dry-matter-base to Carbon base */
#define cTdm 0.4545 /* from dry-matter-base to Carbon base */
#define lTs (3600.0*12.0/100000000.0) /* from micro-mol m-2 s-1 to Mg C ha-1 day-1*/
#define ZAT 273.15 /* zero degree centigrade in absolute temperature */

#define STCIR (0.0111/0.9889) /* standard stable carbon isotope ratio */

#define YES 1
#define NO 0

/***** INITIALIZATION *****/
void initSim(struct Grid *grid);
void initG(FILE *fp_r[26], struct Grid *grid); 
void initC(struct Grid *grid, struct Loct *loct);
void initL(struct Grid *grid, struct Loct *loct, struct Mass *mass, 
	struct Flux *flux, struct Echar *echar);
void dynmcL(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Echar *echar);

void cd_trend(struct Grid *grid);
void read_gcm_clim(struct Grid *grid);
void set_gcm_clim(struct Grid *grid);
void set_cru_clim(struct Grid *grid);
void read_cru_clim(FILE *fp_c[40], struct Grid *grid);
void f_cult_luc(struct Grid *grid);
long basin_id_trip(long original);

/***** CLEARANCE *****/
void clear(struct Grid *grid, struct Loct *loct, struct Echar *echar, 
		struct Mass *mass, struct Flux *flux);
void plant_flux_zero(long month, struct Pflx *flux);
void vanish(struct Mass *mass, struct Flux *flux);
void vlzero(struct Grid *grid, struct Pmas *mass, struct Pflx *flux);

/* EXPERIMENTAL STEPS ***********************/
void cal_stable(struct Grid *grid, struct Loct *loct, 
		struct Echar *echar, struct Mass *mass, struct Flux *flux, FILE *fp_carbon, 
		FILE *fp_nitrogen, FILE *fp_ersn, FILE *fp_ghg, FILE *fp_bioburn, FILE *fp_voc);
void cal_cruclim(struct Grid *grid, struct Loct *loct, 
		struct Echar *echar, struct Mass *mass, struct Flux *flux, FILE *fp_carbon, 
		FILE *fp_nitrogen, FILE *fp_ersn, FILE *fp_ghg, FILE *fp_bioburn, FILE *fp_voc);
void cal_gcmclim2(struct Grid *grid, struct Loct *loct, 
		struct Echar *echar, struct Mass *mass, struct Flux *flux, FILE *fp_carbon, 
		FILE *fp_nitrogen, FILE *fp_ersn, FILE *fp_ghg, FILE *fp_bioburn, FILE *fp_voc);

/***** RADIATION *****/
double sl_dec(struct Grid *grid);
double sl_hgt(struct Grid *grid);
double dlen(struct Grid *grid);
double top_rad(struct Grid *grid);
double gl_rad(struct Grid *grid);
double par(struct Grid *grid);
void net_rad(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Echar *echar);
double albedo_soil(struct Loct *loct, struct Schar *schar);

/***** MASS & PARAMETERS INITIALIZATION *****/
void initVS(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux, 
	struct Echar *echar);
void parameterC3(struct Grid *grid, struct Pchar *C3);
void parameterC4(struct Grid *grid, struct Pchar *C4);
void parameterSoil(struct Grid *grid, struct Schar *Soil);

/***** PHYTOGEOGRAPHY *****/
void c34composition(struct Grid *grid, struct Loct *loct);
void alternate_c34(struct Grid *grid, struct Loct *loct, struct Echar *echar, 
	struct Mass *mass, struct Flux *flux);
void phenol_tran_c34(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, 
	struct Pmas *mass, struct Pflx *flux);

/***** HYDROLOGY *****/
void waterbudget(struct Grid *grid, struct Loct *loct, struct Echar *echar);
double air_density(struct Grid *grid, struct Loct *loct);
double vap_pre_sat(struct Grid *grid);
double slope_vps(struct Grid *grid);
double r_aero(struct Grid *grid);
double pm_evaporation(struct Grid *grid, struct Loct *loct);
double pm_transpiration(struct Grid *grid, struct Loct *loct);
double pm_interception(struct Grid *grid, struct Loct *loct);

/***** SCHEMES ***************************************************************/
void biome_processes(struct Grid *grid, struct Loct *loct, 
		struct Echar *echar, struct Mass *mass, struct Flux *flux);
void forest_process(struct Grid *grid, struct Loct* loct, struct Pflx *flux,
		struct Pchar *pchar, struct Pmas *mass);
void grass_process(struct Grid *grid, struct Loct *loct, struct Pflx *flux, 
		struct Pchar *c34char, struct Pmas *mass);
void agri_process(struct Grid *grid, struct Loct *loct, struct Pflx *flux, 
		struct Pchar *c3char, struct Pmas *mass);
void soil_processes(struct Grid *grid, struct Loct *loct, struct Schar *schar, 
		struct Smas *mass, struct Sflx *flux);

/***** ECOPHYSIOLOGY *****/
void ecophysiology(struct Grid *grid, struct Loct *loct, struct Pchar *veg, struct Pmas *mass);
void incel_cdc(struct Grid *grid, struct Loct *loct, struct Pchar *plant);
void quantum_yield(struct Grid *grid, struct Pchar *plant);
void opt_lai(struct Grid *grid,struct Loct *loct, struct Pchar *plant);
double lai_mass(struct Pmas *mass, struct Pchar *veg);
double irr_attn(struct Grid *grid, struct Loct *loct, struct Pchar *veg);
void qten_ar(struct Grid *grid, struct Pchar *plant);
void spcfc_res_mass(struct Pchar *plant, struct Pmas *mass);
void mortality(struct Grid *grid, struct Pchar *plant);
void stom_cond(struct Grid *grid,struct Loct *loct, struct Pchar *plant);
double canopy_cond(struct Grid *grid,struct Loct *loct, struct Pchar *plant, struct Pmas *mass);
void pc_sat(struct Grid *grid, struct Loct *loct, struct Pchar *veg);

/***** PHENOLOGICAL CYCLES *****/
void growthperiod(struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void phenology(struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void phenology_evergreen(struct Grid *grid, struct Pchar *pchar);
void phenology_colddeciduous(struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void phenology_ariddeciduous(struct Grid *grid, struct Pchar *pchar);
void phenology_bareland(struct Grid *grid, struct Pchar *pchar);
void phenology_grass(struct Grid *grid, struct Loct *loct, struct Pchar *pchar);
void phenology_agriculture(struct Grid *grid, struct Loct *loct, struct Pchar *pchar);

/***** C FLUXES *********************************************************/
/* photosynthesis */
double fgpp(struct Grid *grid, struct Loct *loct, struct Pchar *veg, struct Pmas *mas);
/* respiration */
double frfm(struct Grid *grid, struct Pchar *veg, struct Pmas *mas);
double frfg(struct Grid *grid, struct Pchar *veg, struct Pflx *flux);
double frcm(struct Grid *grid, struct Pchar *veg, struct Pmas *mas);
double frcg(struct Grid *grid, struct Pchar *veg, struct Pflx *flux);
double frrm(struct Grid *grid, struct Pchar *veg, struct Pmas *mas);
double frrg(struct Grid *grid, struct Pchar *veg, struct Pflx *flux);
/* litterfall */
double flf(struct Grid *grid, struct Pchar *veg, struct Pmas *mas);
double flc(struct Grid *grid, struct Pchar *veg, struct Pmas *mas);
double flr(struct Grid *grid, struct Pchar *veg, struct Pmas *mas);
/* allocation */
void allocation(struct Grid *grid, struct Pchar *pchar, struct Pmas *mass, struct Pflx *flux);
void recluit(struct Grid *grid, struct Pchar *veg, struct Pmas *mass);
void reallocation_survival(struct Grid *grid, struct Pchar *veg, struct Pmas *mass);
/* soil */
double frl(struct Grid *grid, struct Loct *loct,struct Schar *soil, struct Smas *mas);
double frh(struct Grid *grid, struct Loct *loct,struct Schar *soil, struct Smas *mas);
double fsf(struct Grid *grid, struct Schar *soil, struct Sflx *flux);

void npp_empirical(struct Grid *grid, struct Loct *loct, struct Flux *flux);

/***** DEALINGS *****/
void beforedeal(struct Grid *grid, struct Pflx *flux);
void afterdeal(struct Grid *grid, struct Pchar *veg, struct Pmas *mass, struct Pflx *flux);
void plant_stand(struct Grid *grid,struct Loct *loct,struct Mass *mass, struct Flux *flux);
double grid_area(double lat1, double lat2, double lon1, double lon2);

/***** CROP CYCLES *****/
void planting(struct Grid *grid, struct Loct *loct, struct Pchar *veg, struct Pmas *mass, struct Pflx *flux);
void harvesting(struct Grid *grid, struct Loct *loct, struct Pchar *veg, struct Pmas *mass, struct Pflx *flux);
void interval(struct Grid *grid, struct Loct *loct, struct Pchar *veg, struct Pmas *mass, struct Pflx *flux);

/***** CARBON DYNAMICS *****/
void greenperiod(struct Grid *grid, struct Loct *loct, struct Pchar *veg, struct Pmas *mass, struct Pflx *flux);					
void leaffall(struct Grid *grid, struct Loct *loct, struct Pchar *veg, struct Pmas *mass, struct Pflx *flux);
void noleafperiod(struct Grid *grid, struct Loct *loct, struct Pchar *veg, struct Pmas *mass, struct Pflx *flux);
void leafemergence(struct Grid *grid, struct Loct *loct, struct Pchar *veg, struct Pmas *mass, struct Pflx *flux);

/***** STABLE CARBON ISOTOPE *****/
void init_d13c(struct Grid *grid, struct Flux *flux, struct Echar *echar, struct Mass *mass);
void photo_13c_frac(struct Grid *grid, struct Loct *loct, struct Pchar *plant);
double deltaTratio(double delta);
double ratioTdelta(double ratio);
double d13c_addition(double d13c_a, double mass_a, double d13c_b, double mass_b);
double d13c_addition3(double d13c_a, double mass_a, double d13c_b, double mass_b, double d13c_c, double mass_c);
void co2_in_canopy(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);	
void d13c_efflux(struct Grid *grid, struct Loct *loct, struct Flux *flux);

/* ADDITIONAL PROCESSES ************/
void f_erosion(struct Grid *grid, struct Loct *loct, struct Echar *echar, struct Mass *mass, struct Flux *flux);
void f_biomassburning(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
void f_luc_emit(struct Grid *grid, struct Mass *mass, struct Flux *flux);
void f_voc_emit_guenther97(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);

void f_casa_mositure(struct Grid *grid, struct Loct *loct);
void f_ch4oxy_ridgewell(struct Grid *grid, struct Loct *loct, struct Flux *flux);
void f_ch4oxy_casa(struct Grid *grid, struct Loct *loct, struct Flux *flux);
void f_ch4oxy_delgrosso(struct Grid *grid, struct Loct *loct, struct Flux *flux);
void f_ch4_emit_cao(struct Grid *grid, struct Loct *loct, struct Flux *flux);
void f_n2o_emit_ngas(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
void f_n2o_emit_casa(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
void f_ch4_emit_veg(struct Grid *grid, struct Loct *loct, struct Echar *echar, struct Mass *mass, struct Flux *flux);

void n_flux_zero(long month, struct Flux *flux);
void n_budget(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
void f_n_deposit(struct Grid *grid, struct Loct *loct);
void n_fertilizer_in(struct Grid *grid, struct Loct *loct);
void f_biolfix(struct Grid *grid, struct Loct *loct, struct Flux *flux);
void f_n_uptake(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux);
void f_n_abandon_salvage(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mass, struct Pflx *flux);
void f_n_alloc(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mass, struct Pflx *flux);
void f_n_realloc(struct Grid *grid, struct Loct *loct, struct Pchar *pchar, struct Pmas *mass, struct Pflx *flux);
void f_nh3_volatilization(struct Grid *grid, struct Loct *loct, struct Smas *mass, struct Sflx *flux);
void f_n_leaching(struct Grid *grid, struct Loct *loct, struct Smas *mass, struct Sflx *flux);
void f_n_mineralz(struct Grid *grid, struct Loct *loct, struct Smas *mass, struct Sflx *flux);
void f_n_leaf_conc(struct Grid *grid, struct Pchar *pchar, struct Pmas *mass);
void f_n_mcrb_abdn(struct Grid *grid, struct Loct *loct, struct Schar *schar, struct Smas *mass, struct Sflx *flux);
void f_n_immoblz(struct Grid *grid, struct Loct *loct, struct Schar *schar, struct Smas *mass, struct Sflx *flux);

/*** OUTPUT ***/
void screenshow(struct Grid *grid, struct Loct *loct, struct Mass *mass, struct Flux *flux, struct Echar *echar);
void publish_cbud(struct Grid*grid, struct Loct *loct, struct Echar *echar, 
		struct Mass *mass, struct Flux *flux, FILE *result);
