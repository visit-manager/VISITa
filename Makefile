# Makefile for VISIT-a by Akihiko Ito

CC = cc
LIB = -lm
OBJS = allocation.o assimilation.o atm_co2.o bioburn.o biome.o cal_cruclim.o cal_gcmclim2.o cal_stable.o canopy.o ch4emit.o ch4ox.o clear.o clim_cru.o clim_gcm.o clim_series.o cultpro.o doc.o ecophysiology.o erosion.o forestpro.o grasspro.o hydrology.o init_grid.o init_sim.o init_vs.o isotope.o landuse.o litterfall.o location.o main.o n2oemit.o n_fertilizer.o n_flux.o open_input.o output_file.o output_screen.o parameter_olson.o parameter_sage.o phenology.o phytgeo_c34.o radiation.o respiration.o soilpro.o vegetdeal.o vocemit.o npp_empirical.o parameter_crop.o

visita: $(OBJS)
	$(CC) -o $@ $(OBJS) $(LIB)

.c.o:
	$(CC) -c $<

allocation.o: structure.h prototype.h setting.h
assimilation.o: structure.h prototype.h setting.h
atm_co2.o: structure.h prototype.h setting.h
bioburn.o: structure.h prototype.h setting.h
biome.o: structure.h prototype.h setting.h
cal_cruclim.o: structure.h prototype.h setting.h
cal_gcmclim2.o: structure.h prototype.h setting.h
cal_stable.o: structure.h prototype.h setting.h
canopy.o: structure.h prototype.h setting.h
ch4emit.o: structure.h prototype.h setting.h
ch4ox.o: structure.h prototype.h setting.h
clear.o: structure.h prototype.h setting.h
clim_cru.o: structure.h prototype.h setting.h
clim_gcm.o: structure.h prototype.h setting.h
clim_series.o: structure.h prototype.h setting.h
cultpro.o: structure.h prototype.h setting.h
doc.o: structure.h prototype.h setting.h
ecophysiology.o: structure.h prototype.h setting.h
erosion.o: structure.h prototype.h setting.h
forestpro.o: structure.h prototype.h setting.h
grasspro.o: structure.h prototype.h setting.h
hydrology.o: structure.h prototype.h setting.h
init_grid.o: structure.h prototype.h setting.h
init_sim.o: structure.h prototype.h setting.h
init_vs.o: structure.h prototype.h setting.h
isotope.o: structure.h prototype.h setting.h
landuse.o: structure.h prototype.h setting.h
litterfall.o: structure.h prototype.h setting.h
location.o: structure.h prototype.h setting.h
main.o: structure.h prototype.h setting.h global_var.h
n2oemit.o: structure.h prototype.h setting.h
n_fertilizer.o: structure.h prototype.h setting.h
n_flux.o: structure.h prototype.h setting.h
npp_empirical.o: structure.h prototype.h setting.h
open_input.o: structure.h prototype.h setting.h
output_file.o: structure.h prototype.h setting.h
output_screen.o: structure.h prototype.h setting.h
parameter_olson.o: structure.h prototype.h setting.h
parameter_sage.o: structure.h prototype.h setting.h
parameter_crop.o: structure.h prototype.h setting.h
phenology.o: structure.h prototype.h setting.h
phytgeo_c34.o: structure.h prototype.h setting.h
radiation.o: structure.h prototype.h setting.h
respiration.o: structure.h prototype.h setting.h
soilpro.o: structure.h prototype.h setting.h
vegetdeal.o: structure.h prototype.h setting.h
vocemit.o: structure.h prototype.h setting.h

clean: 
	rm -f $(OBJS)
