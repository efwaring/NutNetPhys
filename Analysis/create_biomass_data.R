# create site trait dataset for analyses

#########################################
# packages necessary
#########################################

library(lme4)
library(car)
library(emmeans)
library(dplyr)
library(plyr)
library(RColorBrewer)
library(multcompView)
library(ggplot2)

#########################################
# read in and modify site trait dataset 
#########################################

full_biomass<-read.csv("./Data/full-biomass-nutrients-06-December-2017.csv")


full_biomass[full_biomass == 'NULL'] <- NA
full_biomass$pct_N <- as.numeric(full_biomass$pct_N)


full_biomass$Ntrt = 0
full_biomass$Ntrt[full_biomass$trt == 'N' | full_biomass$trt == 'NP' | full_biomass$trt == 'NK' | full_biomass$trt == 'NPK' | full_biomass$trt == 'NPK+Fence'] = 1
full_biomass$Ptrt = 0
full_biomass$Ptrt[full_biomass$trt == 'P' | full_biomass$trt == 'NP' | full_biomass$trt == 'PK' | full_biomass$trt == 'NPK' | full_biomass$trt == 'NPK+Fence'] = 1
full_biomass$Ktrt = 0
full_biomass$Ktrt[full_biomass$trt == 'K' | full_biomass$trt == 'NK' | full_biomass$trt == 'PK' | full_biomass$trt == 'NPK' | full_biomass$trt == 'NPK+Fence'] = 1

full_biomass$Ntrt_fac = as.factor(full_biomass$Ntrt)
full_biomass$Ptrt_fac = as.factor(full_biomass$Ptrt)
full_biomass$Ktrt_fac = as.factor(full_biomass$Ktrt)
## AGN = above ground N!
full_biomass$AGN <- full_biomass$mass*(full_biomass$pct_N*0.01)

full_biomass$Nfix = 'no'
full_biomass$Nfix[full_biomass$category == 'LEGUME'] = 'yes'

# add core data to the dataset
core = read.csv('./Data/comb-by-plot-09-April-2018.csv')
core[core == 'NULL'] <- NA

core[, 26:45] <- sapply(core[, 26:45], as.character)
core[, 26:45] <- sapply(core[, 26:45], as.numeric)
core$par_rat = core$Ground_PAR / core$Ambient_PAR

core_data = core[, c(2, 25, 15, 16, 22, 26:46)]
core_info = core[, c(2, 15, 16, 22, 1, 3:14, 17:21, 23, 24)]

biomass_core_data = join(full_biomass, core_data, by = c("site_code", "year", "block", "plot", "trt"), type = 'left', match = 'first')
biomass_core = join(biomass_core_data, core_info, by = c("site_code", "block", "plot", "trt"), type = 'left', match = 'first')

# add SPEI to "leaf_core" data
spei = read.csv('./Data/CRU-annual_2018-07-06.csv')
spei$p_pet = spei$precip / spei$PET

biomass_core_spei = join(biomass_core, spei, by = c("site_code", "year"), type = 'left', match = 'first')

# check
nrow(full_biomass)
nrow(biomass_core) 
nrow(biomass_core_spei
    )

# read in gs climate (>0°C)
tmp_globe = read.csv('./Analysis/climate_gridded/cru_tmp_climExtract_growingseason_globe.csv')
par_globe = read.csv('./Analysis/climate_gridded/cru_par_climExtract_growingseason_globe.csv')
vpd_globe = read.csv('./Analysis/climate_gridded/cru_vpd_climExtract_growingseason_globe.csv')
z_globe =  read.csv('./Analysis/climate_gridded/z_globe.csv')

biomass_core_spei$tmp = NA
biomass_core_spei$par = NA
biomass_core_spei$vpd = NA
biomass_core_spei$z = NA

climate_df = c()
for (i in 1:nrow(biomass_core_spei
                )){
  
  currentLat = biomass_core_spei$latitude[i]
  currentLon = biomass_core_spei$longitude[i]
  
  clim_comb = tmp_globe
  latClim = clim_comb[ , 3]
  lonClim = clim_comb[ , 2]
  
  best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
  best_lat=latClim[best_lat_pos[1]]
  best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
  best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
  
  tmp = subset(clim_comb, lat==best_lat & lon==best_lon)$tmp
  
  clim_comb = par_globe
  latClim = clim_comb[ , 3]
  lonClim = clim_comb[ , 2]
  
  best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
  best_lat=latClim[best_lat_pos[1]]
  best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
  best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
  
  par = subset(clim_comb, lat==best_lat & lon==best_lon)$par
  
  clim_comb = vpd_globe
  latClim = clim_comb[ , 3]
  lonClim = clim_comb[ , 2]
  
  best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
  best_lat=latClim[best_lat_pos[1]]
  best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
  best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
  
  vpd = subset(clim_comb, lat==best_lat & lon==best_lon)$vpd
  
  clim_comb = z_globe
  latClim = clim_comb[ , 3]
  lonClim = clim_comb[ , 2]
  
  best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
  best_lat=latClim[best_lat_pos[1]]
  best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
  best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
  
  z = subset(clim_comb, lat==best_lat & lon==best_lon)$z
  
  climate_df = rbind(climate_df, c(tmp, par, vpd, z, currentLat, currentLon, best_lat, best_lon, i))
  
  
}

# plot(climate_df[,5] ~ climate_df[,7])
# plot(climate_df[,6] ~ climate_df[,8])

biomass_core_spei$tmp = climate_df[, 1]
biomass_core_spei$par = climate_df[, 2]
biomass_core_spei$vpd = climate_df[, 3]
biomass_core_spei$z = climate_df[, 4]

