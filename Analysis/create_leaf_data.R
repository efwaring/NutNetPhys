# create leaf trait dataset for analyses

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
# read in and modify foliar trait dataset (Jen Firn)
#########################################

leaf = read.csv('./Data/NutNet-foliar-traits-7JAN2017.csv')

leaf[leaf == 'NULL'] <- NA

leaf[, 8:30] <- sapply(leaf[, 8:30], as.character)
leaf[, 8:30] <- sapply(leaf[, 8:30], as.numeric)

leaf$Ntrt = 0
leaf$Ntrt[leaf$trt == 'N' | leaf$trt == 'NP' | leaf$trt == 'NK' | leaf$trt == 'NPK' | leaf$trt == 'NPK+Fence'] = 1
leaf$Ptrt = 0
leaf$Ptrt[leaf$trt == 'P' | leaf$trt == 'NP' | leaf$trt == 'PK' | leaf$trt == 'NPK' | leaf$trt == 'NPK+Fence'] = 1
leaf$Ktrt = 0
leaf$Ktrt[leaf$trt == 'K' | leaf$trt == 'NK' | leaf$trt == 'PK' | leaf$trt == 'NPK' | leaf$trt == 'NPK+Fence'] = 1

leaf$Ntrt_fac = as.factor(leaf$Ntrt)
leaf$Ptrt_fac = as.factor(leaf$Ptrt)
leaf$Ktrt_fac = as.factor(leaf$Ktrt)

leaf$Narea = leaf$leaf_pct_N * (1/leaf$SLA)

leaf$Nfix = 'no'
leaf$Nfix[leaf$Family == 'Fabaceae'] = 'yes'

# add core data to the dataset
core = read.csv('./Data/comb-by-plot-09-April-2018.csv')
core[core == 'NULL'] <- NA

core[, 26:45] <- sapply(core[, 26:45], as.character)
core[, 26:45] <- sapply(core[, 26:45], as.numeric)
core$par_rat = core$Ground_PAR / core$Ambient_PAR

core_data = core[, c(2, 25, 15, 16, 22, 26:46)]
core_info = core[, c(2, 15, 16, 22, 1, 3:14, 17:21, 23, 24)]

leaf_core_data = join(leaf, core_data, by = c("site_code", "year", "block", "plot", "trt"), type = 'left', match = 'first')
leaf_core = join(leaf_core_data, core_info, by = c("site_code", "block", "plot", "trt"), type = 'left', match = 'first')

# add SPEI to "leaf_core" data
spei = read.csv('./Data/CRU-annual_2018-07-06.csv')
spei$p_pet = spei$precip / spei$PET

leaf_core_spei = join(leaf_core, spei, by = c("site_code", "year"), type = 'left', match = 'first')

# check
nrow(leaf)
nrow(leaf_core) # some core data matches the leaf data twice...
nrow(leaf_core_spei)

# read in gs climate (>0°C)
tmp_globe = read.csv('/Users/nicksmith/Dropbox/0Main/Research/Colimitation/Spatial_Maps/cru_tmp_climExtract_growingseason_globe.csv')
par_globe = read.csv('/Users/nicksmith/Dropbox/0Main/Research/Colimitation/Spatial_Maps/cru_par_climExtract_growingseason_globe.csv')
vpd_globe = read.csv('/Users/nicksmith/Dropbox/0Main/Research/Colimitation/Spatial_Maps/cru_vpd_climExtract_growingseason_globe.csv')
z_globe =  read.csv('/Users/nicksmith/Dropbox/0Main/Research/Colimitation/Spatial_Maps/z_globe.csv')

leaf_core_spei$tmp = NA
leaf_core_spei$par = NA
leaf_core_spei$vpd = NA
leaf_core_spei$z = NA

climate_df = c()
for (i in 1:nrow(leaf_core_spei)){
  
  currentLat = leaf_core_spei$latitude[i]
  currentLon = leaf_core_spei$longitude[i]
  
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

leaf_core_spei$tmp = climate_df[, 1]
leaf_core_spei$par = climate_df[, 2]
leaf_core_spei$vpd = climate_df[, 3]
leaf_core_spei$z = climate_df[, 4]

########################################################
# add species information
########################################################
info = read.csv('../Data/full-species-info-25-January-2019.csv')
code_taxon_data = paste(toupper(leaf_core_spei$site_code), toupper(leaf_core_spei$Taxon), sep = ' ')

info_caps = mutate_all(info, .funs = toupper)
info_caps$code_taxon = paste(info_caps$site_code, info_caps$standard_taxon, sep = ' ')

n_info = NULL
for (i in 1:length(code_taxon_data)){
  
  ancillary_data = subset(info_caps, code_taxon == code_taxon_data[i])[, c(6:10)]
  
  ancillary_data$n = i
  
  n_info = rbind(n_info, ancillary_data)
  
}

leaf_core_spei_info = cbind(leaf_core_spei, n_info)

# write.csv(leaf_core_spei_info, "../Data/leaf_plus.csv")