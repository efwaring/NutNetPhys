# create the growing season data for each site in the leaf trait dataset

## libraries
library(R.matlab)
library(tidyverse)

## read in leaf data
leaf_data <- read.csv('../../Data/foliar_cover_updated_2.csv')
colnames(leaf_data)
leaf_data_site <- leaf_data[, c(1, 3, 19, 20)]
leaf_data_site_group_by <- group_by(leaf_data_site, site_code)
leaf_data_site_year <- summarise(leaf_data_site_group_by, 
                                 year_leaf = mean(year, na.rm = T),
                                 latitude = mean(latitude, na.rm = T),
                                 longitude = mean(longitude, na.rm = T))

## read in growing season information data and combine with site_year information
gs_information <- read.csv('Weather_site_inventory_20200921.csv')
head(gs_information)

leaf_data_site_year_gs <- left_join(leaf_data_site_year, gs_information)
head(leaf_data_site_year_gs)

## set path to cru data
cru_path <- '~/Documents/Data/CRU3.24'

## create start and end CRU columns for month and year for each site
leaf_data_site_year_gs$start_mo_CRU_gs <- leaf_data_site_year_gs$gs_start_month
leaf_data_site_year_gs$end_mo_CRU_gs <- leaf_data_site_year_gs$gs_start_month + leaf_data_site_year_gs$gs_len
leaf_data_site_year_gs$end_yr_CRU_gs <- leaf_data_site_year_gs$year_leaf
leaf_data_site_year_gs$start_yr_CRU_gs <- leaf_data_site_year_gs$end_yr_CRU_gs - 29

## add start and end month data for sites without it (fix later?)
leaf_data_site_year_gs[leaf_data_site_year_gs$site_code == 'gilb.za',]$start_mo_CRU_gs <- 10
leaf_data_site_year_gs[leaf_data_site_year_gs$site_code == 'gilb.za',]$end_mo_CRU_gs <- 4
leaf_data_site_year_gs[leaf_data_site_year_gs$site_code == 'summ.za',]$start_mo_CRU_gs <- 9
leaf_data_site_year_gs[leaf_data_site_year_gs$site_code == 'summ.za',]$end_mo_CRU_gs <- 5

## fix 

## loop to create climatologies for each site
vars <- c('pre2','tmp2','vpd2','par','alpha')
varsFileN <- c('pre','tmp','vpd','par','alpha')
sites <- leaf_data_site_year_gs$site_code

latlonClim=readMat('~/Documents/Data/CRU3.24/pre2/cru_pre_1901.mat')[[1]][,1:2]
latClim=latlonClim[,2]
lonClim=latlonClim[,1]
latlonClim_df <- as.data.frame(latlonClim)
colnames(latlonClim_df)=c('lon','lat')

gs_climate <- c()
for(sitei in 1:length(sites)){
  
  site_temp <- sites[sitei]
  site_startyear <- leaf_data_site_year_gs[leaf_data_site_year_gs$site_code == site_temp,]$start_yr_CRU_gs
  site_endyear <- leaf_data_site_year_gs[leaf_data_site_year_gs$site_code == site_temp,]$end_yr_CRU_gs
  site_startmonth <- leaf_data_site_year_gs[leaf_data_site_year_gs$site_code == site_temp,]$start_mo_CRU_gs
  site_endmonth <- leaf_data_site_year_gs[leaf_data_site_year_gs$site_code == site_temp,]$end_mo_CRU_gs
  site_lat <- leaf_data_site_year_gs[leaf_data_site_year_gs$site_code == site_temp,]$latitude
  site_lon <- leaf_data_site_year_gs[leaf_data_site_year_gs$site_code == site_temp,]$longitude
  
  best_lat_pos <- which(abs(latClim - site_lat)==min(abs(latClim - site_lat)))
  best_lat <- latClim[best_lat_pos[1]]
  best_lon_pos <- which(abs(subset(latlonClim_df, lat == best_lat)$lon - site_lon)==min(abs(subset(latlonClim_df, lat == best_lat)$lon - site_lon)))
  best_lon <- subset(latlonClim_df, lat== best_lat)$lon[best_lon_pos[1]]
  
  site_gs_climate <- c()
  for(vari in 1:length(vars)){
    
    site_gs_climate_year <- c()
    for(yeari in c(site_startyear:site_endyear)){
      
      latlon_var <- subset(as.data.frame(readMat(paste('~/Documents/Data/CRU3.24/',vars[vari],'/cru_',varsFileN[vari],'_',yeari,'.mat',sep=''))[[1]]),
                           V1 == best_lon & V2 == best_lat)
      latlon_var_plus <- cbind(latlon_var, latlon_var[,3:14]) # add in second year for southern hemisphere sites
      latlon_var_gs <- latlon_var_plus[,(3 + site_startmonth):(3 + site_endmonth)]
      latlon_var_mean <- rowMeans(latlon_var_gs, dim=1, na.rm=T)
      
      site_gs_climate_year <- c(site_gs_climate_year, latlon_var_mean)
      
    }
    
    site_gs_climate_var <- mean(site_gs_climate_year, na.rm = T)
    site_gs_climate <- c(site_gs_climate, site_gs_climate_var)
    
  }
  gs_climate <- rbind(gs_climate, site_gs_climate)
  
}

gs_climate

gs_climate_df <- as.data.frame(gs_climate)
colnames(gs_climate_df) <- c('pre2_gs', 'tmp2_gs', 'vpd2_gs', 'par2_gs', 'alpha')
gs_climate_df$site_code <- sites
gs_climate_df

write.csv(gs_climate_df,'cru_gs_climate.csv')




