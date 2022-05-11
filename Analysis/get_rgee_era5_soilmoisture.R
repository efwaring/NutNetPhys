# in data ####
rm(list=ls()) # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
dev.off()     # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~
cat("\f")     # ~ ~ ~ ~ ~ ~ ~ ~ ~ ~ ~

# library(doBy)
# library(raster)
# library(sf)
# library(rgee)
# library(tmap)
# library(tidyverse)
# library(lubridate)
# ee_check()
# ee_check_python()
# ee_check_credentials()
# ee_check_python_packages()
# ee_Initialize(drive = T, gcs = T)

# target data ####
# https://developers.google.com/earth-engine/datasets/catalog/ECMWF_ERA5_LAND_HOURLY

setwd("~/Git/NutNetPhys/Analysis/")
env<-read.csv("growing_season_data/Weather_site_inventory_20200921_rev.csv")
df<-read.csv("../Data/processed/traits_v3.csv")

df_year<-df[,c(3:4)]
df_year_groupby<-group_by(df_year, site_code)
df_year_summarize<-summarize(df_year_groupby, year = mean(year))

new_no_coords<-left_join(env, df_year_summarize)
new_no_coords$site_long[new_no_coords$site_code=="smith.us"]=-122.545400
new_no_coords$site_lat[new_no_coords$site_code=="smith.us"]=48.046259
new = st_as_sf(new_no_coords, coords = c("site_long", "site_lat"), crs = 4326)

new$start_mo_with_length<-new$gs_start+new$gs_len
new$end_month<-ifelse(new$start_mo_with_length>12, new$start_mo_with_length-12,new$start_mo_with_length)
new$start_year<-new$year
new$end_year<-ifelse(new$start_mo_with_length>12, new$start_year+1, new$start_year)
new$start_month<-new$gs_start

new$start_day<-1
new$end_day<-28

new$start_date<-with(new, paste(start_year, start_day, start_month, sep = "-"))
new$start_date<-as.Date(new$start_date, format = "%Y-%d-%m")

new$end_date<-with(new, paste(end_year, end_day, end_month, sep = "-"))
new$end_date<-as.Date(new$end_date, format = "%Y-%d-%m")

new<-subset(new, year != "NA")
# tmap_mode("view")
# tm_shape(new)+
#   tm_dots(shape = 16, size = 0.1, col = "black")

# output<-data.frame(site_code = as.numeric(),
#                    daterange = as.numeric(),
#                    soil_moisture = as.numeric())

# this variable for looping is a little inefficient
# but it prevents errors resulting from multiple sites with
# the same growing season start and end dates
new$ID<-with(new, interaction(start_date,site_code,end_date))

# new$end_date<-new$start_date+3 # for testing short time periods
# new<-subset(new, site_code == "shps.us") # for testing a single site

rm(list=setdiff(ls(), "new"))

# start getting sites data ####

# bldr.us ####

levels(as.factor(new$site_code))
bldr.us<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "bldr.us")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>%
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>%
  as.data.frame()

# rbind dummy data to output
bldr.us<-rbind(bldr.us, ee_extracted_values_long)
write.csv(bldr.us, "bldr_us.csv", row.names = F)

# Sys.sleep(300)

# bnch.us ####

levels(as.factor(new$site_code))
bnch.us<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "bnch.us")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>%
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>%
  as.data.frame()

# rbind dummy data to output
bnch.us<-rbind(bnch.us, ee_extracted_values_long)
write.csv(bnch.us, "bnch_us.csv", row.names = F)

# Sys.sleep(300)

# bogong.au ####

levels(as.factor(new$site_code))
bogong.au<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "bogong.au")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>%
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>%
  as.data.frame()

# rbind dummy data to output
bogong.au<-rbind(bogong.au, ee_extracted_values_long)
write.csv(bogong.au, "bogong_au.csv", row.names = F)

# Sys.sleep(300)

# burrawan.au ####

levels(as.factor(new$site_code))
burrawan.au<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "burrawan.au")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>%
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>%
  as.data.frame()

# rbind dummy data to output
burrawan.au<-rbind(burrawan.au, ee_extracted_values_long)
write.csv(burrawan.au, "burrawan_au.csv", row.names = F)

# Sys.sleep(300)

# cbgb.us ####

levels(as.factor(new$site_code))
cbgb.us<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "cbgb.us")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
cbgb.us<-rbind(cbgb.us, ee_extracted_values_long)
write.csv(cbgb.us, "cbgb_us.csv", row.names = F)

# Sys.sleep(300)

# comp.pt ####

levels(as.factor(new$site_code))
comp.pt<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "comp.pt")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
comp.pt<-rbind(comp.pt, ee_extracted_values_long)
write.csv(comp.pt, "comp_pt.csv", row.names = F)
# Sys.sleep(300)
# cowi.ca ####

levels(as.factor(new$site_code))
cowi.ca<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "cowi.ca")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
cowi.ca<-rbind(cowi.ca, ee_extracted_values_long)
write.csv(cowi.ca, "cowi_ca.csv", row.names = F)

# Sys.sleep(300)
# elliot.us ####

levels(as.factor(new$site_code))
elliot.us<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "elliot.us")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
elliot.us<-rbind(elliot.us, ee_extracted_values_long)
write.csv(elliot.us, "elliot_us.csv", row.names = F)

# Sys.sleep(300)
# frue.ch ####

levels(as.factor(new$site_code))
frue.ch<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "frue.ch")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
frue.ch<-rbind(frue.ch, ee_extracted_values_long)
write.csv(frue.ch, "frue_ch.csv", row.names = F)

# Sys.sleep(300)
# gilb.za ####

levels(as.factor(new$site_code))
gilb.za<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "gilb.za")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
gilb.za<-rbind(gilb.za, ee_extracted_values_long)
write.csv(gilb.za, "gilb_za.csv", row.names = F)

# Sys.sleep(300)
# hopl.us ####

levels(as.factor(new$site_code))
hopl.us<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "hopl.us")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
hopl.us<-rbind(hopl.us, ee_extracted_values_long)
write.csv(hopl.us, "hopl_us.csv", row.names = F)
# Sys.sleep(300)
# jena.de ####

levels(as.factor(new$site_code))
jena.de<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "jena.de")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
jena.de<-rbind(jena.de, ee_extracted_values_long)
write.csv(jena.de, "jena_de.csv", row.names = F)
# Sys.sleep(300)
# kiny.au ####

levels(as.factor(new$site_code))
kiny.au<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "kiny.au")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
kiny.au<-rbind(kiny.au, ee_extracted_values_long)
write.csv(kiny.au, "kiny_au.csv", row.names = F)
# Sys.sleep(300)
# konz.us ####

levels(as.factor(new$site_code))
konz.us<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "konz.us")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
konz.us<-rbind(konz.us, ee_extracted_values_long)
write.csv(konz.us, "konz_us.csv", row.names = F)
# Sys.sleep(300)
# lancaster.uk ####

levels(as.factor(new$site_code))
lancaster.uk<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "lancaster.uk")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
lancaster.uk<-rbind(lancaster.uk, ee_extracted_values_long)
write.csv(lancaster.uk, "lancaster_uk.csv", row.names = F)
# Sys.sleep(300)
# look.us ####

levels(as.factor(new$site_code))
look.us<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "look.us")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
look.us<-rbind(look.us, ee_extracted_values_long)
write.csv(look.us, "look_us.csv", row.names = F)
# Sys.sleep(300)
# mcla.us ####

levels(as.factor(new$site_code))
mcla.us<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "mcla.us")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
mcla.us<-rbind(mcla.us, ee_extracted_values_long)
write.csv(mcla.us, "mcla_us.csv", row.names = F)
# Sys.sleep(300)
# mtca.au ####

levels(as.factor(new$site_code))
mtca.au<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "mtca.au")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
mtca.au<-rbind(mtca.au, ee_extracted_values_long)
write.csv(mtca.au, "mtca_au.csv", row.names = F)
# Sys.sleep(300)
# sage.us ####

levels(as.factor(new$site_code))
sage.us<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "sage.us")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
sage.us<-rbind(sage.us, ee_extracted_values_long)
write.csv(sage.us, "sage_us.csv", row.names = F)
# Sys.sleep(300)
# saline.us ####

levels(as.factor(new$site_code))
saline.us<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "saline.us")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
saline.us<-rbind(saline.us, ee_extracted_values_long)
write.csv(saline.us, "saline_us.csv", row.names = F)
# Sys.sleep(300)
# sgs.us ####

levels(as.factor(new$site_code))
sgs.us<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "sgs.us")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
sgs.us<-rbind(sgs.us, ee_extracted_values_long)
write.csv(sgs.us, "sgs_us.csv", row.names = F)
# Sys.sleep(300)
# shps.us ####

levels(as.factor(new$site_code))
shps.us<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "shps.us")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
shps.us<-rbind(shps.us, ee_extracted_values_long)
write.csv(shps.us, "shps_us.csv", row.names = F)
# Sys.sleep(300)
# sier.us ####

levels(as.factor(new$site_code))
sier.us<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "sier.us")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
sier.us<-rbind(sier.us, ee_extracted_values_long)
write.csv(sier.us, "sier_us.csv", row.names = F)
# Sys.sleep(300)
# smith.us ####

levels(as.factor(new$site_code))
smith.us<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "smith.us")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
smith.us<-rbind(smith.us, ee_extracted_values_long)
write.csv(smith.us, "smith_us.csv", row.names = F)
# Sys.sleep(300)
# summ.za ####

levels(as.factor(new$site_code))
summ.za<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "summ.za")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
summ.za<-rbind(summ.za, ee_extracted_values_long)
write.csv(summ.za, "summ_za.csv", row.names = F)
# Sys.sleep(300)
# unc.us ####

levels(as.factor(new$site_code))
unc.us<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "unc.us")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
unc.us<-rbind(unc.us, ee_extracted_values_long)
write.csv(unc.us, "unc_us.csv", row.names = F)
# Sys.sleep(300)
# valm.ch ####

levels(as.factor(new$site_code))
valm.ch<-data.frame(site_code = as.numeric(),
                    daterange = as.numeric(),
                    soil_moisture = as.numeric())
df_temp<-subset(new, site_code == "valm.ch")

# create function for extracting SMAP data
IC <- ee$ImageCollection("ECMWF/ERA5_LAND/HOURLY") %>% # select SMAP
  ee$ImageCollection$filterDate(as.character(unique(df_temp$start_date)), as.character(unique(df_temp$end_date))) %>% # index start and end dates from dummy data
  ee$ImageCollection$map(function(x) x$select("volumetric_soil_water_layer_1")) %>% # select band
  ee$ImageCollection$toBands() # from imagecollection to image

# extract values for sites
ee_extracted_values <- ee_extract(x = IC, y = df_temp["site_code"], sf = FALSE, fun = ee$Reducer$mean(), via = "drive") ## fun = mean

# convert from wide format to long format
ee_extracted_values_long = ee_extracted_values %>% 
  pivot_longer(-site_code, names_to = "daterange", values_to = "soil_moisture") %>% 
  as.data.frame()

# rbind dummy data to output
valm.ch<-rbind(valm.ch, ee_extracted_values_long)
write.csv(valm.ch, "valm_ch.csv", row.names = F)
# Sys.sleep(300)
