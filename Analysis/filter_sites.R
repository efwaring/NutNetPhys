## File to determine overlapping sites between the two data sets
library(dplyr)
library(tidyr)

source("./Analysis/create_biomass_data.R")
source("./Analysis/create_leaf_data.R")


overlap <- intersect(biomass_core_spei$site_code, leaf_core_spei$site_code)

biomass_core_spei <- filter(biomass_core_spei, site_code %in% overlap) 

leaf_core_spei <- filter(biomass_core_spei, site_code %in% overlap) 

full_biomass <- filter(full_biomass, site_code %in% overlap)

leaf <- filter(leaf, site_code %in% overlap)

#test that they all match

biomass_core_spei$site_code %in% overlap

