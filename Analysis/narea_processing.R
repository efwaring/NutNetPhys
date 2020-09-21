# Narea preprocessing

#########################################
# packages necessary
#########################################

library(lme4)
library(car)
library(emmeans)
library(dplyr)
library(RColorBrewer)
library(multcompView)
library(ggplot2)

#########################################
# read in and modify foliar trait dataset (Jen Firn)
#########################################
## load observational data
traits = read.csv('../Data/leaf_plus.csv')
head(traits)

## fiter by overlapping sites
# source('filter_sites.R')
# traits = filter(traits, site_code %in% overlap)

## make sure factors are correctly defined
traits$Ntrt_fac = as.factor(traits$Ntrt)
traits$Ptrt_fac = as.factor(traits$Ptrt)
traits$Ktrt_fac = as.factor(traits$Ktrt)

## add in photosynthetic pathway information
levels(traits$Family) # check Amaranthaceae, Asteraceae, Boraginaceae, Caryophyllaceae, Cyperaceae, Euphorbiaceae,
# Polygonaceae, Poaceae, Scrophulariaceae
# only one!
# C4
traits$photosynthetic_pathway[traits$photosynthetic_pathway == 'NULL'
                              & traits$Family == 'Cyperaceae' & traits$Taxon == 'FIMBRISTYLIS DICHOTOMA'] <- 'C4'

traits$photosynthetic_pathway[traits$photosynthetic_pathway == 'NULL'] <- 'C3'


## calculate lma (in g m-2) and narea
hist(traits$leaf_area_mm2)
traits$la_m2 = traits$leaf_area_mm2 * (1/1000000)
traits$lma = traits$leaf_dry_mass_g / traits$la_m2
traits$narea = (traits$leaf_pct_N / 100) * (traits$lma)
# hist(traits$lma) # some extremely high values
# hist(traits$narea) # some extremely high values
# hist(traits$SLA)
# hist(traits$la_m2)
# hist(traits$leaf_dry_mass_g) # some very hig values
# hist(subset(traits, leaf_dry_mass_g < 100)$leaf_dry_mass_g)

## some light calculations
traits$lai = -log(traits$Ground_PAR / traits$Ambient_PAR) / 0.86 # from: http://manuals.decagon.com/Manuals/10242_Accupar%20LP80_Web.pdf page 41
# hist(traits$lai)
### calculate par per leaf area to assume par absorbed is reduced in dense canopies
traits$par_per_leaf_area = traits$par * ((1 - exp(-0.5 * traits$lai)) / traits$lai) # from Dong et al. (2007) eqn 2
# hist(traits$par_per_leaf_area)
# hist(traits$Ambient_PAR)
# hist(traits$par)
# hist(((1 - exp(-0.5 * traits$lai)) / traits$lai))

## remove C4 and values without d13c
traits$delta = ((-0.008 - traits$leaf_C13_delta_PDB * 0.001) / (1 + traits$leaf_C13_delta_PDB * 0.001)) * 1000
traits$chi[traits$photosynthetic_pathway == 'C3'] = 
  (traits$delta[traits$photosynthetic_pathway == 'C3'] * 0.001 - 0.0044) / (0.027 - 0.0044)
traits$chi[traits$photosynthetic_pathway == 'C4'] = 
  (traits$delta[traits$photosynthetic_pathway == 'C4'] * 0.001 - 0.0044) / ((-0.0057 + 0.03*0.4) - 0.0044)
traits$chi[traits$chi > 1] = 1
traits$chi[traits$chi < 0] = 0
# hist(traits$lma) # some extremely high values
# hist(traits$narea) # some extremely high values
# hist(traits$la_m2)
# hist(traits$leaf_dry_mass_g) # 
# hist(traits$chi) # looks good
# traits = subset(traits, chi > 0.2 & chi < 0.95)

traits_group_by_site_for_count = group_by(traits, site_code)
traits_site_count = summarise(traits_group_by_site_for_count, n = n())
subset(traits_site_count, n <=30)

# subset(traits, site_code == 'bldr.us')
# traits = subset(traits, site_code != 'bldr.us' & 
#                   site_code != 'gilb.za' & 
#                   site_code != 'konz.us' &
#                   site_code != 'saline.us') # remove sites without much leaf data (mostly due to bad chi values)

## add in species composition data
cover = read.csv('../Data/full-cover-31-August-2020.csv')
head(cover)

cover_select = select(cover, site_code, plot, year, Taxon, max_cover)

traits_w_cover = left_join(traits, cover_select)
traits_w_cover$spp_lai = traits_w_cover$lai * (traits_w_cover$max_cover/100)

# write.csv(traits_w_cover, '../Data/processed/traits.csv')

traits_group_by_site = group_by(traits, 
                                    site_code, plot, block, trt, Ntrt, Ptrt, Ktrt)
traits_site = summarise(traits_group_by_site, 
                            narea_mean = mean(narea, na.rm = T), 
                            chi_mean = mean(chi, na.rm = T),
                            tmp_mean = mean(tmp, na.rm = T),
                            vpd_mean = mean(vpd, na.rm = T),
                            z_mean = mean(z, na.rm = T),
                            lma_mean = mean(lma, na.rm = T),
                            par_per_leaf_area_mean = mean(par_per_leaf_area, na.rm = T),
                            par_mean = mean(par, na.rm = T),
                            lai_mean = mean(lai, na.rm = T),
                            p_pet_mean = mean(p_pet, na.rm = T),
                            live_mass_mean = mean(live_mass, na.rm = T))

# write.csv(traits_site, '../Data/processed/traits_site.csv')
