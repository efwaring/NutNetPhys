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
traits$Ntrt_fac = as.factor(traits$Ntrt_fac)
traits$Ptrt_fac = as.factor(traits$Ptrt_fac)
traits$Ktrt_fac = as.factor(traits$Ktrt_fac)

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
traits$delta = (-8 - traits$leaf_C13_delta_PDB) / (1 + traits$leaf_C13_delta_PDB * 0.001)
traits$chi = (traits$delta - 4.4) / (27 - 4.4)
traits$chi[traits$chi > 1] = 1
traits$chi[traits$chi < 0] = 0
# hist(traits$lma) # some extremely high values
# hist(traits$narea) # some extremely high values
# hist(traits$la_m2)
# hist(traits$leaf_dry_mass_g) # 
# hist(traits$chi) # looks good
traits = subset(traits, chi > 0.4 & chi < 1)

write.csv(traits, '../Data/processed/traits.csv')

traits_group_by_site = group_by(traits, 
                                    site_code, plot, block, trt, Ntrt_fac, Ptrt_fac, Ktrt_fac)
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

write.csv(traits_sub_site, '../Data/processed/traits_site.csv')