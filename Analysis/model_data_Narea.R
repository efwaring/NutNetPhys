# Script to compare Narea observed across the Nutrient Network to that predicted
# from climate and LMA

## Note: only works for C3 species

# setwd("/Users/nicksmith/Documents/Git/NutNetPhys/Analysis")

## load libraries
library(R.utils)
library(lme4)
library(car)
# install.packages("devtools")
library(devtools)
devtools::install_github("hohenstein/remef")
library(remef)
library(emmeans)
library(tidyverse)
library(caret)
library(relaimpo)

## load optimal vcmax script (Smith et al., Ecology Letters)

source('optimal_vcmax_R/calc_optimal_vcmax_nutnet.R')
sourceDirectory('optimal_vcmax_R/functions')

## load observational data
traits = read.csv('../Data/leaf_plus.csv')
head(traits)

## fiter by overlapping sites
source('filer_sites.R')

traits = filter(traits, site_code %in% overlap)

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
hist(traits$lma) # some extremely high values
hist(traits$narea) # some extremely high values
hist(traits$SLA)
hist(traits$la_m2)
hist(traits$leaf_dry_mass_g) # some very hig values
hist(subset(traits, leaf_dry_mass_g < 100)$leaf_dry_mass_g)

## some light calculations
traits$lai = -log(traits$Ground_PAR / traits$Ambient_PAR) / 0.86 # from: http://manuals.decagon.com/Manuals/10242_Accupar%20LP80_Web.pdf page 41
hist(traits$lai)
### calculate par per leaf area to assume par absorbed is reduced in dense canopies
traits$par_per_leaf_area = traits$par * ((1 - exp(-0.5 * traits$lai)) / traits$lai) # from Dong et al. (2007) eqn 2
hist(traits$par_per_leaf_area)
hist(traits$Ambient_PAR)
hist(traits$par)
hist(((1 - exp(-0.5 * traits$lai)) / traits$lai))

## remove C4 and values without d13c
traits_sub = subset(traits, lma > 0 & leaf_dry_mass_g < 100 & photosynthetic_pathway == 'C3' & leaf_C13_delta_PDB != 'NA'
                    & leaf_C13_delta_PDB < -20 & trt != 'Fence' & trt != 'NPK+Fence')
traits_sub$delta = (-8 - traits_sub$leaf_C13_delta_PDB) / (1 + traits_sub$leaf_C13_delta_PDB * 0.001)
traits_sub$chi = (traits_sub$delta - 4.4) / (27 - 4.4)
hist(traits_sub$lma) # some extremely high values
hist(traits_sub$narea) # some extremely high values
hist(traits_sub$la_m2)
hist(traits_sub$leaf_dry_mass_g) # 
hist(traits_sub$chi) # looks good

## calculate optimal trait values
traits_photo = calc_optimal_vcmax(tg_c = traits_sub$tmp, paro = traits_sub$par_per_leaf_area, cao = 400, 
                                  vpdo = traits_sub$vpd, z = traits_sub$z,chi = traits_sub$chi)
traits_sub$vcmax25_mod = traits_photo$vcmax25
traits_sub$vcmax25_mod[traits_sub$vcmax25_mod < 0] <- NA
plot(traits_sub$vcmax25_mod ~ traits_sub$tmp)

## calculate N from rubisco and structure
### from Ning Dong and Harrison et al. (2009): Nrubisco (g per m2) = ((vcmax25 * Mr * nr) / (kcat * Nr)) * Mn
### vcmax25 = vcmax at 25C (umol m-2 s-1)
### Nr = number of catalytic sites per mole of Rubisco = 8 mol Rubisco site per mol Rubisco
### kcat = catalytic turnover at 25C = 3.5 CO2 per mol Rubisco site per second
### Mr = molecular mass of Rubisco = 0.55 g Rubisco per umol Rubisco (550000 g Rubisco per mol Rubisco)
### nr = N concentration of Rubisco = 11400 mol N per g Rubisco
### Mn = molecular mass of N = 14 gN per mol N
# https://onlinelibrary.wiley.com/doi/full/10.1111/j.1365-3040.2008.01918.x
### lma conversion from Dong et al. (2017) Biogeosciences
traits_sub$n_rubisco_mod = traits_sub$vcmax25_mod * (1/3500000) * (1/8) * 550000 * .0114 * 14
traits_sub$n_structure_mod = (10^-2.67) * (traits_sub$lma ^ 0.99)
hist(traits_sub$n_rubisco_mod)
hist(traits_sub$n_structure_mod)

## fit some models for N following Dong et al. (2017) but also include fertilization
### all individuals
n_pred_all_lm = lm(log(narea) ~ chi + log(par_per_leaf_area) + tmp + log(lma) + Nfix + trt, data = traits_sub)
summary(n_pred_all_lm) # not much of a chi response
anova(n_pred_all_lm)
plot(resid(n_pred_all_lm) ~ fitted(n_pred_all_lm))
# plot(log(narea) ~ chi + log(par_per_leaf_area) + tmp + log(lma) + Nfix + trt, data = traits_sub)
emmeans(n_pred_all_lm, ~Nfix)
cld(emmeans(n_pred_all_lm, ~trt)) # N is highest
calc.relimp(n_pred_all_lm, rela = T) # lma = 60%, tmp = 26%, par = 1% chi = 12%, trt = 0.2%, Nfix = 2%

### community mean
traits_sub_group_by_site = group_by(subset(traits_sub, Nfix == 'no'), site_code, trt)
traits_sub_site = summarise(traits_sub_group_by_site, 
                            narea_mean = mean(narea, na.rm = T), chi_mean = mean(chi, na.rm = T),
                            tmp_mean = mean(tmp, na.rm = T), lma_mean = mean(lma, na.rm = T),
                            par_per_leaf_area_mean = mean(par_per_leaf_area, na.rm = T),
                            n_rubisco_mod_mean = mean(n_rubisco_mod, na.rm = T),
                            n_structure_mod_mean = mean(n_structure_mod, na.rm = T))
n_pred_com_lm = lm(log(narea_mean) ~ chi_mean + log(par_per_leaf_area_mean) + tmp_mean + log(lma_mean) + trt,
                   data = traits_sub_site)
summary(n_pred_com_lm) # weak chi and par response, but others as expected
anova(n_pred_com_lm)
plot(resid(n_pred_com_lm) ~ fitted(n_pred_com_lm))
# plot(log(narea_mean) ~ chi_mean + log(par_per_leaf_area_mean) + tmp_mean + log(lma_mean) + trt,
#    data = traits_sub_site)
cld(emmeans(n_pred_com_lm, ~trt))
calc.relimp(n_pred_com_lm, rela = T) # lma = 73%, chi = 16%, temp = 9%, par = 2%, treatment = 0.4%

### examining the different components
# n_pred_model_lm = lm(log(narea) ~ (log(n_rubisco_mod) + log(n_structure_mod)) * Nfix, data = traits_sub)
n_pred_model_lm = lm(log(narea) ~ log(n_rubisco_mod) + log(n_structure_mod) + trt, data = subset(traits_sub, Nfix == 'no'))
Anova(n_pred_model_lm)
summary(n_pred_model_lm)
plot(resid(n_pred_model_lm) ~ fitted(n_pred_model_lm))
# plot(log(narea) ~ log(n_rubisco_mod) + log(n_structure_mod) + trt, data = subset(traits_sub, Nfix == 'no'))
# test(emtrends(n_pred_model_lm, ~1, var = 'log(n_rubisco_mod)', at = list(Nfix = 'no')))
# test(emtrends(n_pred_model_lm, ~1, var = 'log(n_rubisco_mod)', at = list(Nfix = 'yes')))
# test(emtrends(n_pred_model_lm, ~1, var = 'log(n_structure_mod)', at = list(Nfix = 'no')))
# test(emtrends(n_pred_model_lm, ~1, var = 'log(n_structure_mod)', at = list(Nfix = 'yes')))
calc.relimp(n_pred_model_lm, rela = T) # 81% structure, 19% rubisco, treatment = 0.03%
plot(log(narea) ~ log(n_rubisco_mod) , data = subset(traits_sub, Nfix == 'no'))
plot(log(narea) ~ log(n_structure_mod) , data = subset(traits_sub, Nfix == 'no'))

### components at the community level
n_pred_model_com_lm = lm(log(narea_mean) ~ log(n_rubisco_mod_mean) + log(n_structure_mod_mean) + trt, data = traits_sub_site)
Anova(n_pred_model_com_lm)
summary(n_pred_model_com_lm)
plot(resid(n_pred_model_com_lm) ~ fitted(n_pred_model_com_lm))
# plot(log(n_pred_model_com_lm) ~ (log(n_pred_model_com_lm) + log(n_structure_mod_mean)), data = traits_sub_site)
calc.relimp(n_pred_model_com_lm, rela = T) # 95% structure, 5% rubisco, 0.5% treatment

#############
# okay do the same thing for %N >> stronger treatment effect
############
nper_pred_all_lm = lm(log(leaf_pct_N) ~ chi + log(par_per_leaf_area) + tmp + log(lma) + Nfix + trt, data = traits_sub)
summary(nper_pred_all_lm) # 
anova(nper_pred_all_lm)
plot(resid(nper_pred_all_lm) ~ fitted(nper_pred_all_lm))
calc.relimp(nper_pred_all_lm, rela = T) # N fixer most important, followed by temperature, treatment = 10%


### components
traits_sub$sla_m2 = traits_sub$la_m2 / traits_sub$leaf_dry_mass_g
traits_sub$n_rubisco_mod_pct = traits_sub$n_rubisco_mod * traits_sub$sla_m2
traits_sub$n_structure_mod_pct = traits_sub$n_structure_mod * traits_sub$sla_m2

n_pred_model_lm_pct = lm(log(leaf_pct_N) ~ log(n_rubisco_mod_pct) + log(n_structure_mod_pct) + trt, data = subset(traits_sub, Nfix == 'no'))
Anova(n_pred_model_lm_pct)
summary(n_pred_model_lm_pct)
plot(resid(n_pred_model_lm_pct) ~ fitted(n_pred_model_lm_pct))
calc.relimp(n_pred_model_lm_pct, rela = T) # 
plot(log(leaf_pct_N) ~ log(n_rubisco_mod_pct) , data = subset(traits_sub, Nfix == 'no'))
plot(log(leaf_pct_N) ~ log(n_structure_mod_pct) , data = subset(traits_sub, Nfix == 'no'))


