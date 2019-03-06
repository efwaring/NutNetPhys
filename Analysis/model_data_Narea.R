# Script to compare Narea observed across the Nutrient Network to that predicted
# from climate and LMA

## Note: only works for C3 species

## load libraries
library(R.utils)

## load optimal vcmax script (Smith et al., Ecology Letters)
sourceDirectory('optimal_vcmax_R/functions')
source('optimal_vcmax_R/calc_optimal_vcmax_nutnet.R')

## load observational data
traits = read.csv('../Data/leaf_plus.csv')
head(traits)

## calculate lma (in g m-2) and narea
traits$lma = 1 / (traits$SLA * (1 / 1000000)) # conver mm2 g-1 to g m-2
traits$narea = (traits$leaf_pct_N / 100) * (traits$lma)
hist(traits$lma) # some extremely high values
hist(traits$narea) # some extremely high values

## remove C4 and values without d13c
traits_sub = subset(traits, lma < 100000 & photosynthetic_pathway == 'C3' & leaf_C13_delta_PDB != 'NA'
                    & leaf_C13_delta_PDB < -20)

## calculate optimal trait values
traits_photo = calc_optimal_vcmax(tg_c = traits_sub$tmp, paro = traits_sub$par, cao = 400, 
                                  vpdo = traits_sub$vpd, z = traits_sub$z, dleaf = traits_sub$leaf_C13_delta_PDB)
traits_sub$vcmax_mod = as.numeric(as.character(traits_photo$vcmax_prime))
traits_sub$vcmax25_mod = traits_sub$vcmax_mod / calc_tresp_mult(traits_sub$tmp, traits_sub$tmp, 25)
traits_sub$vcmax25_mod[traits_sub$vcmax25_mod < 0] <- NA

## calculate N from rubisco and structure
### specific activity of rubisco is 60 umol CO2 gRubisco-1 s-1
### mass ratio of rubisco is 7.16 gRubisco gNRubisco-1
### lma conversion from Dong et al. (2017) Biogeosciences
traits_sub$n_rubisco_mod = traits_sub$vcmax25_mod * (1/60) * (1/7.16)
traits_sub$n_structure_mod = (10^-2.67) * (traits_sub$lma ^ 0.99)
hist(traits_sub$n_rubisco_mod)
hist(traits_sub$n_structure_mod)

## compare modeled and observed N
### first remove very high lma (>1000 g m-2)
n_pred_lmer = lmer(log(narea) ~ log(n_rubisco_mod) + log(n_structure_mod) + (1|site_code), data = traits_sub)
Anova(n_pred_lmer)
summary(n_pred_lmer)
plot(resid(n_pred_lmer) ~ fitted(n_pred_lmer))
# cld(emmeans(n_pred_lmer, ~N))
test(emtrends(n_pred_lmer, ~1, var = "n_rubisco_mod"))
test(emtrends(n_pred_lmer, ~1, var = "n_structure_mod"))

## make a plot or two
traits_sub$n_rubisco_structure_mod = traits_sub$n_rubisco_mod + traits_sub$n_structure_mod
plot(traits_sub$narea ~ traits_sub$n_rubisco_structure_mod) # 
plot(traits_sub$narea ~ traits_sub$n_rubisco_mod) # not a great relationship
plot(traits_sub$narea ~ traits_sub$n_structure_mod) # coming mostly from the LMA results

plot(traits_sub$n_rubisco_mod ~ traits_sub$tmp) # some odd responses
plot(traits_sub$n_rubisco_mod ~ traits_sub$par)
plot(traits_sub$n_rubisco_mod ~ traits_sub$vpd)
plot(traits_sub$n_rubisco_mod ~ traits_sub$leaf_C13_delta_PDB)


