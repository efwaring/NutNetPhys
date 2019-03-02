# Script to compare Narea observed across the Nutrient Network to that predicted
# from climate and LMA

## load libraries
library(R.utils)

## load optimal vcmax script (Smith et al., Ecology Letters)
source('optimal_vcmax_R/calc_optimal_vcmax.R')

## load observational data
traits = read.csv('../Data/leaf_plus.csv')
head(traits)

## calculate lma (in g m-2) and narea
traits$lma = 1 / (traits$SLA * (1 / 1000000)) # conver mm2 g-1 to g m-2
traits$narea = (traits$leaf_pct_N / 100) * (traits$lma)
hist(traits$lma) # some extremely high values
hist(traits$narea) # some extremely high values

## calculate optimal trait values
traits$vcmax_mod = calc_optimal_vcmax(tg_c = traits$tmp, paro = traits$par, cao = 400, 
                                      vpdo = traits$vpd, z = traits$z)$vcmax_prime
traits$vcmax25_mod = traits$vcmax_mod / calc_tresp_mult(traits$tmp, traits$tmp,
                                                        25)

## calculate N from rubisco and structure
### specific activity of rubisco is 60 umol CO2 gRubisco-1 s-1
### mass ratio of rubisco is 7.16 gRubisco gNRubisco-1
### lma conversion from Dong et al. (2017) Biogeosciences
traits$n_rubisco_mod = traits$vcmax25_mod * (1/60) * (1/7.16)
traits$n_structure_mod = (10^-2.67) * (traits$lma ^ 0.99)
hist(traits$n_rubisco_mod)
hist(traits$n_structure_mod)

## compare modeled and observed N
### first remove very high lma (>1000 g m-2)
traits_sub = subset(traits, lma < 1000)
n_pred_lm = lm(narea ~ n_rubisco_mod + n_structure_mod, data = traits_sub)
anova(n_pred_lm)
summary(n_pred_lm)

## make a plot or two
traits_sub$n_rubisco_structure_mod = traits_sub$n_rubisco_mod + traits_sub$n_structure_mod
plot(traits_sub$narea ~ traits_sub$n_rubisco_structure_mod) # 
plot(traits_sub$narea ~ traits_sub$n_rubisco_mod) # not a great relationship
plot(traits_sub$narea ~ traits_sub$n_structure_mod) # coming mostly from the LMA results




