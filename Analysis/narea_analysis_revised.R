# Script using least-cost to explore Narea response at nutnet
## generally follows the analysis structure of Dong et al. (2017) Biogeosciences
## but includes the addition of soil N as a predictor

############################
#### load packages ####
############################
library(tidyverse)
library(lme4)
library(car)
library(r2glmm)
library(treemapify)
library(emmeans)
library(relaimpo)
library(patchwork)
library(multcomp)
library(ggplot2)

############################
#### load functions ####
############################
### functions to calculate vcmax and jmax 
source('optimal_vcmax_R/calc_optimal_vcmax.R')
sourceDirectory('optimal_vcmax_R/functions', modifiedOnly = FALSE)

### function to calculate relative importance for mixed models 
### from https://gist.github.com/BERENZ/e9b581a4b7160357934e
calc.relip.mm <- function(model,type = 'lmg') {
  if (!isLMM(model) & !isGLMM(model)) {
    stop('Currently supports only lmer/glmer objects', call. = FALSE)
  }
  require(lme4)
  X <- getME(model,'X')
  X <- X[ , -1]
  Y <- getME(model, 'y')
  s_resid <- sigma(model)
  s_effect <- getME(model, 'theta') * s_resid
  s2 <- sum(s_resid^2, s_effect^2)
  V <- Diagonal(x = s2, n = nrow(X))
  YX <- cbind(Y, X)
  cov_XY <- solve(t(YX) %*% solve(V) %*% as.matrix(YX))
  colnames(cov_XY) <- rownames(cov_XY) <- colnames(YX)
  importances <- calc.relimp(as.matrix(cov_XY), rela = F, type = type)
  return(importances)
}

calc.relip.boot.mm <- function(model,type = 'lmg') {
  if (!isLMM(model) & !isGLMM(model)) {
    stop('Currently supports only lmer/glmer objects', call. = FALSE)
  }
  require(lme4)
  X <- getME(model,'X')
  X <- X[ , -1]
  Y <- getME(model, 'y')
  s_resid <- sigma(model)
  s_effect <- getME(model, 'theta') * s_resid
  s2 <- sum(s_resid^2, s_effect^2)
  V <- Diagonal(x = s2, n = nrow(X))
  YX <- cbind(Y, X)
  cov_XY <- solve(t(YX) %*% solve(V) %*% as.matrix(YX))
  colnames(cov_XY) <- rownames(cov_XY) <- colnames(YX)
  bootresults <- boot.relimp(as.matrix(cov_XY), b=1000, rela = F, type = type)
  importances <- booteval.relimp(bootresults, norank=T)
  return(importances)
}

############################
#### load and manipulate data ####
############################

### leaf data
leaf_raw <- read.csv('../Data/processed/traits_v2.csv')
leaf = leaf_raw

leaf$Ntrt_fac <- as.factor(leaf$Ntrt)
leaf$Ptrt_fac <- as.factor(leaf$Ptrt)
leaf$Ktrt_fac <- as.factor(leaf$Ktrt)
leaf$block_fac <- as.factor(leaf$block)
leaf$trt_fac <- as.factor(leaf$trt)

# leaf$tmp_scaled <- leaf$tmp - 25
# 
# leaf$pft <- leaf$functional_group
# leaf$pft[leaf$pft == 'NULL'] <- NA
# leaf$grass <- 'no'
# leaf$grass[leaf$pft == 'GRASS'] <- 'yes'
leaf$fence <- 'no'
leaf$fence[leaf$trt == 'Fence' | leaf$trt == 'NPK+Fence'] <- 'yes'

leaf$logpar <- log(leaf$par2_gs)
leaf$logpar_per_leaf_area <- log(leaf$par_per_leaf_area)
leaf$logvpd <- log(leaf$vpd2_gs)
leaf$loglma <- log(leaf$lma)

leaf$nstar <- calc_nstar(leaf$tmp2_gs, leaf$z)
leaf$gammastar <- calc_gammastar_pa(leaf$tmp2_gs, leaf$z)
leaf$km <- calc_km_pa(leaf$tmp2_gs, leaf$z)
leaf$vpd_adj <- calc_vpd(leaf$tmp2_gs, leaf$z, leaf$vpd2_gs)
leaf$patm <- calc_patm(leaf$z)
leaf$ca <- 400 * 1e-6 * leaf$patm
Kp_25 <- 16 # Pa
Ea_Kp <- 36300 # J mol^-1, for the Pa parameter ## Boyd et al 2015
leaf$tmpK <- leaf$tmp2_gs + 273.15
leaf$kp <- Kp_25 * exp((Ea_Kp * (leaf$tmpK - 298.15))/(298.15 * 8.3145 * leaf$tmpK))

leaf$beta_num <- 1.6 * leaf$nstar * leaf$vpd_adj * 1000 * ((leaf$chi) ^ 2)
leaf$beta_denom <- ((1 - leaf$chi)^2) * (leaf$km + leaf$gammastar)
leaf$beta_denom[leaf$photosynthetic_pathway == 'C4'] <- ((1 - leaf$chi[leaf$photosynthetic_pathway == 'C4'])^2) * 
  (leaf$kp[leaf$photosynthetic_pathway == 'C4'] + leaf$gammastar[leaf$photosynthetic_pathway == 'C4'])

leaf$beta <- leaf$beta_num / leaf$beta_denom

## subset outliers
leaf_sub <- subset(leaf, chi < 0.95 & beta < 2000 & chi > 0.45)

# leaf$xi <- sqrt(leaf$beta * (leaf$km + leaf$gammastar) * (1/(1.6 * leaf$nstar)))
# leaf$xi[leaf$photosynthetic_pathway == 'C4'] <- sqrt(leaf$beta[leaf$photosynthetic_pathway == 'C4'] * 
#                                                   (leaf$kp[leaf$photosynthetic_pathway == 'C4'] + 
#                                                      leaf$gammastar[leaf$photosynthetic_pathway == 'C4']) * 
#                                                   (1/(1.6 * leaf$nstar[leaf$photosynthetic_pathway == 'C4'])))
# leaf$chi_check <- (leaf$xi / (leaf$xi + sqrt(leaf$vpd_adj * 1000)))

# plot(leaf$chi_check ~ leaf$chi)
# abline(0, 1)
# 
# hist(leaf$beta_num)
# hist(leaf$beta_denom)
# hist((1 - leaf$chi)^2)
# plot(leaf$beta_num ~ leaf$chi)
# plot(leaf$beta_denom ~ leaf$chi)
# plot(leaf$beta ~ leaf$chi)
# hist(leaf$beta)

leaf_n <- subset(leaf_sub, trt == 'N' | trt == 'Control')
leaf_n$trt_fac

### calculate treatment type averages
leaf_site_N <- leaf_n %>%
  group_by(site_code, Ntrt_fac,
           block_fac, 
           trt_fac, trt,
           Taxon 
           # grass, Nfix, photosynthetic_pathway
  ) %>%
  summarise_at(vars(narea, spp_lai, chi, spp_live_mass, spp_mass_N, max_cover, lma, p_pet, 
                    vpd2_gs, tmp2_gs, alpha,
                    Ambient_PAR, Ground_PAR, par_rat, beta),
               mean, na.rm = TRUE)

### subset by low N treatment
leaf_site_lowN <- subset(leaf_site_N, Ntrt_fac == '0')
nrow(leaf_site_lowN)

### subset by high N treatment
leaf_site_highN <- subset(leaf_site_N, Ntrt_fac == '1')
nrow(leaf_site_highN)

### combine low N and high N subsets
leaf_site_trt <- left_join(leaf_site_lowN, leaf_site_highN, 
                           by = c('site_code', 
                                  'block_fac',
                                  'Taxon' 
                                  # 'grass', 'Nfix', 'photosynthetic_pathway'
                           ))

### calculate percent change from the high N plots to the low N plots
leaf_site_trt$delta_narea <- ((leaf_site_trt$narea.y - 
                                 leaf_site_trt$narea.x) / leaf_site_trt$narea.x) * 100
leaf_site_trt$delta_live_mass <- ((leaf_site_trt$spp_live_mass.y - 
                                     leaf_site_trt$spp_live_mass.x) / leaf_site_trt$spp_live_mass.x) * 100
leaf_site_trt$delta_N_mass <- ((leaf_site_trt$spp_mass_N.y - 
                                  leaf_site_trt$spp_mass_N.x) / leaf_site_trt$spp_mass_N.x) * 100
leaf_site_trt$delta_chi <- ((leaf_site_trt$chi.y - 
                               leaf_site_trt$chi.x) / leaf_site_trt$chi.x) * 100
leaf_site_trt$delta_lma <- ((leaf_site_trt$lma.y - 
                               leaf_site_trt$lma.x) / leaf_site_trt$lma.x) * 100
leaf_site_trt$delta_beta <- ((leaf_site_trt$beta.y - 
                               leaf_site_trt$beta.x) / leaf_site_trt$beta.x) * 100
delta_narea_mad <- mad(leaf_site_trt$delta_narea, na.rm = T)
delta_lma_mad <- mad(leaf_site_trt$delta_lma, na.rm = T)
delta_live_mass_mad <- mad(leaf_site_trt$delta_live_mass, na.rm = T)
delta_N_mass_mad <- mad(leaf_site_trt$delta_N_mass, na.rm = T)
delta_chi_mad <- mad(leaf_site_trt$delta_chi, na.rm = T)
delta_beta_mad <- mad(leaf_site_trt$delta_beta, na.rm = T)

delta_beta_data <- subset(leaf_site_trt,
                               delta_beta < 2 * delta_beta_mad &
                                 delta_beta > 2 * -delta_beta_mad &
                            # delta_lma < 1 * delta_lma_mad &
                            # delta_lma > 1 * -delta_lma_mad &
                                 delta_live_mass < 2 * delta_live_mass_mad &
                                 delta_live_mass > 2 * -delta_live_mass_mad)

############################
#### run analyses ####
############################

#### Hyp 1: soil N reduces beta, but this is reduced when water availability is low
beta_lmer <- lmer(log(beta) ~ trt_fac + alpha +
                         Nfix + # photosynthetic_pathway +
                         (1|Taxon) + (1|Taxon:site_code) + (1|Taxon:site_code:block_fac), 
                       data = leaf_n)
plot(resid(beta_lmer) ~ fitted(beta_lmer))
summary(beta_lmer) # N = 
Anova(beta_lmer, type = 'III')
emmeans(beta_lmer, ~trt_fac)
emtrends(beta_lmer, ~1, var = 'alpha')
# emmeans(beta_lmer, ~trt_fac, at = list(alpha = 1.24))
# emmeans(beta_lmer, ~trt_fac, at = list(alpha = 0.35))

#### Hyp 1 supported, soil N reduces beta, but greater effect with lower water availability

#### Hyp 1b: the beta response is less negative when AGB increases
hist(delta_beta_data$delta_beta)
hist(delta_beta_data$delta_live_mass)
hist(delta_beta_data$delta_lma)
delta_beta_lmer <- lmer(delta_beta ~ delta_live_mass * alpha.x +
                          (1|Taxon) + (1|Taxon:site_code), data = delta_beta_data)
plot(resid(delta_beta_lmer) ~ fitted(delta_beta_lmer))
summary(delta_beta_lmer)
Anova(delta_beta_lmer, type = 'III')
emtrends(delta_beta_lmer, ~1, var = 'delta_live_mass')
emtrends(delta_beta_lmer, ~1, var = 'alpha.x')

plot(delta_beta_data$delta_beta ~ delta_beta_data$delta_live_mass)

#### Hyp 1b: no support that the AGB response matters

#### Hyp 2: chi is negatively related to soil N and vpd, positively related to soil water and T, and lower in C4
chi_lmer <- lmer(log(chi) ~ trt_fac + alpha +
                   photosynthetic_pathway +
                   tmp2_gs + vpd2_gs +
                   (1|Taxon) + (1|Taxon:site_code) + (1|Taxon:site_code:block_fac), 
                 data = leaf_n)
plot(resid(chi_lmer) ~ fitted(chi_lmer))
summary(chi_lmer) # N = 
Anova(chi_lmer, type = 'III')
emmeans(chi_lmer, ~trt_fac) # negative effect
emmeans(chi_lmer, ~photosynthetic_pathway) # lower in C4
emtrends(chi_lmer, ~1, var = 'tmp2_gs') # positive

#### Hyp 2: soil N (+), C4 (-), and T (+) confirmed...no impact of alpha or vpd
