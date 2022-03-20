# Script using least-cost to explore Narea response at nutnet
## generally follows the analysis structure of Dong et al. (2017) Biogeosciences
## but includes the addition of soil N as a predictor

############################
#### load packages ####
############################
library(tidyverse)
library(lme4)
library(car)
# library(r2glmm)
# library(treemapify)
library(marginaleffects)
library(emmeans)
# library(relaimpo)
# library(patchwork)
# library(multcomp)
# library(ggplot2)
library(piecewiseSEM)
library(lavaan)

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
leaf_raw <- read.csv('../Data/processed/traits_v3.csv')
leaf = leaf_raw

leaf$Ntrt_fac <- as.factor(leaf$Ntrt)
leaf$Ptrt_fac <- as.factor(leaf$Ptrt)
leaf$Ktrt_fac <- as.factor(leaf$Ktrt)
leaf$block_fac <- as.factor(leaf$block)
leaf$trt_fac <- as.factor(leaf$trt)

### add growing season length data
gs_information <- read.csv('growing_season_data/Weather_site_inventory_20200921.csv')
head(gs_information)
colnames(gs_information)

gs_length <- gs_information[,c(2, 20)]

leaf <- left_join(leaf, gs_length)
leaf$site_code[is.na(leaf$gs_len)]
leaf$gs_len[leaf$site_code == 'gilb.za'] <- 6
leaf$gs_len[leaf$site_code == 'summ.za'] <- 8

leaf$gs_frac <- leaf$gs_len / 12


## remove c4 (not enough, only 18 out of 313 in final)
# leaf <- subset(leaf, photosynthetic_pathway == 'C3')

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
# leaf$beta_denom[leaf$photosynthetic_pathway == 'C4'] <- ((1 - leaf$chi[leaf$photosynthetic_pathway == 'C4'])^2) * 
#   (leaf$kp[leaf$photosynthetic_pathway == 'C4'] + leaf$gammastar[leaf$photosynthetic_pathway == 'C4'])

leaf$beta <- leaf$beta_num / leaf$beta_denom

leaf$nmass <- leaf$leaf_pct_N

## subset outliers using MAD
nrow(leaf)
beta_median <- median(leaf$beta, na.rm = T)
chi_median <- median(leaf$chi, na.rm = T)
nmass_median <- median(leaf$nmass, na.rm = T)
lma_median <- median(leaf$lma, na.rm = T)
narea_median <- median(leaf$narea, na.rm = T)

beta_mad <- mad(leaf$beta, na.rm = T)
chi_mad <- mad(leaf$chi, na.rm = T)
nmass_mad <- mad(leaf$nmass, na.rm = T)
lma_mad <- mad(leaf$lma, na.rm = T)
narea_mad <- mad(leaf$narea, na.rm = T)

#### remove values more than 3 MAD (VERY CONSERVATIVE: https://www.sciencedirect.com/science/article/pii/S0022103113000668)
# leaf_sub <- subset(leaf, beta < beta_median + 6 * beta_mad &
#                      beta > beta_median - 6 * beta_mad &
#                      nmass < nmass_median + 6 * nmass_mad &
#                      nmass > nmass_median - 6 * nmass_mad &
#                      lma < lma_median + 6 * lma_mad &
#                      lma > lma_median - 6 * lma_mad &
#                      chi < chi_median + 6 * chi_mad &
#                      chi > chi_median - 6 * chi_mad &
#                      narea < narea_median + 6 * narea_mad &
#                      narea > narea_median - 6 * narea_mad)



# leaf_sub <- subset(leaf, chi < 0.95 & beta < 2000 & beta >1 & chi > 0.05)
leaf_sub <- subset(leaf, chi < 0.95 & chi > 0.05)
# leaf_sub <- leaf_sub[!(leaf_sub$photosynthetic_pathway == "C3" & leaf_sub$chi < 0.2),]

nrow(leaf_sub)

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

leaf_n <- leaf_sub
# leaf_n <- subset(leaf_sub, trt == 'N' | trt == 'Control')
# leaf_n$trt_fac

nrow(leaf_n)
levels(as.factor(leaf_n$site_code))
nrow(subset(leaf_n, photosynthetic_pathway == 'C4'))
nrow(subset(leaf_n, Nfix == 'yes'))

### calculate treatment type averages
leaf_site_N <- leaf_n %>%
  group_by(site_code, Ntrt_fac, Ptrt_fac, Ktrt_fac,
           block_fac, 
           trt_fac, trt
           # Taxon 
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
                                  'block_fac'
                                  # 'Taxon' 
                                  # 'grass', 'Nfix', 'photosynthetic_pathway'
                           ))

### calculate percent change from the high N plots to the low N plots
# leaf_site_trt$delta_narea <- ((leaf_site_trt$narea.y - 
#                                  leaf_site_trt$narea.x) / leaf_site_trt$narea.x) * 100
leaf_site_trt$delta_live_mass <- ((leaf_site_trt$spp_live_mass.y - 
                                    leaf_site_trt$spp_live_mass.x) / leaf_site_trt$spp_live_mass.x) * 100
# leaf_site_trt$delta_N_mass <- ((leaf_site_trt$spp_mass_N.y - 
#                                   leaf_site_trt$spp_mass_N.x) / leaf_site_trt$spp_mass_N.x) * 100
# leaf_site_trt$delta_chi <- ((leaf_site_trt$chi.y - 
#                                leaf_site_trt$chi.x) / leaf_site_trt$chi.x) * 100
# leaf_site_trt$delta_lma <- ((leaf_site_trt$lma.y - 
#                                leaf_site_trt$lma.x) / leaf_site_trt$lma.x) * 100
leaf_site_trt$delta_beta <- ((leaf_site_trt$beta.y - 
                              leaf_site_trt$beta.x) / leaf_site_trt$beta.x) * 100
# delta_narea_mad <- mad(leaf_site_trt$delta_narea, na.rm = T)
# delta_lma_mad <- mad(leaf_site_trt$delta_lma, na.rm = T)
# delta_N_mass_mad <- mad(leaf_site_trt$delta_N_mass, na.rm = T)
# delta_chi_mad <- mad(leaf_site_trt$delta_chi, na.rm = T)
delta_live_mass_median <- median(leaf_site_trt$delta_live_mass, na.rm = T)
delta_beta_median <- median(leaf_site_trt$delta_beta, na.rm = T)
delta_live_mass_mad <- mad(leaf_site_trt$delta_live_mass, na.rm = T)
delta_beta_mad <- mad(leaf_site_trt$delta_beta, na.rm = T)

delta_beta_data <- subset(leaf_site_trt,
                              delta_beta < delta_beta_median + 3 * delta_beta_mad &
                              delta_beta > delta_beta_median - 3 * delta_beta_mad &
                              delta_live_mass < delta_live_mass_median + 3 * delta_live_mass_mad &
                              delta_live_mass > delta_live_mass_median - 3 * delta_live_mass_mad)

# delta_beta_data <- leaf_site_trt

############################
#### run analyses ####
############################

#### Hyp 1: soil N reduces beta, water availability increases beta, Nfix increases beta
beta_lmer <- lmer(log(beta) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac + alpha +
                         Nfix + photosynthetic_pathway +
                         (1|Taxon) + (1|site_code) + (1|site_code:block_fac), 
                       data = leaf_n)
plot(resid(beta_lmer) ~ fitted(beta_lmer))
summary(beta_lmer) # N = 2118
rsquared(beta_lmer)
Anova(beta_lmer, type = 'III')
beta_mfx <- marginaleffects(beta_lmer)
summary(beta_mfx)
beta_lmer_Ntrt_emmeans <- summary(emmeans(beta_lmer, ~Ntrt_fac))
# emmeans(beta_lmer, ~Nfix)
beta_lmer_c4_emmeans <- summary(emmeans(beta_lmer, ~photosynthetic_pathway))
beta_change <- ((exp(beta_lmer_Ntrt_emmeans[2, 2]) - exp(beta_lmer_Ntrt_emmeans[1, 2]))/
  exp(beta_lmer_Ntrt_emmeans[1, 2])) * 100
beta_change_c4 <- ((exp(beta_lmer_c4_emmeans[2, 2]) - exp(beta_lmer_c4_emmeans[1, 2]))/
                  exp(beta_lmer_c4_emmeans[1, 2])) * 100

#### H1 support: Soil N reduces beta, C4 reduces beta

#### Hyp 1b: the beta response is less negative when AGB increases
beta_agb_lmer <- lmer(log(beta) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac * spp_live_mass + 
                        spp_live_mass * alpha +
                    Nfix + photosynthetic_pathway +
                    (1|Taxon) + (1|site_code) + (1|site_code:block_fac), 
                  data = leaf_n)
plot(resid(beta_agb_lmer) ~ fitted(beta_agb_lmer))
summary(beta_agb_lmer) # N = 1763
rsquared(beta_agb_lmer)
Anova(beta_agb_lmer, type = 'III')
summary(emmeans(beta_agb_lmer, ~Ntrt_fac))
summary(emmeans(beta_agb_lmer, ~Ptrt_fac))
summary(emmeans(beta_agb_lmer, ~photosynthetic_pathway))

hist(delta_beta_data$delta_beta)
hist(log(delta_beta_data$delta_beta + 100))
hist(delta_beta_data$delta_live_mass)
hist(log(delta_beta_data$delta_live_mass + 100))
# hist(delta_beta_data$delta_lma)
# hist(delta_beta_data$delta_N_mass)
delta_beta_lmer <- lmer(delta_beta ~ delta_live_mass +
                          # (1|Taxon) + 
                          (1|site_code), data = delta_beta_data)
                          # (1|site_code), data = delta_beta_data)
plot(resid(delta_beta_lmer) ~ fitted(delta_beta_lmer))
summary(delta_beta_lmer)
rsquared(delta_beta_lmer)
Anova(delta_beta_lmer, type = 'III')
delta_beta_mfx <- marginaleffects(delta_beta_lmer)
summary(delta_beta_mfx)
test(emtrends(delta_beta_lmer, ~1, var = 'delta_live_mass'))

plot(delta_beta_data$delta_beta ~ delta_beta_data$delta_live_mass)
# plot(log(delta_beta_data$delta_beta + 100) ~ log(delta_beta_data$delta_live_mass + 100))

#### Hyp 1b: support that the AGB response matters (lowers the beta response)

#### Hyp 2: chi is negatively related to soil N and vpd, positively related to soil water and T, and lower in C4
##### note no hypothesized effcect of Nfix as this didn't come out in the beta impacts
chi_lmer <- lmer(logit(chi) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac + # alpha +
                   Nfix + photosynthetic_pathway +
                   tmp2_gs + vpd2_gs +
                   (1|Taxon) + (1|site_code) + (1|site_code:block_fac), 
                 data = leaf_n)
plot(resid(chi_lmer) ~ fitted(chi_lmer))
summary(chi_lmer) # N = 
rsquared(chi_lmer)
Anova(chi_lmer, type = 'III')
chi_mfx <- marginaleffects(chi_lmer)
summary(chi_mfx)
emmeans(chi_lmer, ~Ntrt_fac) # negative effect
emmeans(chi_lmer, ~photosynthetic_pathway) # lower in C4
emtrends(chi_lmer, ~1, var = 'alpha') # positive

#### Hyp 2: soil N (-), alpha (+), and C4 (-) confirmed...no impact of Nfix, T, or vpd

#### Hyp 3: Nmass positively related to light, vpd, soil N, C4, and Nfix
#### negatively related to temperature and soil water
# nmass_lmer <- lmer(log(nmass) ~ trt_fac + alpha + 
#                      Nfix + #photosynthetic_pathway +
#                    tmp2_gs + vpd2_gs +
#                      par2_gs +
#                    (1|Taxon) + (1|site_code) + (1|site_code:block_fac), 
#                  data = leaf_n)
# plot(resid(nmass_lmer) ~ fitted(nmass_lmer))
# summary(nmass_lmer) # N = 
# Anova(nmass_lmer, type = 'III')
# nmass_mfx <- marginaleffects(nmass_lmer)
# summary(nmass_mfx)

#### Hyp 3: soil N (+), light (+), Nfix (+) confirmed
#### no temperature vpd or alpha effect

#### Hyp 4: Marea positively related to light
#### negatively related to temperature
# marea_lmer <- lmer(log(lma) ~ tmp2_gs + par2_gs +
#                      (1|Taxon) + (1|site_code) + (1|site_code:block_fac), 
#                    data = leaf_n)
# plot(resid(marea_lmer) ~ fitted(marea_lmer))
# summary(marea_lmer) # N = 
# Anova(marea_lmer, type = 'III')
# marea_mfx <- marginaleffects(marea_lmer)
# summary(marea_mfx)

#### Hyp 4: no effects, so we should expect same impact on Nmass as to Narea

#### Hyp 5: Narea positively related to light, vpd, soil N, C4, and Nfix
#### negatively related to temperature and soil water
# narea_lmer <- lmer(log(narea) ~ trt_fac + alpha + 
#                      Nfix + #photosynthetic_pathway +
#                      tmp2_gs + vpd2_gs +
#                      par2_gs +
#                      (1|Taxon) + (1|site_code) + (1|site_code:block_fac), 
#                    data = leaf_n)
# plot(resid(narea_lmer) ~ fitted(narea_lmer))
# summary(narea_lmer) # N = 
# Anova(narea_lmer, type = 'III')
# narea_mfx <- marginaleffects(narea_lmer)
# summary(narea_mfx)

#### Hyp 5: only a positive soil N effect, but others washed out

#### Hyp 6: we can estimate leaf chi from optimality
beta_values <- summary(emmeans(beta_lmer, ~ Ntrt_fac * photosynthetic_pathway))
beta_lowN_c3 <- exp(beta_values[1, 3])
beta_highN_c3 <- exp(beta_values[2, 3])
beta_lowN_c4 <- exp(beta_values[3, 3])
beta_highN_c4 <- exp(beta_values[4, 3])
leaf_n$beta_model[leaf_n$Ntrt_fac == '0' & leaf_n$photosynthetic_pathway == 'C3'] <- beta_lowN_c3
leaf_n$beta_model[leaf_n$Ntrt_fac == '0' & leaf_n$photosynthetic_pathway == 'C4'] <- beta_lowN_c4
leaf_n$beta_model[leaf_n$Ntrt_fac == '1' & leaf_n$photosynthetic_pathway == 'C3'] <- beta_highN_c3
leaf_n$beta_model[leaf_n$Ntrt_fac == '1' & leaf_n$photosynthetic_pathway == 'C4'] <- beta_highN_c4
leaf_n$beta_average <- mean(leaf_n$beta, na.rm = T)

pred_c3 <- calc_optimal_vcmax(pathway = "C3",
                           tg_c = leaf_n$tmp2_gs, 
                           paro = leaf_n$par2_gs, 
                           cao = 400, 
                           vpdo = leaf_n$vpd2_gs, 
                           z = leaf_n$z,
                           q0_resp = "no",
                           f = leaf_n$gs_frac,
                           beta = leaf_n$beta_average)

pred_c4 <- calc_optimal_vcmax(pathway = "C4",
                              tg_c = leaf_n$tmp2_gs, 
                              paro = leaf_n$par2_gs, 
                              cao = 400, 
                              vpdo = leaf_n$vpd2_gs, 
                              z = leaf_n$z,
                              q0_resp = "no",
                              f = leaf_n$gs_frac,
                              beta = leaf_n$beta_average)

## add model results to leaf dataset
pred_c3$photosynthetic_pathway <- leaf_n$photosynthetic_pathway
pred_c3$chi_pred <- pred_c3$chi
pred_c3$chi_pred[pred_c3$photosynthetic_pathway == 'C4'] <- 0
pred_c4$photosynthetic_pathway <- leaf_n$photosynthetic_pathway
pred_c4$chi_pred <- pred_c4$chi
pred_c4$chi_pred[pred_c3$photosynthetic_pathway == 'C3'] <- 0

leaf_n$chi_pred <- pred_c3$chi_pred + pred_c4$chi_pred

## predicting chi from optimality
chi_pred_lm <- lm(chi ~ chi_pred + Ntrt_fac, data = leaf_n)
plot(resid(chi_pred_lm) ~ fitted(chi_pred_lm))
summary(chi_pred_lm) # N = 
Anova(chi_pred_lm, type = 'III')
chi_mfx <- marginaleffects(chi_pred_lm)
summary(chi_mfx)

plot(chi ~ chi_pred, data = leaf_n)

## now do the same with the beta adjustments
pred_c3_beta <- calc_optimal_vcmax(pathway = "C3",
                              tg_c = leaf_n$tmp2_gs, 
                              paro = leaf_n$par2_gs, 
                              cao = 400, 
                              vpdo = leaf_n$vpd2_gs, 
                              z = leaf_n$z,
                              q0_resp = "no",
                              f = leaf_n$gs_frac,
                              beta = leaf_n$beta_model)

pred_c4_beta <- calc_optimal_vcmax(pathway = "C4",
                              tg_c = leaf_n$tmp2_gs, 
                              paro = leaf_n$par2_gs, 
                              cao = 400, 
                              vpdo = leaf_n$vpd2_gs, 
                              z = leaf_n$z,
                              q0_resp = "no",
                              f = leaf_n$gs_frac,
                              beta = leaf_n$beta_model)

## add model results to leaf dataset
pred_c3_beta$photosynthetic_pathway <- leaf_n$photosynthetic_pathway
pred_c3_beta$chi_pred_beta <- pred_c3_beta$chi
pred_c3_beta$chi_pred_beta[pred_c3_beta$photosynthetic_pathway == 'C4'] <- 0
pred_c4_beta$photosynthetic_pathway <- leaf_n$photosynthetic_pathway
pred_c4_beta$chi_pred_beta <- pred_c4_beta$chi
pred_c4_beta$chi_pred_beta[pred_c3_beta$photosynthetic_pathway == 'C3'] <- 0

leaf_n$chi_pred_beta <- pred_c3_beta$chi_pred_beta + pred_c4_beta$chi_pred_beta
hist(leaf_n$chi_pred_beta)
hist(leaf_n$chi_pred)

## make dataset summarized by site and plant type
leaf_n_group_by <- group_by(leaf_n, site_code, photosynthetic_pathway)
leaf_n_site_c4 <- summarise(leaf_n_group_by, chi_mean = mean(chi, na.rm = T),
                            chi_pred_beta_mean = mean(chi_pred_beta, na.rm = T))

## predicting chi from optimality
chi_pred_beta_lm <- lm(chi ~ chi_pred_beta, data = leaf_n)
plot(resid(chi_pred_beta_lm) ~ fitted(chi_pred_beta_lm))
summary(chi_pred_beta_lm) # N = 
Anova(chi_pred_beta_lm, type = 'III')
chi_beta_mfx <- marginaleffects(chi_pred_beta_lm)
summary(chi_beta_mfx)
emtrends(chi_pred_beta_lm, ~1, var = 'chi_pred_beta_mean')

plot(chi ~ chi_pred_beta, data = leaf_n)
abline(0,1)

## prediction plot
chi_beta_pred_plot <- ggplot(data = leaf_n, aes(x = chi_pred_beta, y = chi)) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(2), colour = 'black'),
        axis.title.x = element_text(size = rel(2), colour = 'black'),
        axis.text.x = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_abline(slope=1, intercept=0, lty = 2) +
  geom_point(alpha = 0.2) +
  geom_smooth(method = 'lm') + #, formula=y~0+x)
  xlab('Predicted χ') +
  ylab('Observed χ') +
  xlim(c(0.5, 0.9)) +
  ylim(c(0,1))

## SEM prediction of chi
leaf_n$Ntrt_cont[leaf_n$Ntrt_fac == '0'] <- 0
leaf_n$Ptrt_cont[leaf_n$Ptrt_fac == '0'] <- 0
leaf_n$Ntrt_cont[leaf_n$Ntrt_fac == '1'] <- 1
leaf_n$Ptrt_cont[leaf_n$Ptrt_fac == '1'] <- 1
leaf_n$Nfix_cont[leaf_n$Nfix == 'no'] <- 0
leaf_n$Nfix_cont[leaf_n$Nfix == 'yes'] <- 1
leaf_n$photosynthetic_pathway_cont[leaf_n$photosynthetic_pathway == 'C3'] <- 0
leaf_n$photosynthetic_pathway_cont[leaf_n$photosynthetic_pathway == 'C4'] <- 1

leaf_n$scale.beta <- scale(leaf_n$beta)
leaf_n$scale.Ntrt_cont <- scale(leaf_n$Ntrt_cont)
leaf_n$scale.Ptrt_cont <- scale(leaf_n$Ptrt_cont)
leaf_n$scale.Nfix_cont <- scale(leaf_n$Nfix_cont)
leaf_n$scale.photosynthetic_pathway_cont <- scale(leaf_n$photosynthetic_pathway_cont)
leaf_n$scale.alpha <- scale(leaf_n$alpha)
leaf_n$scale.chi <- scale(leaf_n$chi)
leaf_n$scale.vpd2_gs <- scale(leaf_n$vpd2_gs)
leaf_n$scale.tmp2_gs <- scale(leaf_n$tmp2_gs)

chi_path <- 'scale.beta ~ scale.Ntrt_cont + scale.Ptrt_cont + scale.Nfix_cont + scale.photosynthetic_pathway_cont + scale.alpha
             scale.chi ~ scale.beta + scale.vpd2_gs + scale.tmp2_gs'

fit.chi_path <- sem(chi_path, data = leaf_n)
summary(fit.chi_path, standardized = T, rsq = T)











## predicting nmass from optimality
leaf_n_pred$nmass.y <- (leaf_n_pred$nall / leaf_n_pred$lma.y)*100
nmass_pred_lm <- lm(nmass ~ nmass.y + trt_fac, data = leaf_n_pred)
plot(resid(nmass_pred_lm) ~ fitted(nmass_pred_lm))
summary(nmass_pred_lm) # N = 
Anova(nmass_pred_lm, type = 'III')
nmass_mfx <- marginaleffects(nmass_pred_lm)
summary(nmass_mfx)

## predicting narea from optimality
narea_pred_lm <- lm(narea ~ nall + trt_fac, data = leaf_n_pred)
plot(resid(narea_pred_lm) ~ fitted(narea_pred_lm))
summary(narea_pred_lm) # N = 
Anova(narea_pred_lm, type = 'III')
narea_mfx <- marginaleffects(narea_pred_lm)
summary(narea_mfx)



