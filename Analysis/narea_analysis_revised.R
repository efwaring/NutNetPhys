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
library(multcomp)
library(ggplot2)
library(piecewiseSEM)
# library(lavaan)
# library(boot)

############################
#### load functions ####
############################
### functions to calculate vcmax and jmax 
source('../../optimal_vcmax_R/calc_optimal_vcmax.R')
sourceDirectory('../../optimal_vcmax_R/functions', modifiedOnly = FALSE)

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

multiplot <- function(..., plotlist=NULL, cols) {
  require(grid)
  
  # Make a list from the ... arguments and plotlist
  plots <- c(list(...), plotlist)
  
  numPlots = length(plots)
  
  # Make the panel
  plotCols = cols                          # Number of columns of plots
  plotRows = ceiling(numPlots/plotCols) # Number of rows needed, calculated from # of cols
  
  # Set up the page
  grid.newpage()
  pushViewport(viewport(layout = grid.layout(plotRows, plotCols)))
  vplayout <- function(x, y)
    viewport(layout.pos.row = x, layout.pos.col = y)
  
  # Make each plot, in the correct location
  for (i in 1:numPlots) {
    curRow = ceiling(i/plotCols)
    curCol = (i-1) %% plotCols + 1
    print(plots[[i]], vp = vplayout(curRow, curCol ))
  }
  
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

## add soil moisture data
soil <- read.csv('soil_data.csv')

leaf <- left_join(leaf, soil)

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

leaf$pft[leaf$photosynthetic_pathway == 'C4'] <- 'c4'
leaf$pft[leaf$photosynthetic_pathway == 'C3' & leaf$Nfix == 'no'] <- 'c3_noNfix'
leaf$pft[leaf$photosynthetic_pathway == 'C3' & leaf$Nfix == 'yes'] <- 'c3_yesNfix'

## subset outliers using MAD
nrow(leaf)
# beta_median <- median(leaf$beta, na.rm = T)
# chi_median <- median(leaf$chi, na.rm = T)
# nmass_median <- median(leaf$nmass, na.rm = T)
# lma_median <- median(leaf$lma, na.rm = T)
# narea_median <- median(leaf$narea, na.rm = T)
# 
# beta_mad <- mad(leaf$beta, na.rm = T)
# chi_mad <- mad(leaf$chi, na.rm = T)
# nmass_mad <- mad(leaf$nmass, na.rm = T)
# lma_mad <- mad(leaf$lma, na.rm = T)
# narea_mad <- mad(leaf$narea, na.rm = T)

#### remove values more than 2.5 MAD (VERY CONSERVATIVE: https://www.sciencedirect.com/science/article/pii/S0022103113000668)
# leaf_sub <- subset(leaf, beta < beta_median + 2.5 * beta_mad &
#                     beta > beta_median - 2.5 * beta_mad &
#                     chi < chi_median + 2.5 * chi_mad &
#                     chi > chi_median - 2.5 * chi_mad)



# leaf_sub <- subset(leaf, chi < 0.95 & beta < 2000 & beta >1 & chi > 0.05)
leaf_sub <- subset(leaf, chi < 0.95 & chi > 0.3 & pft != 'c4') # c4s are really throwing things off!
# leaf_sub <- leaf_sub[!(leaf_sub$photosynthetic_pathway == "C3" & leaf_sub$chi < 0.3),]
# leaf_sub <- subset(leaf_sub, beta < 1000)

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
# hist(leaf_sub$beta)

leaf_n <- leaf_sub
leaf_n$PKtrt[leaf_n$trt == 'Control' | leaf_n$trt == 'N' | leaf_n$trt == 'Fence'] <- 'Control'
leaf_n$PKtrt[leaf_n$trt == 'P' | leaf_n$trt == 'NP'] <- 'P'
leaf_n$PKtrt[leaf_n$trt == 'K' | leaf_n$trt == 'NK'] <- 'K'
leaf_n$PKtrt[leaf_n$trt == 'PK' | leaf_n$trt == 'NPK' | leaf_n$trt == 'NPK+Fence'] <- 'PK'
# leaf_n <- subset(leaf_sub, trt == 'N' | trt == 'Control')
leaf_n$trt_fac

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
# leaf_site_lowN <- subset(leaf_site_N, Ntrt_fac == '0')
# nrow(leaf_site_lowN)
# 
# ### subset by high N treatment
# leaf_site_highN <- subset(leaf_site_N, Ntrt_fac == '1')
# nrow(leaf_site_highN)

### combine low N and high N subsets
# leaf_site_trt <- left_join(leaf_site_lowN, leaf_site_highN, 
#                            by = c('site_code', 
#                                   'block_fac'
#                                   # 'Taxon' 
#                                   # 'grass', 'Nfix', 'photosynthetic_pathway'
#                            ))
# 
# ### calculate percent change from the high N plots to the low N plots
# # leaf_site_trt$delta_narea <- ((leaf_site_trt$narea.y - 
# #                                  leaf_site_trt$narea.x) / leaf_site_trt$narea.x) * 100
# leaf_site_trt$delta_live_mass <- ((leaf_site_trt$spp_live_mass.y - 
#                                     leaf_site_trt$spp_live_mass.x) / leaf_site_trt$spp_live_mass.x) * 100
# # leaf_site_trt$delta_N_mass <- ((leaf_site_trt$spp_mass_N.y - 
# #                                   leaf_site_trt$spp_mass_N.x) / leaf_site_trt$spp_mass_N.x) * 100
# # leaf_site_trt$delta_chi <- ((leaf_site_trt$chi.y - 
# #                                leaf_site_trt$chi.x) / leaf_site_trt$chi.x) * 100
# # leaf_site_trt$delta_lma <- ((leaf_site_trt$lma.y - 
# #                                leaf_site_trt$lma.x) / leaf_site_trt$lma.x) * 100
# leaf_site_trt$delta_beta <- ((leaf_site_trt$beta.y - 
#                               leaf_site_trt$beta.x) / leaf_site_trt$beta.x) * 100
# # delta_narea_mad <- mad(leaf_site_trt$delta_narea, na.rm = T)
# # delta_lma_mad <- mad(leaf_site_trt$delta_lma, na.rm = T)
# # delta_N_mass_mad <- mad(leaf_site_trt$delta_N_mass, na.rm = T)
# # delta_chi_mad <- mad(leaf_site_trt$delta_chi, na.rm = T)
# delta_live_mass_median <- median(leaf_site_trt$delta_live_mass, na.rm = T)
# delta_beta_median <- median(leaf_site_trt$delta_beta, na.rm = T)
# delta_live_mass_mad <- mad(leaf_site_trt$delta_live_mass, na.rm = T)
# delta_beta_mad <- mad(leaf_site_trt$delta_beta, na.rm = T)
# 
# delta_beta_data <- subset(leaf_site_trt,
#                               delta_beta < delta_beta_median + 3 * delta_beta_mad &
#                               delta_beta > delta_beta_median - 3 * delta_beta_mad &
#                               delta_live_mass < delta_live_mass_median + 3 * delta_live_mass_mad &
#                               delta_live_mass > delta_live_mass_median - 3 * delta_live_mass_mad)

# delta_beta_data <- leaf_site_trt

############################
#### run analyses ####
############################

#### Hyp 1: soil N reduces beta, water availability increases beta, Nfix increases beta
beta_lmer <- lmer(log(beta) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac *
                         pft +
                         p_pet +
                         (1|Taxon) + 
                         #(1|site_code) + # removed because correlated with p_pet 
                         (1|site_code:block_fac) + (1|site_code:Taxon), # +
                         # (1|PKtrt) + (1|site_code:PKtrt), 
                       data = leaf_n)
plot(resid(beta_lmer) ~ fitted(beta_lmer))
summary(beta_lmer) # N = XXX
rsquared(beta_lmer)
Anova(beta_lmer, type = 'II')
# beta_mfx <- marginaleffects(beta_lmer)
# summary(beta_mfx)
beta_lmer_Ntrt_emmeans <- summary(emmeans(beta_lmer, ~Ntrt_fac))
# emmeans(beta_lmer, ~Nfix)
beta_lmer_pft_emmeans <- summary(emmeans(beta_lmer, ~pft))
beta_lmer_Ntrt_pft_emmeans <- summary(emmeans(beta_lmer, ~Ntrt_fac * pft))
cld(emmeans(beta_lmer, ~Ntrt_fac * pft))
beta_change <- ((exp(beta_lmer_Ntrt_emmeans[2, 2]) - exp(beta_lmer_Ntrt_emmeans[1, 2]))/
  exp(beta_lmer_Ntrt_emmeans[1, 2])) * 100
# beta_change_N_c4 <- ((exp(beta_lmer_Ntrt_pft_emmeans[6, 3]) - exp(beta_lmer_Ntrt_pft_emmeans[5, 3]))/
#                        exp(beta_lmer_Ntrt_pft_emmeans[5, 3])) * 100
beta_change_N_c3nofix <- ((exp(beta_lmer_Ntrt_pft_emmeans[2, 3]) - exp(beta_lmer_Ntrt_pft_emmeans[1, 3]))/
                       exp(beta_lmer_Ntrt_pft_emmeans[2, 3])) * 100
beta_change_N_c3yesfix <- ((exp(beta_lmer_Ntrt_pft_emmeans[4, 3]) - exp(beta_lmer_Ntrt_pft_emmeans[3, 3]))/
                            exp(beta_lmer_Ntrt_pft_emmeans[3, 3])) * 100

#### H1 support: Soil N reduces beta but not in N fixers


  

#### Hyp 1b: the beta response is less negative when AGB increases
# beta_agb_lmer <- lmer(log(beta) ~ Ntrt_fac * spp_live_mass + 
#                        p_pet +
#                    Nfix + photosynthetic_pathway +
#                    (1|Taxon) + (1|site_code) + (1|site_code:block_fac) +
#                    (1|site_code:block_fac:PKtrt), 
#                  data = leaf_n)
# plot(resid(beta_agb_lmer) ~ fitted(beta_agb_lmer))
# summary(beta_agb_lmer) # N = 1656
# rsquared(beta_agb_lmer)
# Anova(beta_agb_lmer, type = 'II')
# summary(emmeans(beta_agb_lmer, ~Ntrt_fac))
# summary(emmeans(beta_agb_lmer, ~photosynthetic_pathway))

# hist(delta_beta_data$delta_beta)
# hist(log(delta_beta_data$delta_beta + 100))
# hist(delta_beta_data$delta_live_mass)
# hist(log(delta_beta_data$delta_live_mass + 100))
# # hist(delta_beta_data$delta_lma)
# # hist(delta_beta_data$delta_N_mass)
# delta_beta_lmer <- lmer(delta_beta ~ delta_live_mass +
#                           # (1|Taxon) + 
#                           (1|site_code), data = delta_beta_data)
#                           # (1|site_code), data = delta_beta_data)
# plot(resid(delta_beta_lmer) ~ fitted(delta_beta_lmer))
# summary(delta_beta_lmer)
# rsquared(delta_beta_lmer)
# Anova(delta_beta_lmer, type = 'III')
# delta_beta_mfx <- marginaleffects(delta_beta_lmer)
# summary(delta_beta_mfx)
# test(emtrends(delta_beta_lmer, ~1, var = 'delta_live_mass'))
# 
# plot(delta_beta_data$delta_beta ~ delta_beta_data$delta_live_mass)
# plot(log(delta_beta_data$delta_beta + 100) ~ log(delta_beta_data$delta_live_mass + 100))

#### Hyp 1b: support that the AGB response matters (lowers the beta response)

#### Hyp 2: chi is negatively related to soil N and vpd, positively related to soil water and T, and lower in C4
##### note no hypothesized effcect of Nfix as this didn't come out in the beta impacts
chi_lmer <- lmer(logit(chi) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac *
                   pft +
                   p_pet +
                   (1|Taxon) + 
                   #(1|site_code) + # removed because correlated with p_pet 
                   (1|site_code:block_fac) + (1|site_code:Taxon), # +
                 # (1|PKtrt) + (1|site_code:PKtrt), 
                 data = leaf_n)
plot(resid(chi_lmer) ~ fitted(chi_lmer))
summary(chi_lmer) # N = XXXX
rsquared(chi_lmer)
Anova(chi_lmer, type = 'II')
# chi_mfx <- marginaleffects(chi_lmer)
# summary(chi_mfx)
chi_lmer_Ntrt_emmeans <- summary(emmeans(chi_lmer, ~Ntrt_fac))
# emmeans(chi_lmer, ~Nfix)
chi_lmer_pft_emmeans <- summary(emmeans(chi_lmer, ~pft))
chi_lmer_Ntrt_pft_emmeans <- summary(emmeans(chi_lmer, ~Ntrt_fac * pft))
cld(emmeans(chi_lmer, ~Ntrt_fac * pft))
chi_change <- ((exp(chi_lmer_Ntrt_emmeans[2, 2]) - exp(chi_lmer_Ntrt_emmeans[1, 2]))/
                  exp(chi_lmer_Ntrt_emmeans[1, 2])) * 100
# chi_change_N_c4 <- ((exp(chi_lmer_Ntrt_pft_emmeans[6, 3]) - exp(chi_lmer_Ntrt_pft_emmeans[5, 3]))/
#                        exp(chi_lmer_Ntrt_pft_emmeans[5, 3])) * 100
chi_change_N_c3nofix <- ((exp(chi_lmer_Ntrt_pft_emmeans[2, 3]) - exp(chi_lmer_Ntrt_pft_emmeans[1, 3]))/
                            exp(chi_lmer_Ntrt_pft_emmeans[2, 3])) * 100
chi_change_N_c3yesfix <- ((exp(chi_lmer_Ntrt_pft_emmeans[4, 3]) - exp(chi_lmer_Ntrt_pft_emmeans[3, 3]))/
                             exp(chi_lmer_Ntrt_pft_emmeans[3, 3])) * 100

#### Hyp 2: soil N (-) and C4 (-) confirmed...no impact of Nfix, T, or vpd

cor.test(log(leaf_n$beta), log(leaf_n$nmass), use = "pairwise.complete.obs")
cor.test(logit(leaf_n$chi), log(leaf_n$nmass), use = "pairwise.complete.obs")

#### Hyp 3: Nmass positively related to light, vpd, soil N, C4, and Nfix
#### negatively related to temperature and soil water
# nmass_lmer <- lmer(log(nmass) ~ Ntrt_fac + 
#                      p_pet +
#                      Nfix + photosynthetic_pathway +
#                      tmp2_gs + vpd2_gs +
#                      par2_gs +
#                      (1|Taxon) + (1|site_code) + (1|site_code:block_fac) +
#                      (1|site_code:block_fac:PKtrt), 
#                   data = leaf_n)
# plot(resid(nmass_lmer) ~ fitted(nmass_lmer))
# summary(nmass_lmer) # N = 2009
# Anova(nmass_lmer, type = 'II')
# nmass_mfx <- marginaleffects(nmass_lmer)
# summary(nmass_mfx)

#### Hyp 3: soil N (+), light (+), Nfix (+; separate from effect on beta) confirmed
#### C4 effect is negative (separate from effect on beta)
#### no temperature vpd or soil moisture effect

#### Hyp 4: Marea positively related to light
#### negatively related to temperature
# marea_lmer <- lmer(log(lma) ~ Ntrt_fac + 
#                      p_pet +
#                      tmp2_gs + par2_gs +
#                      (1|Taxon) + (1|site_code) + (1|site_code:block_fac) +
#                      (1|site_code:block_fac:trt), # note that trt (not PKtrt) used here
#                   data = leaf_n)
# plot(resid(marea_lmer) ~ fitted(marea_lmer))
# summary(marea_lmer) # N = 1768
# Anova(marea_lmer, type = 'II')
# marea_mfx <- marginaleffects(marea_lmer)
# summary(marea_mfx)

#### Hyp 4: no effects except separate negative effect of soil N?, so we should expect same impact on Nmass as to Narea

#### Hyp 5: Narea positively related to light, vpd, soil N, C4, and Nfix
#### negatively related to temperature and soil water
# narea_lmer <- lmer(log(narea) ~ Ntrt_fac + 
#                      p_pet +
#                      Nfix + photosynthetic_pathway +
#                      tmp2_gs + vpd2_gs + par2_gs +
#                      (1|Taxon) + (1|site_code) + (1|site_code:block_fac) +
#                      (1|site_code:block_fac:trt), 
#                   data = leaf_n)
# plot(resid(narea_lmer) ~ fitted(narea_lmer))
# summary(narea_lmer) # N = 1766
# Anova(narea_lmer, type = 'II')
# narea_mfx <- marginaleffects(narea_lmer)
# summary(narea_mfx)

#### Hyp 5: only a positive soil N effect and reduction in C4, but others washed out

#### hypothesis 1 and 2 plot

# beta_lmer_Ntrt_emmeans <- summary(emmeans(beta_lmer, ~Ntrt_fac))
# beta_lmer_Ntrt_emmeans_df <- as.data.frame(beta_lmer_Ntrt_emmeans)
# 
# chi_lmer_Ntrt_emmeans <- summary(emmeans(chi_lmer, ~Ntrt_fac))
# chi_lmer_Ntrt_emmeans_df <- as.data.frame(chi_lmer_Ntrt_emmeans)
# 
# beta_lmer_photosynthetic_pathway_emmeans <- summary(emmeans(beta_lmer, ~photosynthetic_pathway))
# beta_lmer_photosynthetic_pathway_emmeans_df <- as.data.frame(beta_lmer_photosynthetic_pathway_emmeans)
# 
# chi_lmer_photosynthetic_pathway_emmeans <- summary(emmeans(chi_lmer, ~photosynthetic_pathway))
# chi_lmer_photosynthetic_pathway_emmeans_df <- as.data.frame(chi_lmer_photosynthetic_pathway_emmeans)

beta_Ntrt_plot <- ggplot(aes(x = pft, y = log(beta), fill = Ntrt_fac), data = leaf_n) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(2), colour = 'black'),
        axis.title.x = element_text(size = rel(2), colour = 'black'),
        axis.text.x = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_violin(col = 'black', adjust = 2) +
  scale_fill_manual(values = c('grey', 'orange')) +
  guides(fill = "none") +
  # geom_boxplot(width = 0.5, fill = 'white') +
  # geom_point(data = beta_lmer_Ntrt_pft_emmeans, aes(y = emmean), size = 6) +
  # geom_errorbar(data = beta_lmer_Ntrt_emmeans, 
  #               aes(y= emmean, ymin = emmean - SE, ymax = emmean + SE, x = Ntrt_fac), 
  #               width = 0.2) + 
  scale_x_discrete(labels = c(expression('non-N-fixer'), 
                              expression('N-fixer'))) +
  ylim(c(4, 8)) +
  ylab('ln(β) (unitless)') +
  xlab('')

chi_Ntrt_plot <- ggplot(aes(x = pft, y = logit(chi), fill = Ntrt_fac), data = leaf_n) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(2), colour = 'black'),
        axis.title.x = element_text(size = rel(2), colour = 'black'),
        axis.text.x = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_violin(col = 'black', adjust = 2) +
  scale_fill_manual(values = c('grey', 'orange')) +
  guides(fill = "none") +
  # geom_boxplot(width = 0.5, fill = 'white', outlier.color = NA) +
  scale_x_discrete(labels = c(expression('C'[3] * ' non-N-fixer'), 
                              expression('C'[3] * ' N-fixer'))) +
  # ylim(c(0, 2)) +
  ylab('logit(χ)') +
  xlab('')

# beta_photosynthetic_pathway_plot <- ggplot(aes(x = photosynthetic_pathway, y = log(beta), fill = photosynthetic_pathway), data = leaf_n) +
#   theme(legend.position = NULL,
#         axis.title.y = element_text(size = rel(2), colour = 'black'),
#         axis.title.x = element_text(size = rel(2), colour = 'black'),
#         axis.text.x = element_text(size = rel(2)),
#         axis.text.y = element_text(size = rel(2)),
#         panel.background = element_rect(fill = 'white', colour = 'black'),
#         panel.grid.major = element_line(colour = "white"),
#         legend.background = element_blank(),
#         legend.box.background = element_rect(colour = "black")) +
#   geom_violin(col = 'black', adjust = 2) +
#   scale_fill_manual(values = c('grey', 'orange')) +
#   guides(fill = "none") +
#   geom_boxplot(width = 0.1, fill = 'black', outlier.color = NA) +
#   scale_x_discrete(labels = c(expression('C'[3]), expression('C'[4]))) +
#   ylim(c(-2, 10)) +
#   ylab('ln(β) (unitless)') +
#   xlab('Photosynthetic pathway')
# 
# chi_photosynthetic_pathway_plot <- ggplot(aes(x = photosynthetic_pathway, y = logit(chi), fill = photosynthetic_pathway), data = leaf_n) +
#   theme(legend.position = NULL,
#         axis.title.y = element_text(size = rel(2), colour = 'black'),
#         axis.title.x = element_text(size = rel(2), colour = 'black'),
#         axis.text.x = element_text(size = rel(2)),
#         axis.text.y = element_text(size = rel(2)),
#         panel.background = element_rect(fill = 'white', colour = 'black'),
#         panel.grid.major = element_line(colour = "white"),
#         legend.background = element_blank(),
#         legend.box.background = element_rect(colour = "black")) +
#   geom_violin(col = 'black', adjust = 2) +
#   scale_fill_manual(values = c('grey', 'orange')) +
#   guides(fill = "none") +
#   geom_boxplot(width = 0.1, fill = 'black', outlier.color = NA) +
#   scale_x_discrete(labels = c(expression('C'[3]), expression('C'[4]))) +
#   ylim(c(-4, 4)) +
#   ylab('logit(χ)') +
#   xlab('Photosynthetic pathway')
# 
# jpeg(filename = "plots/beta_chi.jpeg", width = 14, height = 14, units = 'in', res = 600)
# multiplot(beta_Ntrt_plot, beta_photosynthetic_pathway_plot, 
#           chi_Ntrt_plot, chi_photosynthetic_pathway_plot, cols=2)
# dev.off()
# 
# jpeg(filename = "plots/beta_chi.jpeg", width = 10, height = 6, units = 'in', res = 600)
# multiplot(beta_Ntrt_plot, 
#           chi_Ntrt_plot, cols=2)
# dev.off()

# make plots of changes by sites
leaf_n_site_Ntrt_pft_groupby <- group_by(leaf_n,
                                         site_code,
                                         pft,
                                         Ntrt_fac)

leaf_n_site_Ntrt_pft_summary <- summarise(leaf_n_site_Ntrt_pft_groupby, 
                                   beta_mean = mean(beta, na.rm = T),
                                   chi_mean = mean(chi, na.rm = T),
                                   Ntrt = mean(Ntrt, na.rm = T))

beta_Ntrt_site_plot_c3_noNfix <- ggplot(aes(x = Ntrt, y = (beta_mean)), 
                              data = subset(leaf_n_site_Ntrt_pft_summary, pft == 'c3_noNfix')) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(2), colour = 'black'),
        axis.title.x = element_text(size = rel(2), colour = 'black'),
        axis.text.x = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_point(pch = 16, color="black", alpha = 0.3, size = 3) +
  geom_line(aes(group = site_code), color="black", alpha = 0.3, lwd = 1) +
  geom_point(aes(x = as.numeric(as.character(Ntrt_fac)), y = exp(emmean)), 
             data = subset(beta_lmer_Ntrt_pft_emmeans, pft == 'c3_noNfix'),
             color="black", alpha = 0.7, size = 10) +
  geom_line(aes(x = as.numeric(as.character(Ntrt_fac)), y = exp(emmean)), 
            data = subset(beta_lmer_Ntrt_pft_emmeans, pft == 'c3_noNfix'),
            color="black", alpha = 0.7, lwd = 3) +
  ylab('β (unitless)') +
  xlab('N treatment') +
  scale_x_continuous(breaks = c(0, 1), labels=c('0' = 'Ambient', '1' = 'Added'))

beta_Ntrt_site_plot_c3_yesNfix <- ggplot(aes(x = Ntrt, y = (beta_mean)), 
                              data = subset(leaf_n_site_Ntrt_pft_summary, pft == 'c3_yesNfix')) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(2), colour = 'black'),
        axis.title.x = element_text(size = rel(2), colour = 'black'),
        axis.text.x = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_point(pch = 16, color="black", alpha = 0.3, size = 3) +
  geom_line(aes(group = site_code), color="black", alpha = 0.3, lwd = 1) +
  geom_point(aes(x = as.numeric(as.character(Ntrt_fac)), y = exp(emmean)), 
             data = subset(beta_lmer_Ntrt_pft_emmeans, pft == 'c3_yesNfix'),
             color="black", alpha = 0.7, size = 10) +
  geom_line(aes(x = as.numeric(as.character(Ntrt_fac)), y = exp(emmean)), 
            data = subset(beta_lmer_Ntrt_pft_emmeans, pft == 'c3_yesNfix'),
            color="black", alpha = 0.7, lwd = 3) +
  ylab('β (unitless)') +
  xlab('N treatment') +
  scale_x_continuous(breaks = c(0, 1), labels=c('0' = 'Ambient', '1' = 'Added'))


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
chi_pred_lm <- lm(chi ~ chi_pred + Ntrt_fac + p_pet, data = leaf_n)
plot(resid(chi_pred_lm) ~ fitted(chi_pred_lm))
summary(chi_pred_lm) # N = 
Anova(chi_pred_lm, type = 'II')
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
emtrends(chi_pred_beta_lm, ~1, var = 'chi_pred_beta')

plot(chi ~ chi_pred_beta, data = leaf_n)
abline(0,1)

## prediction plot
chi_beta_pred_plot <- ggplot(data = leaf_n, aes(x = chi_pred_beta, y = chi, color = Ntrt_fac)) +
  theme(legend.position = NULL,
        axis.title.y = element_text(size = rel(2), colour = 'black'),
        axis.title.x = element_text(size = rel(2), colour = 'black'),
        axis.text.x = element_text(size = rel(2)),
        axis.text.y = element_text(size = rel(2)),
        legend.title = element_text(size = rel(1)),
        legend.text = element_text(size = rel(1)),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_point(alpha = 0.2) +
  scale_color_manual(values = c('grey', 'orange'), labels = c('Ambient', 'Added')) +
  geom_abline(slope=1, intercept=0, lty = 2) +
  geom_smooth(data = leaf_n, aes(x = chi_pred_beta, y = chi), 
              method = 'lm', se= T, inherit.aes = F, color = 'black') +
  labs(color = 'N treatment') +
  xlab('Predicted χ') +
  ylab('Observed χ') +
  xlim(c(0.5, 0.9)) +
  ylim(c(0,1))

jpeg(filename = "plots/chi_pred.jpeg", width = 7, height = 7, units = 'in', res = 600)
plot(chi_beta_pred_plot)
dev.off()






## SEM prediction of chi
leaf_n$Ntrt_cont[leaf_n$Ntrt_fac == '0'] <- 0
leaf_n$Ptrt_cont[leaf_n$Ptrt_fac == '0'] <- 0
leaf_n$Ntrt_cont[leaf_n$Ntrt_fac == '1'] <- 1
leaf_n$Ptrt_cont[leaf_n$Ptrt_fac == '1'] <- 1
leaf_n$Nfix_cont[leaf_n$Nfix == 'no'] <- 0
leaf_n$Nfix_cont[leaf_n$Nfix == 'yes'] <- 1
leaf_n$photosynthetic_pathway_cont[leaf_n$photosynthetic_pathway == 'C3'] <- 0
leaf_n$photosynthetic_pathway_cont[leaf_n$photosynthetic_pathway == 'C4'] <- 1

leaf_n$scale.beta <- scale(log(leaf_n$beta))
leaf_n$scale.Ntrt_cont <- scale(leaf_n$Ntrt_cont)
leaf_n$scale.Ptrt_cont <- scale(leaf_n$Ptrt_cont)
leaf_n$scale.Nfix_cont <- scale(leaf_n$Nfix_cont)
leaf_n$scale.photosynthetic_pathway_cont <- scale(leaf_n$photosynthetic_pathway_cont)
leaf_n$scale.p_pet <- scale(leaf_n$p_pet)
leaf_n$scale.chi <- scale(logit(leaf_n$chi))
leaf_n$scale.vpd2_gs <- scale(leaf_n$vpd2_gs)
leaf_n$scale.tmp2_gs <- scale(leaf_n$tmp2_gs)
leaf_n$scale.par2_gs <- scale(leaf_n$par2_gs)
leaf_n$scale.narea <- scale(leaf_n$narea)
leaf_n$scale.nmass <- scale(leaf_n$nmass)
leaf_n$scale.lma <- scale(leaf_n$lma)

chi_path <- 'scale.beta ~ scale.Ntrt_cont + scale.p_pet + scale.photosynthetic_pathway_cont
             scale.chi ~ scale.beta'

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



