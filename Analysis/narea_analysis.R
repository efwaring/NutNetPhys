# Script using least-cost to explore Narea response at nutnet
## generally follows the analysis structure of Dong et al. (2017) Biogeosciences
## but includes the addition of soil N as a predictor


#### load packages ####
library(tidyverse)
library(lme4)
library(car)
library(r2glmm)
library(treemapify)
library(emmeans)
library(relaimpo)
library(patchwork)
library(multcomp)

#### load functions ####
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


#### create hypothesis figure ####
### leaf N and chi by N supply with different differences in N demand

n_supply_trend <- calc_optimal_vcmax(beta = seq(100, 300, 10))
# n_supply_trend$photo_n <- fvcmax25_nrubisco(n_supply_trend$vcmax) + fjmax25_nbioe(n_supply_trend$jmax)
n_supply_trend$photo_n_nobetachange <- n_supply_trend$nphoto[10]

hypothesis_data <- data.frame(cbind(c(n_supply_trend$beta, n_supply_trend$beta), 
                        c(n_supply_trend$nphoto, n_supply_trend$photo_n_nobetachange),
                        c(rep('no_change', 21), rep('change', 21))))
colnames(hypothesis_data) <- c('beta', 'Narea', 'demand')
hypothesis_data$beta <- as.numeric(as.character(hypothesis_data$beta))
hypothesis_data$Narea <- as.numeric(as.character(hypothesis_data$Narea))

(hypothesis_plot <- ggplot(data = hypothesis_data, 
                         aes(x = (1/beta), y = Narea, col = demand)) +
  theme(legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.text = element_text(size = 25),
        axis.title.y = element_text(size = rel(4), colour = 'black'),
        axis.title.x = element_text(size = rel(4), colour = 'black'),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line(size = 4, aes(linetype = demand, color = demand)) +
  scale_linetype_manual(values = c(1, 2), guide = F) +
  scale_colour_manual(values = c('black', 'grey'), 
                      labels = c('∆AGB = high', '∆AGB = low')) +
  guides(color = guide_legend(title = NULL)) +
  ylab(expression(italic('N')['area'])) +
  xlab(expression(italic('N')['availability'])))

#### soil N effects on leaf traits (Narea and chi) ####
### load data
leaf <- read.csv('../Data/processed/traits_v2.csv')

## turn treatment numbers into factors
leaf$Ntrt_fac <- as.factor(leaf$Ntrt)
leaf$Ptrt_fac <- as.factor(leaf$Ptrt)
leaf$Ktrt_fac <- as.factor(leaf$Ktrt)
leaf$block_fac <- as.factor(leaf$block)

## scale temperature data
leaf$tmp_scaled <- leaf$tmp - 25

## assign plant functional groups
leaf$pft <- leaf$functional_group
leaf$pft[leaf$pft == 'NULL'] <- NA
leaf$grass <- 'no'
leaf$grass[leaf$pft == 'GRASS'] <- 'yes'
leaf$fence <- 'no'
leaf$fence[leaf$trt == 'Fence' | leaf$trt == 'NPK+Fence'] <- 'yes'

## create subset of data where chi is greater than 0 and less than 1
leaf_chi_subset = subset(leaf, chi > 0 & chi < 1) # lose 659 points
# leaf_chi_subset1$chi_mad <- abs(leaf_chi_subset1$chi - median(leaf_chi_subset1$chi)) /
#   mad(leaf_chi_subset1$chi)
# hist(leaf_chi_subset1$chi_mad)
# nrow(leaf_chi_subset1)
# leaf_chi_subset <- subset(leaf_chi_subset1, chi_mad < 5)
# nrow(leaf_chi_subset)

### linear mixed effects model
leaf_chi_subset$logpar <- log(leaf_chi_subset$par)
leaf_chi_subset$logpar_per_leaf_area <- log(leaf_chi_subset$par_per_leaf_area)
leaf_chi_subset$logvpd <- log(leaf_chi_subset$vpd)
leaf_chi_subset$loglma <- log(leaf_chi_subset$lma)
leaf$logpar <- log(leaf$par)
leaf$logpar_per_leaf_area <- log(leaf$par_per_leaf_area)
leaf$logvpd <- log(leaf$vpd)
leaf$loglma <- log(leaf$lma)
leafNarea_lmer <- lmer(log(narea) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac + tmp + 
                         logpar_per_leaf_area + loglma + chi + Nfix + photosynthetic_pathway +
                        (1|Taxon) + (1|Taxon:site_code) + (1|Taxon:site_code:block_fac), 
                      data = leaf_chi_subset)
# plot(resid(leafNarea_lmer) ~ fitted(leafNarea_lmer))
summary(leafNarea_lmer) # N = 1,561
Anova(leafNarea_lmer)

# pwpm(emmeans(leafNarea_lmer, ~Ntrt_fac*Ptrt_fac))
# cld(emmeans(leafNarea_lmer, ~Ntrt_fac*Ptrt_fac))

### get some stats
#### soil nitrogen effect
# percentage increase of Narea in plots receiving N compared to plots not receiving N
(summary(emmeans(leafNarea_lmer, ~Ntrt_fac))[2,2] - summary(emmeans(leafNarea_lmer, ~Ntrt_fac))[1,2])/
  summary(emmeans(leafNarea_lmer, ~Ntrt_fac))[1,2]
# 0.286

# percentage increase of Narea in plots receiving N but not P compared to plots not receiving N or P
(summary(emmeans(leafNarea_lmer, ~Ntrt_fac * Ptrt_fac))[2,3] - summary(emmeans(leafNarea_lmer, ~Ntrt_fac * Ptrt_fac))[1,3])/
  summary(emmeans(leafNarea_lmer, ~Ntrt_fac * Ptrt_fac))[1,3]
# 0.352

# percentage increase of Narea in plots receiving N and P compared to plots receiving P but not N
(summary(emmeans(leafNarea_lmer, ~Ntrt_fac * Ptrt_fac))[4,3] - summary(emmeans(leafNarea_lmer, ~Ntrt_fac * Ptrt_fac))[3,3])/
  summary(emmeans(leafNarea_lmer, ~Ntrt_fac * Ptrt_fac))[3,3]
# 0.225

# percentage increase of Narea in N fixers compared to non-N fixers
(summary(emmeans(leafNarea_lmer, ~Nfix))[2,2] - summary(emmeans(leafNarea_lmer, ~Nfix))[1,2])/
  summary(emmeans(leafNarea_lmer, ~Nfix))[1,2]
# 1.022

# percentage increase of Narea in C3s compared to C4s
(summary(emmeans(leafNarea_lmer, ~photosynthetic_pathway))[1,2] - summary(emmeans(leafNarea_lmer, ~photosynthetic_pathway))[2,2])/
  summary(emmeans(leafNarea_lmer, ~photosynthetic_pathway))[1,2]
# 0.516

### make figures
## find slope and intercept from mixed effects model
tmp_slope <- summary(emtrends(leafNarea_lmer, ~tmp, var = "tmp"))[1, 2] # slope = -0.0276
tmp_intercept <- summary(emmeans(leafNarea_lmer, ~tmp, at = list(tmp = 0)))[1, 2] # intercept = 1.75
tmp_seq <- seq(min(leaf$tmp, na.rm = T), max(leaf$tmp, na.rm = T), 0.01)
tmp_trend <- tmp_intercept + tmp_seq * tmp_slope
tmp_trend <- as.data.frame(cbind(tmp_seq, tmp_trend))

(tmp_plot <- ggplot(data = leaf_chi_subset, aes(x = tmp, y = log(narea))) + 
    geom_jitter(pch = 21, fill = "black", alpha = 0.8) + 
    geom_line(data = tmp_trend, aes(x = tmp_seq, y = tmp_trend), 
              col = 'black', lwd = 2, alpha = 0.8) +
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 30, colour = 'black'),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 20, colour = 'black'),
          axis.text.y = element_text(size = 20, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "grey")) +
    ylab(expression('ln ' * italic('N')['area'])) +
    xlab('Mean Annual Growing Season Temperature (°C)'))

## find slope and intercept from mixed effects model
lma_slope <- summary(emtrends(leafNarea_lmer, ~loglma, var = "loglma"))[1, 2] # slope = 0.936
lma_intercept <- summary(emmeans(leafNarea_lmer, ~loglma, at = list(loglma = 0)))[1, 2] # intercept = -3.32
lma_seq <- seq(min(leaf_chi_subset$loglma, na.rm = T), max(leaf_chi_subset$loglma, na.rm = T), 0.01)
lma_trend <- lma_intercept + lma_seq * lma_slope
lma_trend <- as.data.frame(cbind(lma_seq, lma_trend))

(lma_plot <- ggplot(data = leaf_chi_subset, aes(x = log(lma), y = log(narea))) + 
    geom_jitter(pch = 21, fill = "black", alpha = 0.8) + 
    geom_line(data = lma_trend, aes(x = lma_seq, y = lma_trend), 
              col = 'black', lwd = 2, alpha = 0.8) +
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 30, colour = 'black'),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 20, colour = 'black'),
          axis.text.y = element_text(size = 20, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "grey")) +
    ylab(expression('ln ' * italic('N')['area'])) +
    xlab(expression('ln ' * italic('M')['area'])))

## find slope and intercept from mixed effects model
par_slope <- summary(emtrends(leafNarea_lmer, ~logpar_per_leaf_area, var = "logpar_per_leaf_area"))[1, 2] # slope = 0.936
par_intercept <- summary(emmeans(leafNarea_lmer, ~logpar_per_leaf_area, at = list(logpar_per_leaf_area = 0)))[1, 2] # intercept = -3.32
par_seq <- seq(min(leaf_chi_subset$logpar_per_leaf_area, na.rm = T), max(leaf_chi_subset$logpar_per_leaf_area, na.rm = T), 0.01)
par_trend <- par_intercept + par_seq * par_slope
par_trend <- as.data.frame(cbind(par_seq, par_trend))

(par_plot <- ggplot(data = leaf_chi_subset, aes(x = logpar_per_leaf_area, y = log(narea))) + 
    geom_jitter(pch = 21, fill = "black", alpha = 0.8) + 
    geom_line(data = par_trend, aes(x = par_seq, y = par_trend), 
              col = 'black', lwd = 2, alpha = 0.8) +
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 30, colour = 'black'),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 20, colour = 'black'),
          axis.text.y = element_text(size = 20, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "grey")) +
    ylab(expression('ln ' * italic('N')['area'])) +
    xlab(expression('ln ' * italic('I')['g'])))

### assign treatment group labels
leaf_chi_subset$PKgroup[leaf_chi_subset$Ptrt_fac == '0' & leaf_chi_subset$Ktrt_fac == '0'] <- '-P, -K'
leaf_chi_subset$PKgroup[leaf_chi_subset$Ptrt_fac == '1' & leaf_chi_subset$Ktrt_fac == '0'] <- '+P, -K'
leaf_chi_subset$PKgroup[leaf_chi_subset$Ptrt_fac == '0' & leaf_chi_subset$Ktrt_fac == '1'] <- '-P, +K'
leaf_chi_subset$PKgroup[leaf_chi_subset$Ptrt_fac == '1' & leaf_chi_subset$Ktrt_fac == '1'] <- '+P, +K'
    
cld.emmGrid(emmeans(leafNarea_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))
leafNarea_letters <- data.frame(x = c(0.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2),
                          NPgroup = c('-P, -K', '-P, -K', '+P, -K', '+P, -K', 
                                      '-P, +K', '-P, +K', '+P, +K', '+P, +K'),
                          Ntrt_fac = c(0, 1, 0, 1, 0, 1, 0, 1),
                          y = c(3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5, 3.5), 
                          letter = c('a', 'b', 'a', 'b', 
                                     'a', 'b', 'a', 'b'))
leafNarea_letters$Ntrt_fac <- as.factor(leafNarea_letters$Ntrt_fac)

(narea_plot <- ggplot(data = leaf_chi_subset, 
                         aes(x = PKgroup, y = log(narea), fill = Ntrt_fac)) +
  theme(legend.position = "right",
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 15),
        legend.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.y = element_text(size = 30, colour = 'black'),
        axis.title.x = element_text(size = 30, colour = 'black'),
        axis.text.x = element_text(size = 20, colour = 'black'),
        axis.text.y = element_text(size = 20, colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "white")) +
    geom_boxplot(outlier.color = NA) +
    geom_text(data = leafNarea_letters, aes(x = x, y = y, label = letter), size = 6) +
    scale_fill_manual(values = c("gray40", "burlywood1"), labels = c("Ambient", "Added N")) +
    # geom_jitter(size = 0.1) +
    labs(fill = "Soil N") +
    ylab(expression('ln ' * italic('N')['area'])) +
    xlab('P x K treatment'))

## find relative importance for each factor from model
relimp_leafn <- NULL
relimp_leafn$Factor <- c('Soil N', 'Soil P', 'Soil K+µ', 'Temperature', 'PAR', 'LMA', 'χ',
                         'N fixer', 'C3/C4', 'Soil Interactions', 'Unexplained')
relimp_leafn$Importance <- as.numeric(as.character(c(calc.relip.mm(leafNarea_lmer)$lmg[1:9], 
                                                        sum(calc.relip.mm(leafNarea_lmer)$lmg[10:13]),
                                                        1 - sum(calc.relip.mm(leafNarea_lmer)$lmg))))

relimp_leafn_df <- as.data.frame(relimp_leafn)

tm <- treemapify(data = relimp_leafn_df,
                 area = "Importance", start = "topleft")
tm$x <- (tm$xmax + tm$xmin) / 2
tm$y <- (tm$ymax + tm$ymin) / 2

narea_tm <- full_join(relimp_leafn_df, tm, by = "Factor")
narea_tm$name <- c('Soil~N', 'Soil~P', 'Soil~K[+µ]', 'italic(T)[g]', 'italic(I)[g]',
                   'italic(M)[area]', 'χ', 'N~fixer', 
                   'C[3]/C[4]', 'Soil~Interactions', 'Unexplained')

(narea_treemap <- ggplot(narea_tm, 
                         aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                             label = name)) +
    geom_rect(aes(fill = Importance), color = "black") +
    theme(legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "right",
          panel.background = element_rect(fill = 'white'),
          axis.title = element_text(colour = 'white'),
          axis.text = element_text(colour = 'white'),
          axis.ticks = element_line(colour = "white")) + 
    scale_fill_gradient(low = "lightcyan", high = "lightcyan4") +
    geom_text(data = filter(narea_tm, Factor == 'LMA' | Factor == 'PAR'), 
              aes(x = x, y = y), parse = T, size = 14) +
    # geom_text(data = filter(narea_tm, ), 
    #           aes(x = x, y = y), parse = T, size = 8) +
    geom_text(data = filter(narea_tm, Factor == 'χ'), 
              aes(x = x, y = y), parse = T, size = 10, family = 'Times') +
    geom_text(data = filter(narea_tm, Factor == 'Unexplained' | Factor == 'Temperature' | 
                              Factor == 'N fixer' | Factor == 'C3/C4'), 
              aes(x = x, y = y), parse = T, size = 7) +
    geom_text(data = filter(narea_tm, Factor == 'Soil N'), 
              aes(x = x, y = y), parse = T, size = 4) +
    ggrepel::geom_text_repel(data = filter(narea_tm, Factor == 'Soil Interactions' | Factor == 'Soil P' |
                                             Factor == 'Soil K+µ'), 
                             aes(x = x, y = y), parse = T, size = 4, 
                             direction = "y", xlim = c(1.01, NA)) +
    scale_x_continuous(limits = c(0, 1.2), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)))

# (narea_plot_treemap <- narea_plot + narea_treemap +
#   plot_annotation(tag_levels = 'A') & 
#   theme(plot.tag = element_text(size = 24)))
  
### table with model results
Narea_model <- data.frame(Var = c('Soil N', 'Soil P', 'Soil K+µ', 'Temperature', 
                             'ln PAR', 'ln LMA', 'χ', 'N fixer', 'C3/C4',
                             'Soil N x Soil P', 'Soil N x Soil K', 'Soil P x Soil K',
                             'Soil N x Soil P x Soil K'))
Narea_model$df <- as.matrix(Anova(leafNarea_lmer))[1:13, 2]
Narea_model$Slope <- c(NA, NA, NA,
                      summary(emtrends(leafNarea_lmer, ~tmp, var = "tmp"))[1, 2],
                      summary(emtrends(leafNarea_lmer, ~logpar_per_leaf_area, var = "logpar_per_leaf_area"))[1, 2],
                      summary(emtrends(leafNarea_lmer, ~loglma, var = "loglma"))[1, 2],
                      summary(emtrends(leafNarea_lmer, ~chi, var = "chi"))[1, 2],
                      NA, NA, NA, NA, NA, NA)
Narea_model$SE <- c(NA, NA, NA,
                       summary(emtrends(leafNarea_lmer, ~tmp, var = "tmp"))[1, 3],
                       summary(emtrends(leafNarea_lmer, ~logpar_per_leaf_area, var = "logpar_per_leaf_area"))[1, 3],
                       summary(emtrends(leafNarea_lmer, ~loglma, var = "loglma"))[1, 3],
                       summary(emtrends(leafNarea_lmer, ~chi, var = "chi"))[1, 3],
                       NA, NA, NA, NA, NA, NA)
Narea_model$p <- as.matrix(Anova(leafNarea_lmer))[1:13, 3]
Narea_model$RelImp <- as.matrix(calc.relip.mm(leafNarea_lmer)$lmg)[1:13]
Narea_model$RelImp <- Narea_model$RelImp * 100

write.csv(Narea_model, 'tables/Narea_model.csv')


#### Narea predictions ####
### calculate vcmax, jmax, and vpmax with known chi for C3 plants
leaf_chi_subset_c3 <- subset(leaf_chi_subset, photosynthetic_pathway == 'C3')
gas_exchange_pred_c3 <- calc_optimal_vcmax(pathway = "C3",
                                           tg_c = leaf_chi_subset_c3$tmp, 
                                           paro = leaf_chi_subset_c3$par_per_leaf_area, 
                                           cao = 400, 
                                           vpdo = leaf_chi_subset_c3$vpd, 
                                           z = leaf_chi_subset_c3$z,
                                           q0_resp = "no",
                                           chi = leaf_chi_subset_c3$chi,
                                           lma = leaf_chi_subset_c3$lma)

### calculate vcmax, jmax, and vpmax with known chi for C4 plants
leaf_chi_subset_c4 <- subset(leaf_chi_subset, photosynthetic_pathway == 'C4')
gas_exchange_pred_c4 <- calc_optimal_vcmax(pathway = "C4",
                                           tg_c = leaf_chi_subset_c4$tmp, 
                                           paro = leaf_chi_subset_c4$par_per_leaf_area, 
                                           cao = 400, 
                                           vpdo = leaf_chi_subset_c4$vpd, 
                                           z = leaf_chi_subset_c4$z,
                                           q0_resp = "no",
                                           chi = leaf_chi_subset_c4$chi,
                                           lma = leaf_chi_subset_c4$lma)

## add C3 model results to leaf dataset
npred_c3 <- bind_cols(leaf_chi_subset_c3, gas_exchange_pred_c3[ ,39:51])
npred_c3$model_lma <- gas_exchange_pred_c3$lma

## add C4 model results to leaf dataset
npred_c4 <- bind_cols(leaf_chi_subset_c4, gas_exchange_pred_c4[ ,39:51])
npred_c4$model_lma <- gas_exchange_pred_c4$lma

# combine C3 and C4 subsets
leaf_chi_subset_npred <- bind_rows(npred_c3, npred_c4)

leaf_chi_subset_npred$lognphoto <- log(leaf_chi_subset_npred$nphoto)
leaf_chi_subset_npred$lognstructure <- log(leaf_chi_subset_npred$nstructure)

### fit linear mixed effects model
npred_soil_lmer <- lmer(log(narea) ~ lognphoto + lognstructure +
                          Ntrt_fac * Ptrt_fac * Ktrt_fac +
                          Nfix + photosynthetic_pathway + 
                          (1|Taxon) + (1|Taxon:site_code) + (1|Taxon:site_code:block_fac),
                        data = leaf_chi_subset_npred)
# plot(resid(npred_soil_lmer) ~ fitted(npred_soil_lmer))
summary(npred_soil_lmer) # N = 1,548
Anova(npred_soil_lmer)
# cld(emmeans(npred_soil_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))

### make figures
## find slope and intercept from mixed effects model
nphoto_slope <- summary(emtrends(npred_soil_lmer, ~lognphoto, var = "lognphoto"))[1, 2]
nphoto_intercept_lowN <- summary(emmeans(npred_soil_lmer, ~Ntrt_fac, at = list(lognphoto = 0)))[1, 2]
nphoto_intercept_highN <- summary(emmeans(npred_soil_lmer, ~Ntrt_fac, at = list(lognphoto = 0)))[2, 2]
lognphoto_seq <- seq(min(leaf_chi_subset_npred$lognphoto, na.rm = T), 
                     max(leaf_chi_subset_npred$lognphoto, na.rm = T), 0.01)
nphoto_trend_lowN <- nphoto_intercept_lowN + lognphoto_seq * nphoto_slope
nphoto_trend_highN <- nphoto_intercept_highN + lognphoto_seq * nphoto_slope
nphoto_trend <- as.data.frame(cbind(lognphoto_seq, nphoto_trend_lowN, nphoto_trend_highN))

(npred_photo_plot <- ggplot(data = leaf_chi_subset_npred, 
                            aes(x = lognphoto, y = log(narea), color = Ntrt_fac)) +
    theme(legend.position = "none", 
          axis.title.y = element_text(size = 30, colour = 'black'),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 20, colour = 'black'),
          axis.text.y = element_text(size = 20, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "grey")) +
    geom_point(shape = 16, size = 3, alpha = 0.8) +
    scale_color_manual(values = c("black", "burlywood1"), labels = c("Ambient", "Added N")) +
    labs(color = 'Soil N') +
    geom_line(data = nphoto_trend, aes(x = lognphoto_seq, y = nphoto_trend_lowN), 
              col = 'black', lwd = 2, alpha = 0.8) +
    geom_line(data = nphoto_trend, aes(x = lognphoto_seq, y = nphoto_trend_highN), 
              col = 'burlywood1', lwd = 3, alpha = 0.8) +
    scale_y_continuous(limits = c(-2.5, 5.5)) +
    scale_x_continuous(limits = c(-2, 0)) +
    ylab(expression('ln ' * italic('N')['area'])) +
    xlab(expression('ln ' * italic('N')['photo'])))
  
## find slope and intercept from mixed effects model
nstructure_slope <- summary(emtrends(npred_soil_lmer, ~lognstructure, var = "lognstructure"))[1, 2]
nstructure_intercept_lowN <- summary(emmeans(npred_soil_lmer, ~Ntrt_fac, at = list(lognstructure = 0)))[1, 2]
nstructure_intercept_highN <- summary(emmeans(npred_soil_lmer, ~Ntrt_fac, at = list(lognstructure = 0)))[2, 2]
lognstructure_seq <- seq(min(leaf_chi_subset_npred$lognstructure, na.rm = T), 
                         max(leaf_chi_subset_npred$lognstructure, na.rm = T), 0.1)
nstructure_trend_lowN <- nstructure_intercept_lowN + lognstructure_seq * nstructure_slope
nstructure_trend_highN <- nstructure_intercept_highN + lognstructure_seq * nstructure_slope
nstructure_trend <- as.data.frame(cbind(lognstructure_seq, nstructure_trend_lowN, nstructure_trend_highN))

(npred_structure_plot <- ggplot(data = leaf_chi_subset_npred, 
                              aes(x = lognstructure, y = log(narea), color = Ntrt_fac)) +
    theme(legend.position = c(0, 1),
          legend.justification = c(0, 1),
          legend.title = element_text(size = 22),
          legend.text = element_text(size = 20),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          legend.key=element_blank(),
          axis.title.y = element_text(size = 30, colour = 'black'),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 20, colour = 'black'),
          axis.text.y = element_text(size = 20, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "grey")) +
    geom_point(shape = 16, size = 2, alpha = 0.8) +
    scale_color_manual(values = c("black", "burlywood1"), labels = c("Ambient", "Added N")) +
    labs(color = 'Soil N') +
    geom_line(data = nstructure_trend, aes(x = lognstructure_seq, y = nstructure_trend_lowN), 
              col = 'black', lwd = 3, alpha = 0.8) +
    geom_line(data = nstructure_trend, aes(x = lognstructure_seq, y = nstructure_trend_highN), 
              col = 'burlywood1', lwd = 3, alpha = 0.8) +
    scale_y_continuous(limits = c(-2.5, 5.5)) +
    scale_x_continuous(limits = c(-5, 2.75)) +
    ylab(expression('ln ' * italic('N')['area'])) +
    xlab(expression('ln ' * italic('N')['structure'])))

(npred_plot <- npred_structure_plot + npred_photo_plot +
    plot_annotation(tag_levels = 'A') & 
    theme(plot.tag = element_text(size = 18)))

## find relative importance for each factor from model
relimp_leafn_pred <- NULL
relimp_leafn_pred$Factor <- c('N photo', 'N structure', 'Soil N', 'Soil P', 'Soil K+µ', 
                              'N fixer', 'C3/C4', 'Soil Interactions', 'Unexplained')
relimp_leafn_pred$Importance <- as.numeric(as.character(c(calc.relip.mm(npred_soil_lmer)$lmg[1:7], 
                                                     sum(calc.relip.mm(npred_soil_lmer)$lmg[8:11]),
                                                     1 - sum(calc.relip.mm(npred_soil_lmer)$lmg))))
relimp_leafn_pred_df <- as.data.frame(relimp_leafn_pred)

sum(calc.relip.mm(npred_soil_lmer)$lmg[c(3:5, 8:11)])

tm_pred <- treemapify(data = relimp_leafn_pred_df,
                 area = "Importance", start = "topleft")
tm_pred$x <- (tm_pred$xmax + tm_pred$xmin) / 2
tm_pred$y <- (tm_pred$ymax + tm_pred$ymin) / 2

narea_tm_pred <- full_join(relimp_leafn_pred_df, tm_pred, by = "Factor")
narea_tm_pred$name <- c('italic(N)[photo]', 'italic(N)[structure]', 'Soil~N', 'Soil~P', 'Soil~K[+µ]', 
                        'N~fixer', 'C[3]/C[4]', 'Soil~Interactions', 'Unexplained')

(narea_pred_treemap <- ggplot(narea_tm_pred, 
                              aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, 
                                  label = name)) +
    geom_rect(aes(fill = Importance), color = "black") +
    theme(legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.position = "right",
          panel.background = element_rect(fill = 'white'),
          axis.title = element_text(colour = 'white'),
          axis.text = element_text(colour = 'white'),
          axis.ticks = element_line(colour = "white")) + 
    scale_fill_gradient(low = "lightcyan", high = "lightcyan4") +
    geom_text(data = filter(narea_tm_pred, Factor == 'Unexplained'),
              aes(x = x, y = y), parse = T, size = 14) +
    geom_text(data = filter(narea_tm_pred, Factor == 'N structure' | Factor == 'N photo' | 
                              Factor == 'Soil N'),
              aes(x = x, y = y), parse = T, size = 10) +
    geom_text(data = filter(narea_tm_pred, Factor == 'Soil K+µ' | Factor == 'Soil P' | 
                              Factor == 'N fixer'),
              aes(x = x, y = y), parse = T, size = 8) +
    geom_text(data = filter(narea_tm_pred, Factor == 'Soil Interactions'),
              aes(x = x, y = y), parse = T, size = 7) +
    ggrepel::geom_text_repel(data = filter(narea_tm_pred, Factor == 'C3/C4'), 
                             aes(x = x, y = y), parse = T, size = 4, 
                             direction = "y", xlim = c(1.01, NA)) +
    scale_x_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
    scale_y_continuous(expand = c(0, 0)))

### table with model results
Narea_pred_model <- data.frame(Var = c('ln Nphoto', 'ln Nstructure', 'Soil N', 'Soil P', 
                                       'Soil K+µ', 'N fixer', 'C3/C4', 'Soil N x Soil P', 
                                       'Soil N x Soil K', 'Soil P x Soil K', 
                                       'Soil N x Soil P x Soil K')) 
Narea_pred_model$df <- as.matrix(Anova(npred_soil_lmer))[1:11, 2]
Narea_pred_model$Slope <- c(summary(emtrends(npred_soil_lmer, ~lognphoto, var = "lognphoto"))[1, 2],
                            summary(emtrends(npred_soil_lmer, ~lognstructure, var = "lognstructure"))[1, 2],
                            NA, NA, NA, NA, NA, NA, NA, NA, NA)
Narea_pred_model$SE <- c(summary(emtrends(npred_soil_lmer, ~lognphoto, var = "lognphoto"))[1, 3],
                         summary(emtrends(npred_soil_lmer, ~lognstructure, var = "lognstructure"))[1, 3],
                         NA, NA, NA, NA, NA, NA, NA, NA, NA)
Narea_pred_model$p <- as.matrix(Anova(npred_soil_lmer))[1:11, 3]
Narea_pred_model$RelImp <- as.matrix(calc.relip.mm(npred_soil_lmer)$lmg)[1:11]
Narea_pred_model$RelImp <- Narea_pred_model$RelImp * 100

write.csv(Narea_pred_model, 'tables/Narea_pred_model.csv')

#### soil N effects on AGB and LAI ####
### load data
leaf_site <- read.csv('../Data/processed/traits_site.csv')

## turn treatment numbers into factors
leaf_site$Ntrt_fac <- as.factor(leaf_site$Ntrt)
leaf_site$Ptrt_fac <- as.factor(leaf_site$Ptrt)
leaf_site$Ktrt_fac <- as.factor(leaf_site$Ktrt)
leaf_site$block_fac <- as.factor(leaf_site$block)

## linear mixed effects model for leaf area index (lai)
lai_lmer <- lmer(log(lai_mean) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac +
                  (1|site_code) + (1|site_code:block), 
                data = leaf_site)
## plot(resid(lai_lmer) ~ fitted(lai_lmer))
summary(lai_lmer)
Anova(lai_lmer)

lai_model <- data.frame(Var = c('Soil N', 'Soil P', 'Soil K+µ', 'Soil N x Soil P', 
                                      'Soil N x Soil K', 'Soil P x Soil K', 
                                      'Soil N x Soil P x Soil K')) 
lai_model$chisq <- as.matrix(Anova(lai_lmer)[1:7, 1])
lai_model$df <- as.matrix(Anova(lai_lmer)[1:7, 2])
lai_model$p <- as.matrix(Anova(lai_lmer)[1:7, 3])
write.csv(lai_model, 'tables/lai_model.csv')

### linear mixed effects model for mean live mass
live_mass_lmer <- lmer(log(live_mass_mean) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac +
                        (1|site_code) + (1|site_code:block), 
                      data = leaf_site)
# plot(resid(live_mass_lmer) ~ fitted(live_mass_lmer))
summary(live_mass_lmer) # N = 763
Anova(live_mass_lmer)

# pwpm(emmeans(live_mass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))

# percentage increase of live mass in plots receiving N compared to plots not receiving N
(summary(emmeans(live_mass_lmer, ~Ntrt_fac))[2,2] - summary(emmeans(live_mass_lmer, ~Ntrt_fac))[1,2])/
  summary(emmeans(live_mass_lmer, ~Ntrt_fac))[1,2]
# 0.0411

# percentage increase of live mass in plots receiving N and K compared to plots receiving K but not N
(summary(emmeans(live_mass_lmer, ~Ntrt_fac*Ktrt_fac))[4,3] - summary(emmeans(live_mass_lmer, ~Ntrt_fac*Ktrt_fac))[3,3])/
  summary(emmeans(live_mass_lmer, ~Ntrt_fac*Ktrt_fac))[3,3]
# 0.055

live_mass_model <- data.frame(Var = c('Soil N', 'Soil P', 'Soil K+µ', 'Soil N x Soil P', 
                                       'Soil N x Soil K', 'Soil P x Soil K', 
                                       'Soil N x Soil P x Soil K')) 
live_mass_model$chisq <- as.matrix(Anova(live_mass_lmer)[1:7, 1])
live_mass_model$df <- as.matrix(Anova(live_mass_lmer)[1:7, 2])
live_mass_model$p <- as.matrix(Anova(live_mass_lmer)[1:7, 3])
write.csv(live_mass_model, 'tables/live_mass_model.csv')

### assign treatment group labels
leaf_site$PKgroup[leaf_site$Ptrt_fac == '0' & leaf_site$Ktrt_fac == '0'] <- '-P, -K'
leaf_site$PKgroup[leaf_site$Ptrt_fac == '1' & leaf_site$Ktrt_fac == '0'] <- '+P, -K'
leaf_site$PKgroup[leaf_site$Ptrt_fac == '0' & leaf_site$Ktrt_fac == '1'] <- '-P, +K'
leaf_site$PKgroup[leaf_site$Ptrt_fac == '1' & leaf_site$Ktrt_fac == '1'] <- '+P, +K'

### make figures
# cld.emmGrid(emmeans(live_mass_lmer, ~Ntrt_fac * Ptrt_fac * Ktrt_fac))
live_mass_letters <- data.frame(x = c(0.8, 1.2, 1.8, 2.2, 2.8, 3.2, 3.8, 4.2),
                                PKgroup = c('-P, -K', '-P, -K', '-P, +K', '-P, +K',
                                            '+P, -K', '+P, -K', '+P, +K', '+P, +K'),
                                Ntrt_fac = c(0, 1, 0, 1, 0, 1, 0, 1),
                                y = c(7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8, 7.8), 
                                letter = c("ab", "abc", "a", "c", 
                                           "bc", "cd", "bc", "d"))
live_mass_letters$Ntrt_fac <- as.factor(live_mass_letters$Ntrt_fac)

(live_mass_plot <- ggplot(data = leaf_site, 
                          aes(x = PKgroup, y = log(live_mass_mean), fill = Ntrt_fac)) +
    theme(legend.position = "right",
          legend.title = element_text(size = 20),
          legend.text = element_text(size = 15),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          axis.title.y = element_text(size = 30, colour = 'black'),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 20, colour = 'black'),
          axis.text.y = element_text(size = 20, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "white")) +
    geom_boxplot(outlier.color = NA) +
    scale_fill_manual(values = c("gray40", "burlywood1"), labels = c("Ambient", "Added N")) +
    geom_text(data = live_mass_letters, aes(x = x, y = y, label = letter), size = 6) +
    scale_y_continuous(limits = c(3.5, 8)) +
    scale_x_discrete(limits = c('-P, -K', '-P, +K', '+P, -K', '+P, +K')) +
    labs(fill = "Soil N") +
    xlab('P x K treatment') +
    ylab(expression('ln AGB (g)')))

#### supply vs demand effects on leaf N ####
### calculate treatment type averages
leaf_site_N <- leaf_chi_subset_npred %>%
  group_by(site_code, Ntrt_fac, Ptrt_fac, Ktrt_fac, 
           block_fac, 
           fence, trt,
           # plot, 
           Taxon 
           # grass, Nfix, photosynthetic_pathway
           ) %>%
  summarise_at(vars(narea, spp_lai, chi, spp_live_mass, spp_mass_N, max_cover, lma, p_pet, vpd, tmp,
                    Ambient_PAR, Ground_PAR, par_rat),
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
                                       'Ptrt_fac', 'Ktrt_fac',
                                       'fence',
                                       'Taxon' 
                                       # 'grass', 'Nfix', 'photosynthetic_pathway'
                                       ))

### calculate percent change from the high N plots to the low N plots
leaf_site_trt$delta_narea <- ((leaf_site_trt$narea.y - 
                                leaf_site_trt$narea.x) / leaf_site_trt$narea.x) * 100
leaf_site_trt$delta_lai <- ((leaf_site_trt$spp_lai.y - 
                              leaf_site_trt$spp_lai.x) / leaf_site_trt$spp_lai.x) * 100
leaf_site_trt$delta_live_mass <- ((leaf_site_trt$spp_live_mass.y - 
                                    leaf_site_trt$spp_live_mass.x) / leaf_site_trt$spp_live_mass.x) * 100
leaf_site_trt$delta_N_mass <- ((leaf_site_trt$spp_mass_N.y - 
                                    leaf_site_trt$spp_mass_N.x) / leaf_site_trt$spp_mass_N.x) * 100
leaf_site_trt$delta_chi <- ((leaf_site_trt$chi.y - 
                              leaf_site_trt$chi.x) / leaf_site_trt$chi.x) * 100
leaf_site_trt$delta_lma <- ((leaf_site_trt$lma.y - 
                              leaf_site_trt$lma.x) / leaf_site_trt$lma.x) * 100

nrow(subset(leaf_site_trt, delta_narea > -1000))

### find mean absolute deviation (Leys et al., 2013)
delta_narea_mad <- mad(leaf_site_trt$delta_narea, na.rm = T)
delta_lai_mad <- mad(leaf_site_trt$delta_lai, na.rm = T)
delta_lma_mad <- mad(leaf_site_trt$delta_lma, na.rm = T)
delta_live_mass_mad <- mad(leaf_site_trt$delta_live_mass, na.rm = T)
delta_N_mass_mad <- mad(leaf_site_trt$delta_N_mass, na.rm = T)
delta_chi_mad <- mad(leaf_site_trt$delta_chi, na.rm = T)

## remove instances where any ∆ values are 3 times higher than the MAD
# delta_lai_data <- subset(leaf_site_trt, 
#                          delta_narea < 3 * delta_narea_mad & 
#                           delta_narea > 3 * -delta_narea_mad & 
#                           delta_lma < 3 * delta_lma_mad &
#                           delta_lma > 3 * -delta_lma_mad &
#                           delta_lai < 3 * delta_lai_mad & 
#                           delta_lai > 3 * -delta_lai_mad &
#                           delta_chi < 3 * delta_chi_mad & 
#                           delta_chi > 3 * -delta_chi_mad)

delta_live_mass_data <- subset(leaf_site_trt,
                             delta_narea < 3 * delta_narea_mad &
                               delta_narea > 3 * -delta_narea_mad &
                               delta_lma < 3 * delta_lma_mad &
                               delta_lma > 3 * -delta_lma_mad &
                               delta_live_mass < 3 * delta_live_mass_mad &
                               delta_live_mass > 3 * -delta_live_mass_mad) #&
                               # delta_chi < 5 * delta_chi_mad &
                               # delta_chi > 5 * -delta_chi_mad)

# delta_live_mass_data <- leaf_site_trt

# delta_N_mass_data <- subset(leaf_site_trt, 
#                            delta_narea < 5 * delta_narea_mad & 
#                              delta_narea > 5 * -delta_narea_mad &
#                              delta_lma < 5 * delta_lma_mad &
#                              delta_lma > 5 * -delta_lma_mad &
#                              delta_N_mass < 5 * delta_N_mass_mad & 
#                              delta_N_mass > 5 * -delta_N_mass_mad) #&
#                              # delta_chi < 3 * delta_chi_mad &
#                              # delta_chi > 3 * -delta_chi_mad)

### linear mixed effects model for delta Narea by delta live mass
delta_live_mass_lm <- lmer(delta_narea ~ 
                            delta_live_mass + 
                            Ptrt_fac * Ktrt_fac +
                            # photosynthetic_pathway + Nfix +
                            delta_lma +
                            delta_live_mass : delta_lma +
                            # delta_live_mass * delta_chi +
                            # delta_live_mass * delta_lma * delta_chi +
                            (1|Taxon) + (1|Taxon:site_code) +
                            (1|Taxon:block_fac:site_code)
                           , 
                        data = delta_live_mass_data)
# plot(resid(delta_live_mass_lm) ~ fitted(delta_live_mass_lm))
summary(delta_live_mass_lm) # N = 328
Anova(delta_live_mass_lm)

test(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass'))
test(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass', 
                at = list(delta_lma = -25)))
test(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass', 
              at = list(delta_lma = 0)))
test(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass', 
             at = list(delta_lma = 25)))

emmeans(delta_live_mass_lm, ~Ptrt_fac)
# 
# test(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass', 
#               at = list(delta_lma = -25, delta_chi = -10)))
# test(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass', 
#               at = list(delta_lma = 0, delta_chi = -10)))
# test(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass', 
#               at = list(delta_lma = 25, delta_chi = -10)))
# 
# test(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass', 
#               at = list(delta_lma = -25, delta_chi = 0)))
# test(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass', 
#               at = list(delta_lma = 0, delta_chi = 0)))
# test(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass', 
#               at = list(delta_lma = 25, delta_chi = 0)))
# 
# test(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass', 
#               at = list(delta_lma = -25, delta_chi = 10)))
# test(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass', 
#               at = list(delta_lma = 0, delta_chi = 10)))
# test(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass', 
#               at = list(delta_lma = 25, delta_chi = 10)))
# 
# emmeans(delta_live_mass_lm, ~Ptrt_fac)

delta_live_mass_model <- data.frame(Var = c('delta AGB', 'Soil P', 'Soil K', 'C3/C4', 
                                            'N fixer', 'delta LMA', 'delta Chi', 'Soil P x Soil K',
                                            'delta AGB x delta LMA', 'delta AGB x delta Chi',
                                            'delta LMA x delta Chi', 'delta AGB x delta LMA x delta Chi')) 
delta_live_mass_model$chisq <- as.matrix(Anova(delta_live_mass_lm)[1:12, 1])
delta_live_mass_model$df <- as.matrix(Anova(delta_live_mass_lm)[1:12, 2])
delta_live_mass_model$p <- as.matrix(Anova(delta_live_mass_lm)[1:12, 3])
write.csv(delta_live_mass_model, 'tables/delta_live_mass_model.csv')

### make figures
## dataset
delta_live_mass_plot_data <- delta_live_mass_data

## trendline information
# low chi
dlm_lowlma_lowchi_intercept <- summary(emmeans(delta_live_mass_lm, ~delta_live_mass,
                                              at = list(delta_live_mass = 0, delta_lma = -25, delta_chi = -10)))[1, 2]
dlm_lowlma_lowchi_slope <- summary(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass',
                                           at = list(delta_lma = -25, delta_chi = -10)))[1, 2]
dlm_midlma_lowchi_intercept <- summary(emmeans(delta_live_mass_lm, ~delta_live_mass,
                                              at = list(delta_live_mass = 0, delta_lma = 0, delta_chi = -10)))[1, 2]
dlm_midlma_lowchi_slope <- summary(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass',
                                           at = list(delta_lma = 0, delta_chi = -10)))[1, 2]
dlm_highlma_lowchi_intercept <- summary(emmeans(delta_live_mass_lm, ~delta_live_mass,
                                              at = list(delta_live_mass = 0, delta_lma = 25, delta_chi = -10)))[1, 2]
dlm_highlma_lowchi_slope <- summary(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass',
                                           at = list(delta_lma = 25, delta_chi = -10)))[1, 2]

dlm_lowlma_lowchi_trend <- dlm_lowlma_lowchi_slope *
  seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
  dlm_lowlma_lowchi_intercept

dlm_midlma_lowchi_trend <- dlm_midlma_lowchi_slope *
  seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
  dlm_midlma_lowchi_intercept

dlm_highlma_lowchi_trend <- dlm_highlma_lowchi_slope *
  seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
  dlm_highlma_lowchi_intercept

dlm_lowchi_trend_df <- data.frame(seq(min(delta_live_mass_plot_data$delta_live_mass),
                                               max(delta_live_mass_plot_data$delta_live_mass), 1),
                                dlm_lowlma_lowchi_trend,
                                dlm_midlma_lowchi_trend,
                                dlm_highlma_lowchi_trend)
colnames(dlm_lowchi_trend_df) <- c('delta_live_mass', 'delta_narea_lowlma', 'delta_narea_midlma', 'delta_narea_highlma')

# mid chi
dlm_lowlma_midchi_intercept <- summary(emmeans(delta_live_mass_lm, ~delta_live_mass,
                                               at = list(delta_live_mass = 0, delta_lma = -25, delta_chi = 0)))[1, 2]
dlm_lowlma_midchi_slope <- summary(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass',
                                            at = list(delta_lma = -25, delta_chi = 0)))[1, 2]
dlm_midlma_midchi_intercept <- summary(emmeans(delta_live_mass_lm, ~delta_live_mass,
                                               at = list(delta_live_mass = 0, delta_lma = 0, delta_chi = 0)))[1, 2]
dlm_midlma_midchi_slope <- summary(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass',
                                            at = list(delta_lma = 0, delta_chi = 0)))[1, 2]
dlm_highlma_midchi_intercept <- summary(emmeans(delta_live_mass_lm, ~delta_live_mass,
                                                at = list(delta_live_mass = 0, delta_lma = 25, delta_chi = 0)))[1, 2]
dlm_highlma_midchi_slope <- summary(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass',
                                             at = list(delta_lma = 25, delta_chi = 0)))[1, 2]

dlm_lowlma_midchi_trend <- dlm_lowlma_midchi_slope *
  seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
  dlm_lowlma_midchi_intercept

dlm_midlma_midchi_trend <- dlm_midlma_midchi_slope *
  seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
  dlm_midlma_midchi_intercept

dlm_highlma_midchi_trend <- dlm_highlma_midchi_slope *
  seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
  dlm_highlma_midchi_intercept

dlm_midchi_trend_df <- data.frame(seq(min(delta_live_mass_plot_data$delta_live_mass),
                                    max(delta_live_mass_plot_data$delta_live_mass), 1),
                                dlm_lowlma_midchi_trend,
                                dlm_midlma_midchi_trend,
                                dlm_highlma_midchi_trend)
colnames(dlm_midchi_trend_df) <- c('delta_live_mass', 'delta_narea_lowlma', 'delta_narea_midlma', 'delta_narea_highlma')

# high chi
dlm_lowlma_highchi_intercept <- summary(emmeans(delta_live_mass_lm, ~delta_live_mass,
                                               at = list(delta_live_mass = 0, delta_lma = -25, delta_chi = 10)))[1, 2]
dlm_lowlma_highchi_slope <- summary(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass',
                                            at = list(delta_lma = -25, delta_chi = 10)))[1, 2]
dlm_midlma_highchi_intercept <- summary(emmeans(delta_live_mass_lm, ~delta_live_mass,
                                               at = list(delta_live_mass = 0, delta_lma = 0, delta_chi = 10)))[1, 2]
dlm_midlma_highchi_slope <- summary(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass',
                                            at = list(delta_lma = 0, delta_chi = 10)))[1, 2]
dlm_highlma_highchi_intercept <- summary(emmeans(delta_live_mass_lm, ~delta_live_mass,
                                                at = list(delta_live_mass = 0, delta_lma = 25, delta_chi = 10)))[1, 2]
dlm_highlma_highchi_slope <- summary(emtrends(delta_live_mass_lm, ~delta_live_mass, var = 'delta_live_mass',
                                             at = list(delta_lma = 25, delta_chi = 10)))[1, 2]

dlm_lowlma_highchi_trend <- dlm_lowlma_highchi_slope *
  seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
  dlm_lowlma_highchi_intercept

dlm_midlma_highchi_trend <- dlm_midlma_highchi_slope *
  seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
  dlm_midlma_highchi_intercept

dlm_highlma_highchi_trend <- dlm_highlma_highchi_slope *
  seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
  dlm_highlma_highchi_intercept

dlm_highchi_trend_df <- data.frame(seq(min(delta_live_mass_plot_data$delta_live_mass),
                                      max(delta_live_mass_plot_data$delta_live_mass), 1),
                                  dlm_lowlma_highchi_trend,
                                  dlm_midlma_highchi_trend,
                                  dlm_highlma_highchi_trend)
colnames(dlm_highchi_trend_df) <- c('delta_live_mass', 'delta_narea_lowlma', 'delta_narea_midlma', 'delta_narea_highlma')

(dlm_lowchi_plot <- ggplot(data = delta_live_mass_plot_data, 
                                aes(x = delta_live_mass, y = delta_narea, fill = delta_lma)) +
    theme(legend.position = c(1, 1),
          legend.justification = c(1, 1),
          legend.title = element_text(size = 28),
          legend.text = element_text(size = 20),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          plot.title = element_text(size = 30, colour = 'black', hjust = 0.5),
          axis.title.y = element_text(size = 30, colour = 'black'),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 20, colour = 'black'),
          axis.text.y = element_text(size = 20, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "grey")) +
    # geom_point(shape = 21, colour = 'black', stroke = 0.5, size = 3) +
    # scale_size_continuous(range = c(1, 5)) +
    # scale_fill_gradient(low = 'grey80', high = 'grey0') +
    geom_line(data = dlm_lowchi_trend_df,
              aes(x = delta_live_mass, y = delta_narea_lowlma, fill = NULL),
              size = 5, colour = 'palegreen', alpha = 1, lty = 5) +
    geom_line(data = dlm_lowchi_trend_df,
              aes(x = delta_live_mass, y = delta_narea_midlma, fill = NULL),
              size = 5, colour = 'palegreen3', alpha = 1, lty = 3) +
    geom_line(data = dlm_lowchi_trend_df,
              aes(x = delta_live_mass, y = delta_narea_highlma, fill = NULL),
              size = 5, colour = 'palegreen4', alpha = 1, lty = 3) +
    scale_y_continuous(name = expression('∆' * italic('N')['area'] * ' (%)'), 
                       limits = c(-60, 120),
                       breaks = c(-60, 0, 60, 120)) +
    scale_x_continuous(name = expression('∆' *'AGB' * ' (%)')) +
    labs(title = expression("Low ∆"*chi)))

(dlm_midchi_plot <- ggplot(data = delta_live_mass_plot_data, 
                           aes(x = delta_live_mass, y = delta_narea)) +
    theme(legend.position = c(1, 1),
          legend.justification = c(1, 1),
          legend.title = element_text(size = 28),
          legend.text = element_text(size = 20),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          plot.title = element_text(size = 30, colour = 'black', hjust = 0.5),
          # axis.title.y = element_text(size = 30, colour = 'black'),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 20, colour = 'black'),
          axis.text.y = element_text(size = 20, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "grey")) +
    # geom_point(shape = 21, colour = 'black', fill = 'grey', stroke = 0.5, size = 3) +
    # scale_size_continuous(range = c(1, 5)) +
    # scale_fill_gradient(low = 'grey80', high = 'grey0') +
    geom_line(data = dlm_midchi_trend_df,
              aes(x = delta_live_mass, y = delta_narea_lowlma, fill = NULL),
              size = 5, colour = 'palegreen', alpha = 1, lty = 3) +
    geom_line(data = dlm_midchi_trend_df,
              aes(x = delta_live_mass, y = delta_narea_midlma, fill = NULL),
              size = 5, colour = 'palegreen3', alpha = 1, lty = 3) +
    geom_line(data = dlm_midchi_trend_df,
              aes(x = delta_live_mass, y = delta_narea_highlma, fill = NULL),
              size = 5, colour = 'palegreen4', alpha = 1, lty = 3) +
    scale_y_continuous(
      # name = expression('∆' * italic('N')['area'] * ' (%)'), 
                       limits = c(-60, 120),
                       breaks = c(-60, 0, 60, 120)) +
    scale_x_continuous(name = expression('∆' *'AGB' * ' (%)')) +
    labs(title = expression("Medium ∆"*chi)))

(dlm_highchi_plot <- ggplot(data = delta_live_mass_plot_data, 
                           aes(x = delta_live_mass, y = delta_narea)) +
    theme(legend.position = c(1, 1),
          legend.justification = c(1, 1),
          legend.title = element_text(size = 28),
          legend.text = element_text(size = 20),
          legend.background = element_rect(fill = 'white', colour = 'black'),
          plot.title = element_text(size = 30, colour = 'black', hjust = 0.5),
          # axis.title.y = element_text(size = 30, colour = 'black'),
          axis.title.y = element_blank(),
          axis.title.x = element_text(size = 30, colour = 'black'),
          axis.text.x = element_text(size = 20, colour = 'black'),
          axis.text.y = element_text(size = 20, colour = 'black'),
          panel.background = element_rect(fill = 'white', colour = 'black'),
          panel.grid.major = element_line(colour = "grey")) +
    # geom_point(shape = 21, colour = 'black', fill = 'grey', stroke = 0.5, size = 3) +
    # scale_size_continuous(range = c(1, 5)) +
    # scale_fill_gradient(low = 'grey80', high = 'grey0') +
    geom_line(data = dlm_highchi_trend_df,
              aes(x = delta_live_mass, y = delta_narea_lowlma, fill = NULL),
              size = 5, colour = 'palegreen', alpha = 1, lty = 5) +
    geom_line(data = dlm_highchi_trend_df,
              aes(x = delta_live_mass, y = delta_narea_midlma, fill = NULL),
              size = 5, colour = 'palegreen3', alpha = 1, lty = 3) +
    geom_line(data = dlm_highchi_trend_df,
              aes(x = delta_live_mass, y = delta_narea_highlma, fill = NULL),
              size = 5, colour = 'palegreen4', alpha = 1, lty = 1) +
    scale_y_continuous(
      #name = expression('∆' * italic('N')['area'] * ' (%)'), 
                       limits = c(-60, 120),
                       breaks = c(-60, 0, 60, 120)) +
    scale_x_continuous(name = expression('∆' *'AGB' * ' (%)')) +
    labs(title = expression("High ∆"*chi)))

(delta_live_mass_plot <- dlm_lowchi_plot + dlm_midchi_plot + dlm_highchi_plot + 
  plot_layout(guides = 'collect') +
  plot_annotation(tag_levels = 'A') & 
  theme(plot.tag = element_text(size = 16)))

#### what determines the biomass responses??
delta_live_mass_response_lm <- lmer(delta_live_mass ~ 
                                      spp_live_mass.x * Ptrt_fac * Ktrt_fac +
                                      tmp.x + p_pet.x +
                             (1|Taxon) + (1|Taxon:site_code) +
                             (1|Taxon:site_code:block_fac), 
                           data = delta_live_mass_data)
# plot(resid(delta_live_mass_response_lm) ~ fitted(delta_live_mass_response_lm))
summary(delta_live_mass_response_lm) # N = 309
Anova(delta_live_mass_response_lm)

plot(delta_live_mass ~ (spp_live_mass.x), data = delta_live_mass_data, xlim = c(0, 400))

#### discussion things
leaf_site_climate <- leaf_chi_subset_npred %>%
  group_by(site_code) %>%
  summarise_at(vars(par, p_pet, vpd, tmp),
               mean, na.rm = TRUE)
cor(leaf_site_climate$par, leaf_site_climate$tmp)^2

#### save plots ####

ggsave("plots/hypothesis_plot.jpeg", plot = hypothesis_plot,
       width = 20, height = 25, units = "cm")

ggsave("plots/tmp_plot.jpeg", plot = tmp_plot, 
       width = 29, height = 18, units = "cm")
ggsave("plots/lma_plot.jpeg", plot = lma_plot, 
       width = 29, height = 18, units = "cm")
ggsave("plots/narea_plot.jpeg", plot = narea_plot, 
       width = 29, height = 18, units = "cm")
ggsave("plots/narea_treemap.jpeg", plot = narea_treemap, 
       width = 29, height = 18, units = "cm")
# ggsave("plots/narea_plot_treemap.jpeg", plot = narea_plot_treemap, 
#        width = 38, height = 18, units = "cm")

ggsave("plots/npred_plot.jpeg", plot = npred_plot, 
       width = 29, height = 18, units = "cm")
ggsave("plots/narea_pred_treemap.jpeg", plot = narea_pred_treemap, 
       width = 29, height = 18, units = "cm")

ggsave("plots/live_mass_plot.jpeg", plot = live_mass_plot, 
       width = 29, height = 18, units = "cm")

ggsave("plots/delta_live_mass_plot.jpeg", plot = delta_live_mass_plot, 
       width = 60, height = 25, units = "cm")

