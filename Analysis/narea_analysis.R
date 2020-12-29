# Script using least-cost to explore Narea response at nutnet
## generally follows the analysis structure of Dong et al. (2017) Biogeosciences
## but includes the addition of soil N as a predictor

#####################################################################
## load packages
#####################################################################
library(tidyverse)
library(lme4)
library(car)
library(r2glmm)
library(treemapify)
library(emmeans)
library(relaimpo)

#####################################################################
## load functions
#####################################################################
# source('optimal_vcmax_R/calc_optimal_vcmax.R')
source('optimal_vcmax_R/calc_optimal_vcmax_knownchi.R')
sourceDirectory('optimal_vcmax_R/functions')
source('C4model/C4model_knownchi.R')
sourceDirectory('C4model/functions')
source('n_from_gas_exchange/n_from_gas_exchange.R')

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

calc_vcmax_tresp_mult = function(tleaf, tmean, tref){
  
  temp = tleaf + 273.15
  Ha= 71513
  Hd= 200000
  adelS= 668.39
  bdelS= -1.07
  tmeanK=tmean+273.15
  trefK=tref+273.15
  R=8.314
  kbeg=exp(Ha*(temp-trefK)/(trefK*R*temp))
  kend=((1+exp((trefK*(adelS+bdelS*tmean)-Hd)/(trefK*R)))/(1+exp((temp*(adelS+bdelS*tmean)-Hd)/(temp*R))))
  kbeg*kend
  
}

calc_jmax_tresp_mult = function(tleaf, tmean, tref){
  
  temp = tleaf + 273.15
  Ha= 49884
  Hd= 200000
  adelS= 659.7
  bdelS= -0.75
  tmeanK=tmean+273.15
  trefK=tref+273.15
  R=8.314
  kbeg=exp(Ha*(temp-trefK)/(trefK*R*temp))
  kend=((1+exp((trefK*(adelS+bdelS*tmean)-Hd)/(trefK*R)))/(1+exp((temp*(adelS+bdelS*tmean)-Hd)/(temp*R))))
  kbeg*kend
  
}

# relative importance for mixed models from https://gist.github.com/BERENZ/e9b581a4b7160357934e
calc.relip.mm <- function(model,type='lmg') {
  if (!isLMM(model) & !isGLMM(model)) {
    stop('Currently supports only lmer/glmer objects', call. = FALSE)
  }
  require(lme4)
  X <- getME(model,'X')
  X <- X[,-1]
  Y <- getME(model,'y')
  s_resid <- sigma(model)
  s_effect <- getME(model,'theta')*s_resid
  s2 <- sum(s_resid^2,s_effect^2)
  V <- Diagonal(x = s2,n=nrow(X))
  YX <- cbind(Y,X)
  cov_XY <- solve( t(YX) %*% solve(V) %*% as.matrix(YX))
  colnames(cov_XY) <- rownames(cov_XY) <- colnames(YX)
  importances <- calc.relimp(as.matrix(cov_XY),rela=F,type=type)
  return(importances)
}

#####################################################################
## create hypothesis figure
### leaf N and chi by N supply with different differences in N demand
#####################################################################
n_supply_trend = calc_optimal_vcmax(beta = seq(100, 300, 10))
n_supply_trend$photo_n = fvcmax25_nrubisco(n_supply_trend$vcmax) + fjmax25_nbioe(n_supply_trend$jmax)
n_supply_trend$photo_n_nobetachange = n_supply_trend$photo_n[10]

hypothesis_data = data.frame(cbind(c(n_supply_trend$beta, n_supply_trend$beta), 
                        c(n_supply_trend$photo_n, n_supply_trend$photo_n_nobetachange),
                        c(rep('no_change', 21), rep('change', 21))))
colnames(hypothesis_data) = c('beta', 'Narea', 'demand')
hypothesis_data$beta = as.numeric(as.character(hypothesis_data$beta))
hypothesis_data$Narea = as.numeric(as.character(hypothesis_data$Narea))

hypothesis_plot = ggplot(data = hypothesis_data, 
                         aes(x = (1/beta), y = Narea, col = demand)) +
  theme(legend.position = c(0.3, 0.9),
        legend.text = element_text(size = 25),
        axis.title.y=element_text(size=rel(4), colour = 'black'),
        axis.title.x=element_text(size=rel(4), colour = 'black'),
        axis.text.x=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.x=element_blank(),
        axis.ticks.y=element_blank(),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey"),
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black")) +
  geom_line(size = 4, aes(linetype = demand, color = demand)) +
  scale_linetype_manual(values = c(1, 2), guide = F) +
  scale_colour_manual(values = c('black', 'grey'), 
                      labels = c('∆N demand = ∆N supply', '∆N demand = 0')) +
  guides(color = guide_legend(title = NULL)) +
  ylab(expression('Leaf ' * italic('N')['area'])) +
  xlab('N supply')

jpeg(filename = "plots/hypothesis_plot.jpeg", width = 600, height = 700, units = 'px')
plot(hypothesis_plot)
dev.off()

#####################################################################
## soil N effects on leaf traits (Narea and chi)
#####################################################################
### load data
leaf = read.csv('../Data/processed/traits.csv')
leaf$Ntrt_fac = as.factor(leaf$Ntrt)
leaf$Ptrt_fac = as.factor(leaf$Ptrt)
leaf$Ktrt_fac = as.factor(leaf$Ktrt)
leaf$block_fac = as.factor(leaf$block)
leaf$tmp_scaled = leaf$tmp - 25
leaf$pft <- leaf$functional_group
leaf$pft[leaf$pft == 'NULL'] <- NA
leaf$grass <- 'no'
leaf$grass[leaf$pft == 'GRASS'] <- 'yes'
leaf$fence <- 'no'
leaf$fence[leaf$trt == 'Fence' | leaf$trt == 'NPK+Fence'] <- 'yes'
# leaf_Nonly = subset(leaf, Ptrt_fac == '0' & Ktrt_fac == '0')
leaf_chi_subset = subset(leaf, chi > 0.2 & chi < 0.95) # lose 700 points

# hist(leaf$chi)
# hist(log(leaf$Narea))
leafchi_lmer = lmer(logit(chi) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac + tmp + log(vpd) + z +
                      Nfix + photosynthetic_pathway +
                      (1|Taxon) + (1|Taxon:site_code) + (1|Taxon:site_code:block_fac), 
                    data = leaf_chi_subset)
# plot(resid(leafchi_lmer) ~ fitted(leafchi_lmer))
summary(leafchi_lmer)
Anova(leafchi_lmer)
cld(emmeans(leafchi_lmer, ~Ntrt_fac))
(summary(emmeans(leafchi_lmer, ~Ntrt_fac))[2,2] - summary(emmeans(leafchi_lmer, ~Ntrt_fac))[1,2]) / 
  abs(summary(emmeans(leafchi_lmer, ~Ntrt_fac))[1,2])

leaflma_lmer = lmer(log(lma) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac + tmp + log(vpd) + z +
                      Nfix + photosynthetic_pathway +
                      (1|Taxon) + (1|Taxon:site_code) + (1|Taxon:site_code:block_fac), 
                    data = leaf_chi_subset)
# plot(resid(leaflma_lmer) ~ fitted(leaflma_lmer))
summary(leaflma_lmer)
Anova(leaflma_lmer)
cld(emmeans(leaflma_lmer, ~Ntrt_fac))
(summary(emmeans(leaflma_lmer, ~Ntrt_fac))[2,2] - summary(emmeans(leaflma_lmer, ~Ntrt_fac))[1,2]) / 
  abs(summary(emmeans(leaflma_lmer, ~Ntrt_fac))[1,2])

### run analyses
# hist(leaf$narea)
# hist(log(leaf$narea))
leafNarea_lmer = lmer(log(narea) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac + chi + 
                        tmp + log(par) +  log(vpd) + z +
                        log(lma) + Nfix + photosynthetic_pathway +
                        (1|Taxon) + (1|Taxon:site_code) + (1|Taxon:site_code:block_fac), 
                      data = leaf_chi_subset)
# plot(resid(leafNarea_lmer) ~ fitted(leafNarea_lmer))
Anova(leafNarea_lmer)
summary(leafNarea_lmer)

cld(emmeans(leafNarea_lmer, ~Ntrt_fac))
(exp(summary(emmeans(leafNarea_lmer, ~Ntrt_fac))[2,2]) - exp(summary(emmeans(leafNarea_lmer, ~Ntrt_fac))[1,2])) / 
  exp(abs(summary(emmeans(leafNarea_lmer, ~Ntrt_fac))[1,2]))

cld(emmeans(leafNarea_lmer, ~Ntrt_fac * Ptrt_fac))

cld(emmeans(leafNarea_lmer, ~Nfix))
cld(emmeans(leafNarea_lmer, ~photosynthetic_pathway))
test(emtrends(leafNarea_lmer, ~1, var = 'chi'))
test(emtrends(leafNarea_lmer, ~1, var = 'log(lma)'))
test(emtrends(leafNarea_lmer, ~1, var = 'tmp'))
test(emtrends(leafNarea_lmer, ~1, var = 'log(par)'))

calc.relip.mm(leafNarea_lmer)

### make figure
narea_plot = ggplot(data = leaf_chi_subset, 
                         aes(x = Ntrt_fac, y = log(narea))) +
  theme(legend.position = "none", 
        axis.title.y=element_text(size=rel(4), colour = 'black'),
        axis.title.x=element_text(size=rel(4), colour = 'black'),
        axis.text.x=element_text(size=rel(3), colour = 'black'),
        axis.text.y=element_text(size=rel(3), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey"),
        plot.tag = element_text(size = 30)) +
  geom_boxplot(outlier.color = NA, fill = 'white', lwd = 5) +
  geom_dotplot(binaxis = 'y', binwidth = 0.1, stackdir = 'center', 
               fill = 'burlywood1', alpha = 0.8) +
  scale_x_discrete(labels = c('Ambient', 'Added N')) +
  xlab('') +
  ylab(expression('Log leaf ' * italic('N')['area'])) +
  xlab('Soil N treatment') +
  labs(tag = 'A')

# jpeg(filename = "plots/narea_plot.jpeg", width = 800, height = 900, units = 'px')
# plot(narea_plot)
# dev.off()

# chi_plot = ggplot(data = leaf_Nonly, 
#                   aes(x = Ntrt_fac, y = chi)) +
#   theme(legend.position = "none", 
#         axis.title.y=element_text(size=rel(2.5), colour = 'black'),
#         axis.title.x=element_text(size=rel(2.5), colour = 'black'),
#         axis.text.x=element_text(size=rel(2), colour = 'black'),
#         axis.text.y=element_text(size=rel(2), colour = 'black'),
#         panel.background = element_rect(fill = 'white', colour = 'black'),
#         panel.grid.major = element_line(colour = "grey")) +
#   geom_boxplot(outlier.color = NA, fill = 'white') +
#   geom_dotplot(binaxis = 'y', binwidth = 0.005, stackdir = 'center', alpha = 0.5) +
#   scale_x_discrete(labels = c('Ambient', 'Added N')) +
#   xlab('Nitrogen treatment') +
#   ylab(expression('C'[i] * '/C'[a])) +
#   annotate("text", x = 1.5, y = 0.93, label = "p = 0.05", size = 8)

relimp_leafn = calc.relip.mm(leafNarea_lmer)$lmg
relimp_leafn_df = NULL
relimp_leafn_df$Factor = c('Soil N', 'Soil P', 'Soil K+µ', 'χ', 'Temperature', 'PAR', 'LMA', 'N fixer', 'C3/C4', 'Soil interactions')
relimp_leafn_df$Importance = as.numeric(as.character(c(relimp_leafn[1:9], sum(relimp_leafn[10:13]))))
sum(relimp_leafn[c(1:3, 10:13)]) # importance of soil
relimp_leafn_df = as.data.frame(relimp_leafn_df)
unexplained = 1 - sum(relimp_leafn_df$Importance)
unexplained_df = data.frame("Unexplained", unexplained)
names(unexplained_df) = c("Factor", "Importance")
relimp_leafn_df = rbind(relimp_leafn_df, unexplained_df)

narea_treemap = ggplot(relimp_leafn_df, aes(area = Importance, label = Factor)) +
  theme(plot.tag = element_text(size = 30)) +
  geom_treemap(fill = c('burlywood1', 'burlywood2', 'burlywood3', 
                        'blue', 'yellow', 'red', 'grey', 'orange', 'green', 
                        'burlywood4', 'white'), colour = 'black') +
  geom_treemap_text(colour = "black", place = "centre",
                    grow = TRUE) +
    labs(tag = 'B')

# jpeg(filename = "plots/narea_treemap.jpeg", width = 900, height = 900, units = 'px')
# plot(narea_treemap)
# dev.off()

jpeg(filename = "plots/narea_plot_treemap.jpeg", width = 1700, height = 900, units = 'px')
multiplot(narea_plot, narea_treemap, cols = 2)
dev.off()

#####################################################################
## Narea predictions
#####################################################################
gas_exchange_pred_c3 = calc_optimal_vcmax_knownchi(tg_c = leaf$tmp, 
                                       paro = leaf$par, 
                                       cao = 400, 
                                       vpdo = leaf$vpd, 
                                       z = leaf$z,
                                       chi = leaf$chi,
                                       q0 = 0.1)

gas_exchange_pred_c4 = C4model_knownchi(tg_c = leaf$tmp, 
                                                   paro = leaf$par, 
                                                   cao = 400, 
                                                   vpdo = leaf$vpd, 
                                                   z = leaf$z,
                                                   chi = leaf$chi,
                                        q0 = 0.1)

leaf$vcmax25[leaf$photosynthetic_pathway == 'C3'] = gas_exchange_pred_c3$vcmax[leaf$photosynthetic_pathway == 'C3'] /
  calc_vcmax_tresp_mult(gas_exchange_pred_c3$tg_c[leaf$photosynthetic_pathway == 'C3'], 
                        gas_exchange_pred_c3$tg_c[leaf$photosynthetic_pathway == 'C3'], 25)
leaf$jmax25[leaf$photosynthetic_pathway == 'C3'] = gas_exchange_pred_c3$jmax[leaf$photosynthetic_pathway == 'C3'] /
  calc_jmax_tresp_mult(gas_exchange_pred_c3$tg_c[leaf$photosynthetic_pathway == 'C3'], 
                       gas_exchange_pred_c3$tg_c[leaf$photosynthetic_pathway == 'C3'], 25)
leaf$vpmax25[leaf$photosynthetic_pathway == 'C3'] = 0
leaf$vcmax25[leaf$photosynthetic_pathway == 'C4'] = gas_exchange_pred_c4$vcmax[leaf$photosynthetic_pathway == 'C4'] /
  calc_vcmax_tresp_mult(gas_exchange_pred_c4$tg_c[leaf$photosynthetic_pathway == 'C4'], 
                        gas_exchange_pred_c4$tg_c[leaf$photosynthetic_pathway == 'C4'], 25)
leaf$jmax25[leaf$photosynthetic_pathway == 'C4'] = gas_exchange_pred_c4$jmax[leaf$photosynthetic_pathway == 'C4'] /
  calc_jmax_tresp_mult(gas_exchange_pred_c4$tg_c[leaf$photosynthetic_pathway == 'C4'], 
                       gas_exchange_pred_c4$tg_c[leaf$photosynthetic_pathway == 'C4'], 25)
leaf$vpmax25[leaf$photosynthetic_pathway == 'C4'] = gas_exchange_pred_c4$vpmax[leaf$photosynthetic_pathway == 'C4'] /
  calc_vcmax_tresp_mult(gas_exchange_pred_c4$tg_c[leaf$photosynthetic_pathway == 'C4'], 
                        gas_exchange_pred_c4$tg_c[leaf$photosynthetic_pathway == 'C4'], 25)

leaf$nrubisco = fvcmax25_nrubisco(leaf$vcmax25)
leaf$nbioe = fjmax25_nbioe(leaf$jmax25)
leaf$npep = fvpmax25_npep(leaf$vpmax25)
leaf$nstructure = flma_nstructure(leaf$lma)
leaf$nall = leaf$nrubisco + leaf$nbioe + leaf$nstructure + leaf$npep
leaf$nphoto = leaf$nrubisco + leaf$nbioe + leaf$npep
leaf$nrubisco_frac = leaf$nrubisco / leaf$nall
leaf$nphoto_frac = leaf$nphoto / leaf$nall

# npred_lmer = lmer(log(narea) ~ nphoto + nstructure +
#                     (1|site_code) + (1|site_code:block_fac),
#                   data = leaf)
# Anova(npred_lmer)
# summary(npred_lmer)
# calc.relip.mm(npred_lmer)

leaf$lognphoto = log(leaf$nphoto)
leaf$lognstructure = log(leaf$nstructure)
npred_soil_lmer = lmer(log(narea) ~ lognphoto + lognstructure + 
                         Ntrt_fac * Ptrt_fac * Ktrt_fac +
                         Nfix + photosynthetic_pathway +
                         (1|Taxon) + (1|Taxon:site_code) + (1|Taxon:site_code:block_fac),
                  data = leaf_chi_subset)
plot(resid(npred_soil_lmer) ~ fitted(npred_soil_lmer))
Anova(npred_soil_lmer)
summary(npred_soil_lmer)
calc.relip.mm(npred_soil_lmer)

nphoto_slope = summary(emtrends(npred_soil_lmer, ~1, var = "lognphoto"))[1, 2]
nphoto_intercept_lowN = summary(emmeans(npred_soil_lmer, ~Ntrt_fac, at = list(lognphoto = 0)))[1, 2]
nphoto_intercept_highN = summary(emmeans(npred_soil_lmer, ~Ntrt_fac, at = list(lognphoto = 0)))[2, 2]
lognphoto_seq = seq(min(leaf$lognphoto), max(leaf$lognphoto), 0.01)
nphoto_trend_lowN = nphoto_intercept_lowN + lognphoto_seq * nphoto_slope
nphoto_trend_highN = nphoto_intercept_highN + lognphoto_seq * nphoto_slope
nphoto_trend = as.data.frame(cbind(lognphoto_seq, nphoto_trend_lowN, nphoto_trend_highN))

npred_photo_plot = ggplot(data = leaf_chi_subset, 
                    aes(x = lognphoto, y = log(narea), color = Ntrt_fac)) +
  theme(legend.position = "none", 
        axis.title.y=element_text(size=rel(2.5), colour = 'black'),
        axis.title.x=element_text(size=rel(2.5), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_point(shape = 16, size = 3, alpha = 0.5) +
  scale_color_manual(values = c("black", "burlywood1"), labels = c("Ambient", "Added N")) +
  labs(color = 'Soil N') +
  geom_line(data = nphoto_trend, aes(x = lognphoto_seq, y = nphoto_trend_lowN), 
            col = 'black', lwd = 3, alpha = 0.8) +
  geom_line(data = nphoto_trend, aes(x = lognphoto_seq, y = nphoto_trend_highN), 
            col = 'burlywood1', lwd = 3, alpha = 0.8) +
  ylab(expression('Log leaf ' * italic('N')['area'])) +
  xlab(expression('Log leaf ' * italic('N')['photo']))

nstructure_slope = summary(emtrends(npred_soil_lmer, ~1, var = "lognstructure"))[1, 2]
nstructure_intercept_lowN = summary(emmeans(npred_soil_lmer, ~Ntrt_fac, at = list(lognstructure = 0)))[1, 2]
nstructure_intercept_highN = summary(emmeans(npred_soil_lmer, ~Ntrt_fac, at = list(lognstructure = 0)))[2, 2]
lognstructure_seq = seq(min(leaf$lognstructure, na.rm = T), max(leaf$lognstructure, na.rm = T), 0.1)
nstructure_trend_lowN = nstructure_intercept_lowN + lognstructure_seq * nstructure_slope
nstructure_trend_highN = nstructure_intercept_highN + lognstructure_seq * nstructure_slope
nstructure_trend = as.data.frame(cbind(lognstructure_seq, nstructure_trend_lowN, nstructure_trend_highN))

npred_structure_plot = ggplot(data = leaf_chi_subset, 
                              aes(x = lognstructure, y = log(narea), color = Ntrt_fac)) +
  theme(legend.position = c(0, 1),
        legend.justification = c(0, 1),
        legend.title = element_text(size = 20),
        legend.text = element_text(size = 20),
        legend.background = element_rect(fill = 'white', colour = 'black'),
        axis.title.y=element_text(size=rel(2.5), colour = 'black'),
        axis.title.x=element_text(size=rel(2.5), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_point(shape = 16, size = 3, alpha = 0.5) +
  scale_color_manual(values = c("black", "burlywood1"), labels = c("Ambient", "Added N")) +
  labs(color = 'Soil N') +
  geom_line(data = nstructure_trend, aes(x = lognstructure_seq, y = nstructure_trend_lowN), 
            col = 'black', lwd = 3, alpha = 0.8) +
  geom_line(data = nstructure_trend, aes(x = lognstructure_seq, y = nstructure_trend_highN), 
            col = 'burlywood1', lwd = 3, alpha = 0.8) +
  ylab(expression('Log leaf ' * italic('N')['area'])) +
  xlab(expression('Log leaf ' * italic('N')['structure']))

jpeg(filename = "plots/npred_plot.jpeg", width = 800, height = 500, units = 'px')
multiplot(npred_structure_plot, npred_photo_plot, cols=2)
dev.off()

relimp_leafn_pred = calc.relip.mm(npred_soil_lmer)$lmg
relimp_leafn_pred_df = NULL
relimp_leafn_pred_df$Factor = c('Nphoto', 'Nstructure', 'Soil N', 'Soil P', 'Soil K+µ', 'N fixer', 'C3/C4', 'Soil interactions')
relimp_leafn_pred_df$Importance = as.numeric(as.character(c(relimp_leafn_pred[1:7], sum(relimp_leafn_pred[8:11]))))
sum(relimp_leafn_pred[c(3:5, 8:11)]) # importance of soil
relimp_leafn_pred_df = as.data.frame(relimp_leafn_pred_df)
unexplained = 1 - sum(relimp_leafn_pred_df$Importance)
unexplained_pred_df = data.frame("Unexplained", unexplained)
names(unexplained_pred_df) = c("Factor", "Importance")
relimp_leafn_pred_df = rbind(relimp_leafn_pred_df, unexplained_pred_df)

narea_pred_treemap = ggplot(relimp_leafn_pred_df, aes(area = Importance, label = Factor)) +
  theme(plot.tag = element_text(size = 30)) +
  geom_treemap(fill = c('darkgreen', 'grey', 'burlywood1', 'burlywood2', 'burlywood3', 
                        'orange', 'green', 
                        'burlywood4', 'white'), colour = 'black') +
  geom_treemap_text(colour = "black", place = "centre",
                    grow = TRUE)

jpeg(filename = "plots/narea_pred_treemap.jpeg", width = 850, height = 900, units = 'px')
plot(narea_pred_treemap)
dev.off()

#####################################################################
## soil N effects on LAI
#####################################################################
leaf_site = read.csv('../Data/processed/traits_site.csv')
leaf_site$Ntrt_fac = as.factor(leaf_site$Ntrt)
leaf_site$Ptrt_fac = as.factor(leaf_site$Ptrt)
leaf_site$Ktrt_fac = as.factor(leaf_site$Ktrt)
leaf_site$block_fac = as.factor(leaf_site$block)
# leaf_site_Nonly = subset(leaf_site, Ptrt_fac == '0' & Ktrt_fac == '0')

# hist(leaf_site_Nonly$lai_mean)
# hist(log(leaf_site_Nonly$lai_mean))
lai_lmer = lmer(log(lai_mean) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac +
                  # tmp_mean + log(par_mean) +  log(vpd_mean) + z_mean +
                  (1|site_code) + (1|site_code:block), 
                data = leaf_site)
# plot(resid(lai_lmer) ~ fitted(lai_lmer))
Anova(lai_lmer)
summary(lai_lmer)
cld(emmeans(lai_lmer, ~Ntrt_fac))
cld(emmeans(lai_lmer, ~Ntrt_fac*Ptrt_fac))
(exp(summary(emmeans(lai_lmer, ~Ntrt_fac))[2,2]) - exp(summary(emmeans(lai_lmer, ~Ntrt_fac))[1,2])) / 
  abs(exp(summary(emmeans(lai_lmer, ~Ntrt_fac))[1,2]))

calc.relip.mm(lai_lmer)

# hist(leaf_site_Nonly$live_mass_mean)
# hist(log(leaf_site_Nonly$live_mass_mean))
live_mass_lmer = lmer(log(live_mass_mean) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac +
                        # tmp_mean + log(par_mean) +  log(vpd_mean) + z_mean +
                        (1|site_code) + (1|site_code:block), 
                      data = leaf_site)
# # plot(resid(live_mass_lmer) ~ fitted(live_mass_lmer))
Anova(live_mass_lmer)
cld(emmeans(live_mass_lmer, ~Ntrt_fac))
(exp(summary(emmeans(live_mass_lmer, ~Ntrt_fac))[2,2]) - exp(summary(emmeans(live_mass_lmer, ~Ntrt_fac))[1,2])) / 
 abs(exp(summary(emmeans(live_mass_lmer, ~Ntrt_fac))[1,2]))

leaf_site$NPgroup[leaf_site$Ntrt_fac == '0' & leaf_site$Ptrt_fac == '0'] = '-N, -P'
leaf_site$NPgroup[leaf_site$Ntrt_fac == '1' & leaf_site$Ptrt_fac == '0'] = '+N, -P'
leaf_site$NPgroup[leaf_site$Ntrt_fac == '0' & leaf_site$Ptrt_fac == '1'] = '-N, +P'
leaf_site$NPgroup[leaf_site$Ntrt_fac == '1' & leaf_site$Ptrt_fac == '1'] = '+N, +P'

lai_plot = ggplot(data = leaf_site, 
                  aes(x = NPgroup, y = log(lai_mean))) +
  theme(legend.position = "none", 
        axis.title.y=element_text(size=rel(2.5), colour = 'black'),
        axis.title.x=element_text(size=rel(2.5), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_boxplot(outlier.color = NA, fill = 'white') +
  geom_dotplot(aes(fill = NPgroup), binaxis = 'y', binwidth = 0.1, stackdir = 'center', alpha = 0.5) +
  scale_fill_manual(values = c('grey', 'blue', 'red', 'purple')) +
  # scale_x_discrete(labels = c('Ambient', 'Added N')) +
  xlab('N x P treatment') +
  ylab(expression('Log LAI (m' ^ '2' *' m' ^'-2' * ')'))

live_mass_plot = ggplot(data = leaf_site, 
                  aes(x = NPgroup, y = log(live_mass_mean))) +
  theme(legend.position = "none", 
        axis.title.y=element_text(size=rel(2.5), colour = 'black'),
        axis.title.x=element_text(size=rel(2.5), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_boxplot(outlier.color = NA, fill = 'white') +
  geom_dotplot(aes(fill = NPgroup), binaxis = 'y', binwidth = 0.1, stackdir = 'center', alpha = 0.5) +
  scale_fill_manual(values = c('grey', 'blue', 'red', 'purple')) +
  # scale_x_discrete(labels = c('Ambient', 'Added N')) +
  xlab('N x P treatment') +
  ylab(expression('Log Biomass (g m' ^ '2' * ')'))

jpeg(filename = "plots/plant_plot.jpeg", width = 550, height = 900, units = 'px')
multiplot(lai_plot, live_mass_plot, cols=1)
dev.off()

# jpeg(filename = "plots/lai_plot.jpeg", width = 500, height = 550, units = 'px')
# plot(lai_plot)
# dev.off()

#####################################################################
## supply vs demand effects on leaf N
#####################################################################

### non-delta response
#### note: not really sure what this is telling us!
# narea_supdem_lmer = lmer(log(narea) ~ Ntrt_fac * spp_lai +
#                            (1|Taxon) + (1|Taxon:site_code) + (1|Taxon:site_code:block_fac), 
#                          data = leaf)
# plot(resid(narea_supdem_lmer) ~ fitted(narea_supdem_lmer))
# Anova(narea_supdem_lmer)
# test(emtrends(narea_supdem_lmer, ~Ntrt_fac, var = 'spp_lai'))
# cld(emtrends(narea_supdem_lmer, ~Ntrt_fac, var = 'spp_lai'))

### calculate deltas (per species per treatment type per block per site)
#### questions: should we include block into the grouping? this yields 2692 observations
#### removing block (and plot) yields 1267 obs

leaf_site_N_group_by = group_by(leaf_chi_subset, 
                                site_code, Ntrt_fac, Ptrt_fac, Ktrt_fac,
                                block_fac,
                                fence, trt, plot,
                                Taxon, grass, Nfix, photosynthetic_pathway)
leaf_site_N = summarise(leaf_site_N_group_by,
                        n = n(),
                        narea_mean = mean(narea, na.rm = T),
                        spp_lai_mean = mean(spp_lai, na.rm = T),
                        chi_mean = mean(chi, na.rm = T),
                        spp_live_mass_mean = mean(spp_live_mass, na.rm = T),
                        spp_N_mass_mean = mean(spp_mass_N, na.rm = T),
                        max_cover_mean = mean(max_cover, na.rm = T),
                        lma_mean = mean(lma, na.rm = T),
                        p_pet_mean = mean(p_pet, na.rm = T),
                        vpd_mean = mean(vpd, na.rm = T),
                        tmp_mean = mean(tmp, na.rm = T))
# leaf_site_N$PKtrt_fac = as.factor(paste(leaf_site_N$Ptrt_fac, leaf_site_N$Ktrt_fac, sep = ''))
leaf_site_lowN = subset(leaf_site_N, Ntrt_fac == '0')
leaf_site_highN = subset(leaf_site_N, Ntrt_fac == '1')
nrow(leaf_site_lowN)
nrow(leaf_site_highN)
leaf_site_trt = left_join(leaf_site_lowN, leaf_site_highN, 
                                by = c('site_code', 
                                       'block_fac', 
                                       'Ptrt_fac', 'Ktrt_fac',
                                       'fence',
                                       'Taxon', 'grass', 'Nfix', 'photosynthetic_pathway'))
nrow(leaf_site_trt)

leaf_site_trt$delta_narea = ((leaf_site_trt$narea_mean.y - 
                                leaf_site_trt$narea_mean.x) / leaf_site_trt$narea_mean.x) * 100
leaf_site_trt$delta_lai = ((leaf_site_trt$spp_lai_mean.y - 
                              leaf_site_trt$spp_lai_mean.x) / leaf_site_trt$spp_lai_mean.x) * 100
leaf_site_trt$delta_live_mass = ((leaf_site_trt$spp_live_mass_mean.y - 
                                    leaf_site_trt$spp_live_mass_mean.x) / leaf_site_trt$spp_live_mass_mean.x) * 100
leaf_site_trt$delta_N_mass = ((leaf_site_trt$spp_N_mass_mean.y - 
                                    leaf_site_trt$spp_N_mass_mean.x) / leaf_site_trt$spp_N_mass_mean.x) * 100
leaf_site_trt$delta_chi = ((leaf_site_trt$chi_mean.y - 
                              leaf_site_trt$chi_mean.x) / leaf_site_trt$chi_mean.x) * 100
leaf_site_trt$delta_lma = ((leaf_site_trt$lma_mean.y - 
                              leaf_site_trt$lma_mean.x) / leaf_site_trt$lma_mean.x) * 100

## check variables
### lots of outliers!
### remove using mean absolute deviation (Leys et al., 2013); nix anything over 2.5x times the MAD
### mad calculations
delta_narea_mad = mad(leaf_site_trt$delta_narea, na.rm = T)
delta_lai_mad = mad(leaf_site_trt$delta_lai, na.rm = T)
delta_lma_mad = mad(leaf_site_trt$delta_lma, na.rm = T)
delta_live_mass_mad = mad(leaf_site_trt$delta_live_mass, na.rm = T)
delta_N_mass_mad = mad(leaf_site_trt$delta_N_mass, na.rm = T)
delta_chi_mad = mad(leaf_site_trt$delta_chi, na.rm = T)


hist(leaf_site_trt$delta_narea)
sort(leaf_site_trt$delta_narea) # 4 measurements > 1000% change, 2 > 3000%
hist(subset(leaf_site_trt, delta_narea < 3 * delta_narea_mad & delta_narea > 3 * -delta_narea_mad)$delta_narea)

hist(leaf_site_trt$delta_lai)
sort(leaf_site_trt$delta_lai) # 13 measurements > 1000%, 4 measurements > 4000%
hist(subset(leaf_site_trt, delta_lai < 3 * delta_lai_mad & delta_lai > 3 * -delta_lai_mad)$delta_lai)

hist(leaf_site_trt$delta_live_mass)
sort(leaf_site_trt$delta_live_mass)
hist(subset(leaf_site_trt, delta_live_mass < 3 * delta_live_mass_mad & 
              delta_live_mass > 3 * -delta_live_mass_mad)$delta_live_mass)
hist(log(subset(leaf_site_trt, delta_live_mass < 3 * delta_live_mass_mad & 
              delta_live_mass > 3 * -delta_live_mass_mad)$delta_live_mass + 100))

hist(leaf_site_trt$delta_N_mass)
sort(leaf_site_trt$delta_N_mass)
hist(subset(leaf_site_trt, delta_N_mass < 3 * delta_N_mass_mad & 
              delta_N_mass > 3 * -delta_N_mass_mad)$delta_N_mass)

hist(leaf_site_trt$delta_lma)
sort(leaf_site_trt$delta_lma)
hist(subset(leaf_site_trt, delta_lma < 3 * delta_lma_mad &
              delta_lma > 3 * -delta_lma_mad)$delta_lma)

hist(leaf_site_trt$delta_chi)
sort(leaf_site_trt$delta_chi)
hist(subset(leaf_site_trt, delta_chi < 3 * delta_chi_mad & delta_chi > 3 * -delta_chi_mad)$delta_chi)


### model response
delta_lai_data = subset(leaf_site_trt, 
                        delta_narea < 3 * delta_narea_mad & 
                          delta_narea > 3 * -delta_narea_mad & 
                          delta_lma < 3 * delta_lma_mad &
                          delta_lma > 3 * -delta_lma_mad &
                          delta_lai < 3 * delta_lai_mad & 
                          delta_lai > 3 * -delta_lai_mad &
                          delta_chi < 3 * delta_chi_mad & 
                          delta_chi > 3 * -delta_chi_mad)
sort(delta_lai_data$delta_narea)

delta_lai_lm = lmer(delta_narea ~ delta_lai + 
                      Ptrt_fac * Ktrt_fac +
                      photosynthetic_pathway + Nfix +
                      delta_lai * delta_lma * delta_chi +
                      (1|Taxon) + (1|Taxon:site_code) +
                      (1|Taxon:site_code:block_fac), 
                  data = delta_lai_data)
plot(resid(delta_lai_lm) ~ fitted(delta_lai_lm))
Anova(delta_lai_lm)
summary(delta_lai_lm)
test(emtrends(delta_lai_lm, ~1, var = 'delta_lma'))
test(emtrends(delta_lai_lm, ~1, var = 'delta_lai'))
test(emtrends(delta_lai_lm, ~1, var = 'delta_lai', at = list(delta_lma = -25, delta_chi = -15)))
test(emtrends(delta_lai_lm, ~1, var = 'delta_lai', at = list(delta_lma = 0, delta_chi = -15)))
test(emtrends(delta_lai_lm, ~1, var = 'delta_lai', at = list(delta_lma = 25, delta_chi = -15)))
test(emtrends(delta_lai_lm, ~1, var = 'delta_lai', at = list(delta_lma = -25, delta_chi = 0)))
test(emtrends(delta_lai_lm, ~1, var = 'delta_lai', at = list(delta_lma = 0, delta_chi = 0)))
test(emtrends(delta_lai_lm, ~1, var = 'delta_lai', at = list(delta_lma = 25, delta_chi = 0)))
test(emtrends(delta_lai_lm, ~1, var = 'delta_lai', at = list(delta_lma = -25, delta_chi = 15)))
test(emtrends(delta_lai_lm, ~1, var = 'delta_lai', at = list(delta_lma = 0, delta_chi = 15)))
test(emtrends(delta_lai_lm, ~1, var = 'delta_lai', at = list(delta_lma = 25, delta_chi = 15))) ## only one!

delta_live_mass_data = subset(leaf_site_trt, 
                              delta_narea < 3 * delta_narea_mad & 
                                delta_narea > 3 * -delta_narea_mad &
                                delta_lma < 3 * delta_lma_mad &
                                delta_lma > 3 * -delta_lma_mad &
                                delta_live_mass < 3 * delta_live_mass_mad & 
                                delta_live_mass > 3 * -delta_live_mass_mad &
                                delta_chi < 3 * delta_chi_mad &
                                delta_chi > 3 * -delta_chi_mad)

delta_live_mass_lm = lmer(delta_narea ~ delta_live_mass + 
                            Ptrt_fac * Ktrt_fac +
                            photosynthetic_pathway + Nfix +
                            delta_live_mass * delta_lma * delta_chi +
                            (1|Taxon) + (1|Taxon:site_code) +
                            (1|Taxon:site_code:block_fac), 
                        data = delta_live_mass_data)
plot(resid(delta_live_mass_lm) ~ fitted(delta_live_mass_lm))
Anova(delta_live_mass_lm)
summary(delta_live_mass_lm)
test(emtrends(delta_live_mass_lm, ~1, var = 'delta_lma'))
test(emtrends(delta_live_mass_lm, ~1, var = 'delta_live_mass'))
test(emtrends(delta_live_mass_lm, ~1, var = 'delta_live_mass', at = list(delta_lma = -25)))
test(emtrends(delta_live_mass_lm, ~1, var = 'delta_live_mass', at = list(delta_lma = 0))) # marginal
test(emtrends(delta_live_mass_lm, ~1, var = 'delta_live_mass', at = list(delta_lma = 25))) # ding!
# test(emtrends(delta_live_mass_lm, ~1, var = 'delta_live_mass', at = list(delta_lma = -25, delta_chi = -15)))
# test(emtrends(delta_live_mass_lm, ~1, var = 'delta_live_mass', at = list(delta_lma = 0, delta_chi = -15)))
# test(emtrends(delta_live_mass_lm, ~1, var = 'delta_live_mass', at = list(delta_lma = 25, delta_chi = -15)))
# test(emtrends(delta_live_mass_lm, ~1, var = 'delta_live_mass', at = list(delta_lma = -25, delta_chi = 0)))
# test(emtrends(delta_live_mass_lm, ~1, var = 'delta_live_mass', at = list(delta_lma = 0, delta_chi = 0)))
# test(emtrends(delta_live_mass_lm, ~1, var = 'delta_live_mass', at = list(delta_lma = 25, delta_chi = 0))) # ding!
# test(emtrends(delta_live_mass_lm, ~1, var = 'delta_live_mass', at = list(delta_lma = -25, delta_chi = 15)))
# test(emtrends(delta_live_mass_lm, ~1, var = 'delta_live_mass', at = list(delta_lma = 0, delta_chi = 15)))
# test(emtrends(delta_live_mass_lm, ~1, var = 'delta_live_mass', at = list(delta_lma = 25, delta_chi = 15)))

delta_N_mass_data = subset(leaf_site_trt, 
                              delta_narea < 3 * delta_narea_mad & 
                                delta_narea > 3 * -delta_narea_mad &
                                delta_lma < 3 * delta_lma_mad &
                                delta_lma > 3 * -delta_lma_mad &
                                delta_N_mass < 3 * delta_N_mass_mad & 
                                delta_N_mass > 3 * -delta_N_mass_mad &
                                delta_chi < 3 * delta_chi_mad &
                                delta_chi > 3 * -delta_chi_mad)

delta_N_mass_lm = lmer(delta_narea ~ delta_N_mass + 
                         Ptrt_fac * Ktrt_fac +
                         Nfix + # photosynthetic_pathway + Nfix + # no C4
                         delta_N_mass * delta_lma * delta_chi +
                         (1|Taxon) + (1|Taxon:site_code) +
                         (1|Taxon:site_code:block_fac), 
                          data = delta_N_mass_data)
plot(resid(delta_N_mass_lm) ~ fitted(delta_N_mass_lm))
Anova(delta_N_mass_lm)
summary(delta_N_mass_lm)
test(emtrends(delta_N_mass_lm, ~1, var = 'delta_lma'))
test(emtrends(delta_N_mass_lm, ~1, var = 'delta_N_mass'))
test(emtrends(delta_N_mass_lm, ~1, var = 'delta_N_mass', at = list(delta_chi = -25)))
test(emtrends(delta_N_mass_lm, ~1, var = 'delta_N_mass', at = list(delta_chi = 0)))
test(emtrends(delta_N_mass_lm, ~1, var = 'delta_N_mass', at = list(delta_chi = 25)))

## make delta_plot
### dataset
delta_live_mass_plot_data =  delta_live_mass_data

### trendline information
# delta_live_mass_plot_intercept = summary(emmeans(delta_live_mass_lm, ~1, at = list(delta_live_mass = 0)))[1, 2]
# delta_live_mass_plot_intercept_Nfix_no = summary(emmeans(delta_live_mass_lm, ~Nfix, at = list(delta_live_mass = 0)))[1, 2]
# delta_live_mass_plot_intercept_Nfix_yes = summary(emmeans(delta_live_mass_lm, ~Nfix, at = list(delta_live_mass = 0)))[2, 2]
# delta_live_mass_plot_slope = summary(emtrends(delta_live_mass_lm, ~1, var = 'delta_live_mass'))[1, 2]
# delta_live_mass_plot_trend = delta_live_mass_plot_slope * seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
#   delta_live_mass_plot_intercept
# delta_live_mass_plot_trend_Nfix_no = delta_live_mass_plot_slope * seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
#   delta_live_mass_plot_intercept_Nfix_no
# delta_live_mass_plot_trend_Nfix_yes = delta_live_mass_plot_slope * seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
#   delta_live_mass_plot_intercept_Nfix_yes
# 
# delta_live_mass_plot_trend_df = data.frame(seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1),
#                                     delta_live_mass_plot_trend, delta_live_mass_plot_trend_Nfix_no, delta_live_mass_plot_trend_Nfix_yes)
# colnames(delta_live_mass_plot_trend_df) = c('delta_live_mass', 'delta_narea', 'delta_narea_Nfix_no', 'delta_narea_Nfix_yes')

delta_live_mass_plot_intercept_lowlma = summary(emmeans(delta_live_mass_lm, ~1, 
                                              at = list(delta_live_mass = 0, delta_lma = -25)))[1, 2]
delta_live_mass_plot_slope_lowlma = summary(emtrends(delta_live_mass_lm, ~1, var = 'delta_live_mass', 
                                           at = list(delta_lma = -25)))[1, 2]
delta_live_mass_plot_intercept_midlma = summary(emmeans(delta_live_mass_lm, ~1, 
                                              at = list(delta_live_mass = 0, delta_lma = 0)))[1, 2]
delta_live_mass_plot_slope_midlma = summary(emtrends(delta_live_mass_lm, ~1, var = 'delta_live_mass', 
                                           at = list(delta_lma = 0)))[1, 2]
delta_live_mass_plot_intercept_highlma = summary(emmeans(delta_live_mass_lm, ~1, 
                                              at = list(delta_live_mass = 0, delta_lma = 25)))[1, 2]
delta_live_mass_plot_slope_highlma = summary(emtrends(delta_live_mass_lm, ~1, var = 'delta_live_mass', 
                                           at = list(delta_lma = 25)))[1, 2]
delta_live_mass_plot_trend_lowlma = delta_live_mass_plot_slope_lowlma * 
  seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
  delta_live_mass_plot_intercept_lowlma
delta_live_mass_plot_trend_midlma = delta_live_mass_plot_slope_midlma * 
  seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
  delta_live_mass_plot_intercept_midlma
delta_live_mass_plot_trend_highlma = delta_live_mass_plot_slope_highlma * 
  seq(min(delta_live_mass_plot_data$delta_live_mass), max(delta_live_mass_plot_data$delta_live_mass), 1) +
  delta_live_mass_plot_intercept_highlma

delta_live_mass_plot_trend_df = data.frame(seq(min(delta_live_mass_plot_data$delta_live_mass), 
                                               max(delta_live_mass_plot_data$delta_live_mass), 1),
                                    delta_live_mass_plot_trend_lowlma, 
                                    delta_live_mass_plot_trend_midlma, 
                                    delta_live_mass_plot_trend_highlma)
colnames(delta_live_mass_plot_trend_df) = c('delta_live_mass', 
                                            'delta_narea_lowlma', 
                                            'delta_narea_midlma', 
                                            'delta_narea_highlma')


delta_live_mass_plot = ggplot(data = delta_live_mass_plot_data, 
       aes(x = delta_live_mass, y = delta_narea, fill = delta_lma, size = delta_lma)) +
  theme(legend.position = "none", 
        axis.title.y=element_text(size=rel(2.5), colour = 'black'),
        axis.title.x=element_text(size=rel(2.5), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_point(shape = 21, colour = 'black', stroke = 0.5, alpha = 0.7) +
  scale_size_continuous(range = c(1, 3)) +
  scale_fill_gradient(low = 'blue', high = 'red') +
  geom_line(data = delta_live_mass_plot_trend_df, 
            aes(x = delta_live_mass, y = delta_narea_lowlma, fill = NULL), 
            size = 3, colour = 'blue', alpha = 1, lty = 2) +
  geom_line(data = delta_live_mass_plot_trend_df, 
            aes(x = delta_live_mass, y = delta_narea_midlma, fill = NULL), 
            size = 3, colour = 'purple', alpha = 1, lty = 2) +
  geom_line(data = delta_live_mass_plot_trend_df, 
            aes(x = delta_live_mass, y = delta_narea_highlma, fill = NULL), 
            size = 3, colour = 'red', alpha = 1, lty = 1) +
  # geom_line(data = delta_live_mass_plot_trend_df, aes(x = delta_live_mass, y = delta_narea_Nfix_yes), 
  #           size = 4, colour = 'green', alpha = 0.7) +
  # ylim(c(-50, 50)) +
  ylab(expression('∆' * italic('N')['area'] * ' (%)')) +
  xlab(expression('∆' *'AG Biomass' * ' (%)'))

jpeg(filename = "plots/delta_live_mass_plot.jpeg", width = 400, height = 450, units = 'px')
plot(delta_live_mass_plot)
dev.off()


