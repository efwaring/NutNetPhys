# Script using least-cost to explore Narea response at nutnet

#####################################################################
## load packages
#####################################################################
library(tidyverse)
library(lme4)
library(car)
library(r2glmm)
library(treemapify)

#####################################################################
## load functions
#####################################################################
source('optimal_vcmax_R/calc_optimal_vcmax.R')
source('optimal_vcmax_R/calc_optimal_vcmax_knownchi.R')
sourceDirectory('optimal_vcmax_R/functions')
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
leaf_Nonly = subset(leaf, Ptrt_fac == '0' & Ktrt_fac == '0')

### run analyses
# hist(leaf$narea)
# hist(log(leaf$narea))
leafNarea_lmer = lmer(log(narea) ~ Ntrt_fac + 
                        (1|site_code), 
                      data = leaf_Nonly)
# plot(resid(leafNarea_lmer) ~ fitted(leafNarea_lmer))
Anova(leafNarea_lmer)
cld(emmeans(leafNarea_lmer, ~Ntrt_fac))
(exp(summary(emmeans(leafNarea_lmer, ~Ntrt_fac))[2,2]) - exp(summary(emmeans(leafNarea_lmer, ~Ntrt_fac))[1,2])) / 
  exp(abs(summary(emmeans(leafNarea_lmer, ~Ntrt_fac))[1,2]))
r2beta(leafNarea_lmer, partial = T)

# hist(leaf$chi)
# hist(log(leaf$Narea))
leafchi_lmer = lmer(chi ~ Ntrt_fac + 
                        (1|site_code), 
                      data = leaf_Nonly)
# plot(resid(leafNarea_lmer) ~ fitted(leafNarea_lmer))
Anova(leafchi_lmer)
cld(emmeans(leafchi_lmer, ~Ntrt_fac))
(summary(emmeans(leafchi_lmer, ~Ntrt_fac))[2,2] - summary(emmeans(leafchi_lmer, ~Ntrt_fac))[1,2]) / 
  abs(summary(emmeans(leafchi_lmer, ~Ntrt_fac))[1,2])
r2beta(leafchi_lmer, partial = T)

### make figure
narea_plot = ggplot(data = leaf_Nonly, 
                         aes(x = Ntrt_fac, y = log(narea))) +
  theme(legend.position = "none", 
        axis.title.y=element_text(size=rel(2.5), colour = 'black'),
        axis.title.x=element_text(size=rel(2.5), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_boxplot(outlier.color = NA, fill = 'white') +
  geom_dotplot(binaxis = 'y', binwidth = 0.1, stackdir = 'center', alpha = 0.5) +
  scale_x_discrete(labels = c('Ambient', 'Added N')) +
  xlab('') +
  ylab(expression('Log leaf ' * italic('N')['area'])) +
  annotate("text", x = 1.5, y = 6.7, label = "p = 0.82", size = 8)

chi_plot = ggplot(data = leaf_Nonly, 
                  aes(x = Ntrt_fac, y = chi)) +
  theme(legend.position = "none", 
        axis.title.y=element_text(size=rel(2.5), colour = 'black'),
        axis.title.x=element_text(size=rel(2.5), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_boxplot(outlier.color = NA, fill = 'white') +
  geom_dotplot(binaxis = 'y', binwidth = 0.005, stackdir = 'center', alpha = 0.5) +
  scale_x_discrete(labels = c('Ambient', 'Added N')) +
  xlab('Nitrogen treatment') +
  ylab(expression('C'[i] * '/C'[a])) +
  annotate("text", x = 1.5, y = 0.93, label = "p = 0.05", size = 8)

jpeg(filename = "plots/leaf_plot.jpeg", width = 550, height = 900, units = 'px')
multiplot(narea_plot, chi_plot, cols=1)
dev.off()

#####################################################################
## soil N effects on LAI
#####################################################################
leaf_site = read.csv('../Data/processed/traits_site.csv')
leaf_site$Ntrt_fac = as.factor(leaf_site$Ntrt)
leaf_site$Ptrt_fac = as.factor(leaf_site$Ptrt)
leaf_site$Ktrt_fac = as.factor(leaf_site$Ktrt)
leaf_site$block_fac = as.factor(leaf_site$block)
leaf_site_Nonly = subset(leaf_site, Ptrt_fac == '0' & Ktrt_fac == '0')

# hist(leaf_site_Nonly$lai_mean)
# hist(log(leaf_site_Nonly$lai_mean))
lai_lmer = lmer(log(lai_mean) ~ Ntrt_fac + 
                      (1|site_code) + (1|site_code:block), 
                    data = leaf_site_Nonly)
# plot(resid(lai_lmer) ~ fitted(lai_lmer))
Anova(lai_lmer)
cld(emmeans(lai_lmer, ~Ntrt_fac))
(exp(summary(emmeans(lai_lmer, ~Ntrt_fac))[2,2]) - exp(summary(emmeans(lai_lmer, ~Ntrt_fac))[1,2])) / 
  abs(exp(summary(emmeans(lai_lmer, ~Ntrt_fac))[1,2]))

# hist(leaf_site_Nonly$live_mass_mean)
# hist(log(leaf_site_Nonly$live_mass_mean))
live_mass_lmer = lmer(log(live_mass_mean) ~ Ntrt_fac + 
                  (1|site_code) + (1|site_code:block), 
                data = leaf_site_Nonly)
# plot(resid(live_mass_lmer) ~ fitted(live_mass_lmer))
Anova(live_mass_lmer)
cld(emmeans(live_mass_lmer, ~Ntrt_fac))
(exp(summary(emmeans(live_mass_lmer, ~Ntrt_fac))[2,2]) - exp(summary(emmeans(live_mass_lmer, ~Ntrt_fac))[1,2])) / 
  abs(exp(summary(emmeans(live_mass_lmer, ~Ntrt_fac))[1,2]))

lai_plot = ggplot(data = leaf_site_Nonly, 
                  aes(x = Ntrt_fac, y = log(lai_mean))) +
  theme(legend.position = "none", 
        axis.title.y=element_text(size=rel(2.5), colour = 'black'),
        axis.title.x=element_text(size=rel(2.5), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_boxplot(outlier.color = NA, fill = 'white') +
  geom_dotplot(binaxis = 'y', binwidth = 0.1, stackdir = 'center', alpha = 0.5) +
  scale_x_discrete(labels = c('Ambient', 'Added N')) +
  xlab('') +
  ylab(expression('Log LAI (m' ^ '2' *' m' ^'-2' * ')')) +
  annotate("text", x = 1.5, y = 1.3, label = "p < 0.05", size = 8)

live_mass_plot = ggplot(data = leaf_site_Nonly, 
                  aes(x = Ntrt_fac, y = log(live_mass_mean))) +
  theme(legend.position = "none", 
        axis.title.y=element_text(size=rel(2.5), colour = 'black'),
        axis.title.x=element_text(size=rel(2.5), colour = 'black'),
        axis.text.x=element_text(size=rel(2), colour = 'black'),
        axis.text.y=element_text(size=rel(2), colour = 'black'),
        panel.background = element_rect(fill = 'white', colour = 'black'),
        panel.grid.major = element_line(colour = "grey")) +
  geom_boxplot(outlier.color = NA, fill = 'white') +
  geom_dotplot(binaxis = 'y', binwidth = 0.07, stackdir = 'center', alpha = 0.5) +
  scale_x_discrete(labels = c('Ambient', 'Added N')) +
  xlab('Nitrogen treatment') +
  ylab(expression('Log Biomass (g m' ^ '2' * ')')) +
  annotate("text", x = 1.5, y = 7.3, label = "p < 0.01", size = 8)

jpeg(filename = "plots/plant_plot.jpeg", width = 550, height = 900, units = 'px')
multiplot(lai_plot, live_mass_plot, cols=1)
dev.off()

#####################################################################
## Narea predictions
#####################################################################
gas_exchange_pred = calc_optimal_vcmax_knownchi(tg_c = leaf_Nonly$tmp, 
                                       paro = leaf_Nonly$par_per_leaf_area, 
                                       cao = 400, 
                                       vpdo = leaf_Nonly$vpd, 
                                       z = leaf_Nonly$z,
                                       chi = leaf_Nonly$chi)

leaf_Nonly$vcmax25 = gas_exchange_pred$vcmax /
  calc_vcmax_tresp_mult(gas_exchange_pred$tg_c, gas_exchange_pred$tg_c, 25)
leaf_Nonly$jmax25 = gas_exchange_pred$jmax /
  calc_jmax_tresp_mult(gas_exchange_pred$tg_c, gas_exchange_pred$tg_c, 25)

leaf_Nonly$nrubisco = fvcmax25_nrubisco(leaf_Nonly$vcmax25)
leaf_Nonly$nbioe = fjmax25_nbioe(leaf_Nonly$jmax25)
leaf_Nonly$nstructure = flma_nstructure(leaf_Nonly$lma)
leaf_Nonly$nall = leaf_Nonly$nrubisco + leaf_Nonly$nbioe + leaf_Nonly$nstructure

npred_lm = lm(narea ~ nall + Ntrt_fac, data = leaf_Nonly)
Anova(npred_lm)
r2beta(npred_lm, partial = T)

leafNarea_traits_lmer = lmer(log(narea) ~ 
                   chi + lma +
                   (1|site_code), 
                 data = leaf_Nonly)
Anova(leafNarea_traits_lmer)
plot(resid(leafNarea_traits_lmer) ~ fitted(leafNarea_traits_lmer))
emtrends(leafNarea_traits_lmer, ~1, var = 'chi')
emtrends(leafNarea_traits_lmer, ~1, var = 'lma')
r2beta(leafNarea_traits_lmer, partial = T)


