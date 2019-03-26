# Goal: examine above ground N data and test hypotheses related to optimal N allocation
#
# Lizz Waring
#
# General strategy: ask questions sequentially from least to most complex

#########################################
# packages necessary
#########################################

library(lme4)
library(car)
library(emmeans)
library(dplyr)
library(plyr)
library(RColorBrewer)
library(afex)

#########################################
# read in and modify biomass data
#########################################
source("./Analysis/create_biomass_data.R")

#########################################
# Q1: does Above ground N differ amongst treatments
#########################################

hist(biomass_core_spei$AGN)
AGN_mixed = lmer(AGN ~ Ntrt_fac * Nfix+ p_pet + (1|site_code) + (1|category) ,
                data = biomass_core_spei)
plot(resid(AGN_mixed) ~ fitted(AGN_mixed))
Anova(AGN_mixed)
emmeans(AGN_mixed, ~Ntrt_fac)
cld(emmeans(AGN_mixed, ~Ntrt_fac * Nfix))

hist(full_biomass$AGN)
AGN_full_lmer = lmer(AGN ~ Ntrt_fac * Ptrt_fac * Ktrt_fac * Nfix + p_pet + (1|site_code) + (1|category) + (1|site_code:category), data = biomass_core_spei)
plot(resid(AGN_full_lmer) ~ fitted(AGN_full_lmer))
Anova(AGN_full_lmer)
cld(emmeans(AGN_full_lmer, ~Ntrt_fac))
cld(emmeans(AGN_full_lmer, ~Ntrt_fac * Nfix)) 



#########################################
# A1: Yes, unsurpringinly AGN increases with treatments, no effect of legumes
# 
#########################################

#########################################
# Q2: does the change in AGNdiffer across climates?
#########################################
#
#
# calculate ∆AGN

# aggregate by site, N treatment
# remove fixers

biomass_core_spei_nofix = subset(biomass_core_spei, Nfix == 'no')

biomass_core_spei_nofix_group = group_by(biomass_core_spei_nofix, site_code, year, Ntrt_fac)
biomass_core_spei_nofix_mean = summarise_all(biomass_core_spei_nofix_group, .funs = funs(Mean = mean(., na.rm = T), SD = sd(., na.rm = T)))

biomass_core_spei_nofix_mean_yN = subset(biomass_core_spei_nofix_mean, Ntrt_fac == '1')
biomass_core_spei_nofix_mean_nN = subset(biomass_core_spei_nofix_mean, Ntrt_fac == '0')

nrow(biomass_core_spei_nofix_mean_yN)
nrow(biomass_core_spei_nofix_mean_nN)

deltaN = ((biomass_core_spei_nofix_mean_yN$AGN_Mean - biomass_core_spei_nofix_mean_nN$AGN_Mean) / biomass_core_spei_nofix_mean_nN$AGN_Mean) * 100


summary(lm(deltaN ~ biomass_core_spei_nofix_mean_nN$p_pet_Mean)) # nothing
summary(lm(deltaN ~ biomass_core_spei_nofix_mean_yN$p_pet_Mean)) # nothing
summary(lm(deltaN ~ biomass_core_spei_nofix_mean_nN$par_rat_Mean)) # nothing
summary(lm(deltaN ~ biomass_core_spei_nofix_mean_yN$par_rat_Mean)) 

plot(deltaN ~ biomass_core_spei_nofix_mean_nN$p_pet_Mean)
plot(deltaN ~ biomass_core_spei_nofix_mean_nN$par_rat_Mean)

summary(lm(deltaN ~ biomass_core_spei_nofix_mean_nN$precip.mean_Mean)) # nothing
summary(lm(deltaN ~ biomass_core_spei_nofix_mean_yN$precip.mean_Mean)) # nothing

plot(deltaN ~ biomass_core_spei_nofix_mean_nN$precip.mean_Mean)
plot(deltaN ~ biomass_core_spei_nofix_mean_nN$precip.mean_Mean)

lm_deltaN = lm(deltaN ~ biomass_core_spei_nofix_mean_nN$p_pet_Mean * biomass_core_spei_nofix_mean_nN$par_rat_Mean)
anova(lm_deltaN)
summary(lm_deltaN)

#########################################
# A2: does the change in AGN differ across climates? nope, not really! basically no change regardless of site/climate
#########################################
#


#########################################
# Let's make some plots
#########################################

# N fert effect on Narea
par(mfrow = c(1, 1), oma = c(4, 4, 1, 1))
leafNarea_lmer_sum = summary(emmeans(leafNarea_lmer, ~Ntrt_fac * Nfix))
leafNarea_boxplot = boxplot(log(subset(biomass_core_spei, Ntrt_fac == '0' & Nfix == 'no')$Narea), log(subset(biomass_core_spei, Ntrt_fac == '1' & Nfix == 'no')$Narea), log(subset(biomass_core_spei, Ntrt_fac == '0' & Nfix == 'yes')$Narea), log(subset(biomass_core_spei, Ntrt_fac == '1' & Nfix == 'yes')$Narea), ylim = c(-20, 0), yaxt = 'n', xaxt = 'n', ylab = '', xlab = '', col = c('white', 'white', 'yellow', 'yellow'))
points(leafNarea_lmer_sum[,3], pch = 16, col = c('purple', 'purple', 'red', 'red'), cex = 2)
segments(seq(1, 4, 1), leafNarea_lmer_sum[,3] + leafNarea_lmer_sum[,4] * 1.96, seq(1, 4, 1), leafNarea_lmer_sum[,3] - leafNarea_lmer_sum[,4] * 1.96, col = c('purple', 'purple', 'red', 'red'), lwd = 6)
axis(1, at = seq(1, 4, 1), c("-N", '+N', '-N', '+N'), cex.axis = 2)
text(1.5, 0, 'non-fixers', cex = 2)
text(3.5, 0, 'N-fixers', cex = 2)
axis(2, seq(-20, 0, 5), cex.axis = 2, las = 1)
mtext(side = 2, 'log(Narea)', cex = 3, line = 5)

# climate effect on ∆N
par(mfrow = c(1, 1), oma = c(4, 4, 1, 1))
palette = colorRampPalette(brewer.pal(9,'Blues'))
p_pet_cols <- palette(9)[as.numeric(cut(leaf_core_spei_nofix_mean_nN$p_pet_Mean,breaks = 9))]
plot(deltaNarea ~ leaf_core_spei_nofix_mean_nN$par_rat_Mean, pch =21, cex = 3, bg = p_pet_cols,  yaxt = 'n', xaxt = 'n', ylab = '', xlab = '', xlim = c(0, 1), ylim = c(-100, 100))
axis(1, seq(0, 1, 0.2), cex.axis = 2)
axis(2, seq(-100, 100, 50), cex.axis = 2, las = 1)
mtext(side = 1, 'Light interception', line = 4, cex = 3)
mtext(side = 2, '∆Narea', line = 5, cex = 3)

# mean(deltaNarea, na.rm = T) # 12.1
# (sd(deltaNarea, na.rm = T) / sqrt(25)) * 1.96 # 10.8
# text(0.5, 90, 'Mean ± 95% CI: 12.1 ± 10.8', cex = 1.5)

#########################################
# Q3: can we predict AGN?
#########################################

# does N area increase with abs(latitude)?

plot(biomass_core_spei$AGN ~ abs(biomass_core_spei$latitude))
summary(lm(log(biomass_core_spei$AGN) ~ abs(biomass_core_spei$latitude))) # increases with latitude (temperature effect)

# what about with climate?

# read in gs climate (>0°C)

biomass_core_spei$tmp = NA
biomass_core_spei$par = NA
biomass_core_spei$vpd = NA
biomass_core_spei$z = NA

climate_df = c()
for (i in 1:nrow(biomass_core_spei)){
  
  currentLat = biomass_core_spei$latitude[i]
  currentLon = biomass_core_spei$longitude[i]
  
  clim_comb = tmp_globe
  latClim = clim_comb[ , 3]
  lonClim = clim_comb[ , 2]
  
  best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
  best_lat=latClim[best_lat_pos[1]]
  best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
  best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
  
  tmp = subset(clim_comb, lat==best_lat & lon==best_lon)$tmp
  
  clim_comb = par_globe
  latClim = clim_comb[ , 3]
  lonClim = clim_comb[ , 2]
  
  best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
  best_lat=latClim[best_lat_pos[1]]
  best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
  best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
  
  par = subset(clim_comb, lat==best_lat & lon==best_lon)$par
  
  clim_comb = vpd_globe
  latClim = clim_comb[ , 3]
  lonClim = clim_comb[ , 2]
  
  best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
  best_lat=latClim[best_lat_pos[1]]
  best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
  best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
  
  vpd = subset(clim_comb, lat==best_lat & lon==best_lon)$vpd
  
  clim_comb = z_globe
  latClim = clim_comb[ , 3]
  lonClim = clim_comb[ , 2]
  
  best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
  best_lat=latClim[best_lat_pos[1]]
  best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
  best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
  
  z = subset(clim_comb, lat==best_lat & lon==best_lon)$z
  
  climate_df = rbind(climate_df, c(tmp, par, vpd, z, currentLat, currentLon, best_lat, best_lon, i))
  
  
}

plot(climate_df[,5] ~ climate_df[,7])
plot(climate_df[,6] ~ climate_df[,8])

biomass_core_spei$tmp = climate_df[, 1]
biomass_core_spei$par = climate_df[, 2]
biomass_core_spei$vpd = climate_df[, 3]
biomass_core_spei$z = climate_df[, 4]

biomass_core_spei_by_site_Ntrt = group_by(biomass_core_spei, site_code, Ntrt_fac)
biomass_core_spei_by_site_Ntrt_mean = summarise_all(biomass_core_spei_by_site_Ntrt, .funs = funs(Mean = mean(., na.rm = T), SD = sd(., na.rm = T)))

### ok let's try to compare predictions
AGN_lmer_clim = lm(log(AGN_Mean) ~ tmp_Mean + par_Mean + Ntrt_fac, data = biomass_core_spei_by_site_Ntrt_mean)
plot(resid(AGN_lmer_clim) ~ fitted(AGN_lmer_clim))
Anova(AGN_lmer_clim) # just Light!
emtrends(AGN_lmer_clim, ~1, var = 'par_Mean')
emmeans(AGN_lmer_clim, ~ par_Mean, at = list(par_Mean = 0))

# try to predict using optimal Vcmax


plot(log(AGN_Mean) ~ par_Mean, data = biomass_core_spei_by_site_Ntrt_mean) # denser canopy more AGN



