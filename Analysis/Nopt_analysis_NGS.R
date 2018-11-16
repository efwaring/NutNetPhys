# Goal: examine leaf N data and test hypotheses related to optimal N allocation
#
# Nick Smith
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

#########################################
# read in and modify foliar trait dataset (Jen Firn)
#########################################

leaf = read.csv('/Users/nicksmith/Dropbox/NutNet_meeting_2018/Nopt_analyses/Data/NutNet-foliar-traits-7JAN2017.csv')

leaf[leaf == 'NULL'] <- NA

leaf[, 8:30] <- sapply(leaf[, 8:30], as.character)
leaf[, 8:30] <- sapply(leaf[, 8:30], as.numeric)

leaf$Ntrt = 0
leaf$Ntrt[leaf$trt == 'N' | leaf$trt == 'NP' | leaf$trt == 'NK' | leaf$trt == 'NPK' | leaf$trt == 'NPK+Fence'] = 1
leaf$Ptrt = 0
leaf$Ptrt[leaf$trt == 'P' | leaf$trt == 'NP' | leaf$trt == 'PK' | leaf$trt == 'NPK' | leaf$trt == 'NPK+Fence'] = 1
leaf$Ktrt = 0
leaf$Ktrt[leaf$trt == 'K' | leaf$trt == 'NK' | leaf$trt == 'PK' | leaf$trt == 'NPK' | leaf$trt == 'NPK+Fence'] = 1

leaf$Ntrt_fac = as.factor(leaf$Ntrt)
leaf$Ptrt_fac = as.factor(leaf$Ptrt)
leaf$Ktrt_fac = as.factor(leaf$Ktrt)

leaf$Narea = leaf$leaf_pct_N * (1/leaf$SLA)

leaf$Nfix = 'no'
leaf$Nfix[leaf$Family == 'Fabaceae'] = 'yes'

# add core data to the dataset
core = read.csv('/Users/nicksmith/Dropbox/NutNet_meeting_2018/Nopt_analyses/Data/comb-by-plot-09-April-2018.csv')
core[core == 'NULL'] <- NA

core[, 26:45] <- sapply(core[, 26:45], as.character)
core[, 26:45] <- sapply(core[, 26:45], as.numeric)

core$par_rat = core$Ground_PAR / core$Ambient_PAR
leaf_core = join(leaf, core, by = c("site_code", "year", "block", "plot", "trt"), type = 'left', match = 'first')

# add SPEI to "leaf_core" data
spei = read.csv('/Users/nicksmith/Dropbox/NutNet_meeting_2018/Nopt_analyses/Data/CRU-annual_2018-07-06.csv')
spei$p_pet = spei$precip / spei$PET

leaf_core_spei = join(leaf_core, spei, by = c("site_code", "year"), type = 'left', match = 'first')

# check
nrow(leaf)
nrow(leaf_core) # some core data matches the leaf data twice...
nrow(leaf_core_spei)

#########################################
# Q1: does leaf N differ amongst treatments
#########################################

hist(leaf_core_spei $leaf_pct_N)
leafNper_lmer = lmer(leaf_pct_N ~ Ntrt_fac * Ptrt_fac * Ktrt_fac * Nfix + p_pet + (1|site_code) + (1|Family) + (1|Taxon) + (1|site_code:Family) + (1|site_code:Taxon), data = leaf_core_spei)
plot(resid(leafNper_lmer) ~ fitted(leafNper_lmer))
Anova(leafNper_lmer)
emmeans(leafNper_lmer, ~Ntrt_fac)
cld(emmeans(leafNper_lmer, ~Ntrt_fac * Ptrt_fac))
cld(emmeans(leafNper_lmer, ~Ntrt_fac * Nfix))

hist(leaf$Narea)
hist(log(leaf$Narea))
leafNarea_lmer = lmer(log(Narea) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac * Nfix + p_pet + (1|site_code) + (1|Family) + (1|Taxon) + (1|site_code:Family) + (1|site_code:Taxon), data = leaf_core_spei)
plot(resid(leafNarea_lmer) ~ fitted(leafNarea_lmer))
Anova(leafNarea_lmer)
cld(emmeans(leafNarea_lmer, ~Ntrt_fac))
cld(emmeans(leafNarea_lmer, ~Ntrt_fac * Nfix)) # ~2% increase in non-N fixers

hist(leaf$SLA_num)
hist(log(leaf$SLA_num))
SLA_lmer = lmer(log(leaf$SLA) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac * Nfix + p_pet + (1|site_code) + (1|Family) + (1|Taxon) + (1|site_code:Family) + (1|site_code:Taxon), data = leaf_core_spei)
plot(resid(SLA_lmer) ~ fitted(SLA_lmer))
Anova(SLA_lmer)
cld(emmeans(SLA_lmer, ~Ntrt_fac)) # (~1% increase)

#########################################
# A1: yes, N is a little higher (P < 0.001), however the actual effect is quite small (2%). No effect with N fixers. SLA is slightly increased with soil N as well (1%). PLants seem to be "moderately" building larger, higher N leaves with N...not a big effect
# 
#########################################

#########################################
# Q2: does the change in leaf N differ across climates?
#########################################
#
#
# calculate ∆leafN

# aggregate by site, N treatment
# remove fixers

leaf_core_spei_nofix = subset(leaf_core_spei, Nfix == 'no')

leaf_core_spei_nofix_group = group_by(leaf_core_spei_nofix, site_code, year, Ntrt_fac)
leaf_core_spei_nofix_mean = summarise_all(leaf_core_spei_nofix_group, .funs = funs(Mean = mean(., na.rm = T), SD = sd(., na.rm = T)))

leaf_core_spei_nofix_mean_yN = subset(leaf_core_spei_nofix_mean, Ntrt_fac == '1')
leaf_core_spei_nofix_mean_nN = subset(leaf_core_spei_nofix_mean, Ntrt_fac == '0')

nrow(leaf_core_spei_nofix_mean_yN)
nrow(leaf_core_spei_nofix_mean_nN)

deltaNper = ((leaf_core_spei_nofix_mean_yN$leaf_pct_N_Mean - leaf_core_spei_nofix_mean_nN$leaf_pct_N_Mean) / leaf_core_spei_nofix_mean_nN$leaf_pct_N_Mean) * 100

deltaNarea = ((leaf_core_spei_nofix_mean_yN$Narea_Mean - leaf_core_spei_nofix_mean_nN$Narea_Mean) / leaf_core_spei_nofix_mean_nN$Narea_Mean) * 100

summary(lm(deltaNper ~ leaf_core_spei_nofix_mean_nN$p_pet_Mean)) # nothing
summary(lm(deltaNarea ~ leaf_core_spei_nofix_mean_nN$p_pet_Mean)) # nothing
summary(lm(deltaNper ~ leaf_core_spei_nofix_mean_nN$p_pet_Mean)) # nothing
summary(lm(deltaNarea ~ leaf_core_spei_nofix_mean_nN$p_pet_Mean)) # nothing

plot(deltaNper ~ leaf_core_spei_nofix_mean_nN$p_pet_Mean)
plot(deltaNarea ~ leaf_core_spei_nofix_mean_nN$p_pet_Mean)
plot(deltaNper ~ leaf_core_spei_nofix_mean_nN$par_rat_Mean)
plot(deltaNarea ~ leaf_core_spei_nofix_mean_nN$par_rat_Mean)

summary(lm(deltaNper ~ leaf_core_spei_nofix_mean_nN$precip.mean_Mean)) # nothing
summary(lm(deltaNarea ~ leaf_core_spei_nofix_mean_nN$precip.mean_Mean)) # nothing

plot(deltaNper ~ leaf_core_spei_nofix_mean_nN$precip.mean_Mean)
plot(deltaNarea ~ leaf_core_spei_nofix_mean_nN$precip.mean_Mean)

lm_deltaNarea = lm(deltaNarea ~ leaf_core_spei_nofix_mean_nN$p_pet_Mean * leaf_core_spei_nofix_mean_nN$par_rat_Mean)
anova(lm_deltaNarea)
summary(lm_deltaNarea)

#########################################
# A2: does the change in leaf N differ across climates? nope, not really! basically no change regardless of site/climate
#########################################
#


#########################################
# Let's make some plots
#########################################

# N fert effect on Narea
par(mfrow = c(1, 1), oma = c(4, 4, 1, 1))
leafNarea_lmer_sum = summary(emmeans(leafNarea_lmer, ~Ntrt_fac * Nfix))
leafNarea_boxplot = boxplot(log(subset(leaf_core_spei, Ntrt_fac == '0' & Nfix == 'no')$Narea), log(subset(leaf_core_spei, Ntrt_fac == '1' & Nfix == 'no')$Narea), log(subset(leaf_core_spei, Ntrt_fac == '0' & Nfix == 'yes')$Narea), log(subset(leaf_core_spei, Ntrt_fac == '1' & Nfix == 'yes')$Narea), ylim = c(-20, 0), yaxt = 'n', xaxt = 'n', ylab = '', xlab = '', col = c('white', 'white', 'yellow', 'yellow'))
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
# Q3: can we predict Narea?
#########################################

# does N area increase with abs(latitude)?

plot(leaf_core_spei$Narea ~ abs(leaf_core_spei$latitude))
summary(lm(log(leaf_core_spei$Narea) ~ abs(leaf_core_spei$latitude))) # increases with latitude (temperature effect)

# what about with climate?

# read in gs climate (>0°C)
tmp_globe = read.csv('/Users/nicksmith/Dropbox/0Main/Research/Colimitation/Spatial_Maps/cru_tmp_climExtract_growingseason_globe.csv')
par_globe = read.csv('/Users/nicksmith/Dropbox/0Main/Research/Colimitation/Spatial_Maps/cru_par_climExtract_growingseason_globe.csv')
vpd_globe = read.csv('/Users/nicksmith/Dropbox/0Main/Research/Colimitation/Spatial_Maps/cru_vpd_climExtract_growingseason_globe.csv')
z_globe =  read.csv('/Users/nicksmith/Dropbox/0Main/Research/Colimitation/Spatial_Maps/z_globe.csv')

leaf_core_spei$tmp = NA
leaf_core_spei$par = NA
leaf_core_spei$vpd = NA
leaf_core_spei$z = NA

climate_df = c()
for (i in 1:nrow(leaf_core_spei)){
	
	currentLat = leaf_core_spei$latitude[i]
	currentLon = leaf_core_spei$longitude[i]
	
	clim_comb = tmp_globe
	latClim = climb_comb[ , 3]
	lonClim = climb_comb[ , 2]
	
	best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
    best_lat=latClim[best_lat_pos[1]]
    best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
    best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
    
    tmp = subset(clim_comb, lat==best_lat & lon==best_lon)$tmp
    
    clim_comb = par_globe
	latClim = climb_comb[ , 3]
	lonClim = climb_comb[ , 2]
	
	best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
    best_lat=latClim[best_lat_pos[1]]
    best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
    best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
    
    par = subset(clim_comb, lat==best_lat & lon==best_lon)$par
    
    clim_comb = vpd_globe
	latClim = climb_comb[ , 3]
	lonClim = climb_comb[ , 2]
	
	best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
    best_lat=latClim[best_lat_pos[1]]
    best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
    best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
    
    vpd = subset(clim_comb, lat==best_lat & lon==best_lon)$vpd
    
    clim_comb = z_globe
	latClim = climb_comb[ , 3]
	lonClim = climb_comb[ , 2]
	
	best_lat_pos=which(abs(latClim-currentLat)==min(abs(latClim-currentLat)))
    best_lat=latClim[best_lat_pos[1]]
    best_lon_pos=which(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)==min(abs(subset(clim_comb, lat== best_lat)$lon-currentLon)))
    best_lon=subset(clim_comb, lat== best_lat)$lon[best_lon_pos[1]]
    
    z = subset(clim_comb, lat==best_lat & lon==best_lon)$z
    
    climate_df = rbind(climate_df, c(tmp, par, vpd, z, currentLat, currentLon, best_lat, best_lon, i))
	
	
}

plot(climate_df[,5] ~ climate_df[,7])
plot(climate_df[,6] ~ climate_df[,8])

leaf_core_spei$tmp = climate_df[, 1]
leaf_core_spei$par = climate_df[, 2]
leaf_core_spei$vpd = climate_df[, 3]
leaf_core_spei$z = climate_df[, 4]

leaf_core_spei_by_site_Ntrt = group_by(leaf_core_spei, site_code, Ntrt_fac)
leaf_core_spei_by_site_Ntrt_mean = summarise_all(leaf_core_spei_by_site_Ntrt, .funs = funs(Mean = mean(., na.rm = T), SD = sd(., na.rm = T)))

### ok let's try to compare predictions
leafNarea_lmer_clim = lm(log(Narea_Mean) ~ tmp_Mean + par_Mean + Ntrt_fac, data = leaf_core_spei_by_site_Ntrt_mean)
plot(resid(leafNarea_lmer_clim) ~ fitted(leafNarea_lmer_clim))
Anova(leafNarea_lmer_clim) # just temperature
emtrends(leafNarea_lmer_clim, ~1, var = 'tmp_Mean')
emmeans(leafNarea_lmer_clim, ~ tmp_Mean, at = list(tmp_Mean = 0))

# try to predict using optimal Vcmax


plot(log(Narea_Mean) ~ tmp_Mean, data = leaf_core_spei_by_site_Ntrt_mean)
abline(a = -0.2185934, b = -4.55581)



