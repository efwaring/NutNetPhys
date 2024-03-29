# Goal: examine leaf N data and test hypotheses related to optimal N allocation
#
# Nick Smith
#
# General strategy: ask questions sequentially from least to most complex

# setwd("/Users/nicksmith/Documents/Git/NutNetPhys/Analysis")

#########################################
# packages necessary
#########################################

library(lme4)
library(car)
library(emmeans)
library(dplyr)
library(plyr)
library(RColorBrewer)
library(multcompView)
library(ggplot2)

#########################################
# read in and modify foliar trait dataset (Jen Firn)
#########################################

leaf = read.csv('../Data/NutNet-foliar-traits-7JAN2017.csv')

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
core = read.csv('../Data/comb-by-plot-09-April-2018.csv')
core[core == 'NULL'] <- NA

core[, 26:45] <- sapply(core[, 26:45], as.character)
core[, 26:45] <- sapply(core[, 26:45], as.numeric)

core$par_rat = core$Ground_PAR / core$Ambient_PAR
leaf_core = join(leaf, core, by = c("site_code", "year", "block", "plot", "trt"), type = 'left', match = 'first')

# add SPEI to "leaf_core" data
spei = read.csv('../Data/CRU-annual_2018-07-06.csv')
spei$p_pet = spei$precip / spei$PET

leaf_core_spei = join(leaf_core, spei, by = c("site_code", "year"), type = 'left', match = 'first')

# check
nrow(leaf)
nrow(leaf_core) # some core data matches the leaf data twice...
nrow(leaf_core_spei)

#########################################
# Q1: does leaf N differ amongst treatments
#########################################

# hist(leaf_core_spei $leaf_pct_N)
# leafNper_lmer = lmer(leaf_pct_N ~ Ntrt_fac + 
#                        (1|site_code) + (1|Taxon) + (1|site_code:Taxon), 
#                      data = leaf_core_spei)
# plot(resid(leafNper_lmer) ~ fitted(leafNper_lmer))
# Anova(leafNper_lmer)
# emmeans(leafNper_lmer, ~Ntrt_fac)
# cld(emmeans(leafNper_lmer, ~Ntrt_fac)) # percent N is higher in N addition by ~12%

hist(leaf$Narea)
hist(log(leaf$Narea))
leafNarea_lmer = lmer(log(Narea) ~ Ntrt_fac * p_pet * Ambient_PAR + 
                        (1|site_code), 
                      data = leaf_core_spei)
# plot(resid(leafNarea_lmer) ~ fitted(leafNarea_lmer))
Anova(leafNarea_lmer)
cld(emmeans(leafNarea_lmer, ~Ntrt_fac))
(summary(emmeans(leafNarea_lmer, ~Ntrt_fac))[2,2] - summary(emmeans(leafNarea_lmer, ~Ntrt_fac))[1,2]) / abs(summary(emmeans(leafNarea_lmer, ~Ntrt_fac))[1,2])
test(emtrends(leafNarea_lmer, ~1, var = 'p_pet'))
test(emtrends(leafNarea_lmer, ~1, var = 'Ambient_PAR'))

hist(leaf$SLA_num)
hist(log(leaf$SLA_num))
SLA_lmer = lmer(log(leaf$SLA) ~ Ntrt_fac * p_pet * Ambient_PAR + 
                  (1|site_code), 
                data = leaf_core_spei)
# plot(resid(leafNarea_lmer) ~ fitted(leafNarea_lmer))
# plot(resid(SLA_lmer) ~ fitted(SLA_lmer))
Anova(SLA_lmer)
cld(emmeans(SLA_lmer, ~Ntrt_fac))
(summary(emmeans(SLA_lmer, ~Ntrt_fac))[2,2] - summary(emmeans(SLA_lmer, ~Ntrt_fac))[1,2]) / abs(summary(emmeans(SLA_lmer, ~Ntrt_fac))[1,2])
test(emtrends(SLA_lmer, ~1, var = 'p_pet'))
test(emtrends(SLA_lmer, ~1, var = 'Ambient_PAR'))

#########################################
# A1: yes, Narea is higher (P < 0.001), 
# however the actual effect is quite small (<2% increase)
# Plants seem to be "moderately" building larger, higher N leaves with N...not a big effect
# no interaction with climate or light availability
#########################################
# 
# #########################################
# # Q2: does the change in leaf N differ across climates?
# #########################################
# #
# #
# # calculate ∆leafN
# 
# # aggregate by site, N treatment
# 
# # leaf_core_spei_nofix = subset(leaf_core_spei, Nfix == 'no')
# 
# leaf_core_spei_group = group_by(leaf_core_spei, site_code, year, Ntrt_fac)
# leaf_core_spei_mean = summarise_all(leaf_core_spei_group, 
#                                     .funs = funs(Mean = mean(., na.rm = T), SD = sd(., na.rm = T)))
# 
# leaf_core_spei_mean_yN = subset(leaf_core_spei_mean, Ntrt_fac == '1')
# leaf_core_spei_mean_nN = subset(leaf_core_spei_mean, Ntrt_fac == '0')
# 
# # check
# nrow(leaf_core_spei_mean_yN)
# nrow(leaf_core_spei_mean_nN)
# 
# # deltaNper = ((leaf_core_spei_mean_yN$leaf_pct_N_Mean - leaf_core_spei_mean_nN$leaf_pct_N_Mean) / 
# #                leaf_core_spei_mean_nN$leaf_pct_N_Mean) * 100
# 
# deltaNarea = ((leaf_core_spei_mean_yN$Narea_Mean - leaf_core_spei_mean_nN$Narea_Mean) / 
#                 leaf_core_spei_mean_nN$Narea_Mean) * 100
# 
# # summary(lm(deltaNper ~ leaf_core_spei_mean_nN$p_pet_Mean)) # nothing
# anova(lm(deltaNarea ~ leaf_core_spei_mean_nN$p_pet_Mean * leaf_core_spei_mean_nN$par_rat_Mean)) # nothing
# 
# 
# #########################################
# # A2: nope, not really! basically no change regardless of site/climate
# #########################################

# deltaNarea_line = summary(lm_deltaNarea)$coefficients[1, 1] + summary(lm_deltaNarea)$coefficients[2, 1] * seq(0, 1, 0.1)
# 
# par(mfrow = c(1, 1), oma = c(4, 4, 1, 1))
# palette = colorRampPalette(brewer.pal(9,'Blues'))
# p_pet_cols <- palette(9)[as.numeric(cut(leaf_core_spei_mean_nN$p_pet_Mean,breaks = 9))]
# plot(deltaNarea ~ leaf_core_spei_mean_nN$par_rat_Mean, pch =21, cex = 3, bg = 'blue',  yaxt = 'n', xaxt = 'n', ylab = '', xlab = '', xlim = c(0, 1), ylim = c(-100, 300))
# axis(1, seq(0, 1, 0.2), cex.axis = 2)
# axis(2, seq(-100, 300, 100), cex.axis = 2, las = 1)
# mtext(side = 1, 'Light interception', line = 4, cex = 2)
# mtext(side = 2, '∆N per leaf area (%)', line = 5, cex = 2)
# lines(deltaNarea_line ~ seq(0, 1, 0.1), lwd = 4)
# text(0.5, 290, 'p = 0.35; mean = +12%', cex = 1.8)

#########################################
# Q3: does plot N differ among treatments?
#########################################

plot = read.csv('../Data/full-biomass-nutrients-06-December-2017.csv')

plot$Ntrt = 0
plot$Ntrt[plot$trt == 'N' | plot$trt == 'NP' | plot$trt == 'NK' | plot$trt == 'NPK' | plot$trt == 'NPK+Fence'] = 1
plot$Ptrt = 0
plot$Ptrt[plot$trt == 'P' | plot$trt == 'NP' | plot$trt == 'PK' | plot$trt == 'NPK' | plot$trt == 'NPK+Fence'] = 1
plot$Ktrt = 0
plot$Ktrt[plot$trt == 'K' | plot$trt == 'NK' | plot$trt == 'PK' | plot$trt == 'NPK' | plot$trt == 'NPK+Fence'] = 1

plot$pct_N_num = as.numeric(as.character(plot$pct_N))
plot$N_tot = as.numeric(as.character(plot$pct_N)) * as.numeric(as.character(plot$mass))

plot$Ntrt_fac = as.factor(plot$Ntrt)
plot$Ptrt_fac = as.factor(plot$Ptrt)
plot$Ktrt_fac = as.factor(plot$Ktrt)

plot$Nfix = 'no'
plot$Nfix[plot$category == "LEGUME"] = 'yes'
plot$Nfix[plot$category == "BRYOPHYTE"] = 'bryo'

# plot_pct_N_num_lmer = lmer(pct_N_num ~ Ntrt_fac * Ptrt_fac * Ktrt_fac * Nfix + (1|site_code), data = subset(plot, live == 1 & category != 'ANNUAL' & category != 'PERENNIAL' & category != 'LIVE'))
# Anova(plot_pct_N_num_lmer)
# plot(resid(plot_pct_N_num_lmer) ~ fitted(plot_pct_N_num_lmer))
# cld(emmeans(plot_pct_N_num_lmer, ~Ntrt_fac * Nfix)) # 16% increase in non-N fixing non-bryophytes
# cld(emmeans(plot_pct_N_num_lmer, ~Ntrt_fac))

# plot_N_tot_lmer = lmer(log(N_tot) ~ Ntrt_fac * Ptrt_fac * Ktrt_fac * Nfix + (1|site_code), data = subset(plot, live == 1 & category != 'ANNUAL' & category != 'PERENNIAL' & category != 'LIVE'))
# Anova(plot_N_tot_lmer)
# plot(resid(plot_N_tot_lmer) ~ fitted(plot_N_tot_lmer))
# cld(emmeans(plot_N_tot_lmer, ~Ntrt_fac * Nfix)) # 9% increase in non-N fixing non-bryophytes
# cld(emmeans(plot_N_tot_lmer, ~Ntrt_fac))
# 
# plot_N_tot_plot <- ggplot(subset(plot, live == 1 & category != 'ANNUAL' & category != 'PERENNIAL' & category != 'LIVE' & Nfix == 'no'), xlab = c('-N', '+N'), aes(x = Ntrt_fac, y = log(N_tot))) +
# 	geom_boxplot()
# 
# plot_N_tot_plot + scale_x_discrete(labels = c('-N', '+N'))
# 
# plot(log(N_tot) ~ Ntrt_fac,  data = subset(plot, live == 1 & category != 'ANNUAL' & category != 'PERENNIAL' & category != 'LIVE' & Nfix == 'no'), xlab = c('-N', '+N'))

#########################################
# A3: does plot N differ among treatments? yes, it seems that the N treatment increases plot %N by ~ 16% in the normal species. No change in N fixers and bryophytes. Interestingly, the total plot N only goes up by ~9% in the "normal" species. This is still quite a bit bigger than the change in per leaf area N!
#########################################


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
# Q4: can we predict Narea?
#########################################

# does N area increase with abs(latitude)?

plot(leaf_core_spei$Narea ~ abs(leaf_core_spei$latitude))
summary(lm(log(leaf_core_spei$Narea) ~ abs(leaf_core_spei$latitude))) # increases with latitude (temperature effect)

# what about with climate?

# read in gs climate (>0°C)
tmp_globe = read.csv('climate_gridded/cru_tmp_climExtract_growingseason_globe.csv')
par_globe = read.csv('climate_gridded/cru_par_climExtract_growingseason_globe.csv')
vpd_globe = read.csv('climate_gridded/cru_vpd_climExtract_growingseason_globe.csv')
z_globe =  read.csv('climate_gridded/z_globe.csv')

leaf_core_spei$tmp = NA
leaf_core_spei$par = NA
leaf_core_spei$vpd = NA
leaf_core_spei$z = NA

climate_df = c()
for (i in 1:nrow(leaf_core_spei)){
	
	currentLat = leaf_core_spei$latitude[i]
	currentLon = leaf_core_spei$longitude[i]
	
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

# plot(climate_df[,5] ~ climate_df[,7])
# plot(climate_df[,6] ~ climate_df[,8])

leaf_core_spei$tmp = climate_df[, 1]
leaf_core_spei$par = climate_df[, 2]
leaf_core_spei$vpd = climate_df[, 3]
leaf_core_spei$z = climate_df[, 4]

# write.csv(leaf_core_spei, "../Data/leaf_plus.csv")

leaf_core_spei_by_site_Ntrt = group_by(leaf_core_spei, site_code, Ntrt_fac)
leaf_core_spei_by_site_Ntrt_mean = summarise_all(leaf_core_spei_by_site_Ntrt, .funs = funs(Mean = mean(., na.rm = T), SD = sd(., na.rm = T)))

### ok let's try to compare predictions
leafNarea_lmer_clim = lm(log(Narea_Mean) ~ tmp_Mean + par_Mean + Ntrt_fac, data = leaf_core_spei_by_site_Ntrt_mean)
plot(resid(leafNarea_lmer_clim) ~ fitted(leafNarea_lmer_clim))
Anova(leafNarea_lmer_clim) # just temperature
emtrends(leafNarea_lmer_clim, ~1, var = 'tmp_Mean')
emmeans(leafNarea_lmer_clim, ~ tmp_Mean, at = list(tmp_Mean = 0))

# try to predict using optimal Vcmax


ggplot(leaf_core_spei_by_site_Ntrt_mean, aes(x = tmp_Mean, y = log(Narea_Mean))) + 
  geom_point() +
  stat_smooth(method = "lm", col = "red")



