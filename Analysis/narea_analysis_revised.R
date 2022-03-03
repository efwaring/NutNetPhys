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
leaf = subset(leaf_raw, chi > 0 & chi < 1)

leaf$Ntrt_fac <- as.factor(leaf$Ntrt)
leaf$Ptrt_fac <- as.factor(leaf$Ptrt)
leaf$Ktrt_fac <- as.factor(leaf$Ktrt)
leaf$block_fac <- as.factor(leaf$block)

# leaf$tmp_scaled <- leaf$tmp - 25
# 
# leaf$pft <- leaf$functional_group
# leaf$pft[leaf$pft == 'NULL'] <- NA
# leaf$grass <- 'no'
# leaf$grass[leaf$pft == 'GRASS'] <- 'yes'
leaf$fence <- 'no'
leaf$fence[leaf$trt == 'Fence' | leaf$trt == 'NPK+Fence'] <- 'yes'

leaf$logpar <- log(leaf$par)
leaf$logpar_per_leaf_area <- log(leaf$par_per_leaf_area)
leaf$logvpd <- log(leaf$vpd)
leaf$loglma <- log(leaf$lma)

leaf$nstar <- calc_nstar(leaf$tmp, leaf$z)
leaf$gammastar <- calc_gammastar_pa(leaf$tmp, leaf$z)
leaf$km <- calc_km_pa(leaf$tmp, leaf$z)
## need to calculate adjusted vpd and ca!!
leaf$ca <- calc_
leaf$beta <- 1.6 * leaf$nstar * ((leaf$vpd * 1000) / ((1 - leaf$chi)^2)) * (leaf$chi - (leaf$gammastar))



