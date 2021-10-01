##%######################################################%##
#                                                          #
####     Additional tests - mcmc version of models      ####
#                                                          #
##%######################################################%##

# in this script I will rerun the models using an mcmc approach to help
# test whether the models are appropriately specified. 

rm(list = ls())

# load libraries
library(MCMCglmm)
library(parallel)
library(coda)
library(snow)

# sort directories
datadir <- "6_RunLUClimateModels/"
outdir <- "14_Additional_Tests/"

# read in predicts data
predictsSites <- readRDS(file = paste0(datadir,"PREDICTSSiteData.rds"))


# Specify the model

#### 1. abundance model, mean anomaly ####

model_data <- predictsSites[!is.na(predictsSites$LogAbund), ]
model_data <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ]


abund_mean_mod <- MCMCglmm(LogAbund ~ UI2 + poly(StdTmeanAnomalyRS,1) + UI2:poly(StdTmeanAnomalyRS,1), 
         random = ~ SS + SSB,
         family = "gaussian",
         data = model_data,
         verbose = T,
         nitt = 23000,
         thin = 20)

# take a look at the summary
summary.MCMCglmm(abund_mean_mod)
plot(abund_mean_mod) # this plots the diagnostic plots

save(abund_mean_mod, file = paste0(outdir, "/MeanAnomalyModelAbund_mcmc.rdata"))
# load(file = paste0(outdir, "/MeanAnomalyModelAbund_mcmc.rdata"))



# try using snow
# run model for 4 chains to get diagnostics
st1 <- Sys.time()

cl <- snow::makeCluster(4)

snow::clusterExport(
  cl = cl,
  list = c('model_data', 'MCMCglmm'),envir = environment())


ab_mean_mod <- parLapply(cl, 1:4, function(i){
  MCMCglmm(LogAbund ~ UI2 + poly(StdTmeanAnomalyRS,1) + UI2:poly(StdTmeanAnomalyRS,1), 
           random = ~ SS + SSB,
           family = "gaussian",
           data = model_data,
           verbose = T,
           nitt = 23000,
           thin = 20)
}) 


snow::stopCluster(cl)

st2 <- Sys.time()

print(st2 - st1) # Time difference of 48.81089 secs

# combines outputs from each chain for use with gelman functions next
ab1 <- lapply(ab_mean_mod, function(m) m$Sol)
ab1 <- do.call(mcmc.list, ab1)

# model checks
par(mfrow=c(4,2), mar=c(2,2,1,2))
gelman.plot(ab1, auto.layout=F) # looking at when chains converge

# table of Rhat values
gelman.diag(ab1) # rhat values

par(mfrow=c(8,2), mar=c(2, 1, 1, 1))
plot(ab1, ask=F, auto.layout=F) # more trace plots but with multiple chains, check for mixing


# load REML version of the model and compare outputs
load(paste0(datadir, "/MeanAnomalyModelAbund.rdata"))

summary(MeanAnomalyModelAbund$model)
summary(abund_mean_mod)
# values are very similar

#### 2. Richness, mean anomaly ####

model_data <- predictsSites[!is.na(predictsSites$StdTmeanAnomaly), ]


rich_mean_model <- MCMCglmm(Species_richness ~ UI2 + poly(StdTmeanAnomalyRS,1) + UI2:poly(StdTmeanAnomalyRS,1), 
                   random = ~ SS + SSB+SSBS, # removed SSBS, much better diagnostics
                   family = "poisson", 
                   data = model_data, 
                   verbose = T, 
                   nitt = 23000, 
                   thin = 20)

# take a look at the model summary
summary.MCMCglmm(rich_mean_model)
plot(rich_mean_model) # this plots the diagnostic plots including random effects

# trace plots look better when SSBS is removed

# save model output
save(rich_mean_model, file = paste0(outdir, "/MeanAnomalyModelRich_mcmc.rdata"))
# load(file = paste0(outdir, "/MeanAnomalyModelRich_mcmc.rdata"))



st1 <- Sys.time()

cl <- snow::makeCluster(4)

snow::clusterExport(
  cl = cl,
  list = c('model_data', 'MCMCglmm'),envir = environment())


sr_mean_mod <- parLapply(cl, 1:4, function(i){
                MCMCglmm(Species_richness ~ UI2 + poly(StdTmeanAnomalyRS,1) + UI2:poly(StdTmeanAnomalyRS,1), 
                        random = ~ SS + SSB, # removed SSBS, much better diagnostics
                        family = "poisson", 
                        data = model_data, 
                        verbose = T, 
                        nitt = 23000, 
                        thin = 20)})

snow::stopCluster(cl)

st2 <- Sys.time()

print(st2 - st1) # Time difference of 2.690455 mins

# combines outputs from each chain for use with gelman functions next
sr1 <- lapply(sr_mean_mod, function(m) m$Sol)
sr1 <- do.call(mcmc.list, sr1)

# model checks
par(mfrow=c(4,2), mar=c(2,2,1,2))
gelman.plot(sr1, auto.layout=F) # looking at when chains converge

# table of Rhat values
gelman.diag(sr1) # rhat values

par(mfrow=c(8,2), mar=c(2, 1, 1, 1))
plot(sr1, ask=F, auto.layout=F) # more trace plots but with multiple chains, check for mixing


# load REML version of the model and compare outputs
load(paste0(datadir, "/MeanAnomalyModelRich.rdata"))

summary(MeanAnomalyModelRich$model)
summary.MCMCglmm(rich_mean_model)
# some slight differences between some variables


#### 3. Abundance, max anomaly ####

model_data <- predictsSites[!is.na(predictsSites$LogAbund), ]
model_data <- predictsSites[!is.na(predictsSites$StdTmaxAnomalyRS), ]


Abun_max_model <- MCMCglmm(LogAbund ~ UI2 + poly(StdTmaxAnomalyRS,1) + UI2:poly(StdTmaxAnomalyRS,1), 
                           random = ~ SS + SSB, 
                           family = "gaussian", 
                           data = model_data, 
                           verbose = T, 
                           nitt = 23000, 
                           thin = 20)

# take a look at the summary
summary.MCMCglmm(Abun_max_model)
plot(Abun_max_model) # this plots the diagnostic plots

save(Abun_max_model, file = paste0(outdir, "/MaxAnomalyModelAbund_mcmc.rdata"))
# load(file = paste0(outdir, "/MaxAnomalyModelAbund_mcmc.rdata"))


# run in parallel for 4 chains
st1 <- Sys.time()

cl <- snow::makeCluster(4)

snow::clusterExport(
  cl = cl,
  list = c('model_data', 'MCMCglmm'),envir = environment())


ab_max_mod <- parLapply(cl, 1:4, function(i){
MCMCglmm(LogAbund ~ UI2 + poly(StdTmaxAnomalyRS,1) + UI2:poly(StdTmaxAnomalyRS,1), 
                        random = ~ SS + SSB, 
                        family = "gaussian", 
                        data = model_data, 
                        verbose = T, 
                        nitt = 23000, 
                        thin = 20)})

snow::stopCluster(cl)

st2 <- Sys.time()

print(st2 - st1) # Time difference of 28.93384 secs

# combines outputs from each chain for use with gelman functions next
ab2 <- lapply(ab_max_mod, function(m) m$Sol)
ab2 <- do.call(mcmc.list, ab2)

# model checks
par(mfrow=c(4,2), mar=c(2,2,1,2))
gelman.plot(ab2, auto.layout=F) # looking at when chains converge

# table of Rhat values
gelman.diag(ab2) # rhat values

par(mfrow=c(8,2), mar=c(2, 1, 1, 1))
plot(ab2, ask=F, auto.layout=F) # more trace plots but with multiple chains, check for mixing



# load REML version of the model and compare outputs
load(paste0(datadir, "/MaxAnomalyModelAbund.rdata"))

summary(MaxAnomalyModelAbund$model)
summary.MCMCglmm(Abun_max_model)



#### 4. Richness, max anomaly ####

model_data <- predictsSites[!is.na(predictsSites$StdTmaxAnomaly), ]


rich_max_model <- MCMCglmm(Species_richness ~ UI2 + poly(StdTmaxAnomalyRS,1) + UI2:poly(StdTmaxAnomalyRS,1), 
                           random = ~ SS + SSB + SSBS, # removed SSBS for better trace plots
                           family = "poisson", 
                           data = model_data, 
                           verbose = T, 
                           nitt = 23000, 
                           thin = 20)


# take a look at the summary
summary.MCMCglmm(rich_max_model)
plot(rich_max_model) # this plots the diagnostic plots


# save model output
save(rich_max_model, file = paste0(outdir, "/MaxAnomalyModelRich_mcmc.rdata"))
# load(file = paste0(outdir, "/MaxAnomalyModelRich_mcmc.rdata"))


# run in parallel for 4 chains
st1 <- Sys.time()

cl <- snow::makeCluster(4)

snow::clusterExport(
  cl = cl,
  list = c('model_data', 'MCMCglmm'),envir = environment())


sr_max_mod <- parLapply(cl, 1:4, function(i){
MCMCglmm(Species_richness ~ UI2 + poly(StdTmaxAnomalyRS,1) + UI2:poly(StdTmaxAnomalyRS,1), 
                        random = ~ SS + SSB, 
                        family = "poisson", 
                        data = model_data, 
                        verbose = T, 
                        nitt = 23000, 
                        thin = 20)})

snow::stopCluster(cl)

st2 <- Sys.time()

print(st2 - st1) # Time difference of 

# combines outputs from each chain for use with gelman functions next
sr2 <- lapply(sr_max_mod, function(m) m$Sol)
sr2 <- do.call(mcmc.list, sr2)

# model checks
par(mfrow=c(4,2), mar=c(2,2,1,2))
gelman.plot(sr2, auto.layout=F) # looking at when chains converge

# table of Rhat values
gelman.diag(sr2) # rhat values

par(mfrow=c(8,2), mar=c(2, 1, 1, 1))
plot(sr2, ask=F, auto.layout=F) # more trace plots but with multiple chains, check for mixing


# load REML version of the model and compare outputs
load(paste0(datadir, "/MaxAnomalyModelRich.rdata"))

summary(MaxAnomalyModelRich$model)
summary.MCMCglmm(rich_max_model)




#### save output tables ####

# trying out functions from here:
# https://gkhajduk.github.io/2017-10-25-cleanMCMCglmm/

library(dplyr)

source("clean.MCMC.R")

# create a list of models to get the outputs for
mod_list <- list(abund_mean_mod, rich_mean_model, Abun_max_model, rich_max_model)

# set a list of names for the models to be identified in the table
mod_list_names <- list("Abun_mean_anomaly", "Rich_mean_anomaly", "Abun_max_anomaly", "Rich_max_anomaly")

# apply the clean.mcmc function 
readyList <- mapply(cbind, lapply(mod_list, clean.MCMC), "modelName" = mod_list_names, SIMPLIFY = F)

# turn list of data frames into one dataframe
mcmcOutputs <- as.data.frame(do.call(rbind, readyList), stringsAsFactors = FALSE)

# now can use another package to create a nice summary table
library(stargazer)
stargazer(mcmcOutputs, type = "text", summary = FALSE, 
          out = paste0(outdir, "MCMC_output_tables.htm"))

