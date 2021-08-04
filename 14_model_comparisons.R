##%######################################################%##
#                                                          #
####       testing variations on the main models        ####
#                                                          #
##%######################################################%##

# in this script I will test models using average temperature, the anomaly 
# (difference in temperature) and the standardused anomaly and compare them.

rm(list = ls())

# load libraries
library(devtools)
library(StatisticalModels)
library(predictsFunctions)
source("Functions.R")
library(sjPlot)

# directories 
predictsDataDir <- "6_RunLUClimateModels/"
outDir <- "14_Additional_Tests/"

# load in the data
predictsSites <- readRDS(file = paste0(predictsDataDir,"PREDICTSSiteData.rds"))


# rescale the other variables
predictsSites$avg_tempRS <- StdCenterPredictor(predictsSites$avg_temp)
predictsSites$TmeanAnomalyRS <- StdCenterPredictor(predictsSites$TmeanAnomaly)
  

#### first run the abundance models ####
# running the model selection process


# 1. Abundance, average temperature
avgtempModelAbund <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(avg_tempRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(avg_tempRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(avgtempModelAbund, file = paste0(outDir, "avgtempModelAbund.rdata"))



# 2. Abundance,  mean anomaly
AnomalyModelAbund <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(TmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(TmeanAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(AnomalyModelAbund, file = paste0(outDir, "AnomalyModelAbund.rdata"))



# 3. Abundance, standardised mean anomaly
StdAnomalyModelAbund <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(StdAnomalyModelAbund, file = paste0(outDir, "StdAnomalyModelAbund.rdata"))


### compare AIC values ###

print(AIC(avgtempModelAbund$model,AnomalyModelAbund$model,StdAnomalyModelAbund$model))


#df      AIC
#avgtempModelAbund$model    11 15107.56
#AnomalyModelAbund$model    11 15122.35
#StdAnomalyModelAbund$model 11 15092.59


#### species richness models ####

# 1. Richness, average temp
avgtempModelRich <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
                                    fitFamily = "poisson",fixedFactors = "UI2",
                                    fixedTerms = list(avg_tempRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(avg_tempRS,1)"),
                                    saveVars = c("Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(avgtempModelRich, file = paste0(outDir, "avgtempModelRich.rdata"))

# 2. Richness, mean anomaly
AnomalyModelRich <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
                                    fitFamily = "poisson",fixedFactors = "UI2",
                                    fixedTerms = list(TmeanAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(TmeanAnomalyRS,1)"),
                                    saveVars = c("Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(AnomalyModelRich, file = paste0(outDir, "AnomalyModelRich.rdata"))

# 3. Richness, Std mean anomaly
StdAnomalyModelRich <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
                                    fitFamily = "poisson",fixedFactors = "UI2",
                                    fixedTerms = list(StdTmeanAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                    saveVars = c("Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(StdAnomalyModelRich, file = paste0(outDir, "StdAnomalyModelRich.rdata"))


### compare AIC values ###

print(AIC(avgtempModelRich$model,AnomalyModelRich$model,StdAnomalyModelRich$model))

# 
# df      AIC
# avgtempModelRich$model    11 33352.27
# AnomalyModelRich$model    11 33364.64
# StdAnomalyModelRich$model 11 33277.79