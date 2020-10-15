##%######################################################%##
#                                                          #
####                    Model Checks                    ####
#                                                          #
##%######################################################%##

# run checks for models used in paper

rm(list = ls())

# load libraries


# directories
LUclimmod <- "6_RunLUClimateModels/"
nathabmod <- "7_RunLUClimateNHModels/"


# read in model files
load(paste0(LUclimmod, "MeanAnomalyModelAbund.rdata"))
load(paste0(LUclimmod, "MeanAnomalyModelRich.rdata"))
load(paste0(LUclimmod, "MaxAnomalyModelAbund.rdata"))
load(paste0(LUclimmod, "MaxAnomalyModelRich.rdata"))

load(paste0(nathabmod, "MeanAnomalyModelAbun_NH.rdata"))
load(paste0(nathabmod, "RichMeanAnomalyModel_NH.rdata"))
load(paste0(nathabmod, "AbundMaxAnomalyModel_NH.rdata"))
load(paste0(nathabmod, "RichMaxAnomalyModel_NH.rdata"))


# run checks

## 1. Checking the fitted vs residuals relationship



## 2. Normality of Residuals



## 3. Check for spatial autocorrelation
