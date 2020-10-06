##%######################################################%##
#                                                          #
####            Testing different baselines             ####
#                                                          #
##%######################################################%##

# In this script, different baselines for the climate anomaly are tested.
# In the main analysis the baseline years are 1901-1905. 

# alternatives to test:
# longer time period for baseline 1901-1920
# different time period 1920-1925, 1940-1945

# calculate the SCA and SCAM metrics using these baselines
# rerun the model with LU and SCA/SCAM interaction 
# look at Figure 1 equivalents. 

rm(list = ls())

# load libraries
library(raster)
library(sp)
library(dismo)
library(rgdal)
library(RColorBrewer)
library(ncdf4)
library(rasterVis)
library(gridExtra)
library(cowplot)
library(ggplot2)
library(viridis)

# directories
dataDir <- "0_data/"
outDir <- "10_SCA_Baseline_testing/"
dir.create(outDir)


# extract data for dif baseline periods

# load in the mean temperature data from CRU
tmp <- stack(paste0(dataDir,"cru_ts4.03.1901.2018.tmp.dat.nc"),varname="tmp")

# take names of values for 1901 to 1905
tmp1901_1905 <- tmp[[names(tmp)[1:60]]]

# calculate the mean and sd of the baseline values
tmp1901_1905mean <- calc(tmp1901_1905, base::mean)
tmp1901_1905sd <- calc(tmp1901_1905, stats::sd)





