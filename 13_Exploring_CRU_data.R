##%######################################################%##
#                                                          #
####        13. Exploring data bias in CRU data         ####
#                                                          #
##%######################################################%##

# In this script, I will take a look at the data coverage of the CRU data
# using the associated stn data 
# (https://crudata.uea.ac.uk/cru/data/hrg/cru_ts_4.03/cruts.1905011326.v4.03/tmn/)



rm(list = ls())

# load libraries
library(raster)

# organise directories
datadir <- "0_data"
outdir <- "13_Exploring_CRU_data"
dir.create(outdir)

# load in the stn data
stn <- stack(paste0(datadir,"/cru_ts4.03.1901.2018.tmp.dat.nc"),varname = "stn")

# take the baseline years and take a look
stn_0130 <- stn[[1:360]]

stn_0130_mean <- mean(stn_0130)
plot(stn_0130_mean, main = "mean number of stations \nacross baseline 1901-1930")

stn_0130_range <- range(stn_0130)
plot(stn_0130_range)

# load in the predicts data and do some summaries for the sites
predictsDir <- "1_PreparePREDICTSData/"

# read in the predicts sites - insect subset
predicts_sites <- readRDS(paste0(predictsDir,"PREDICTSSiteData.rds"))

# Remove sites without coordinates
predicts_sites2 <- predicts_sites[!is.na(predicts_sites$Latitude), ] # this has already been done elsewhere

wgs84 <- crs(stn)
# Create spatial map of PREDICTS sites
predicts_sp <- SpatialPointsDataFrame(
  coords = cbind(predicts_sites2$Longitude, predicts_sites2$Latitude), 
  data = predicts_sites2, proj4string = wgs84)

# add the sites onto the plots
plot(stn_0130_mean, main = "mean number of stations \nacross baseline 1901-1930")
plot(predicts_sp, add = T)


# extract just from the mean number of station?
nstn <- extract(stn_0130_mean, predicts_sp)
  
nstn <- round(nstn, digits = 2)

table(nstn)
