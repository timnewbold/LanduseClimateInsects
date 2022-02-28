##%######################################################%##
#                                                          #
####      Organise Natural habitat info for sites       ####
#                                                          #
##%######################################################%##


# This script calculates the proportion of NH surrounding each site 
# in PREDICTS subset at multiple spatial scales: 1,3,5 and 10km.


# directories
dataDir <- "0_data/"
predictsDir <- "4_PREDICTSMatchClimateIndex/"
outDir <- "5_PREDICTSMatchPropNatHab/"

if(!dir.exists(outDir)) dir.create(outDir)

sink(paste0(outDir,"log.txt"))

t.start <- Sys.time()

print(t.start)

# load libraries
library(raster)
library(sp)
library(predictsFunctions)
library(snow)

# required crs
wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"

# read in the predicts data with climate info
predicts_sp <- readRDS(paste0(predictsDir,"PREDICTSSitesWithClimateData_update.rds"))

# load the NH raster layers from the Hoskins et al dataset
PV <- raster(paste0(
  dataDir,"PRI_1km_2005_0ice.bil"), 
  crs=wgs84)
SV <- raster(paste0(
  dataDir,"SEC_1km_2005_0ice.bil"), 
  crs=wgs84)

# run in parallel
nCores <- parallel::detectCores()

st1 <- Sys.time()

cl <- snow::makeCluster(nCores-1)

# export data to cores
snow::clusterExport(
  cl = cl,
  list = c('predicts_sp','PV','SV','buffer','crop','cellStats'),envir = environment())

# Time difference of 31.47568 mins

# for each site, get the proportion of NH surrounding the site
natHabitat <- data.frame(t(parSapply(cl = cl,X = (1:nrow(predicts_sp)),FUN = function(i){
  cat(paste0("Processing site ",i," of ",nrow(predicts_sp),"\r"))
  
  buff1 <- buffer(predicts_sp[i, ], width=1000)
  PV_crop1 <- crop(PV, buff1)
  PV_mean1 <- cellStats(PV_crop1, stat="mean")
  SV_crop1 <- crop(SV, buff1)
  SV_mean1 <- cellStats(SV_crop1, stat="mean")
  # predicts_sp$PV_1000[i] <- PV_mean1
  # predicts_sp$SV_1000[i] <-  SV_mean1
  
  
  buff3 <- buffer(predicts_sp[i, ], width=3000)
  PV_crop3 <- crop(PV, buff3)
  PV_mean3 <- cellStats(PV_crop3, stat="mean")
  SV_crop3 <- crop(SV, buff3)
  SV_mean3 <- cellStats(SV_crop3, stat="mean")
  # predicts_sp$PV_3000[i] <- PV_mean3
  # predicts_sp$SV_3000[i] <-  SV_mean3
  
  
  buff5 <- buffer(predicts_sp[i, ], width=5000)
  PV_crop5 <- crop(PV, buff5)
  PV_mean5 <- cellStats(PV_crop5, stat="mean")
  SV_crop5 <- crop(SV, buff5)
  SV_mean5 <- cellStats(SV_crop5, stat="mean")
  # predicts_sp$PV_5000[i] <- PV_mean5
  # predicts_sp$SV_5000[i] <-  SV_mean5
  
  
  buff10 <- buffer(predicts_sp[i, ], width=10000)
  PV_crop10 <- crop(PV, buff10)
  PV_mean10 <- cellStats(PV_crop10, stat="mean")
  SV_crop10 <- crop(SV, buff1)
  SV_mean10 <- cellStats(SV_crop10, stat="mean")
  # predicts_sp$PV_10000[i] <- PV_mean10
  # predicts_sp$SV_10000[i] <-  SV_mean10
  
  return(c(PV1=PV_mean1,SV1=SV_mean1,
           PV3=PV_mean3,SV3=SV_mean3,
           PV5=PV_mean5,SV5=SV_mean5,
           PV10=PV_mean10,SV10=SV_mean10))
  
})))

snow::stopCluster(cl)

st2 <- Sys.time()

print(st2 - st1) # Time difference of 8.622247 mins

predicts_sp$PV_1000 <- natHabitat$PV1
predicts_sp$SV_1000 <- natHabitat$SV1
predicts_sp$PV_3000 <- natHabitat$PV3
predicts_sp$SV_3000 <- natHabitat$SV3
predicts_sp$PV_5000 <- natHabitat$PV5
predicts_sp$SV_5000 <- natHabitat$SV5
predicts_sp$PV_10000 <- natHabitat$PV10
predicts_sp$SV_10000 <- natHabitat$SV10

saveRDS(object = predicts_sp,file = paste0(outDir,"PREDICTSSitesWithClimateAndNatHab.rds"))



t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()

