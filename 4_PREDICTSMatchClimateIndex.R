##Read in packages
library(raster)
library(sp)
library(dismo)
library(rgdal)
library(predictsFunctions)
library(Rfast)
library(snow)
source("Functions.R")

dataDir <- "0_data/"
predictsDir <- "1_PreparePREDICTSData/"

outDir <- "4_PREDICTSMatchClimateIndex/"

##Path for monthly mean temperature from CRUv4.03
tmp.path <- "Data/cru_ts4.03.1901.2018.tmp.dat.nc"

##Path for mean monthly maximum from CRUv4.03
tmx.path <- "Data/cru_ts4.03.1901.2018.tmx.dat.nc"

output <- "Outputs/predicts_climate_info.rds"

##Read in average temperature data from CRU v4.03
tmp <- stack(paste0(dataDir,"cru_ts4.03.1901.2018.tmp.dat.nc"),varname = "tmp")

##Both cru data and predicts is in WGS84
wgs84 <- crs(tmp)

predicts_sites <- readRDS(paste0(predictsDir,"PREDICTSSiteData.rds"))
# Remove sites without coordinates
predicts_sites2 <- predicts_sites[!is.na(predicts_sites$Latitude), ]

# Create spatial map of PREDICTS sites
predicts_sp <- SpatialPointsDataFrame(
  coords = cbind(predicts_sites2$Longitude, predicts_sites2$Latitude), 
  data = predicts_sites2, proj4string = wgs84)

##Create a layer of mean monthly temperatures from 1901 to 1905
tmp1901_1905 <- tmp[[names(tmp)[1:60]]]

##Calculate the mean of 1901 to 1905 mean monthly temperatUres
tmp1901_1905mean <- raster::calc(tmp1901_1905, base::mean)

##Names of tmp layer, needed for subsettting
names <- names(tmp)

##Create a list of a all names of tmp layers, that will be used for matching later on
names_sub <- substr(names, 2, 8) 

##Spatial points for rasterizing
SP <- SpatialPoints(predicts_sp, proj4string=wgs84)

nCores <- parallel::detectCores()

st1 <- Sys.time()

cl <- snow::makeCluster(nCores-1)

snow::clusterExport(
  cl = cl,
  list = c('predicts_sp','names_sub','names',
           'tmp','SP','rasterize','crop','trim',
           'cellStats'),envir = environment())

temperatureVars <- data.frame(t(parSapply(
  cl = cl,X = (1:nrow(predicts_sp)),FUN = function(i){
  
  #Get end sample date for sample in predicts
  sampDate <- predicts_sp$Sample_end_latest[i]

  #Reformat date for string matching
  sampDate <- substr(sampDate,1, 7)
  sampDate <- gsub("-", ".", sampDate, fixed = TRUE)

  #Match date in predicts with month in CRU climate data
  month_match <- which(names_sub==sampDate)
  surrounding_months <- names[(month_match-11):(month_match)]

  #Create a layer for average temperature in the year preceding end sample date
  temp <- tmp[[surrounding_months]]

  ## Mask to improve speed

  mask <- trim(rasterize(SP[i, ], temp[[1]]))
  mapCrop <- crop(temp, mask)

  avg_temp <- mean(cellStats(mapCrop, stat = "mean"))
  sd_temp <- sd(cellStats(mapCrop, stat = "mean"))

  return(c(avg_temp=avg_temp,sd_temp=sd_temp))
  
})))

snow::stopCluster(cl)

st2 <- Sys.time()

print(st2 - st1)

#Calculate climate anomaly from sample compared with 1901 to 1905 mean
predicts_sp$avg_temp <- temperatureVars$avg_temp
predicts_sp$avg_temp_sd <- temperatureVars$sd_temp
predicts_sp$climate_anomaly <- temperatureVars$avg_temp-extract(tmp1901_1905mean, predicts_sp)

#####TMAX

##Read in tmax data
tmx <- stack(paste0(dataDir,"cru_ts4.03.1901.2018.tmx.dat.nc"))

#Create raster stack for 1901 to 1905
tmx1901_1905 <- tmx[[names(tmx)[1:60]]]
names(tmx1901_1905)

#Function to calculate mean of hottest quarter
max_quarter_fast <- function(x) {
  
  mean(Rfast::Sort(x, TRUE)[1:3])
  
}

#Hottest quarter for each year
tmx1901max <- calc(tmx1901_1905[[(names(tmx)[1:12])]], max_quarter_fast)
tmx1902max <- calc(tmx1901_1905[[(names(tmx)[13:24])]], max_quarter_fast)
tmx1903max <- calc(tmx1901_1905[[(names(tmx)[25:36])]], max_quarter_fast)
tmx1904max <- calc(tmx1901_1905[[(names(tmx)[37:48])]], max_quarter_fast)
tmx1905max <- calc(tmx1901_1905[[(names(tmx)[49:60])]], max_quarter_fast)

#Stack hottest quarters
tmx1901_1905max <- stack(tmx1901max, tmx1902max, tmx1903max, tmx1904max, tmx1905max)

#Mean of hottest quarter
tmx1901_1905max_mean <- calc(tmx1901_1905max, mean) 

#Baseline tmax
predicts_sp$tmax_baseline <- extract(tmx1901_1905max_mean, predicts_sp)

##names of tmax layer needed for subsetting
names_tmx <- names(tmx)

##For matching
names_sub_tmx <- substr(names_tmx, 2, 8)

st1 <- Sys.time()

cl <- snow::makeCluster(nCores-1)

snow::clusterExport(
  cl = cl,
  list = c('predicts_sp','names_sub_tmx','names_tmx','tmx',
           'trim','rasterize','crop','cellStats','SP',
           'max_quarter_fast'))

temperatureVarsTmax <- data.frame(t(parSapply(
  cl = cl,X = 1:nrow(predicts_sp),FUN = function(i){
    
    sampDate <- predicts_sp$Sample_end_latest[i]
    sampDate<- substr(sampDate,1, 7)
    sampDate <- gsub("-", ".", sampDate, fixed = TRUE)
    
    month_match <- which(names_sub_tmx==sampDate)
    surrounding_months <- names_tmx[(month_match-11):(month_match)]
    
    max_temp <- tmx[[surrounding_months]]
    
    ## Mask to improve speed
    
    mask <- trim(rasterize(SP[i, ], max_temp[[1]]))
    crop <- crop(max_temp, mask)
    
    tmax <- max(cellStats(crop, stat = "mean"))
    tmax_quarter <- max_quarter_fast(cellStats(crop, stat = "mean"))
    
    return(c(tmax=tmax,tmax_quarter=tmax_quarter))
    
    
  })))

snow::stopCluster(cl)

st2 <- Sys.time()

print(st2 - st1)

predicts_sp$tmax <- temperatureVarsTmax$tmax
predicts_sp$tmax_quarter <- temperatureVarsTmax$tmax_quarter

predicts_sp$tmax_anomaly <- predicts_sp$tmax - predicts_sp$tmax_baseline
predicts_sp$tmax_quarter_anomaly <- predicts_sp$tmax_quarter - predicts_sp$tmax_baseline

####### Historic variability in temperatures (Standard deviation of 1901 to 1905 Mean Monthly T )

tmp1901_1905sd <- calc(tmp1901_1905, stats::sd)
predicts_sp$historic_sd <- extract(tmp1901_1905sd, predicts_sp)
tmx1901_1905_sd <- calc(tmx1901_1905, stats::sd)
predicts_sp$historic_sd_tmax <- extract(tmx1901_1905_sd, predicts_sp)

##Save results

saveRDS(object = predicts_sp, file = paste0(outDir,"PREDICTSSitesWithClimateData.rds"))




