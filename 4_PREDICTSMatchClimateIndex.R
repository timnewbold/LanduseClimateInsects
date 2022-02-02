##%######################################################%##
#                                                          #
####      Match climate values with PREDICTS sites      ####
#                                                          #
##%######################################################%##

# This script organises the climate data associated with each PREDICTS site

# updated version with mean anomaly based on active months and edited max anomaly
# July 2021


# load libraries
library(raster)
library(sp)
library(dismo)
library(rgdal)
library(predictsFunctions)
library(Rfast)
library(snow)
source("Functions.R")

# directories
dataDir <- "0_data/"
predictsDir <- "1_PreparePREDICTSData/"
outDir <- "4_PREDICTSMatchClimateIndex/"

if(!dir.exists(outDir)) dir.create(outDir)


##Path for monthly mean temperature from CRUv4.03
tmp.path <- "Data/cru_ts4.03.1901.2018.tmp.dat.nc"

##Path for mean monthly maximum from CRUv4.03
tmx.path <- "Data/cru_ts4.03.1901.2018.tmx.dat.nc"

output <- "Outputs/predicts_climate_info.rds"

##Read in average temperature data from CRU v4.03
tmp <- stack(paste0(dataDir,"cru_ts4.03.1901.2018.tmp.dat.nc"),varname = "tmp")

##Both cru data and predicts is in WGS84
wgs84 <- crs(tmp)

# read in the predicts sites - insect subset
predicts_sites <- readRDS(paste0(predictsDir,"PREDICTSSiteData.rds"))

# Remove sites without coordinates
predicts_sites2 <- predicts_sites[!is.na(predicts_sites$Latitude), ]

# Create spatial map of PREDICTS sites
predicts_sp <- SpatialPointsDataFrame(
  coords = cbind(predicts_sites2$Longitude, predicts_sites2$Latitude), 
  data = predicts_sites2, proj4string = wgs84)

##Create a layer of mean monthly temperatures from 1901 to 1930
tmp1901_1930 <- tmp[[names(tmp)[1:360]]]

##Read in tmax data
tmx <- stack(paste0(dataDir,"cru_ts4.03.1901.2018.tmx.dat.nc"),varname = "tmx")

#Create raster stack for 1901 to 1930
tmx1901_1930 <- tmx[[names(tmx)[1:360]]]

##Names of tmp layer, needed for subsettting
names_tmp <- names(tmp)

##Create a list of a all names of tmp layers, that will be used for matching later on
names_sub <- substr(names_tmp, 2, 8) 

##Spatial points for rasterizing
SP <- SpatialPoints(predicts_sp, proj4string=wgs84)

### set the threshold temperature for insect activity
thresh <- 10


#### calculating the mean based anomaly for each site ####

# this now looks at the 5 years preceding each sample data
# assesses which months are above the temp threshold
# then uses these months to calculate the present and baseline temperature means
# and the baseline temperature sd. 

# Time difference of 7.857316 mins

nCores <- parallel::detectCores()

st1 <- Sys.time()

cl <- snow::makeCluster(nCores-1)

# export to clusters
snow::clusterExport(
  cl = cl,
  list = c('predicts_sp','names_sub','names_tmp', 'values', 'names', 'length', 'mean', 'sd',
           'tmp', 'tmx', 'SP','rasterize','crop','trim', 'grep', 'sapply', 'strsplit',
           'cellStats', 'thresh', 'tmp1901_1930', 'tmx1901_1930'),envir = environment())

temperatureVars <- data.frame(t(parSapply(
  cl = cl,X = (1:nrow(predicts_sp)),FUN = function(i){

    #for testing
    #temperatureVars <- NULL
    #for(i in 1538:nrow(predicts_sp)){
    #for(i in 1538:1580){
    #print(i)
    
    #Get end sample date for sample in predicts
    sampDate <- predicts_sp$Sample_end_latest[i]
    
    #Reformat date for string matching
    sampDate <- substr(sampDate,1, 7)
    sampDate <- gsub("-", ".", sampDate, fixed = TRUE)
    
    #Match date in predicts with month in CRU climate data
    month_match <- which(names_sub==sampDate)
    
    # edit: use months from 5 year pre-sample, rather than 1 year
    surrounding_months <- names_tmp[(month_match-59):(month_match)]
    
    #Create a layer for average temperature in the year preceding end sample date
    temp <- tmp[[surrounding_months]]
    temp_mx <- tmx[[surrounding_months]]
    
    ## Mask to improve speed
    mask <- trim(rasterize(SP[i, ], temp[[1]]))
    mapCrop <- crop(temp, mask)
    

    # there are instances where there are no months above the threshold and
    # other instances where points do not line up with the tmp layers
    # so this if statement is necessary to avoid errors in those instances.
    if(!length(names(mapCrop)[values(mapCrop) >= thresh]) == 0 & length(values(mapCrop)[!is.na(values(mapCrop))]) > 0 ){
      
      # Get the average temperature for each month across 5 years
      vals <- NULL
      
      # for each month, get the average temp over the 5 years
      for(j in 1:12){
        
        if(j < 10){ mon <- paste0(0, j) }else {mon <- j}
        
        monthmean <- values(mean(mapCrop[[grep(mon, sapply(strsplit(names(mapCrop), "[.]"), "[[", 2))  ]]))
        
        vals <- rbind(vals, c(mon, monthmean))
        
      }
      

      vals <- as.data.frame(vals)
      vals$V2 <- as.numeric(as.character(vals$V2))
      
      # which months are the 5 year average >= the threshold
      vals <- vals[vals$V2 >= thresh, ]
      
      
      # which are the good months
      months <- vals$V1
      
      # how many months are at or above the threshold?
      n_months <- length(months)
      
      # calculate the "present day" mean and sd
      avg_temp <- mean(vals$V2)
      #sd_temp <- sd(vals$V2) # don't actually use this
      
      
      
      # get 3 hottest month vals for those months that meet the threshold
      
      # loop through each year, get the 3 hottest months
      
      mask_mx <- trim(rasterize(SP[i, ], temp_mx[[1]]))
      mapCrop_mx <- crop(temp_mx, mask_mx)
      
      # for each year get the 3 hottest monts then take the mean of these
      yr1 <- mean(sort(values(mapCrop_mx[[1:12]]), decreasing = T)[1:3])
      yr2 <- mean(sort(values(mapCrop_mx[[13:24]]), decreasing = T)[1:3])
      yr3 <- mean(sort(values(mapCrop_mx[[25:36]]), decreasing = T)[1:3])
      yr4 <- mean(sort(values(mapCrop_mx[[37:48]]), decreasing = T)[1:3])
      yr5 <- mean(sort(values(mapCrop_mx[[49:60]]), decreasing = T)[1:3])
      
      # then take the mean across years to get one max value for the present
      # this is the mean for the "present day" max (mean of 3 hottest months)
      max_temp <- mean(yr1, yr2, yr3, yr4, yr5)
      

      ### now work out the baseline mean and sd for the active months and hottest months ###

      # get the values for that grid cell across all years
      baseline <- crop(tmp1901_1930, mask)
      
      # subset the baseline to just the required months
      baseline <-  baseline[[names(baseline)[sapply(strsplit(names(baseline), "[.]"), "[[", 2) %in% months]]]
      
      # get the mean and sd
      mean_baseline <- mean(values(baseline))
      sd_mean_baseline <- sd(values(baseline))
      
      
      
      ### now max temps, get 3 hottest years for each year ###
      baseline_max <- crop(tmx1901_1930, mask_mx)
      
      # for each year, get the 3 hottest years #
      # get the year names
      yrs <- unique(gsub("\\..*", "", names(baseline_max)))
      
      max_vals <- NULL
      
      for(yr in yrs){
        
        mx3 <- sort(values(baseline_max[[grep(yr, names(baseline_max))]]), decreasing = T)[1:3]
        
        max_vals <- c(max_vals, mx3)
        
      }
    
      # get the baseline values for this site
      mean_baseline_mx <- mean(max_vals)
      sd_baseline_mx <- sd(max_vals)
    
    
      ### now calc the anomaly for that site, using the site specific baselines ###
      
      Anom <- avg_temp - mean_baseline
      StdAnom <-  Anom/sd_mean_baseline
      
      tmax_anomaly <- max_temp - mean_baseline_mx
      StdTmaxAnomaly <- tmax_anomaly/sd_baseline_mx
      
      
      return(c(avg_temp=avg_temp, Anom = Anom, StdAnom = StdAnom, n_months = n_months,
               max_temp = max_temp, tmax_anomaly = tmax_anomaly, StdTmaxAnomaly = StdTmaxAnomaly))
      
      
      #temperatureVars <- rbind(temperatureVars, c(avg_temp=avg_temp, Anom = Anom, StdAnom = StdAnom, n_months = n_months,
      #max_temp = max_temp, tmax_anomaly = tmax_anomaly, StdTmaxAnomaly = StdTmaxAnomaly))
      
    }else{
      avg_temp <- NA
      Anom <- NA
      StdAnom <- NA
      n_months = 0
      max_temp <- NA
      tmax_anomaly <- NA
      StdTmaxAnomaly <- NA
      
      #temperatureVars <- rbind(temperatureVars, c(avg_temp=avg_temp, Anom = Anom, StdAnom = StdAnom, n_months = n_months,
      #max_temp = max_temp, tmax_anomaly = tmax_anomaly, StdTmaxAnomaly = StdTmaxAnomaly))
      
      return(c(avg_temp=avg_temp, Anom = Anom, StdAnom = StdAnom, n_months = n_months,
               max_temp = max_temp, tmax_anomaly = tmax_anomaly, StdTmaxAnomaly = StdTmaxAnomaly))      
    }
    
    
  }
  
)))


snow::stopCluster(cl)

st2 <- Sys.time()

print(st2 - st1) # Time difference of 7.857316 mins


# organise the anomaly info along with the predicts data
temperatureVars <- as.data.frame(temperatureVars)


## add new values in temperatureVars into predicts dataset
predicts_sp$avg_temp <- temperatureVars$avg_temp
predicts_sp$TmeanAnomaly <- temperatureVars$Anom
predicts_sp$StdTmeanAnomaly <- temperatureVars$StdAnom
predicts_sp$n_months <- temperatureVars$n_months
predicts_sp$max_temp <- temperatureVars$max_temp
predicts_sp$TmaxAnomaly <- temperatureVars$tmax_anomaly
predicts_sp$StdTmaxAnomaly <- temperatureVars$StdTmaxAnomaly


# save
saveRDS(object = predicts_sp, file = paste0(outDir,"PREDICTSSitesWithClimateData_update.rds"))




