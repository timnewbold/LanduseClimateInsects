##Read in packages
library(raster)
library(sp)
library(dismo)
library(rgdal)
library(predictsFunctions)
library(Rfast)
source("r/functions.r")


##The following file paths are required for this script to run:
##Path for predicts database
predicts.path <-"Data/database.rds"

##Path for monthly mean temperature from CRUv4.03
tmp.path <- "Data/cru_ts4.03.1901.2018.tmp.dat.nc"

##Path for mean monthly maximum from CRUv4.03
tmx.path <- "Data/cru_ts4.03.1901.2018.tmx.dat.nc"

output <- "Outputs/predicts_climate_info.rds"

##Read in average temperature data from CRU v4.03
tmp <- stack(tmp.path)

##Both cru data and predicts is in WGS84
wgs84 <- crs(tmp)


#Read in predicts
predicts <- readRDS(predicts.path)

##Merge sites and create Spatial points dataframe
predicts_sites <- create_sites(predicts)
predicts_sites2 <- predicts_sites[!is.na(predicts_sites$Latitude), ]
predicts_sp <- SpatialPointsDataFrame(coords = cbind(predicts_sites2$Longitude, predicts_sites2$Latitude), data = predicts_sites2, proj4string = wgs84)




##Create a layer of mean monthly temperatures from 1901 to 1905
tmp1901_1905 <- tmp[[names(tmp)[1:60]]]
names(tmp1901_1905)

##Calculate the mean of 1901 to 1905 mean monthly temperatUres
tmp1901_1905mean <- raster::calc(tmp1901_1905, base::mean)


##Names of tmp layer, needed for subsettting
names <- names(tmp)

##Create a list of a all names of tmp layers, that will be used for matching later on
names_sub <- substr(names, 2, 8) 
names_sub
##Empty vector for storing matching months
month_match <-vector()

predicts_sp$climate_anomaly <- rep(0, 1)
predicts_sp$avg_temp <- rep(0, 1)
predicts_sp$avg_temp_sd <- rep(0, 1)

##Spatial points for rasterizing
SP <- SpatialPoints(predicts_sp, proj4string=wgs84)
###For each site in predicts find climate anomaly, match with CRU climate data
system.time(for(i in 1:nrow(predicts_sp)) {
  #Get end sample date for sample in predicts
  date <- predicts_sp$Sample_end_latest[i]
  
  #Reformat date for string matching
  date<- substr(date,1, 7)
  date <- gsub("-", ".", date, fixed = TRUE)
  
  #Match date in predicts with month in CRU climate data
  month_match <- which(names_sub==date)
  surrounding_months <- names[(month_match-11):(month_match)]
  
  #Create a layer for average temperature in the year preceding end sample date
  temp <- tmp[[surrounding_months]]
  
  ## Mask to improve speed
  
  mask <- trim(rasterize(SP[i, ], temp[[1]]))
  crop <- crop(temp, mask)
  
  
  #Extract average temperature 
  predicts_sp$avg_temp[i] <-  mean(cellStats(crop, stat = "mean"))
  predicts_sp$avg_temp_sd[i] <- sd(cellStats(crop, stat = "mean"))
  
})



#Calculate climate anomaly from sample compared with 1901 to 1905 mean
predicts_sp$climate_anomaly <- predicts_sp$avg_temp-extract(tmp1901_1905mean, predicts_sp)


#####TMAX

##Read in tmax data
tmx <- stack(tmx.path)

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


predicts_sp$tmax <- rep(0, 1)
predicts_sp$tmax_quarter <- rep(0, 1)
##As above but for tmax

system.time(for(i in 1:nrow(predicts_sp)) {
  date <- predicts_sp$Sample_end_latest[i]
  date<- substr(date,1, 7)
  date <- gsub("-", ".", date, fixed = TRUE)
  month_match <- which(names_sub_tmx==date)
  surrounding_months <- names_tmx[(month_match-11):(month_match)]
  
  max_temp <- tmx[[surrounding_months]]
  
  ## Mask to improve speed
  
  mask <- trim(rasterize(SP[i, ], max_temp[[1]]))
  crop <- crop(max_temp, mask)
  
  
  #Extract average temperature 
  predicts_sp$tmax[i] <-  max(cellStats(crop, stat = "mean"))
  predicts_sp$tmax_quarter[i] <- max_quarter_fast(cellStats(crop, stat = "mean"))
  
  
  
})

head(predicts_sp)


predicts_sp$tmax_anomaly <- predicts_sp$tmax - predicts_sp$tmax_baseline
predicts_sp$tmax_quarter_anomaly <- predicts_sp$tmax_quarter - predicts_sp$tmax_baseline

####### Historic variability in temperatures (Standard deviation of 1901 to 1905 Mean Monthly T )

tmp1901_1905sd <- calc(tmp1901_1905, stats::sd)
predicts_sp$historic_sd <- extract(tmp1901_1905sd, predicts_sp)
tmx1901_1905_sd <- calc(tmx1901_1905, stats::sd)
predicts_sp$historic_sd_tmax <- extract(tmx1901_1905_sd, predicts_sp)



##Save results


saveRDS(predicts_sp, output)




