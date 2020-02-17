##This script will obtain climate data for all sites in the predicts database
##For the script to work, you will need to download and unzip the CRUv4.03 climate data
## and provide put the file path in this script

##Read in packages
library(raster)
library(sp)
library(dismo)
library(rgdal)
library(predictsFunctions)

##Also makes sure the convenience_functions.r file has been run

##The following file paths are required for this script to run:
##Path for predicts database
predicts.path <-"Your path here"

##Path for monthly mean temperature from CRUv4.03
tmp.path <- "Your path here"

##Path for mean monthly maximum from CRUv4.03
tmx.path <- "Your path here"

##
percNH.path <- "Your path here"


##Location for where you'd like the climate data to be saved (as a csv file)
output <- "Your path here"


#Read in predicts
predicts <- readRDS(predicts.path)

##Merge sites
predicts_sites <- create_sites(predicts)


###Climate


##Predicts coordinantes
coords_predicts_sites <- cbind(predicts_sites$Longitude, predicts_sites$Latitude)


##Read in average temperature data from CRU v4.03
tmp <- stack(tmp.path)

##Create a layer of mean monthly temperatures from 1901 to 1905
tmp1901_1905 <- tmp[[names(tmp)[1:60]]]
names(tmp1901_1905)

##Calculate the mean of 1901 to 1905 mean monthly temperatyres
tmp1901_1905mean <- raster::calc(tmp1901_1905, base::mean)

##Calculate the standard deviation of 1901 to 1905 mean monthly temperatures 
tmp1901_1905sd <- calc(tmp1901_1905, stats::sd)  


nrow(predicts_sites) 
##Names of tmp layer, needed for subsettting
names <- names(tmp)
nrow(predicts_sites)

##Create a list of a all names of tmp layers, that will be used for matching later on
names_sub <- substr(names, 2, 8) 
names_sub
##Empty vector for storing matching months
month_match <-vector()

predicts_sites$rate <- rep(0, 1)
predicts_sites$avg_temp_sd <- rep(0, 1)
predicts_sites$avg_temp <- rep(0, 1)
head(predicts_sites)



###For each site in predicts find climate anomaly, match with CRU climate data
for(i in 1:nrow(predicts_sites)) {
  #Get end sample date for sample in predicts
  date <- predicts_sites$Sample_end_latest[i]
  #Reformat date for string matching
  date<- substr(date,1, 7)
  date <- gsub("-", ".", date, fixed = TRUE)
  #Match date in predicts with month in CRU climate data
  month_match <- which(names_sub==date)
  surrounding_months <- names[(month_match-11):(month_match)]
  #Create a layer for average temperature in the year preceding end sample date
  avg_temp <- calc(tmp[[surrounding_months]], base::mean)
  #Get coordinates for predicts sites
  coords <- matrix(coords_predicts_sites[i, ], ncol=2)
  #Extract average temperature 
  predicts_sites$avg_temp[i] <- extract(avg_temp, coords)
  #Calculate climate anomaly from sample compared with 1901 to 1905 mean
  rate <- avg_temp-tmp1901_1905mean
  predicts_sites$rate[i] <- extract(rate, coords)
  
}

month_match <-vector()

####As above but for sd of average temperature in year preceding sample end date
for(i in 1:nrow(predicts_sites)) {
  date <- predicts_sites$Sample_end_latest[i]
  date<- substr(date,1, 7)
  date <- gsub("-", ".", date, fixed = TRUE)
  month_match <- which(names_sub==date)
  surrounding_months <- names[(month_match-11):(month_match)]
  avg_temp_sd <- calc(tmp[[surrounding_months]], stats::sd)
  coords <- matrix(coords_predicts_sites[i, ], ncol=2)
  predicts_sites$avg_temp_sd[i] <- extract(avg_temp_sd, coords)
  
}


#####TMAX

##Read in tmax data
tmx <- stack(tmx.path)

#Create raster stack for 1901 to 1905
tmx1901_1905 <- tmx[[names(tmx)[1:60]]]
names(tmx1901_1905)



###Calculate maximum mean monthly temperature for each year from 1901 to 1905
tmx1901max <- calc(tmx1901_1905[[(names(tmx)[1:12])]], base::max)
tmx1902max <- calc(tmx1901_1905[[(names(tmx)[13:24])]], base::max)
tmx1903max <- calc(tmx1901_1905[[(names(tmx)[25:36])]], base::max)
tmx1904max <- calc(tmx1901_1905[[(names(tmx)[37:48])]], base::max)
tmx1905max <- calc(tmx1901_1905[[(names(tmx)[49:60])]], base::max)

plot(tmx1901max)

##Stack maximum temp
tmx1901_1905max <- stack(tmx1901max, tmx1902max, tmx1903max, tmx1904max, tmx1905max)

##Calculate mean of maximum mean monthy temperatures between 1901 and 1905
tmx1901_1905mean.max <- calc(tmx1901_1905max, base::mean)
plot(tmx1901_1905mean.max) 


##names of tmax layer needed for subsetting
names_tmx <- names(tmx)

##For matching
names_sub_tmx <- substr(names_tmx, 2, 8)


##To store values
predicts_sites$tmax <- rep(0,1)
predicts_sites$tmax_anomaly <- rep(0,1)
nrow(predicts_sites)

##empty month match
month_match <-vector()

##Same as above but maxing
for(i in 1:nrow(predicts_sites)) {
  date <- predicts_sites$Sample_end_latest[i]
  date<- substr(date,1, 7)
  date <- gsub("-", ".", date, fixed = TRUE)
  month_match <- which(names_sub_tmx==date)
  surrounding_months <- names_tmx[(month_match-11):(month_match)]
  max_temp <- calc(tmx[[surrounding_months]], base::max)
  coords <- matrix(coords_predicts_sites[i, ], ncol=2)
  predicts_sites$tmax[i] <- extract(max_temp, coords)
  tmax_anomaly <- max_temp-tmx1901_1905mean.max
  predicts_sites$tmax_anomaly[i] <- extract(tmax_anomaly, coords)
  
}



####### Historic variability in temperatures (Standard deviation of 1901 to 1905 Mean Monthly T )

tmp1901_1905sd <- calc(tmp1901_1905, stats::sd)
predicts_sites$historic_sd <- extract(tmp1901_1905sd, coords_predicts_sites)


##Extract percent natural; habitat data

##Create a versions of predicts sites with no NAs
predicts_sites_no_nas <- predicts_sites[!is.na(predicts_sites$Longitude) & !is.na(predicts_sites$Latitude), ]
## WGS84 CRS
crs <- CRS("+init=epsg:4326")
crs
##Create spatial points object for predicts sites
coords_predicts_sites <- SpatialPoints(cbind(predicts_sites_no_nas$Longitude, predicts_sites_no_nas$Latitude), proj4string =crs)

##Read in percent natural data
percent_natural <- raster(percNH.path)

##Extract percent natural data
predicts_sites_no_nas$percNH <- extract(percent_natural, coords_predicts_sites)

##Now we need recombine with the orginal frame

#We need to make on site vector a character string in order to match factors
predicts_sites_no_nas$SSBS <- as.character(predicts_sites_no_nas$SSBS)


##Create a vector to store percNH values. After storage, values that are 2000 will be NAs
predicts_sites$percNH <- rep(2000,1)


##Slot data into orginal predicts data frame, matching by site
for(i in 1:nrow(predicts_sites_no_nas)) {
  s <- unique(predicts_sites_no_nas$SSBS)[i]
  predicts_sites[predicts_sites$SSBS==s, ]$percNH <- predicts_sites_no_nas$percNH[predicts_sites_no_nas$SSBS==s]
  print(i)
}

##Quick check to make the data has been stored in the right place
length(unique(predicts_sites$SSBS))
length(unique(predicts_sites_no_nas$SSBS))
s <- unique(predicts_sites$SSBS)[22008]

#This should be true
predicts_sites_no_nas[predicts_sites_no_nas$SSBS==s, ]$percNH == predicts_sites[predicts_sites$SSBS==s, ]$percNH



##Change values of 2000 to NAs

a <- which(predicts_sites$percNH==2000)
predicts_sites[a, ]$percNH <- NA

##Save results
write.csv(predicts_sites, output)



