##%######################################################%##
#                                                          #
####            Testing different thresholds            ####
####              for active temperatures               ####
#                                                          #
##%######################################################%##

# testing some lower temperature thresholds for insect active months. 
# There are going to be some insects active at temps lower than 10 degrees. 

rm(list = ls())


# directories
dataDir <- "0_data/"
predictsDir <- "1_PreparePREDICTSData/"
outDir <- "14_Additional_Tests/"
if(!dir.exists(outDir)) dir.create(outDir)

sink(paste0(outDir,"log.txt"))

t.start <- Sys.time()

print(t.start)

##Read in packages
library(raster)
library(sp)
library(dismo)
library(rgdal)
library(predictsFunctions)
library(Rfast)
library(snow)
source("Functions.R")

##Path for monthly mean temperature from CRUv4.03
tmp.path <- "Data/cru_ts4.03.1901.2018.tmp.dat.nc"

##Path for mean monthly maximum from CRUv4.03
tmx.path <- "Data/cru_ts4.03.1901.2018.tmx.dat.nc"

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
#thresh <- 8
thresh <- 6


#### calculating the mean based anomaly for each site ####

# this now looks at the 5 years preceding each sample data
# assesses which months are above the temp threshold
# then uses these months to calculate the present and baseline temperature means
# and the baseline temperature sd. 

# Time difference of 13.14845 mins 

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
    # other instances where points do not line up with the tmp layers (in the sea?)
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

print(st2 - st1) # Time difference of 12.9704 mins


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

# organise the data
predictsSites <- predicts_sp@data


# remove the Urban sites
#predictsSites$UI2[(predictsSites$UI2=="Urban")] <- NA
predictsSites$UI2 <- factor(predictsSites$UI2)
predictsSites$UI2 <- relevel(predictsSites$UI2,ref="Primary vegetation")

# organise the climate anomaly data
#predictsSites$TmeanAnomaly <- predictsSites$climate_anomaly
#predictsSites$StdTmeanAnomaly <- predictsSites$TmeanAnomaly / predictsSites$historic_sd
predictsSites$StdTmeanAnomalyRS <- StdCenterPredictor(predictsSites$StdTmeanAnomaly)

#predictsSites$TmaxAnomaly <- predictsSites$tmax_quarter_anomaly
#predictsSites$StdTmaxAnomaly <- predictsSites$TmaxAnomaly / predictsSites$historic_sd_tmax
# rescale the variable
predictsSites$StdTmaxAnomalyRS <- StdCenterPredictor(predictsSites$StdTmaxAnomaly)

# # organise NH data - combine primary and secondary vegetation at each scale
# predictsSites$NH_1000 <- predictsSites$PV_1000 + predictsSites$SV_1000
# predictsSites$NH_3000 <- predictsSites$PV_3000 + predictsSites$SV_3000
# predictsSites$NH_5000 <- predictsSites$PV_5000 + predictsSites$SV_5000
# predictsSites$NH_10000 <- predictsSites$PV_10000 + predictsSites$SV_10000
# 
# 
# # rescale the variables
# predictsSites$NH_1000.rs <- StdCenterPredictor(predictsSites$NH_1000)
# predictsSites$NH_3000.rs <- StdCenterPredictor(predictsSites$NH_3000)
# predictsSites$NH_5000.rs <- StdCenterPredictor(predictsSites$NH_5000)
# predictsSites$NH_10000.rs <- StdCenterPredictor(predictsSites$NH_10000)

# rescaling abundance and log values
predictsSites <- RescaleAbundance(predictsSites)

# charlie added this line as later bits were throwing errors
predictsSites <- droplevels(predictsSites)


# some of the climate values are NA since they do not meet the thresholds
nrow(predictsSites[is.na(predictsSites$avg_temp), ]) # 19
predictsSites <- predictsSites[!is.na(predictsSites$avg_temp), ]


# # take a look at possible correlations between variables
# plot(predictsSites$avg_temp, predictsSites$TmeanAnomaly,
#      xlab = "Average temperature", 
#      ylab = "Anomaly (difference between present and baseline)")
# cor(predictsSites$avg_temp, predictsSites$TmeanAnomaly)
# 
# 
# plot(predictsSites$avg_temp, predictsSites$StdTmeanAnomaly,
#      xlab = "Average temperature", 
#      ylab = "Standardised climate anomaly")
# cor(predictsSites$avg_temp, predictsSites$StdTmeanAnomaly)


# save
saveRDS(object = predictsSites, file = paste0(outDir,"/PREDICTSSitesWithClimateData_thresh", thresh, ".rds"))

#predictsSites <- readRDS(file = paste0(outDir,"PREDICTSSiteData.rds"))

# running the model selection process

# 1. Abundance, mean anomaly
MeanAnomalyModelAbund <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS"))

# save the model output
save(MeanAnomalyModelAbund, file = paste0(outDir, "/MeanAnomalyModelAbund_thresh_", thresh, ".rdata"))

# 2. Richness, mean anomaly
MeanAnomalyModelRich <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
                                    fitFamily = "poisson",fixedFactors = "UI2",
                                    fixedTerms = list(StdTmeanAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                    saveVars = c("Total_abundance", "SSBS"))

# save model output
save(MeanAnomalyModelRich, file = paste0(outDir, "/MeanAnomalyModelRich_thresh_", thresh, ".rdata"))


#### plot the abundance and richenss mean anom models ####


# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)


pdf(file = paste0(outDir, "/Figure2_MeanAnom_Abun_Rich_thresh_", thresh, ".pdf"), width = 8, height = 4)

par(mfrow=c(1,2))
par(las=1)
par(mgp=c(1.6,0.3,0)) # 
par(mar=c(2.6,2.6,1,1)) # margins around plot
par(tck=-0.01) # tick mark size
#par(pty="s") # set plot type to be square

nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelAbund$data$UI2)))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund$model,data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# plot #

# set up plotting window
plot(-9e99,-9e99,xlim=c(min(nd$StdTmeanAnomaly),2),
     ylim=c(-70,80),
     xlab="Standardised Temperature Anomaly",ylab="Change in total abundance (%)", cex.lab = 0.8, cex.axis = 0.8)

title("a", adj = 0, cex.main = 1)


invisible(mapply(FUN = function(preds,col){
  
  preds <- na.omit(preds)
  
  X.Vec <- c(preds$StdTmeanAnomaly, max(preds$StdTmeanAnomaly), 
             rev(preds$StdTmeanAnomaly), min(preds$StdTmeanAnomaly))
  Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
             rev(preds$PredUpper), (preds$PredLower)[1])
  
  polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
  
  points(x = preds$StdTmeanAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
  
},split(nd,nd$UI2),c("#009E73", "#D55E00", "#E69F00", "#0072B2")))

# add some gridlines
abline(h=150,lty=1,col="#0000000C")
abline(h=100,lty=1,col="#0000000C")
abline(h=50,lty=1,col="#0000000C")
abline(h=0,lty=2,col="#030303")
abline(h=-50,lty=1,col="#0000000C")
abline(v=0,lty=1,col="#0000000C")
abline(v=0.5,lty=1,col="#0000000C")
abline(v=1,lty=1,col="#0000000C")
abline(v=1.5,lty=1,col="#0000000C")
abline(v=2,lty=1,col="#0000000C")

# add legend
legend(
  x = -0.1, y = 80,bty="n",
  legend = c("Primary","Secondary",
             "Agriculture_Low",
             "Agriculture_High"),
  col = c("#009E73", "#0072B2",
          "#E69F00", "#D55E00"),
  lty=1,lwd=2, cex = 0.8)


## now the species richness plot 
nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelRich$data$UI2)))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  MeanAnomalyModelRich$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelRich$model,data = nd, nIters = 10000)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# plot #

# set up plotting window
plot(-9e99,-9e99,xlim=c(min(nd$StdTmeanAnomaly),2),
     ylim=c(-70,80),
     xlab="Standardised Temperature Anomaly",ylab="Change in species richness (%)", cex.lab = 0.8, cex.axis = 0.8)

title("b", adj = 0, cex.main = 1)


invisible(mapply(FUN = function(preds,col){
  
  preds <- na.omit(preds)
  
  X.Vec <- c(preds$StdTmeanAnomaly, max(preds$StdTmeanAnomaly), 
             rev(preds$StdTmeanAnomaly), min(preds$StdTmeanAnomaly))
  Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
             rev(preds$PredUpper), (preds$PredLower)[1])
  
  polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
  
  points(x = preds$StdTmeanAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
  
},split(nd,nd$UI2),c("#009E73", "#D55E00", "#E69F00", "#0072B2")))

# add some gridlines
abline(h=150,lty=1,col="#0000000C")
abline(h=100,lty=1,col="#0000000C")
abline(h=50,lty=1,col="#0000000C")
abline(h=0,lty=2,col="#030303")
abline(h=-50,lty=1,col="#0000000C")
abline(v=0,lty=1,col="#0000000C")
abline(v=0.5,lty=1,col="#0000000C")
abline(v=1,lty=1,col="#0000000C")
abline(v=1.5,lty=1,col="#0000000C")
abline(v=2,lty=1,col="#0000000C")

# add legend
# legend(
#   x = 0,y = 90,bty="n",
#   legend = c("Primary","Secondary",
#              "Agriculture_Low",
#              "Agriculture_High"),
#   col = c("#009E73", "#0072B2",
#           "#E69F00", "#D55E00"),
#   lty=1,lwd=2, cex = 0.8)


dev.off()



t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()
