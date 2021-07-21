##%######################################################%##
#                                                          #
####   Anomaly based on seasonal activity of insects    ####
#                                                          #
##%######################################################%##

# This script will recalculate the anomaly based on just those months
# where insects can be considered active. 

# A paper by Johansson et al (2020) in Scientific Reports looks into 
# a similar issues.  They use months where the temperature is above 10 degrees C 
# as a threshold for insect activity. https://doi.org/10.1038/s41598-020-65608-7

# questions/issues to consider:
# 1. Do we only want consecutive months where temp is > 10?
# 2. Calculate all temperature using this threshold, i.e. the sample year temp
# and the baseline temp and the baseline sd?


## issues encountered:
# 1. when using mean temp >= 10 as the threshold, there are
# instances where there are no months which meet this criteria.  The same occurs
# when a threshold of >= 8 is used. 

rm(list = ls())


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
outDir <- "12_Anomaly_Recalc_active_months/"
dir.create(outDir)

##Path for monthly mean temperature from CRUv4.03
tmp.path <- "Data/cru_ts4.03.1901.2018.tmp.dat.nc"

##Path for mean monthly maximum from CRUv4.03
#tmx.path <- "Data/cru_ts4.03.1901.2018.tmx.dat.nc"

output <- "Outputs/predicts_climate_info.rds"

##Read in average temperature data from CRU v4.03
tmp <- stack(paste0(dataDir,"cru_ts4.03.1901.2018.tmp.dat.nc"),varname = "tmp")
#tmn <- stack(paste0(dataDir,"cru_ts4.03.1901.2018.tmn.dat.nc"),varname = "tmn")

##Both cru data and predicts is in WGS84
wgs84 <- crs(tmp)

# read in the predicts sites - insect subset
predicts_sites <- readRDS(paste0(predictsDir,"PREDICTSSiteData.rds"))

# Remove sites without coordinates
predicts_sites2 <- predicts_sites[!is.na(predicts_sites$Latitude), ] # this has already been done elsewhere

# Create spatial map of PREDICTS sites
predicts_sp <- SpatialPointsDataFrame(
  coords = cbind(predicts_sites2$Longitude, predicts_sites2$Latitude), 
  data = predicts_sites2, proj4string = wgs84)

##Create a layer of mean monthly temperatures from 1901 to 1930
tmp1901_1930 <- tmp[[names(tmp)[1:360]]]


##Calculate the mean of 1901 to 19030 mean monthly temperatUres
#tmp1901_1930mean <- raster::calc(tmp1901_1930, base::mean)
# we now want to calculate the mean just across months that are hotter than the threshold

##Names of tmp layer, needed for subsettting
names_tmp <- names(tmp)

##Create a list of a all names of tmp layers, that will be used for matching later on
names_sub <- substr(names_tmp, 2, 8) 

##Spatial points for rasterizing
SP <- SpatialPoints(predicts_sp, proj4string=wgs84)

### set the threshold for insect activity
thresh <- 10



nCores <- parallel::detectCores()

st1 <- Sys.time()

cl <- snow::makeCluster(nCores-1)

# Time difference 
snow::clusterExport(
  cl = cl,
  list = c('predicts_sp','names_sub','names_tmp', 'values', 'names', 'length', 'mean', 'sd',
           'tmp','SP','rasterize','crop','trim', 'grep', 'sapply', 'strsplit',
           'cellStats', 'thresh', 'tmp1901_1930'),envir = environment())

temperatureVars <- data.frame(t(parSapply(
  cl = cl,X = (1:nrow(predicts_sp)),FUN = function(i){
    #cl = cl,X = (1:10),FUN = function(i){
      
    #temperatureVars <- NULL
    
    #for(i in 1538:nrow(predicts_sp)){
      
    #print(i)
    
      #Get end sample date for sample in predicts
      sampDate <- predicts_sp$Sample_end_latest[i]
      
      #Reformat date for string matching
      sampDate <- substr(sampDate,1, 7)
      sampDate <- gsub("-", ".", sampDate, fixed = TRUE)
      
      #Match date in predicts with month in CRU climate data
      month_match <- which(names_sub==sampDate)
      
      # surrounding_months <- names[(month_match-11):(month_match)]
      # edit: use months from 5 year pre-sample, rather than 1 year
      surrounding_months <- names_tmp[(month_match-59):(month_match)]
      
      
      #Create a layer for average temperature in the year preceding end sample date
      temp <- tmp[[surrounding_months]]
      
      ## Mask to improve speed
      
      mask <- trim(rasterize(SP[i, ], temp[[1]]))
      mapCrop <- crop(temp, mask)
      
      
      if(!length(names(mapCrop)[values(mapCrop) >= thresh]) == 0 & length(values(mapCrop)[!is.na(values(mapCrop))]) > 0 ){
        

        
        #if(length(values(mapCrop)[!is.na(values(mapCrop))]) > 0 ){
          
          # Get the average temperature for each month across 5 years
          
          vals <- NULL
          
          for(j in 1:12){
            
            if(j < 10){ mon <- paste0(0, j) }else {mon <- j}
            
            monthmean <- values(mean(mapCrop[[grep(mon, sapply(strsplit(names(mapCrop), "[.]"), "[[", 2))  ]]))
            
            vals <- rbind(vals, c(mon, monthmean))

          }
          
          
          vals <- as.data.frame(vals)
          vals$V2 <- as.numeric(as.character(vals$V2))
          
          # determine months that meet the threshold
          
          #months <- mapCrop[[names(mapCrop)[values(mapCrop) >= thresh]]]
          
          # which months are the 5 year average >= the threshold

          vals <- vals[vals$V2 >= thresh, ]
          
          # which are the good months
          months <- vals$V1
          
          n_months <- length(months)
            
          avg_temp <- mean(vals$V2)
          sd_temp <- sd(vals$V2)

          #avg_temp <- mean(values(months))
          #sd_temp <- sd(values(months))
          
          #avg_temp <- mean(cellStats(mapCrop, stat = "mean")) # I don't think the two means are necessary here
          #sd_temp <- sd(cellStats(mapCrop, stat = "mean"))
          
          
          ### now work out the baseline for those months ###
          # extract month numbers
          
          #months_n <- sapply(strsplit(names(months), "[.]"), "[[", 2)    
          
          baseline <- crop(tmp1901_1930, mask)
          
          # subset the baseline to just the required months
          baseline <-  baseline[[names(baseline)[sapply(strsplit(names(baseline), "[.]"), "[[", 2) %in% months]]]
          
          mean_baseline <- mean(values(baseline))
          sd_baseline <- sd(values(baseline))
          
          ### now calc the anomaly for that site, using the site specific baselines ###
          
          Anom <- avg_temp - mean_baseline
          StdAnom <-  Anom/sd_baseline
          
          
          return(c(avg_temp=avg_temp, sd_temp=sd_temp, Anom = Anom, StdAnom = StdAnom, n_months = n_months))
         # temperatureVars <- rbind(temperatureVars, c(avg_temp=avg_temp, sd_temp=sd_temp, Anom = Anom, StdAnom = StdAnom))
          
        }else{
          mean_baseline <- NA
          sd_baseline <- NA
          Anom <- NA
          StdAnom <-  NA
          avg_temp <- NA
          sd_temp <- NA
          n_months = 0
          
         # temperatureVars <- rbind(temperatureVars, c(avg_temp=avg_temp, sd_temp=sd_temp, Anom = Anom, StdAnom = StdAnom))
          return(c(avg_temp=avg_temp, sd_temp=sd_temp, Anom = Anom, StdAnom = StdAnom, n_months = n_months))
          
        }
        
        
      }
    
)))

snow::stopCluster(cl)

st2 <- Sys.time()

print(st2 - st1) # 6.851037 mins: time on Charlie's desktop for 6095 sites

# there seems to be 26 NA using this method. Some of these are from points being in 
# the sea and others are where there are no months in the present year that are above the threshold.
temperatureVars <- as.data.frame(temperatureVars)


## add new values in temperatureVars into predicts dataset
predicts_sp$avg_temp <- temperatureVars$avg_temp
predicts_sp$avg_temp_sd <- temperatureVars$sd_temp
predicts_sp$climate_anomaly <- temperatureVars$Anom
predicts_sp$Std_climate_anomaly <- temperatureVars$StdAnom
predicts_sp$n_months <- temperatureVars$n_months



saveRDS(object = predicts_sp, file = paste0(outDir,"PREDICTSSitesWithClimateData_seasonTest_5.rds"))


table(temperatureVars$n_months)

#  0    3    4    5    6    7    8    9   10   11   12 
# 26  192 1060  811 1034  496   52  155   52   62 2155



#### organise the data ####


# ran the NH script 5 inbetween this #


# read in the predicts data
predictsSites <- readRDS(paste0(outDir,"PREDICTSSitesWithClimateData_seasonTest_5.rds"))
predictsSites <- predictsSites@data

# remove the Urban sites
#predictsSites$UI2[(predictsSites$UI2=="Urban")] <- NA
predictsSites$UI2 <- factor(predictsSites$UI2)
predictsSites$UI2 <- relevel(predictsSites$UI2,ref="Primary vegetation")

# rescale the anomaly
predictsSites$StdTmeanAnomalyRS <- StdCenterPredictor(predictsSites$Std_climate_anomaly)

# organise NH data - combine primary and secondary vegetation at each scale
predictsSites$NH_1000 <- predictsSites$PV_1000 + predictsSites$SV_1000
predictsSites$NH_3000 <- predictsSites$PV_3000 + predictsSites$SV_3000
predictsSites$NH_5000 <- predictsSites$PV_5000 + predictsSites$SV_5000
predictsSites$NH_10000 <- predictsSites$PV_10000 + predictsSites$SV_10000


# rescale the variables
predictsSites$NH_1000.rs <- StdCenterPredictor(predictsSites$NH_1000)
predictsSites$NH_3000.rs <- StdCenterPredictor(predictsSites$NH_3000)
predictsSites$NH_5000.rs <- StdCenterPredictor(predictsSites$NH_5000)
predictsSites$NH_10000.rs <- StdCenterPredictor(predictsSites$NH_10000)

# rescaling abundance and log values
predictsSites <- RescaleAbundance(predictsSites)

# charlie added this line as later bits were throwing errors
predictsSites <- droplevels(predictsSites)

# save the dataset
saveRDS(object = predictsSites,file = paste0(outDir,"PREDICTSSiteData_seasonTest_5.rds"))



##### model selection ####

# just using mean anomaly now. 

# 1. Abundance, mean anomaly
MeanAnomalyModelAbund <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyModelAbund, file = paste0(outDir, "MeanAnomalyModelAbund_5.rdata"))

# 2. Richness, mean anomaly
MeanAnomalyModelRich <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
                                    fitFamily = "poisson",fixedFactors = "UI2",
                                    fixedTerms = list(StdTmeanAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                    saveVars = c("Total_abundance", "SSBS", "NH_3000"))

# save model output
save(MeanAnomalyModelRich, file = paste0(outDir, "/MeanAnomalyModelRich_5.rdata"))



#### plots ####

exclQuantiles <- c(0.025,0.975)


pdf(file = paste0(outDir, "/MeanAnom_Abun_5.pdf"), width = 4, height = 4)

par(mfrow=c(1,1))
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
plot(-9e99,-9e99,xlim=c(min(nd$StdTmeanAnomaly),max(nd$StdTmeanAnomaly)),
     ylim=c(min(nd$PredLower,na.rm = TRUE),max(nd$PredUpper,na.rm = TRUE)),
     xlab="Standardised Temperature Anomaly",ylab="Abundance (%)", cex.lab = 0.8, cex.axis = 0.8)

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
abline(h=150,lty=1,col="#00000022")
abline(h=100,lty=1,col="#00000022")
abline(h=50,lty=1,col="#00000022")
abline(h=0,lty=1,col="#00000022")
abline(h=-50,lty=1,col="#00000022")
abline(v=0,lty=1,col="#00000022")
abline(v=1,lty=1,col="#00000022")
abline(v=2,lty=1,col="#00000022")

# add legend
legend(
  x = -1,y = 45,bty="n",
  legend = c("Primary","Secondary",
             "Agriculture_Low",
             "Agriculture_High"),
  col = c("#009E73", "#0072B2",
          "#E69F00", "#D55E00"),
  lty=1,lwd=2, cex = 0.8)


dev.off()


######## Extended data fig 2 - Mean Anom - SR ############

pdf(file = paste0(outDir, "/MeanAnom_Rich_5.pdf"), width = 4, height = 4)

par(mfrow=c(1,1))
par(las=1)
par(mgp=c(1.6,0.3,0)) # 
par(mar=c(2.6,2.6,1,1)) # margins around plot
par(tck=-0.01) # tick mark size
#par(pty="s") # set plot type to be square

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
plot(-9e99,-9e99,xlim=c(min(nd$StdTmeanAnomaly),max(nd$StdTmeanAnomaly)),
     ylim=c(min(nd$PredLower,na.rm = TRUE),max(nd$PredUpper,na.rm = TRUE)),
     xlab="Standardised Temperature Anomaly",ylab="Species Richness (%)", cex.lab = 0.8, cex.axis = 0.8)

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
abline(h=150,lty=1,col="#00000022")
abline(h=100,lty=1,col="#00000022")
abline(h=50,lty=1,col="#00000022")
abline(h=0,lty=1,col="#00000022")
abline(h=-50,lty=1,col="#00000022")
abline(v=0,lty=1,col="#00000022")
abline(v=1,lty=1,col="#00000022")
abline(v=2,lty=1,col="#00000022")

# add legend
legend(
  x = -1,y = 50,bty="n",
  legend = c("Primary","Secondary",
             "Agriculture_Low",
             "Agriculture_High"),
  col = c("#009E73", "#0072B2",
          "#E69F00", "#D55E00"),
  lty=1,lwd=2, cex = 0.8)


dev.off()





