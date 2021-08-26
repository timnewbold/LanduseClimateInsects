##%######################################################%##
#                                                          #
####            Testing different baselines             ####
#                                                          #
##%######################################################%##

# In this script, different baselines for the climate anomaly are tested.
# In the main analysis the baseline years are 1901-1905. 

# alternatives to test:
# 5 year, 10 year and 20 year baselines

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
source("Functions.R")
library(reshape)
library(predictsFunctions)
library(Rfast)
library(snow)

# directories
dataDir <- "0_data/"
outDir <- "10_SCA_Baseline_testing/"
dir.create(outDir)
predictsDataDir <- "6_RunLUClimateModels/"


# data required for the plots is created in script 3_PrepareClimateIndecMaps_update where 
# the map data can be generated for each of the different baseline lengths.

# load in the datasets 

load(paste0("3_PrepareClimateIndexMaps", "/Map_data_tempvars_2004_06.rdata"))
plot_data_0130 <- temperatureVars
load(paste0(outDir, "Map_data_tempvars_2004_06_thresh_10_baseline5.rdata"))
plot_data_0105 <- temperatureVars
load(paste0(outDir, "Map_data_tempvars_2004_06_thresh_10_baseline10.rdata"))
plot_data_0110 <- temperatureVars
load(paste0(outDir, "Map_data_tempvars_2004_06_thresh_10_baseline20.rdata"))
plot_data_0120 <- temperatureVars

# load in the mean temperature data from CRU
tmp <- stack(paste0(dataDir,"cru_ts4.03.1901.2018.tmp.dat.nc"),varname="tmp")


# take the current predicts time data
tmp2004_6 <- tmp[[names(tmp)[1237:1272]]]


# which raster to use as present
pre_ras <- tmp2004_6 
#pre_ras <- tmp2016_18

# get vector of positions in values of raster that are not NA
ras <- pre_ras[[1]]

# create a set of points based on the non NA cells

wgs84 <- crs(tmp)

pnts <- rasterToPoints(ras, spatial = T)

SP <- SpatialPoints(pnts, proj4string=wgs84)  


SP_df <- as.data.frame(SP)

# add data to lat/longs of global cells
all_data <- cbind(SP_df, plot_data_0105[, 4], plot_data_0110[, 4], plot_data_0120[, 4], plot_data_0130[, 4])

names(all_data)[3:6] <- c("1901-1905", "1901-1910", "1901-1920", "1901-1930")

all_data <- melt(all_data, id.vars = c("x", "y"))


# create maps of each anomaly

all_plot <- all_data[!is.na(all_data$value), ]

# organise breaks, colours and labels
brks <- c(-2,-0.5,-0.2,-0.1,0,0.1,0.5,0.75,1,1.5,2,5,50)
cols <- c(rev(brewer.pal(n = 8,name = "Greens"))[4:8],
          (brewer.pal(n = 8,name = "Purples"))[4:6],
          (brewer.pal(n = 8,name = "Oranges"))[5:8])
labs <- c("-2 : -0.5","-0.5 : -0.2","-0.2 : -0.1","-0.1 : 0",
          "0 : 0.1","0.1 : 0.5","0.5 : 0.75","0.75 : 1","1 : 1.5","1.5 : 2","2 : 5","> 5")

# assign values into bins
all_plot$bins <- cut(all_plot$value, 
                     breaks = brks, 
                     labels = labs,
                     include.lowest = TRUE)


names(all_plot)[3] <- "year"

year_labs <- c("A. 1901-1905", "B. 1901-1910", "C. 1901-1920", "D. 1901-1930")
names(year_labs) <- unique(all_plot$year)

all_plot$year <- as.factor(all_plot$year)

View(all_plot[is.na(all_plot$bins), ])

all_plot <- all_plot[!is.na(all_plot$bins), ]


# plot
ggplot(all_plot) + 
  geom_raster(aes(x = x, y = y, fill = bins), alpha = 0.9) +
  #scale_fill_viridis_c(option = "magma", values = c(0, 0.2, 1)) + 
  scale_fill_manual(values = cols) +
  facet_wrap(~ year, nrow = 2, labeller = labeller(year = year_labs)) +
  xlab("") +
  ylab("") +
  labs(fill = "Standardised\nTemperature\nAnomaly") +
  theme_bw() +
  theme(legend.position = 'bottom', 
        #panel.border = element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_blank(),
        #legend.key.width = unit(3, "cm"),
        axis.ticks = element_blank(), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), 
        #legend.key.size = unit(0.2,"cm",
        panel.border = element_rect(size = 0.2),
        strip.background = element_rect(size = 0.2, fill = "transparent"),
        strip.text = element_text(size = 10))


# save plot as pdf
ggsave(filename = paste0(outDir, "/Supplementary_baseline_test_anomaly_maps.pdf"), width = 12, height = 8)



##%######################################################%##
#                                                          #
####               Get the values for the               ####
####        dif baselines for each predicts site        ####
#                                                          #
##%######################################################%##


### set the threshold temperature for insect activity
thresh <- 10

# read in the predicts sites - insect subset
predicts_sites <- readRDS(paste0(predictsDataDir,"PREDICTSSiteData.rds"))

# Remove sites without coordinates
predicts_sites2 <- predicts_sites[!is.na(predicts_sites$Latitude), ]

# Create spatial map of PREDICTS sites
predicts_sp <- SpatialPointsDataFrame(
  coords = cbind(predicts_sites2$Longitude, predicts_sites2$Latitude), 
  data = predicts_sites2, proj4string = wgs84)

##Spatial points for rasterizing
SP <- SpatialPoints(predicts_sp, proj4string=wgs84)


##Names of tmp layer, needed for subsettting
names_tmp <- names(tmp)

##Create a list of a all names of tmp layers, that will be used for matching later on
names_sub <- substr(names_tmp, 2, 8) 


#### calculating the mean based anomaly for each site ####

# this now looks at the 5 years preceding each sample data
# assesses which months are above the temp threshold
# then uses these months to calculate the present and baseline temperature means
# and the baseline temperature sd. 

# Time difference of 7.857316 mins

# run the code for each baseline

nCores <- parallel::detectCores()

st1 <- Sys.time()

cl <- snow::makeCluster(nCores-1)

# export to clusters
snow::clusterExport(
  cl = cl,
  list = c('predicts_sp','names_sub','names_tmp', 'values', 'names', 'length', 'mean', 'sd',
           'tmp', 'SP','rasterize','crop','trim', 'grep', 'sapply', 'strsplit',
           'cellStats', 'thresh', 'tmp1901_1930', 'tmp1901_1905', 'tmp1901_1910', 'tmp1901_1920'),envir = environment())

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
    #temp_mx <- tmx[[surrounding_months]]
    
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
      #sd_temp <- sd(vals$V2) # don't actually use this
      
      
      
      # get 3 hottest month vals for those months that meet the threshold
      
      # loop through each year, get the 3 hottest months
      
      # mask_mx <- trim(rasterize(SP[i, ], temp_mx[[1]]))
      # mapCrop_mx <- crop(temp_mx, mask_mx)
      # 
      # # for each year get the 3 hottest monts then take the mean of these
      # yr1 <- mean(sort(values(mapCrop_mx[[1:12]]), decreasing = T)[1:3])
      # yr2 <- mean(sort(values(mapCrop_mx[[13:24]]), decreasing = T)[1:3])
      # yr3 <- mean(sort(values(mapCrop_mx[[25:36]]), decreasing = T)[1:3])
      # yr4 <- mean(sort(values(mapCrop_mx[[37:48]]), decreasing = T)[1:3])
      # yr5 <- mean(sort(values(mapCrop_mx[[49:60]]), decreasing = T)[1:3])
      # 
      # # then take the mean across years to get one max value for the present
      # # this is the mean for the "present day" max (mean of 3 hottest months)
      # max_temp <- mean(yr1, yr2, yr3, yr4, yr5)
      
      
      ### now work out the baseline mean and sd for the active months and hottest months ###
      
      # get the values for that grid cell across all years
      #baseline <- crop(tmp1901_1905, mask)
      #baseline <- crop(tmp1901_1910, mask)
      #baseline <- crop(tmp1901_1920, mask)
      baseline <- crop(tmp1901_1930, mask)
      
      # subset the baseline to just the required months
      baseline <-  baseline[[names(baseline)[sapply(strsplit(names(baseline), "[.]"), "[[", 2) %in% months]]]
      
      # get the mean and sd
      mean_baseline <- mean(values(baseline))
      sd_mean_baseline <- sd(values(baseline))
      
      
      
      ### now max temps, get 3 hottest years for each year ###
      #baseline_max <- crop(tmx1901_1930, mask_mx)
      
      # for each year, get the 3 hottest years #
      # get the year names
      # yrs <- unique(gsub("\\..*", "", names(baseline_max)))
      # 
      # max_vals <- NULL
      # 
      # for(yr in yrs){
      #   
      #   mx3 <- sort(values(baseline_max[[grep(yr, names(baseline_max))]]), decreasing = T)[1:3]
      #   
      #   max_vals <- c(max_vals, mx3)
      #   
      # }
      # 
      # # get the baseline values for this site
      # mean_baseline_mx <- mean(max_vals)
      # sd_baseline_mx <- sd(max_vals)
      # 
      
      ### now calc the anomaly for that site, using the site specific baselines ###
      
      Anom <- avg_temp - mean_baseline
      StdAnom <-  Anom/sd_mean_baseline
      
      # tmax_anomaly <- max_temp - mean_baseline_mx
      # StdTmaxAnomaly <- tmax_anomaly/sd_baseline_mx
      # 
      
      return(c(avg_temp=avg_temp, Anom = Anom, StdAnom = StdAnom, n_months = n_months#,
               #max_temp = max_temp, tmax_anomaly = tmax_anomaly, StdTmaxAnomaly = StdTmaxAnomaly
               ))
      
      
      #temperatureVars <- rbind(temperatureVars, c(avg_temp=avg_temp, Anom = Anom, StdAnom = StdAnom, n_months = n_months,
      #max_temp = max_temp, tmax_anomaly = tmax_anomaly, StdTmaxAnomaly = StdTmaxAnomaly))
      
    }else{
      avg_temp <- NA
      Anom <- NA
      StdAnom <- NA
      n_months = 0
      #max_temp <- NA
      #tmax_anomaly <- NA
      #StdTmaxAnomaly <- NA
      
      #temperatureVars <- rbind(temperatureVars, c(avg_temp=avg_temp, Anom = Anom, StdAnom = StdAnom, n_months = n_months,
      #max_temp = max_temp, tmax_anomaly = tmax_anomaly, StdTmaxAnomaly = StdTmaxAnomaly))
      
      return(c(avg_temp=avg_temp, Anom = Anom, StdAnom = StdAnom, n_months = n_months#,
               #max_temp = max_temp, tmax_anomaly = tmax_anomaly, StdTmaxAnomaly = StdTmaxAnomaly
               ))      
    }
    
    
  }
  
)))


snow::stopCluster(cl)

st2 <- Sys.time()

print(st2 - st1) # Time difference of 7.857316 mins for longest baseline

# save each iteration separately
temperatureVars_1901_1905 <- temperatureVars
temperatureVars_1901_1910 <- temperatureVars
temperatureVars_1901_1920 <- temperatureVars
temperatureVars_1901_1930 <- temperatureVars

# # organise the anomaly info along with the predicts data
# temperatureVars <- as.data.frame(temperatureVars)
# 

## add new values in temperatureVars into predicts dataset
predicts_sp$StdTmeanAnomaly_1901_1905 <- temperatureVars_1901_1905$StdAnom
predicts_sp$StdTmeanAnomaly_1901_1910 <- temperatureVars_1901_1910$StdAnom
predicts_sp$StdTmeanAnomaly_1901_1920 <- temperatureVars_1901_1920$StdAnom
predicts_sp$StdTmeanAnomaly_1901_1930 <- temperatureVars_1901_1930$StdAnom



# save
saveRDS(object = predicts_sp, file = paste0(outDir,"PREDICTSSitesWithClimateData_baselines.rds"))


#predicts_sp <- readRDS(file = paste0(outDir,"PREDICTSSitesWithClimateData_baselines.rds"))
#names(predicts_sp)


##%######################################################%##
#                                                          #
####    Run models using new baselines for SCA         ####
#                                                          #
##%######################################################%##

# run the climate:LU model for abundance with each new baseline to compare results


# centre the predictors
predicts_sp$anomStd_0105RS <- StdCenterPredictor(predicts_sp$StdTmeanAnomaly_1901_1905)
predicts_sp$anomStd_0110RS <- StdCenterPredictor(predicts_sp$StdTmeanAnomaly_1901_1910)
predicts_sp$anomStd_0120RS <- StdCenterPredictor(predicts_sp$StdTmeanAnomaly_1901_1920)
predicts_sp$anomStd_0130RS <- StdCenterPredictor(predicts_sp$StdTmeanAnomaly_1901_1930)

predictsSites <- as.data.frame(predicts_sp)

# due to using active months there are some NAs for the anomalies, remove NA
predictsSites <- predictsSites[!is.na(predictsSites$anomStd_0130RS), ] # 6069 rows
 
# run abundance model using each of the anomalies

# 1. Abundance, mean anomaly 1901-1905 baseline
MeanAnomalyAbund_0105 <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(anomStd_0105RS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(anomStd_0105RS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyAbund_0105, file = paste0(outDir, "/MeanAnomalyAbund_0105.rdata"))

summary(MeanAnomalyAbund_0105$model)
# LogAbund ~ UI2 + poly(anomStd_0105RS, 1) + UI2:poly(anomStd_0105RS,      1) + (1 | SS) + (1 | SSB)


# 2. Abundance, mean anomaly 1901-1910 baseline
MeanAnomalyAbund_0110 <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(anomStd_0110RS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(anomStd_0110RS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyAbund_0110, file = paste0(outDir, "/MeanAnomalyAbund_0110.rdata"))


summary(MeanAnomalyAbund_0110$model)
# LogAbund ~ UI2 + UI2:poly(anomStd_0110RS, 1) + poly(anomStd_0110RS,      1) + (1 | SS) + (1 | SSB)

# 3. Abundance, mean anomaly  1901-1920 baseline
MeanAnomalyAbund_0120 <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(anomStd_0120RS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(anomStd_0120RS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyAbund_0120, file = paste0(outDir, "/MeanAnomalyAbund_0120.rdata"))


summary(MeanAnomalyAbund_0120$model)
#  LogAbund ~ UI2 + UI2:poly(anomStd_0120RS, 1) + poly(anomStd_0120RS,      1) + (1 | SS) + (1 | SSB)


# 4. Abundance, mean anomaly  1901-1930 baseline
MeanAnomalyAbund_0130 <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(anomStd_0130RS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(anomStd_0130RS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyAbund_0130, file = paste0(outDir, "/MeanAnomalyAbund_0130.rdata"))

summary(MeanAnomalyAbund_0130$model)
#  LogAbund ~ UI2 + UI2:poly(anomStd_0130RS, 1) + poly(anomStd_0130RS,      1) + (1 | SS) + (1 | SSB)


######## plot figures ######## 
  
#load(paste0(outDir, "/MeanAnomalyAbund_0105.rdata"))
#load(paste0(outDir, "/MeanAnomalyAbund_0110.rdata"))
#load(paste0(outDir, "/MeanAnomalyAbund_0120.rdata"))
#load(paste0(outDir, "/MeanAnomalyAbund_0130.rdata"))




# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)


nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyAbund_0105$data$anomStd_0105RS),
                        to = max(MeanAnomalyAbund_0105$data$anomStd_0105RS),
                        length.out = 200),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyAbund_0105$data$UI2)))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly_1901_1905)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyAbund_0105$data$anomStd_0105RS[
  MeanAnomalyAbund_0105$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyAbund_0105$data$anomStd_0105RS[
  MeanAnomalyAbund_0105$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyAbund_0105$data$anomStd_0105RS[
  MeanAnomalyAbund_0105$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyAbund_0105$data$anomStd_0105RS[
  MeanAnomalyAbund_0105$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# rename variable to match model output
names(nd)[1] <-"anomStd_0105RS"

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyAbund_0105$model,data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$anomStd_0105RS < QPV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$anomStd_0105RS > QPV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$anomStd_0105RS < QSV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$anomStd_0105RS > QSV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$anomStd_0105RS < QAL[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$anomStd_0105RS > QAL[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$anomStd_0105RS < QAH[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$anomStd_0105RS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


p1 <- ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  theme_bw() + 
  #labs(fill = "% NH", col = "% NH") + 
  ylab("Total Abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 120)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = c(0.3, 0.85)) +
  ggtitle("A. baseline: 1901-1905")



# 2.

nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyAbund_0110$data$anomStd_0110RS),
                        to = max(MeanAnomalyAbund_0110$data$anomStd_0110RS),
                        length.out = 200),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyAbund_0110$data$UI2)))

# back transform the predictors
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly_1901_1910)

# set richness and abundance to 0 - to be predicted
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyAbund_0110$data$anomStd_0110RS[
  MeanAnomalyAbund_0110$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyAbund_0110$data$anomStd_0110RS[
  MeanAnomalyAbund_0110$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyAbund_0110$data$anomStd_0110RS[
  MeanAnomalyAbund_0110$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyAbund_0110$data$anomStd_0110RS[
  MeanAnomalyAbund_0110$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# rename variable to match model output
names(nd2)[1] <-"anomStd_0110RS"

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyAbund_0110$model,data = nd2)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$anomStd_0110RS < QPV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$anomStd_0110RS > QPV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$anomStd_0110RS < QSV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$anomStd_0110RS > QSV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$anomStd_0110RS < QAL[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$anomStd_0110RS > QAL[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$anomStd_0110RS < QAH[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$anomStd_0110RS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd2$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


p2 <- ggplot(data = nd2, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd2$PredLower, ymax = nd2$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  theme_bw() + 
  #labs(fill = "% NH", col = "% NH") + 
  ylab("Total Abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 120)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = c(0.3, 0.85)) +
  ggtitle("B. baseline: 1901-1910")




# 3.

nd3 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyAbund_0120$data$anomStd_0120RS),
                        to = max(MeanAnomalyAbund_0120$data$anomStd_0120RS),
                        length.out = 300),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyAbund_0120$data$UI2)))

# back transform the predictors
nd3$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly_1901_1920)

# set richness and abundance to 0 - to be predicted
nd3$LogAbund <- 0
nd3$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd3$UI2=="Primary vegetation") & (nd3$StdTmeanAnomaly==min(abs(nd3$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyAbund_0120$data$anomStd_0120RS[
  MeanAnomalyAbund_0120$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyAbund_0120$data$anomStd_0120RS[
  MeanAnomalyAbund_0120$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyAbund_0120$data$anomStd_0120RS[
  MeanAnomalyAbund_0120$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyAbund_0120$data$anomStd_0120RS[
  MeanAnomalyAbund_0120$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# rename variable to match model output
names(nd3)[1] <-"anomStd_0120RS"

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyAbund_0120$model,data = nd3)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$anomStd_0120RS < QPV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$anomStd_0120RS > QPV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$anomStd_0120RS < QSV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$anomStd_0120RS > QSV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$anomStd_0120RS < QAL[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$anomStd_0120RS > QAL[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$anomStd_0120RS < QAH[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$anomStd_0120RS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd3$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd3$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3$UI2 <- factor(nd3$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


p3 <- ggplot(data = nd3, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd3$PredLower, ymax = nd3$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  theme_bw() + 
  #labs(fill = "% NH", col = "% NH") + 
  ylab("Total Abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 120)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = c(0.3, 0.85)) +
  ggtitle("C. baseline: 1901-1920")






# 4.

nd4 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyAbund_0130$data$anomStd_0130RS),
                        to = max(MeanAnomalyAbund_0130$data$anomStd_0130RS),
                        length.out = 300),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyAbund_0130$data$UI2)))

# back transform the predictors
nd4$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd4$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly_1901_1930)

# set richness and abundance to 0 - to be predicted
nd4$LogAbund <- 0
nd4$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd4$UI2=="Primary vegetation") & (nd4$StdTmeanAnomaly==min(abs(nd4$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyAbund_0130$data$anomStd_0130RS[
  MeanAnomalyAbund_0130$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyAbund_0130$data$anomStd_0130RS[
  MeanAnomalyAbund_0130$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyAbund_0130$data$anomStd_0130RS[
  MeanAnomalyAbund_0130$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyAbund_0130$data$anomStd_0130RS[
  MeanAnomalyAbund_0130$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# rename variable to match model output
names(nd4)[1] <-"anomStd_0130RS"

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyAbund_0130$model,data = nd4)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$anomStd_0130RS < QPV[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$anomStd_0130RS > QPV[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$anomStd_0130RS < QSV[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$anomStd_0130RS > QSV[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$anomStd_0130RS < QAL[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$anomStd_0130RS > QAL[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$anomStd_0130RS < QAH[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$anomStd_0130RS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd4$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd4$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4$UI2 <- factor(nd4$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


p4 <- ggplot(data = nd4, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd4$PredLower, ymax = nd4$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  theme_bw() + 
  #labs(fill = "% NH", col = "% NH") + 
  ylab("Total Abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 120)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = c(0.3, 0.85)) +
  ggtitle("D. baseline: 1901-1930")




# organise plots
cowplot::plot_grid(p1, p2, p3, p4)

ggsave(filename = paste0(outDir, "Extended_Data_Baselines_plots_Abun.pdf"), width = 10, height = 10, units = "in")





#### species richness models ####

# 1. SR, mean anomaly 1901-1905 baseline
MeanAnomalyRich_0105 <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
                                     fitFamily = "poisson",fixedFactors = "UI2",
                                     fixedTerms = list(anomStd_0105RS=1),
                                     randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                     fixedInteractions = c("UI2:poly(anomStd_0105RS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyRich_0105, file = paste0(outDir, "/MeanAnomalyRich_0105.rdata"))




# 2. SR, mean anomaly 1901-1910 baseline
MeanAnomalyRich_0110 <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
                                     fitFamily = "poisson",fixedFactors = "UI2",
                                     fixedTerms = list(anomStd_0110RS=1),
                                     randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                     fixedInteractions = c("UI2:poly(anomStd_0110RS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyRich_0110, file = paste0(outDir, "/MeanAnomalyRich_0110.rdata"))


# 3. SR, mean anomaly  1901-1920 baseline
MeanAnomalyRich_0120 <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
                                     fitFamily = "poisson",fixedFactors = "UI2",
                                     fixedTerms = list(anomStd_0120RS=1),
                                     randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                     fixedInteractions = c("UI2:poly(anomStd_0120RS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyRich_0120, file = paste0(outDir, "/MeanAnomalyRich_0120.rdata"))


# 4. SR, mean anomaly  1901-1920 baseline
MeanAnomalyRich_0130 <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
                                    fitFamily = "poisson",fixedFactors = "UI2",
                                    fixedTerms = list(anomStd_0130RS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(anomStd_0130RS,1)"),
                                    saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyRich_0130, file = paste0(outDir, "/MeanAnomalyRich_0130.rdata"))





# take a look at model outputs
summary(MeanAnomalyRich_0105$model)
summary(MeanAnomalyRich_0110$model)
summary(MeanAnomalyRich_0120$model)
summary(MeanAnomalyRich_0130$model)



######## plot figures ######## 

#load(paste0(outDir, "/MeanAnomalyRich_0105.rdata"))
#load(paste0(outDir, "/MeanAnomalyRich_0110.rdata"))
#load(paste0(outDir, "/MeanAnomalyRich_0120.rdata"))
#load(paste0(outDir, "/MeanAnomalyRich_0130.rdata"))




# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)


nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyRich_0105$data$anomStd_0105RS),
                        to = max(MeanAnomalyRich_0105$data$anomStd_0105RS),
                        length.out = 200),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyRich_0105$data$UI2)))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly_1901_1905)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyRich_0105$data$anomStd_0105RS[
  MeanAnomalyRich_0105$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyRich_0105$data$anomStd_0105RS[
  MeanAnomalyRich_0105$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyRich_0105$data$anomStd_0105RS[
  MeanAnomalyRich_0105$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyRich_0105$data$anomStd_0105RS[
  MeanAnomalyRich_0105$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# rename variable to match model output
names(nd)[1] <-"anomStd_0105RS"

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyRich_0105$model,data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$anomStd_0105RS < QPV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$anomStd_0105RS > QPV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$anomStd_0105RS < QSV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$anomStd_0105RS > QSV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$anomStd_0105RS < QAL[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$anomStd_0105RS > QAL[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$anomStd_0105RS < QAH[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$anomStd_0105RS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


p1 <- ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  theme_bw() + 
  #labs(fill = "% NH", col = "% NH") + 
  ylab("Species Richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-75, 150)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = c(0.3, 0.85)) +
  ggtitle("A. baseline: 1901-1905")



# 2.

nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyRich_0110$data$anomStd_0110RS),
                        to = max(MeanAnomalyRich_0110$data$anomStd_0110RS),
                        length.out = 200),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyRich_0110$data$UI2)))

# back transform the predictors
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly_1901_1910)

# set richness and abundance to 0 - to be predicted
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyRich_0110$data$anomStd_0110RS[
  MeanAnomalyRich_0110$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyRich_0110$data$anomStd_0110RS[
  MeanAnomalyRich_0110$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyRich_0110$data$anomStd_0110RS[
  MeanAnomalyRich_0110$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyRich_0110$data$anomStd_0110RS[
  MeanAnomalyRich_0110$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# rename variable to match model output
names(nd2)[1] <-"anomStd_0110RS"

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyRich_0110$model,data = nd2)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$anomStd_0110RS < QPV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$anomStd_0110RS > QPV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$anomStd_0110RS < QSV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$anomStd_0110RS > QSV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$anomStd_0110RS < QAL[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$anomStd_0110RS > QAL[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$anomStd_0110RS < QAH[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$anomStd_0110RS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd2$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


p2 <- ggplot(data = nd2, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd2$PredLower, ymax = nd2$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  theme_bw() + 
  #labs(fill = "% NH", col = "% NH") + 
  ylab("Species Richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-75, 150)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = c(0.3, 0.85)) +
  ggtitle("B. baseline: 1901-1910")




# 3.

nd3 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyRich_0120$data$anomStd_0120RS),
                        to = max(MeanAnomalyRich_0120$data$anomStd_0120RS),
                        length.out = 300),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyRich_0120$data$UI2)))

# back transform the predictors
nd3$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly_1901_1920)

# set richness and abundance to 0 - to be predicted
nd3$LogAbund <- 0
nd3$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd3$UI2=="Primary vegetation") & (nd3$StdTmeanAnomaly==min(abs(nd3$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyRich_0120$data$anomStd_0120RS[
  MeanAnomalyRich_0120$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyRich_0120$data$anomStd_0120RS[
  MeanAnomalyRich_0120$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyRich_0120$data$anomStd_0120RS[
  MeanAnomalyRich_0120$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyRich_0120$data$anomStd_0120RS[
  MeanAnomalyRich_0120$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# rename variable to match model output
names(nd3)[1] <-"anomStd_0120RS"

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyRich_0120$model,data = nd3)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$anomStd_0120RS < QPV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$anomStd_0120RS > QPV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$anomStd_0120RS < QSV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$anomStd_0120RS > QSV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$anomStd_0120RS < QAL[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$anomStd_0120RS > QAL[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$anomStd_0120RS < QAH[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$anomStd_0120RS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd3$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd3$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3$UI2 <- factor(nd3$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


p3 <- ggplot(data = nd3, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd3$PredLower, ymax = nd3$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  theme_bw() + 
  #labs(fill = "% NH", col = "% NH") + 
  ylab("Species Richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-75, 150)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = c(0.3, 0.85)) +
  ggtitle("C. baseline: 1901-1920")




# 4.

nd4 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyRich_0130$data$anomStd_0130RS),
                        to = max(MeanAnomalyRich_0130$data$anomStd_0130RS),
                        length.out = 300),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyRich_0130$data$UI2)))

# back transform the predictors
nd4$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd4$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly_1901_1930)

# set richness and abundance to 0 - to be predicted
nd4$LogAbund <- 0
nd4$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd4$UI2=="Primary vegetation") & (nd4$StdTmeanAnomaly==min(abs(nd4$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyRich_0130$data$anomStd_0130RS[
  MeanAnomalyRich_0130$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyRich_0130$data$anomStd_0130RS[
  MeanAnomalyRich_0130$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyRich_0130$data$anomStd_0130RS[
  MeanAnomalyRich_0130$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyRich_0130$data$anomStd_0130RS[
  MeanAnomalyRich_0130$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# rename variable to match model output
names(nd4)[1] <-"anomStd_0130RS"

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyRich_0130$model,data = nd4)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$anomStd_0130RS < QPV[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$anomStd_0130RS > QPV[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$anomStd_0130RS < QSV[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$anomStd_0130RS > QSV[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$anomStd_0130RS < QAL[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$anomStd_0130RS > QAL[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$anomStd_0130RS < QAH[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$anomStd_0130RS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd4$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd4$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd4$UI2 <- factor(nd4$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


p4 <- ggplot(data = nd4, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd4$PredLower, ymax = nd4$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  theme_bw() + 
  #labs(fill = "% NH", col = "% NH") + 
  ylab("Species Richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-75, 150)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = c(0.3, 0.85)) +
  ggtitle("D. baseline: 1901-1930")



# organise plots
cowplot::plot_grid(p1, p2, p3, p4)

ggsave(filename = paste0(outDir, "Extended_Data_Baselines_plots_Rich.pdf"), width = 10, height = 10, units = "in")





##%######################################################%##
#                                                          #
####             Testing anomaly including              ####
####          current temperature variability           ####
#                                                          #
##%######################################################%##


# need to create new anomaly vals based on current data
# usual temp different but divided by recent variabilty

# read in the data for the anomaly maps
load("3_PrepareClimateIndexMaps/Map_data_tempvars_2004_06.rdata") # temperatureVars

# create global points of each cell, extract the new standardising variable
tmp <- stack(paste0(dataDir,"cru_ts4.03.1901.2018.tmp.dat.nc"),varname="tmp")
# take the current predicts time data
tmp2004_6 <- tmp[[names(tmp)[1237:1272]]]
# which raster to use as present
pre_ras <- tmp2004_6
# get vector of positions in values of raster that are not NA
ras <- pre_ras[[1]]
wgs84 <- crs(tmp)
pnts <- rasterToPoints(ras, spatial = T)
SP <- SpatialPoints(pnts, proj4string=wgs84)

# load in current varibility data from WorldClim
# BIO4 = Temperature Seasonality (standard deviation ×100)

rdsdir <- "Z:/Datasets/wc2.1_10m_bio"

BIO4 <- raster(paste0(rdsdir,"/wc2.1_10m_bio_4.tif"))


# divide by 100
BIO4 <- BIO4/100

plot(BIO4)
# this is already the average variability across the years (I think)

# get the new data to the same resolution as the CRU data
bio4_agg <- aggregate(BIO4, fact = 3, fun = mean)

plot(bio4_agg)

# a lot of 0s which leads to Inf values, adding a low number to compensate
bio4_agg[bio4_agg == 0, ] <- 0.01


# extract the bio4 value for each terrestrial cell

temperatureVars$seas <- extract(bio4_agg, SP)

# recalculate the anomaly for each cell
temperatureVars$new_anom <- temperatureVars$Anom/temperatureVars$seas

SP_df <- as.data.frame(SP)

plot_data <- cbind(SP_df, temperatureVars)



###### maps of new anomaly ###### 

# remove any NAs
plot_data <- plot_data[!is.na(plot_data$new_anom), ]

# organise breaks, colours and labels
brks <- c(-2,-0.5,-0.2,-0.1,0,0.1,0.5,0.75,1,1.5,2,5,50)
cols <- c(rev(brewer.pal(n = 8,name = "Greens"))[4:8],
          (brewer.pal(n = 8,name = "Purples"))[4:6],
          (brewer.pal(n = 8,name = "Oranges"))[5:8])
labs <- c("-2 : -0.5","-0.5 : -0.2","-0.2 : -0.1","-0.1 : 0",
          "0 : 0.1","0.1 : 0.5","0.5 : 0.75","0.75 : 1","1 : 1.5","1.5 : 2","2 : 5","> 5")

# assign values into bins
plot_data$bins <- cut(plot_data$new_anom, 
                     breaks = brks, 
                     labels = labs,
                     include.lowest = TRUE)


# plot
ggplot(plot_data) + 
  geom_raster(aes(x = x, y = y, fill = bins), alpha = 0.9) +
  #scale_fill_viridis_c(option = "magma", values = c(0, 0.2, 1)) + 
  scale_fill_manual(values = cols) +
  xlab("") +
  ylab("") +
  labs(fill = "Standardised\nTemperature\nAnomaly") +
  theme_bw() +
  theme(legend.position = 'bottom', 
        #panel.border = element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_blank(),
        #legend.key.width = unit(3, "cm"),
        axis.ticks = element_blank(), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), 
        #legend.key.size = unit(0.2,"cm",
        panel.border = element_rect(size = 0.2),
        strip.background = element_rect(size = 0.2, fill = "transparent"),
        strip.text = element_text(size = 10))


# save plot as pdf
ggsave(filename = paste0(outDir, "/Extended_Data_anomaly_map_recentvar.pdf"), width = 6, height = 4)







#### Get the anomalies for each predicts site ####



# read in predicts data
predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSiteData.rds"))

wgs84 <- crs(tmp)

# Create spatial map of PREDICTS sites
predicts_sp <- SpatialPointsDataFrame(
  coords = cbind(predictsSites$Longitude, predictsSites$Latitude), 
  data = predictsSites, proj4string = wgs84)



# extract baseline mean and current variation
tmpSd_pres <- extract(bio4_agg, predicts_sp)


# get new anomaly values based on present variation anomaly
predictsSites$anomStd_new <- predictsSites$TmeanAnomaly / tmpSd_pres


## now prep for modelling

# centre the predictors
predictsSites$anomStd_newRS <- StdCenterPredictor(predictsSites$anomStd_new)



#### Run abundance models ####


# 1. Abundance, new mean anomaly
MeanAnomalyAbund_newvar <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(anomStd_newRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(anomStd_newRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyAbund_newvar, file = paste0(outDir, "/MeanAnomalyAbund_newvar.rdata"))


summary(MeanAnomalyAbund_newvar$model)
#LogAbund ~ UI2 + UI2:poly(anomStd_newRS, 1) + poly(anomStd_newRS, 1) + (1 | SS) + (1 | SSB)


# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)


nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyAbund_newvar$data$anomStd_newRS),
                        to = max(MeanAnomalyAbund_newvar$data$anomStd_newRS),
                        length.out = 200),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyAbund_newvar$data$UI2)))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$anomStd_new)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyAbund_newvar$data$anomStd_newRS[
  MeanAnomalyAbund_newvar$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyAbund_newvar$data$anomStd_newRS[
  MeanAnomalyAbund_newvar$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyAbund_newvar$data$anomStd_newRS[
  MeanAnomalyAbund_newvar$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyAbund_newvar$data$anomStd_newRS[
  MeanAnomalyAbund_newvar$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# rename variable to match model output
names(nd)[1] <-"anomStd_newRS"

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyAbund_newvar$model,data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$anomStd_newRS < QPV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$anomStd_newRS > QPV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$anomStd_newRS < QSV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$anomStd_newRS > QSV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$anomStd_newRS < QAL[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$anomStd_newRS > QAL[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$anomStd_newRS < QAH[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$anomStd_newRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


p1 <- ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  theme_bw() + 
  #labs(fill = "% NH", col = "% NH") + 
  ylab("Total Abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 3)) +
  ylim(c(-75, 170)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = c(0.3, 0.85)) +
  ggtitle("A.")







#### species richness models ####

# 1. SR, mean anomaly 1901-1905 baseline
MeanAnomalyRich_new <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
                                    fitFamily = "poisson",fixedFactors = "UI2",
                                    fixedTerms = list(anomStd_newRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(anomStd_newRS,1)"),
                                    saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyRich_new, file = paste0(outDir, "/MeanAnomalyRich_new.rdata"))

summary(MeanAnomalyRich_new$model)
# Species_richness ~ UI2 + poly(anomStd_newRS, 1) + UI2:poly(anomStd_newRS,      1) + (1 | SS) + (1 | SSB) + (1 | SSBS)



# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)


nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyRich_new$data$anomStd_newRS),
                        to = max(MeanAnomalyRich_new$data$anomStd_newRS),
                        length.out = 200),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyRich_new$data$UI2)))

# back transform the predictors
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$anomStd_new)

# set richness and abundance to 0 - to be predicted
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyRich_new$data$anomStd_newRS[
  MeanAnomalyRich_new$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyRich_new$data$anomStd_newRS[
  MeanAnomalyRich_new$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyRich_new$data$anomStd_newRS[
  MeanAnomalyRich_new$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyRich_new$data$anomStd_newRS[
  MeanAnomalyRich_new$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# rename variable to match model output
names(nd2)[1] <-"anomStd_newRS"

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyRich_new$model,data = nd2)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$anomStd_newRS < QPV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$anomStd_newRS > QPV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$anomStd_newRS < QSV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$anomStd_newRS > QSV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$anomStd_newRS < QAL[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$anomStd_newRS > QAL[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$anomStd_newRS < QAH[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$anomStd_newRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd2$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


p2 <- ggplot(data = nd2, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd2$PredLower, ymax = nd2$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  theme_bw() + 
  #labs(fill = "% NH", col = "% NH") + 
  ylab("Species Richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-75, 200)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = c(0.3, 0.85)) +
  ggtitle("B.")



# save plots in one pdf


plot_grid(p1, p2, ncol = 2)

ggsave(filename = paste0(outDir, "SupplementaryFig_Baselines_recentvar.pdf"), width = 10, height = 5, units = "in")
