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
source("Functions.R")

# directories
dataDir <- "0_data/"
outDir <- "10_SCA_Baseline_testing/"
dir.create(outDir)
predictsDataDir <- "6_RunLUClimateModels/"


# extract data for dif baseline periods

# load in the mean temperature data from CRU
tmp <- stack(paste0(dataDir,"cru_ts4.03.1901.2018.tmp.dat.nc"),varname="tmp")

# currently used baseline
# take names of values for 1901 to 1905
tmp1901_1905 <- tmp[[names(tmp)[1:60]]]
tmp1901_1905mean <- calc(tmp1901_1905, base::mean)
tmp1901_1905sd <- calc(tmp1901_1905, stats::sd)

# take names of values for additional time frames
tmp1901_1920 <- tmp[[names(tmp)[1:240]]]
tmp1920_1925 <- tmp[[names(tmp)[229:300]]]
tmp1940_1945 <- tmp[[names(tmp)[469:540]]]


# calculate the mean and sd of the baseline values
tmp1901_1920mean <- calc(tmp1901_1920, base::mean)
tmp1901_1920sd <- calc(tmp1901_1920, stats::sd)

tmp1920_1925mean <- calc(tmp1920_1925, base::mean)
tmp1920_1925sd <- calc(tmp1920_1925, stats::sd)

tmp1940_1945mean <- calc(tmp1940_1945, base::mean)
tmp1940_1945sd <- calc(tmp1940_1945, stats::sd)



### Calculate the standardised anomaly ###


# get the present day values
tmp2004_6 <- tmp[[names(tmp)[1237:1272]]]

# calc the mean for present time period
tmp2004_6mean <- calc(tmp[[names(tmp)[1237:1272]]], base::mean)


# standardise the baseline
anom_01_05 <- (calc(tmp2004_6, base::mean)-tmp1901_1905mean)  / tmp1901_1905sd
anom_01_20 <- (calc(tmp2004_6, base::mean)-tmp1901_1920mean)  / tmp1901_1920sd
anom_20_25 <- (calc(tmp2004_6, base::mean)-tmp1920_1925mean)  / tmp1920_1925sd
anom_40_45 <- (calc(tmp2004_6, base::mean)-tmp1940_1945mean)  / tmp1940_1945sd


###### maps of each anomaly ###### 
plot_data_0105 <- as.data.frame(anom_01_05, xy = TRUE)
plot_data_0120 <- as.data.frame(anom_01_20, xy = TRUE)
plot_data_2025 <- as.data.frame(anom_20_25, xy = TRUE)
plot_data_4045 <- as.data.frame(anom_40_45, xy = TRUE)


# add details to data table
plot_data_0105$year <- "1901-1905"
plot_data_0120$year <- "1901-1920"
plot_data_2025$year <- "1920-1925"
plot_data_4045$year <- "1940-1945"

# combine the two together
all_plot <- rbind(plot_data_0105, plot_data_0120, plot_data_2025, plot_data_4045)

all_plot <- all_plot[!is.na(all_plot$layer), ]

# organise breaks, colours and labels
brks <- c(-1,-0.5,-0.2,-0.1,0,0.1,0.5,0.75,1,1.5,2,5,50)
cols <- c(rev(brewer.pal(n = 8,name = "Greens"))[4:8],
          (brewer.pal(n = 8,name = "Purples"))[4:6],
          (brewer.pal(n = 8,name = "Oranges"))[5:8])
labs <- c("-1 : -0.5","-0.5 : -0.2","-0.2 : -0.1","-0.1 : 0",
          "0 : 0.1","0.1 : 0.5","0.5 : 0.75","0.75 : 1","1 : 1.5","1.5 : 2","2 : 5","> 5")

# assign values into bins
all_plot$bins <- cut(all_plot$layer, 
                     breaks = brks, 
                     labels = labs,
                     include.lowest = TRUE)

# plot
ggplot(all_plot) + 
  geom_raster(aes(x = x, y = y, fill = bins), alpha = 0.9) +
  #scale_fill_viridis_c(option = "magma", values = c(0, 0.2, 1)) + 
  scale_fill_manual(values = cols) +
  facet_wrap(~ year, nrow = 2) +
  xlab("") +
  ylab("") +
  labs(fill = "Standardised\nClimate\nAnomaly") +
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
ggsave(filename = paste0(outDir, "/Extended_Data_anomaly_maps.pdf"), width = 12, height = 8)





##%######################################################%##
#                                                          #
###     Run models using new baselines for SCA           ###
#                                                          #
##%######################################################%##

# run the clumate:LU model for abundance with each new baseline to compare results

# read in predicts data
predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSiteData.rds"))

wgs84 <- crs(tmp)

# Create spatial map of PREDICTS sites
predicts_sp <- SpatialPointsDataFrame(
  coords = cbind(predictsSites$Longitude, predictsSites$Latitude), 
  data = predictsSites, proj4string = wgs84)


# extract mean and sd from baseline datasets for each predicts site
tmpMean_0120 <- extract(tmp1901_1920mean, predicts_sp)
tmpSd_0120 <- extract(tmp1901_1920sd, predicts_sp)

tmpMean_2025 <- extract(tmp1920_1925mean, predicts_sp)
tmpSd_2025 <- extract(tmp1920_1925sd, predicts_sp)

tmpMean_4045 <- extract(tmp1940_1945mean, predicts_sp)
tmpSd_4045 <- extract(tmp1940_1945sd, predicts_sp)


# calc anomaly using the avgtmp column in predictsSites
predictsSites$anomStd_0120 <- (predictsSites$avg_temp - tmpMean_0120) / tmpSd_0120
predictsSites$anomStd_2025 <- (predictsSites$avg_temp - tmpMean_2025) / tmpSd_2025
predictsSites$anomStd_4045 <- (predictsSites$avg_temp - tmpMean_4045) / tmpSd_4045
  
  
  
# centre the predictors
predictsSites$anomStd_0120RS <- StdCenterPredictor(predictsSites$anomStd_0120)
predictsSites$anomStd_2025RS <- StdCenterPredictor(predictsSites$anomStd_2025)
predictsSites$anomStd_4045RS <- StdCenterPredictor(predictsSites$anomStd_4045)



# run abundance model using each of the anomalies

# 1. Abundance, mean anomaly 1901-1920 baseline
MeanAnomalyAbund_0120 <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(anomStd_0120RS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(anomStd_0120RS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyAbund_0120, file = paste0(outDir, "/MeanAnomalyAbund_0120.rdata"))




# 2. Abundance, mean anomaly 1920-1925 baseline
MeanAnomalyAbund_2025 <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(anomStd_2025RS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(anomStd_2025RS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyAbund_2025, file = paste0(outDir, "/MeanAnomalyAbund_2025.rdata"))


# 3. Abundance, mean anomaly  1940-1941 baseline
MeanAnomalyAbund_4045 <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(anomStd_4045RS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(anomStd_4045RS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyAbund_4045, file = paste0(outDir, "/MeanAnomalyAbund_4045.rdata"))



# take a look at model outputs
summary(MeanAnomalyAbund_0120$model)
summary(MeanAnomalyAbund_2025$model)
summary(MeanAnomalyAbund_4045$model)



######## plot figures ######## 
  
# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)


nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyAbund_0120$data$anomStd_0120RS),
                        to = max(MeanAnomalyAbund_0120$data$anomStd_0120RS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyAbund_0120$data$UI2)))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$anomStd_0120)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))

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
names(nd)[1] <-"anomStd_0120RS"

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyAbund_0120$model,data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$anomStd_0120RS < QPV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$anomStd_0120RS > QPV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$anomStd_0120RS < QSV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$anomStd_0120RS > QSV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$anomStd_0120RS < QAL[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$anomStd_0120RS > QAL[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$anomStd_0120RS < QAH[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$anomStd_0120RS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


p1 <- ggplot(data = nd, aes(x = anomStd_0120RS, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  theme_bw() + 
  #labs(fill = "% NH", col = "% NH") + 
  ylab("Total Abundance (%)") +
  xlab("Standardised Climate Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-75, 100)) + 
  theme(aspect.ratio = 1, text = element_text(size = 9),
        legend.title = element_blank(), 
        legend.position = c(0.2, 0.85)) +
  ggtitle("baseline: 1901-1920")



# 2.

nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyAbund_2025$data$anomStd_2025RS),
                        to = max(MeanAnomalyAbund_2025$data$anomStd_2025RS),
                        length.out = 400),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyAbund_2025$data$UI2)))

# back transform the predictors
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$anomStd_2025)

# set richness and abundance to 0 - to be predicted
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyAbund_2025$data$anomStd_2025RS[
  MeanAnomalyAbund_2025$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyAbund_2025$data$anomStd_2025RS[
  MeanAnomalyAbund_2025$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyAbund_2025$data$anomStd_2025RS[
  MeanAnomalyAbund_2025$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyAbund_2025$data$anomStd_2025RS[
  MeanAnomalyAbund_2025$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# rename variable to match model output
names(nd2)[1] <-"anomStd_2025RS"

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyAbund_2025$model,data = nd2)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$anomStd_2025RS < QPV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$anomStd_2025RS > QPV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$anomStd_2025RS < QSV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$anomStd_2025RS > QSV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$anomStd_2025RS < QAL[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$anomStd_2025RS > QAL[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$anomStd_2025RS < QAH[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$anomStd_2025RS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd2$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


p2 <- ggplot(data = nd2, aes(x = anomStd_2025RS, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd2$PredLower, ymax = nd2$PredUpper, fill = UI2), alpha = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  theme_bw() + 
  #labs(fill = "% NH", col = "% NH") + 
  ylab("Total Abundance (%)") +
  xlab("Standardised Climate Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-75, 100)) + 
  theme(aspect.ratio = 1, text = element_text(size = 9),
        legend.title = element_blank(), 
        legend.position = c(0.2, 0.85)) +
  ggtitle("baseline: 1920-1925")




# 3.

nd3 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyAbund_4045$data$anomStd_4045RS),
                        to = max(MeanAnomalyAbund_4045$data$anomStd_4045RS),
                        length.out = 300),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyAbund_4045$data$UI2)))

# back transform the predictors
nd3$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3$StdTmeanAnomalyRS,
  originalX = predictsSites$anomStd_4045)

# set richness and abundance to 0 - to be predicted
nd3$LogAbund <- 0
nd3$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd3$UI2=="Primary vegetation") & (nd3$StdTmeanAnomaly==min(abs(nd3$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyAbund_4045$data$anomStd_4045RS[
  MeanAnomalyAbund_4045$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyAbund_4045$data$anomStd_4045RS[
  MeanAnomalyAbund_4045$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyAbund_4045$data$anomStd_4045RS[
  MeanAnomalyAbund_4045$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyAbund_4045$data$anomStd_4045RS[
  MeanAnomalyAbund_4045$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# rename variable to match model output
names(nd3)[1] <-"anomStd_4045RS"

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyAbund_4045$model,data = nd3)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$anomStd_4045RS < QPV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$anomStd_4045RS > QPV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$anomStd_4045RS < QSV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$anomStd_4045RS > QSV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$anomStd_4045RS < QAL[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$anomStd_4045RS > QAL[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$anomStd_4045RS < QAH[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$anomStd_4045RS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd3$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd3$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd3$UI2 <- factor(nd3$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


p3 <- ggplot(data = nd3, aes(x = anomStd_4045RS, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd3$PredLower, ymax = nd3$PredUpper, fill = UI2), alpha = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  theme_bw() + 
  #labs(fill = "% NH", col = "% NH") + 
  ylab("Total Abundance (%)") +
  xlab("Standardised Climate Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-75, 100)) + 
  theme(aspect.ratio = 1, text = element_text(size = 9),
        legend.title = element_blank(), 
        legend.position = c(0.2, 0.85)) +
  ggtitle("baseline: 1940-1945")


# organise plots
plot_grid(p1, p2, p3)

ggsave(filename = paste0(outDir, "Extended_Data_Baselines_plots.pdf"), width = 10, height = 10, units = "in")






