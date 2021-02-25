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
tmp1901_1930 <- tmp[[names(tmp)[1:360]]]
tmp1901_1930mean <- calc(tmp1901_1930, base::mean)
tmp1901_1930sd <- calc(tmp1901_1930, stats::sd)


# alternative baselines
# 5 years
tmp1901_1905 <- tmp[[names(tmp)[1:60]]]
tmp1901_1905mean <- calc(tmp1901_1905, base::mean)
tmp1901_1905sd <- calc(tmp1901_1905, stats::sd)

# 10 years
tmp1901_1910 <- tmp[[names(tmp)[1:120]]]
tmp1901_1910mean <- calc(tmp1901_1910, base::mean)
tmp1901_1910sd <- calc(tmp1901_1910, stats::sd)

# 20 years
tmp1901_1920 <- tmp[[names(tmp)[1:240]]]
tmp1901_1920mean <- calc(tmp1901_1920, base::mean)
tmp1901_1920sd <- calc(tmp1901_1920, stats::sd)



### Calculate the standardised anomaly ###


# get the present day values
tmp2004_6 <- tmp[[names(tmp)[1237:1272]]]

# calc the mean for present time period
tmp2004_6mean <- calc(tmp[[names(tmp)[1237:1272]]], base::mean)


# standardise the baselines
anom_01_30 <- (calc(tmp2004_6, base::mean)-tmp1901_1930mean)  / tmp1901_1930sd

anom_01_05 <- (calc(tmp2004_6, base::mean)-tmp1901_1905mean)  / tmp1901_1905sd
anom_01_10 <- (calc(tmp2004_6, base::mean)-tmp1901_1910mean)  / tmp1901_1910sd
anom_01_20 <- (calc(tmp2004_6, base::mean)-tmp1901_1920mean)  / tmp1901_1920sd


###### maps of each anomaly ###### 
plot_data_0130 <- as.data.frame(anom_01_30, xy = TRUE)
plot_data_0105 <- as.data.frame(anom_01_05, xy = TRUE)
plot_data_0110 <- as.data.frame(anom_01_10, xy = TRUE)
plot_data_0120 <- as.data.frame(anom_01_20, xy = TRUE)


# add details to data table
plot_data_0130$year <- "1901-1930"
plot_data_0105$year <- "1901-1905"
plot_data_0110$year <- "1901-1910"
plot_data_0120$year <- "1901-1920"

# combine the two together
all_plot <- rbind(plot_data_0130, plot_data_0105, plot_data_0110, plot_data_0120)

all_plot <- all_plot[!is.na(all_plot$layer), ]

# organise breaks, colours and labels
brks <- c(-2,-0.5,-0.2,-0.1,0,0.1,0.5,0.75,1,1.5,2,5,50)
cols <- c(rev(brewer.pal(n = 8,name = "Greens"))[4:8],
          (brewer.pal(n = 8,name = "Purples"))[4:6],
          (brewer.pal(n = 8,name = "Oranges"))[5:8])
labs <- c("-2 : -0.5","-0.5 : -0.2","-0.2 : -0.1","-0.1 : 0",
          "0 : 0.1","0.1 : 0.5","0.5 : 0.75","0.75 : 1","1 : 1.5","1.5 : 2","2 : 5","> 5")

# assign values into bins
all_plot$bins <- cut(all_plot$layer, 
                     breaks = brks, 
                     labels = labs,
                     include.lowest = TRUE)


year_labs <- c("D. 1901-1930", "A. 1901-1905", "B. 1901-1910", "C. 1901-1920")
names(year_labs) <- unique(all_plot$year)

all_plot$year <- as.factor(all_plot$year)

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
ggsave(filename = paste0(outDir, "/Extended_Data_anomaly_maps.pdf"), width = 12, height = 8)





##%######################################################%##
#                                                          #
####    Run models using new baselines for SCA         ####
#                                                          #
##%######################################################%##

# run the climate:LU model for abundance with each new baseline to compare results

# read in predicts data
predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSiteData.rds"))

wgs84 <- crs(tmp)

# Create spatial map of PREDICTS sites
predicts_sp <- SpatialPointsDataFrame(
  coords = cbind(predictsSites$Longitude, predictsSites$Latitude), 
  data = predictsSites, proj4string = wgs84)


# extract mean and sd from baseline datasets for each predicts site
tmpMean_0105 <- extract(tmp1901_1905mean, predicts_sp)
tmpSd_0105 <- extract(tmp1901_1905sd, predicts_sp)

tmpMean_0110 <- extract(tmp1901_1910mean, predicts_sp)
tmpSd_0110 <- extract(tmp1901_1910sd, predicts_sp)

tmpMean_0120 <- extract(tmp1901_1920mean, predicts_sp)
tmpSd_0120 <- extract(tmp1901_1920sd, predicts_sp)


# calc anomaly using the avgtmp column in predictsSites
predictsSites$anomStd_0105 <- (predictsSites$avg_temp - tmpMean_0105) / tmpSd_0105
predictsSites$anomStd_0110 <- (predictsSites$avg_temp - tmpMean_0110) / tmpSd_0110
predictsSites$anomStd_0120 <- (predictsSites$avg_temp - tmpMean_0120) / tmpSd_0120
  
  
  
# centre the predictors
predictsSites$anomStd_0105RS <- StdCenterPredictor(predictsSites$anomStd_0105)
predictsSites$anomStd_0110RS <- StdCenterPredictor(predictsSites$anomStd_0110)
predictsSites$anomStd_0120RS <- StdCenterPredictor(predictsSites$anomStd_0120)



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




# 2. Abundance, mean anomaly 1901-1910 baseline
MeanAnomalyAbund_0110 <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(anomStd_0110RS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(anomStd_0110RS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyAbund_0110, file = paste0(outDir, "/MeanAnomalyAbund_0110.rdata"))


# 3. Abundance, mean anomaly  1901-1920 baseline
MeanAnomalyAbund_0120 <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(anomStd_0120RS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(anomStd_0120RS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyAbund_0120, file = paste0(outDir, "/MeanAnomalyAbund_0120.rdata"))



# take a look at model outputs
summary(MeanAnomalyAbund_0105$model)
summary(MeanAnomalyAbund_0110$model)
summary(MeanAnomalyAbund_0120$model)



######## plot figures ######## 
  
#load(paste0(outDir, "/MeanAnomalyAbund_0105.rdata"))
#load(paste0(outDir, "/MeanAnomalyAbund_0110.rdata"))
#load(paste0(outDir, "/MeanAnomalyAbund_0120.rdata"))




# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)


nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyAbund_0105$data$anomStd_0105RS),
                        to = max(MeanAnomalyAbund_0105$data$anomStd_0105RS),
                        length.out = 300),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyAbund_0105$data$UI2)))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$anomStd_0105)

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
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  theme_bw() + 
  #labs(fill = "% NH", col = "% NH") + 
  ylab("Total Abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-75, 120)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = c(0.3, 0.85)) +
  ggtitle("A. baseline: 1901-1905")



# 2.

nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyAbund_0110$data$anomStd_0110RS),
                        to = max(MeanAnomalyAbund_0110$data$anomStd_0110RS),
                        length.out = 400),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyAbund_0110$data$UI2)))

# back transform the predictors
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$anomStd_0110)

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
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  theme_bw() + 
  #labs(fill = "% NH", col = "% NH") + 
  ylab("Total Abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-75, 120)) + 
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
  originalX = predictsSites$anomStd_0120)

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
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  theme_bw() + 
  #labs(fill = "% NH", col = "% NH") + 
  ylab("Total Abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-75, 120)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        legend.title = element_blank(), 
        legend.position = c(0.3, 0.85)) +
  ggtitle("C. baseline: 1901-1920")


# organise plots
cowplot::plot_grid(p1, p2, p3)

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



# take a look at model outputs
summary(MeanAnomalyRich_0105$model)
summary(MeanAnomalyRich_0110$model)
summary(MeanAnomalyRich_0120$model)



######## plot figures ######## 

#load(paste0(outDir, "/MeanAnomalyRich_0105.rdata"))
#load(paste0(outDir, "/MeanAnomalyRich_0110.rdata"))
#load(paste0(outDir, "/MeanAnomalyRich_0120.rdata"))




# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)


nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyRich_0105$data$anomStd_0105RS),
                        to = max(MeanAnomalyRich_0105$data$anomStd_0105RS),
                        length.out = 300),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyRich_0105$data$UI2)))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$anomStd_0105)

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
                        length.out = 400),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyRich_0110$data$UI2)))

# back transform the predictors
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$anomStd_0110)

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
  originalX = predictsSites$anomStd_0120)

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


# organise plots
cowplot::plot_grid(p1, p2, p3)

ggsave(filename = paste0(outDir, "Extended_Data_Baselines_plots_Rich.pdf"), width = 10, height = 10, units = "in")





##%######################################################%##
#                                                          #
####             Testing anomaly including              ####
####          current temperature variability           ####
#                                                          #
##%######################################################%##


# need to create new anomaly vals based on current data
# usual temp different but divided by recent variabilty


# get the CRU data
tmp <- stack(paste0(dataDir,"cru_ts4.03.1901.2018.tmp.dat.nc"),varname="tmp")

# currently used baseline, just get the mean
# take names of values for 1901 to 1930
tmp1901_1930 <- tmp[[names(tmp)[1:360]]]
tmp1901_1930mean <- calc(tmp1901_1930, base::mean)


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

# calculate the new anomaly for plotting only
dif <- (calc(tmp2004_6, base::mean)-tmp1901_1930mean)
plot(dif)

anom_new <-  dif / bio4_agg

plot(anom_new)

# a very few very large values for the anomaly, so remove anything greater than 10
anom_new[anom_new > 10] <- NA

plot(anom_new)


###### maps of each anomaly ###### 
plot_data <- as.data.frame(anom_new, xy = TRUE)

# remove any NAs
plot_data <- plot_data[!is.na(plot_data$layer), ]

# organise breaks, colours and labels
brks <- c(-2,-0.5,-0.2,-0.1,0,0.1,0.5,0.75,1,1.5,2,5,50)
cols <- c(rev(brewer.pal(n = 8,name = "Greens"))[4:8],
          (brewer.pal(n = 8,name = "Purples"))[4:6],
          (brewer.pal(n = 8,name = "Oranges"))[5:8])
labs <- c("-2 : -0.5","-0.5 : -0.2","-0.2 : -0.1","-0.1 : 0",
          "0 : 0.1","0.1 : 0.5","0.5 : 0.75","0.75 : 1","1 : 1.5","1.5 : 2","2 : 5","> 5")

# assign values into bins
plot_data$bins <- cut(plot_data$layer, 
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
tmpMean_0130 <- extract(tmp1901_1930mean, predicts_sp)
tmpSd_pres <- extract(bio4_agg, predicts_sp)


# get new anomaly values based on present variation anomaly
predictsSites$anomStd_new <- (predictsSites$avg_temp - tmpMean_0130) / tmpSd_pres


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
                        length.out = 300),
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
                        length.out = 300),
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

ggsave(filename = paste0(outDir, "Extended_Data_Baselines_recentvar.pdf"), width = 10, height = 5, units = "in")
