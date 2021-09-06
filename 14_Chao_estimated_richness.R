##%######################################################%##
#                                                          #
#### Additional tests: Chao-estimated species richness  ####
#                                                          #
##%######################################################%##


rm(list = ls())

#library(devtools)
#install_github(repo = "timnewbold/predicts-demo",subdir = "predictsFunctions")
#install_github(repo = "timnewbold/StatisticalModels")

# load required libraries
library(predictsFunctions)
library(ggplot2)
library(StatisticalModels)
library(yarg)
source("Functions.R")

# directories
dataDir <- "0_data/"
outDir <- "14_Additional_Tests/"

# Set the path to your local copy of the database
predicts.path <- paste0(dataDir,"database.rds")

# Read in the PREDICTS data
predicts <- ReadPREDICTS(predicts.path)

# Select only data for insects
predicts <- predicts[(predicts$Class=="Insecta"),]

# hosts <- read.csv(paste0(dataDir,"HostSpecies_PREDICTS.csv"))
# predicts <- predicts[(predicts$Best_guess_binomial %in% hosts$Host.binomial),]

# Correct effort-sensitive abundance measures (assumes linear relationship between effort and recorded abundance)
predicts <- CorrectSamplingEffort(diversity = predicts)
# Correcting 870378 values for sensitivity to sampling effort (most of these are 0s, 19184 non-zero)

table(predicts$Diversity_metric)


# insects should not have diversity metric "percent cover", this is a mistake in the database
# remove those sites that are the problem
predicts <- predicts[!predicts$Diversity_metric == "percent cover", ]


# Merge sites that have the same coordinates (e.g. multiple traps on a single transect)
predicts <- predictsFunctions::MergeSites(diversity = predicts)


# remove rows where land use or use intensity info is missing
predicts.complete <- droplevels(predicts[(
  predicts$Predominant_land_use!="Cannot decide"),])
predicts.complete <- droplevels(predicts.complete[(
  predicts.complete$Use_intensity!="Cannot decide"),])

nrow(predicts.complete)
# 779,912 records

# Calculate site metrics of diversity
# this needs Tim's yarg package
sites <- SiteMetrics(diversity = predicts,
                     extra.cols = c("Predominant_land_use",
                                    "SSB","SSBS", "Biome"),
                     srEstimators = "Chao")
# 7567 rows

# Chao cannot be estimated for all sites, so drop those that it can't and see what is left
sites <- sites[!is.na(sites$ChaoR),] # 4782 rows

# simple plot of estimated species richness against sampled species richness
par(tck=-0.01,mgp=c(1.6,0.2,0),mar=c(2.7,2.7,0.2,0.2),las=1)
plot(sites$Species_richness+1,sites$ChaoR+1,log="xy",pch=16,
     xlab="Sampled species richness",ylab="Estimated species richness")
abline(0,1,lwd=2,col="#ff0000")


# histogram of the ratios of estimated to sampled species richness:
par(tck=-0.01,mgp=c(1.6,0.2,0),mar=c(2.7,2.7,0.2,0.2),las=1)
hist(x = log10((sites$ChaoR-sites$Species_richness)+1),
     xaxt="n",xlab="Estimated - Sampled richness",main=NULL)
axis(1,at=log10(c(0,1,5,10,100,1000)+1),labels=c(0,1,5,10,100,1000))


# Test whether completeness of sampling is related to sampled or estimated species richness:
par(mfrow=c(1,2),tck=-0.01,mgp=c(1.6,0.2,0),mar=c(2.7,2.7,0.2,0.2),las=1)
plot(x = sites$Species_richness+1,y = (sites$ChaoR - sites$Species_richness)+1,log="xy",pch=16,
     xlab="Sampled species richness",ylab="Estimated - Sampled richness",
     xaxt="n",yaxt="n")
axis(1,at=c(1,11,51,101,501),labels=c(0,10,50,100,500))
axis(2,at=c(1,11,51,101,501),labels=c(0,10,50,100,500))
plot(x = sites$ChaoR+1,y = (sites$ChaoR - sites$Species_richness)+1,log="xy",pch=16,
     xlab="Estimated species richness",ylab="Estimated - Sampled richness",
     xaxt="n",yaxt="n")
axis(1,at=c(1,11,51,101,501,1001),labels=c(0,10,50,100,500,1000))
axis(2,at=c(1,11,51,101,501),labels=c(0,10,50,100,500))

# round the estimated species richness values to integers.
sites$ChaoR <- round(sites$ChaoR,0)

# First, we will rearrange the land-use classification a bit
sites$LandUse <- paste(sites$Predominant_land_use)

# Drop classification where land use could not be identified
sites$LandUse[(sites$LandUse=="Cannot decide")] <- NA

# Drop classification where the stage of recovery of secondary vegetation is unknown
# sites$LandUse[(sites$LandUse=="Secondary vegetation (indeterminate age)")] <- NA

# Now make the variable a factor, and set the reference level to primary vegetation
sites$LandUse <- factor(sites$LandUse)
sites$LandUse <- relevel(sites$LandUse,ref="Primary vegetation")

sites$Use_intensity[sites$Use_intensity=="Cannot decide"] <- NA

# combine LU and UI 
sites$UI <- paste0(sites$LandUse,'_',sites$Use_intensity)
sites$UI[grep("NA",sites$UI)] <- NA

# recode according to land use and use intensity combinations
sites$UI2 <- dplyr::recode(sites$UI,
                           'Primary vegetation_Minimal use' = 'Primary vegetation',
                           'Cropland_Light use' = 'Agriculture_High',
                           'Secondary vegetation (indeterminate age)_Minimal use' = 'Secondary vegetation',
                           'Urban_Light use' = 'Urban',
                           'Secondary vegetation (indeterminate age)_Light use' = 'Secondary vegetation',
                           'Cropland_Intense use' = 'Agriculture_High',
                           'Cropland_Minimal use' = 'Agriculture_Low',
                           'Pasture_Light use' = 'Agriculture_Low',
                           'Pasture_Minimal use' = 'Agriculture_Low',
                           'Intermediate secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Mature secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Secondary vegetation (indeterminate age)_Intense use' = 'Secondary vegetation',
                           'Pasture_Intense use' = 'Agriculture_High',
                           'Urban_Minimal use' = 'Urban',
                           'Primary vegetation_Light use' = 'Primary vegetation',
                           'Young secondary vegetation_Light use' = 'Secondary vegetation',
                           'Urban_Intense use' = 'Urban',
                           'Primary vegetation_Intense use' = 'Primary vegetation',
                           'Young secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Mature secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Plantation forest_Minimal use' = 'Agriculture_Low',
                           'Plantation forest_Intense use' = 'Agriculture_High',
                           'Young secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Plantation forest_Light use' = 'Agriculture_High',
                           'Mature secondary vegetation_Light use' = 'Secondary vegetation',
                           'Intermediate secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Intermediate secondary vegetation_Light use' = 'Secondary vegetation')

# 
sites$Use_intensity[((sites$LandUse=="Mature secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"
sites$Use_intensity[((sites$LandUse=="Intermediate secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"
sites$Use_intensity[((sites$LandUse=="Young secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"

# remove the urban sites and NA in UI2
sites <- sites[!sites$UI2 == "Urban", ]
sites <- sites[!is.na(sites$UI2), ]


sites <- droplevels(sites)

# Remove sites without coordinates
sites <- sites[!is.na(sites$Latitude), ]

# nsites with chao richness estimated 4268


# load in original dataset and match climate anomaly data to sites
predictsSites <- readRDS("5_PREDICTSMatchPropNatHab/PREDICTSSitesWithClimateAndNatHab.rds")
predictsSites <- predictsSites@data

predictsSites <- predictsSites[ , c("SSBS", "StdTmeanAnomaly", "StdTmaxAnomaly")]

# combine to get the anomaly info
sites2 <- merge(sites, predictsSites, by = "SSBS")


sites2$UI2 <- factor(sites2$UI2)
sites2$UI2 <- relevel(sites2$UI2,ref="Primary vegetation")

# organise the climate anomaly data
sites2$StdTmeanAnomalyRS <- StdCenterPredictor(sites2$StdTmeanAnomaly)

# rescale the variable
sites2$StdTmaxAnomalyRS <- StdCenterPredictor(sites2$StdTmaxAnomaly)


# charlie added this line as later bits were throwing errors
sites2 <- droplevels(sites2)


save(sites2, file = paste0(outDir, "PREDICTS_sites_ChaoR.rdata"))

#### run the model with ChaoR as the response ####

# 1. Chao Richness, mean anomaly
MeanAnomalyModelChaoR <- GLMERSelect(modelData = sites2,responseVar = "ChaoR",
                                    fitFamily = "poisson", fixedFactors = "UI2",
                                    fixedTerms = list(StdTmeanAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                    saveVars = c("Total_abundance", "SSBS"))

save(MeanAnomalyModelChaoR, file = paste0(outDir, "/MeanAnomalyModelChaoR.rdata"))


summary(MeanAnomalyModelChaoR$model)
# ChaoR ~ UI2 + poly(StdTmeanAnomalyRS, 1) + UI2:poly(StdTmeanAnomalyRS,      1) + (1 | SS) + (1 | SSB) + (1 | SSBS)


# 1. Chao Richness, max anomaly
MaxAnomalyModelChaoR <- GLMERSelect(modelData = sites2,responseVar = "ChaoR",
                                     fitFamily = "poisson", fixedFactors = "UI2",
                                     fixedTerms = list(StdTmaxAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                     fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)"),
                                     saveVars = c("Total_abundance", "SSBS"))


save(MaxAnomalyModelChaoR, file = paste0(outDir, "/MaxAnomalyModelChaoR.rdata"))

summary(MaxAnomalyModelChaoR$model)
# ChaoR ~ UI2 + poly(StdTmaxAnomalyRS, 1) + UI2:poly(StdTmaxAnomalyRS,      1) + (1 | SS) + (1 | SSB) + (1 | SSBS)


#### plot the results ####

#load(paste0(outDir, "MeanAnomalyModelChaoR.rdata"))
#load(paste0(outDir, "MaxAnomalyModelChaoR.rdata"))


#### 1. mean anom ChaoR ####

# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)

nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelChaoR$data$UI2)))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = sites2$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$ChaoR <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyModelChaoR$data$StdTmeanAnomalyRS[
  MeanAnomalyModelChaoR$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelChaoR$model,data = nd)

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

nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))


p1 <- ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  theme_bw() + 
  ylab("Change in Chao estimated richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.25, 2)) +
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"), 
        legend.position = "none") +
  ggtitle("a.")


#### maximum anomaly ####

nd2 <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelChaoR$data$StdTmaxAnomalyRS),
                        to = max(MaxAnomalyModelChaoR$data$StdTmaxAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxAnomalyModelChaoR$data$UI2)))

# back transform the predictors
nd2$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmaxAnomalyRS,
  originalX = sites2$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd2$ChaoR <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomaly==min(abs(nd2$StdTmaxAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MaxAnomalyModelChaoR$data$StdTmaxAnomalyRS[
  MaxAnomalyModelChaoR$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxAnomalyModelChaoR$data$StdTmaxAnomalyRS[
  MaxAnomalyModelChaoR$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxAnomalyModelChaoR$data$StdTmaxAnomalyRS[
  MaxAnomalyModelChaoR$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxAnomalyModelChaoR$data$StdTmaxAnomalyRS[
  MaxAnomalyModelChaoR$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MaxAnomalyModelChaoR$model,data = nd2)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmaxAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmaxAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmaxAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmaxAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmaxAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmaxAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd2$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))


p2 <- ggplot(data = nd2, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd2$PredLower, ymax = nd2$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
  theme_bw() + 
  ylab("Change in Chao estimated richness (%)") +
  xlab("Standardised Maximum Temperature Anomaly") +
  xlim(c(-0.25, 2)) +
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"), 
        legend.position = c(0.3, 0.8), legend.title = element_blank()) +
  ggtitle("b.")



library(cowplot)

plot_grid(p1, p2)

ggsave2(filename = paste0(outDir, "Chao_Richness_LUCC_plots.pdf"), width = 8, height = 4)


#### create a map of the sites ####


# plot the raster in ggplot
map.world <- map_data('world')

# map of sites
ggplot() +
  geom_map(data=map.world, map=map.world,
           aes(x=long, y=lat, group=group, map_id=region),
           fill= "grey", colour="grey", size=0.2) +
  geom_point(data = sites2, aes(x = Longitude, y = Latitude), col = c("#1E90FF"), fill = c("#104E8B"), shape = 21) +
  theme(axis.title = element_blank(), 
        plot.background = element_blank(), 
        panel.background = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank()) 

# take a look at how many in each realm
sites2$Tropical <- NA

sites2[sites2$Latitude > -23.44 & sites2$Latitude < 23.44, 'Tropical'] <- "Tropical"
sites2[is.na(sites2$Tropical), 'Tropical'] <- "Temperate"

table(sites2$Tropical)

#Temperate  Tropical 
#3287       981
