##%######################################################%##
#                                                          #
####    Investigating "habitat as described" column     ####
#                                                          #
##%######################################################%##

# one reviewer was concerned about the classification of primary vegetation, investigating here. 

# load required libraries
library(predictsFunctions)
library(StatisticalModels)
library(ggplot2)
source("Functions.R")

# directories
dataDir <- "0_data/"
outDir <- "14_Additional_Tests/"
dir.create(outDir)

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
predicts <- MergeSites(diversity = predicts)


# remove rows where land use or use intensity info is missing
predicts.complete <- droplevels(predicts[(
  predicts$Predominant_land_use!="Cannot decide"),])
predicts.complete <- droplevels(predicts.complete[(
  predicts.complete$Use_intensity!="Cannot decide"),])

nrow(predicts.complete)
# 779,912 records

# get counts of n species in major groups
species <- unique(predicts.complete[,c('Order','Taxon_name_entered')])

order.counts <- tapply(X = species$Taxon_name_entered,
                       INDEX = species$Order,
                       FUN = function(sp) length(unique(sp)))



# Calculate site metrics of diversity
sites <- SiteMetrics(diversity = predicts,
                     extra.cols = c("Predominant_land_use",
                                    "SSB","SSBS", "Habitat_as_described"))

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

# transform abundance values 
sites$LogAbund <- log(sites$Total_abundance+1)


# Remove sites without coordinates
sites <- sites[!is.na(sites$Latitude), ]


##### investigate "habitat as described" for the primary vegetation sites #####

sites_prim <- sites[sites$UI2 == "Primary vegetation", ] # 1516 sites

# take a look and make a list of the ones that look a bit weird
View(unique(sites_prim$Habitat_as_described))

View(sites[sites$Habitat_as_described == "Traditional Bedouin agricultural garden, diversity of fruit and vegetable crops, irrigated, no herbicide or pesticide, goat manure for fertiliser", ])

# make a list of potentially dodgy descriptions inc site ids

query <- NULL

query <- rbind(query, sites[sites$Habitat_as_described == "Urban Park", c("SSBS", "Habitat_as_described")])
query <- rbind(query, sites[sites$Habitat_as_described == "Savanna farmland", c("SSBS", "Habitat_as_described")])
query <- rbind(query, sites[sites$Habitat_as_described == "Healthland/Park, Rural", c("SSBS", "Habitat_as_described")])
query <- rbind(query, sites[sites$Habitat_as_described == "Canal Side, Rural with Agriculture", c("SSBS", "Habitat_as_described")]) # potential reference site in this case
query <- rbind(query, sites[sites$Habitat_as_described == "beside road, people present", c("SSBS", "Habitat_as_described")]) # potential reference site in this case
query <- rbind(query, sites[sites$Habitat_as_described == "landslide", c("SSBS", "Habitat_as_described")]) # potential reference site in this case
query <- rbind(query, sites[sites$Habitat_as_described == "Mainly agriculture with some forest within 2Km of site", c("SSBS", "Habitat_as_described")]) # potential reference site in this case
query <- rbind(query, sites[sites$Habitat_as_described == "Mainly agriculure, some forest and urban development within 2Km of site", c("SSBS", "Habitat_as_described")]) # potential reference site in this case
query <- rbind(query, sites[sites$Habitat_as_described == "Mostly agriculture within 2Km of site", c("SSBS", "Habitat_as_described")]) # potential reference site in this case
query <- rbind(query, sites[sites$Habitat_as_described == "on road bank", c("SSBS", "Habitat_as_described")]) # potential reference site in this case
query <- rbind(query, sites[sites$Habitat_as_described == "Roadside vegetation", c("SSBS", "Habitat_as_described")]) # potential reference site in this case
query <- rbind(query, sites[sites$Habitat_as_described == "Traditional Bedouin agricultural garden, diversity of fruit and vegetable crops, irrigated, no herbicide or pesticide, goat manure for fertiliser", c("SSBS", "Habitat_as_described")]) # potential reference site in this case


# 95 sites with possibly dodgy assignment of primary vegetation

table(sites$UI2)

# Agriculture_High      Agriculture_Low   Primary vegetation Secondary vegetation 
# 1779                 1317                 1516                 1483

##%######################################################%##
#                                                          #
####         remove potentially dodgy sites and         ####
####            rerun climate and LU models             ####
#                                                          #
##%######################################################%##

# read in a later dataset to get climate info

preddir <-"6_RunLUClimateModels"

sites2 <- readRDS(paste0(preddir, "/PREDICTSSiteData.rds"))

sites2 <- sites2[!sites2$SSBS %in% query$SSBS, ] # 6000 rows


# 1. Abundance, mean anomaly
MeanAnomalyModelAbund <- GLMERSelect(modelData = sites2,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyModelAbund, file = paste0(outDir, "/MeanAnomalyModelAbund.rdata"))

# plot the results
exclQuantiles <- c(0.025,0.975)


pdf(file = paste0(outDir, "/Figure2_MeanAnom_Abun_removedsites.pdf"), width = 4, height = 4)

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
  originalX = sites2$StdTmeanAnomaly)

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
     ylim=c(min(nd$PredLower,na.rm = TRUE),max(nd$PredUpper,na.rm = TRUE)),
     xlab="Standardised Temperature Anomaly",ylab="Change in Abundance (%)", cex.lab = 0.8, cex.axis = 0.8)

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
abline(h=60,lty=1,col="#00000022")
abline(h=40,lty=1,col="#00000022")
abline(h=20,lty=1,col="#00000022")
abline(h=0,lty=1,col="#00000022")
abline(h=-20,lty=1,col="#00000022")
abline(h=-40,lty=1,col="#00000022")
abline(h=-60,lty=1,col="#00000022")
abline(v=0,lty=1,col="#00000022")
abline(v=0.5,lty=1,col="#00000022")
abline(v=1,lty=1,col="#00000022")
abline(v=1.5,lty=1,col="#00000022")
abline(v=2,lty=1,col="#00000022")

# add legend
legend(
  x = 0, y = 70,bty="n",
  legend = c("Primary","Secondary",
             "Agriculture_Low",
             "Agriculture_High"),
  col = c("#009E73", "#0072B2",
          "#E69F00", "#D55E00"),
  lty=1,lwd=2, cex = 0.8)


dev.off()


