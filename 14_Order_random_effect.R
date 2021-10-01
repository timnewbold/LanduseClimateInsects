##%######################################################%##
#                                                          #
####       Additional tests, taxon random effect        ####
#                                                          #
##%######################################################%##

rm(list = ls())

# load required libraries
library(predictsFunctions)
library(ggplot2)
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




# take a look at how many sites have more than one order recorded.
test <- predicts[, c("SSBS", "Order")]

head(test)

result <- NULL

for(i in unique(test$SSBS)){
  
  n <- length(unique(test[test$SSBS == i, 2]))
  
  out <- c(as.character(i), n)
  
  result <- rbind(result, out)
}

result <- as.data.frame(result)
result[, 2] <- as.numeric(as.character(result[, 2] ))

# how many record more than one Order?
nrow(result[result$V2 > 1, ]) # 1365 sites

table(result$V2)

#    1    2    3    4    5    6    7    8    9   10   11   12   13   15   16   19 
# 6202  330  257   67   34    8   21   43   65  154  152   75   10   64   59   26 

# create categories based on orders, and then a separate group for those that 
# look at more than one order


# create a new order column
predicts$Order_new <- predicts$Order
predicts$Order_new <- as.character(predicts$Order_new)


# for those sites in the table above that record more than one taxon, create 
# a new group
predicts[predicts$SSBS %in% result[result$V2 > 1, "V1"], "Order_new"] <- "Multiple"


table(predicts$Order_new)
#         Coleoptera     Diptera   Hemiptera Hymenoptera    Isoptera Lepidoptera 
# 30      315636        4115       22477      139581         486      120217 
# Multiple  Orthoptera 
# 222020        1730 

# remove the one study where no orders where known
predicts <- predicts[!predicts$Order_new == "", ] # 826262 rows

# set as a factor level
predicts$Order_new <- as.factor(predicts$Order_new)

predicts <- droplevels(predicts)

# Calculate site metrics of diversity
sites <- SiteMetrics(diversity = predicts,
                     extra.cols = c("Predominant_land_use",
                                    "SSB","SSBS", "Biome", "Order_new"))

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


# save the prepared dataset
saveRDS(object = sites,file = paste0(outDir,"PREDICTSSiteData_newOrder.rds"))


# load in version with the climate data, merge
predictsSites <- readRDS("5_PREDICTSMatchPropNatHab/PREDICTSSitesWithClimateAndNatHab.rds")
predictsSites <- predictsSites@data

predictsSites <- predictsSites[ , c("SSBS", "StdTmeanAnomaly")]

predictsSites <- merge(sites, predictsSites, by = "SSBS")


predictsSites$UI2 <- factor(predictsSites$UI2)
predictsSites$UI2 <- relevel(predictsSites$UI2,ref="Primary vegetation")


#rescale the climate data
predictsSites$StdTmeanAnomalyRS <- StdCenterPredictor(predictsSites$StdTmeanAnomaly)


# rescaling abundance and log values
predictsSites <- RescaleAbundance(predictsSites)

# charlie added this line as later bits were throwing errors
predictsSites <- droplevels(predictsSites)


# some of the climate values are NA since they do not meet the thresholds
predictsSites <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ]

saveRDS(object = predictsSites,file = paste0(outDir,"PREDICTSSiteData_newOrder.rds"))
#predictsSites <- readRDS(file = paste0(outDir,"PREDICTSSiteData_newOrder.rds"))

#### try the model with the order classifications as a random effect

# 1. Abundance, mean anomaly
am1 <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)+(1|Order_new)",
                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS"))

# save the model output
save(am1, file = paste0(outDir, "/MeanAnomalyModelAbund_Order.rdata"))
#load(paste0(outDir, "/MeanAnomalyModelAbund_Order.rdata"))

summary(am1$model)
R2GLMER(am1$model) # check the R2 values 


# 2. Abundance, mean anomaly, original
am2 <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                   fitFamily = "gaussian",fixedFactors = "UI2",
                   fixedTerms = list(StdTmeanAnomalyRS=1),
                   randomStruct = "(1|SS)+(1|SSB)",
                   fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                   saveVars = c("Species_richness", "Total_abundance", "SSBS"))

# save the model output
save(am2, file = paste0(outDir, "/MeanAnomalyModelAbund_NoOrder.rdata"))
#load(paste0(outDir, "/MeanAnomalyModelAbund_NoOrder.rdata"))

summary(am2$model)
R2GLMER(am2$model) # check the R2 values 

AIC(am1$model, am2$model)

am1$stats
am2$stats

plot(am1$model)
plot(am2$model)


#### take a look at the plots ####



# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)


pdf(file = paste0(outDir, "/Figure2_MeanAnom_Abun_Rich_incOrder.pdf"), width = 8, height = 4)

par(mfrow=c(1,2))
par(las=1)
par(mgp=c(1.6,0.3,0)) # 
par(mar=c(2.6,2.6,1,1)) # margins around plot
par(tck=-0.01) # tick mark size
#par(pty="s") # set plot type to be square

nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(am1$data$StdTmeanAnomalyRS),
                        to = max(am1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(am1$data$UI2)))

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

QPV <- quantile(x = am1$data$StdTmeanAnomalyRS[
  am1$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = am1$data$StdTmeanAnomalyRS[
  am1$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = am1$data$StdTmeanAnomalyRS[
  am1$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = am1$data$StdTmeanAnomalyRS[
  am1$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = am1$model,data = nd)

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

title("inc Order random effect", adj = 0, cex.main = 1)


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


nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(am2$data$StdTmeanAnomalyRS),
                        to = max(am2$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(am2$data$UI2)))

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

QPV <- quantile(x = am2$data$StdTmeanAnomalyRS[
  am2$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = am2$data$StdTmeanAnomalyRS[
  am2$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = am2$data$StdTmeanAnomalyRS[
  am2$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = am2$data$StdTmeanAnomalyRS[
  am2$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = am2$model,data = nd)

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

title("Original model", adj = 0, cex.main = 1)


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


dev.off()
