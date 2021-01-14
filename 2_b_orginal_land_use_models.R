##%######################################################%##
#                                                          #
#### SCA models with orginal PREDICTS land use classes  ####
#                                                          #
##%######################################################%##



###  notes;

# I think we can present a graph in the supplementary information
# that shows PV, Cropland, pasture and plantation forest. The graph will show
# steep declines in cropland and pasture, but pasture will be relatively flat. 
# If we then provide a second graph of pasture split by land use intensity 
# (minimal, light and intense) it will show that  SCA has minimal impact in
# minimal and light but a large impact on intese. 

# load libraries
library(StatisticalModels)
library(predictsFunctions)
source("Functions.R")


# directories 
predictsDataDir <- "6_RunLUClimateModels/"
outDir <- "2_RunSimpleLUIModel/"

# load dataset
predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSiteData.rds"))


predictsSites <- predictsSites[!is.na(predictsSites$StdTmeanAnomaly), ] # 6308


# Combine secondary land uses
predictsSites$LandUse2 <- predictsSites$LandUse
predictsSites$LandUse2 <- as.character(predictsSites$LandUse2)
predictsSites$LandUse2[grep("secondary", predictsSites$LandUse2)] <- "Secondary vegetation"
predictsSites$LandUse2[grep("Secondary", predictsSites$LandUse2)] <- "Secondary vegetation"

unique(predictsSites$LandUse2)

predictsSites$LandUse2 <- factor(predictsSites$LandUse2, 
                                    levels = c("Primary vegetation", "Secondary vegetation", "Plantation forest", "Pasture", "Cropland"))


#### first run model with original land use classifications ####

# 1. Abundance, mean anomaly
MeanAnomalyModelAbund_LandUse <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "LandUse2",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("LandUse2:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyModelAbund_LandUse, file = paste0(outDir, "/MeanAnomalyModelAbund_LandUse_Test.rdata"))

# take a look at the results
summary(MeanAnomalyModelAbund_LandUse$model)




#### then run model with land use and use intensity combined classification ####

# Combine secondary land uses
predictsSites$UI_test <- predictsSites$UI
predictsSites$UI_test <- as.character(predictsSites$UI_test)

predictsSites$UI_test <- sub("Mature secondary vegetation", "Secondary vegetation", predictsSites$UI_test)
predictsSites$UI_test <- sub("Intermediate secondary vegetation", "Secondary vegetation", predictsSites$UI_test)
predictsSites$UI_test <- sub("Young secondary vegetation", "Secondary vegetation", predictsSites$UI_test)

predictsSites$UI_test <- sub("\\)", "", predictsSites$UI_test)
predictsSites$UI_test <- sub("\\(", "", predictsSites$UI_test)

predictsSites$UI_test <- sub("Secondary vegetation indeterminate age", "Secondary vegetation", predictsSites$UI_test)

unique(predictsSites$UI_test)

predictsSites$UI_test <- factor(predictsSites$UI_test, 
                                 levels = c("Primary vegetation_Minimal use", "Primary vegetation_Light use","Primary vegetation_Intense use",
                                            "Secondary vegetation_Minimal use", "Secondary vegetation_Light use", "Secondary vegetation_Intense use",
                                            "Plantation forest_Minimal use", "Plantation forest_Light use", "Plantation forest_Intense use",
                                            "Pasture_Minimal use", "Pasture_Light use", "Pasture_Intense use",
                                            "Cropland_Minimal use","Cropland_Light use", "Cropland_Intense use"))



# 1. Abundance, mean anomaly
MeanAnomalyModelAbund_UI <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI_test",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI_test:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyModelAbund_UI, file = paste0(outDir, "/MeanAnomalyModelAbund_UI_Test.rdata"))

# take a look at the results
summary(MeanAnomalyModelAbund_UI$model)



##### Plot the results #####

exclQuantiles <- c(0.025,0.975)

pdf(file = paste0(outDir, "/Test_LandUse.pdf"), width = 4, height = 4)

par(mfrow=c(1,1))
par(las=1)
par(mgp=c(1.6,0.3,0)) # 
par(mar=c(2.6,2.6,1,1)) # margins around plot
par(tck=-0.01) # tick mark size
#par(pty="s") # set plot type to be square

nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund_LandUse$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelAbund_LandUse$data$StdTmeanAnomalyRS),
                        length.out = 100),
  LandUse2=factor(c("Primary vegetation",  "Secondary vegetation", "Plantation forest", "Pasture", "Cropland"),
             levels = levels(MeanAnomalyModelAbund_LandUse$data$LandUse2)))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$LandUse2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyModelAbund_LandUse$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_LandUse$data$LandUse2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelAbund_LandUse$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_LandUse$data$LandUse2=="Secondary vegetation"],
  probs = exclQuantiles)
QPL <- quantile(x = MeanAnomalyModelAbund_LandUse$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_LandUse$data$LandUse2=="Plantation forest"],
  probs = exclQuantiles)
QPA <- quantile(x = MeanAnomalyModelAbund_LandUse$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_LandUse$data$LandUse2=="Pasture"],
  probs = exclQuantiles)
QCR <- quantile(x = MeanAnomalyModelAbund_LandUse$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_LandUse$data$LandUse2=="Cropland"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund_LandUse$model,data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd$LandUse2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd$LandUse2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd$LandUse2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd$LandUse2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd$LandUse2=="Plantation forest") & (nd$StdTmeanAnomalyRS < QPL[1])),] <- NA
a.preds.tmean[which((nd$LandUse2=="Plantation forest") & (nd$StdTmeanAnomalyRS > QPL[2])),] <- NA
a.preds.tmean[which((nd$LandUse2=="Pasture") & (nd$StdTmeanAnomalyRS < QPA[1])),] <- NA
a.preds.tmean[which((nd$LandUse2=="Pasture") & (nd$StdTmeanAnomalyRS > QPA[2])),] <- NA
a.preds.tmean[which((nd$LandUse2=="Cropland") & (nd$StdTmeanAnomalyRS < QCR[1])),] <- NA
a.preds.tmean[which((nd$LandUse2=="Cropland") & (nd$StdTmeanAnomalyRS > QCR[2])),] <- NA

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
     xlab="Standardised Climate Anomaly",ylab="Abundance (%)", cex.lab = 0.8, cex.axis = 0.8)

invisible(mapply(FUN = function(preds,col){
  
  preds <- na.omit(preds)
  
  X.Vec <- c(preds$StdTmeanAnomaly, max(preds$StdTmeanAnomaly), 
             rev(preds$StdTmeanAnomaly), min(preds$StdTmeanAnomaly))
  Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
             rev(preds$PredUpper), (preds$PredLower)[1])
  
  polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
  
  points(x = preds$StdTmeanAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
  
},split(nd,nd$LandUse2),c("#009E73", "#0072B2", "#6959CD", "#E69F00", "#D55E00"))) 

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
  x = -0.6,y = 120,bty="n",
  legend = c("Primary","Secondary", "Plantation forest",
             "Pasture",
             "Cropland"),
  col = c("#009E73", "#0072B2", "#6959CD",
          "#E69F00", "#D55E00"),
  lty=1,lwd=2, cex = 0.8)


dev.off()






nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund_UI$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelAbund_UI$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI_test=factor(c("Primary vegetation_Minimal use", "Primary vegetation_Light use", "Primary vegetation_Intense use",
                 "Secondary vegetation_Minimal use", "Secondary vegetation_Light use", "Secondary vegetation_Intense use",
                 "Plantation forest_Minimal use", "Plantation forest_Light use", "Plantation forest_Intense use", 
                 "Pasture_Minimal use", "Pasture_Light use" , "Pasture_Intense use", "Cropland_Minimal use", "Cropland_Light use", "Cropland_Intense use" ),
                  levels = levels(MeanAnomalyModelAbund_UI$data$UI_test)))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$UI_test=="Primary vegetation_Minimal use") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPVM <- quantile(x = MeanAnomalyModelAbund_UI$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_UI$data$UI_test=="Primary vegetation_Minimal use"],
  probs = exclQuantiles)
QPVL <- quantile(x = MeanAnomalyModelAbund_UI$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_UI$data$UI_test=="Primary vegetation_Light use"],
  probs = exclQuantiles)
QPVI <- quantile(x = MeanAnomalyModelAbund_UI$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_UI$data$UI_test=="Primary vegetation_Intense use"],
  probs = exclQuantiles)



QSVM <- quantile(x = MeanAnomalyModelAbund_UI$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_UI$data$UI_test=="Secondary vegetation_Minimal use"],
  probs = exclQuantiles)
QSVL <- quantile(x = MeanAnomalyModelAbund_UI$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_UI$data$UI_test=="Secondary vegetation_Light use"],
  probs = exclQuantiles)
QSVI <- quantile(x = MeanAnomalyModelAbund_UI$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_UI$data$UI_test=="Secondary vegetation_Intense use"],
  probs = exclQuantiles)


QPLM <- quantile(x = MeanAnomalyModelAbund_UI$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_UI$data$UI_test=="Plantation forest_Minimal use"],
  probs = exclQuantiles)
QPLL <- quantile(x = MeanAnomalyModelAbund_UI$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_UI$data$UI_test=="Plantation forest_Light use"],
  probs = exclQuantiles)
QPLI <- quantile(x = MeanAnomalyModelAbund_UI$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_UI$data$UI_test=="Plantation forest_Intense use"],
  probs = exclQuantiles)


QPAM <- quantile(x = MeanAnomalyModelAbund_UI$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_UI$data$UI_test=="Pasture_Minimal use"],
  probs = exclQuantiles)
QPAL <- quantile(x = MeanAnomalyModelAbund_UI$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_UI$data$UI_test=="Pasture_Light use"],
  probs = exclQuantiles)
QPAI <- quantile(x = MeanAnomalyModelAbund_UI$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_UI$data$UI_test=="Pasture_Intense use"],
  probs = exclQuantiles)


QCRM <- quantile(x = MeanAnomalyModelAbund_UI$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_UI$data$UI_test=="Cropland_Minimal use"],
  probs = exclQuantiles)
QCRL <- quantile(x = MeanAnomalyModelAbund_UI$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_UI$data$UI_test=="Cropland_Light use"],
  probs = exclQuantiles)
QCRI <- quantile(x = MeanAnomalyModelAbund_UI$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_UI$data$UI_test=="Cropland_Intense use"],
  probs = exclQuantiles)


# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund_UI$model,data = nd)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd$UI_test=="Primary vegetation_Minimal use") & (nd$StdTmeanAnomalyRS < QPVM[1])),] <- NA
a.preds.tmean[which((nd$UI_test=="Primary vegetation_Minimal use") & (nd$StdTmeanAnomalyRS > QPVM[2])),] <- NA
a.preds.tmean[which((nd$UI_test=="Primary vegetation_Light use") & (nd$StdTmeanAnomalyRS < QPVL[1])),] <- NA
a.preds.tmean[which((nd$UI_test=="Primary vegetation_Light use") & (nd$StdTmeanAnomalyRS > QPVL[2])),] <- NA
a.preds.tmean[which((nd$UI_test=="Primary vegetation_Intense use") & (nd$StdTmeanAnomalyRS < QPVI[1])),] <- NA
a.preds.tmean[which((nd$UI_test=="Primary vegetation_Intense use") & (nd$StdTmeanAnomalyRS > QPVI[2])),] <- NA


a.preds.tmean[which((nd$UI_test=="Secondary vegetation_Minimal use") & (nd$StdTmeanAnomalyRS < QSVM[1])),] <- NA
a.preds.tmean[which((nd$UI_test=="Secondary vegetation_Minimal use") & (nd$StdTmeanAnomalyRS > QSVM[2])),] <- NA
a.preds.tmean[which((nd$UI_test=="Secondary vegetation_Light use") & (nd$StdTmeanAnomalyRS < QSVL[1])),] <- NA
a.preds.tmean[which((nd$UI_test=="Secondary vegetation_Light use") & (nd$StdTmeanAnomalyRS > QSVL[2])),] <- NA
a.preds.tmean[which((nd$UI_test=="Secondary vegetation_Intense use") & (nd$StdTmeanAnomalyRS < QSVI[1])),] <- NA
a.preds.tmean[which((nd$UI_test=="Secondary vegetation_Intense use") & (nd$StdTmeanAnomalyRS > QSVI[2])),] <- NA


a.preds.tmean[which((nd$UI_test=="Plantation forest_Minimal use") & (nd$StdTmeanAnomalyRS < QPLM[1])),] <- NA
a.preds.tmean[which((nd$UI_test=="Plantation forest_Minimal use") & (nd$StdTmeanAnomalyRS > QPLM[2])),] <- NA
a.preds.tmean[which((nd$UI_test=="Plantation forest_Light use") & (nd$StdTmeanAnomalyRS < QPLL[1])),] <- NA
a.preds.tmean[which((nd$UI_test=="Plantation forest_Light use") & (nd$StdTmeanAnomalyRS > QPLL[2])),] <- NA
a.preds.tmean[which((nd$UI_test=="Plantation forest_Intense use") & (nd$StdTmeanAnomalyRS < QPLI[1])),] <- NA
a.preds.tmean[which((nd$UI_test=="Plantation forest_Intense use") & (nd$StdTmeanAnomalyRS > QPLI[2])),] <- NA


a.preds.tmean[which((nd$UI_test=="Pasture_Minimal use") & (nd$StdTmeanAnomalyRS < QPAM[1])),] <- NA
a.preds.tmean[which((nd$UI_test=="Pasture_Minimal use") & (nd$StdTmeanAnomalyRS > QPAM[2])),] <- NA
a.preds.tmean[which((nd$UI_test=="Pasture_Light use") & (nd$StdTmeanAnomalyRS < QPAL[1])),] <- NA
a.preds.tmean[which((nd$UI_test=="Pasture_Light use") & (nd$StdTmeanAnomalyRS > QPAL[2])),] <- NA
a.preds.tmean[which((nd$UI_test=="Pasture_Intense use") & (nd$StdTmeanAnomalyRS < QPAI[1])),] <- NA
a.preds.tmean[which((nd$UI_test=="Pasture_Intense use") & (nd$StdTmeanAnomalyRS > QPAI[2])),] <- NA


a.preds.tmean[which((nd$UI_test=="Cropland_Minimal use") & (nd$StdTmeanAnomalyRS < QCRM[1])),] <- NA
a.preds.tmean[which((nd$UI_test=="Cropland_Minimal use") & (nd$StdTmeanAnomalyRS > QCRM[2])),] <- NA
a.preds.tmean[which((nd$UI_test=="Cropland_Light use") & (nd$StdTmeanAnomalyRS < QCRL[1])),] <- NA
a.preds.tmean[which((nd$UI_test=="Cropland_Light use") & (nd$StdTmeanAnomalyRS > QCRL[2])),] <- NA
a.preds.tmean[which((nd$UI_test=="Cropland_Intense use") & (nd$StdTmeanAnomalyRS < QCRI[1])),] <- NA
a.preds.tmean[which((nd$UI_test=="Cropland_Intense use") & (nd$StdTmeanAnomalyRS > QCRI[2])),] <- NA


# Get the median, upper and lower quants for the plot
nd$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


nd$Landuse <- sub("_.*", "", nd$UI_test)
nd$UI <- sub(".*_", "", nd$UI_test)


ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI_test, linetype = UI), size = 1) +
  geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI_test), alpha = 0.2) +
  facet_wrap(~Landuse) + 
  scale_fill_manual(values = c("#458B00", "#66CD00", "#76EE00", "#27408B", "#3A5FCD", "#436EEE", "#473C8B", "#6959CD", "#836FFF", "#CD6600", "#EE7600", "#FF8C00", "#CD9B1D", "#EEB422", "#FFC125")) +
  scale_colour_manual(values = c("#458B00", "#66CD00", "#76EE00", "#27408B", "#3A5FCD", "#436EEE", "#473C8B", "#6959CD", "#836FFF", "#CD6600", "#EE7600", "#FF8C00", "#CD9B1D", "#EEB422", "#FFC125")) +
  xlim(c(-0.5,1.5)) +
  ylim(c(-100, 200))
  
ggsave(filename = paste0(outDir, "/Test_UI.pdf"), width = 8, height = 6)
