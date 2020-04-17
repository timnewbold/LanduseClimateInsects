
library(StatisticalModels)
library(predictsFunctions)
source("Functions.R")

###Create Models for all insects in predicts for standardised climate anomaly and Land interactions

predictsDataDir <- "5_PREDICTSMatchPropNatHab/"

outDir <- ""

predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSitesWithClimateAndNatHab.rds"))
predictsSites <- predictsSites@data

predictsSites$UI2[(predictsSites$UI2=="Urban")] <- NA
predictsSites$UI2 <- factor(predictsSites$UI2)
predictsSites$UI2 <- relevel(predictsSites$UI2,ref="Primary vegetation")

predictsSites$TmeanAnomaly <- predictsSites$climate_anomaly
predictsSites$StdTmeanAnomaly <- predictsSites$TmeanAnomaly / predictsSites$historic_sd
predictsSites$StdTmeanAnomalyRS <- StdCenterPredictor(predictsSites$StdTmeanAnomaly)

predictsSites$TmaxAnomaly <- predictsSites$tmax_anomaly
predictsSites$StdTmaxAnomaly <- predictsSites$TmaxAnomaly / predictsSites$historic_sd_tmax
predictsSites$StdTmaxAnomalyRS <- StdCenterPredictor(predictsSites$StdTmaxAnomaly)

predictsSites$NH_1000 <- predictsSites$PV_1000 + predictsSites$SV_1000
predictsSites$NH_3000 <- predictsSites$PV_3000 + predictsSites$SV_3000
predictsSites$NH_5000 <- predictsSites$PV_5000 + predictsSites$SV_5000
predictsSites$NH_10000 <- predictsSites$PV_10000 + predictsSites$SV_10000

predictsSites <- RescaleAbundance(predictsSites)

MeanAnomalyModelAbund <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                fitFamily = "gaussian",fixedFactors = "UI2",
                                fixedTerms = list(StdTmeanAnomalyRS=1),
                                randomStruct = "(1|SS)+(1|SSB)",
                                fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))
MeanAnomalyModelRich <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
                                     fitFamily = "poisson",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("Total_abundance", "SSBS", "NH_3000"))

nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelAbund$data$UI2)))
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)
nd$LogAbund <- 0
nd$Species_richness <- 0

refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))

exclQuantiles <- c(0.01,0.99)


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

a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund$model,data = nd)
a.preds.tmean <- exp(a.preds.tmean)-0.01

a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS > QAH[2])),] <- NA

nd$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                                FUN = median,na.rm=TRUE))*100)-100
nd$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                               FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                               FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

plot(-9e99,-9e99,xlim=c(min(nd$StdTmeanAnomaly),max(nd$StdTmeanAnomaly)),
     ylim=c(min(nd$PredLower,na.rm = TRUE),max(nd$PredUpper,na.rm = TRUE)),
     xlab="Climate anomaly",ylab="Abundance (%)")

invisible(mapply(FUN = function(preds,col){
  
  preds <- na.omit(preds)
  
  X.Vec <- c(preds$StdTmeanAnomaly, max(preds$StdTmeanAnomaly), 
             rev(preds$StdTmeanAnomaly), min(preds$StdTmeanAnomaly))
  Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
             rev(preds$PredUpper), (preds$PredLower)[1])
  
  polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
  
  points(x = preds$StdTmeanAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
  
},split(nd,nd$UI2),c("#009E73", "#D55E00", "#E69F00", "#0072B2")))

abline(h=0,lty=2,col="#00000044")
abline(h=-50,lty=2,col="#00000044")
abline(v=0,lty=2,col="#00000044")
abline(v=1,lty=2,col="#00000044")





climate <- GLMER(predictsSites, 
                 responseVar = "LogAbund", fitFamily = "gaussian",
                 fixedStruct = "UI2*StdTmeanAnomalyRS",
                 randomStruct = "(1|SS) + (1|SSB)", 
                 saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))
summary(climate$model)
vif.mer(climate$model)

cols <- c("#009E73", "#D55E00", "#E69F00", "#0072B2")
PlotGLMERContinuous(model = climate$model,data = climate$data,
                    effects = "StdTmeanAnomalyRS", byFactor = "UI2",
                    xlab = "Standardised centered climate anomaly",ylab = "Scaled Abundance",
                    logLink = "e", ylim = c(0,0.5), main="All insects Land*std_climate_anomaly",
                    line.cols =cols, plotRug = T,seMultiplier = 1)



tmax <- GLMER(predictsSites, 
              responseVar = "LogAbund", fitFamily = "gaussian",
              fixedStruct = "UI2*TmaxAnomalyRS",
              randomStruct = "(1|SS) + (1|SSB)", saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))
summary(tmax$model)
vif.mer(tmax$model)
PlotGLMERContinuous(model = tmax$model,data = tmax$data,
                    effects = "TmaxAnomalyRS", byFactor = "UI2",
                    xlab = "Standardised centered tmax anomaly",ylab = "Scaled Abundance",
                    logLink = "e", main ="All insects Land*std_tmax_anomaly", ylim = c(0,0.5),
                    line.cols =cols, plotRug = TRUE,seMultiplier = 1)
legend(-2, 0.2, legend = levels(climate$data$Land), col = cols, pch=19, cex=1.4)


climate_rich <- GLMER(insect_sites, 
                      responseVar = "Species_richness", fitFamily = "poisson",
                      fixedStruct = "Land*std_climate_anomaly.s",
                      randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))
summary(climate_rich$model)
vif.mer(climate_rich$model)
par(mfrow=c(1,1))
PlotGLMERContinuous(model = climate_rich$model,data = climate_rich$data,
                    effects = "std_climate_anomaly.s", byFactor = "Land",
                    xlab = "Standardised centered climate anomaly",ylab = "Species Richness",
                    logLink = "e", main ="All insects Land*std_climate_anomaly", ylim = c(0,20),
                    line.cols =cols, plotRug = TRUE,seMultiplier = 1)


legend(-2, 6, legend = levels(climate$data$Land), col = cols, pch=19, cex=1.4)

tmax_rich <- GLMER(insect_sites, 
                   responseVar = "Species_richness", fitFamily = "poisson",
                   fixedStruct = "Land*std_tmax_anomaly.s",
                   randomStruct = "(1|SS) + (1|SSB)+ (1|SSBS)", saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))
summary(tmax_rich$model)
vif.mer(tmax_rich$model)
PlotGLMERContinuous(model = tmax_rich$model,data = tmax_rich$data,
                    effects = "std_tmax_anomaly.s", byFactor = "Land",
                    xlab = "Standardised centered tmax anomaly",ylab = "Species richness",
                    logLink = "e", main ="All insects Land*std_tmax_anomaly", ylim = c(0,20),
                    line.cols =cols, plotRug = TRUE,seMultiplier = 1)
legend(-2,6,  legend = levels(climate$data$Land), col = cols, pch=19, cex=1.4)

