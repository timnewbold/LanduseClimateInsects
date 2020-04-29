
library(StatisticalModels)
library(predictsFunctions)
source("Functions.R")

###Create Models for all insects in predicts for standardised climate anomaly and Land interactions

predictsDataDir <- "5_PREDICTSMatchPropNatHab/"

outDir <- "6_RunLUClimateModels/"

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

predictsSites$NH_1000.rs <- StdCenterPredictor(predictsSites$NH_1000)
predictsSites$NH_3000.rs <- StdCenterPredictor(predictsSites$NH_3000)
predictsSites$NH_5000.rs <- StdCenterPredictor(predictsSites$NH_5000)
predictsSites$NH_10000.rs <- StdCenterPredictor(predictsSites$NH_10000)

predictsSites <- RescaleAbundance(predictsSites)

saveRDS(object = predictsSites,file = paste0(outDir,"PREDICTSSiteData.rds"))

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
MaxAnomalyModelAbund <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmaxAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))
MaxAnomalyModelRich <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
                                    fitFamily = "poisson",fixedFactors = "UI2",
                                    fixedTerms = list(StdTmaxAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)"),
                                    saveVars = c("Total_abundance", "SSBS", "NH_3000"))


exclQuantiles <- c(0.025,0.975)

pdf(file = paste0(outDir,"LUClimateAnomalyInteractions.pdf"),width = 17.5/2.54,height = 16/2.54)

par(mfrow=c(2,2))

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

if(!is.null(MeanAnomalyModelAbund$model)){
  
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
       xlab="Mean temperature anomaly",ylab="Abundance (%)")
  
  invisible(mapply(FUN = function(preds,col){
    
    preds <- na.omit(preds)
    
    X.Vec <- c(preds$StdTmeanAnomaly, max(preds$StdTmeanAnomaly), 
               rev(preds$StdTmeanAnomaly), min(preds$StdTmeanAnomaly))
    Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
               rev(preds$PredUpper), (preds$PredLower)[1])
    
    polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
    
    points(x = preds$StdTmeanAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
    
  },split(nd,nd$UI2),c("#009E73", "#D55E00", "#E69F00", "#0072B2")))
  
  abline(h=150,lty=1,col="#00000022")
  abline(h=100,lty=1,col="#00000022")
  abline(h=50,lty=1,col="#00000022")
  abline(h=0,lty=1,col="#00000022")
  abline(h=-50,lty=1,col="#00000022")
  abline(v=0,lty=1,col="#00000022")
  abline(v=1,lty=1,col="#00000022")
  abline(v=2,lty=1,col="#00000022")
  
  
} else {
  frame()
}


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

if(!is.null(MeanAnomalyModelRich$model)){
  
  s.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelRich$model,data = nd)
  s.preds.tmean <- exp(s.preds.tmean)
  
  s.preds.tmean <- sweep(x = s.preds.tmean,MARGIN = 2,STATS = s.preds.tmean[refRow,],FUN = '/')
  
  s.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1])),] <- NA
  s.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2])),] <- NA
  s.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS < QSV[1])),] <- NA
  s.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS > QSV[2])),] <- NA
  s.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS < QAL[1])),] <- NA
  s.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS > QAL[2])),] <- NA
  s.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS < QAH[1])),] <- NA
  s.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS > QAH[2])),] <- NA
  
  nd$PredMedian <- ((apply(X = s.preds.tmean,MARGIN = 1,
                           FUN = median,na.rm=TRUE))*100)-100
  nd$PredUpper <- ((apply(X = s.preds.tmean,MARGIN = 1,
                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
  nd$PredLower <- ((apply(X = s.preds.tmean,MARGIN = 1,
                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
  
  plot(-9e99,-9e99,xlim=c(min(nd$StdTmeanAnomaly),max(nd$StdTmeanAnomaly)),
       ylim=c(min(nd$PredLower,na.rm = TRUE),max(nd$PredUpper,na.rm = TRUE)),
       xlab="Mean temperature anomaly",ylab="Richness (%)")
  
  invisible(mapply(FUN = function(preds,col){
    
    preds <- na.omit(preds)
    
    X.Vec <- c(preds$StdTmeanAnomaly, max(preds$StdTmeanAnomaly), 
               rev(preds$StdTmeanAnomaly), min(preds$StdTmeanAnomaly))
    Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
               rev(preds$PredUpper), (preds$PredLower)[1])
    
    polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
    
    points(x = preds$StdTmeanAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
    
  },split(nd,nd$UI2),c("#009E73", "#D55E00", "#E69F00", "#0072B2")))
  
  
  abline(h=150,lty=1,col="#00000022")
  abline(h=100,lty=1,col="#00000022")
  abline(h=50,lty=1,col="#00000022")
  abline(h=0,lty=1,col="#00000022")
  abline(h=-50,lty=1,col="#00000022")
  abline(v=0,lty=1,col="#00000022")
  abline(v=1,lty=1,col="#00000022")
  abline(v=2,lty=1,col="#00000022")
  
  
}

frame()


nd <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
                        to = max(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxAnomalyModelAbund$data$UI2)))
nd$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)
nd$LogAbund <- 0
nd$Species_richness <- 0

refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmaxAnomaly==min(abs(nd$StdTmaxAnomaly))))

QPV <- quantile(x = MaxAnomalyModelAbund$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxAnomalyModelAbund$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxAnomalyModelAbund$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxAnomalyModelAbund$data$StdTmaxAnomalyRS[
  MaxAnomalyModelAbund$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

if(!is.null(MaxAnomalyModelAbund$model)){
  
  a.preds.tmax <- PredictGLMERRandIter(model = MaxAnomalyModelAbund$model,data = nd)
  a.preds.tmax <- exp(a.preds.tmax)-0.01
  
  a.preds.tmax <- sweep(x = a.preds.tmax,MARGIN = 2,STATS = a.preds.tmax[refRow,],FUN = '/')
  
  a.preds.tmax[which((nd$UI2=="Primary vegetation") & (nd$StdTmaxAnomalyRS < QPV[1])),] <- NA
  a.preds.tmax[which((nd$UI2=="Primary vegetation") & (nd$StdTmaxAnomalyRS > QPV[2])),] <- NA
  a.preds.tmax[which((nd$UI2=="Secondary vegetation") & (nd$StdTmaxAnomalyRS < QSV[1])),] <- NA
  a.preds.tmax[which((nd$UI2=="Secondary vegetation") & (nd$StdTmaxAnomalyRS > QSV[2])),] <- NA
  a.preds.tmax[which((nd$UI2=="Agriculture_Low") & (nd$StdTmaxAnomalyRS < QAL[1])),] <- NA
  a.preds.tmax[which((nd$UI2=="Agriculture_Low") & (nd$StdTmaxAnomalyRS > QAL[2])),] <- NA
  a.preds.tmax[which((nd$UI2=="Agriculture_High") & (nd$StdTmaxAnomalyRS < QAH[1])),] <- NA
  a.preds.tmax[which((nd$UI2=="Agriculture_High") & (nd$StdTmaxAnomalyRS > QAH[2])),] <- NA
  
  nd$PredMedian <- ((apply(X = a.preds.tmax,MARGIN = 1,
                           FUN = median,na.rm=TRUE))*100)-100
  nd$PredUpper <- ((apply(X = a.preds.tmax,MARGIN = 1,
                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
  nd$PredLower <- ((apply(X = a.preds.tmax,MARGIN = 1,
                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
  
  plot(-9e99,-9e99,xlim=c(min(nd$StdTmaxAnomaly),max(nd$StdTmaxAnomaly)),
       ylim=c(min(nd$PredLower,na.rm = TRUE),max(nd$PredUpper,na.rm = TRUE)),
       xlab="Maximum temperature anomaly",ylab="Abundance (%)")
  
  invisible(mapply(FUN = function(preds,col){
    
    preds <- na.omit(preds)
    
    X.Vec <- c(preds$StdTmaxAnomaly, max(preds$StdTmaxAnomaly), 
               rev(preds$StdTmaxAnomaly), min(preds$StdTmaxAnomaly))
    Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
               rev(preds$PredUpper), (preds$PredLower)[1])
    
    polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
    
    points(x = preds$StdTmaxAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
    
  },split(nd,nd$UI2),c("#009E73", "#D55E00", "#E69F00", "#0072B2")))
  
  abline(h=150,lty=1,col="#00000022")
  abline(h=100,lty=1,col="#00000022")
  abline(h=50,lty=1,col="#00000022")
  abline(h=0,lty=1,col="#00000022")
  abline(h=-50,lty=1,col="#00000022")
  abline(v=0,lty=1,col="#00000022")
  abline(v=1,lty=1,col="#00000022")
  abline(v=2,lty=1,col="#00000022")
  
} else {
  frame()
}



QPV <- quantile(x = MaxAnomalyModelRich$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MaxAnomalyModelRich$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MaxAnomalyModelRich$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MaxAnomalyModelRich$data$StdTmaxAnomalyRS[
  MaxAnomalyModelRich$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

if(!is.null(MaxAnomalyModelRich$model)){
  
  s.preds.tmax <- PredictGLMERRandIter(model = MaxAnomalyModelRich$model,data = nd)
  s.preds.tmax <- exp(s.preds.tmax)
  
  s.preds.tmax <- sweep(x = s.preds.tmax,MARGIN = 2,STATS = s.preds.tmax[refRow,],FUN = '/')
  
  s.preds.tmax[which((nd$UI2=="Primary vegetation") & (nd$StdTmaxAnomalyRS < QPV[1])),] <- NA
  s.preds.tmax[which((nd$UI2=="Primary vegetation") & (nd$StdTmaxAnomalyRS > QPV[2])),] <- NA
  s.preds.tmax[which((nd$UI2=="Secondary vegetation") & (nd$StdTmaxAnomalyRS < QSV[1])),] <- NA
  s.preds.tmax[which((nd$UI2=="Secondary vegetation") & (nd$StdTmaxAnomalyRS > QSV[2])),] <- NA
  s.preds.tmax[which((nd$UI2=="Agriculture_Low") & (nd$StdTmaxAnomalyRS < QAL[1])),] <- NA
  s.preds.tmax[which((nd$UI2=="Agriculture_Low") & (nd$StdTmaxAnomalyRS > QAL[2])),] <- NA
  s.preds.tmax[which((nd$UI2=="Agriculture_High") & (nd$StdTmaxAnomalyRS < QAH[1])),] <- NA
  s.preds.tmax[which((nd$UI2=="Agriculture_High") & (nd$StdTmaxAnomalyRS > QAH[2])),] <- NA
  
  nd$PredMedian <- ((apply(X = s.preds.tmax,MARGIN = 1,
                           FUN = median,na.rm=TRUE))*100)-100
  nd$PredUpper <- ((apply(X = s.preds.tmax,MARGIN = 1,
                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
  nd$PredLower <- ((apply(X = s.preds.tmax,MARGIN = 1,
                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
  
  plot(-9e99,-9e99,xlim=c(min(nd$StdTmaxAnomaly),max(nd$StdTmaxAnomaly)),
       ylim=c(min(nd$PredLower,na.rm = TRUE),max(nd$PredUpper,na.rm = TRUE)),
       xlab="Maximum temperature anomaly",ylab="Richness (%)")
  
  invisible(mapply(FUN = function(preds,col){
    
    preds <- na.omit(preds)
    
    X.Vec <- c(preds$StdTmaxAnomaly, max(preds$StdTmaxAnomaly), 
               rev(preds$StdTmaxAnomaly), min(preds$StdTmaxAnomaly))
    Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
               rev(preds$PredUpper), (preds$PredLower)[1])
    
    polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
    
    points(x = preds$StdTmaxAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
    
  },split(nd,nd$UI2),c("#009E73", "#D55E00", "#E69F00", "#0072B2")))
  
  abline(h=150,lty=1,col="#00000022")
  abline(h=100,lty=1,col="#00000022")
  abline(h=50,lty=1,col="#00000022")
  abline(h=0,lty=1,col="#00000022")
  abline(h=-50,lty=1,col="#00000022")
  abline(v=0,lty=1,col="#00000022")
  abline(v=1,lty=1,col="#00000022")
  abline(v=2,lty=1,col="#00000022")
  
} else {
  frame()
}

invisible(dev.off())
