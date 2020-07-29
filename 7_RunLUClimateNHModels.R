##%######################################################%##
#                                                          #
####         Run models including percentage NH         ####
#                                                          #
##%######################################################%##

# This script runs the more complex models looking at the buffering effect
# of natural habitat on climate effects across land uses.

# load libraries
library(StatisticalModels)
library(predictsFunctions)
source("Functions.R")

# directories
predictsDataDir <- "6_RunLUClimateModels/"
outDir <- "7_RunLUClimateNHModels/"


###Create Models for all insects in predicts for standardised climate anomaly and Land interactions

predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSiteData.rds"))

# remove any NAs
modelData <- na.omit(predictsSites[,c(
  'LogAbund','UI2','StdTmeanAnomalyRS','SS','SSB','SSBS','NH_5000.rs')])

# run models with and without NH interaction, abundance models
AbundMeanAnomalyModel0 <- GLMER(modelData = modelData,responseVar = "LogAbund",fitFamily = "gaussian",
                fixedStruct = "UI2 * StdTmeanAnomalyRS",
                randomStruct = "(1|SS)+(1|SSB)")
AbundMeanAnomalyModel1 <- GLMER(modelData = modelData,responseVar = "LogAbund",fitFamily = "gaussian",
                fixedStruct = "UI2 * StdTmeanAnomalyRS * NH_5000.rs",
                randomStruct = "(1|SS)+(1|SSB)")

print(anova(AbundMeanAnomalyModel0$model,AbundMeanAnomalyModel1$model))

# run models with and without NH interaction, species richness models
modelData <- na.omit(predictsSites[,c(
  'Species_richness','UI2','StdTmeanAnomalyRS','SS','SSB','SSBS','NH_5000.rs')])

RichMeanAnomalyModel0 <- GLMER(modelData = modelData,responseVar = "Species_richness",fitFamily = "poisson",
                               fixedStruct = "UI2 * StdTmeanAnomalyRS",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)")
RichMeanAnomalyModel1 <- GLMER(modelData = modelData,responseVar = "Species_richness",fitFamily = "poisson",
                               fixedStruct = "UI2 * StdTmeanAnomalyRS * NH_5000.rs",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)")

print(anova(RichMeanAnomalyModel0$model,RichMeanAnomalyModel1$model))

# run models with and without NH interaction, abundance models, max anomaly
modelData <- na.omit(predictsSites[,c(
  'LogAbund','UI2','StdTmaxAnomalyRS','SS','SSB','SSBS','NH_5000.rs')])

AbundMaxAnomalyModel0 <- GLMER(modelData = modelData,responseVar = "LogAbund",fitFamily = "gaussian",
                                fixedStruct = "UI2 * StdTmaxAnomalyRS",
                                randomStruct = "(1|SS)+(1|SSB)")
AbundMaxAnomalyModel1 <- GLMER(modelData = modelData,responseVar = "LogAbund",fitFamily = "gaussian",
                                fixedStruct = "UI2 * StdTmaxAnomalyRS * NH_5000.rs",
                                randomStruct = "(1|SS)+(1|SSB)")

print(anova(AbundMaxAnomalyModel0$model,AbundMaxAnomalyModel1$model))

# run models with and without NH interaction, species richness models, max anomaly
modelData <- na.omit(predictsSites[,c(
  'Species_richness','UI2','StdTmaxAnomalyRS','SS','SSB','SSBS','NH_5000.rs')])

RichMaxAnomalyModel0 <- GLMER(modelData = modelData,responseVar = "Species_richness",fitFamily = "poisson",
                               fixedStruct = "UI2 * StdTmaxAnomalyRS",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)")
RichMaxAnomalyModel1 <- GLMER(modelData = modelData,responseVar = "Species_richness",fitFamily = "poisson",
                               fixedStruct = "UI2 * StdTmaxAnomalyRS * NH_5000.rs",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)")

print(anova(RichMaxAnomalyModel0$model,RichMaxAnomalyModel1$model))

# For insects

# 25% NH cover = -1.015469 in rescaled values
# 50% NH cover = 0.01493712 in rescaled values
# 75% NH cover = 1.045653 in rescaled values
# 100% NH cover = 2.073849 in rescaled values

# For disease hosts

# 25% NH cover = -1.122627 in rescaled values
# 50% NH cover = -0.1351849 in rescaled values
# 75% NH cover = 0.8498366 in rescaled values
# 100% NH cover = 1.834235 in rescaled values

# create matrix for predictions
nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(AbundMeanAnomalyModel1$data$UI2)),
  NH_5000.rs=c(-1.015469,0.01493712,1.045653,2.073849))
  # NH_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform NH data range
nd$NH_5000 <- round(BackTransformCentreredPredictor(
  transformedX = nd$NH_5000.rs,originalX = predictsSites$NH_5000)*100,0)

# set values for richness and abundance
nd$LogAbund <- 0
nd$Species_richness <- 0

# set the reference row
refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))) & 
                  (nd$NH_5000==100))

# quantiles for presenting results
exclQuantiles <- c(0.025,0.975)


QPV <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  AbundMeanAnomalyModel1$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict results 
a.preds.tmean <- PredictGLMERRandIter(model = AbundMeanAnomalyModel1$model,data = nd)

# transform results
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to percentage of reference row
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS > QAH[2])),] <- NA

# get the median and upper/lower intervals for plots
nd$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

pdf(file = paste0(outDir,"AbundanceMeanAnomaly.pdf"),width = 20/2.54,height = 10/2.54)

par(mfrow=c(1,2))


nd_ab_mean <- nd

## 1. plot for abundance response to mean anomaly wuth 

ylims <- with(nd[nd$UI2 %in% c("Agriculture_Low","Agriculture_High"),],
              c(min(PredLower,na.rm = TRUE),max(PredUpper,na.rm = TRUE)))
xlims <- with(nd[nd$UI2 %in% c("Agriculture_Low","Agriculture_High"),],
              c(min(StdTmeanAnomaly,na.rm = TRUE),max(StdTmeanAnomaly,na.rm = TRUE)))

invisible(lapply(X = split(x = nd,f = nd$UI2)[c(3,2)],FUN = function(preds.lu){
  
  plot(-9e99,-9e99,xlim=xlims,ylim=ylims,
       xlab="Mean temperature anomaly",ylab="Abundance (%)", cex.lab = 0.8, cex.axis = 0.8)
  
  invisible(mapply(FUN = function(preds,col){
    
    preds <- na.omit(preds)
    
    X.Vec <- c(preds$StdTmeanAnomaly, max(preds$StdTmeanAnomaly), 
               rev(preds$StdTmeanAnomaly), min(preds$StdTmeanAnomaly))
    Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
               rev(preds$PredUpper), (preds$PredLower)[1])
    
    polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
    
    points(x = preds$StdTmeanAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
    
    
  },split(x = preds.lu,f = preds.lu$NH_5000),c("#a50026","#f46d43","#74add1","#313695")))
  
  
  abline(h=150,lty=1,col="#00000022")
  abline(h=100,lty=1,col="#00000022")
  abline(h=50,lty=1,col="#00000022")
  abline(h=0,lty=1,col="#00000022")
  abline(h=-50,lty=1,col="#00000022")
  abline(v=0,lty=1,col="#00000022")
  abline(v=1,lty=1,col="#00000022")
  abline(v=2,lty=1,col="#00000022")
  
  

}))

# add legend to righthand plot
legend(
  x = 1.5,y = 112, bty="n",
  legend = c("25%", "50%", "75%", "100%"),
  col = c("#a50026","#f46d43","#74add1","#313695"),
  lty=1,lwd=2, cex = 0.7, title = "% NH", title.adj = 0.2)


invisible(dev.off())



pdf(file = paste0(outDir,"RichnessMeanAnomaly.pdf"),width = 17.5/2.54,height = 8/2.54)

QPV <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
  RichMeanAnomalyModel1$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

s.preds.tmean <- PredictGLMERRandIter(model = RichMeanAnomalyModel1$model,data = nd)
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

par(mfrow=c(1,2))

ylims <- with(nd[nd$UI2 %in% c("Agriculture_Low","Agriculture_High"),],
              c(min(PredLower,na.rm = TRUE),max(PredUpper,na.rm = TRUE)))
xlims <- with(nd[nd$UI2 %in% c("Agriculture_Low","Agriculture_High"),],
              c(min(StdTmeanAnomaly,na.rm = TRUE),max(StdTmeanAnomaly,na.rm = TRUE)))

invisible(lapply(X = split(x = nd,f = nd$UI2)[c(3,2)],FUN = function(preds.lu){
  
  plot(-9e99,-9e99,xlim=xlims,ylim=ylims,
       xlab="Mean temperature anomaly",ylab="Richness (%)")
  
  invisible(mapply(FUN = function(preds,col){
    
    preds <- na.omit(preds)
    
    X.Vec <- c(preds$StdTmeanAnomaly, max(preds$StdTmeanAnomaly), 
               rev(preds$StdTmeanAnomaly), min(preds$StdTmeanAnomaly))
    Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
               rev(preds$PredUpper), (preds$PredLower)[1])
    
    polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
    
    points(x = preds$StdTmeanAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
    
  },split(x = preds.lu,f = preds.lu$NH_5000),c("#a50026","#f46d43","#74add1","#313695")))
  
  
  abline(h=150,lty=1,col="#00000022")
  abline(h=100,lty=1,col="#00000022")
  abline(h=50,lty=1,col="#00000022")
  abline(h=0,lty=1,col="#00000022")
  abline(h=-50,lty=1,col="#00000022")
  abline(v=0,lty=1,col="#00000022")
  abline(v=1,lty=1,col="#00000022")
  abline(v=2,lty=1,col="#00000022")
  
  
}))

invisible(dev.off())


nd <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS),
                        to = max(AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(AbundMaxAnomalyModel1$data$UI2)),
  NH_5000.rs=c(-1.015469,0.01493712,1.045653,2.073849))
  # NH_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))
nd$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)
nd$NH_5000 <- round(BackTransformCentreredPredictor(
  transformedX = nd$NH_5000.rs,originalX = predictsSites$NH_5000)*100,0)
nd$LogAbund <- 0
nd$Species_richness <- 0

refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmaxAnomaly==min(abs(nd$StdTmaxAnomaly))) & 
                  (nd$NH_5000==100))

exclQuantiles <- c(0.025,0.975)

pdf(file = paste0(outDir,"AbundanceMaxAnomaly.pdf"),width = 17.5/2.54,height = 8/2.54)

QPV <- quantile(x = AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS[
  AbundMaxAnomalyModel1$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS[
  AbundMaxAnomalyModel1$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS[
  AbundMaxAnomalyModel1$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS[
  AbundMaxAnomalyModel1$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

a.preds.tmax <- PredictGLMERRandIter(model = AbundMaxAnomalyModel1$model,data = nd)
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

par(mfrow=c(1,2))

ylims <- with(nd[nd$UI2 %in% c("Agriculture_Low","Agriculture_High"),],
              c(min(PredLower,na.rm = TRUE),max(PredUpper,na.rm = TRUE)))
xlims <- with(nd[nd$UI2 %in% c("Agriculture_Low","Agriculture_High"),],
              c(min(StdTmaxAnomaly,na.rm = TRUE),max(StdTmaxAnomaly,na.rm = TRUE)))

invisible(lapply(X = split(x = nd,f = nd$UI2)[c(3,2)],FUN = function(preds.lu){
  
  plot(-9e99,-9e99,xlim=xlims,ylim=ylims,
       xlab="Maximum temperature anomaly",ylab="Abundance (%)")
  
  invisible(mapply(FUN = function(preds,col){
    
    preds <- na.omit(preds)
    
    X.Vec <- c(preds$StdTmaxAnomaly, max(preds$StdTmaxAnomaly), 
               rev(preds$StdTmaxAnomaly), min(preds$StdTmaxAnomaly))
    Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
               rev(preds$PredUpper), (preds$PredLower)[1])
    
    polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
    
    points(x = preds$StdTmaxAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
    
  },split(x = preds.lu,f = preds.lu$NH_5000),c("#a50026","#f46d43","#74add1","#313695")))
  
  
  abline(h=150,lty=1,col="#00000022")
  abline(h=100,lty=1,col="#00000022")
  abline(h=50,lty=1,col="#00000022")
  abline(h=0,lty=1,col="#00000022")
  abline(h=-50,lty=1,col="#00000022")
  abline(v=0,lty=1,col="#00000022")
  abline(v=1,lty=1,col="#00000022")
  abline(v=2,lty=1,col="#00000022")
  
  
}))

invisible(dev.off())



pdf(file = paste0(outDir,"RichnessMaxAnomaly.pdf"),width = 17.5/2.54,height = 8/2.54)

QPV <- quantile(x = RichMaxAnomalyModel1$data$StdTmaxAnomalyRS[
  RichMaxAnomalyModel1$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = RichMaxAnomalyModel1$data$StdTmaxAnomalyRS[
  RichMaxAnomalyModel1$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = RichMaxAnomalyModel1$data$StdTmaxAnomalyRS[
  RichMaxAnomalyModel1$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = RichMaxAnomalyModel1$data$StdTmaxAnomalyRS[
  RichMaxAnomalyModel1$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

s.preds.tmax <- PredictGLMERRandIter(model = RichMaxAnomalyModel1$model,data = nd)
s.preds.tmax <- exp(s.preds.tmax)-0.01

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

par(mfrow=c(1,2))

ylims <- with(nd[nd$UI2 %in% c("Agriculture_Low","Agriculture_High"),],
              c(min(PredLower,na.rm = TRUE),max(PredUpper,na.rm = TRUE)))
xlims <- with(nd[nd$UI2 %in% c("Agriculture_Low","Agriculture_High"),],
              c(min(StdTmaxAnomaly,na.rm = TRUE),max(StdTmaxAnomaly,na.rm = TRUE)))

invisible(lapply(X = split(x = nd,f = nd$UI2)[c(3,2)],FUN = function(preds.lu){
  
  plot(-9e99,-9e99,xlim=xlims,ylim=ylims,
       xlab="Maximum temperature anomaly",ylab="Richness (%)")
  
  invisible(mapply(FUN = function(preds,col){
    
    preds <- na.omit(preds)
    
    X.Vec <- c(preds$StdTmaxAnomaly, max(preds$StdTmaxAnomaly), 
               rev(preds$StdTmaxAnomaly), min(preds$StdTmaxAnomaly))
    Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
               rev(preds$PredUpper), (preds$PredLower)[1])
    
    polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
    
    points(x = preds$StdTmaxAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
    
  },split(x = preds.lu,f = preds.lu$NH_5000),c("#a50026","#f46d43","#74add1","#313695")))
  
  
  abline(h=150,lty=1,col="#00000022")
  abline(h=100,lty=1,col="#00000022")
  abline(h=50,lty=1,col="#00000022")
  abline(h=0,lty=1,col="#00000022")
  abline(h=-50,lty=1,col="#00000022")
  abline(v=0,lty=1,col="#00000022")
  abline(v=1,lty=1,col="#00000022")
  abline(v=2,lty=1,col="#00000022")
  
  
}))

invisible(dev.off())




###### figures for manuscript ######

# abundance response to mean anomaly only, all LUs


library(ggplot2)

# set factor levels
nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))

nd$NH_5000 <- factor(nd$NH_5000, levels = c("100", "75", "50", "25"))

# plot
ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = NH_5000), size = 1) +
  geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = NH_5000), alpha = 0.2) +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 2) + 
  theme_bw() + 
  labs(fill = "% NH", col = "% NH") + 
  ylab("Abundance (%)") +
  xlab("Mean temperature anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 200)) + 
  theme(aspect.ratio = 1, text = element_text(size = 15))

# save
ggsave(filename = paste0(outDir, "Figure_3_ab_mean.pdf"), height = 8, width = 8)
