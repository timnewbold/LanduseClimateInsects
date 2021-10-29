##%######################################################%##
#                                                          #
#### 14. Testing squareroot transformation of abundance ####
#                                                          #
##%######################################################%##

# Here I test whether a squareroot transformation of scaled abundance, rather
# than log transformed scaled abundance has a impact on the behaviour of the
# model residuals. 

rm(list = ls())


# load libraries
library(devtools)
#install_github("timnewbold/StatisticalModels")
library(StatisticalModels)
#install_github(repo = "timnewbold/predicts-demo",subdir = "predictsFunctions", force = T)
library(predictsFunctions)
source("Functions.R")
library(sjPlot)

# directories 
predictsDataDir <- "6_RunLUClimateModels/"
outDir <- "14_Additional_Tests/"


# load data
predictsSites <- readRDS(file = paste0(predictsDataDir,"PREDICTSSiteData.rds"))

# squareroot transform rescaled abundance rather than log transform
predictsSites$SqrtAbund <- sqrt(predictsSites$Total_abundance_RS)

summary(predictsSites$SqrtAbund)

hist(predictsSites$SqrtAbund)



# 1. Abundance, mean anomaly

model_data <- predictsSites[!is.na(predictsSites$SqrtAbund), ] # 5759
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ] # 5735


MeanAnomalyModelAbund_sqrtAb <- GLMERSelect(modelData = model_data,responseVar = "SqrtAbund",
                                               fitFamily = "gaussian",fixedFactors = "UI2",
                                               fixedTerms = list(StdTmeanAnomalyRS=1),
                                               randomStruct = "(1|SS)+(1|SSB)",
                                               fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                               saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

summary(MeanAnomalyModelAbund_sqrtAb$model)

# save the model output
save(MeanAnomalyModelAbund_sqrtAb, file = paste0(outDir, "/MeanAnomalyModelAbund_sqrtAb.rdata"))





#### Plots ####


# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)


pdf(file = paste0(outDir, "/Figure2_MeanAnom_sqrtAb.pdf"), width = 8, height = 4)

par(mfrow=c(1,1))
par(las=1)
par(mgp=c(1.6,0.3,0)) # 
par(mar=c(2.6,2.6,1,1)) # margins around plot
par(tck=-0.01) # tick mark size
#par(pty="s") # set plot type to be square

nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund_sqrtAb$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelAbund_sqrtAb$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelAbund_sqrtAb$data$UI2)))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$SqrtAbund <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyModelAbund_sqrtAb$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_sqrtAb$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelAbund_sqrtAb$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_sqrtAb$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyModelAbund_sqrtAb$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_sqrtAb$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyModelAbund_sqrtAb$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_sqrtAb$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund_sqrtAb$model,data = nd)

# back transform the abundance values
a.preds.tmean <- a.preds.tmean^2

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

title("a", adj = 0, cex.main = 1)


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



#### model check plots ####


pdf(NULL)
dev.control(displaylist="enable")
plot(fitted(MeanAnomalyModelAbund_sqrtAb$model), residuals(MeanAnomalyModelAbund_sqrtAb$model))
p1 <- recordPlot()
invisible(dev.off())


pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(MeanAnomalyModelAbund_sqrtAb$model), main = "")
qqline(resid(MeanAnomalyModelAbund_sqrtAb$model))
p2 <- recordPlot()
invisible(dev.off())

pdf(NULL)
dev.control(displaylist="enable")
plot(model_data$SqrtAbund,fitted(MeanAnomalyModelAbund_sqrtAb$model), 
     xlab = "Observed values", ylab = "Fitted values") 
abline(a = 0, b = 1, col = "red", lwd = 2)
p3 <- recordPlot()
invisible(dev.off())



cowplot::plot_grid(p1,p2,p3,
                   labels = c("a.", "b.", "c."))#

ggsave(file = paste0(outDir, "MeanAnomalyModelAbund_sqrtAb_modelchecks.pdf"), height = 10, width = 10) 

