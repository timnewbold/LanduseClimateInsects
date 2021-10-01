##%######################################################%##
#                                                          #
####           Testing abundance models when            ####
####             abundance is not rescaled              ####
#                                                          #
##%######################################################%##

# in this script I will rerun the models in script 6 but without 
# rescaling abundance as the response variable to see if this improves
# model checks. 

# load libraries
library(devtools)
#install_github("timnewbold/StatisticalModels")
library(StatisticalModels)
#install_github(repo = "timnewbold/predicts-demo",subdir = "predictsFunctions", force = T)
library(predictsFunctions)
source("Functions.R")
library(sjPlot)
library(cowplot)
library(ggplot2)
library(roquefort)

# directories 
predictsDataDir <- "6_RunLUClimateModels/"
outDir <- "14_Additional_Tests/"

# load the dataset
predictsSites <- readRDS(file = paste0(predictsDataDir,"PREDICTSSiteData.rds"))


model_data <- predictsSites[!is.na(predictsSites$Total_abundance), ] # 5759 rows

model_data$LogAbund_new <- log(model_data$Total_abundance + 0.01)

# 1. Abundance, mean anomaly
MeanAnomalyModelAbund_test <- GLMERSelect(modelData = model_data,responseVar = "LogAbund_new",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyModelAbund_test, file = paste0(outDir, "MeanAnomalyModelAbund_notRSAbun.rdata"))



## compare summaries with original model

# load in the original model
load(paste0(predictsDataDir, "MeanAnomalyModelAbund.rdata"))

summary(MeanAnomalyModelAbund_test$model)

summary(MeanAnomalyModelAbund$model)




#### plot results from both models ####

# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)


pdf(file = paste0(outDir, "/Figure2_Abundance_notrescaled.pdf"), width = 8, height = 4)

par(mfrow=c(1,2))
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
  originalX = predictsSites$StdTmeanAnomaly)

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
     ylim=c(-70,80),
     xlab="Standardised Temperature Anomaly",ylab="Change in total abundance (%)", cex.lab = 0.8, cex.axis = 0.8)

title("a. original model", adj = 0, cex.main = 1)


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

### and the second plot for where abundance has not been rescaled

nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund_test$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelAbund_test$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelAbund_test$data$UI2)))

# back transform the predictors
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd2$LogAbund_new <- 0
nd2$Species_richness <- 0

# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))))

# adjust plot 1: mean anomaly and abundance

QPV <- quantile(x = MeanAnomalyModelAbund_test$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_test$data$UI2=="Primary vegetation"],
  probs = exclQuantiles)
QSV <- quantile(x = MeanAnomalyModelAbund_test$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_test$data$UI2=="Secondary vegetation"],
  probs = exclQuantiles)
QAL <- quantile(x = MeanAnomalyModelAbund_test$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_test$data$UI2=="Agriculture_Low"],
  probs = exclQuantiles)
QAH <- quantile(x = MeanAnomalyModelAbund_test$data$StdTmeanAnomalyRS[
  MeanAnomalyModelAbund_test$data$UI2=="Agriculture_High"],
  probs = exclQuantiles)

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund_test$model,data = nd2)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd2$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# plot #

# set up plotting window
plot(-9e99,-9e99,xlim=c(min(nd2$StdTmeanAnomaly),2),
     ylim=c(-70,200),
     xlab="Standardised Temperature Anomaly",ylab="Change in total abundance (%)", cex.lab = 0.8, cex.axis = 0.8)

title("b. Abundance just log(x+0.01) transformed", adj = 0, cex.main = 1)


invisible(mapply(FUN = function(preds,col){
  
  preds <- na.omit(preds)
  
  X.Vec <- c(preds$StdTmeanAnomaly, max(preds$StdTmeanAnomaly), 
             rev(preds$StdTmeanAnomaly), min(preds$StdTmeanAnomaly))
  Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
             rev(preds$PredUpper), (preds$PredLower)[1])
  
  polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
  
  points(x = preds$StdTmeanAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
  
},split(nd2,nd2$UI2),c("#009E73", "#D55E00", "#E69F00", "#0072B2")))

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
  x = -0.1, y = 200,bty="n",
  legend = c("Primary","Secondary",
             "Agriculture_Low",
             "Agriculture_High"),
  col = c("#009E73", "#0072B2",
          "#E69F00", "#D55E00"),
  lty=1,lwd=2, cex = 0.8)


dev.off()

##### model checks for both models ####




mod_list <- c("MeanAnomalyModelAbund", "MeanAnomalyModelAbund_test")

#x <- mod_list[1]

for(x in mod_list){
  
  # only do this one if not SR model
  if(grepl("Abun", x) ==1) {
    
    
    ## 1. Checking the fitted vs residuals relationship
    p1 <- plot(get(x)$model)
  }
  
  ## 2. Normality of Residuals
  pdf(NULL)
  dev.control(displaylist="enable")
  qqnorm(resid(get(x)$model), main = "")
  qqline(resid(get(x)$model))
  p2 <- recordPlot()
  invisible(dev.off())
  
  
  
  ## 3. Check for spatial autocorrelation
  
  sa_test<-roquefort::SpatialAutocorrelationTest(model=get(x), all.data=predictsSites)
  
  #summary(sa_test)
  
  # percentage of studies that show spatial autocorrelation?
  perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100
  
  
  sa_test_vals <- as.data.frame(sa_test$P)
  sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)
  
  label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")
  
  p3 <- ggplot(data = sa_test_vals ) +
    geom_histogram(aes(x = sa_test_vals$`sa_test$P`)) +
    geom_vline(xintercept = 0.05, col = "red") +
    geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
    theme_bw() +
    ylim(c(0, 100)) +
    xlab("P-value") +
    ylab("Frequency") +
    theme(panel.grid = element_blank(), 
          aspect.ratio = 1)
  
  
  
  if(grepl("Abun", x) == 1) {
    
    cowplot::plot_grid(p1,p2,p3,
                       labels = c(paste0("A.", x), "B.", "C."))
    ggsave(file = paste0(outDir, x, "_model_checks.pdf"), height = 10, width = 10)
    
    rm(p1, p2, p3)
    rm(perc_auto)
    
    
  }else{
    cowplot::plot_grid(p2,p3,
                       labels = c("A.", "B."))#
    
    ggsave(file = paste0(outDir, x, "_model_checks.pdf"), height = 5, width = 10) 
    
    rm(p2, p3)
    rm(perc_auto)
    
  }
  
}




