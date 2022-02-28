##%######################################################%##
#                                                          #
####             Run the climate/LUI models             ####
#                                                          #
##%######################################################%##

# This script runs the mixed effects models looking at the effect of climate
# change and landuse/use intensity.
rm(list = ls())

# directories 
predictsDataDir <- "5_PREDICTSMatchPropNatHab/"
outDir <- "6_RunLUClimateModels/"

if(!dir.exists(outDir)) dir.create(outDir)

sink(paste0(outDir,"log.txt"))

t.start <- Sys.time()

print(t.start)



# load libraries
library(devtools)
#install_github("timnewbold/StatisticalModels")
library(StatisticalModels)
#install_github(repo = "timnewbold/predicts-demo",subdir = "predictsFunctions", force = T)
library(predictsFunctions)
source("Functions.R")
library(sjPlot)
library(cowplot)


###Create Models for all insects in predicts for standardised climate anomaly and Land interactions

# read in the predicts data
predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSitesWithClimateAndNatHab.rds"))
predictsSites <- predictsSites@data

# remove the Urban sites
predictsSites$UI2 <- factor(predictsSites$UI2)
predictsSites$UI2 <- relevel(predictsSites$UI2,ref="Primary vegetation")

# rescale the variable
predictsSites$StdTmeanAnomalyRS <- StdCenterPredictor(predictsSites$StdTmeanAnomaly)
predictsSites$StdTmaxAnomalyRS <- StdCenterPredictor(predictsSites$StdTmaxAnomaly)

# organise NH data - combine primary and secondary vegetation at each scale
predictsSites$NH_1000 <- predictsSites$PV_1000 + predictsSites$SV_1000
predictsSites$NH_3000 <- predictsSites$PV_3000 + predictsSites$SV_3000
predictsSites$NH_5000 <- predictsSites$PV_5000 + predictsSites$SV_5000
predictsSites$NH_10000 <- predictsSites$PV_10000 + predictsSites$SV_10000


# rescale the variables
predictsSites$NH_1000.rs <- StdCenterPredictor(predictsSites$NH_1000)
predictsSites$NH_3000.rs <- StdCenterPredictor(predictsSites$NH_3000)
predictsSites$NH_5000.rs <- StdCenterPredictor(predictsSites$NH_5000)
predictsSites$NH_10000.rs <- StdCenterPredictor(predictsSites$NH_10000)

# rescaling abundance and log values
predictsSites <- RescaleAbundance(predictsSites)

# charlie added this line as later bits were throwing errors
predictsSites <- droplevels(predictsSites)


# some of the climate values are NA since they do not meet the thresholds
predictsSites <- predictsSites[!is.na(predictsSites$avg_temp), ]

# take a look at possible correlations between variables
plot(predictsSites$avg_temp, predictsSites$TmeanAnomaly,
     xlab = "Average temperature", 
     ylab = "Anomaly (difference between present and baseline)")
cor(predictsSites$avg_temp, predictsSites$TmeanAnomaly) # -0.42


plot(predictsSites$avg_temp, predictsSites$StdTmeanAnomaly,
     xlab = "Average temperature", 
     ylab = "Standardised climate anomaly")
cor(predictsSites$avg_temp, predictsSites$StdTmeanAnomaly) # 0.15

cor(predictsSites$TmeanAnomaly, predictsSites$StdTmeanAnomaly) # 0.21


# save the dataset
saveRDS(object = predictsSites,file = paste0(outDir,"PREDICTSSiteData.rds"))
#predictsSites <- readRDS(file = paste0(outDir,"PREDICTSSiteData.rds"))

##%######################################################%##
#                                                          #
####        Running the model selection process         ####
#                                                          #
##%######################################################%##




# 1. Abundance, mean anomaly

model_data <- predictsSites[!is.na(predictsSites$LogAbund), ] # 5735
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ] # 5735


MeanAnomalyModelAbund <- GLMERSelect(modelData = model_data,responseVar = "LogAbund",
                                fitFamily = "gaussian",fixedFactors = "UI2",
                                fixedTerms = list(StdTmeanAnomalyRS=1),
                                randomStruct = "(1|SS)+(1|SSB)",
                                fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyModelAbund, file = paste0(outDir, "/MeanAnomalyModelAbund.rdata"))



# 2. Richness, mean anomaly

model_data <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ] # 6069


MeanAnomalyModelRich <- GLMERSelect(modelData = model_data,responseVar = "Species_richness",
                                     fitFamily = "poisson",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("SSBS"))

# save model output
save(MeanAnomalyModelRich, file = paste0(outDir, "/MeanAnomalyModelRich.rdata"))

# 3. Abundance, max anomaly

model_data <- predictsSites[!is.na(predictsSites$LogAbund), ] # 5735
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ] # 5735


MaxAnomalyModelAbund <- GLMERSelect(modelData = model_data,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmaxAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save model output
save(MaxAnomalyModelAbund, file = paste0(outDir, "/MaxAnomalyModelAbund.rdata"))

# 4. Richness, max anomaly

model_data <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ] # 6069

MaxAnomalyModelRich <- GLMERSelect(modelData = model_data,responseVar = "Species_richness",
                                    fitFamily = "poisson",fixedFactors = "UI2",
                                    fixedTerms = list(StdTmaxAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)"),
                                    saveVars = c("SSBS"))

# save model output
save(MaxAnomalyModelRich, file = paste0(outDir, "/MaxAnomalyModelRich.rdata"))


#### save model output tables for use in supplementary information ####
# use function from sjPlot library to save neat versions of model output table.

tab_model(MeanAnomalyModelAbund$model, transform = NULL, file = paste0(outDir, "/AbunMeanAnom_output_table.html"))
summary(MeanAnomalyModelAbund$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelAbund$model) # check the R2 values 

tab_model(MeanAnomalyModelRich$model, transform = NULL, file = paste0(outDir, "/RichMeanAnom_output_table.html"))
summary(MeanAnomalyModelRich$model) # check the table against the outputs
R2GLMER(MeanAnomalyModelRich$model) # check the R2 values (function has issues with the richness models, so use output from Tim's formula)

#$conditional
#[1] 0.6385927

#$marginal
#[1] 0.01041383

tab_model(MaxAnomalyModelAbund$model, transform = NULL, file = paste0(outDir, "/AbunMaxAnom_output_table.html"))
summary(MaxAnomalyModelAbund$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelAbund$model) # check the R2 values 

tab_model(MaxAnomalyModelRich$model, transform = NULL, file = paste0(outDir, "/RichMaxAnom_output_table.html"))
summary(MaxAnomalyModelRich$model) # check the table against the outputs
R2GLMER(MaxAnomalyModelRich$model) # check the R2 values (function has issues with the richness models, so use output from Tim's formula)
#$conditional
#[1] 0.636048

#$marginal
#[1] 0.009140814

##%######################################################%##
#                                                          #
####                  save model stats                  ####
#                                                          #
##%######################################################%##

# save the stats tables

MeanAbunStats <- as.data.frame(MeanAnomalyModelAbund$stats)
MeanRichStats <- as.data.frame(MeanAnomalyModelRich$stats)
MaxAbunStats <- as.data.frame(MaxAnomalyModelAbund$stats)
MaxRichStats <- as.data.frame(MaxAnomalyModelRich$stats)


MeanAbunStats$significant <- NA
MeanRichStats$significant <- NA
MaxAbunStats$significant <- NA
MaxRichStats$significant <- NA


# function to check significance
checksig <- function(x){
  if(x <= 0.05){ 
    res <- "Yes" 
  } else { 
    res <- "No" }
  return(res)}

# add values to table
MeanAbunStats$significant <- sapply(X = MeanAbunStats$P, FUN = checksig)
MeanRichStats$significant <- sapply(X = MeanRichStats$P, FUN = checksig)
MaxAbunStats$significant <- sapply(X = MaxAbunStats$P, FUN = checksig)
MaxRichStats$significant <- sapply(X = MaxRichStats$P, FUN = checksig)


# save the stats tables
write.csv(MeanAbunStats, file = paste0(outDir, "/MeanAnomAbun_Stats.csv"), row.names = FALSE)
write.csv(MeanRichStats, file = paste0(outDir, "/MeanAnomRich_Stats.csv"), row.names = FALSE)
write.csv(MaxAbunStats, file = paste0(outDir, "/MaxAnomAbun_Stats.csv"), row.names = FALSE)
write.csv(MaxRichStats, file = paste0(outDir, "/MaxAnomRich_Stats.csv"), row.names = FALSE)




##%######################################################%##
#                                                          #
####                 Manuscript Figures                 ####
#                                                          #
##%######################################################%##

#### Plot results ####

#predictsSites <- readRDS(file = paste0(outDir,"PREDICTSSiteData.rds"))
#load(paste0(outDir, "/MeanAnomalyModelAbund.rdata"))
#load(paste0(outDir, "/MeanAnomalyModelRich.rdata"))
#load(paste0(outDir, "/MaxAnomalyModelAbund.rdata"))
#load(paste0(outDir, "/MaxAnomalyModelRich.rdata"))


# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)

# 
# pdf(file = paste0(outDir, "/Figure2_MeanAnom_Abun_Rich.pdf"), width = 8, height = 4)
# 
# par(mfrow=c(1,2))
# par(las=1)
# par(mgp=c(1.6,0.3,0)) # 
# par(mar=c(2.6,2.6,1,1)) # margins around plot
# par(tck=-0.01) # tick mark size
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
  
  
  # set factor levels
  nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
  
  
  # plot
  p1 <- ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
    geom_line(aes(col = UI2), size = 0.75) +
    geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
    geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
    scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
    scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
    theme_bw() + 
    scale_x_continuous(breaks = c(0,0.5, 1,1.5, 2), limits = c(0, 2)) +
    scale_y_continuous(breaks = c(-75, -50, -25, 0, 25, 50, 75, 100), limits = c(-75, 100)) +
    ylab("Change in total abundance (%)") +
    xlab("Standardised Temperature Anomaly") +
    #xlim(c(-1, 5)) +
    #ylim(c(-65, 60)) + 
    theme(aspect.ratio = 1, 
          title = element_text(size = 8, face = "bold"),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 7),
          legend.position = c(0.2, 0.8),
          legend.background = element_blank(), 
          legend.text = element_text(size = 6), 
          legend.title = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.2),
          panel.border = element_rect(size = 0.2), 
          axis.ticks = element_line(size = 0.2)) + 
    ggtitle("a")
  
  
  # # plottin with base R #
  # 
  # # set up plotting window
  # plot(-9e99,-9e99,xlim=c(min(nd$StdTmeanAnomaly),2),
  #      ylim=c(-70,80),
  #      xlab="Standardised Temperature Anomaly",ylab="Change in total abundance (%)", cex.lab = 0.8, cex.axis = 0.8)
  # 
  # title("a", adj = 0, cex.main = 1)
  # 
  # 
  # invisible(mapply(FUN = function(preds,col){
  #   
  #   preds <- na.omit(preds)
  #   
  #   X.Vec <- c(preds$StdTmeanAnomaly, max(preds$StdTmeanAnomaly), 
  #              rev(preds$StdTmeanAnomaly), min(preds$StdTmeanAnomaly))
  #   Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
  #              rev(preds$PredUpper), (preds$PredLower)[1])
  #   
  #   polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
  #   
  #   points(x = preds$StdTmeanAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
  #   
  # },split(nd,nd$UI2),c("#009E73", "#D55E00", "#E69F00", "#0072B2")))
  # 
  # # add some gridlines
  # abline(h=150,lty=1,col="#0000000C")
  # abline(h=100,lty=1,col="#0000000C")
  # abline(h=50,lty=1,col="#0000000C")
  # abline(h=0,lty=2,col="#030303")
  # abline(h=-50,lty=1,col="#0000000C")
  # abline(v=0,lty=1,col="#0000000C")
  # abline(v=0.5,lty=1,col="#0000000C")
  # abline(v=1,lty=1,col="#0000000C")
  # abline(v=1.5,lty=1,col="#0000000C")
  # abline(v=2,lty=1,col="#0000000C")
  # 
  # # add legend
  # legend(
  #   x = -0.1, y = 80,bty="n",
  #   legend = c("Primary","Secondary",
  #              "Agriculture_Low",
  #              "Agriculture_High"),
  #   col = c("#009E73", "#0072B2",
  #           "#E69F00", "#D55E00"),
  #   lty=1,lwd=2, cex = 0.8)
  # 

## now the species richness plot 
  nd2 <- expand.grid(
    StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
                          to = max(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
                          length.out = 100),
    UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
               levels = levels(MeanAnomalyModelRich$data$UI2)))
  
  # back transform the predictors
  nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
    transformedX = nd2$StdTmeanAnomalyRS,
    originalX = predictsSites$StdTmeanAnomaly)
  
  # set richness and abundance to 0 - to be predicted
  nd2$LogAbund <- 0
  nd2$Species_richness <- 0
  
  # reference for % difference = primary vegetation and positive anomaly closest to 0
  refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))))
  
  # adjust plot 1: mean anomaly and abundance
  
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
  
  # predict the results
  a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelRich$model,data = nd2, nIters = 10000)
  
  # back transform the abundance values
  a.preds.tmean <- exp(a.preds.tmean)
  
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
  
  # set factor levels
  nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
  
  # plot 
  p2 <- ggplot(data = nd2, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
    geom_line(aes(col = UI2), size = 0.75) +
    geom_ribbon(aes(ymin = nd2$PredLower, ymax = nd2$PredUpper, fill = UI2), alpha = 0.2) +
    geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
    scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
    scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
    scale_x_continuous(breaks = c(0,0.5, 1,1.5, 2), limits = c(0, 2)) +
    scale_y_continuous(breaks = c(-75, -50, -25, 0, 25, 50, 75, 100), limits = c(-75, 100)) +
    theme_bw() + 
    ylab("Change in species richness (%)") +
    xlab("Standardised Temperature Anomaly") +
    # xlim(c(-1, 5)) +
    # ylim(c(-65, 60)) + 
    theme(aspect.ratio = 1, 
          title = element_text(size = 8, face = "bold"),
          axis.text = element_text(size = 7),
          axis.title = element_text(size = 7),
          legend.position = "none",
          legend.text = element_text(size = 6), 
          legend.title = element_blank(), 
          panel.grid.minor = element_blank(),
          panel.grid.major = element_line(size = 0.2),
          panel.border = element_rect(size = 0.2), 
          axis.ticks = element_line(size = 0.2)) + 
    ggtitle("b")
  
 
  # combine plots
  cowplot::plot_grid(p1, p2)
  
  ggsave(filename = paste0(outDir, "Figure2_MeanAnom_Abun_Rich.pdf"), plot = last_plot(), width = 183, height = 100, units = "mm", dpi = 300)

  
  
  #  # plotting with base R#
  # 
  # # set up plotting window
  # plot(-9e99,-9e99,xlim=c(min(nd$StdTmeanAnomaly),2),
  #      ylim=c(-70,80),
  #      xlab="Standardised Temperature Anomaly",ylab="Change in species richness (%)", cex.lab = 0.8, cex.axis = 0.8)
  # 
  # title("b", adj = 0, cex.main = 1)
  # 
  # 
  # invisible(mapply(FUN = function(preds,col){
  #   
  #   preds <- na.omit(preds)
  #   
  #   X.Vec <- c(preds$StdTmeanAnomaly, max(preds$StdTmeanAnomaly), 
  #              rev(preds$StdTmeanAnomaly), min(preds$StdTmeanAnomaly))
  #   Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
  #              rev(preds$PredUpper), (preds$PredLower)[1])
  #   
  #   polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
  #   
  #   points(x = preds$StdTmeanAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
  #   
  # },split(nd,nd$UI2),c("#009E73", "#D55E00", "#E69F00", "#0072B2")))
  # 
  # # add some gridlines
  # abline(h=150,lty=1,col="#0000000C")
  # abline(h=100,lty=1,col="#0000000C")
  # abline(h=50,lty=1,col="#0000000C")
  # abline(h=0,lty=2,col="#030303")
  # abline(h=-50,lty=1,col="#0000000C")
  # abline(v=0,lty=1,col="#0000000C")
  # abline(v=0.5,lty=1,col="#0000000C")
  # abline(v=1,lty=1,col="#0000000C")
  # abline(v=1.5,lty=1,col="#0000000C")
  # abline(v=2,lty=1,col="#0000000C")
  # 
  # # add legend
  # # legend(
  # #   x = 0,y = 90,bty="n",
  # #   legend = c("Primary","Secondary",
  # #              "Agriculture_Low",
  # #              "Agriculture_High"),
  # #   col = c("#009E73", "#0072B2",
  # #           "#E69F00", "#D55E00"),
  # #   lty=1,lwd=2, cex = 0.8)
  # 
  # 
  # dev.off()
  # 
  # 

  #### Extended data fig 2 - Max Anom ####
  
  #pdf(file = paste0(outDir,"LUClimateAnomalyInteractions_Mean_anom.pdf"),width = 17.5/2.54,height = 16/2.54)
  #pdf(file = paste0(outDir,"LUClimateAnomalyInteractions_Mean_anom.pdf"),width = 17.5/2.54,height = 16/2.54)
  
  # #pdf(file = paste0(outDir,"Extended_Data2_MaxAnom.pdf"),width = 8,height = 4)
  # jpeg(file = paste0(outDir,"Extended_Data2_MaxAnom.jpeg"),width = 8,height = 4, units = "mm")
  # 
  # par(mfrow=c(1,2))
  # par(las=1)
  # par(mgp=c(1.6,0.3,0)) # 
  # par(mar=c(2.6,2.6,1,1)) # margins around plot
  # par(tck=-0.01) # tick mark size
  # #par(pty="s") # set plot type to be square
  # 
  
  # # create matrix for predictions
  # nd <- expand.grid(
  #   StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
  #                         to = max(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
  #                         length.out = 100),
  #   UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
  #              levels = levels(MeanAnomalyModelAbund$data$UI2)))
  # 
  # # back transform the predictors
  # nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  #   transformedX = nd$StdTmeanAnomalyRS,
  #   originalX = predictsSites$StdTmeanAnomaly)
  # 
  # # set richness and abundance to 0 - to be predicted
  # nd$LogAbund <- 0
  # nd$Species_richness <- 0
  # 
  # # reference for % difference = primary vegetation and positive anomaly closest to 0
  # refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))
  # 
  # # Quantiles for each land use
  # QPV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  #   MeanAnomalyModelAbund$data$UI2=="Primary vegetation"],
  #   probs = exclQuantiles)
  # QSV <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  #   MeanAnomalyModelAbund$data$UI2=="Secondary vegetation"],
  #   probs = exclQuantiles)
  # QAL <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  #   MeanAnomalyModelAbund$data$UI2=="Agriculture_Low"],
  #   probs = exclQuantiles)
  # QAH <- quantile(x = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS[
  #   MeanAnomalyModelAbund$data$UI2=="Agriculture_High"],
  #   probs = exclQuantiles)
  # 
  # if(!is.null(MeanAnomalyModelAbund$model)){
  #   
  #   # predict the results
  #   a.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelAbund$model,data = nd)
  #   
  #   # back transform the abundance values
  #   a.preds.tmean <- exp(a.preds.tmean)-0.01
  #   
  #   # convert to relative to reference
  #   a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')
  #   
  #   # remove anything above and below the quantiles
  #   a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1])),] <- NA
  #   a.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2])),] <- NA
  #   a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS < QSV[1])),] <- NA
  #   a.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS > QSV[2])),] <- NA
  #   a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS < QAL[1])),] <- NA
  #   a.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS > QAL[2])),] <- NA
  #   a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS < QAH[1])),] <- NA
  #   a.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS > QAH[2])),] <- NA
  #   
  #   # Get the median, upper and lower quants for the plot
  #   nd$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
  #                            FUN = median,na.rm=TRUE))*100)-100
  #   nd$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
  #                           FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
  #   nd$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
  #                           FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
  #   
  #   # plot #
  #   
  #   # set up plotting window
  #   plot(-9e99,-9e99,xlim=c(min(nd$StdTmeanAnomaly),max(nd$StdTmeanAnomaly)),
  #        ylim=c(min(nd$PredLower,na.rm = TRUE),max(nd$PredUpper,na.rm = TRUE)),
  #        xlab="Mean temperature anomaly",ylab="Abundance (%)")
  #   
  #   invisible(mapply(FUN = function(preds,col){
  #     
  #     preds <- na.omit(preds)
  #     
  #     X.Vec <- c(preds$StdTmeanAnomaly, max(preds$StdTmeanAnomaly), 
  #                rev(preds$StdTmeanAnomaly), min(preds$StdTmeanAnomaly))
  #     Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
  #                rev(preds$PredUpper), (preds$PredLower)[1])
  #     
  #     polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
  #     
  #     points(x = preds$StdTmeanAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
  #     
  #   },split(nd,nd$UI2),c("#009E73", "#D55E00", "#E69F00", "#0072B2")))
  #   
  #   # add some gridlines
  #   abline(h=150,lty=1,col="#00000022")
  #   abline(h=100,lty=1,col="#00000022")
  #   abline(h=50,lty=1,col="#00000022")
  #   abline(h=0,lty=1,col="#00000022")
  #   abline(h=-50,lty=1,col="#00000022")
  #   abline(v=0,lty=1,col="#00000022")
  #   abline(v=1,lty=1,col="#00000022")
  #   abline(v=2,lty=1,col="#00000022")
  #   
  #   # add legend
  #   legend(
  #     x = -0.6,y = 115,bty="n",
  #     legend = c("Primary","Secondary",
  #                "Agriculture_extensive",
  #                "Agriculture_intensive"),
  #     col = c("#009E73", "#0072B2",
  #             "#E69F00", "#D55E00"),
  #     lty=1,lwd=2)
  #   
  #   # add title
  #   #title("c.", adj = 0)
  #   
  #   p1 <- recordPlot()
  #   
  # } else {
  #   frame()
  # }
  # 
  # 
  # QPV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  #   MeanAnomalyModelRich$data$UI2=="Primary vegetation"],
  #   probs = exclQuantiles)
  # QSV <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  #   MeanAnomalyModelRich$data$UI2=="Secondary vegetation"],
  #   probs = exclQuantiles)
  # QAL <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  #   MeanAnomalyModelRich$data$UI2=="Agriculture_Low"],
  #   probs = exclQuantiles)
  # QAH <- quantile(x = MeanAnomalyModelRich$data$StdTmeanAnomalyRS[
  #   MeanAnomalyModelRich$data$UI2=="Agriculture_High"],
  #   probs = exclQuantiles)
  # 
  # if(!is.null(MeanAnomalyModelRich$model)){
  #   
  #   s.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelRich$model,data = nd)
  #   s.preds.tmean <- exp(s.preds.tmean)
  #   
  #   s.preds.tmean <- sweep(x = s.preds.tmean,MARGIN = 2,STATS = s.preds.tmean[refRow,],FUN = '/')
  #   
  #   s.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1])),] <- NA
  #   s.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2])),] <- NA
  #   s.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS < QSV[1])),] <- NA
  #   s.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS > QSV[2])),] <- NA
  #   s.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS < QAL[1])),] <- NA
  #   s.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS > QAL[2])),] <- NA
  #   s.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS < QAH[1])),] <- NA
  #   s.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS > QAH[2])),] <- NA
  #   
  #   nd$PredMedian <- ((apply(X = s.preds.tmean,MARGIN = 1,
  #                            FUN = median,na.rm=TRUE))*100)-100
  #   nd$PredUpper <- ((apply(X = s.preds.tmean,MARGIN = 1,
  #                           FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
  #   nd$PredLower <- ((apply(X = s.preds.tmean,MARGIN = 1,
  #                           FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
  #   
  #   plot(-9e99,-9e99,xlim=c(min(nd$StdTmeanAnomaly),max(nd$StdTmeanAnomaly)),
  #        ylim=c(min(nd$PredLower,na.rm = TRUE),max(nd$PredUpper,na.rm = TRUE)),
  #        xlab="Mean temperature anomaly",ylab="Richness (%)")
  #   
  #   invisible(mapply(FUN = function(preds,col){
  #     
  #     preds <- na.omit(preds)
  #     
  #     X.Vec <- c(preds$StdTmeanAnomaly, max(preds$StdTmeanAnomaly), 
  #                rev(preds$StdTmeanAnomaly), min(preds$StdTmeanAnomaly))
  #     Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
  #                rev(preds$PredUpper), (preds$PredLower)[1])
  #     
  #     polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
  #     
  #     points(x = preds$StdTmeanAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
  #     
  #   },split(nd,nd$UI2),c("#009E73", "#D55E00", "#E69F00", "#0072B2")))
  #   
  #   
  #   abline(h=150,lty=1,col="#00000022")
  #   abline(h=100,lty=1,col="#00000022")
  #   abline(h=50,lty=1,col="#00000022")
  #   abline(h=0,lty=1,col="#00000022")
  #   abline(h=-50,lty=1,col="#00000022")
  #   abline(v=0,lty=1,col="#00000022")
  #   abline(v=1,lty=1,col="#00000022")
  #   abline(v=2,lty=1,col="#00000022")
  #   
  #   
  #   p2 <- recordPlot()
  #   
  # } else {
  #   frame()
  # }
  
  
  ### Extended Data 2 - maximum anomaly ###
  
  nd <- expand.grid(
    StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
                         to = max(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
                         length.out = 200),
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
  
 # if(!is.null(MaxAnomalyModelAbund$model)){
    
    a.preds.tmax <- PredictGLMERRandIter(model = MaxAnomalyModelAbund$model,data = nd, nIters = 10000)
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
    
    
    
    # set factor levels
    nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
    
    
    # plot
    p1 <- ggplot(data = nd, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
      geom_line(aes(col = UI2), size = 0.75) +
      geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
      geom_hline(yintercept = 0, lty = "dashed") +
      scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
      scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
      theme_bw() + 
      scale_x_continuous(breaks = c(-1,0,1,2,3,4,5), limits = c(-1, 5)) +
      scale_y_continuous(breaks = c(-60, -40, -20, 0, 20, 40, 60), limits = c(-65, 60)) +
      labs(fill = "% NH", col = "% NH") + 
      ylab("Change in total abundance (%)") +
      xlab("Maximum Temperature Anomaly") +
      #xlim(c(-1, 5)) +
      #ylim(c(-65, 60)) + 
      theme(aspect.ratio = 1, 
            title = element_text(size = 8, face = "bold"),
            axis.text = element_text(size = 7),
            legend.position = c(0.2, 0.8),
            legend.text = element_text(size = 6), 
            legend.title = element_blank(), 
            legend.background = element_blank()) + 
      ggtitle("a")
    
    
    # 
    # 
    # 
    # plot(-9e99,-9e99,xlim=c(-1,max(nd$StdTmaxAnomaly)),
    #      ylim=c(-60,60),
    #      xlab="Maximum temperature anomaly",ylab="Change in total abundance (%)", cex.lab = 0.8, cex.axis = 0.8)
    # 
    # title("a", adj = 0, cex.main = 1)
    # 
    # invisible(mapply(FUN = function(preds,col){
    #   
    #   preds <- na.omit(preds)
    #   
    #   X.Vec <- c(preds$StdTmaxAnomaly, max(preds$StdTmaxAnomaly), 
    #              rev(preds$StdTmaxAnomaly), min(preds$StdTmaxAnomaly))
    #   Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
    #              rev(preds$PredUpper), (preds$PredLower)[1])
    #   
    #   polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
    #   
    #   points(x = preds$StdTmaxAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
    #   
    # },split(nd,nd$UI2),c("#009E73", "#D55E00", "#E69F00", "#0072B2")))
    # 
    # abline(h=40,lty=1,col="#0000000C")
    # abline(h=20,lty=1,col="#0000000C")
    # abline(h=0,lty=2,col="#030303")
    # abline(h=-20,lty=1,col="#0000000C")
    # abline(h=-40,lty=1,col="#0000000C")
    # abline(h=-60,lty=1,col="#0000000C")
    # abline(v=-1,lty=1,col="#0000000C")
    # abline(v=0,lty=1,col="#0000000C")
    # abline(v=1,lty=1,col="#0000000C")
    # abline(v=2,lty=1,col="#0000000C")
    # abline(v=3,lty=1,col="#0000000C")
    # abline(v=4,lty=1,col="#0000000C")
    # abline(v=5,lty=1,col="#0000000C")
    # 
    # legend(
    #   x = -1,y = 60 ,bty="n",
    #   legend = c("Primary","Secondary",
    #              "Agriculture_Low",
    #              "Agriculture_High"),
    #   col = c("#009E73", "#0072B2",
    #           "#E69F00", "#D55E00"),
    #   lty=1,lwd=2, cex = 0.8)
    # 
    # #p3 <- recordPlot()
    # 
    
    
    # Organise table to values to be predicted
    
    nd2 <- expand.grid(
      StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
                           to = max(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
                           length.out = 200),
      UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
                 levels = levels(MaxAnomalyModelAbund$data$UI2)))
    nd2$StdTmaxAnomaly <- BackTransformCentreredPredictor(
      transformedX = nd2$StdTmaxAnomalyRS,
      originalX = predictsSites$StdTmaxAnomaly)
    nd2$LogAbund <- 0
    nd2$Species_richness <- 0
    
    # set a reference row for percentage change
    refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomaly==min(abs(nd2$StdTmaxAnomaly))))
  
  # determine the quantiles of the available data
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
  
    # predict values
    s.preds.tmax <- PredictGLMERRandIter(model = MaxAnomalyModelRich$model,data = nd2, nIters = 10000)
    s.preds.tmax <- exp(s.preds.tmax)
    
    # convert to percentage change relative to reference row
    s.preds.tmax <- sweep(x = s.preds.tmax,MARGIN = 2,STATS = s.preds.tmax[refRow,],FUN = '/')
    
    # exclude extreme ends of data
    s.preds.tmax[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomalyRS < QPV[1])),] <- NA
    s.preds.tmax[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomalyRS > QPV[2])),] <- NA
    s.preds.tmax[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmaxAnomalyRS < QSV[1])),] <- NA
    s.preds.tmax[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmaxAnomalyRS > QSV[2])),] <- NA
    s.preds.tmax[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmaxAnomalyRS < QAL[1])),] <- NA
    s.preds.tmax[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmaxAnomalyRS > QAL[2])),] <- NA
    s.preds.tmax[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmaxAnomalyRS < QAH[1])),] <- NA
    s.preds.tmax[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmaxAnomalyRS > QAH[2])),] <- NA
    
    # get the median and quantiles for plotting
    nd2$PredMedian <- ((apply(X = s.preds.tmax,MARGIN = 1,
                             FUN = median,na.rm=TRUE))*100)-100
    nd2$PredUpper <- ((apply(X = s.preds.tmax,MARGIN = 1,
                            FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
    nd2$PredLower <- ((apply(X = s.preds.tmax,MARGIN = 1,
                            FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
    
    # set factor levels
    nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))
    
    # plot 
    p2 <- ggplot(data = nd2, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
      geom_line(aes(col = UI2), size = 0.75) +
      geom_ribbon(aes(ymin = nd2$PredLower, ymax = nd2$PredUpper, fill = UI2), alpha = 0.2) +
      geom_hline(yintercept = 0, lty = "dashed") +
      scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
      scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
      scale_x_continuous(breaks = c(-1,0,1,2,3,4,5), limits = c(-1, 5)) +
      scale_y_continuous(breaks = c(-60, -40, -20, 0, 20, 40, 60), limits = c(-65, 60)) +
      theme_bw() + 
      labs(fill = "% NH", col = "% NH") + 
      ylab("Change in species richness (%)") +
      xlab("Maximum Temperature Anomaly") +
      # xlim(c(-1, 5)) +
      # ylim(c(-65, 60)) + 
      theme(aspect.ratio = 1, 
            title = element_text(size = 8, face = "bold"),
            axis.text = element_text(size = 7),
            legend.position = "none",
            legend.text = element_text(size = 6), 
            legend.title = element_blank()) + 
      ggtitle("b")
    
    # combine plots
    cowplot::plot_grid(p1, p2)
    
    ggsave(filename = paste0(outDir, "Extended_Data2_MaxAnom.jpeg"), plot = last_plot(), width = 183, height = 150, units = "mm", dpi = 300)
    
    
  #   
  #   
  #   plot(-9e99,-9e99,xlim=c(-1,max(nd$StdTmaxAnomaly)),
  #        ylim=c(-60,60),
  #        xlab="Maximum temperature anomaly",ylab="Change in species richness (%)", cex.lab = 0.8, cex.axis = 0.8)
  #   
  #   title("b", adj = 0, cex.main = 1)
  #   
  #   invisible(mapply(FUN = function(preds,col){
  #     
  #     preds <- na.omit(preds)
  #     
  #     X.Vec <- c(preds$StdTmaxAnomaly, max(preds$StdTmaxAnomaly), 
  #                rev(preds$StdTmaxAnomaly), min(preds$StdTmaxAnomaly))
  #     Y.Vec <- c(preds$PredLower, tail(preds$PredUpper, 1), 
  #                rev(preds$PredUpper), (preds$PredLower)[1])
  #     
  #     polygon(x = X.Vec,y = Y.Vec,col=paste0(col,"33"),border=NA)
  #     
  #     points(x = preds$StdTmaxAnomaly,y = preds$PredMedian,type="l",lwd=2,col=paste0(col))
  #     
  #   },split(nd,nd$UI2),c("#009E73", "#D55E00", "#E69F00", "#0072B2")))
  #   
  #   abline(h=40,lty=1,col="#0000000C")
  #   abline(h=20,lty=1,col="#0000000C")
  #   abline(h=0,lty=2,col="#030303")
  #   abline(h=-20,lty=1,col="#0000000C")
  #   abline(h=-40,lty=1,col="#0000000C")
  #   abline(h=-60,lty=1,col="#0000000C")
  #   abline(v=-1,lty=1,col="#0000000C")
  #   abline(v=0,lty=1,col="#0000000C")
  #   abline(v=1,lty=1,col="#0000000C")
  #   abline(v=2,lty=1,col="#0000000C")
  #   abline(v=3,lty=1,col="#0000000C")
  #   abline(v=4,lty=1,col="#0000000C")
  #   abline(v=5,lty=1,col="#0000000C")
  #   #p4 <- recordPlot()
  #   
  #   
  # 
  # 
  # invisible(dev.off())
    
    
    t.end <- Sys.time()
    
    print(round(t.end - t.start,0))
    
    sink()
    