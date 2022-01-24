##%######################################################%##
#                                                          #
####        14. Additional test: Remove outliers        ####
#                                                          #
##%######################################################%##


# in this script, I will check for outliers in the models, remove
# these from the dataset and rerun the models to see if they have
# an effect on the results.

 
### The check_outliers function, from the performance package using cook's distance, doesn't find any outliers. ###


rm(list = ls())

# libraries
library(performance)
source("Functions.R")
library(influence.ME)
library(StatisticalModels)


# directories
moddir <- "6_RunLUClimateModels/"
outdir <- "14_Additional_Tests/"



# load in models
load(paste0(moddir, "MeanAnomalyModelAbund.rdata")) # MeanAnomalyModelAbund
load(paste0(moddir, "MeanAnomalyModelRich.rdata"))  # MeanAnomalyModelRich
load(paste0(moddir, "MaxAnomalyModelAbund.rdata"))  # MaxAnomalyModelAbund
load(paste0(moddir, "MaxAnomalyModelRich.rdata"))   # MaxAnomalyModelRich


# 1. Abundance Mean anom

# load the datasets in to help the functions
modelData <- MeanAnomalyModelAbund$data

# check for outliers
check_outliers(MeanAnomalyModelAbund$model, method = "cook")

# OK: No outliers detected.


# 2. Richness Mean anom

# load the datasets in to help the functions
modelData <- MeanAnomalyModelRich$data

# check for outliers
check_outliers(MeanAnomalyModelRich$model, method = "cook")

# OK: No outliers detected.
#Warning messages:the hat matrix may not make sense for GLMMs

# 3. Abundance Max anom

# load the datasets in to help the functions
modelData <- MaxAnomalyModelAbund$data

# check for outliers
check_outliers(MaxAnomalyModelAbund$model, method = "cook")

# OK: No outliers detected.


# 4. Richness Mean anom

# load the datasets in to help the functions
modelData <- MaxAnomalyModelRich$data

# check for outliers
check_outliers(MaxAnomalyModelRich$model, method = "cook")

# OK: No outliers detected.
#Warning messages:the hat matrix may not make sense for GLMMs



#### using an alternative package to check outliers, just in case ####

# 1. abundance mean anomaly

# result1 <- cooks.distance(MeanAnomalyModelAbund$model, sort = T)
# range(result1)
# summary(result1 < 1)# all true
# result1[result1 > 0.4]

modelData <- MeanAnomalyModelAbund$data

length(unique(modelData$SS)) # 243

maxIters=10000
optimizer="bobyqa"
alt.est.a <- influence(MeanAnomalyModelAbund$model, "SS")
plot(alt.est.a, which = "cook", sort = T)
result1 <- cooks.distance(alt.est.a, sort = T)
# looks like there are some with larger distances
# 2 with higher values (previewed the plot and made it large to see the names)
# SC1_2011__Meijer 1, CC1_2007__Ewers 1

# 2. Richness mean anomaly

# result2 <- cooks.distance(MeanAnomalyModelRich$model)
# range(result2)
# summary(result2 < 1)# all true
# result2[result2 > 0.4]

modelData <- MeanAnomalyModelRich$data

length(unique(modelData$SS)) # 263

# aditional info needed for the influence function
maxIters=10000
optimizer="bobyqa"
fitFamily = "poisson"
alt.est.a <- influence(MeanAnomalyModelRich$model, "SS")
plot(alt.est.a, which = "cook", sort = T)
# looks like there are some with larger distances
# 3 with higher values (previewed the plot and made it large to see the names)
# SC1_2011__Meijer 1, CC1_2007__Ewers 1, HP1_2014__Gray 1


# 3. Abundance Max Anomaly

  
# result3 <- cooks.distance(MaxAnomalyModelAbund$model)
# range(result3)
# summary(result3 < 1)# all true
# result3[result3 > 0.4]

modelData <- MaxAnomalyModelAbund$data

length(unique(modelData$SS)) # 243
  
# aditional info needed for the influence function
maxIters=10000
optimizer="bobyqa"
#fitFamily = "poisson"
alt.est.a <- influence(MaxAnomalyModelAbund$model, "SS")
plot(alt.est.a, which = "cook", sort = T)
# looks like there are some with larger distances
# 3 with higher values (previewed the plot and made it large to see the names)
# DI1_2005__Tylianakis 1, CC1_2007__Ewers 1, SC1_2011__Meijer 1


# 4. Richness Max Anomaly

# result4 <- cooks.distance(MaxAnomalyModelRich$model)
# range(result4)
# summary(result4 < 1)# all true
# result4[result4 > 0.4]

modelData <- MaxAnomalyModelRich$data

length(unique(modelData$SS)) # 263

# aditional info needed for the influence function
maxIters=10000
optimizer="bobyqa"
fitFamily = "poisson"
alt.est.a <- influence(MaxAnomalyModelRich$model, "SS")
plot(alt.est.a, which = "cook", sort = T)
# looks like there are some with larger distances
# 4 with higher values (previewed the plot and made it large to see the names)
# CC1_2007__Ewers 1, SC1_2011__Meijer 1, HP1_2014__Gray 1, DI1_2005__Tylianakis 1




##%######################################################%##
#                                                          #
####      Rerun the models removing these outliers      ####
#                                                          #
##%######################################################%##

datadir <- "6_RunLUClimateModels/"

predictsSites <- readRDS(file = paste0(datadir,"PREDICTSSiteData.rds"))



# 1. Abundance, mean anomaly

# list of studies that had larger Cook's distance values
outliers <- c("SC1_2011__Meijer 1", "CC1_2007__Ewers 1")

# subset the data to exclude the potential outliers
model_data <- predictsSites[!predictsSites$SS %in% outliers, ] # 5710 rows

model_data <- model_data[!is.na(model_data$LogAbund), ] #5374 rows

length(unique(model_data$SS)) # 242

MeanAnomalyModelAbund <- GLMERSelect(modelData = model_data,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(MeanAnomalyModelAbund, file = paste0(outdir, "/MeanAnomalyModelAbund_rmout.rdata"))
#load(file = paste0(outdir, "/MeanAnomalyModelAbund_rmout.rdata"))



# 2. Richness, mean anomaly


# list of studies that had larger Cook's distance values
outliers <- c("SC1_2011__Meijer 1", "CC1_2007__Ewers 1", "HP1_2014__Gray 1")


# subset the data to exclude the potential outliers
model_data <- predictsSites[!predictsSites$SS %in% outliers, ] # 5419 rows

model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ]

nrow(predictsSites[is.na(predictsSites$StdTmeanAnomalyRS), ])

length(unique(model_data$SS)) # 260

MeanAnomalyModelRich <- GLMERSelect(modelData = model_data,responseVar = "Species_richness",
                                    fitFamily = "poisson",fixedFactors = "UI2",
                                    fixedTerms = list(StdTmeanAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                    saveVars = c("SSBS"))

# save model output
save(MeanAnomalyModelRich, file = paste0(outdir, "/MeanAnomalyModelRich_rmout.rdata"))
#load(file = paste0(outdir, "/MeanAnomalyModelRich_rmout.rdata"))

MeanAnomalyModelRich$model


# 3. Abundance, max anomaly

# list of studies that had larger Cook's distance values
outliers <- c("DI1_2005__Tylianakis 1", "CC1_2007__Ewers 1", "SC1_2011__Meijer 1")

# subset the data to exclude the potential outliers
model_data <- predictsSites[!predictsSites$SS %in% outliers, ] # 5662 rows

model_data <- model_data[!is.na(model_data$LogAbund), ] # 5326 rows

length(unique(model_data$SS)) # 24


MaxAnomalyModelAbund <- GLMERSelect(modelData = model_data,responseVar = "LogAbund",
                                    fitFamily = "gaussian",fixedFactors = "UI2",
                                    fixedTerms = list(StdTmaxAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)",
                                    fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)"),
                                    saveVars = c("SSBS"))

# save model output
save(MaxAnomalyModelAbund, file = paste0(outdir, "/MaxAnomalyModelAbund_rmout.rdata"))
#load(file = paste0(outdir, "/MaxAnomalyModelAbund_rmout.rdata"))


# 4. Richness, max anomaly
# list of studies that had larger Cook's distance values
outliers <- c("CC1_2007__Ewers 1", "SC1_2011__Meijer 1", "HP1_2014__Gray 1", "DI1_2005__Tylianakis 1")


model_data <- predictsSites[!predictsSites$SS %in% outliers, ] # 5397 rows

model_data <- model_data[!is.na(model_data$StdTmaxAnomalyRS), ] # 5371

length(unique(model_data$SS)) # 259

model_data <- model_data[, c("Species_richness", "UI2", "StdTmaxAnomalyRS", "SS", "SSB", "SSBS")]

MaxAnomalyModelRich <- GLMERSelect(modelData = model_data,responseVar = "Species_richness",
                                    fitFamily = "poisson",fixedFactors = "UI2",
                                    fixedTerms = list(StdTmaxAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)"))

# save model output
save(MaxAnomalyModelRich, file = paste0(outdir, "/MaxAnomalyModelRich_rmout.rdata"))
#load(file = paste0(outdir, "/MaxAnomalyModelRich_rmout.rdata"))

MaxAnomalyModelRich$model


##%######################################################%##
#                                                          #
####                  Replot Figure 2                   ####
#                                                          #
##%######################################################%##


# set quantiles of predicted result to be presented in the plots
exclQuantiles <- c(0.025,0.975)

# 
# pdf(file = paste0(outdir, "/Figure2_MeanAnom_Abun_Rich_outliersrem.pdf"), width = 8, height = 4)
# 
# par(mfrow=c(1,2))
# par(las=1)
# par(mgp=c(1.6,0.3,0)) # 
# par(mar=c(2.6,2.6,1,1)) # margins around plot
# par(tck=-0.01) # tick mark size
# #par(pty="s") # set plot type to be square

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
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0, 0.5,1, 1.5,2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-80, -60, -40, -20, 0, 20, 40, 60), limits = c(-80, 70)) +
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
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
# # plot #
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


nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))


# plot
p2 <- ggplot(data = nd2, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 0.75) +
  geom_ribbon(aes(ymin = nd2$PredLower, ymax = nd2$PredUpper, fill = UI2), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00","#D55E00")) +
  theme_bw() + 
  scale_x_continuous(breaks = c(0, 0.5,1, 1.5,2), limits = c(0, 2)) +
  scale_y_continuous(breaks = c(-50, -25, 0, 25, 50, 75, 100), limits = c(-50, 110)) +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank(), 
        legend.background = element_blank()) + 
  ggtitle("b")

# combine plots
plot_grid(p1, p2)

ggsave(filename = paste0(outdir, "Extended_Data3_MeanAnom_outrem.jpeg"), plot = last_plot(), width = 183, height = 150, units = "mm", dpi = 300)


# 
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



#### max anom figure ####
# 
# pdf(file = paste0(outdir,"Extended_Data3_MaxAnom_outliersrem.pdf"),width = 8,height = 4)
# 
# par(mfrow=c(1,2))
# par(las=1)
# par(mgp=c(1.6,0.3,0)) # 
# par(mar=c(2.6,2.6,1,1)) # margins around plot
# par(tck=-0.01) # tick mark size
# #par(pty="s") # set plot type to be square


nd <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
                       to = max(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
                       length.out = 250),
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
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5), limits = c(-1, 6)) +
  scale_y_continuous(breaks = c(-80, -60, -40, -20, 0, 20, 40, 60, 80), limits = c(-80, 85)) +
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



nd2 <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
                       to = max(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
                       length.out = 250),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxAnomalyModelAbund$data$UI2)))

nd2$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)

nd2$LogAbund <- 0
nd2$Species_richness <- 0

refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomaly==min(abs(nd2$StdTmaxAnomaly))))



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


s.preds.tmax <- PredictGLMERRandIter(model = MaxAnomalyModelRich$model,data = nd2, nIters = 10000)
s.preds.tmax <- exp(s.preds.tmax)

s.preds.tmax <- sweep(x = s.preds.tmax,MARGIN = 2,STATS = s.preds.tmax[refRow,],FUN = '/')

s.preds.tmax[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomalyRS < QPV[1])),] <- NA
s.preds.tmax[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomalyRS > QPV[2])),] <- NA
s.preds.tmax[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmaxAnomalyRS < QSV[1])),] <- NA
s.preds.tmax[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmaxAnomalyRS > QSV[2])),] <- NA
s.preds.tmax[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmaxAnomalyRS < QAL[1])),] <- NA
s.preds.tmax[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmaxAnomalyRS > QAL[2])),] <- NA
s.preds.tmax[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmaxAnomalyRS < QAH[1])),] <- NA
s.preds.tmax[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmaxAnomalyRS > QAH[2])),] <- NA

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
  theme_bw() + 
  scale_x_continuous(breaks = c(-1, 0, 1, 2, 3, 4, 5), limits = c(-1, 6)) +
  scale_y_continuous(breaks = c(-80, -60, -40, -20, 0, 20, 40, 60, 80), limits = c(-80, 85)) +
  ylab("Change in species richness (%)") +
  xlab("Maximum Temperature Anomaly") +
  #xlim(c(-1, 5)) +
  #ylim(c(-65, 60)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        legend.position = "none",
        legend.text = element_text(size = 6), 
        legend.title = element_blank(), 
        legend.background = element_blank()) + 
  ggtitle("b")


# combine plots
plot_grid(p1, p2)

ggsave(filename = paste0(outdir, "Extended_Data4_MaxAnom_outrem.jpeg"), plot = last_plot(), width = 183, height = 150, units = "mm", dpi = 300)


# 
# plot(-9e99,-9e99,xlim=c(-1,max(nd$StdTmaxAnomaly)),
#      ylim=c(-60,60),
#      xlab="Maximum temperature anomaly",ylab="Change in species richness (%)", cex.lab = 0.8, cex.axis = 0.8)
# 
# title("b", adj = 0, cex.main = 1)
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
# #p4 <- recordPlot()
# 
# 
# 
# 
# invisible(dev.off())
# 



