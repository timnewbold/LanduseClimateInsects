##%######################################################%##
#                                                          #
####                Run models inc Realm                ####
#                                                          #
##%######################################################%##

# running models from script 6 and 7 but with the inclusion of realm to 
# check whether responses differ between tropical and temperate regions


# load libraries
library(StatisticalModels)
library(predictsFunctions)
source("Functions.R")
library(ggplot2)
library(cowplot)


# directories
predictsDataDir <- "6_RunLUClimateModels/"
outDir <- "9_RunModels_Tropical/"
dir.create(outDir)

# read in PREDICTS data
predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSiteData.rds"))

# categorise sites into temperate/tropical
# get the tropical values
predictsSites$Tropical <- NA

predictsSites[predictsSites$Latitude > -23.44 & predictsSites$Latitude < 23.44, 'Tropical'] <- "Tropical"

# label the remaining as temperate
predictsSites[is.na(predictsSites$Tropical), 'Tropical'] <- "Temperate"

# set as a factor
predictsSites$Tropical <- as.factor(predictsSites$Tropical)
# levels: Temperate Tropical

table(predictsSites$Tropical)

# Temperate  Tropical 
#    5455      2339 

table(predictsSites$UI2, predictsSites$Tropical)

#                        Temperate Tropical
# Primary vegetation         951      613
# Agriculture_High          1350      429
# Agriculture_Low           1237      265
# Secondary vegetation      1048      435


# 1. Abundance, mean anomaly
MeanAnomalyModelAbund <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",
                                     fixedFactors = c("UI2", "Tropical"),
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)",
                                                           "UI2:Tropical",
                                                           "Tropical:poly(StdTmeanAnomalyRS,1)",
                                                           "Tropical:poly(StdTmeanAnomalyRS,1):UI2"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

summary(MeanAnomalyModelAbund$model)


# selected model:
# LogAbund ~ UI2 + + poly(StdTmeanAnomalyRS, 1) + Tropical +
# UI2:poly(StdTmeanAnomalyRS, 1) + UI2:Tropical +
# Tropical:poly(StdTmeanAnomalyRS, 1):UI2 + 
# (1 | SS) + (1 | SSB)


# save the model output
save(MeanAnomalyModelAbund, file = paste0(outDir, "/MeanAnomalyModelAbund_3way.rdata"))

# 2. Richness, mean anomaly
MeanAnomalyModelRich <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
                                    fitFamily = "poisson",
                                    fixedFactors = c("UI2", "Tropical"),
                                    fixedTerms = list(StdTmeanAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)",
                                                          "UI2:Tropical",
                                                          "Tropical:poly(StdTmeanAnomalyRS,1)",
                                                          "Tropical:poly(StdTmeanAnomalyRS,1):UI2"),
                                    saveVars = c("Total_abundance", "SSBS", "NH_3000"))

summary(MeanAnomalyModelRich$model)

# selected model:
# Species_richness ~ UI2 + Tropical + poly(StdTmeanAnomalyRS, 1) +  
# UI2:poly(StdTmeanAnomalyRS, 1) + UI2:Tropical + Tropical:poly(StdTmeanAnomalyRS, 1) + 
# Tropical:poly(StdTmeanAnomalyRS, 1):UI2 + 
# (1 | SS) +  (1 | SSB) + (1 | SSBS)

# Model failed to converge with max|grad| = 0.0187175 (tol = 0.001, component 1)


# save model output
save(MeanAnomalyModelRich, file = paste0(outDir, "/MeanAnomalyModelRich_3way.rdata"))


# 3. Abundance, max anomaly
MaxAnomalyModelAbund <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                    fitFamily = "gaussian",
                                    fixedFactors = c("UI2", "Tropical"),
                                    fixedTerms = list(StdTmaxAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)",
                                    fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)",
                                                          "UI2:Tropical",
                                                          "Tropical:poly(StdTmaxAnomalyRS,1)",
                                                          "Tropical:poly(StdTmaxAnomalyRS,1):UI2"),
                                    saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

summary(MaxAnomalyModelAbund$model)

# selected model:
# LogAbund ~ UI2 + UI2:Tropical + Tropical + (1 | SS) + (1 | SSB)

# save model output
save(MaxAnomalyModelAbund, file = paste0(outDir, "/MaxAnomalyModelAbund_3way.rdata"))


# 4. Richness, max anomaly
MaxAnomalyModelRich <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
                                   fitFamily = "poisson",
                                   fixedFactors = c("UI2", "Tropical"),
                                   fixedTerms = list(StdTmaxAnomalyRS=1),
                                   randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                   fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1)",
                                                         "UI2:Tropical",
                                                         "Tropical:poly(StdTmaxAnomalyRS,1)",
                                                         "Tropical:poly(StdTmaxAnomalyRS,1):UI2"),  
                                   saveVars = c("Total_abundance", "SSBS", "NH_3000"))

summary(MaxAnomalyModelRich$model)

# selected model:
# Species_richness ~ UI2 + Tropical + poly(StdTmaxAnomalyRS, 1) + 
# UI2:poly(StdTmaxAnomalyRS, 1) + UI2:Tropical + Tropical:poly(StdTmaxAnomalyRS, 1):UI2 +  
# (1 | SS) + (1 | SSB) + (1 | SSBS)


# Model failed to converge with max|grad| = 0.0118007 (tol = 0.001, component 1)


# save model output
save(MaxAnomalyModelRich, file = paste0(outDir, "/MaxAnomalyModelRich_3way.rdata"))


#load(paste0(outDir, "/MeanAnomalyModelAbund.rdata"))
#load(paste0(outDir, "/MeanAnomalyModelRich.rdata"))
#load(paste0(outDir, "/MaxAnomalyModelAbund.rdata"))
#load(paste0(outDir, "/MaxAnomalyModelRich.rdata"))


summary(MeanAnomalyModelAbund$model)
summary(MeanAnomalyModelRich$model)
summary(MaxAnomalyModelAbund$model)
summary(MaxAnomalyModelRich$model)

############## plotting ##############


### 1. Abundance, mean anomaly

nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelAbund$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelAbund$data$UI2)),
  Tropical = factor(c("Tropical", "Temperate"),
                    levels = levels(MeanAnomalyModelAbund$data$Tropical)))

# back transform the predictors
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd$LogAbund <- 0
nd$Species_richness <- 0


# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))) & (nd$Tropical == "Temperate"))

# adjust plot 1: mean anomaly and abundance

exclQuantiles <- c(0.025,0.975)


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


nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

p1 <- ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  facet_wrap(~Tropical, ncol = 2) + 
  theme_bw() + 
  labs(fill = "% NH", col = "% NH") + 
  ylab("Total Abundance (%)") +
  xlab("Standardised Climate Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 500)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12))





### 2. SR mean anomaly

nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
                        to = max(MeanAnomalyModelRich$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MeanAnomalyModelRich$data$UI2)),
  Tropical = factor(c("Tropical", "Temperate"),
                    levels = levels(MeanAnomalyModelRich$data$Tropical)))

# back transform the predictors
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# set richness and abundance to 0 - to be predicted
nd2$LogAbund <- 0
nd2$Species_richness <- 0


# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) & (nd2$Tropical == "Temperate"))

# quantiles of data to show on the plot
exclQuantiles <- c(0.025,0.975)


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
sr.preds.tmean <- PredictGLMERRandIter(model = MeanAnomalyModelRich$model,data = nd2)

# back transform the abundance values
sr.preds.tmean <- exp(sr.preds.tmean)

# convert to relative to reference
sr.preds.tmean <- sweep(x = sr.preds.tmean,MARGIN = 2,STATS = sr.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
sr.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS < QPV[1])),] <- NA
sr.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS > QPV[2])),] <- NA
sr.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS < QSV[1])),] <- NA
sr.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS > QSV[2])),] <- NA
sr.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS < QAL[1])),] <- NA
sr.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS > QAL[2])),] <- NA
sr.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS < QAH[1])),] <- NA
sr.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd2$PredMedian <- ((apply(X = sr.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = sr.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = sr.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

p2 <- ggplot(data = nd2, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd2$PredLower, ymax = nd2$PredUpper, fill = UI2), alpha = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  facet_wrap(~Tropical, ncol = 2) + 
  theme_bw() + 
  labs(fill = "% NH", col = "% NH") + 
  ylab("Species Richness (%)") +
  xlab("Standardised Climate Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 500)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12))






### 3. Abundance, max anomaly


nd3 <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
                        to = max(MaxAnomalyModelAbund$data$StdTmaxAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxAnomalyModelAbund$data$UI2)),
  Tropical = factor(c("Tropical", "Temperate"),
                    levels = levels(MaxAnomalyModelAbund$data$Tropical)))

# back transform the predictors
nd3$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd3$LogAbund <- 0
nd3$Species_richness <- 0


# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomaly==min(abs(nd3$StdTmaxAnomaly))) & (nd3$Tropical == "Temperate"))

# quantiles of data to show on the plot
exclQuantiles <- c(0.025,0.975)


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

# predict the results
a.preds.tmean <- PredictGLMERRandIter(model = MaxAnomalyModelAbund$model,data = nd3)

# back transform the abundance values
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to relative to reference
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmaxAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmaxAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmaxAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmaxAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmaxAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmaxAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd3$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd3$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


nd3$UI2 <- factor(nd3$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

p3 <- ggplot(data = nd3, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd3$PredLower, ymax = nd3$PredUpper, fill = UI2), alpha = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  facet_wrap(~Tropical, ncol = 2) + 
  theme_bw() + 
  labs(fill = "% NH", col = "% NH") + 
  ylab("Total Abundance (%)") +
  xlab("Standardised Climate Anomaly Maximum") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 150)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12))





### 4. Species richness, max anomaly


nd4 <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(MaxAnomalyModelRich$data$StdTmaxAnomalyRS),
                       to = max(MaxAnomalyModelRich$data$StdTmaxAnomalyRS),
                       length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(MaxAnomalyModelRich$data$UI2)),
  Tropical = factor(c("Tropical", "Temperate"),
                    levels = levels(MaxAnomalyModelRich$data$Tropical)))

# back transform the predictors
nd4$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd4$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)

# set richness and abundance to 0 - to be predicted
nd4$LogAbund <- 0
nd4$Species_richness <- 0


# reference for % difference = primary vegetation and positive anomaly closest to 0
refRow <- which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomaly==min(abs(nd4$StdTmaxAnomaly))) & (nd4$Tropical == "Temperate"))

# quantiles of data to show on the plot
exclQuantiles <- c(0.025,0.975)


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

# predict the results
sr.preds.tmean <- PredictGLMERRandIter(model = MaxAnomalyModelRich$model,data = nd4)

# back transform the abundance values
sr.preds.tmean <- exp(sr.preds.tmean)

# convert to relative to reference
sr.preds.tmean <- sweep(x = sr.preds.tmean,MARGIN = 2,STATS = sr.preds.tmean[refRow,],FUN = '/')

# remove anything above and below the quantiles
sr.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomalyRS < QPV[1])),] <- NA
sr.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmaxAnomalyRS > QPV[2])),] <- NA
sr.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmaxAnomalyRS < QSV[1])),] <- NA
sr.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmaxAnomalyRS > QSV[2])),] <- NA
sr.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmaxAnomalyRS < QAL[1])),] <- NA
sr.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmaxAnomalyRS > QAL[2])),] <- NA
sr.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmaxAnomalyRS < QAH[1])),] <- NA
sr.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmaxAnomalyRS > QAH[2])),] <- NA

# Get the median, upper and lower quants for the plot
nd4$PredMedian <- ((apply(X = sr.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd4$PredUpper <- ((apply(X = sr.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4$PredLower <- ((apply(X = sr.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100


nd4$UI2 <- factor(nd4$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High"))

p4 <- ggplot(data = nd4, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = UI2), size = 1) +
  geom_ribbon(aes(ymin = nd4$PredLower, ymax = nd4$PredUpper, fill = UI2), alpha = 0.2) +
  scale_fill_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  scale_colour_manual(values = c("#009E73", "#0072B2","#E69F00", "#D55E00")) +
  facet_wrap(~Tropical, ncol = 2) + 
  theme_bw() + 
  labs(fill = "% NH", col = "% NH") + 
  ylab("Species Richness (%)") +
  xlab("Standardised Climate Anomaly Maximum") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 150)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12))

#"#009E73" - green
#"#0072B2" - blue
#"#E69F00" - yellow
#"#D55E00" - red



# organise plots into one document

plot_grid(p1, p2, p3, p4, ncol = 1)

# save plot
ggsave(filename = paste0(outDir, "Plots_climate_LU_Tropical_ALL_allinteractions.pdf"), width = 8, height = 12, units = "in")


