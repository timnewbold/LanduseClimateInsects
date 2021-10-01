##%######################################################%##
#                                                          #
####             Testing buffer zone sizes              ####
####             for natural habitat metric             ####
#                                                          #
##%######################################################%##

# testing the difference between model results when different 
# buffer sizes are used


rm(list = ls())

library(StatisticalModels)
library(ggplot2)
library(cowplot)
source('Functions.R')


# directories
predictsDataDir <- "6_RunLUClimateModels/"
outDir <- "14_Additional_Tests/"


###Create Models for all insects in predicts for Standardised climate anomaly and Land interactions

predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSiteData.rds"))


#### Abundance models ####

modelData <- predictsSites[!is.na(predictsSites$LogAbund), ] # 5735

AbundMeanAnomalyModel_1000 <- GLMER(modelData = modelData,responseVar = "LogAbund",fitFamily = "gaussian",
                                fixedStruct = "UI2 * StdTmeanAnomalyRS * NH_1000.rs",
                                randomStruct = "(1|SS)+(1|SSB)",
                                saveVars = c("SSBS"))
AbundMeanAnomalyModel_3000 <- GLMER(modelData = modelData,responseVar = "LogAbund",fitFamily = "gaussian",
                                fixedStruct = "UI2 * StdTmeanAnomalyRS * NH_3000.rs",
                                randomStruct = "(1|SS)+(1|SSB)",
                                saveVars = c("SSBS"))
AbundMeanAnomalyModel_5000 <- GLMER(modelData = modelData,responseVar = "LogAbund",fitFamily = "gaussian",
                                fixedStruct = "UI2 * StdTmeanAnomalyRS * NH_5000.rs",
                                randomStruct = "(1|SS)+(1|SSB)",
                                saveVars = c("SSBS"))
AbundMeanAnomalyModel_10000 <- GLMER(modelData = modelData,responseVar = "LogAbund",fitFamily = "gaussian",
                                fixedStruct = "UI2 * StdTmeanAnomalyRS * NH_10000.rs",
                                randomStruct = "(1|SS)+(1|SSB)",
                                saveVars = c("SSBS"))


### models for the stats


MeanAbun_1000 <- GLMERSelect(modelData = modelData,
                        responseVar = "LogAbund",
                        fitFamily = "gaussian",fixedFactors = "UI2",
                        fixedTerms = list(StdTmeanAnomalyRS=1, NH_1000.rs=1),
                        randomStruct = "(1|SS)+(1|SSB)",
                        fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1):poly(NH_1000.rs,1)"),
                        fitInteractions = TRUE,
                        saveVars = c("SSBS"))


test <- GLMERSelect(modelData = modelData,
            responseVar = "LogAbund",
            fitFamily = "gaussian",fixedFactors = "UI2",
            fixedTerms = list(StdTmeanAnomalyRS=1, NH_1000.rs=1),
            randomStruct = "(1|SS)+(1|SSB)",
            fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1):poly(NH_1000.rs,1)"))

test$model


save(MeanAbun_1000, file = paste0(outDir, "MeanAbun_buffer_1000.rdata"))

MeanAbun_3000 <- GLMERSelect(modelData = modelData,
                        responseVar = "LogAbund",
                        fitFamily = "gaussian",fixedFactors = "UI2",
                        fixedTerms = list(StdTmeanAnomalyRS=1, NH_3000.rs=1),
                        randomStruct = "(1|SS)+(1|SSB)",
                        fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1):poly(NH_3000.rs,1)"),
                        fitInteractions = TRUE,
                        saveVars = c("SSBS"))

save(MeanAbun_3000, file = paste0(outDir, "MeanAbun_buffer_3000.rdata"))

MeanAbun_5000 <- GLMERSelect(modelData = modelData,
                        responseVar = "LogAbund",
                        fitFamily = "gaussian",fixedFactors = "UI2",
                        fixedTerms = list(StdTmeanAnomalyRS=1, NH_5000.rs=1),
                        randomStruct = "(1|SS)+(1|SSB)",
                        fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1):poly(NH_5000.rs,1)"),
                        fitInteractions = TRUE,
                        saveVars = c("SSBS"))

save(MeanAbun_5000, file = paste0(outDir, "MeanAbun_buffer_5000.rdata"))


MeanAbun_10000 <- GLMERSelect(modelData = modelData,
                        responseVar = "LogAbund",
                        fitFamily = "gaussian",fixedFactors = "UI2",
                        fixedTerms = list(StdTmeanAnomalyRS=1, NH_10000.rs=1),
                        randomStruct = "(1|SS)+(1|SSB)",
                        fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1):poly(NH_10000.rs,1)"),
                        fitInteractions = TRUE,
                        saveVars = c("SSBS"))

save(MeanAbun_10000, file = paste0(outDir, "MeanAbun_buffer_10000.rdata"))


# take a look at the results

MeanAbun_1000$stats
MeanAbun_3000$stats
MeanAbun_5000$stats
MeanAbun_10000$stats

#### Species richness models ####

modelData <- predictsSites

RichMeanAnomalyModel_1000 <- GLMER(modelData = modelData,responseVar = "Species_richness",fitFamily = "poisson",
                               fixedStruct = "UI2 * StdTmeanAnomalyRS * NH_1000.rs",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                               saveVars = c("SSBS"))

RichMeanAnomalyModel_3000 <- GLMER(modelData = modelData,responseVar = "Species_richness",fitFamily = "poisson",
                               fixedStruct = "UI2 * StdTmeanAnomalyRS * NH_3000.rs",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                               saveVars = c("SSBS"))

RichMeanAnomalyModel_5000 <- GLMER(modelData = modelData,responseVar = "Species_richness",fitFamily = "poisson",
                               fixedStruct = "UI2 * StdTmeanAnomalyRS * NH_5000.rs",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                               saveVars = c("SSBS"))

RichMeanAnomalyModel_10000 <- GLMER(modelData = modelData,responseVar = "Species_richness",fitFamily = "poisson",
                               fixedStruct = "UI2 * StdTmeanAnomalyRS * NH_10000.rs",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                               saveVars = c("SSBS"))


# take a look at the results
fixef(RichMeanAnomalyModel_1000$model)
fixef(RichdMeanAnomalyModel_3000$model)
fixef(RichMeanAnomalyModel_5000$model)
fixef(RichMeanAnomalyModel_10000$model)




#### plots of abundance results ####

### 1000 buffer

AbundMeanAnomalyModel1 <- AbundMeanAnomalyModel_1000

nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(AbundMeanAnomalyModel1$data$UI2)),
  NH_1000.rs=c(-0.96,0.01493712,1.00,2.0))
# NH_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform NH data range
nd$NH_1000 <- round(BackTransformCentreredPredictor(
  transformedX = nd$NH_1000.rs,originalX = predictsSites$NH_1000)*100,0)

nd$NH_1000[nd$NH_1000 == 99] <- 100

# set values for richness and Abundance
nd$LogAbund <- 0
nd$Species_richness <- 0

# set the reference row
refRow <- which((nd$UI2=="Primary vegetation") & 
                  nd$StdTmeanAnomaly== min(nd$StdTmeanAnomaly[nd$StdTmeanAnomaly > 0]) &
                  #(nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) & 
                  (nd$NH_1000==100))

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
a.preds.tmean <- PredictGLMERRandIter(model = AbundMeanAnomalyModel1$model,data = nd, nIters = 10000)

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
# Abundance response to mean anomaly only, all LUs

# set factor levels
nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))

nd$NH_1000 <- factor(nd$NH_1000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd <- nd[nd$UI2 %in% c("Agriculture_Low", "Agriculture_High"), ]

# plot
p1 <- ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = NH_1000), size = 1) +
  geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = NH_1000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 2, labeller = as_labeller(c('Agriculture_Low' = "a              Agriculture_Low", 'Agriculture_High' = "b              Agriculture_High"))) + 
  theme_bw() + 
  labs(fill = "% NH", col = "% NH") + 
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 210)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"))



### 3000 buffer

AbundMeanAnomalyModel1 <- AbundMeanAnomalyModel_3000

nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(AbundMeanAnomalyModel1$data$UI2)),
  NH_3000.rs=c(-0.99,0.01493712,1.02,2.03))
# NH_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform NH data range
nd2$NH_3000 <- round(BackTransformCentreredPredictor(
  transformedX = nd2$NH_3000.rs,originalX = predictsSites$NH_3000)*100,0)

nd2$NH_3000[nd2$NH_3000 == 99] <- 100

# set values for richness and Abundance
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# set the reference row
refRow <- which((nd2$UI2=="Primary vegetation") & 
                  nd2$StdTmeanAnomaly== min(nd2$StdTmeanAnomaly[nd2$StdTmeanAnomaly > 0]) &
                  #(nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) & 
                  (nd2$NH_3000==100))

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
a.preds.tmean <- PredictGLMERRandIter(model = AbundMeanAnomalyModel1$model,data = nd2, nIters = 10000)

# transform results
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to percentage of reference row
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS > QAH[2])),] <- NA

# get the median and upper/lower intervals for plots
nd2$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# Abundance response to mean anomaly only, all LUs

# set factor levels
nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))

nd2$NH_3000 <- factor(nd2$NH_3000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd2 <- nd2[nd2$UI2 %in% c("Agriculture_Low", "Agriculture_High"), ]

# plot
p2 <- ggplot(data = nd2, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = NH_3000), size = 1) +
  geom_ribbon(aes(ymin = nd2$PredLower, ymax = nd2$PredUpper, fill = NH_3000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 2, labeller = as_labeller(c('Agriculture_Low' = "c              Agriculture_Low", 'Agriculture_High' = "d              Agriculture_High"))) + 
  theme_bw() + 
  labs(fill = "% NH", col = "% NH") + 
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 210)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"))





### 5000 buffer

AbundMeanAnomalyModel1 <- AbundMeanAnomalyModel_5000

nd3 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(AbundMeanAnomalyModel1$data$UI2)),
  NH_5000.rs=c(-1.05,0.01493712,1.045653,2.073849))
# NH_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd3$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform NH data range
nd3$NH_5000 <- round(BackTransformCentreredPredictor(
  transformedX = nd3$NH_5000.rs,originalX = predictsSites$NH_5000)*100,0)

nd3$NH_5000[nd3$NH_5000 == 99] <- 100

# set values for richness and Abundance
nd3$LogAbund <- 0
nd3$Species_richness <- 0

# set the reference row
refRow <- which((nd3$UI2=="Primary vegetation") & 
                  nd3$StdTmeanAnomaly== min(nd3$StdTmeanAnomaly[nd3$StdTmeanAnomaly > 0]) &
                  #(nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) & 
                  (nd3$NH_5000==100))

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
a.preds.tmean <- PredictGLMERRandIter(model = AbundMeanAnomalyModel1$model,data = nd3, nIters = 10000)

# transform results
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to percentage of reference row
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmeanAnomalyRS > QAH[2])),] <- NA

# get the median and upper/lower intervals for plots
nd3$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd3$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# Abundance response to mean anomaly only, all LUs

# set factor levels
nd3$UI2 <- factor(nd3$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))

nd3$NH_5000 <- factor(nd3$NH_5000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd3 <- nd3[nd3$UI2 %in% c("Agriculture_Low", "Agriculture_High"), ]

# plot
p3 <- ggplot(data = nd3, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = NH_5000), size = 1) +
  geom_ribbon(aes(ymin = nd3$PredLower, ymax = nd3$PredUpper, fill = NH_5000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 2, labeller = as_labeller(c('Agriculture_Low' = "e              Agriculture_Low", 'Agriculture_High' = "f              Agriculture_High"))) + 
  theme_bw() + 
  labs(fill = "% NH", col = "% NH") + 
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 210)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"))


### 10000 buffer


AbundMeanAnomalyModel1 <- AbundMeanAnomalyModel_10000

nd4 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(AbundMeanAnomalyModel1$data$UI2)),
  NH_10000.rs=c(-1.08,-0.01,1.045653,2.08))
# NH_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd4$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd4$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform NH data range
nd4$NH_10000 <- round(BackTransformCentreredPredictor(
  transformedX = nd4$NH_10000.rs,originalX = predictsSites$NH_10000)*100,0)

nd4$NH_10000[nd4$NH_10000 == 99] <- 100

# set values for richness and Abundance
nd4$LogAbund <- 0
nd4$Species_richness <- 0

# set the reference row
refRow <- which((nd4$UI2=="Primary vegetation") & 
                  nd4$StdTmeanAnomaly== min(nd4$StdTmeanAnomaly[nd4$StdTmeanAnomaly > 0]) &
                  #(nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) & 
                  (nd4$NH_10000==100))

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
a.preds.tmean <- PredictGLMERRandIter(model = AbundMeanAnomalyModel1$model,data = nd4, nIters = 10000)

# transform results
a.preds.tmean <- exp(a.preds.tmean)-0.01

# convert to percentage of reference row
a.preds.tmean <- sweep(x = a.preds.tmean,MARGIN = 2,STATS = a.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
a.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmeanAnomalyRS < QPV[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmeanAnomalyRS > QPV[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmeanAnomalyRS < QSV[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmeanAnomalyRS > QSV[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmeanAnomalyRS < QAL[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmeanAnomalyRS > QAL[2])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmeanAnomalyRS < QAH[1])),] <- NA
a.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmeanAnomalyRS > QAH[2])),] <- NA

# get the median and upper/lower intervals for plots
nd4$PredMedian <- ((apply(X = a.preds.tmean,MARGIN = 1,
                          FUN = median,na.rm=TRUE))*100)-100
nd4$PredUpper <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4$PredLower <- ((apply(X = a.preds.tmean,MARGIN = 1,
                         FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# Abundance response to mean anomaly only, all LUs

# set factor levels
nd4$UI2 <- factor(nd4$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))

nd4$NH_10000 <- factor(nd4$NH_10000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd4 <- nd4[nd4$UI2 %in% c("Agriculture_Low", "Agriculture_High"), ]

# plot
p4 <- ggplot(data = nd4, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = NH_10000), size = 1) +
  geom_ribbon(aes(ymin = nd4$PredLower, ymax = nd4$PredUpper, fill = NH_10000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 2, labeller = as_labeller(c('Agriculture_Low' = "g              Agriculture_Low", 'Agriculture_High' = "h              Agriculture_High"))) + 
  theme_bw() + 
  labs(fill = "% NH", col = "% NH") + 
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 210)) + 
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"))



cowplot::plot_grid(p1, p2, p3, p4, nrow = 4, labels = c("1000m", "3000m", "5000m", "10000m"))

# save
ggsave(filename = paste0(outDir, "Abun_buffers_plots.pdf"), height = 16, width = 9)



#### Plots for Richness results ####

### 1000m buffer

RichMeanAnomalyModel1 <- RichMeanAnomalyModel_1000

### mean anomaly species richness 
nd <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(RichMeanAnomalyModel1$data$UI2)),
  NH_1000.rs=c(-0.96,0.01493712,1.00,2.0))
#NH_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform NH data range
nd$NH_1000 <- round(BackTransformCentreredPredictor(
  transformedX = nd$NH_1000.rs,originalX = predictsSites$NH_1000)*100,0)

nd$NH_1000[nd$NH_1000 == 99] <- 100

# set values for richness and Abundance
nd$LogAbund <- 0
nd$Species_richness <- 0

# set the reference row
refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))) &
                  (nd$NH_1000==100))

# quantiles for presenting results
exclQuantiles <- c(0.025,0.975)


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

# predict results
s.preds.tmean <- PredictGLMERRandIter(model = RichMeanAnomalyModel1$model,data = nd)

# transform results
s.preds.tmean <- exp(s.preds.tmean)

# convert to percentage of reference row
s.preds.tmean <- sweep(x = s.preds.tmean,MARGIN = 2,STATS = s.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
s.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS < QPV[1])),] <- NA
s.preds.tmean[which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomalyRS > QPV[2])),] <- NA
s.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS < QSV[1])),] <- NA
s.preds.tmean[which((nd$UI2=="Secondary vegetation") & (nd$StdTmeanAnomalyRS > QSV[2])),] <- NA
s.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS < QAL[1])),] <- NA
s.preds.tmean[which((nd$UI2=="Agriculture_Low") & (nd$StdTmeanAnomalyRS > QAL[2])),] <- NA
s.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS < QAH[1])),] <- NA
s.preds.tmean[which((nd$UI2=="Agriculture_High") & (nd$StdTmeanAnomalyRS > QAH[2])),] <- NA

# get the median and upper/lower intervals for plots
nd$PredMedian <- ((apply(X = s.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd$PredUpper <- ((apply(X = s.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd$PredLower <- ((apply(X = s.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# Abundance response to mean anomaly only, all LUs

# set factor levels
nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))

nd$NH_1000 <- factor(nd$NH_1000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd <- nd[nd$UI2 %in% c("Agriculture_Low", "Agriculture_High"), ]

# plot
p1 <- ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) +
  geom_line(aes(col = NH_1000), size = 1) +
  geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = NH_1000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 2, labeller = as_labeller(c('Agriculture_Low' = "a              Agriculture_Low", 'Agriculture_High' = "b              Agriculture_High"))) +
  theme_bw() +
  labs(fill = "% NH", col = "% NH") +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 150)) +
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"))



### 3000m buffer

RichMeanAnomalyModel1 <- RichMeanAnomalyModel_3000

### mean anomaly species richness 
nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(RichMeanAnomalyModel1$data$UI2)),
  NH_3000.rs=c(-0.99,0.01493712,1.02,2.03))
#NH_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform NH data range
nd2$NH_3000 <- round(BackTransformCentreredPredictor(
  transformedX = nd2$NH_3000.rs,originalX = predictsSites$NH_3000)*100,0)

nd2$NH_3000[nd2$NH_3000 == 99] <- 100

# set values for richness and Abundance
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# set the reference row
refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) &
                  (nd2$NH_3000==100))

# quantiles for presenting results
exclQuantiles <- c(0.025,0.975)


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

# predict results
s.preds.tmean <- PredictGLMERRandIter(model = RichMeanAnomalyModel1$model,data = nd2)

# transform results
s.preds.tmean <- exp(s.preds.tmean)

# convert to percentage of reference row
s.preds.tmean <- sweep(x = s.preds.tmean,MARGIN = 2,STATS = s.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
s.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS < QPV[1])),] <- NA
s.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS > QPV[2])),] <- NA
s.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS < QSV[1])),] <- NA
s.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS > QSV[2])),] <- NA
s.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS < QAL[1])),] <- NA
s.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS > QAL[2])),] <- NA
s.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS < QAH[1])),] <- NA
s.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS > QAH[2])),] <- NA

# get the median and upper/lower intervals for plots
nd2$PredMedian <- ((apply(X = s.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd2$PredUpper <- ((apply(X = s.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd2$PredLower <- ((apply(X = s.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# Abundance response to mean anomaly only, all LUs

# set factor levels
nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))

nd2$NH_3000 <- factor(nd2$NH_3000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd2 <- nd2[nd2$UI2 %in% c("Agriculture_Low", "Agriculture_High"), ]

# plot
p2 <- ggplot(data = nd2, aes(x = StdTmeanAnomaly, y = PredMedian)) +
  geom_line(aes(col = NH_3000), size = 1) +
  geom_ribbon(aes(ymin = nd2$PredLower, ymax = nd2$PredUpper, fill = NH_3000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 2, labeller = as_labeller(c('Agriculture_Low' = "c              Agriculture_Low", 'Agriculture_High' = "d              Agriculture_High"))) +
  theme_bw() +
  labs(fill = "% NH", col = "% NH") +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 150)) +
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"))




### 5000m buffer

RichMeanAnomalyModel1 <- RichMeanAnomalyModel_5000

### mean anomaly species richness 
nd3 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(RichMeanAnomalyModel1$data$UI2)),
  NH_5000.rs=c(-1.05,0.01493712,1.045653,2.073849))
#NH_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd3$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform NH data range
nd3$NH_5000 <- round(BackTransformCentreredPredictor(
  transformedX = nd3$NH_5000.rs,originalX = predictsSites$NH_5000)*100,0)

nd3$NH_5000[nd3$NH_5000 == 99] <- 100

# set values for richness and Abundance
nd3$LogAbund <- 0
nd3$Species_richness <- 0

# set the reference row
refRow <- which((nd3$UI2=="Primary vegetation") & (nd3$StdTmeanAnomaly==min(abs(nd3$StdTmeanAnomaly))) &
                  (nd3$NH_5000==100))

# quantiles for presenting results
exclQuantiles <- c(0.025,0.975)


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

# predict results
s.preds.tmean <- PredictGLMERRandIter(model = RichMeanAnomalyModel1$model,data = nd3)

# transform results
s.preds.tmean <- exp(s.preds.tmean)

# convert to percentage of reference row
s.preds.tmean <- sweep(x = s.preds.tmean,MARGIN = 2,STATS = s.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
s.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmeanAnomalyRS < QPV[1])),] <- NA
s.preds.tmean[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmeanAnomalyRS > QPV[2])),] <- NA
s.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmeanAnomalyRS < QSV[1])),] <- NA
s.preds.tmean[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmeanAnomalyRS > QSV[2])),] <- NA
s.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmeanAnomalyRS < QAL[1])),] <- NA
s.preds.tmean[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmeanAnomalyRS > QAL[2])),] <- NA
s.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmeanAnomalyRS < QAH[1])),] <- NA
s.preds.tmean[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmeanAnomalyRS > QAH[2])),] <- NA

# get the median and upper/lower intervals for plots
nd3$PredMedian <- ((apply(X = s.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd3$PredUpper <- ((apply(X = s.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3$PredLower <- ((apply(X = s.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# Abundance response to mean anomaly only, all LUs

# set factor levels
nd3$UI2 <- factor(nd3$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))

nd3$NH_5000 <- factor(nd3$NH_5000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd3 <- nd3[nd3$UI2 %in% c("Agriculture_Low", "Agriculture_High"), ]

# plot
p3 <- ggplot(data = nd3, aes(x = StdTmeanAnomaly, y = PredMedian)) +
  geom_line(aes(col = NH_5000), size = 1) +
  geom_ribbon(aes(ymin = nd3$PredLower, ymax = nd3$PredUpper, fill = NH_5000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 2, labeller = as_labeller(c('Agriculture_Low' = "e              Agriculture_Low", 'Agriculture_High' = "f              Agriculture_High"))) +
  theme_bw() +
  labs(fill = "% NH", col = "% NH") +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 150)) +
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"))


### 10000m buffer


RichMeanAnomalyModel1 <- RichMeanAnomalyModel_10000

### mean anomaly species richness 
nd4 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(RichMeanAnomalyModel1$data$UI2)),
  NH_10000.rs=c(-1.08,-0.01,1.045653,2.08))
#NH_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd4$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd4$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform NH data range
nd4$NH_10000 <- round(BackTransformCentreredPredictor(
  transformedX = nd4$NH_10000.rs,originalX = predictsSites$NH_10000)*100,0)

nd4$NH_10000[nd4$NH_10000 == 99] <- 100

# set values for richness and Abundance
nd4$LogAbund <- 0
nd4$Species_richness <- 0

# set the reference row
refRow <- which((nd4$UI2=="Primary vegetation") & (nd4$StdTmeanAnomaly==min(abs(nd4$StdTmeanAnomaly))) &
                  (nd4$NH_10000==100))

# quantiles for presenting results
exclQuantiles <- c(0.025,0.975)


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

# predict results
s.preds.tmean <- PredictGLMERRandIter(model = RichMeanAnomalyModel1$model,data = nd4)

# transform results
s.preds.tmean <- exp(s.preds.tmean)

# convert to percentage of reference row
s.preds.tmean <- sweep(x = s.preds.tmean,MARGIN = 2,STATS = s.preds.tmean[refRow,],FUN = '/')

# set anything outside the desired quantiles to NA
s.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmeanAnomalyRS < QPV[1])),] <- NA
s.preds.tmean[which((nd4$UI2=="Primary vegetation") & (nd4$StdTmeanAnomalyRS > QPV[2])),] <- NA
s.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmeanAnomalyRS < QSV[1])),] <- NA
s.preds.tmean[which((nd4$UI2=="Secondary vegetation") & (nd4$StdTmeanAnomalyRS > QSV[2])),] <- NA
s.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmeanAnomalyRS < QAL[1])),] <- NA
s.preds.tmean[which((nd4$UI2=="Agriculture_Low") & (nd4$StdTmeanAnomalyRS > QAL[2])),] <- NA
s.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmeanAnomalyRS < QAH[1])),] <- NA
s.preds.tmean[which((nd4$UI2=="Agriculture_High") & (nd4$StdTmeanAnomalyRS > QAH[2])),] <- NA

# get the median and upper/lower intervals for plots
nd4$PredMedian <- ((apply(X = s.preds.tmean,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd4$PredUpper <- ((apply(X = s.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd4$PredLower <- ((apply(X = s.preds.tmean,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# Abundance response to mean anomaly only, all LUs

# set factor levels
nd4$UI2 <- factor(nd4$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))

nd4$NH_10000 <- factor(nd4$NH_10000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd4 <- nd4[nd4$UI2 %in% c("Agriculture_Low", "Agriculture_High"), ]

# plot
p4 <- ggplot(data = nd4, aes(x = StdTmeanAnomaly, y = PredMedian)) +
  geom_line(aes(col = NH_10000), size = 1) +
  geom_ribbon(aes(ymin = nd4$PredLower, ymax = nd4$PredUpper, fill = NH_10000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 2, labeller = as_labeller(c('Agriculture_Low' = "g              Agriculture_Low", 'Agriculture_High' = "h              Agriculture_High"))) +
  theme_bw() +
  labs(fill = "% NH", col = "% NH") +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 150)) +
  theme(aspect.ratio = 1, text = element_text(size = 12),
        strip.text.x = element_text(hjust = 0, size = 12, face = "bold"))




cowplot::plot_grid(p1, p2, p3, p4, nrow = 4, labels = c("1000m", "3000m", "5000m", "10000m"))

# save
ggsave(filename = paste0(outDir, "Rich_buffers_plots.pdf"), height = 16, width = 9)


