##%######################################################%##
#                                                          #
####         Run models including percentage NH         ####
#                                                          #
##%######################################################%##

# This script runs the more complex models looking at the buffering effect
# of natural habitat on climate effects across land uses.


# directories
predictsDataDir <- "6_RunLUClimateModels/"
outDir <- "7_RunLUClimateNHModels/"
if(!dir.exists(outDir)) dir.create(outDir)


sink(paste0(outDir,"log.txt"))

t.start <- Sys.time()

print(t.start)



# load libraries
library(StatisticalModels)
library(predictsFunctions)
source("Functions.R")
library(ggplot2)
library(cowplot)
library(sjPlot)


# read in data
predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSiteData.rds"))

### explore the data 

# look at spread of natural habitat availability across land used
predictsSites$UI2 <- factor(predictsSites$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))

ggplot(data = predictsSites) + 
  geom_boxplot(aes(x = UI2, y = NH_5000)) + 
  theme_bw()

ggsave(filename = paste0(outDir, "Boxplot_NH_LU.pdf"))

# look at spread of NH across anomlay values coloured by LU
ggplot(data = predictsSites)+
  geom_point(aes(x = StdTmeanAnomaly, y = NH_5000, col = UI2)) +
  scale_colour_manual(values = c("#009E73", "#0072B2" , "#E69F00", "#D55E00")) + 
  theme_bw()

ggsave(filename = paste0(outDir, "Scatter_NH_LU_Anom.pdf"))


#### Run the models ####

# remove any NAs
modelData <- na.omit(predictsSites[,c(
  'LogAbund','UI2','StdTmeanAnomalyRS','SS','SSB','SSBS','NH_5000.rs')]) # 5735 rows

# run models with and without NH interaction, Abundance models
AbundMeanAnomalyModel0 <- GLMER(modelData = modelData,responseVar = "LogAbund",fitFamily = "gaussian",
                fixedStruct = "UI2 * StdTmeanAnomalyRS",
                randomStruct = "(1|SS)+(1|SSB)")
AbundMeanAnomalyModel1 <- GLMER(modelData = modelData,responseVar = "LogAbund",fitFamily = "gaussian",
                fixedStruct = "UI2 * StdTmeanAnomalyRS * NH_5000.rs",
                randomStruct = "(1|SS)+(1|SSB)",
                saveVars = c("SSBS"))

print(anova(AbundMeanAnomalyModel0$model,AbundMeanAnomalyModel1$model))

# save for use in predictions script
save(AbundMeanAnomalyModel1, file = paste0(outDir, "/MeanAnomalyModelAbun_NH.rdata"))
#load( file = paste0(outDir, "/MeanAnomalyModelAbun_NH.rdata"))


# run models with and without NH interaction, species richness models
modelData <- na.omit(predictsSites[,c(
  'Species_richness','UI2','StdTmeanAnomalyRS','SS','SSB','SSBS','NH_5000.rs')]) # 6069 rows

RichMeanAnomalyModel0 <- GLMER(modelData = modelData,responseVar = "Species_richness",fitFamily = "poisson",
                               fixedStruct = "UI2 * StdTmeanAnomalyRS",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)")
RichMeanAnomalyModel1 <- GLMER(modelData = modelData,responseVar = "Species_richness",fitFamily = "poisson",
                               fixedStruct = "UI2 * StdTmeanAnomalyRS * NH_5000.rs",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                               saveVars = c("SSBS"))

print(anova(RichMeanAnomalyModel0$model,RichMeanAnomalyModel1$model))


save(RichMeanAnomalyModel1, file = paste0(outDir, "/RichMeanAnomalyModel_NH.rdata"))
#load(file = paste0(outDir, "/RichMeanAnomalyModel_NH.rdata"))


# run models with and without NH interaction, Abundance models, max anomaly
modelData <- na.omit(predictsSites[,c(
  'LogAbund','UI2','StdTmaxAnomalyRS','SS','SSB','SSBS','NH_5000.rs')]) # 5735 rows

AbundMaxAnomalyModel0 <- GLMER(modelData = modelData,responseVar = "LogAbund",fitFamily = "gaussian",
                                fixedStruct = "UI2 * StdTmaxAnomalyRS",
                                randomStruct = "(1|SS)+(1|SSB)")
AbundMaxAnomalyModel1 <- GLMER(modelData = modelData,responseVar = "LogAbund",fitFamily = "gaussian",
                                fixedStruct = "UI2 * StdTmaxAnomalyRS * NH_5000.rs",
                                randomStruct = "(1|SS)+(1|SSB)",
                               saveVars = c("SSBS"))

print(anova(AbundMaxAnomalyModel0$model,AbundMaxAnomalyModel1$model))


save(AbundMaxAnomalyModel1, file = paste0(outDir, "/AbundMaxAnomalyModel_NH.rdata"))
#load(file = paste0(outDir, "/AbundMaxAnomalyModel_NH.rdata"))



# run models with and without NH interaction, species richness models, max anomaly
modelData <- na.omit(predictsSites[,c(
  'Species_richness','UI2','StdTmaxAnomalyRS','SS','SSB','SSBS','NH_5000.rs')]) # 6069 rows

RichMaxAnomalyModel0 <- GLMER(modelData = modelData,responseVar = "Species_richness",fitFamily = "poisson",
                               fixedStruct = "UI2 * StdTmaxAnomalyRS",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)")
RichMaxAnomalyModel1 <- GLMER(modelData = modelData,responseVar = "Species_richness",fitFamily = "poisson",
                               fixedStruct = "UI2 * StdTmaxAnomalyRS * NH_5000.rs",
                               randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                              saveVars = c("SSBS"))

print(anova(RichMaxAnomalyModel0$model,RichMaxAnomalyModel1$model))

save(RichMaxAnomalyModel1, file = paste0(outDir, "/RichMaxAnomalyModel_NH.rdata"))
#load(file = paste0(outDir, "/RichMaxAnomalyModel_NH.rdata"))


##%######################################################%##
#                                                          #
####                    Model stats                     ####
#                                                          #
##%######################################################%##

# rerun models using GLMERSelect function to extract stats easily.

# Mean Anom

modelData <- na.omit(predictsSites[,c(
  'LogAbund','UI2','StdTmeanAnomalyRS','SS','SSB','SSBS','NH_5000.rs')]) # 5735 rows

# 1. Abundance, mean anomaly
MeanAbun <- GLMERSelect(modelData = modelData,
                        responseVar = "LogAbund",
                        fitFamily = "gaussian",fixedFactors = "UI2",
                        fixedTerms = list(StdTmeanAnomalyRS=1, NH_5000.rs=1),
                        randomStruct = "(1|SS)+(1|SSB)",
                        fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1):poly(NH_5000.rs,1)"),
                        fitInteractions = TRUE,
                        saveVars = c("SSBS"))


MeanAbun_stats <- as.data.frame(MeanAbun$stats)
MeanAbun_stats$significant <- NA

# function to check significance
checksig <- function(x){
  if(x <= 0.05){ 
    res <- "Yes" 
  } else { 
    res <- "No" }
  return(res)}

MeanAbun_stats$significant <- sapply(X = MeanAbun_stats$P, FUN = checksig)
write.csv(MeanAbun_stats, file = paste0(outDir, "/MeanAnomAbun_Stats.csv"), row.names = FALSE)


# 2. SR, mean anomaly

modelData <- na.omit(predictsSites[,c(
  'Species_richness','UI2','StdTmeanAnomalyRS','SS','SSB','SSBS','NH_5000.rs')]) # 6069 rows

MeanRich <- GLMERSelect(modelData = modelData,
                        responseVar = "Species_richness",
                        fitFamily = "poisson",fixedFactors = "UI2",
                        fixedTerms = list(StdTmeanAnomalyRS=1, NH_5000.rs=1),
                        randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                        fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1):poly(NH_5000.rs,1)"),
                        fitInteractions = TRUE)


MeanRich_stats <- as.data.frame(MeanRich$stats)
MeanRich_stats$significant <- NA
MeanRich_stats$significant <- sapply(X = MeanRich_stats$P, FUN = checksig)
write.csv(MeanRich_stats, file = paste0(outDir, "/MeanAnomRich_Stats.csv"), row.names = FALSE)

summary(RichMeanAnomalyModel1$model)
summary(MeanRich$model)


# Max Anom
modelData <- na.omit(predictsSites[,c(
  'LogAbund','UI2','StdTmaxAnomalyRS','SS','SSB','SSBS','NH_5000.rs')]) # 5735 rows

# 3. Abundance, max anomaly
MaxAbun <- GLMERSelect(modelData = modelData,
                        responseVar = "LogAbund",
                        fitFamily = "gaussian",fixedFactors = "UI2",
                        fixedTerms = list(StdTmaxAnomalyRS=1, NH_5000.rs=1),
                        randomStruct = "(1|SS)+(1|SSB)",
                        fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1):poly(NH_5000.rs,1)"),
                        fitInteractions = TRUE,
                        saveVars = c("SSBS"))

MaxAbun_stats <- as.data.frame(MaxAbun$stats)
MaxAbun_stats$significant <- NA
MaxAbun_stats$significant <- sapply(X = MaxAbun_stats$P, FUN = checksig)
write.csv(MaxAbun_stats, file = paste0(outDir, "/MaxAnomAbun_Stats.csv"), row.names = FALSE)

summary(AbundMaxAnomalyModel1$model)
summary(MaxAbun$model)


# 4. SR, max anomaly
modelData <- na.omit(predictsSites[,c(
  'Species_richness','UI2','StdTmaxAnomalyRS','SS','SSB','SSBS','NH_5000.rs')]) # 6069 rows


MaxRich <- GLMERSelect(modelData = modelData,
                        responseVar = "Species_richness",
                        fitFamily = "poisson",fixedFactors = "UI2",
                        fixedTerms = list(StdTmaxAnomalyRS=1, NH_5000.rs=1),
                        randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                        fixedInteractions = c("UI2:poly(StdTmaxAnomalyRS,1):poly(NH_5000.rs,1)"),
                        fitInteractions = TRUE)

MaxRich_stats <- as.data.frame(MaxRich$stats)
MaxRich_stats$significant <- NA
MaxRich_stats$significant <- sapply(X = MaxRich_stats$P, FUN = checksig)
write.csv(MaxRich_stats, file = paste0(outDir, "/MaxAnomRich_Stats.csv"), row.names = FALSE)

summary(RichMaxAnomalyModel1$model)
summary(MaxRich$model)


#### save model output tables using sjPlot package ####

tab_model(AbundMeanAnomalyModel1$model, transform = NULL, file = paste0(outDir, "/AbunMeanAnomNH_output_table.html"))
summary(AbundMeanAnomalyModel1$model)
R2GLMER(AbundMeanAnomalyModel1$model)

tab_model(RichMeanAnomalyModel1$model, transform = NULL, file = paste0(outDir, "/RichMeanAnomNH_output_table.html"))
summary(RichMeanAnomalyModel1$model)
R2GLMER(RichMeanAnomalyModel1$model) # use these values in the table

#$conditional
#[1] 0.638607

#$marginal
#[1] 0.01279844

tab_model(AbundMaxAnomalyModel1$model, transform = NULL, file = paste0(outDir, "/AbunMaxAnomNH_output_table.html"))
summary(AbundMaxAnomalyModel1$model)
R2GLMER(AbundMaxAnomalyModel1$model)

tab_model(RichMaxAnomalyModel1$model, transform = NULL, file = paste0(outDir, "/RichMaxAnomNH_output_table.html"))
summary(RichMaxAnomalyModel1$model)
R2GLMER(RichMaxAnomalyModel1$model) # use these values in the table

#$conditional
#[1] 0.6361979

#$marginal
#[1] 0.01214215

# For insects

# 25% NH cover = -1.015469 in rescaled values
# 50% NH cover = 0.01493712 in rescaled values
# 75% NH cover = 1.045653 in rescaled values
# 100% NH cover = 2.073849 in rescaled values


# create matrix for predictions
nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(AbundMeanAnomalyModel1$data$UI2)),
  NH_5000.rs=c(-1.015469,0.01493712,1.045653,2.073849))

# back transform the climate data range
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform NH data range
nd2$NH_5000 <- round(BackTransformCentreredPredictor(
  transformedX = nd2$NH_5000.rs,originalX = predictsSites$NH_5000)*100,0)

# set values for richness and Abundance
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# set the reference row
refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) & 
                  (nd2$NH_5000==100))

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
a.preds.tmean <- PredictGLMERRandIter(model = AbundMeanAnomalyModel1$model,data = nd2)

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



##%######################################################%##
#                                                          #
####             Figures for manuscript                 ####
#                                                          #
##%######################################################%##

# load(paste0(outDir, "MeanAnomalyModelAbun_NH.rdata"))
# load(paste0(outDir, "RichMeanAnomalyModel_NH.rdata"))
# load(paste0(outDir, "AbundMaxAnomalyModel_NH.rdata"))
# load(paste0(outDir, "RichMaxAnomalyModel_NH.rdata"))



nd2 <- expand.grid(
  StdTmeanAnomalyRS=seq(from = min(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        to = max(AbundMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                        length.out = 100),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(AbundMeanAnomalyModel1$data$UI2)),
  NH_5000.rs=c(-1.05,0.01493712,1.045653,2.073849))

# back transform the climate data range
nd2$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd2$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform NH data range
nd2$NH_5000 <- round(BackTransformCentreredPredictor(
  transformedX = nd2$NH_5000.rs,originalX = predictsSites$NH_5000)*100,0)

nd2$NH_5000[nd2$NH_5000 == 99] <- 100

# set values for richness and Abundance
nd2$LogAbund <- 0
nd2$Species_richness <- 0

# set the reference row
refRow <- which((nd2$UI2=="Primary vegetation") & 
                  nd2$StdTmeanAnomaly== min(nd2$StdTmeanAnomaly[nd2$StdTmeanAnomaly > 0]) &
                  #(nd2$StdTmeanAnomaly==min(abs(nd2$StdTmeanAnomaly))) & 
                  (nd2$NH_5000==100))

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

nd2$NH_5000 <- factor(nd2$NH_5000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd2 <- nd2[nd2$UI2 %in% c("Agriculture_Low", "Agriculture_High"), ]

# plot
p1 <- ggplot(data = nd2, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
  geom_line(aes(col = NH_5000), size = 0.75) +
  geom_ribbon(aes(ymin = nd2$PredLower, ymax = nd2$PredUpper, fill = NH_5000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 2, labeller = as_labeller(c('Agriculture_Low' = "a                       Agriculture_Low", 'Agriculture_High' = "b                        Agriculture_High"))) + 
  theme_bw() + 
  labs(fill = "% NH", col = "% NH") + 
  ylab("Change in total abundance (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 150)) + 
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        strip.text.x = element_text(hjust = 0, size = 8, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 7), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2), 
        strip.background = element_rect(size = 0.2)) 


### mean anomaly species richness 
nd <- expand.grid(
 StdTmeanAnomalyRS=seq(from = min(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                       to = max(RichMeanAnomalyModel1$data$StdTmeanAnomalyRS),
                       length.out = 100),
 UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
            levels = levels(RichMeanAnomalyModel1$data$UI2)),
 NH_5000.rs=c(-1.05,0.01493712,1.045653,2.073849))
#NH_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

# back transform the climate data range
nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd$StdTmeanAnomalyRS,
  originalX = predictsSites$StdTmeanAnomaly)

# back transform NH data range
nd$NH_5000 <- round(BackTransformCentreredPredictor(
  transformedX = nd$NH_5000.rs,originalX = predictsSites$NH_5000)*100,0)

nd$NH_5000[nd2$NH_5000 == 99] <- 100

# set values for richness and Abundance
nd$LogAbund <- 0
nd$Species_richness <- 0

# set the reference row
refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))) &
                  (nd$NH_5000==100))

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

nd$NH_5000 <- factor(nd$NH_5000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd <- nd[nd$UI2 %in% c("Agriculture_Low", "Agriculture_High"), ]

# plot
p2 <- ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) +
  geom_line(aes(col = NH_5000), size = 0.75) +
  geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = NH_5000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed", size = 0.2) +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 2, labeller = as_labeller(c('Agriculture_Low' = "c                       Agriculture_Low", 'Agriculture_High' = "d                        Agriculture_High"))) +
  theme_bw() +
  labs(fill = "% NH", col = "% NH") +
  ylab("Change in species richness (%)") +
  xlab("Standardised Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 150)) +
  theme(aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        axis.title = element_text(size = 7),
        strip.text.x = element_text(hjust = 0, size = 8, face = "bold"),
        legend.position = "right",
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 7), 
        panel.grid.minor = element_blank(),
        panel.grid.major = element_line(size = 0.2),
        panel.border = element_rect(size = 0.2), 
        axis.ticks = element_line(size = 0.2), 
        strip.background = element_rect(size = 0.2))



cowplot::plot_grid(p1, p2, nrow = 2)

# save
ggsave(filename = paste0(outDir, "Figure_4_absr_NH.pdf"), plot = last_plot(), width = 180, height = 170, units = "mm", dpi = 300)



#### extended data - max anom

# The interaction between NH, climate and LU is not significant so this interaction is not plotted. 

# 
# nd2 <- expand.grid(
#   StdTmaxAnomalyRS=seq(from = min(AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS),
#                        to = max(AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS),
#                        length.out = 200),
#   UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
#              levels = levels(AbundMaxAnomalyModel1$data$UI2)),
#   NH_5000.rs=c(-1.05,0.01493712,1.045653,2.073849))
# # NH_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))
# nd2$StdTmaxAnomaly <- BackTransformCentreredPredictor(
#   transformedX = nd2$StdTmaxAnomalyRS,
#   originalX = predictsSites$StdTmaxAnomaly)
# nd2$NH_5000 <- round(BackTransformCentreredPredictor(
#   transformedX = nd2$NH_5000.rs,originalX = predictsSites$NH_5000)*100,0)
# 
# nd2$NH_5000[nd2$NH_5000  == 99] <- 100
# nd2$LogAbund <- 0
# nd2$Species_richness <- 0
# 
# refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomaly==min(abs(nd2$StdTmaxAnomaly))) & 
#                   (nd2$NH_5000==100))
# 
# exclQuantiles <- c(0.025,0.975)
# 
# QPV <- quantile(x = AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS[
#   AbundMaxAnomalyModel1$data$UI2=="Primary vegetation"],
#   probs = exclQuantiles)
# QSV <- quantile(x = AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS[
#   AbundMaxAnomalyModel1$data$UI2=="Secondary vegetation"],
#   probs = exclQuantiles)
# QAL <- quantile(x = AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS[
#   AbundMaxAnomalyModel1$data$UI2=="Agriculture_Low"],
#   probs = exclQuantiles)
# QAH <- quantile(x = AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS[
#   AbundMaxAnomalyModel1$data$UI2=="Agriculture_High"],
#   probs = exclQuantiles)
# 
# a.preds.tmax <- PredictGLMERRandIter(model = AbundMaxAnomalyModel1$model,data = nd2, nIters = 10000)
# 
# a.preds.tmax <- exp(a.preds.tmax)-0.01
# 
# a.preds.tmax <- sweep(x = a.preds.tmax,MARGIN = 2,STATS = a.preds.tmax[refRow,],FUN = '/')
# 
# a.preds.tmax[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomalyRS < QPV[1])),] <- NA
# a.preds.tmax[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomalyRS > QPV[2])),] <- NA
# a.preds.tmax[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmaxAnomalyRS < QSV[1])),] <- NA
# a.preds.tmax[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmaxAnomalyRS > QSV[2])),] <- NA
# a.preds.tmax[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmaxAnomalyRS < QAL[1])),] <- NA
# a.preds.tmax[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmaxAnomalyRS > QAL[2])),] <- NA
# a.preds.tmax[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmaxAnomalyRS < QAH[1])),] <- NA
# a.preds.tmax[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmaxAnomalyRS > QAH[2])),] <- NA
# 
# nd2$PredMedian <- ((apply(X = a.preds.tmax,MARGIN = 1,
#                           FUN = median,na.rm=TRUE))*100)-100
# nd2$PredUpper <- ((apply(X = a.preds.tmax,MARGIN = 1,
#                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd2$PredLower <- ((apply(X = a.preds.tmax,MARGIN = 1,
#                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# 
# # set factor levels
# nd2$UI2 <- factor(nd2$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))
# 
# nd2$NH_5000 <- factor(nd2$NH_5000, levels = c("100", "75", "50", "25"))
# 
# # just take agriculture values
# nd2 <- nd2[nd2$UI2 %in% c("Agriculture_Low", "Agriculture_High"), ]
# 
# # plot
# p1 <- ggplot(data = nd2, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
#   geom_line(aes(col = NH_5000), size = 1) +
#   geom_ribbon(aes(ymin = nd2$PredLower, ymax = nd2$PredUpper, fill = NH_5000), alpha = 0.2) +
#   geom_hline(yintercept =  0, lty = "dashed") + 
#   scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
#   scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
#   facet_wrap(~UI2, ncol = 2, labeller = as_labeller(c('Agriculture_Low' = "a              Agriculture_Low", 'Agriculture_High' = "b              Agriculture_High"))) + 
#   theme_bw() + 
#   labs(fill = "% NH", col = "% NH") + 
#   ylab("Change in total abundance (%)") +
#   xlab("Maximum Temperature Anomaly") +
#   xlim(c(-0.5, 2)) +
#   ylim(c(-100, 150)) + 
#   theme(aspect.ratio = 1, text = element_text(size = 12),
#         strip.text.x = element_text(hjust = 0, size = 12, face = "bold"))
# 
# 
# #ggsave(filename = paste0(outDir, "Extended_Data6_MaxAnomAbun_NH.pdf"), height = 4, width = 8)


# richness and max anomaly

# richness model

nd3 <- expand.grid(
  StdTmaxAnomalyRS=seq(from = min(RichMaxAnomalyModel1$data$StdTmaxAnomalyRS),
                       to = max(RichMaxAnomalyModel1$data$StdTmaxAnomalyRS),
                       length.out = 200),
  UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
             levels = levels(RichMaxAnomalyModel1$data$UI2)),
  NH_5000.rs=c(-1.05,0.01493712,1.045653,2.073849))
# NH_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))

nd3$StdTmaxAnomaly <- BackTransformCentreredPredictor(
  transformedX = nd3$StdTmaxAnomalyRS,
  originalX = predictsSites$StdTmaxAnomaly)

nd3$NH_5000 <- round(BackTransformCentreredPredictor(
  transformedX = nd3$NH_5000.rs,originalX = predictsSites$NH_5000)*100,0)

nd3$NH_5000[nd3$NH_5000 == 99] <- 100

nd3$LogAbund <- 0
nd3$Species_richness <- 0

refRow <- which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomaly==min(abs(nd3$StdTmaxAnomaly))) & 
                  (nd3$NH_5000==100))

exclQuantiles <- c(0.025,0.975)

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

s.preds.tmax <- PredictGLMERRandIter(model = RichMaxAnomalyModel1$model,data = nd3, nIters = 10000)
s.preds.tmax <- exp(s.preds.tmax)

s.preds.tmax <- sweep(x = s.preds.tmax,MARGIN = 2,STATS = s.preds.tmax[refRow,],FUN = '/')

s.preds.tmax[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomalyRS < QPV[1])),] <- NA
s.preds.tmax[which((nd3$UI2=="Primary vegetation") & (nd3$StdTmaxAnomalyRS > QPV[2])),] <- NA
s.preds.tmax[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmaxAnomalyRS < QSV[1])),] <- NA
s.preds.tmax[which((nd3$UI2=="Secondary vegetation") & (nd3$StdTmaxAnomalyRS > QSV[2])),] <- NA
s.preds.tmax[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmaxAnomalyRS < QAL[1])),] <- NA
s.preds.tmax[which((nd3$UI2=="Agriculture_Low") & (nd3$StdTmaxAnomalyRS > QAL[2])),] <- NA
s.preds.tmax[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmaxAnomalyRS < QAH[1])),] <- NA
s.preds.tmax[which((nd3$UI2=="Agriculture_High") & (nd3$StdTmaxAnomalyRS > QAH[2])),] <- NA

nd3$PredMedian <- ((apply(X = s.preds.tmax,MARGIN = 1,
                         FUN = median,na.rm=TRUE))*100)-100
nd3$PredUpper <- ((apply(X = s.preds.tmax,MARGIN = 1,
                        FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
nd3$PredLower <- ((apply(X = s.preds.tmax,MARGIN = 1,
                        FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100

# set factor levels
nd3$UI2 <- factor(nd3$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))

nd3$NH_5000 <- factor(nd3$NH_5000, levels = c("100", "75", "50", "25"))

# just take agriculture values
nd3 <- nd3[nd3$UI2 %in% c("Agriculture_Low", "Agriculture_High"), ]

# plot
p2 <-ggplot(data = nd3, aes(x = StdTmaxAnomaly, y = PredMedian)) + 
  geom_line(aes(col = NH_5000), size = 0.75) +
  geom_ribbon(aes(ymin = nd3$PredLower, ymax = nd3$PredUpper, fill = NH_5000), alpha = 0.2) +
  geom_hline(yintercept = 0, lty = "dashed") +
  scale_fill_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  scale_colour_manual(values = rev(c("#a50026","#f46d43","#74add1","#313695"))) +
  facet_wrap(~UI2, ncol = 2, labeller = as_labeller(c('Agriculture_Low' = "a                                Agriculture_Low", 'Agriculture_High' = "b                                Agriculture_High"))) + 
  theme_bw() + 
  labs(fill = "% NH", col = "% NH") + 
  ylab("Change in species richness (%)") +
  xlab("Maximum Temperature Anomaly") +
  xlim(c(-0.5, 2)) +
  ylim(c(-100, 150)) + 
  theme(aspect.ratio = 1, 
        #title = element_text(size = 8, face = "bold"),
        strip.text.x = element_text(hjust = 0, size = 8, face = "bold"),
        axis.text = element_text(size = 7),
        text = element_text(size = 7),
        legend.position = "right",
        legend.text = element_text(size = 7), 
        legend.title = element_text(size = 7))


# abundance plot non-significant so not including in figure. 
# organise plots together
#cowplot::plot_grid(p1, p2, nrow = 2)

# save
ggsave(filename = paste0(outDir, "Extended_Data6_sr_MAX.pdf"), height = 4, width = 8)
ggsave(filename = paste0(outDir, "Extended_Data6_sr_MAX.jpeg"), plot = last_plot(), width = 183, height = 150, units = "mm", dpi = 300)


#### plotting using base R if preferred


# pdf(file = paste0(outDir,"AbundanceMeanAnomaly.pdf"),width = 20/2.54,height = 10/2.54)
# 
# par(mfrow=c(1,2))
# 
# 
# nd2_ab_mean <- nd2
# 
# ## 1. plot for Abundance response to mean anomaly wuth 
# 
# ylims <- with(nd2[nd2$UI2 %in% c("Agriculture_Low","Agriculture_High"),],
#               c(min(PredLower,na.rm = TRUE),max(PredUpper,na.rm = TRUE)))
# xlims <- with(nd2[nd2$UI2 %in% c("Agriculture_Low","Agriculture_High"),],
#               c(min(StdTmeanAnomaly,na.rm = TRUE),max(StdTmeanAnomaly,na.rm = TRUE)))
# 
# invisible(lapply(X = split(x = nd2,f = nd2$UI2)[c(3,2)],FUN = function(preds.lu){
#   
#   plot(-9e99,-9e99,xlim=xlims,ylim=ylims,
#        xlab="Mean temperature anomaly",ylab="Abundance (%)", cex.lab = 0.8, cex.axis = 0.8)
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
#     
#   },split(x = preds.lu,f = preds.lu$NH_5000),c("#a50026","#f46d43","#74add1","#313695")))
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
#   
# }))
# 
# # add legend to righthand plot
# legend(
#   x = 1.5,y = 112, bty="n",
#   legend = c("25%", "50%", "75%", "100%"),
#   col = c("#a50026","#f46d43","#74add1","#313695"),
#   lty=1,lwd=2, cex = 0.7, title = "% NH", title.adj = 0.2)
# 
# 
# invisible(dev.off())
# 
# 
# 
# pdf(file = paste0(outDir,"RichnessMeanAnomaly.pdf"),width = 17.5/2.54,height = 8/2.54)
# 
# QPV <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
#   RichMeanAnomalyModel1$data$UI2=="Primary vegetation"],
#   probs = exclQuantiles)
# QSV <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
#   RichMeanAnomalyModel1$data$UI2=="Secondary vegetation"],
#   probs = exclQuantiles)
# QAL <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
#   RichMeanAnomalyModel1$data$UI2=="Agriculture_Low"],
#   probs = exclQuantiles)
# QAH <- quantile(x = RichMeanAnomalyModel1$data$StdTmeanAnomalyRS[
#   RichMeanAnomalyModel1$data$UI2=="Agriculture_High"],
#   probs = exclQuantiles)
# 
# s.preds.tmean <- PredictGLMERRandIter(model = RichMeanAnomalyModel1$model,data = nd2)
# s.preds.tmean <- exp(s.preds.tmean)
# 
# s.preds.tmean <- sweep(x = s.preds.tmean,MARGIN = 2,STATS = s.preds.tmean[refRow,],FUN = '/')
# 
# s.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS < QPV[1])),] <- NA
# s.preds.tmean[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmeanAnomalyRS > QPV[2])),] <- NA
# s.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS < QSV[1])),] <- NA
# s.preds.tmean[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmeanAnomalyRS > QSV[2])),] <- NA
# s.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS < QAL[1])),] <- NA
# s.preds.tmean[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmeanAnomalyRS > QAL[2])),] <- NA
# s.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS < QAH[1])),] <- NA
# s.preds.tmean[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmeanAnomalyRS > QAH[2])),] <- NA
# 
# nd2$PredMedian <- ((apply(X = s.preds.tmean,MARGIN = 1,
#                           FUN = median,na.rm=TRUE))*100)-100
# nd2$PredUpper <- ((apply(X = s.preds.tmean,MARGIN = 1,
#                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd2$PredLower <- ((apply(X = s.preds.tmean,MARGIN = 1,
#                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# par(mfrow=c(1,2))
# 
# ylims <- with(nd2[nd2$UI2 %in% c("Agriculture_Low","Agriculture_High"),],
#               c(min(PredLower,na.rm = TRUE),max(PredUpper,na.rm = TRUE)))
# xlims <- with(nd2[nd2$UI2 %in% c("Agriculture_Low","Agriculture_High"),],
#               c(min(StdTmeanAnomaly,na.rm = TRUE),max(StdTmeanAnomaly,na.rm = TRUE)))
# 
# invisible(lapply(X = split(x = nd2,f = nd2$UI2)[c(3,2)],FUN = function(preds.lu){
#   
#   plot(-9e99,-9e99,xlim=xlims,ylim=ylims,
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
#   },split(x = preds.lu,f = preds.lu$NH_5000),c("#a50026","#f46d43","#74add1","#313695")))
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
# }))
# 
# invisible(dev.off())
# 
# 
# nd2 <- expand.grid(
#   StdTmaxAnomalyRS=seq(from = min(AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS),
#                        to = max(AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS),
#                        length.out = 100),
#   UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
#              levels = levels(AbundMaxAnomalyModel1$data$UI2)),
#   NH_5000.rs=c(-1.015469,0.01493712,1.045653,2.073849))
# # NH_5000.rs=c(-1.122627,-0.1351849,0.8498366,1.834235))
# nd2$StdTmaxAnomaly <- BackTransformCentreredPredictor(
#   transformedX = nd2$StdTmaxAnomalyRS,
#   originalX = predictsSites$StdTmaxAnomaly)
# nd2$NH_5000 <- round(BackTransformCentreredPredictor(
#   transformedX = nd2$NH_5000.rs,originalX = predictsSites$NH_5000)*100,0)
# nd2$LogAbund <- 0
# nd2$Species_richness <- 0
# 
# refRow <- which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomaly==min(abs(nd2$StdTmaxAnomaly))) & 
#                   (nd2$NH_5000==100))
# 
# exclQuantiles <- c(0.025,0.975)
# 
# pdf(file = paste0(outDir,"AbundanceMaxAnomaly.pdf"),width = 17.5/2.54,height = 8/2.54)
# 
# QPV <- quantile(x = AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS[
#   AbundMaxAnomalyModel1$data$UI2=="Primary vegetation"],
#   probs = exclQuantiles)
# QSV <- quantile(x = AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS[
#   AbundMaxAnomalyModel1$data$UI2=="Secondary vegetation"],
#   probs = exclQuantiles)
# QAL <- quantile(x = AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS[
#   AbundMaxAnomalyModel1$data$UI2=="Agriculture_Low"],
#   probs = exclQuantiles)
# QAH <- quantile(x = AbundMaxAnomalyModel1$data$StdTmaxAnomalyRS[
#   AbundMaxAnomalyModel1$data$UI2=="Agriculture_High"],
#   probs = exclQuantiles)
# 
# a.preds.tmax <- PredictGLMERRandIter(model = AbundMaxAnomalyModel1$model,data = nd2)
# a.preds.tmax <- exp(a.preds.tmax)-0.01
# 
# a.preds.tmax <- sweep(x = a.preds.tmax,MARGIN = 2,STATS = a.preds.tmax[refRow,],FUN = '/')
# 
# a.preds.tmax[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomalyRS < QPV[1])),] <- NA
# a.preds.tmax[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomalyRS > QPV[2])),] <- NA
# a.preds.tmax[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmaxAnomalyRS < QSV[1])),] <- NA
# a.preds.tmax[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmaxAnomalyRS > QSV[2])),] <- NA
# a.preds.tmax[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmaxAnomalyRS < QAL[1])),] <- NA
# a.preds.tmax[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmaxAnomalyRS > QAL[2])),] <- NA
# a.preds.tmax[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmaxAnomalyRS < QAH[1])),] <- NA
# a.preds.tmax[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmaxAnomalyRS > QAH[2])),] <- NA
# 
# nd2$PredMedian <- ((apply(X = a.preds.tmax,MARGIN = 1,
#                           FUN = median,na.rm=TRUE))*100)-100
# nd2$PredUpper <- ((apply(X = a.preds.tmax,MARGIN = 1,
#                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd2$PredLower <- ((apply(X = a.preds.tmax,MARGIN = 1,
#                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# par(mfrow=c(1,2))
# 
# ylims <- with(nd2[nd2$UI2 %in% c("Agriculture_Low","Agriculture_High"),],
#               c(min(PredLower,na.rm = TRUE),max(PredUpper,na.rm = TRUE)))
# xlims <- with(nd2[nd2$UI2 %in% c("Agriculture_Low","Agriculture_High"),],
#               c(min(StdTmaxAnomaly,na.rm = TRUE),max(StdTmaxAnomaly,na.rm = TRUE)))
# 
# invisible(lapply(X = split(x = nd2,f = nd2$UI2)[c(3,2)],FUN = function(preds.lu){
#   
#   plot(-9e99,-9e99,xlim=xlims,ylim=ylims,
#        xlab="Maximum temperature anomaly",ylab="Abundance (%)")
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
#   },split(x = preds.lu,f = preds.lu$NH_5000),c("#a50026","#f46d43","#74add1","#313695")))
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
# }))
# 
# invisible(dev.off())
# 
# 
# 
# pdf(file = paste0(outDir,"RichnessMaxAnomaly.pdf"),width = 17.5/2.54,height = 8/2.54)
# 
# QPV <- quantile(x = RichMaxAnomalyModel1$data$StdTmaxAnomalyRS[
#   RichMaxAnomalyModel1$data$UI2=="Primary vegetation"],
#   probs = exclQuantiles)
# QSV <- quantile(x = RichMaxAnomalyModel1$data$StdTmaxAnomalyRS[
#   RichMaxAnomalyModel1$data$UI2=="Secondary vegetation"],
#   probs = exclQuantiles)
# QAL <- quantile(x = RichMaxAnomalyModel1$data$StdTmaxAnomalyRS[
#   RichMaxAnomalyModel1$data$UI2=="Agriculture_Low"],
#   probs = exclQuantiles)
# QAH <- quantile(x = RichMaxAnomalyModel1$data$StdTmaxAnomalyRS[
#   RichMaxAnomalyModel1$data$UI2=="Agriculture_High"],
#   probs = exclQuantiles)
# 
# s.preds.tmax <- PredictGLMERRandIter(model = RichMaxAnomalyModel1$model,data = nd2)
# s.preds.tmax <- exp(s.preds.tmax)-0.01
# 
# s.preds.tmax <- sweep(x = s.preds.tmax,MARGIN = 2,STATS = s.preds.tmax[refRow,],FUN = '/')
# 
# s.preds.tmax[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomalyRS < QPV[1])),] <- NA
# s.preds.tmax[which((nd2$UI2=="Primary vegetation") & (nd2$StdTmaxAnomalyRS > QPV[2])),] <- NA
# s.preds.tmax[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmaxAnomalyRS < QSV[1])),] <- NA
# s.preds.tmax[which((nd2$UI2=="Secondary vegetation") & (nd2$StdTmaxAnomalyRS > QSV[2])),] <- NA
# s.preds.tmax[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmaxAnomalyRS < QAL[1])),] <- NA
# s.preds.tmax[which((nd2$UI2=="Agriculture_Low") & (nd2$StdTmaxAnomalyRS > QAL[2])),] <- NA
# s.preds.tmax[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmaxAnomalyRS < QAH[1])),] <- NA
# s.preds.tmax[which((nd2$UI2=="Agriculture_High") & (nd2$StdTmaxAnomalyRS > QAH[2])),] <- NA
# 
# nd2$PredMedian <- ((apply(X = s.preds.tmax,MARGIN = 1,
#                           FUN = median,na.rm=TRUE))*100)-100
# nd2$PredUpper <- ((apply(X = s.preds.tmax,MARGIN = 1,
#                          FUN = quantile,probs = 0.975,na.rm=TRUE))*100)-100
# nd2$PredLower <- ((apply(X = s.preds.tmax,MARGIN = 1,
#                          FUN = quantile,probs = 0.025,na.rm=TRUE))*100)-100
# 
# par(mfrow=c(1,2))
# 
# ylims <- with(nd2[nd2$UI2 %in% c("Agriculture_Low","Agriculture_High"),],
#               c(min(PredLower,na.rm = TRUE),max(PredUpper,na.rm = TRUE)))
# xlims <- with(nd2[nd2$UI2 %in% c("Agriculture_Low","Agriculture_High"),],
#               c(min(StdTmaxAnomaly,na.rm = TRUE),max(StdTmaxAnomaly,na.rm = TRUE)))
# 
# invisible(lapply(X = split(x = nd2,f = nd2$UI2)[c(3,2)],FUN = function(preds.lu){
#   
#   plot(-9e99,-9e99,xlim=xlims,ylim=ylims,
#        xlab="Maximum temperature anomaly",ylab="Richness (%)")
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
#   },split(x = preds.lu,f = preds.lu$NH_5000),c("#a50026","#f46d43","#74add1","#313695")))
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
# }))
# 
# invisible(dev.off())
# 

t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()