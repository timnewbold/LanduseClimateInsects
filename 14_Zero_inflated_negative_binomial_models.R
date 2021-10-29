##%######################################################%##
#                                                          #
####  Additional tests: Trialling zero-inflated models  ####
#                                                          #
##%######################################################%##

# look into this package: https://cran.r-project.org/web/packages/glmmTMB/vignettes/glmmTMB.pdf

# start by testing a zero inflated negative binomial model for total abundance
# then try just zero inflated versions of the scaled abundance and richness models


rm(list = ls())

# load libraries
library(devtools)
library(StatisticalModels)
library(predictsFunctions)
library(glmmTMB)
source("Functions.R")
library(DHARMa)
library(cowplot)
library(sjPlot)


# directories 
predictsDataDir <- "6_RunLUClimateModels/"
outDir <- "14_Additional_Tests/"


# load in the data
predictsSites <- readRDS(file = paste0(predictsDataDir,"PREDICTSSiteData.rds"))

predictsSites <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ] # 6069

predictsSites <- predictsSites[!is.na(predictsSites$Total_abundance), ]
#### trying out the abundance model using a zero inflated negative binomial model ####

# total abundacne glmmTMB model
abun_NB_mod <- glmmTMB::glmmTMB(Total_abundance ~ UI2 + poly(StdTmeanAnomalyRS,1) + UI2:poly(StdTmeanAnomalyRS,1) + (1|SS) + (1|SSB), data = predictsSites, ziformula = ~1, family = list(family = "nbinom2", link = "log"))

# save the model output
save(abun_NB_mod, file = paste0(outDir, "AbunMeanAnom_zeroinf_negbinom.rdata"))
#load(paste0(outDir, "AbunMeanAnom_zeroinf_negbinom.rdata"))



#### compare coefficients with original model ####

load("6_RunLUClimateModels/MeanAnomalyModelAbund.rdata") # MeanAnomalyModelAbund

summary(MeanAnomalyModelAbund$model)
summary(abun_NB_mod)

# are they comparable though?
#AIC(MeanAnomalyModelAbund$model, abun_NB_mod, abunRS_NB_mod)


# present model outputs

tab_model(abun_NB_mod, MeanAnomalyModelAbund$model, 
          transform = NULL, 
          show.ci = F,
          show.icc = F,
          show.ngroups = F,
          show.obs =  F, 
          file = paste0(outDir, "Abundance_negbinom_Zinf_results_tab.html"))

#### model checks ####


pdf(NULL)
dev.control(displaylist="enable")
plot(fitted(abun_NB_mod), residuals(abun_NB_mod))
p1 <- recordPlot()
invisible(dev.off())


pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(abun_NB_mod), main = "")
qqline(resid(abun_NB_mod))
p2 <- recordPlot()
invisible(dev.off())

cowplot::plot_grid(p1,p2,
                   labels = c("A.", "B."))#

ggsave(file = paste0(outDir, "Ab_negbinom_Zinf_model_checks.pdf"), height = 5, width = 10) 



# these are generated in each of the functions below, run first to reduce time
abun_NB_mod_simres <- simulateResiduals(abun_NB_mod)

# Q-Q plot and residuals vs predicts
plot(abun_NB_mod_simres)


testOutliers(abun_NB_mod_simres)
testDispersion(abun_NB_mod_simres)
testUniformity(abun_NB_mod_simres)
testQuantiles(abun_NB_mod_simres)
testCategorical(abun_NB_mod_simres, catPred = predictsSites$UI2)


# table a look at these tests for the original model

plot(MeanAnomalyModelAbund$model)
qqnorm(resid(MeanAnomalyModelAbund$model), main = "")
qqline(resid(MeanAnomalyModelAbund$model))


MeanAnomalyModelAbund_simres <- simulateResiduals(MeanAnomalyModelAbund$model)

plot(MeanAnomalyModelAbund_simres)

testZeroInflation(MeanAnomalyModelAbund_simres)

testOutliers(MeanAnomalyModelAbund_simres)
testDispersion(MeanAnomalyModelAbund_simres)
testUniformity(MeanAnomalyModelAbund_simres)
testQuantiles(MeanAnomalyModelAbund_simres)
testCategorical(MeanAnomalyModelAbund_simres, catPred = predictsSites$UI2)




#### Plot results of glmmTMB model ####


pdf(file = paste0(outDir, "/Figure2_MeanAnom_Abun_ZInegbin.pdf"), width = 4, height = 4)

par(mfrow=c(1,1))
par(las=1)
par(mgp=c(1.6,0.3,0)) # 
par(mar=c(2.6,2.6,1,1)) # margins around plot
par(tck=-0.01) # tick mark size
#par(pty="s") # set plot type to be square

exclQuantiles <- c(0.025,0.975)

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
nd$Total_abundance <- 0

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
a.preds.tmean <- PredictGLMERRandIter_TMB(model = abun_NB_mod, data = nd)

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

title("zero-inflated negative binomial model", adj = 0, cex.main = 1)


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


dev.off()




##%######################################################%##
#                                                          #
####      zero inflated model for species richness      ####
#                                                          #
##%######################################################%##

# load in the data
predictsSites <- readRDS(file = paste0(predictsDataDir,"PREDICTSSiteData.rds"))

predictsSites <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ] # 6069


#### richness model, zero inflated poisson ####
SR_Zinf_mod <- glmmTMB::glmmTMB(Species_richness ~ UI2 + poly(StdTmeanAnomalyRS,1) + UI2:poly(StdTmeanAnomalyRS,1) + (1|SS) + (1|SSB) + (1|SSBS), data = predictsSites, ziformula = ~1, family=poisson)

# save the model output
save(SR_Zinf_mod, file = paste0(outDir, "RichnessMeanAnom_zeroinf.rdata"))
#load(paste0(outDir, "RichnessMeanAnom_zeroinf.rdata"))


fixed.effects(SR_Zinf_mod)

#plot(abun_NB_mod)
qqnorm(resid(SR_Zinf_mod), main = "")
qqline(resid(SR_Zinf_mod))


pdf(NULL)
dev.control(displaylist="enable")
plot(fitted(SR_Zinf_mod), residuals(SR_Zinf_mod))
p1 <- recordPlot()
invisible(dev.off())


pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(SR_Zinf_mod), main = "")
qqline(resid(SR_Zinf_mod))
p2 <- recordPlot()
invisible(dev.off())

cowplot::plot_grid(p1,p2,
                   labels = c("A.", "B."))#

ggsave(file = paste0(outDir, "SR_Zinf_model_checks.pdf"), height = 5, width = 10) 


# these are generated in each of the functions below, run first to reduce time
SR_Zinf_mod_simres <- simulateResiduals(SR_Zinf_mod)

# Q-Q plot and residuals vs predicts
plot(SR_Zinf_mod_simres)


# compare to the original model
load("6_RunLUClimateModels/MeanAnomalyModelRich.rdata")


qqnorm(resid(MeanAnomalyModelRich$model), main = "")
qqline(resid(MeanAnomalyModelRich$model))

plot(fitted(MeanAnomalyModelRich$model), residuals(MeanAnomalyModelRich$model))

# zero-inflated model has more positively skewed q-q plot than the original model

# compare coeffs

fixef(SR_Zinf_mod)
fixef(MeanAnomalyModelRich$model)

# coefficents are very similar


# present model outputs
library(sjPlot)
tab_model(SR_Zinf_mod, MeanAnomalyModelRich$model, 
          transform = NULL, 
          show.ci = F,
          show.icc = F,
          show.ngroups = F,
          show.obs =  F, 
          file = paste0(outDir, "Richness_Zinf_results_tab.html"))



##%######################################################%##
#                                                          #
####          abundance model, zero inflated            ####
#                                                          #
##%######################################################%##


model_data <- predictsSites[!is.na(predictsSites$LogAbund), ]
model_data <- droplevels(model_data)

Ab_Zinf_mod <- glmmTMB::glmmTMB(LogAbund ~ UI2 + poly(StdTmeanAnomalyRS,1) + UI2:poly(StdTmeanAnomalyRS,1) + (1|SS) + (1|SSB), data = model_data, ziformula = ~1, family=gaussian)

# save the model output
save(Ab_Zinf_mod, file = paste0(outDir, "AbundanceMeanAnom_zeroinf.rdata"))
#load(paste0(outDir, "AbundanceMeanAnom_zeroinf.rdata"))


fixed.effects(Ab_Zinf_mod)

### cant extract residuals for this model for some reason, so cannot generate plots. Doesn't work using Dharma package either. 

#plot(abun_NB_mod)
qqnorm(residuals(Ab_Zinf_mod, type = "response"), main = "")
qqline(resid(Ab_Zinf_mod))


pdf(NULL)
dev.control(displaylist="enable")
plot(fitted(Ab_Zinf_mod), residuals(Ab_Zinf_mod))
p1 <- recordPlot()
invisible(dev.off())


pdf(NULL)
dev.control(displaylist="enable")
qqnorm(resid(Ab_Zinf_mod), main = "")
qqline(resid(Ab_Zinf_mod))
p2 <- recordPlot()
invisible(dev.off())

cowplot::plot_grid(p1,p2,
                   labels = c("A.", "B."))#

ggsave(file = paste0(outDir, "Ab_Zinf_model_checks.pdf"), height = 5, width = 10) 


# these are generated in each of the functions below, run first to reduce time
Ab_Zinf_mod_simres <- simulateResiduals(Ab_Zinf_mod)

# Q-Q plot and residuals vs predicts
plot(Ab_Zinf_mod_simres)


# compare to the original model
load("6_RunLUClimateModels/MeanAnomalyModelAbund.rdata")


qqnorm(resid(MeanAnomalyModelAbund$model), main = "")
qqline(resid(MeanAnomalyModelAbund$model))

plot(fitted(MeanAnomalyModelAbund$model), residuals(MeanAnomalyModelAbund$model))


# compare coeffs

fixef(Ab_Zinf_mod)
fixef(MeanAnomalyModelAbund$model)

# coefficents are very similar

# present model outputs
tab_model(Ab_Zinf_mod, MeanAnomalyModelAbund$model, 
          transform = NULL, 
          show.ci = F,
          show.icc = F,
          show.ngroups = F,
          show.obs =  F, 
          file = paste0(outDir, "Abundance_Zinf_results_tab.html"))



