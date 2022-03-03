##%######################################################%##
#                                                          #
####        13. Exploring data bias in CRU data         ####
#                                                          #
##%######################################################%##

# In this script, I will take a look at the data coverage of the CRU data
# using the stn variable within the mean temp data

rm(list = ls())


# organise directories
datadir <- "0_data/"
outdir <- "13_Exploring_CRU_data/"
dir.create(outdir)

sink(paste0(outdir,"log.txt"))

t.start <- Sys.time()

print(t.start)

# load libraries
library(raster)
library(predictsFunctions)
library(StatisticalModels)
source("Functions.R")

# load in the stn data
stn <- stack(paste0(datadir,"/cru_ts4.03.1901.2018.tmp.dat.nc"),varname = "stn")

# take the baseline years and take a look
stn_0130 <- stn[[1:360]]

stn_0130_mean <- mean(stn_0130)
plot(stn_0130_mean, main = "mean number of stations \nacross baseline 1901-1930")

stn_0130_range <- range(stn_0130)
plot(stn_0130_range)

# load in the predicts data and do some summaries for the sites
predictsDir <- "6_RunLUClimateModels/"

# read in the predicts sites - insect subset
predicts_sites <- readRDS(paste0(predictsDir,"PREDICTSSiteData.rds"))

# Remove sites without coordinates
predicts_sites2 <- predicts_sites[!is.na(predicts_sites$Latitude), ] # this has already been done elsewhere

predicts_sites2 <- predicts_sites2[!is.na(predicts_sites2$StdTmeanAnomaly), ] # 6069

wgs84 <- crs(stn)
# Create spatial map of PREDICTS sites
predicts_sp <- SpatialPointsDataFrame(
  coords = cbind(predicts_sites2$Longitude, predicts_sites2$Latitude), 
  data = predicts_sites2, proj4string = wgs84)

# add the sites onto the plots
plot(stn_0130_mean, main = "Mean number of stations \nacross baseline 1901-1930")
plot(predicts_sp, add = T)
abline(h = -23.44, lty = "dashed")
abline(h = 23.44, lty = "dashed")


# extract just from the mean number of station?
nstn <- extract(stn_0130_mean, predicts_sp)

predicts_sp$stn <- nstn

predicts_sp_stn <- predicts_sp[predicts_sp$stn >= 6, ]

plot(stn_0130_mean, main = "Mean number of stations \nacross baseline 1901-1930")
plot(predicts_sp_stn, add = T)
abline(h = -23.44, lty = "dashed")
abline(h = 23.44, lty = "dashed")

# round to whole numbers  
#nstn <- round(nstn, digits = 0)
# take a look at spread
#table(nstn)




### look at just those within tropics ###

# add values to predicts data
predicts_sites2$nstn <- round(nstn, digits = 0)

trop_data <- predicts_sites2[predicts_sites2$Latitude > -23.44 & predicts_sites2$Latitude < 23.44, ] # 1742 rows


table(trop_data$nstn)

#   0   1   2   3   4   5   6   7   8 
# 148 178 519 544  39  65   6  29 214 

test <- trop_data[trop_data$nstn >6, ] #243

test <- trop_data[!is.na(trop_data$LogAbund), ]

test <- test[test$nstn >6, ] #231

# species richness data
table(predicts_sites2$nstn)

#   0    1    2    3    4    5    6    7    8 
# 148  178  505  750  147   66   21   66 4188

test <- trop_data[trop_data$nstn >3, ] #353


# abundance data
ab_dat <- predicts_sites2[!is.na(predicts_sites2$LogAbund), ] #5735
table(ab_dat$nstn)

#  0    1    2    3    4    5    6    7    8 
#127  174  413  732   88   66   15   54 4066 

test <- ab_dat[ab_dat$nstn > 6, ]


##%######################################################%##
#                                                          #
####            Rerun models removing those             ####
####             with few stns contributing             ####
#                                                          #
##%######################################################%##

# sequentially remove sites and rerun models, remove those with 0, then 1, then 2, then 3 (just to here?)

# write a for loop to run the model removing rows

# nstn are not whole numbers here but the average values
vals <- c(0, 1, 2, 3, 4, 5, 6)

#### abundance models, all data ####

# i <- 0
for(i in vals){
  
  model_data <- predicts_sites2[predicts_sites2$nstn > i, ]
  
  model_data <- droplevels(model_data)
  
  # model_output <- GLMERSelect(modelData = model_data,responseVar = "LogAbund",
  #                             fitFamily = "gaussian",fixedFactors = "UI2",
  #                             fixedTerms = list(StdTmeanAnomalyRS=1),
  #                             randomStruct = "(1|SS)+(1|SSB)",
  #                             fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
  #                             saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))
  
  #summary(model_output$model)
  
  model_output <- GLMER(modelData = model_data,responseVar = "LogAbund", 
                        fitFamily = "gaussian",
                        fixedStruct = "UI2 * StdTmeanAnomalyRS",
                        randomStruct = "(1|SS)+(1|SSB)")
  
  
  # save model output
  save(model_output, file = paste0(outdir, "/Model_output_removing_", i, "_Abun.rdata"))
  
  # save model stats
  write.csv(model_output$stats, file = paste0(outdir, "/Model_stats_removing_", i, "_Abun.csv"))
  
  
  # plot the results
  # the 0 one doesn't get a refrow with lnegth.out set to 300
  if(i == 0 | i == 2 | i == 3 | i == 4 | i == 5){
    nd <- expand.grid(
      StdTmeanAnomalyRS=seq(from = min(model_output$data$StdTmeanAnomalyRS),
                            to = max(model_output$data$StdTmeanAnomalyRS),
                            length.out = 100),
      UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
                 levels = levels(model_output$data$UI2)))
  }else{
  
  nd <- expand.grid(
    StdTmeanAnomalyRS=seq(from = min(model_output$data$StdTmeanAnomalyRS),
                          to = max(model_output$data$StdTmeanAnomalyRS),
                          length.out = 300),
    UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
               levels = levels(model_output$data$UI2)))
  
  }
  
  # back transform the predictors
  nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
    transformedX = nd$StdTmeanAnomalyRS,
    originalX = predicts_sites2$StdTmeanAnomaly)
  
  # set richness and abundance to 0 - to be predicted
  nd$LogAbund <- 0
  nd$Species_richness <- 0
  
  # reference for % difference = primary vegetation and positive anomaly closest to 0
  refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))
  
  # adjust plot 1: mean anomaly and abundance
  exclQuantiles <- c(0.025,0.975)
  
  QPV <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Primary vegetation"],
    probs = exclQuantiles)
  QSV <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Secondary vegetation"],
    probs = exclQuantiles)
  QAL <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Agriculture_Low"],
    probs = exclQuantiles)
  QAH <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Agriculture_High"],
    probs = exclQuantiles)
  
  # predict the results
  a.preds.tmean <- PredictGLMERRandIter(model = model_output$model,data = nd, nIters = 2000)
  
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
  
  nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))
  
  # plot
  ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
    geom_line(aes(col = UI2), size = 1) +
    geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
    scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
    theme_bw() + 
    ylab("Change in Total Abundance (%)") +
    xlab("Standardised Temperature Anomaly") +
    xlim(c(-0.5, 2)) +
    ylim(c(-100, 150)) + 
    theme(aspect.ratio = 1, text = element_text(size = 12),
          strip.text.x = element_text(hjust = 0, size = 12, face = "bold"), legend.title = element_blank()) +
    ggtitle(paste0("Sites with more than ", i, " stns"))
  
  # save
  ggsave(filename = paste0(outdir, "/Plot_removing_", i, "_Abun.pdf"), height = 5, width = 6)
  

  
 rm(nd, model_data, model_output)
  
}



#### Richness models, all data ####


# i <- 0
for(i in vals){
  
  model_data <- predicts_sites2[predicts_sites2$nstn > i, ]
  
  model_data <- droplevels(model_data)
  

  model_output <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                                 fixedStruct = "UI2 * StdTmeanAnomalyRS",
                                 randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)")
  
  
  # save model output
  save(model_output, file = paste0(outdir, "/Model_output_removing_", i, "_Richness.rdata"))
  
  # save model stats
  write.csv(model_output$stats, file = paste0(outdir, "/Model_stats_removing_", i, "_Richness.csv"))
  
  
  # plot the results
  # the 0 one doesn't get a refrow with lnegth.out set to 300
  if(i == 0 | i == 2 | i == 3 | i == 4| i == 5){
    nd <- expand.grid(
      StdTmeanAnomalyRS=seq(from = min(model_output$data$StdTmeanAnomalyRS),
                            to = max(model_output$data$StdTmeanAnomalyRS),
                            length.out = 100),
      UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
                 levels = levels(model_output$data$UI2)))
  }else{
    
    nd <- expand.grid(
      StdTmeanAnomalyRS=seq(from = min(model_output$data$StdTmeanAnomalyRS),
                            to = max(model_output$data$StdTmeanAnomalyRS),
                            length.out = 300),
      UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
                 levels = levels(model_output$data$UI2)))
    
  }
  
  # back transform the predictors
  nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
    transformedX = nd$StdTmeanAnomalyRS,
    originalX = predicts_sites2$StdTmeanAnomaly)
  
  # set richness and abundance to 0 - to be predicted
  nd$LogAbund <- 0
  nd$Species_richness <- 0
  
  # reference for % difference = primary vegetation and positive anomaly closest to 0
  refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))
  
  # adjust plot 1: mean anomaly and abundance
  exclQuantiles <- c(0.025,0.975)
  
  QPV <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Primary vegetation"],
    probs = exclQuantiles)
  QSV <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Secondary vegetation"],
    probs = exclQuantiles)
  QAL <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Agriculture_Low"],
    probs = exclQuantiles)
  QAH <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Agriculture_High"],
    probs = exclQuantiles)
  
  # predict the results
  a.preds.tmean <- PredictGLMERRandIter(model = model_output$model,data = nd)
  
  # back transform the abundance values
  a.preds.tmean <- exp(a.preds.tmean)
  
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
  
  nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))
  
  # plot
  ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
    geom_line(aes(col = UI2), size = 1) +
    geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
    scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
    theme_bw() + 
    ylab("Change in Species Richness (%)") +
    xlab("Standardised Temperature Anomaly") +
    xlim(c(-0.5, 2)) +
    ylim(c(-100, 150)) + 
    theme(aspect.ratio = 1, text = element_text(size = 12),
          strip.text.x = element_text(hjust = 0, size = 12, face = "bold"), legend.title = element_blank()) +
    ggtitle(paste0("Sites with more than ", i, " stns"))
  
  # save
  ggsave(filename = paste0(outdir, "/Plot_removing_", i, "_Richness.pdf"), height = 5, width = 6)
  
  
  
  rm(nd, model_data, model_output)
  
}




#### Tropical sites, Abundance ####

vals <- c(0,1,2,3,4,5,6)


trop_data <- predicts_sites2[predicts_sites2$Latitude > -23.44 & predicts_sites2$Latitude < 23.44, ] # 1742 rows

table(trop_data[!is.na(trop_data$LogAbund), "nstn"])
# 0   1   2   3   4   5   7   8 
# 127 174 367 586  39  65  17 214 


# i <- 0
for(i in vals){
  
  model_data <- trop_data[trop_data$nstn > i, ]
  
  model_data <- model_data[!is.na(model_data$LogAbund), ]
  
  model_data <- droplevels(model_data)
  

  model_output <- GLMER(modelData = model_data,responseVar = "LogAbund", 
                        fitFamily = "gaussian",
                        fixedStruct = "UI2 * StdTmeanAnomalyRS",
                        randomStruct = "(1|SS)+(1|SSB)")
  # save model output
  save(model_output, file = paste0(outdir, "/Model_output_removing_", i, "_Abun_Tropical.rdata"))
  
  # save model stats
  write.csv(model_output$stats, file = paste0(outdir, "/Model_stats_removing_", i, "_Abun_Tropical.csv"))
  
  
  # plot the results
  # the 0 one doesn't get a refrow with lnegth.out set to 300
  if(i == 0 | i == 5){
    nd <- expand.grid(
      StdTmeanAnomalyRS=seq(from = min(model_output$data$StdTmeanAnomalyRS),
                            to = max(model_output$data$StdTmeanAnomalyRS),
                            length.out = 200),
      UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
                 levels = levels(model_output$data$UI2)))
  }else{
    
    nd <- expand.grid(
      StdTmeanAnomalyRS=seq(from = min(model_output$data$StdTmeanAnomalyRS),
                            to = max(model_output$data$StdTmeanAnomalyRS),
                            length.out = 300),
      UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
                 levels = levels(model_output$data$UI2)))
    
  }
  
  # back transform the predictors
  nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
    transformedX = nd$StdTmeanAnomalyRS,
    originalX = predicts_sites2$StdTmeanAnomaly)
  
  # set richness and abundance to 0 - to be predicted
  nd$LogAbund <- 0
  nd$Species_richness <- 0
  
  # reference for % difference = primary vegetation and positive anomaly closest to 0
  refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))
  
  # adjust plot 1: mean anomaly and abundance
  exclQuantiles <- c(0.025,0.975)
  
  QPV <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Primary vegetation"],
    probs = exclQuantiles)
  QSV <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Secondary vegetation"],
    probs = exclQuantiles)
  QAL <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Agriculture_Low"],
    probs = exclQuantiles)
  QAH <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Agriculture_High"],
    probs = exclQuantiles)
  
  # predict the results
  a.preds.tmean <- PredictGLMERRandIter(model = model_output$model,data = nd)
  
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
  
  nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))
  
  # plot
  ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
    geom_line(aes(col = UI2), size = 1) +
    geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
    scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
    theme_bw() + 
    ylab("Change in Total Abundance (%)") +
    xlab("Standardised Temperature Anomaly") +
    xlim(c(-0.5, 2)) +
    ylim(c(-100, 150)) + 
    theme(aspect.ratio = 1, text = element_text(size = 12),
          strip.text.x = element_text(hjust = 0, size = 12, face = "bold"), legend.title = element_blank()) +
    ggtitle(paste0("Sites with more than ", i, " stns"))
  
  # save
  ggsave(filename = paste0(outdir, "/Plot_removing_", i, "_Abun_Tropical.pdf"), height = 5, width = 6)
  
  
  
  rm(nd, model_data, model_output)
  
}




#### Tropical sites, Richness ####

table(trop_data$nstn)
#   0   1   2   3   4   5   6   7   8 
# 148 178 459 604  39  65   6  29 214


# i <- 0
for(i in vals){
  
  print(i)
  
  model_data <- trop_data[trop_data$nstn > i, ]
  
  model_data <- droplevels(model_data)

  
  model_output <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                        fixedStruct = "UI2 * StdTmeanAnomalyRS",
                        randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)")
  
  
  # save model output
  save(model_output, file = paste0(outdir, "/Model_output_removing_", i, "_Tropical_Rich.rdata"))
  
  # save model stats
  write.csv(model_output$stats, file = paste0(outdir, "/Model_stats_removing_", i, "_Tropical_Rich.csv"))
  
  
  # plot the results
  # the 0 one doesn't get a refrow with lnegth.out set to 300
  if(i == 0){
    nd <- expand.grid(
      StdTmeanAnomalyRS=seq(from = min(model_output$data$StdTmeanAnomalyRS),
                            to = max(model_output$data$StdTmeanAnomalyRS),
                            length.out = 175),
      UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
                 levels = levels(model_output$data$UI2)))
  }else{
    
    nd <- expand.grid(
      StdTmeanAnomalyRS=seq(from = min(model_output$data$StdTmeanAnomalyRS),
                            to = max(model_output$data$StdTmeanAnomalyRS),
                            length.out = 300),
      UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
                 levels = levels(model_output$data$UI2)))
    
  }
  
  # back transform the predictors
  nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
    transformedX = nd$StdTmeanAnomalyRS,
    originalX = predicts_sites2$StdTmeanAnomaly)
  
  # set richness and abundance to 0 - to be predicted
  nd$LogAbund <- 0
  nd$Species_richness <- 0
  
  # reference for % difference = primary vegetation and positive anomaly closest to 0
  refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))
  
  # adjust plot 1: mean anomaly and abundance
  exclQuantiles <- c(0.025,0.975)
  
  QPV <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Primary vegetation"],
    probs = exclQuantiles)
  QSV <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Secondary vegetation"],
    probs = exclQuantiles)
  QAL <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Agriculture_Low"],
    probs = exclQuantiles)
  QAH <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Agriculture_High"],
    probs = exclQuantiles)
  
  # predict the results
  a.preds.tmean <- PredictGLMERRandIter(model = model_output$model,data = nd)
  
  # back transform the abundance values
  a.preds.tmean <- exp(a.preds.tmean)
  
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
  
  nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))
  
  
  if(i %in% c(3,4,5,6)){
    
    # plot
    ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
      geom_line(aes(col = UI2), size = 1) +
      geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
      geom_hline(yintercept = 0, linetype = "dashed") +
      scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
      scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
      theme_bw() + 
      ylab("Change in Species Richness (%)") +
      xlab("Standardised Temperature Anomaly") +
      xlim(c(-0.5, 2)) +
      ylim(c(-100, 450)) + 
      theme(aspect.ratio = 1, text = element_text(size = 12),
            strip.text.x = element_text(hjust = 0, size = 12, face = "bold"), legend.title = element_blank()) +
      ggtitle(paste0("Sites with more than ", i, " stns"))
    
    # save
    ggsave(filename = paste0(outdir, "/Plot_removing_", i, "_Tropical_Rich.pdf"), height = 5, width = 6)
    
    
    
    
  }else{
  
  # plot
  ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
    geom_line(aes(col = UI2), size = 1) +
    geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
    scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
    theme_bw() + 
    ylab("Change in Species Richness (%)") +
    xlab("Standardised Temperature Anomaly") +
    xlim(c(-0.5, 2)) +
    ylim(c(-100, 150)) + 
    theme(aspect.ratio = 1, text = element_text(size = 12),
          strip.text.x = element_text(hjust = 0, size = 12, face = "bold"), legend.title = element_blank()) +
    ggtitle(paste0("Sites with more than ", i, " stns"))
  
  # save
  ggsave(filename = paste0(outdir, "/Plot_removing_", i, "_Tropical_Rich.pdf"), height = 5, width = 6)
  
  
  }
  
  
  rm(nd, model_data, model_output)
  
}



#### Temperate sites, Abundance ####

vals <- c(0,1,2,3,4,5,6)


temp_data <- predicts_sites2[!predicts_sites2$SS %in% trop_data$SS, ] # 4327 rows

table(temp_data[!is.na(temp_data$LogAbund), "nstn"])
#  2    3    4    5    6    7    8 
# 46  146   49    1   15   37 3852 


# i <- 0
for(i in vals){
  
  model_data <- temp_data[temp_data$nstn > i, ]
  
  model_data <- model_data[!is.na(model_data$LogAbund), ]
  
  model_data <- droplevels(model_data)
  

  model_output <- GLMER(modelData = model_data,responseVar = "LogAbund", 
                        fitFamily = "gaussian",
                        fixedStruct = "UI2 * StdTmeanAnomalyRS",
                        randomStruct = "(1|SS)+(1|SSB)")
  # save model output
  save(model_output, file = paste0(outdir, "/Model_output_removing_", i, "_Abun_Temperate.rdata"))
  
  # save model stats
  write.csv(model_output$stats, file = paste0(outdir, "/Model_stats_removing_", i, "_Abun_Temperate.csv"))
  
  
  # plot the results
  # the 0 one doesn't get a refrow with lnegth.out set to 300
  #if(i == 0 | i == 5){
    nd <- expand.grid(
      StdTmeanAnomalyRS=seq(from = min(model_output$data$StdTmeanAnomalyRS),
                            to = max(model_output$data$StdTmeanAnomalyRS),
                            length.out = 200),
      UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
                 levels = levels(model_output$data$UI2)))
  # }else{
  #   
  #   nd <- expand.grid(
  #     StdTmeanAnomalyRS=seq(from = min(model_output$data$StdTmeanAnomalyRS),
  #                           to = max(model_output$data$StdTmeanAnomalyRS),
  #                           length.out = 300),
  #     UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
  #                levels = levels(model_output$data$UI2)))
  #   
  # }
  
  # back transform the predictors
  nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
    transformedX = nd$StdTmeanAnomalyRS,
    originalX = predicts_sites2$StdTmeanAnomaly)
  
  # set richness and abundance to 0 - to be predicted
  nd$LogAbund <- 0
  nd$Species_richness <- 0
  
  # reference for % difference = primary vegetation and positive anomaly closest to 0
  refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))
  
  # adjust plot 1: mean anomaly and abundance
  exclQuantiles <- c(0.025,0.975)
  
  QPV <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Primary vegetation"],
    probs = exclQuantiles)
  QSV <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Secondary vegetation"],
    probs = exclQuantiles)
  QAL <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Agriculture_Low"],
    probs = exclQuantiles)
  QAH <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Agriculture_High"],
    probs = exclQuantiles)
  
  # predict the results
  a.preds.tmean <- PredictGLMERRandIter(model = model_output$model,data = nd, nIters = 2000)
  
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
  
  nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))
  
  # plot
  ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
    geom_line(aes(col = UI2), size = 1) +
    geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
    scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
    theme_bw() + 
    ylab("Change in Total Abundance (%)") +
    xlab("Standardised Temperature Anomaly") +
    xlim(c(-0.5, 2)) +
    ylim(c(-100, 150)) + 
    theme(aspect.ratio = 1, text = element_text(size = 12),
          strip.text.x = element_text(hjust = 0, size = 12, face = "bold"), legend.title = element_blank()) +
    ggtitle(paste0("Sites with more than ", i, " stns"))
  
  # save
  ggsave(filename = paste0(outdir, "/Plot_removing_", i, "_Abun_Temperate.pdf"), height = 5, width = 6)
  
  
  
  rm(nd, model_data, model_output)
  
}




#### Temperate sites, Richness ####

table(temp_data$nstn)
#  2    3    4    5    6    7    8 
# 46  146  108    1   15   37 3974 


# i <- 0
for(i in vals){
  
  print(i)
  
  model_data <- temp_data[temp_data$nstn > i, ]
  
  model_data <- droplevels(model_data)
  

  model_output <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
                        fixedStruct = "UI2 * StdTmeanAnomalyRS",
                        randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)")
  
  
  # save model output
  save(model_output, file = paste0(outdir, "/Model_output_removing_", i, "_Temperate_Rich.rdata"))
  
  # save model stats
  #write.csv(model_output$stats, file = paste0(outdir, "/Model_stats_removing_", i, "_Temperate_Rich.csv"))
  
  
  # plot the results
  # the 0 one doesn't get a refrow with lnegth.out set to 300
  #if(i == 0){
    nd <- expand.grid(
      StdTmeanAnomalyRS=seq(from = min(model_output$data$StdTmeanAnomalyRS),
                            to = max(model_output$data$StdTmeanAnomalyRS),
                            length.out = 275),
      UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
                 levels = levels(model_output$data$UI2)))
  # }else{
  #   
  #   nd <- expand.grid(
  #     StdTmeanAnomalyRS=seq(from = min(model_output$data$StdTmeanAnomalyRS),
  #                           to = max(model_output$data$StdTmeanAnomalyRS),
  #                           length.out = 300),
  #     UI2=factor(c("Primary vegetation","Secondary vegetation","Agriculture_Low","Agriculture_High"),
  #                levels = levels(model_output$data$UI2)))
  #   
  # }
  # 
  # back transform the predictors
  nd$StdTmeanAnomaly <- BackTransformCentreredPredictor(
    transformedX = nd$StdTmeanAnomalyRS,
    originalX = predicts_sites2$StdTmeanAnomaly)
  
  # set richness and abundance to 0 - to be predicted
  nd$LogAbund <- 0
  nd$Species_richness <- 0
  
  # reference for % difference = primary vegetation and positive anomaly closest to 0
  refRow <- which((nd$UI2=="Primary vegetation") & (nd$StdTmeanAnomaly==min(abs(nd$StdTmeanAnomaly))))
  

  # adjust plot 1: mean anomaly and abundance
  exclQuantiles <- c(0.025,0.975)
  
  QPV <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Primary vegetation"],
    probs = exclQuantiles)
  QSV <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Secondary vegetation"],
    probs = exclQuantiles)
  QAL <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Agriculture_Low"],
    probs = exclQuantiles)
  QAH <- quantile(x = model_output$data$StdTmeanAnomalyRS[
    model_output$data$UI2=="Agriculture_High"],
    probs = exclQuantiles)
  
  # predict the results
  a.preds.tmean <- PredictGLMERRandIter(model = model_output$model,data = nd)
  
  # back transform the abundance values
  a.preds.tmean <- exp(a.preds.tmean)
  
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
  
  nd$UI2 <- factor(nd$UI2, levels = c("Primary vegetation", "Secondary vegetation", "Agriculture_Low", "Agriculture_High" ))
  
  # plot
  ggplot(data = nd, aes(x = StdTmeanAnomaly, y = PredMedian)) + 
    geom_line(aes(col = UI2), size = 1) +
    geom_ribbon(aes(ymin = nd$PredLower, ymax = nd$PredUpper, fill = UI2), alpha = 0.2) +
    geom_hline(yintercept = 0, linetype = "dashed") +
    scale_fill_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
    scale_colour_manual(values = c("#009E73", "#0072B2", "#E69F00", "#D55E00")) +
    theme_bw() + 
    ylab("Change in Species Richness (%)") +
    xlab("Standardised Temperature Anomaly") +
    xlim(c(-0.5, 2)) +
    ylim(c(-100, 200)) + 
    theme(aspect.ratio = 1, text = element_text(size = 12),
          strip.text.x = element_text(hjust = 0, size = 12, face = "bold"), legend.title = element_blank()) +
    ggtitle(paste0("Sites with more than ", i, " stns"))
  
  # save
  ggsave(filename = paste0(outdir, "/Plot_removing_", i, "_Temperate_Rich.pdf"), height = 5, width = 6)
  
  
  
  rm(nd, model_data, model_output)
  
}







#### maps of stns ####

pdf(file = paste0(outdir, "/Maps_removing_sites.pdf"), width = 10, height = 6)
par(mfrow = c(2,3))

predicts_sp_1 <- predicts_sp[predicts_sp$stn > 0, ]

plot(stn_0130_mean, main = "Sites supported by 1 or more stations")
plot(predicts_sp_1, add = T)
abline(h = -23.44, lty = "dashed")
abline(h = 23.44, lty = "dashed")

predicts_sp_2 <- predicts_sp[predicts_sp$stn > 1, ]

plot(stn_0130_mean, main = "Sites supported by 2 or more stations")
plot(predicts_sp_2, add = T)
abline(h = -23.44, lty = "dashed")
abline(h = 23.44, lty = "dashed")

predicts_sp_3 <- predicts_sp[predicts_sp$stn > 2, ]

plot(stn_0130_mean, main = "Sites supported by 3 or more stations")
plot(predicts_sp_3, add = T)
abline(h = -23.44, lty = "dashed")
abline(h = 23.44, lty = "dashed")

predicts_sp_4 <- predicts_sp[predicts_sp$stn > 3, ]

plot(stn_0130_mean, main = "Sites supported by 4 or more stations")
plot(predicts_sp_4, add = T)
abline(h = -23.44, lty = "dashed")
abline(h = 23.44, lty = "dashed")

predicts_sp_5 <- predicts_sp[predicts_sp$stn > 4, ]

plot(stn_0130_mean, main = "Sites supported by 5 or more stations")
plot(predicts_sp_5, add = T)
abline(h = -23.44, lty = "dashed")
abline(h = 23.44, lty = "dashed")

predicts_sp_6 <- predicts_sp[predicts_sp$stn > 5, ]

plot(stn_0130_mean, main = "Sites supported by 6 or more stations")
plot(predicts_sp_6, add = T)
abline(h = -23.44, lty = "dashed")
abline(h = 23.44, lty = "dashed")

dev.off()


#### histograms of data ####

pdf(file = paste0(outdir, "/Histograms_anomaly_removing_sites.pdf"), width = 10, height = 8)
par(mfrow = c(2,2))


hist(predicts_sp$StdTmeanAnomaly, main = "All sites", xlab = "Standardised Temperature Anomaly")
hist(predicts_sp_1$StdTmeanAnomaly, main = "Sites supported by 1 or more stations", xlab = "Standardised Temperature Anomaly")
hist(predicts_sp_2$StdTmeanAnomaly, main = "Sites supported by 2 or more stations", xlab = "Standardised Temperature Anomaly")
hist(predicts_sp_3$StdTmeanAnomaly, main = "Sites supported by 3 or more stations", xlab = "Standardised Temperature Anomaly")


dev.off()


pdf(file = paste0(outdir, "/Histograms_anomaly_removing_sites_tropical.pdf"), width = 10, height = 8)
par(mfrow = c(2,3))

#hist(trop_data$StdTmeanAnomaly, main = "All sites", xlab = "Standardised Temperature Anomaly")
hist(trop_data[trop_data$nstn > 0, "StdTmeanAnomaly"], main = "Sites supported by 1 or more stations", xlab = "Standardised Temperature Anomaly")
hist(trop_data[trop_data$nstn > 1, "StdTmeanAnomaly"], main = "Sites supported by 2 or more stations", xlab = "Standardised Temperature Anomaly")
hist(trop_data[trop_data$nstn > 2, "StdTmeanAnomaly"], main = "Sites supported by 3 or more stations", xlab = "Standardised Temperature Anomaly")
hist(trop_data[trop_data$nstn > 3, "StdTmeanAnomaly"], main = "Sites supported by 4 or more stations", xlab = "Standardised Temperature Anomaly")
hist(trop_data[trop_data$nstn > 4, "StdTmeanAnomaly"], main = "Sites supported by 5 or more stations", xlab = "Standardised Temperature Anomaly")
hist(trop_data[trop_data$nstn >5, "StdTmeanAnomaly"], main = "Sites supported by 6 or more stations", xlab = "Standardised Temperature Anomaly")


dev.off()

t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()