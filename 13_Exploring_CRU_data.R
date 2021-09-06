##%######################################################%##
#                                                          #
####        13. Exploring data bias in CRU data         ####
#                                                          #
##%######################################################%##

# In this script, I will take a look at the data coverage of the CRU data
# using the stn variable within the mean temp data

rm(list = ls())

# load libraries
library(raster)
library(predictsFunctions)
library(StatisticalModels)
source("Functions.R")


# organise directories
datadir <- "0_data"
outdir <- "13_Exploring_CRU_data"
dir.create(outdir)

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

# round to whole numbers  
#nstn <- round(nstn, digits = 0)
# take a look at spread
#table(nstn)

#   0    1    2    3    4    5    6    7    8 
# 148  178  565  690  147   66   21   66 4195 


### look at just those within tropics ###

# add values to predicts data
predicts_sites2$nstn <- nstn

trop <- predicts_sites2[predicts_sites2$Latitude > -23.44 & predicts_sites2$Latitude < 23.44, ] # 1742 rows


table(round(trop$nstn, digits = 0))

#   0   1   2   3   4   5   6   7   8 
# 148 178 519 544  39  65   6  29 214 



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


# i <- 0
for(i in vals){
  
  model_data <- predicts_sites2[predicts_sites2$nstn > i, ]
  
  model_data <- droplevels(model_data)
  
  model_output <- GLMERSelect(modelData = model_data,responseVar = "LogAbund",
                              fitFamily = "gaussian",fixedFactors = "UI2",
                              fixedTerms = list(StdTmeanAnomalyRS=1),
                              randomStruct = "(1|SS)+(1|SSB)",
                              fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                              saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))
  
  #summary(model_output$model)
  
  
  # save model output
  save(model_output, file = paste0(outdir, "/Model_output_removing_", i, ".rdata"))
  
  # save model stats
  write.csv(model_output$stats, file = paste0(outdir, "/Model_stats_removing_", i, ".csv"))
  
  
  # plot the results
  # the 0 one doesn't get a refrow with lnegth.out set to 300
  if(i == 0){
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
    originalX = model_data$StdTmeanAnomaly)
  
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
  ggsave(filename = paste0(outdir, "/Plot_removing_", i, ".pdf"), height = 5, width = 6)
  

  
 rm(nd, model_data, model_output)
  
}


