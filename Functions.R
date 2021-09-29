###Functions
library(StatisticalModels)
library(predictsFunctions)
library(dplyr)


##Function for standardizing and centering predictor variables on very different scales, e.g. percNH
StdCenterPredictor <- function(x) {
  variable <- x
  sd <- sd(na.omit(variable))
  mean <- mean(na.omit(variable))
  variable.s <- (variable - mean)/sd
  return(variable.s)
}

BackTransformCentreredPredictor <- function(transformedX,originalX){
  
  sd <- sd(na.omit(originalX))
  mean <- mean(na.omit(originalX))
  backTX <- (transformedX * sd) + mean
  
  return(backTX)
  
}

RescaleAbundance <- function(sites){
  
  sites <- droplevels(sites)
  
  StudyMaxAbund <- suppressWarnings(tapply(
    X = sites$Total_abundance,INDEX = sites$SS,
    FUN = max,na.rm=TRUE))
  
  StudyMaxAbund[StudyMaxAbund == -Inf] <- NA
  
  AllMaxAbund <- StudyMaxAbund[match(sites$SS,names(StudyMaxAbund))]
 
  sites$Total_abundance_RS <- sites$Total_abundance/AllMaxAbund
  
  sites$LogAbund <- log(sites$Total_abundance_RS+0.01)
   
  return(sites)
}

##For a set of merged sites in predicts:
##This functions finds climate data, calculates scaled abundance values, logs abundance values,
##organizes the LandUse class, calculate standardise climate anomalies, and finds the percNH data
##In order to work, you must run the climate_percNH.R script and use the csv from that as the
## predicts_sites_filename arguement

organize_land_use_find_climate <- function(aggregated_sites) {
  
  sites <- aggregated_sites
  predicts_sites <- readRDS("Outputs/predicts_sites_info.rds")
  
  print("Organising Land Use")
  #####Organizing Land use
  sites$LandUse <- paste(sites$Predominant_land_use)
  sites$LandUse[(sites$LandUse=="Cannot decide")] <- NA
  sites$LandUse[(sites$LandUse=="Secondary vegetation (indeterminate age)")] <- NA
  sites$LandUse[(sites$LandUse=="Urban")] <- NA
  sites$LandUse <- factor(sites$LandUse)
  sites$LandUse <- relevel(sites$LandUse,ref="Primary vegetation")
  sites$Use_intensity[sites$Use_intensity=="Cannot decide"] <- NA
  sites$LandUse <- relevel(sites$LandUse, ref="Primary vegetation")
  sites$Use_intensity <- relevel(sites$Use_intensity, ref="Minimal use")
  
  ###Create agriculture intensity land class
  
  sites$Land <- paste(sites$LandUse)
  sites$Use_intensity <- as.character(sites$Use_intensity)
  sites$Use_intensity[is.na(sites$Use_intensity)] <- "Unknown"
  sites$Land <- interaction(sites$Land, sites$Use_intensity)
  sites$Land <- as.character(sites$Land)
  sites$Land[sites$Land=="Cropland.Intense use"] <- "High Agriculture"
  sites$Land[sites$Land=="Cropland.Light use"] <- "Low Agriculture"
  sites$Land[sites$Land=="Cropland.Minimal use"] <- "Low Agriculture"
  sites$Land[sites$Land=="Cropland.Unknown"] <- NA
  sites$Land[sites$Land=="Plantation forest.Light use"] <- "High Agriculture"
  sites$Land[sites$Land=="Plantation forest.Intense use"] <- "High Agriculture"
  sites$Land[sites$Land=="Plantation forest.Minimal use"] <- "Low Agriculture"
  sites$Land[sites$Land=="Plantation forest.Unknown"] <- NA
  sites$Land[sites$Land=="Pasture.Light use"] <- "Low Agriculture"
  sites$Land[sites$Land=="Pasture.Minimal use"] <- "Low Agriculture"
  sites$Land[sites$Land=="Pasture.Intense use"] <- "High Agriculture"
  sites$Land[sites$Land=="Pasture.Unknown"] <- NA
  sites$Land[sites$Land=="Primary vegetation.Intense use"] <- "PV"
  sites$Land[sites$Land=="Primary vegetation.Light use"] <- "PV"
  sites$Land[sites$Land=="Primary vegetation.Minimal use"] <- "PV"
  sites$Land[sites$Land=="Primary vegetation.Unknown"] <- "PV"
  sites$Land[sites$Land=="Intermediate secondary vegetation.Intense use"] <- "SV"
  sites$Land[sites$Land=="Intermediate secondary vegetation.Light use"] <- "SV"
  sites$Land[sites$Land=="Intermediate secondary vegetation.Minimal use"] <- "SV"
  sites$Land[sites$Land=="Intermediate secondary vegetation.Unknown"] <- "SV"
  sites$Land[sites$Land=="Mature secondary vegetation.Intense use"] <- "SV"
  sites$Land[sites$Land=="Mature secondary vegetation.Light use"] <- "SV"
  sites$Land[sites$Land=="Mature secondary vegetation.Minimal use"] <- "SV"
  sites$Land[sites$Land=="Mature secondary vegetation.Unknown"] <- "SV"
  sites$Land[sites$Land=="Young secondary vegetation.Intense use"] <- "SV"
  sites$Land[sites$Land=="Young secondary vegetation.Unknown"] <- "SV"
  sites$Land[sites$Land=="Young secondary vegetation.Minimal use"] <- "SV"
  sites$Land[sites$Land=="Young secondary vegetation.Light use"] <- "SV"
  sites$Land <- factor(sites$Land, levels = c("PV", "SV", "Low Agriculture", "High Agriculture"))
  
  print(table(sites$Land))
  
  ##Calculate log abundance
  sites$log_abundance <- log(sites$Total_abundance+1)
  
  
  ##Calculate Scaled abundance values
  n_sites <- length(unique(sites$SSBS))
  n_studies <- length(unique(sites$SS))
  print(paste("Calculating Scaled Abundance for",n_sites , "sites across",n_studies , "studies"))
  sites$Scaled_Abundance <- rep(3,1)
  for(i in 1:length(unique(sites$SS))) {
    s <- unique(sites$SS)[i]
    g <- sites[sites$SS==s, ]
    max <- max(g$Total_abundance)
    l <- g$Total_abundance / max
    
    sites[sites$SS==s, ]$Scaled_Abundance <- l
  }
  not_abundance <- which(sites$Scaled_Abundance==3)
  if(length(not_abundance)>0) sites[not_abundance, ]$Scaled_Abundance <- NA
  sites$log_scaled_abundance <- log(sites$Scaled_Abundance+0.01)
  
  ###Finding climate data for sites
  
  ##One Needs to be a character vector to match a factor
  predicts_sites$SSBS <- as.character(predicts_sites$SSBS)
  print("Finding Climate Data")
  site_info <- dplyr::select(predicts_sites, SSBS,climate_anomaly:SV_10000)
  sites <- left_join(sites, site_info)
  
  ##Combine PV and SV
  sites$NH_1000 <- sites$PV_1000 + sites$SV_1000
  sites$NH_3000 <- sites$PV_3000 + sites$SV_3000
  sites$NH_5000 <- sites$PV_5000 + sites$SV_5000
  sites$NH_10000 <- sites$PV_10000 + sites$SV_10000
  
  print("Standardising Climate Data")
  
  
  sites$std_climate_anomaly <- sites$climate_anomaly / sites$historic_sd
  sites$std_tmax_anomaly <- sites$tmax_anomaly / sites$historic_sd_tmax
  sites$std_tmax_quarter_anomaly <- sites$tmax_quarter_anomaly / sites$historic_sd_tmax
  
  ##Standardising and centering variables on very different scales
  sites$NH_1000.s <- StdCenterPredictor(sites$NH_1000)
  sites$NH_3000.s <- StdCenterPredictor(sites$NH_3000)
  sites$NH_5000.s <- StdCenterPredictor(sites$NH_5000)
  sites$NH_10000.s <- StdCenterPredictor(sites$NH_10000)
  
  print("Done")
  
  return(sites)
}
max_quarter_fast <- function(x) {
  
  mean(Rfast::Sort(x, TRUE)[1:3])
  
}



vif.mer <- function (fit) {
  ## adapted from rms::vif
  
  v <- vcov(fit)
  nam <- names(fixef(fit))
  
  ## exclude intercepts
  ns <- sum(1 * (nam == "Intercept" | nam == "(Intercept)"))
  if (ns > 0) {
    v <- v[-(1:ns), -(1:ns), drop = FALSE]
    nam <- nam[-(1:ns)]
  }
  
  d <- diag(v)^0.5
  v <- diag(solve(v/(d %o% d)))
  names(v) <- nam
  v
}


# edited version of StatisticalModels function for use with glmmTMB model

PredictGLMERRandIter_TMB <- function(model,data,nIters=1000){
  
  # stopifnot((class(model)[1] == "lmerMod") | (class(model)[1] == "glmerMod"))
  
  preds <- sapply(X = 1:nIters,FUN = function(i){
    
    mm<-model.matrix(terms(model),data)
    
    if(ncol(mm)>length(lme4::fixef(model)[1]$cond)){
      mm <- mm[,-which(!(names(mm[1,]) %in% names(
        lme4::fixef(model))))]
    }
    
    fe.draw <- mvrnorm(n = 1,mu = fixef(model)[1]$cond,Sigma = vcov(model)[1]$cond)
    
    y <- mm %*% fe.draw
    
  })
  
  return(preds)
}



