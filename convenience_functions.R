###Functions
library(StatisticalModels)
library(predictsFunctions)

##Function for standardizing and centering predictor variables on very different scales, e.g. percNH
StdCenterPredictor <- function(x) {
  variable <- x
  sd <- sd(na.omit(variable))
  mean <- mean(na.omit(variable))
  variable.s <- (variable - mean)/sd
  return(variable.s)
}

##For a set of merged sites in predicts:
##This functions finds climate data, calculates scaled abundance values, logs abundance values,
##organizes the LandUse class, calculate standardise climate anomalies, and finds the percNH data
##In order to work, you must run the climate_percNH.R script and use the csv from that as the
## predicts_sites_filename arguement
organize_land_use_find_climate <- function(aggregated_sites, predicts_sites_filename) {
  
  sites <- aggregated_sites
  predicts_sites <- read.csv(predicts_sites_filename)
  
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
  for(i in 1:nrow(sites)) {
    
    sites$avg_temp[i] <-  predicts_sites$avg_temp[predicts_sites$SSBS==sites$SSBS[i]]
    sites$tmax[i] <-  predicts_sites$tmax[predicts_sites$SSBS==sites$SSBS[i]]
    sites$climate_anomaly[i] <-  predicts_sites$rate[predicts_sites$SSBS==sites$SSBS[i]]
    sites$avg_temp_sd[i] <-  predicts_sites$avg_temp_sd[predicts_sites$SSBS==sites$SSBS[i]]
    sites$historic_sd[i] <-  predicts_sites$historic_sd[predicts_sites$SSBS==sites$SSBS[i]]
    sites$tmax_anomaly[i] <-  predicts_sites$tmax_anomaly[predicts_sites$SSBS==sites$SSBS[i]]
    sites$percNH[i] <-  predicts_sites$percNH[predicts_sites$SSBS==sites$SSBS[i]]
  }
  
  
  sites$tropics <- ifelse(abs(sites$Latitude<23.5), "tropics", "temperate")
  print("Standardising Climate Data")
  sites$std_climate_anomaly <- sites$climate_anomaly / sites$historic_sd
  sites$std_tmax_anomaly <- sites$tmax_anomaly / sites$historic_sd
  
  ##Standardising and centering variables on very different scales
  sites$percNH.s <-StdCenterPredictor(sites$percNH)
  print("Done")
  
  return(sites)
}



##This function requires predictsFunctions in order to work.
create_sites <- function(predicts_data) {
  format <- predicts_data
  #Correct Sampling effort
  format <- CorrectSamplingEffort(diversity = format)
  
  # Merge sites that have the same coordinates (e.g. multiple traps on a single transect)
  format <- MergeSites(diversity = format)
  
  # Calculate site metrics of diversity - currently only species richness and total abundance
  format <- SiteMetrics(diversity = format,
                        extra.cols = c("Predominant_land_use",
                                       "SSB","SSBS","Km_to_nearest_edge_of_habitat", "Ecoregion","Biome", "Habitat_patch_area_square_metres"))
  return(format)
}


