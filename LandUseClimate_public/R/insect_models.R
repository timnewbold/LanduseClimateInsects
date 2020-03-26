
library(StatisticalModels)
library(predictsFunctions)
source("r/functions.r")

###Create Models for all insects in predicts for standardised climate anomaly and Land interactions

predicts <- readRDS("Data/database.rds")
insect_predicts <- droplevels(subset(predicts, predicts$Class=="Insecta"))

##Create sites
insect_sites <- create_sites(insect_predicts)
insect_sites <- organize_land_use_find_climate(insect_sites)


#Drop places with no land or climate data
insect_sites <- insect_sites[!is.na(insect_sites$Land), ]
insect_sites <- insect_sites[!is.na(insect_sites$std_climate_anomaly), ]

#Center std climate anomaly
insect_sites$std_climate_anomaly.s <- StdCenterPredictor(insect_sites$std_climate_anomaly)
insect_sites$std_tmax_anomaly.s <- StdCenterPredictor(insect_sites$std_tmax_anomaly)


climate <- GLMER(insect_sites, 
                 responseVar = "log_scaled_abundance", fitFamily = "gaussian",
                 fixedStruct = "Land*std_climate_anomaly.s",
                 randomStruct = "(1|SS) + (1|SSB)", saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))
summary(climate$model)
vif.mer(climate$model)

cols <- c("#003300", "#00cc99", "#ffcc00", "#800000")
PlotGLMERContinuous(model = climate$model,data = climate$data,
                    effects = "std_climate_anomaly.s", byFactor = "Land",
                    xlab = "Standardised centered climate anomaly",ylab = "Scaled Abundance",
                    logLink = "e", ylim = c(0,0.5), main="All insects Land*std_climate_anomaly",
                    line.cols =cols, plotRug = T,seMultiplier = 1)



tmax <- GLMER(insect_sites, 
              responseVar = "log_scaled_abundance", fitFamily = "gaussian",
              fixedStruct = "Land*std_tmax_anomaly.s",
              randomStruct = "(1|SS) + (1|SSB)", saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))
summary(tmax$model)
vif.mer(tmax$model)
PlotGLMERContinuous(model = tmax$model,data = tmax$data,
                    effects = "std_tmax_anomaly.s", byFactor = "Land",
                    xlab = "Standardised centered tmax anomaly",ylab = "Scaled Abundance",
                    logLink = "e", main ="All insects Land*std_tmax_anomaly", ylim = c(0,0.5),
                    line.cols =cols, plotRug = TRUE,seMultiplier = 1)
legend(-2, 0.2, legend = levels(climate$data$Land), col = cols, pch=19, cex=1.4)


climate_rich <- GLMER(insect_sites, 
                      responseVar = "Species_richness", fitFamily = "poisson",
                      fixedStruct = "Land*std_climate_anomaly.s",
                      randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))
summary(climate_rich$model)
vif.mer(climate_rich$model)
par(mfrow=c(1,1))
PlotGLMERContinuous(model = climate_rich$model,data = climate_rich$data,
                    effects = "std_climate_anomaly.s", byFactor = "Land",
                    xlab = "Standardised centered climate anomaly",ylab = "Species Richness",
                    logLink = "e", main ="All insects Land*std_climate_anomaly", ylim = c(0,20),
                    line.cols =cols, plotRug = TRUE,seMultiplier = 1)


legend(-2, 6, legend = levels(climate$data$Land), col = cols, pch=19, cex=1.4)

tmax_rich <- GLMER(insect_sites, 
                   responseVar = "Species_richness", fitFamily = "poisson",
                   fixedStruct = "Land*std_tmax_anomaly.s",
                   randomStruct = "(1|SS) + (1|SSB)+ (1|SSBS)", saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))
summary(tmax_rich$model)
vif.mer(tmax_rich$model)
PlotGLMERContinuous(model = tmax_rich$model,data = tmax_rich$data,
                    effects = "std_tmax_anomaly.s", byFactor = "Land",
                    xlab = "Standardised centered tmax anomaly",ylab = "Species richness",
                    logLink = "e", main ="All insects Land*std_tmax_anomaly", ylim = c(0,20),
                    line.cols =cols, plotRug = TRUE,seMultiplier = 1)
legend(-2,6,  legend = levels(climate$data$Land), col = cols, pch=19, cex=1.4)

