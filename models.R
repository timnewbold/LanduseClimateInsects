#Load in required packages
library(predictsFunctions)
library(StatisticalModels)

predicts.path <-"//ad.ucl.ac.uk/homec/zcbtpmc/Documents/Predicts/R/database.rds"
predicts <- readRDS(predicts.path)


#Subset predicts to include only insects and drop unused factor levels
insects_predicts <- droplevels(subset(predicts, predicts$Class=="Insecta"))

#Correct sampling effort, merge sites and calculate Total abundance and species richness
insect_sites <- create_sites(insects_predicts)


##This should be the file location for output of the climate_percNH script
predicts.climate.data.path <- "//ad.ucl.ac.uk/homec/zcbtpmc/Documents/Predicts/pollinator_taxonomic_geographic_dist_text-analysis-master/pollinator_taxonomic_geographic_dist_text-analysis-master/data_my_project/predicts_sites_6.csv"


#Organise Land use, calculate scaled abundance, find climate data for insect sites
insect_sites <- organize_land_use_find_climate(insect_sites, predicts.climate.data.path)


#Reference for colnames
colnames(insect_sites)

#LandUse is Landcover excluding urban site
#log_abundance is the log(x+1) transformation of Total_abundance
#Scaled_abundance scales abundance within studies in predicts, giving the maximum abundance in a study a value of 1
#and all other abundances with that study relative to that value
#log_scaled_abundance is the log(x+0.01) transformation of Scaled_abundance
#avg_temp is the Mean annual temp of the site in year preceding the end sample date
#tmax is the mean maximum temperature of the warmest month in the year preceding the end sample date
#avg_temp_sd is the standard deviation of mean monthly temperatures in the year preceding the end sample date
#historic_sd is the standard deviation of mean monthly temperatures between 1901 and 1905
#climate_anomaly is avg_temp - the mean temperature for all months between 1901 and 1905
#tmax_anomaly is tmax - mean of the hottest month of each year from 1901 to 1905
#percNH is the percentage natural habitat in the by 5km grid cell in which the site is found
#percNH.s is the rescaled and centered version of percNH
#std_climate anomaly is climate_anomaly / historic_sd
#std_tmax_anomaly is tmax_anomaly / historic_sd



###Main models additive and interactive

mod_add <- GLMER(insect_sites, responseVar = "log_scaled_abundance", fitFamily = "gaussian",
                 fixedStruct = "LandUse*std_climate_anomaly+LandUse*percNH.s",
                 randomStruct = "(1|SS) + (1|SSB)")


mod_int <- GLMER(insect_sites, responseVar = "log_scaled_abundance", fitFamily = "gaussian",
                 fixedStruct = "LandUse*std_climate_anomaly*percNH.s",
                 randomStruct = "(1|SS) + (1|SSB)")



mod_no_percHN <- GLMER(insect_sites, responseVar = "log_scaled_abundance", fitFamily = "gaussian",
                       fixedStruct = "LandUse*std_climate_anomaly",
                       randomStruct = "(1|SS) + (1|SSB)")

mod_just_landuse <- GLMER(insect_sites, responseVar = "log_scaled_abundance", fitFamily = "gaussian",
                       fixedStruct = "LandUse",
                       randomStruct = "(1|SS) + (1|SSB)")




###Refitting models to compare AIC so that they are all fit with the same data
mod_int <- GLMER(insect_sites, responseVar = "log_scaled_abundance", fitFamily = "gaussian",
                 fixedStruct = "LandUse*std_climate_anomaly*percNH.s",
                 randomStruct = "(1|SS) + (1|SSB)")

mod_add_test <- GLMER(mod_int$data, responseVar = "log_scaled_abundance", fitFamily = "gaussian",
                 fixedStruct = "LandUse*std_climate_anomaly+LandUse*percNH.s",
                 randomStruct = "(1|SS) + (1|SSB)")


mod_no_percHN_test <- GLMER(mod_int$data, responseVar = "log_scaled_abundance", fitFamily = "gaussian",
                       fixedStruct = "LandUse*std_climate_anomaly",
                       randomStruct = "(1|SS) + (1|SSB)")

mod_just_landuse_test <- GLMER(mod_int$data, responseVar = "log_scaled_abundance", fitFamily = "gaussian",
                          fixedStruct = "LandUse",
                          randomStruct = "(1|SS) + (1|SSB)")


anova(mod_just_landuse_test$model,mod_no_percHN_test$model, mod_add_test$model,mod_int$model )


##Model checks

#interactive model
qqnorm(residuals(mod_int$model))
qqline(residuals(mod_int$model))
abline(v=2, col="red")
abline(v=-2, col="red")
plot(fitted(mod_int$model), residuals(mod_int$model))


#additive model
qqnorm(residuals(mod_add$model))
qqline(residuals(mod_add$model))
abline(v=2, col="red")
abline(v=-2, col="red")
plot(fitted(mod_add$model), residuals(mod_add$model))


