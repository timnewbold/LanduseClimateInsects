
library(predictsFunctions)

dataDir <- "0_data/"
outDir <- "1_PreparePREDICTSData/"

# Set the path to your local copy of the database
predicts.path <- paste0(dataDir,"database.rds")

# Read in the PREDICTS data
predicts <- ReadPREDICTS(predicts.path)

predicts <- predicts[(predicts$Class=="Insecta"),]

# Correct effort-sensitive abundance measures (assumes linear relationship between effort and recorded abundance)
predicts <- CorrectSamplingEffort(diversity = predicts)

# Merge sites that have the same coordinates (e.g. multiple traps on a single transect)
predicts <- MergeSites(diversity = predicts)

# Calculate site metrics of diversity
sites <- SiteMetrics(diversity = predicts,
                     extra.cols = c("Predominant_land_use",
                                    "SSB","SSBS"))

# First, we will rearrange the land-use classification a bit
sites$LandUse <- paste(sites$Predominant_land_use)
# Drop classification where land use could not be identified
sites$LandUse[(sites$LandUse=="Cannot decide")] <- NA
# Drop classification where the stage of recovery of secondary vegetation is unknown
# sites$LandUse[(sites$LandUse=="Secondary vegetation (indeterminate age)")] <- NA
# Now make the variable a factor, and set the reference level to primary vegetation
sites$LandUse <- factor(sites$LandUse)
sites$LandUse <- relevel(sites$LandUse,ref="Primary vegetation")

sites$Use_intensity[sites$Use_intensity=="Cannot decide"] <- NA

sites$UI <- paste0(sites$LandUse,'_',sites$Use_intensity)
sites$UI[grep("NA",sites$UI)] <- NA

sites$UI2 <- dplyr::recode(sites$UI,
                           'Primary vegetation_Minimal use' = 'Primary vegetation',
                           'Cropland_Light use' = 'Agriculture_High',
                           'Secondary vegetation (indeterminate age)_Minimal use' = 'Secondary vegetation',
                           'Urban_Light use' = 'Urban',
                           'Secondary vegetation (indeterminate age)_Light use' = 'Secondary vegetation',
                           'Cropland_Intense use' = 'Agriculture_High',
                           'Cropland_Minimal use' = 'Agriculture_Low',
                           'Pasture_Light use' = 'Agriculture_Low',
                           'Pasture_Minimal use' = 'Agriculture_Low',
                           'Intermediate secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Mature secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Secondary vegetation (indeterminate age)_Intense use' = 'Secondary vegetation',
                           'Pasture_Intense use' = 'Agriculture_High',
                           'Urban_Minimal use' = 'Urban',
                           'Primary vegetation_Light use' = 'Primary vegetation',
                           'Young secondary vegetation_Light use' = 'Secondary vegetation',
                           'Urban_Intense use' = 'Urban',
                           'Primary vegetation_Intense use' = 'Primary vegetation',
                           'Young secondary vegetation_Minimal use' = 'Secondary vegetation',
                           'Mature secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Plantation forest_Minimal use' = 'Agriculture_Low',
                           'Plantation forest_Intense use' = 'Agriculture_High',
                           'Young secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Plantation forest_Light use' = 'Agriculture_High',
                           'Mature secondary vegetation_Light use' = 'Secondary vegetation',
                           'Intermediate secondary vegetation_Intense use' = 'Secondary vegetation',
                           'Intermediate secondary vegetation_Light use' = 'Secondary vegetation')

sites$Use_intensity[((sites$LandUse=="Mature secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"
sites$Use_intensity[((sites$LandUse=="Intermediate secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"
sites$Use_intensity[((sites$LandUse=="Young secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"

sites <- droplevels(sites)

sites$LogAbund <- log(sites$Total_abundance+1)

saveRDS(object = sites,file = paste0(outDir,"PREDICTSSiteData.rds"))
