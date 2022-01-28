##%######################################################%##
#                                                          #
####         Organise PREDICTS data for insects         ####
#                                                          #
##%######################################################%##

# This script takes the complete PREDICTS database, selects those entries for
# insects, and organises the data for analysis.

# directories
dataDir <- "0_data/"
outDir <- "1_PreparePREDICTSData/"

if(!dir.exists(outDir)) dir.create(outDir)

sink(paste0(outDir,"log.txt"))

t.start <- Sys.time()

print(t.start)

# load required libraries
library(predictsFunctions)
library(ggplot2)

sessionInfo()

# Set the path to your local copy of the database
predicts.path <- paste0(dataDir,"database.rds")

# Read in the PREDICTS data
predicts <- ReadPREDICTS(predicts.path)

# Select only data for insects
predicts <- predicts[(predicts$Class=="Insecta"),]

# Correct effort-sensitive abundance measures (assumes linear relationship between effort and recorded abundance)
predicts <- CorrectSamplingEffort(diversity = predicts)
# Correcting 870378 values for sensitivity to sampling effort (most of these are 0s, 19184 non-zero)

table(predicts$Diversity_metric)


# insects should not have diversity metric "percent cover", this is a mistake in the database
# remove those sites that are the problem
predicts <- predicts[!predicts$Diversity_metric == "percent cover", ]


# Merge sites that have the same coordinates (e.g. multiple traps on a single transect)
predicts <- MergeSites(diversity = predicts)


# remove rows where land use or use intensity info is missing
predicts.complete <- droplevels(predicts[(
  predicts$Predominant_land_use!="Cannot decide"),])
predicts.complete <- droplevels(predicts.complete[(
  predicts.complete$Use_intensity!="Cannot decide"),])

nrow(predicts.complete)
# 779,912 records

# get counts of n species in major groups
species <- unique(predicts.complete[,c('Order','Taxon_name_entered')])

order.counts <- tapply(X = species$Taxon_name_entered,
                       INDEX = species$Order,
                       FUN = function(sp) length(unique(sp)))



# Calculate site metrics of diversity
sites <- SiteMetrics(diversity = predicts,
                     extra.cols = c("Predominant_land_use",
                                    "SSB","SSBS", "Biome"))

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

# combine LU and UI 
sites$UI <- paste0(sites$LandUse,'_',sites$Use_intensity)
sites$UI[grep("NA",sites$UI)] <- NA

# recode according to land use and use intensity combinations
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

# 
sites$Use_intensity[((sites$LandUse=="Mature secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"
sites$Use_intensity[((sites$LandUse=="Intermediate secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"
sites$Use_intensity[((sites$LandUse=="Young secondary vegetation") & 
                       (sites$Use_intensity=="Intense use"))] <- "Light use"

# remove the urban sites and NA in UI2
sites <- sites[!sites$UI2 == "Urban", ]
sites <- sites[!is.na(sites$UI2), ]


sites <- droplevels(sites)

# transform abundance values 
sites$LogAbund <- log(sites$Total_abundance+1)


# Remove sites without coordinates
sites <- sites[!is.na(sites$Latitude), ]


# save the prepared dataset
saveRDS(object = sites,file = paste0(outDir,"PREDICTSSiteData.rds"))



##%######################################################%##
#                                                          #
#### redo dataset summaries after removing other sites  ####
#                                                          #
##%######################################################%##


predicts2 <- predicts[predicts$SSBS %in% sites$SSBS, ]


table(predicts2$Diversity_metric)


nrow(predicts2)
# 756,879 records

# get counts of n species in major groups
species <- unique(predicts2[,c('Order','Taxon_name_entered')])

order.counts <- tapply(X = species$Taxon_name_entered,
                       INDEX = species$Order,
                       FUN = function(sp) length(unique(sp)))

sum(order.counts)
# 17889

#### basic map of PREDICTS sites ####

# plot the raster in ggplot
map.world <- map_data('world')

# map of sites
p1 <-ggplot() +
         geom_map(data=map.world, map=map.world,
           aes(x=long, y=lat, group=group, map_id=region),
           fill= "grey", colour="grey", size=0.2) +
         geom_point(data = sites, aes(x = Longitude, y = Latitude), col = c("#1E90FF"), fill = c("#104E8B"), shape = 21) +
         theme(axis.title = element_blank(), 
               plot.background = element_blank(), 
               panel.background = element_blank(),
               axis.text = element_blank(),
               axis.ticks = element_blank()) +
         ggtitle("a.")
   

# save plot
ggsave(filename = paste0(outDir, "/PREDICTS_points_map.pdf"), height = 4, width = 8)



### Basic summaries ###

# nstudies/ nsites - all
length(unique(sites$SS)) # 264
length(unique(sites$SSBS)) # 6095

# nstudies/nsites - abun
length(unique(sites[!is.na(sites$LogAbund) , 'SS'])) # 244
length(unique(sites[!is.na(sites$LogAbund) , 'SSBS'])) # 5759

# reviewer request, UI2 by Biome

table(sites$Biome, sites$UI2)

#                                                         Agriculture_High Agriculture_Low Primary vegetation Secondary vegetation
# Boreal Forests/Taiga                                                    2               6                172                   13
# Temperate Conifer Forests                                               4             103                 10                   88
# Temperate Broadleaf & Mixed Forests                                  1308             671                315                  787
# Montane Grasslands & Shrublands                                         2             200                247                   33
# Temperate Grasslands, Savannas & Shrublands                            15              11                 15                   27
# Mediterranean Forests, Woodlands & Scrub                               21              32                 96                   58
# Deserts & Xeric Shrublands                                              0              30                 16                   16
# Tropical & Subtropical Grasslands, Savannas & Shrublands               56              47                175                   78
# Tropical & Subtropical Coniferous Forests                               2              26                 32                   43
# Flooded Grasslands & Savannas                                           6               6                  0                    0
# Tropical & Subtropical Dry Broadleaf Forests                           62              66                 13                   50
# Tropical & Subtropical Moist Broadleaf Forests                        293             119                420                  283
# Mangroves                                                               8               0                  5                    7

table(sites[!is.na(sites$LogAbund), 'Biome'], sites[!is.na(sites$LogAbund), 'UI2'])

#                                                          Agriculture_High Agriculture_Low Primary vegetation Secondary vegetation
# Boreal Forests/Taiga                                                    2               6                172                   13
# Temperate Conifer Forests                                               4             103                 10                   88
# Temperate Broadleaf & Mixed Forests                                  1280             669                289                  729
# Montane Grasslands & Shrublands                                         2             200                247                   33
# Temperate Grasslands, Savannas & Shrublands                            15              11                 15                   27
# Mediterranean Forests, Woodlands & Scrub                               21              26                 92                   58
# Deserts & Xeric Shrublands                                              0              30                 16                   16
# Tropical & Subtropical Grasslands, Savannas & Shrublands               56              47                165                   70
# Tropical & Subtropical Coniferous Forests                               2              26                 32                   43
# Flooded Grasslands & Savannas                                           0               0                  0                    0
# Tropical & Subtropical Dry Broadleaf Forests                           62              66                 13                   50
# Tropical & Subtropical Moist Broadleaf Forests                        265             110                354                  204
# Mangroves                                                               8               0                  5                    7

t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()
