##%######################################################%##
#                                                          #
####        Predicting values from model outputs        ####
#                                                          #
##%######################################################%##

# Here, I predict values for certain fixed effect combinations and 
# express as percentage declines for use in the text of the paper. 

rm(list = ls())

# load libraries
library(StatisticalModels)
library(predictsFunctions)


source('Functions.R')


# directories
predictsDataDir <- "6_RunLUClimateModels/"
moddir <- "6_RunLUClimateModels/"
outdir <- "8_Predictions/"
dir.create(outdir)

# read in the predicts data
predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSiteData.rds"))


#### Hyp 1: land use effect only ####

# these values can be retrieved from the info in script 2_RunSimpleLUIModel.R
# used the median values that were used for plotting of Figure 1.




#### Hyp 2: Land use and climate anomaly interaction ####

# first one, looking at SCA of 1

# load in models
load(file = paste0(moddir, "/MeanAnomalyModelAbund.rdata"))
load(file = paste0(moddir, "/MeanAnomalyModelRich.rdata"))
#load(file = paste0(moddir, "/MaxAnomalyModelAbund.rdata"))
#load(file = paste0(moddir, "/MaxAnomalyModelRich.rdata"))

# create matrix for predictions
# Primary, Low, High
# SCA = 1
# abun and richness = 0

# what is the rescaled value of SCA of 1
BackTransformCentreredPredictor(transformedX = 1.68, originalX = predictsSites$StdTmeanAnomaly) # 1.68 gives about 1 

# what is the rescaled value of SCA of 0
BackTransformCentreredPredictor(transformedX = -0.9, originalX = predictsSites$StdTmeanAnomaly) # -0.9 gives about 0 

# reference is primary with 0 climate change so have 0 for that row

data_tab <- data.frame(UI2 = c("Primary vegetation", "Agriculture_Low", "Agriculture_High", "Agriculture_Low", "Agriculture_High"), 
                       StdTmeanAnomalyRS = c(-0.9,1.68,1.68,-0.9,-0.9),
                       LogAbund = 0,
                       Species_richness = 0)

# factor the LU info
data_tab$UI2 <- factor(data_tab$UI2, levels = levels(predictsSites$UI2))

# predict the results
result.sr <- PredictGLMER(model = MeanAnomalyModelRich$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.sr <- exp(result.sr)

# add in the LU info
result.sr$UI2 <- data_tab$UI2

# express as a percentage of primary
result.sr$perc <- ((result.sr$y/result.sr$y[1]) * 100) - 100

# add in SCA vals
result.sr$SCA <- c(0, 1, 1, 0,0)  

# now for the abundance model  
result.ab <- PredictGLMER(model = MeanAnomalyModelAbund$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.ab <- exp(result.ab)-0.01

# add in the LU info
result.ab$UI2 <- data_tab$UI2

# express as a percentage of primary
result.ab$perc <- ((result.ab$y/result.ab$y[1]) * 100) - 100

# add in SCA vals
result.ab$SCA <- c(0, 1, 1, 0,0) 

# combine results into a table for saving
all_res <- rbind(result.ab, result.sr)

all_res$measure <- c(rep("ab", 5), rep("sr", 5))

# save table
write.csv(all_res, file = paste0(outdir, "/percentage_change_LU_CC.csv"))





### Hyp 3:  Land use, climate and natural habitat interactions ###

# 25% NH cover = -1.015469 in rescaled values
# 50% NH cover = 0.01493712 in rescaled values
# 75% NH cover = 1.045653 in rescaled values
# 100% NH cover = 2.073849 in rescaled values

# what is the reference here
# primary vegetation, SCA = 0, NH = 100

# what is the resscaled value of NH
#BackTransformCentreredPredictor(transformedX = 2.073849, originalX = predictsSites$NH_5000) # -0.99 gives about 0 


load(file = "7_RunLUClimateNHModels/MeanAnomalyModelAbun_NH.rdata")



data_tab <- data.frame(UI2 = c("Primary vegetation", rep("Agriculture_Low", 2), rep("Agriculture_High", 2)), 
                       StdTmeanAnomalyRS = c(-0.9,1.68,1.68, 1.68,1.68),
                       NH_5000.rs = c(2.073849, 1.045653, -1.015469, 1.045653, -1.015469 ), #100, 75, 25
                       LogAbund = 0,
                       Species_richness = 0)

# factor the LU info
data_tab$UI2 <- factor(data_tab$UI2, levels = levels(predictsSites$UI2))


# now for the abundance model  
result.ab <- PredictGLMER(model = AbundMeanAnomalyModel1$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

# backtransform
result.ab <- exp(result.ab)-0.01

# add in the LU info
result.ab$UI2 <- data_tab$UI2

result.ab$NH <- c(100, 75, 25, 75, 25)
result.ab$CC <- c(0, 1, 1, 1, 1)

# express as a percentage of primary
result.ab$perc <- ((result.ab$y/result.ab$y[1]) * 100) - 100

result.ab$metric <- "ab"

# save table
write.csv(result.ab, file = paste0(outdir, "/percentage_change_LU_CC_NH.csv"))



