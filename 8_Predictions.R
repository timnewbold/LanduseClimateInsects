##%######################################################%##
#                                                          #
####        Predicting values from model outputs        ####
#                                                          #
##%######################################################%##

# Here, I predict values for certain fixed effect combinations and 
# express as percentage declines for use in the text of the paper. 

rm(list = ls())

# load libraries




### Hyp 1: land use effect only ###

# these values can be retrieved from the info in script 2_RunSimpleLUIModel.R
# used the median values that were used for plotting of Figure 1.


### Hyp 2: Land use and climate anomaly interaction ###


# create matrix for predictions
# Primary, Low, High
# SCA = 1
# abun and richness = 0

# what is the resscaled value of SCA of 1
BackTransformCentreredPredictor(transformedX = 1.4, originalX = predictsSites$StdTmeanAnomaly) # 1.4 gives about 1 

data_tab <- data.frame(UI2 = c("Primary vegetation", "Agriculture_Low", "Agriculture_High"), 
                       StdTmeanAnomalyRS = c(1.4,1.4,1.4),
                       LogAbund = 0,
                       Species_richness = 0)


data_tab$UI2 <- factor(data_tab$UI2, levels = levels(predictsSites$UI2))

result.sr <- PredictGLMER(model = MeanAnomalyModelRich$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

result.sr <- exp(result.sr)

result.sr$UI2 <- data_tab$UI2

# express as a percentage of primary
result$perc <- 

result.ab <- PredictGLMER(model = MeanAnomalyModelAbund$model, data = data_tab, se.fit = TRUE, seMultiplier = 1.96)

