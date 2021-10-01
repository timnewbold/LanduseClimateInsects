##%######################################################%##
#                                                          #
####     Testing different abundance transformation     ####
#                                                          #
##%######################################################%##

rm(list = ls())


# load libraries
library(devtools)
#install_github("timnewbold/StatisticalModels")
library(StatisticalModels)
#install_github(repo = "timnewbold/predicts-demo",subdir = "predictsFunctions", force = T)
library(predictsFunctions)
source("Functions.R")
library(sjPlot)

# directories 
dataDir <- "6_RunLUClimateModels/"
outdir <- "14_Additional_Tests/"

# load the data
predictsSites <- readRDS(file = paste0(dataDir,"PREDICTSSiteData.rds"))


# try arcsine squareroot transformation of abundance
predictsSites$ArcAbun <- asin(sqrt(predictsSites$Total_abundance_RS))

summary(predictsSites$ArcAbun )

hist(predictsSites$LogAbund)
hist(predictsSites$ArcAbun)



model_data <- predictsSites[!is.na(predictsSites$ArcAbun), ] # 5759
model_data <- model_data[!is.na(model_data$StdTmeanAnomalyRS), ] # 5735

#### run model ####
MeanAnomalyModelAbund <- GLMERSelect(modelData = model_data,responseVar = "ArcAbun",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("SSBS"))

# save the model output
save(MeanAnomalyModelAbund, file = paste0(outdir, "/MeanAnomalyModelAbund_arcsinesqrt.rdata"))

summary(MeanAnomalyModelAbund$model)



#### model check plots ####

mod_list <- c("MeanAnomalyModelAbund")

for(x in mod_list){
  
  # only do this one if not SR model
  if(grepl("Abun", x) ==1) {
    
    
    ## 1. Checking the fitted vs residuals relationship
    p1 <- plot(get(x)$model)
  }
  
  ## 2. Normality of Residuals
  pdf(NULL)
  dev.control(displaylist="enable")
  qqnorm(resid(get(x)$model), main = "")
  qqline(resid(get(x)$model))
  p2 <- recordPlot()
  invisible(dev.off())
  
  
  
  ## 3. Check for spatial autocorrelation
  
  sa_test<-roquefort::SpatialAutocorrelationTest(model=get(x), all.data=predictsSites)
  
  #summary(sa_test)
  
  # percentage of studies that show spatial autocorrelation?
  perc_auto <- (length(which(sa_test$P<0.05))/length(sa_test$P))*100
  
  
  sa_test_vals <- as.data.frame(sa_test$P)
  sa_test_vals$`sa_test$P` <- round(sa_test_vals$`sa_test$P`, digits = 4)
  
  label1 <- paste0("P < 0.05 \nin ", round(perc_auto, 1), "% \nof studies")
  
  p3 <- ggplot(data = sa_test_vals ) +
    geom_histogram(aes(x = sa_test_vals$`sa_test$P`)) +
    geom_vline(xintercept = 0.05, col = "red") +
    geom_text(aes(x = 0.9, y = 90, label = label1), size = 4, check_overlap = T) +
    theme_bw() +
    ylim(c(0, 100)) +
    xlab("P-value") +
    ylab("Frequency") +
    theme(panel.grid = element_blank(), 
          aspect.ratio = 1)
  
  
  
  if(grepl("Abun", x) == 1) {
    
    
    predData <- predictsSites[!is.na(predictsSites$ArcAbun), ]
    predData <- predData[!is.na(predData$StdTmeanAnomalyRS), ]
    
    # 4. plot of observed vs fitted values
    
    pdf(NULL)
    dev.control(displaylist="enable")
    plot(predData$ArcAbun,fitted(get(x)$model), 
         xlab = "Observed values", ylab = "Fitted values") 
    abline(a = 0, b = 1, col = "red", lwd = 2)
    p4 <- recordPlot()
    invisible(dev.off())
    
    
    cowplot::plot_grid(p1,p2,p3, p4,
                       labels = c("A.", "B.", "C.", "D."))
    
    ggsave(file = paste0(outdir, x, "_model_checks_arcsinsqrt.pdf"), height = 10, width = 10)
    
    rm(p1, p2, p3, p4)
    rm(perc_auto)
    
    
  }else{
    
    
    # 4. plot of observed vs fitted values
    
    pdf(NULL)
    dev.control(displaylist="enable")
    plot(predictsSites$Species_richness,fitted(get(x)$model), 
         xlab = "Observed values", ylab = "Fitted values") 
    abline(a = 0, b = 1, col = "red", lwd = 2)
    p4 <- recordPlot()
    invisible(dev.off())
    
    cowplot::plot_grid(p2,p3,p4,
                       labels = c("A.", "B.", "C."))#
    
    ggsave(file = paste0(outdir, x, "_model_checks_arcsinsqrt.pdf"), height = 10, width = 10) 
    
    rm(p2, p3, p4)
    rm(perc_auto)
    
  }
  
}
