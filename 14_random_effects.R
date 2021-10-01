##%######################################################%##
#                                                          #
####          testing random effect structures          ####
#                                                          #
##%######################################################%##



# load libraries
library(devtools)
#install_github("timnewbold/StatisticalModels")
library(StatisticalModels)
#install_github(repo = "timnewbold/predicts-demo",subdir = "predictsFunctions", force = T)
library(predictsFunctions)
source("Functions.R")
library(sjPlot)

# directories 
predictsDataDir <- "6_RunLUClimateModels/"
outDir <- "14_Additional_Tests/"

# load in the data
predictsSites <- readRDS(file = paste0(predictsDataDir,"PREDICTSSiteData.rds"))


#### abundance models ####

# 1. original formulation
mod1 <- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                                     fitFamily = "gaussian",fixedFactors = "UI2",
                                     fixedTerms = list(StdTmeanAnomalyRS=1),
                                     randomStruct = "(1|SS)+(1|SSB)",
                                     fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                     saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(mod1, file = paste0(outDir, "/MeanAnomalyModelAbund.rdata"))

mod1 <- lmer(formula = LogAbund ~ UI2*StdTmeanAnomalyRS + (1|SS)+(1|SSB), 
              data = predictsSites)

plot(mod1)
qqnorm(resid(mod1))
qqline(resid(mod1))

mod2 <- lmer(formula = LogAbund ~ UI2*StdTmeanAnomalyRS + (1|SS)+(1|SS:SSB), 
             data = predictsSites)

plot(mod2)
qqnorm(resid(mod2))
qqline(resid(mod2))


# 2. Adding SSS
mod2<- GLMERSelect(modelData = predictsSites,responseVar = "LogAbund",
                    fitFamily = "gaussian",fixedFactors = "UI2",
                    fixedTerms = list(StdTmeanAnomalyRS=1),
                    randomStruct = "(1|SS)+(1|SSB)",
                    fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                    saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000"))

# save the model output
save(mod2, file = paste0(outDir, "/MeanAnomalyModelAbund.rdata"))



#### richness models ####

# 1. Richness, mean anomaly
srmod1 <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
                                    fitFamily = "poisson",fixedFactors = "UI2",
                                    fixedTerms = list(StdTmeanAnomalyRS=1),
                                    randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",
                                    fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                                    saveVars = c("Total_abundance", "SSBS", "NH_3000"))

# save model output
save(srmod1, file = paste0(outDir, "/MeanAnomalyModelRich.rdata"))

# 2. Richness, mean anomaly
srmod2 <- GLMERSelect(modelData = predictsSites,responseVar = "Species_richness",
                      fitFamily = "poisson",fixedFactors = "UI2",
                      fixedTerms = list(StdTmeanAnomalyRS=1),
                      randomStruct = "(1|SS)+(1|SSB)+(1|SSS)+(1|SSBS)",
                      fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                      saveVars = c("Total_abundance", "SSBS", "NH_3000"))

# save model output
save(srmod1, file = paste0(outDir, "/MeanAnomalyModelRich.rdata"))


#### model check plots ####

mod_list <- c("mod1")


#x <- mod_list[1]

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
    
    cowplot::plot_grid(p1,p2,p3,
                       labels = c("A.", "B.", "C."))
    ggsave(file = paste0(outdir, x, "_model_checks.pdf"), height = 10, width = 10)
    
    rm(p1, p2, p3)
    rm(perc_auto)
    
    
  }else{
    cowplot::plot_grid(p2,p3,
                       labels = c("A.", "B."))#
    
    ggsave(file = paste0(outdir, x, "_model_checks.pdf"), height = 5, width = 10) 
    
    rm(p2, p3)
    rm(perc_auto)
    
  }
  
}







