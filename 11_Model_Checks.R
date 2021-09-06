##%######################################################%##
#                                                          #
####                    Model Checks                    ####
#                                                          #
##%######################################################%##

# run checks for models used in paper

rm(list = ls())

# load libraries
library(cowplot)
library(ggplot2)
library(roquefort)


# directories
LUclimmod <- "6_RunLUClimateModels/"
nathabmod <- "7_RunLUClimateNHModels/"
predictsDataDir <- "6_RunLUClimateModels/"
outdir <- "11_Model_Checks/"
dir.create(outdir)

# read in model files
load(paste0(LUclimmod, "MeanAnomalyModelAbund.rdata"))
load(paste0(LUclimmod, "MeanAnomalyModelRich.rdata"))
load(paste0(LUclimmod, "MaxAnomalyModelAbund.rdata"))
load(paste0(LUclimmod, "MaxAnomalyModelRich.rdata"))

load(paste0(nathabmod, "MeanAnomalyModelAbun_NH.rdata"))
load(paste0(nathabmod, "RichMeanAnomalyModel_NH.rdata"))
load(paste0(nathabmod, "AbundMaxAnomalyModel_NH.rdata"))
load(paste0(nathabmod, "RichMaxAnomalyModel_NH.rdata"))


# try to write a function to run these checks on each model output and save a pdf of all figs.

mod_list <- c("MeanAnomalyModelAbund", "MeanAnomalyModelRich", "MaxAnomalyModelAbund", "MaxAnomalyModelRich",
              "AbundMeanAnomalyModel1", "RichMeanAnomalyModel1", "AbundMaxAnomalyModel1", "RichMaxAnomalyModel1")


# load dataset
predictsSites <- readRDS(paste0(predictsDataDir,"PREDICTSSiteData.rds"))


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



##%######################################################%##
#                                                          #
####      Additional checks - performance package       ####
#                                                          #
##%######################################################%##

library(performance)

# first the LU:CC models
check_model(MeanAnomalyModelAbund$model, check = c("vif", "normality", "ncv", "homogeneity", "reqq"))
check_model(MeanAnomalyModelRich$model, check = c("vif",  "ncv", "homogeneity", "reqq")) # 'check_normality()' does not support models of class 'glmerMod'.
check_model(MaxAnomalyModelAbund$model, check = c("vif", "normality", "ncv", "homogeneity", "reqq"))
check_model(MaxAnomalyModelRich$model, check = c("vif",  "ncv", "homogeneity", "reqq"))


# problem with homogeneity of variance in the abundance models (try adding in block?)

modelData <- MeanAnomalyModelAbund$data
check_outliers(MeanAnomalyModelAbund$model)
# Warning: 3 outliers detected (cases 229, 1815, 4874).

modelData <- MeanAnomalyModelRich$data
check_outliers(MeanAnomalyModelRich$model)
# Warning: 1 outliers detected (cases 229).
# Warning message:  In hatvalues.merMod(model) : the hat matrix may not make sense for GLMMs

modelData <- MaxAnomalyModelAbund$data
check_outliers(MaxAnomalyModelAbund$model)
# Warning: 3 outliers detected (cases 229, 1815, 4874).

modelData <- MaxAnomalyModelRich$data
check_outliers(MaxAnomalyModelRich$model)
# Warning: 1 outliers detected (cases 229).
# Warning message: In hatvalues.merMod(model) : the hat matrix may not make sense for GLMMs

# check autocorrelation of the residuals
check_autocorrelation(MeanAnomalyModelAbund$model)
# Warning: Autocorrelated residuals detected (p < .001).
check_autocorrelation(MeanAnomalyModelRich$model)
# Warning: Autocorrelated residuals detected (p < .001).
check_autocorrelation(MaxAnomalyModelAbund$model)
# Warning: Autocorrelated residuals detected (p < .001).
check_autocorrelation(MaxAnomalyModelRich$model)
# Warning: Autocorrelated residuals detected (p < .001).

### then the LU:CC:NH models ###

check_model(AbundMeanAnomalyModel1$model, check = c("vif", "normality", "ncv", "homogeneity", "reqq"))
check_model(RichMeanAnomalyModel1$model, check = c("vif",  "ncv", "homogeneity", "reqq")) # 'check_normality()' does not support models of class 'glmerMod'.
check_model(AbundMaxAnomalyModel1$model, check = c("vif", "normality", "ncv", "homogeneity", "reqq"))
check_model(RichMaxAnomalyModel1$model, check = c("vif",  "ncv", "homogeneity", "reqq"))


# problem with homogeneity of variance in the abundance models (try adding in block?)

modelData <- AbundMeanAnomalyModel1$data
check_outliers(AbundMeanAnomalyModel1$model)
# OK: No outliers detected.

modelData <- RichMeanAnomalyModel1$data
check_outliers(RichMeanAnomalyModel1$model)
# OK: No outliers detected.
# Warning message: In hatvalues.merMod(model) : the hat matrix may not make sense for GLMMs

modelData <- AbundMaxAnomalyModel1$data
check_outliers(AbundMaxAnomalyModel1$model)
# OK: No outliers detected.

modelData <- RichMaxAnomalyModel1$data
check_outliers(RichMaxAnomalyModel1$model)
# OK: No outliers detected.
# Warning message: In hatvalues.merMod(model) : the hat matrix may not make sense for GLMMs


# check autocorrelation of the residuals
check_autocorrelation(AbundMeanAnomalyModel1$model)
# Warning: Autocorrelated residuals detected (p < .001).
check_autocorrelation(RichMeanAnomalyModel1$model)
# Warning: Autocorrelated residuals detected (p < .001).
check_autocorrelation(AbundMaxAnomalyModel1$model)
# Warning: Autocorrelated residuals detected (p < .001).
check_autocorrelation(RichMaxAnomalyModel1$model)
# Warning: Autocorrelated residuals detected (p < .001).


##%######################################################%##
#                                                          #
####         Additional checks - DHARMa package         ####
#                                                          #
##%######################################################%##

library(DHARMa)

# working through using vignette https://cran.r-project.org/web/packages/DHARMa/vignettes/DHARMa.html

simulationOutput <- simulateResiduals(fittedModel = MeanAnomalyModelAbund$model, plot = T)

plot(simulationOutput)

plotResiduals(simulationOutput, form = MeanAnomalyModelAbund$data$UI2)
plotResiduals(simulationOutput, form = MeanAnomalyModelAbund$data$StdTmeanAnomalyRS)


simulationOutput2 <- simulateResiduals(fittedModel = MeanAnomalyModelRich$model, plot = F)

plot(simulationOutput2)


plotResiduals(simulationOutput2, form = MeanAnomalyModelRich$data$UI2)
plotResiduals(simulationOutput2, form = MeanAnomalyModelRich$data$StdTmeanAnomalyRS)


