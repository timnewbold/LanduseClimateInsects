##%######################################################%##
#                                                          #
####                    Model Checks                    ####
#                                                          #
##%######################################################%##

# In this scripts, model checks are carried out and plots presented in 
# supplementary information generated.

rm(list = ls())

# directories
LUclimmod <- "6_RunLUClimateModels/"
nathabmod <- "7_RunLUClimateNHModels/"
predictsDataDir <- "6_RunLUClimateModels/"
outdir <- "11_Model_Checks/"
if(!dir.exists(outdir)) dir.create(outdir)


sink(paste0(outdir,"log.txt"))

t.start <- Sys.time()

print(t.start)

# load libraries
library(cowplot)
library(ggplot2)
library(StatisticalModels)

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
predictsSites <- predictsSites[!is.na(predictsSites$StdTmeanAnomalyRS), ]

#x <- mod_list[1]

# loop through models and generate plots for checking assumptions etc
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

sa_test<-SpatialAutocorrelationTest(model=get(x), all.data=predictsSites)

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
  
  
 predData <- predictsSites[!is.na(predictsSites$LogAbund), ]
 
 
 # 4. plot of observed vs fitted values
  pdf(NULL)
 dev.control(displaylist="enable")
 plot(predData$LogAbund,fitted(get(x)$model), 
      xlab = "Observed values", ylab = "Fitted values") 
 abline(a = 0, b = 1, col = "red", lwd = 2)
 p4 <- recordPlot()
 invisible(dev.off())
  

 cowplot::plot_grid(p1,p2,p3, p4,
          labels = c("A.", "B.", "C.", "D."))

  ggsave(file = paste0(outdir, x, "_model_checks.pdf"), height = 10, width = 10)
  
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

ggsave(file = paste0(outdir, x, "_model_checks.pdf"), height = 10, width = 10) 

rm(p2, p3, p4)
rm(perc_auto)

}

}




##%######################################################%##
#                                                          #
####    Additional plots: residuals with anomaly/LU     ####
#                                                          #
##%######################################################%##

# additional plots for the abundance models to investigate residuals

mod_list <- c("MeanAnomalyModelAbund", "MaxAnomalyModelAbund", 
              "AbundMeanAnomalyModel1",  "AbundMaxAnomalyModel1")


mod_res <- resid(MeanAnomalyModelAbund$model)

plot_data1 <- cbind(mod_res, MeanAnomalyModelAbund$data[, c("UI2", "StdTmeanAnomalyRS")])

p1 <- ggplot(data = plot_data1) +
  geom_point(aes(x = StdTmeanAnomalyRS, y = mod_res, col = UI2), size = 0.5)+
  scale_colour_manual(values = c("#009E73","#D55E00", "#E69F00", "#0072B2")) +
  theme_bw() + 
  ggtitle("Model of total abundance as a function of the STA \nand Land use") + 
  ylab("Residuals")+
  theme(legend.position = "bottom", legend.title = element_blank(), text = element_text(size = 8))



mod_res2 <- resid(MaxAnomalyModelAbund$model)

plot_data2 <- cbind(mod_res2, MaxAnomalyModelAbund$data[, c("UI2", "StdTmaxAnomalyRS")])

p2 <- ggplot(data = plot_data2) +
  geom_point(aes(x = StdTmaxAnomalyRS, y = mod_res, col = UI2), size = 0.5)+
  scale_colour_manual(values = c("#009E73","#D55E00", "#E69F00", "#0072B2")) +
  theme_bw() + 
  ggtitle("Model of total abundance as a function of the STMA \nand Land use") + 
  ylab("Residuals") +
  theme(legend.position = "none", text = element_text(size = 8))


mod_res3 <- resid(AbundMeanAnomalyModel1$model)

plot_data3 <- cbind(mod_res3, AbundMeanAnomalyModel1$data[, c("UI2", "StdTmeanAnomalyRS")])

p3 <- ggplot(data = plot_data3) +
  geom_point(aes(x = StdTmeanAnomalyRS, y = mod_res, col = UI2), size = 0.5)+
  scale_colour_manual(values = c("#009E73","#D55E00", "#E69F00", "#0072B2")) +
  theme_bw() + 
  ggtitle("Model of total abundance as a function of the STA, \nLand use and natural habitat availability") + 
  ylab("Residuals") +
  theme(legend.position = "none", text = element_text(size = 8))


mod_res4 <- resid(AbundMaxAnomalyModel1$model)

plot_data4 <- cbind(mod_res4, AbundMaxAnomalyModel1$data[, c("UI2", "StdTmaxAnomalyRS")])

p4 <- ggplot(data = plot_data4) +
  geom_point(aes(x = StdTmaxAnomalyRS, y = mod_res, col = UI2), size = 0.5)+
  scale_colour_manual(values = c("#009E73","#D55E00", "#E69F00", "#0072B2")) +
  theme_bw() + 
  ggtitle("Model of total abundance as a function of the STMA, \nLand use and natural habitat availability") + 
  ylab("Residuals") +
  theme(legend.position = "none", text = element_text(size = 8))

legend <- get_legend(
  # create some space to the left of the legend
  p1 + theme(legend.box.margin = margin(0, 0, 0, 12), text = element_text(size = 12))
)



plots <- cowplot::plot_grid(p1 + theme(legend.position = "none"), p2, p3, p4, 
                   nrow = 2)


cowplot::plot_grid(plots, legend, 
                   nrow = 2, 
                   rel_heights = c(6, 1))

ggsave(filename = paste0(outdir, "Residual_plots_by_anom_All.pdf"), height = 7, width = 8)


t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()