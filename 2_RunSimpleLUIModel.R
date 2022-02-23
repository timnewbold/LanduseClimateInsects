##%######################################################%##
#                                                          #
####           Land use/use intensity models            ####
#                                                          #
##%######################################################%##

# This script runs simple models of biodiversity responses to
# land use/use intensity only

# directories
inDir <- "1_PreparePREDICTSData/"
outDir <- "2_RunSimpleLUIModel/"

if(!dir.exists(outDir)) dir.create(outDir)

sink(paste0(outDir,"log.txt"))

t.start <- Sys.time()

print(t.start)

# load required libraries
library(StatisticalModels)
library(ggplot2)
library(cowplot)
library(sjPlot)


sessionInfo()

# read in the formatted PREDICTS data
sites <- readRDS(file = paste0(inDir,"PREDICTSSiteData.rds")) # 6095 rows


##%######################################################%##
#                                                          #
####              SPECIES RICHNESS MODELS               ####
#                                                          #
##%######################################################%##


# remove NAs in the specified columns
model_data_sr <- na.omit(sites[,c('Species_richness','LandUse','Use_intensity','UI2','SS','SSB','SSBS')])

# summaries
length(unique(model_data_sr$SS)) # 264
length(unique(model_data_sr$SSBS)) # 6095


# look at the spread of land use/use intensity categories
print(table(model_data_sr$UI2))

# Agriculture_High      Agriculture_Low   Primary vegetation Secondary vegetation 
# 1779                 1317                 1516                 1483

# run set of simple models with different fixed effects structures
sm0 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "1",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

sm1 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "LandUse",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

sm2 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

sm3 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "UI2",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

sm4 <- GLMER(modelData = model_data_sr,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)
# fixed-effect model matrix is rank deficient so dropping 3 columns / coefficients

# take a look at the AICs
print(AIC(sm0$model,sm1$model,sm2$model,sm3$model,sm4$model))
# df      AIC
# sm0$model  4 35851.75
# sm1$model 11 35581.07
# sm2$model 13 35581.86
# sm3$model  7 35614.58
# sm4$model 24 35538.23

##%######################################################%##
#                                                          #
####                 ABUNDANCE MODELS                   ####
#                                                          #
##%######################################################%##


model_data_ab <- na.omit(sites[,c('LogAbund','LandUse','Use_intensity','UI2','SS','SSB','SSBS')])

# summaries
length(unique(model_data_ab$SS)) # 244
length(unique(model_data_ab$SSBS)) # 5759

# look at the spread of land use/use intensity categories
print(table(model_data_ab$UI2))

# Agriculture_High      Agriculture_Low   Primary vegetation Secondary vegetation 
# 1717                 1294                 1410                 1338

# run set of simple models with different fixed effects structures
am0 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "1",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am1 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "LandUse",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am2 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am3 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "UI2",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

am4 <- GLMER(modelData = model_data_ab,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)
#fixed-effect model matrix is rank deficient so dropping 3 columns / coefficients



# take a look at the AICs
print(AIC(am0$model,am1$model,am2$model,am3$model,am4$model))
# df      AIC
# am0$model  4 15504.51
# am1$model 11 15380.29
# am2$model 13 15380.70
# am3$model  7 15378.95
# am4$model 24 15347.04

##%######################################################%##
#                                                          #
####          Predict responses for plotting            ####
#                                                          #
##%######################################################%##


# create dataframe for values to predict response to
nd <- data.frame(UI2=factor(c("Primary vegetation","Secondary vegetation",
                              "Agriculture_Low","Agriculture_High")),
                 Species_richness=0,
                 LogAbund=0)

## species richness predictions ##
s.preds <- PredictGLMERRandIter(model = sm3$model, data = nd)

s.preds <- exp(s.preds)

# convert to percentage difference from primary vegetation
s.preds <- sweep(x = s.preds, MARGIN = 2, STATS = s.preds[1,], FUN = '/')

# get quantiles
s.preds.median <- ((apply(X = s.preds,MARGIN = 1,FUN = median))*100)-100
s.preds.upper <- ((apply(X = s.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
s.preds.lower <- ((apply(X = s.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100


## abundance predictions ##

a.preds <- PredictGLMERRandIter(model = am3$model,data = nd)

a.preds <- exp(a.preds)-1

# convert to percentage difference from primary vegetation
a.preds <- sweep(x = a.preds,MARGIN = 2,STATS = a.preds[1,],FUN = '/')

# get quantiles
a.preds.median <- ((apply(X = a.preds,MARGIN = 1,FUN = median))*100)-100
a.preds.upper <- ((apply(X = a.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
a.preds.lower <- ((apply(X = a.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100



# combine data into one table for plotting
abun_res <- as.data.frame(cbind(a.preds.median, a.preds.lower, a.preds.upper))
rich_res <- as.data.frame(cbind(s.preds.median, s.preds.lower, s.preds.upper))
colnames(abun_res) <- c("median", "lower", "upper")
colnames(rich_res) <- c("median", "lower", "upper")
abun_res$metric <- "abun"
rich_res$metric <- "rich"
abun_res$LU <- factor(c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"), levels = c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"))
rich_res$LU <- factor(c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"), levels = c("Primary vegetation","Secondary vegetation", "Agriculture_Low","Agriculture_High"))

abun_res[abun_res$LU == "Primary vegetation", c("lower", "upper")] <- NA
rich_res[abun_res$LU == "Primary vegetation", c("lower", "upper")] <- NA


##%######################################################%##
#                                                          #
####  plot species richness and abundance predictions   ####
#                                                          #
##%######################################################%##


# Figure 1, includes map of sites:

# load world map
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
        axis.ticks = element_blank(),
        #plot.title = element_text(hjust=0.05, size = 8, face = 'bold')
        title = element_text(size = 8, face = "bold"))+
  ggtitle("a")

#pdf(file = paste0(outDir,"LUI_Plot.pdf"),width = 8.5/2.54,height = 12/2.54, onefile = T)


# point plots

p2 <- ggplot(data = abun_res) +
      geom_point(aes(x = LU, y = median, col = LU), size = 1.2) + 
      geom_errorbar(aes(x = LU, ymin = lower, ymax = upper, col = LU), size = 0.2, width = 0.2)+
      geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
      xlab("") +
      ylab("Change in total abundance (%)") +
      scale_color_manual(values = c("#009E73","#0072B2","#E69F00","#D55E00")) +
      theme(legend.position = "none", 
            aspect.ratio = 1, 
            title = element_text(size = 8, face = "bold"),
            axis.text.y = element_text(size = 7),
            axis.text.x = element_text(size = 7, angle = 45, vjust = 0.5),
            axis.title = element_text(size = 7),
            panel.grid.minor = element_blank(),
            panel.grid.major = element_blank(),
            panel.border = element_blank(), 
            panel.background = element_blank(), 
            axis.ticks = element_line(size = 0.2), 
            axis.line = element_line(size = 0.2)) +
      ggtitle("b")


p3 <- ggplot(data = rich_res) +
  geom_point(aes(x = LU, y = median, col = LU), size = 1.2) + 
  geom_errorbar(aes(x = LU, ymin = lower, ymax = upper, col = LU), size = 0.2, width = 0.2)+
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.2) +
  xlab("") +
  ylab("Change in species richness (%)") +
  scale_color_manual(values = c("#009E73","#0072B2","#E69F00","#D55E00")) +
  ylim(c(-50, 0)) +
  theme(legend.position = "none", 
        aspect.ratio = 1, 
        title = element_text(size = 8, face = "bold"),
        axis.text.y = element_text(size = 7),
        axis.text.x = element_text(size = 7, angle = 45, vjust = 0.5),
        axis.title = element_text(size = 7),
        panel.grid.minor = element_blank(),
        panel.grid.major = element_blank(),
        panel.border = element_blank(), 
        panel.background = element_blank(), 
        axis.ticks = element_line(size = 0.2), 
        axis.line = element_line(size = 0.2)) +
  ggtitle("c")


p5 <- cowplot::plot_grid(p1, cowplot::plot_grid(p2, p3), ncol = 1, rel_heights = c(1, 1))

ggsave(filename = paste0(outDir, "Figure1_map_simplemods.pdf"), plot = last_plot(), width = 120, height = 150, units = "mm", dpi = 300)


### plotting using base plot if preferred ###

# #par(mfrow=c(2,1))
# par(cex=1)
# par(cex.lab=1)
# par(cex.axis=1)
# par(cex.main=1)
# par(ps=10)
# par(las=1)
# par(mgp=c(2.5,1,0))
# par(mar=c(4,3.5,0.1,1.5))
# par(tck=-0.01)
# 
# # set colours
# errbar.cols <- c("#009E73","#0072B2","#E69F00","#D55E00")
# 
# # species richness plot
# 
# 
# errbar(x = 1:4,y = s.preds.median,yplus = s.preds.upper,yminus = s.preds.lower,
#        col=errbar.cols,errbar.col = errbar.cols,ylim=c(-50,0),xaxt="n",
#        ylab="Change in species richness (%)",xlab="",bty="l")
# 
# axis(side = 1,at = 1:4, labels = c("Primary","Secondary","Agriculture\nLow","Agriculture\nHigh"), cex.axis = 1, las = 2)
# 
# 
# abline(h=0,col="#00000077",lty=2)
# 
# #title(main = "b.", adj = 0, cex.main = 1, line = 1)
# 
# p3 <- recordPlot()
# 
# # abundance plot
# 
# errbar(x = 1:4,y = a.preds.median,yplus = a.preds.upper,yminus = a.preds.lower,
#        col=errbar.cols,errbar.col = errbar.cols,ylim=c(-50,0),xaxt="n",
#        ylab="Change in total abundance (%)",xlab="",bty="l")
# 
# axis(side = 1,at = 1:4,labels = c("Primary","Secondary","Agriculture\nLow","Agriculture\nHigh"), cex.axis = 1, las = 2)
# 
# abline(h=0,col="#00000077",lty=2)
# 
# #title(main = "c.", adj = 0, cex.main = 1, line = 1)
# 
# p2 <- recordPlot()
# 
# 
# #invisible(dev.off())
# 
# # Organise the plots
# # 
# # p4 <- plot_grid(p1, p2, p3, ncol = 1, scale = 0.85, labels = c("a", "b", "c"), label_x = 0.1)
# # 
# # # save figure
# # save_plot(paste0(outDir, "Figure_1.pdf"), p4, base_height = 8, base_width = 6)
# # 
# 
# # alternative fig 1 format
# 
# p5 <- plot_grid(p1, plot_grid(p2, p3, scale = 0.78, labels = c("b", "c"), label_x = 0.2), ncol = 1, rel_heights = c(0.7, 1), labels = c("a"), label_x = 0.1)
# 
# save_plot(paste0(outDir, "Figure_1_alt.pdf"), p5, base_height = 9, base_width = 9)
# 


##%######################################################%##
#                                                          #
####                    Model stats                     ####
#                                                          #
##%######################################################%##


# models used = sm3 and am3

summary(sm3$model)

# rerun the model using GLMERSelect function to easily get stats

sm3_2 <- GLMERSelect(modelData = model_data_sr,
                     responseVar = "Species_richness",
                      fitFamily = "poisson",
                      fixedFactors = "UI2",
                      #fixedTerms = list(StdTmeanAnomalyRS=1),
                      randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)"#,
                      #fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                      #saveVars = c("Total_abundance", "SSBS", "NH_3000")
                      )


am3_2 <- GLMERSelect(modelData = model_data_ab,
                     responseVar = "LogAbund",
                     fitFamily = "gaussian",
                     fixedFactors = "UI2",
                     #fixedTerms = list(StdTmeanAnomalyRS=1),
                     randomStruct = "(1|SS)+(1|SSB)"#,
                     #fixedInteractions = c("UI2:poly(StdTmeanAnomalyRS,1)"),
                     #saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_3000")
                     )

summary(sm3_2$model)
summary(sm3$model)
summary(am3_2$model)
summary(am3$model)

# save the stats info
sm3stats <- as.data.frame(sm3_2$stats)
am3stats <- as.data.frame(am3_2$stats)

sm3stats$significant <- NA
am3stats$significant <- NA


# function to check significance
checksig <- function(x){
  if(x <= 0.05){ 
    res <- "Yes" 
  } else { 
    res <- "No" }
  return(res)}

# add values to table
sm3stats$significant <- sapply(X = sm3stats$P, FUN = checksig)
am3stats$significant <- sapply(X = am3stats$P, FUN = checksig)


# save the stats tables
write.csv(sm3stats, file = paste0(outDir, "/SR_Stats.csv"), row.names = FALSE)
write.csv(am3stats, file = paste0(outDir, "/Abun_Stats.csv"), row.names = FALSE)



### save model output tables ###

tab_model(am3_2$model, transform = NULL, file = paste0(outDir, "/AbunLU_output_table.html"))
summary(am3_2$model)
R2GLMER(am3_2$model)

tab_model(sm3_2$model, transform = NULL, file = paste0(outDir, "/SRLU_output_table.html"))
summary(sm3_2$model)
R2GLMER(sm3_2$model)

t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()
