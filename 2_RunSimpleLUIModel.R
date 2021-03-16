##%######################################################%##
#                                                          #
####           Land use/use intensity models            ####
#                                                          #
##%######################################################%##

# This script runs simple models of biodiversity responses to
# land use/use intensity only

# load required libraries
library(StatisticalModels)
library(ggplot2)
library(cowplot)


# directories
inDir <- "1_PreparePREDICTSData/"
outDir <- "2_RunSimpleLUIModel/"

# read in the formatted PREDICTS data
sites <- readRDS(file = paste0(inDir,"PREDICTSSiteData.rds")) # 6095 rows


## SPECIES RICHNESS MODELS ##

# remove NAs in the specified columns
model_data_sr <- na.omit(sites[,c('Species_richness','LandUse','Use_intensity','UI2','SS','SSB','SSBS')])

# summaries
length(unique(model_data_sr$SS)) # 264
length(unique(model_data_sr$SSBS)) # 6095


# look at the spread of land use/use intensity categories
print(table(model_data_sr$UI2))

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


## ABUNDANCE MODELS ##

model_data_ab <- na.omit(sites[,c('LogAbund','LandUse','Use_intensity','UI2','SS','SSB','SSBS')])

# summaries
length(unique(model_data_ab$SS)) # 244
length(unique(model_data_ab$SSBS)) # 5759

# look at the spread of land use/use intensity categories
print(table(model_data_ab$UI2))


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


### Predict responses for plotting ###


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



### plot species richness and abundance predictions ###

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
        plot.title = element_text(hjust=0.05, size = 11, face = 'bold'))

#pdf(file = paste0(outDir,"LUI_Plot.pdf"),width = 8.5/2.54,height = 12/2.54, onefile = T)


#par(mfrow=c(2,1))
par(cex=1)
par(cex.lab=1)
par(cex.axis=1)
par(cex.main=1)
par(ps=10)
par(las=1)
par(mgp=c(1.6,0.2,0))
par(mar=c(2.5,2.5,1,0.5))
par(tck=-0.01)

# set colours
errbar.cols <- c("#009E73","#0072B2","#E69F00","#D55E00")

# species richness plot


errbar(x = 1:4,y = s.preds.median,yplus = s.preds.upper,yminus = s.preds.lower,
       col=errbar.cols,errbar.col = errbar.cols,ylim=c(min(s.preds.lower),max(s.preds.upper)),xaxt="n",
       ylab="Species richness (%)",xlab="Land use",bty="l")

axis(side = 1,at = 1:4,labels = c("PV","SV","AG.Low","AG.Hi"))

abline(h=0,col="#00000077",lty=2)

#title(main = "b.", adj = 0, cex.main = 1, line = 1)

p2 <- recordPlot()

# abundance plot

errbar(x = 1:4,y = a.preds.median,yplus = a.preds.upper,yminus = a.preds.lower,
       col=errbar.cols,errbar.col = errbar.cols,ylim=c(min(a.preds.lower),max(a.preds.upper)),xaxt="n",
       ylab="Total abundance (%)",xlab="Land use",bty="l")

axis(side = 1,at = 1:4,labels = c("PV","SV","AG.Low","AG.Hi"))

abline(h=0,col="#00000077",lty=2)

#title(main = "c.", adj = 0, cex.main = 1, line = 1)

p3 <- recordPlot()


#invisible(dev.off())

# Organise the plots

p4 <- plot_grid(p1, p2, p3, ncol = 1, scale = 0.85, labels = c("a", "b", "c"), label_x = 0.1)

# save figure
save_plot(paste0(outDir, "Figure_1.pdf"), p4, base_height = 8, base_width = 6)


# alternative fig 1 format

p5 <- plot_grid(p1, plot_grid(p2, p3, scale = 0.78, labels = c("b", "c"), label_x = 0.2), ncol = 1, rel_heights = c(0.7, 1), labels = c("a"), label_x = 0.1)

save_plot(paste0(outDir, "Figure_1_alt.pdf"), p5, base_height = 8, base_width = 8)



##%######################################################%##
#                                                          #
####                    Model stats                     ####
#                                                          #
##%######################################################%##


# models used = sm3 and am3

summary(sm3$model)

# rerun the model using GLMERSelect function to get stats

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

library(sjPlot)

tab_model(am3_2$model, transform = NULL, file = paste0(outDir, "/AbunLU_output_table.html"))
summary(am3_2$model)
R2GLMER(am3_2$model)

tab_model(sm3_2$model, transform = NULL, file = paste0(outDir, "/SRLU_output_table.html"))
summary(sm3_2$model)
R2GLMER(sm3_2$model)
