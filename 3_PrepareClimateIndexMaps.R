##%######################################################%##
#                                                          #
####        Organising and mapping climate data         ####
#                                                          #
##%######################################################%##

# This script explores the climate data, produces anomalies and associated maps.

# load required libraries
library(raster)
library(sp)
library(dismo)
library(rgdal)
library(RColorBrewer)
library(ncdf4)
library(rasterVis)
library(gridExtra)
library(cowplot)
library(ggplot2)
library(viridis)

# directories
dataDir <- "0_data/"
outDir <- "3_PrepareClimateIndexMaps/"

# load in the mean temperature data from CRU
tmp <- stack(paste0(dataDir,"cru_ts4.03.1901.2018.tmp.dat.nc"),varname="tmp")

# take names of values for 1901 to 1930
tmp1901_1930 <- tmp[[names(tmp)[1:360]]]

# calculate the mean and sd of the baseline values
tmp1901_1930mean <- calc(tmp1901_1930, base::mean)
tmp1901_1930sd <- calc(tmp1901_1930, stats::sd)

### Note: 2005 is the mean year for insect data

par(mfrow=c(1,1))
par(oma=c(5,5,5,5))

##### 2004 to 2006 maps #####

tmp2004_6 <- tmp[[names(tmp)[1237:1272]]]


### Calculate the standardised anomaly ###

# calc the mean for present time period
tmp2004_6mean <- calc(tmp[[names(tmp)[1237:1272]]], base::mean)

# calc mean for baseline
tmp2004_6_climate_anomaly <- (calc(tmp2004_6, base::mean)-tmp1901_1930mean)

# standardise the baseline
tmp2004_6std_climate_anomaly <- (calc(tmp2004_6, base::mean)-tmp1901_1930mean)  / tmp1901_1930sd

# info for map
breaks <- c(-2,0,0.25, 0.5,0.75, 1,1.5, 2,2.5, 3,4,5, 10)
pallete <- colorRampPalette(c("lightblue","red", "black"))

# plot map
plot(tmp2004_6std_climate_anomaly, breaks=breaks, col=pallete(12), main="Mean Standardized Climate Anomaly 2004 to 2006")



## Climate anomaly (non standardised, for comparison) ##
warming2004_6 <- calc(tmp2004_6, base::mean)-tmp1901_1905mean
breaks2 <- c(0,0.25,0.5,0.75, 1 ,1.5,2,2.5,3, 4)
length(breaks2)
pallete2 <- colorRampPalette(c("lightblue","red", "black"))
plot(warming2004_6, breaks=breaks2, col=pallete2(11), main="Mean warming 2004 to 2006 Since 1901")



## anomaly in 1970 ##

names(tmp)[829:840] ##1970
tmp1970sd <-  calc(tmp[[names(tmp)[829:840]]], stats::sd)
tmp1970mean <-  calc(tmp[[names(tmp)[829:840]]], base::mean)
tmp1970std_climate_anomaly <- (tmp1970mean-tmp1901_1930mean)  / tmp1901_1930sd
breaks <- c(-2,0,0.25, 0.5,0.75, 1,1.5, 2,2.5, 3,4,5)
pallete <- colorRampPalette(c("lightblue","red", "black"))
plot(tmp1970std_climate_anomaly, breaks=breaks, col=pallete(12), main="1970")



## anomaly, present 2016 to 2018 ##

names(tmp)[1381:1416] ## 2016 to 2018 
tmp2016_18sd <-  calc(tmp[[names(tmp)[1381:1416]]], stats::sd)
tmp2016_18mean <-  calc(tmp[[names(tmp)[1381:1416]]], base::mean)
tmp2016_18std_climate_anomaly_2 <-  (tmp2016_18mean - tmp1901_1930mean) / tmp2016_18sd
tmp2016_18std_climate_anomaly <- (tmp2016_18mean-tmp1901_1930mean)  / tmp1901_1930sd
breaks <- c(-2,0,0.25, 0.5,0.75, 1,1.5, 2,2.5, 3,4,5, 15)
pallete <- colorRampPalette(c("lightblue","red", "black"))
plot(tmp2016_18std_climate_anomaly, breaks=breaks, col=pallete(12), main="2016_18")


#### Calculating future anomalies ####

months.1979.2013 <- 937:1356

hist.mean.temp.1979.2013 <- stack(stackApply(x = tmp[[months.1979.2013]],
                                       indices = (rep(1:35,each=12)),fun = mean))
hist.mean.temp.1979.2013 <- stackApply(x = hist.mean.temp.1979.2013,indices = rep(1,35),
                                       fun = mean)

# future temperature estimates from ISIMIP
futureClimateDir <- "0_data/ISIMIPAnomalies"

# selection of years
years <- 2069:2071

# file path for ISIMIP data
all.files <- dir(path = futureClimateDir,recursive = TRUE,full.names = TRUE)

# Behrman projection 
behrCRS <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')

# using RCP 8.5
mean.temp.2069.2071 <- stack(lapply(X = years,FUN = function(yr){
  
  print(yr)
  
  all.model.files <- all.files[grepl("rcp85",all.files) & grepl(yr,all.files)]
  
  # Check that there are the same files for each scenario-year combination
  stopifnot(all(sapply(
    X = gsub("0_data/ISIMIPAnomalies/","",all.model.files),function(f) return(strsplit(x = f,split = "[-_]",fixed = FALSE)[[1]][1]))==
      c("GFDL","HadGEM2","IPSL","MIROC5")))
  
  # what are each of these files?
  
  # 
  meant.anom <- mean(stack(lapply(X = all.model.files,function(f){
    
    ras <- stack(f)$"X0.1"
    
  })),na.rm=TRUE)
  
  meant <- hist.mean.temp.1979.2013 + (meant.anom/10)
  
  return(meant)
  
}))

mean.temp.2069.2071 <- stackApply(x = mean.temp.2069.2071,indices = rep(1,3),fun = mean)

# calc the anomalies for the future years
tmp2069_71_climate_anomaly <- (mean.temp.2069.2071-tmp1901_1930mean)
tmp2069_71std_climate_anomaly <- (mean.temp.2069.2071-tmp1901_1930mean)  / tmp1901_1930sd

brks <- c(-50,-5,-2,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,0.75,1,1.5,2,5,50)
brks2 <- c(-1,-0.75,-0.5,-0.25,-0.1,0,0.1,0.25,0.5,0.75,1,1.5,3.1)
cols <- c(rev(brewer.pal(n = 9,name = "Greens"))[3:9],
          (brewer.pal(n = 9,name = "Purples"))[4:8],
          (brewer.pal(n = 9,name = "Oranges"))[6:9])
cols2 <- c(rev(brewer.pal(n = 9,name = "Blues"))[5:9],
          (brewer.pal(n = 9,name = "Reds"))[3:9])

pdf(file = paste0(outDir,"ClimateIndexMaps.pdf"),width = 8.5/2.54,height = 14.77/2.54)

par(ps=10)
par(cex=1)
par(cex.lab=1)
par(cex.axis=1)
par(cex.main=1)

layout.mat <- matrix(data = 1:5,nrow = 5,ncol = 1,byrow = TRUE)
layout(mat = layout.mat,widths = 8.5,heights = c(2,3.59,3.59,3.59,2))

# par(mfrow=c(3,1))

par(mar=c(0,0,0,0))

# image(tmp1970std_climate_anomaly,breaks=brks,col=cols,xaxt="n",yaxt="n",bty="n")
plot.new()
legend(0,1,c("< -0.75","-0.75 : -0.5","-0.5 : -0.25","-0.25 : -0.1","-0.1 : 0",
       "0 : 0.1","0.1 : 0.25","0.25 : 0.5","0.5 : 0.75","0.75 : 1","1 : 1.5","> 1.5"),
       ncol=4,fill=cols2,bty="n",cex=1)
image(tmp2004_6_climate_anomaly,breaks=brks2,col=cols2,xaxt="n",yaxt="n",bty="n",
      xlim=c(-180,180),ylim=c(-68,84))
text(-175,-20,"2005\n(Absolute)",pos=4)
image(tmp2004_6std_climate_anomaly,breaks=brks,col=cols,xaxt="n",yaxt="n",bty="n")
text(-175,-20,"2005\n(Standardised)",pos=4)
image(tmp2069_71std_climate_anomaly,breaks=brks,col=cols,xaxt="n",yaxt="n",bty="n")
text(-175,-20,"2070 - RCP 8.5\n(Standardised)",pos=4)
plot.new()
legend(0,1,c("< -5","-5 : -2","-2 : -1","-1 : -0.5","-0.5 : -0.2","-0.2 : -0.1","-0.1 : 0",
             "0 : 0.1","0.1 : 0.2","0.2 : 0.5","0.5 : 0.75","0.75 : 1","1 : 1.5","1.5 : 2","2 : 5","> 5"),
       ncol=4,fill=cols,bty="n",cex=1)
# image(tmp2016_18std_climate_anomaly,breaks=brks,col=cols,xaxt="n",yaxt="n",bty="n")

invisible(dev.off())



##%######################################################%##
#                                                          #
####                 Manuscript Figures                 ####
#                                                          #
##%######################################################%##


### first , the absolute change ###

# convert raster to dataframe
plot_data <- as.data.frame(tmp2004_6_climate_anomaly, xy = TRUE)

# organise breaks, colours and labels
brks <- c(-0.5,-0.2,-0.1,0,0.1,0.5,0.75,1,1.5,3,5)
cols <- c(rev(brewer.pal(n = 8,name = "Greens"))[5:6],
          (brewer.pal(n = 8,name = "Purples"))[4:8],
          (brewer.pal(n = 8,name = "Oranges"))[3:5])
labs <- c("-0.5 : -0.2","-0.2 : -0.1","-0.1 : 0",
          "0 : 0.1","0.1 : 0.5","0.5 : 0.75","0.75 : 1","1 : 1.5","1.5 : 3", "3 : 5")

# assign values into bins
plot_data$bins <- cut(plot_data$layer, 
                     breaks = brks, 
                     labels = labs,
                     include.lowest = TRUE)

# plot the raster
p1 <- ggplot(plot_data[!is.na(plot_data$layer),]) + 
          geom_raster(aes(x = x, y = y, fill = bins), na.rm = TRUE) +
          scale_fill_manual(values = cols) + 
          xlab("") +
          ylab("") +
          labs(fill = "Absolute Temperature Change") +
          theme_bw() +
          theme(legend.position = 'bottom', 
                panel.border = element_blank(), 
                panel.grid = element_blank(),
                axis.text = element_blank(),
                #legend.key.width = unit(3, "cm"),
                axis.ticks = element_blank(), 
                legend.text = element_text(size = 10), 
                legend.title = element_text(size = 12), 
                legend.key.size = unit(0.2,"cm")) +
          #guides(fill = guide_colorbar(title.position = "top")) + 
  ggtitle("a.")

# get the mean climate value for each row of the dataset
rows <- init(tmp2004_6_climate_anomaly, v='row')
ravg <- zonal(tmp2004_6_climate_anomaly, rows, fun = 'mean', na.rm = T)
ravg <- as.data.frame(ravg)

# plot the marginal plot
p2 <- ggplot(data = ravg) +
        geom_line( aes(x = zone, y = mean), col = c("#473C8B")) +
        geom_ribbon(aes(ymin = min(ravg$mean, na.rm = T), ymax = mean, x = zone), fill = c("#473C8B"), alpha = 0.7) +
        theme_bw() + 
        scale_x_reverse(limits = c(300, 1), expand = c(0,0)) +
        scale_y_continuous(limits = c(min(ravg$mean, na.rm = T), max(ravg$mean, na.rm = T)), expand = c(0,0)) +
        theme(panel.border = element_blank(), 
              panel.grid = element_blank(),
              axis.text = element_blank(),
              axis.title = element_blank(),
              axis.ticks = element_blank()
              ) + 
              coord_flip()



#### now the standardised anomaly ####

# convert raster to dataframe
plot_data2 <- as.data.frame(tmp2004_6std_climate_anomaly, xy = TRUE)

# organise breaks, colours and labels
brks2 <- c(-0.5,-0.2,-0.1,0,0.1,0.5,0.75,1,1.5,3, 5, 10.4)
cols2 <- c(rev(brewer.pal(n = 8,name = "Greens"))[5:6],
          (brewer.pal(n = 8,name = "Purples"))[4:8],
          (brewer.pal(n = 8,name = "Oranges"))[3:6])
labs2 <- c("-0.5 : -0.2","-0.2 : -0.1","-0.1 : 0",
          "0 : 0.1","0.1 : 0.5","0.5 : 0.75","0.75 : 1","1 : 1.5","1.5 : 3", "3 : 5", "> 5")

# assign values into bins
plot_data2$bins <- cut(plot_data2$layer, 
                      breaks = brks2, 
                      labels = labs2,
                      include.lowest = TRUE)

# plot the raster
p3 <- ggplot(plot_data2[!is.na(plot_data2$layer),]) + 
  geom_raster(aes(x = x, y = y, fill = bins), na.rm = TRUE) +
  scale_fill_manual(values = cols2) + 
  xlab("") +
  ylab("") +
  labs(fill = "Standardised Temperature Anomaly") +
  theme_bw() +
  theme(legend.position = 'bottom', 
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_blank(),
        #legend.key.width = unit(3, "cm"),
        axis.ticks = element_blank(), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 12), 
        legend.key.size = unit(0.2,"cm")) +
  #guides(fill = guide_colorbar(title.position = "top")) +
  ggtitle("b.")

# get the mean climate value for each row of the dataset
rows2 <- init(tmp2004_6std_climate_anomaly, v='row')
ravg2 <- zonal(tmp2004_6std_climate_anomaly, rows2, fun = 'mean', na.rm = T)
ravg2 <- as.data.frame(ravg2)

# plot the marginal plot
p4 <- ggplot(data = ravg2) +
  geom_line( aes(x = zone, y = mean), col = c("#473C8B")) +
  geom_ribbon(aes(ymin = min(ravg$mean, na.rm = T), ymax = mean, x = zone), fill = c("#473C8B"), alpha = 0.7) +
  theme_bw() + 
  scale_x_reverse(limits = c(300, 1), expand = c(0,0)) +
  scale_y_continuous(limits = c(min(ravg$mean, na.rm = T), max(ravg$mean, na.rm = T)), expand = c(0,0)) +
  theme(panel.border = element_blank(), 
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()
    ) + 
  coord_flip()



# setting up a black plot to fill the gap
p0 <- ggplot() +
  theme_bw() +
  theme(panel.border = element_blank(), 
    panel.grid = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    axis.ticks = element_blank()) 
  


# organise the plots and legends into one object


final_plot <- cowplot::plot_grid(
  cowplot::plot_grid(
    cowplot::plot_grid(
    p1 + theme(legend.position = "none")
    , p2
    , nrow = 1
    , align = "hv"
    , rel_widths = c(3,1))
  , cowplot::plot_grid(
    get_legend(p1)
    #, p0
    #, ncol = 2
    , rel_widths = c(2,1))

 , nrow = 2
 , rel_heights = c(3, 1)
),
cowplot::plot_grid(
  cowplot::plot_grid(
    p3 + theme(legend.position = "none")
    , p4
    , nrow = 1
    , align = "hv"
    , rel_widths = c(3,1))
  , cowplot::plot_grid(
    get_legend(p3)
  #  , NULL
  #  , ncol = 2
    , rel_widths = c(1,1))
  
  , nrow = 2
  , rel_heights = c(3, 1)
  
),
nrow = 2
)

# save as a pdf
ggsave(filename = paste0(outDir, "/Extended_Data1_maps.pdf"), plot = last_plot(), width = 9, height = 10)

# save final_plot as an rdata file to be used in later scripts
save(final_plot, file = paste0(outDir, "/abs_and_anom_maps.rdata"))





##### Manuscript, Figure 4, present and future anomaly #####

plot_data_pres <- as.data.frame(tmp2016_18std_climate_anomaly, xy = TRUE)
plot_data_fut <- as.data.frame(tmp2069_71std_climate_anomaly, xy = TRUE)

# add years to data table
plot_data_pres$year <- 2018
plot_data_fut$year <- 2070

# combine the two together
all_plot <- rbind(plot_data_pres, plot_data_fut)

all_plot <- all_plot[!is.na(all_plot$layer), ]

# organise breaks, colours and labels
brks <- c(-0.5,-0.2,-0.1,0,0.1,0.5,0.75,1,1.5,2,5,10,50)
cols <- c(rev(brewer.pal(n = 8,name = "Greens"))[5:8],
          (brewer.pal(n = 8,name = "Purples"))[4:6],
          (brewer.pal(n = 9,name = "Oranges"))[5:9])
labs <- c("-0.5 : -0.2","-0.2 : -0.1","-0.1 : 0",
          "0 : 0.1","0.1 : 0.5","0.5 : 0.75","0.75 : 1","1 : 1.5","1.5 : 2","2 : 5","5 : 10", "> 10")

# assign values into bins
all_plot$bins <- cut(all_plot$layer, 
                     breaks = brks, 
                     labels = labs,
                     include.lowest = TRUE)

# plot
ggplot(all_plot) + 
  geom_raster(aes(x = x, y = y, fill = bins), alpha = 0.9) +
  #scale_fill_viridis_c(option = "magma", values = c(0, 0.2, 1)) + 
  scale_fill_manual(values = cols) +
  facet_grid(~ year) +
  xlab("") +
  ylab("") +
  labs(fill = "Standardised\nTemperature\nAnomaly") +
  theme_bw() +
  theme(legend.position = 'bottom', 
        #panel.border = element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_blank(),
        #legend.key.width = unit(3, "cm"),
        axis.ticks = element_blank(), 
        legend.text = element_text(size = 10), 
        legend.title = element_text(size = 10), 
        #legend.key.size = unit(0.2,"cm",
        panel.border = element_rect(size = 0.2),
        strip.background = element_rect(size = 0.2, fill = "transparent"),
        strip.text = element_text(size = 10))
 

# save plot as pdf
ggsave(filename = paste0(outDir, "/Figure4_Current_Future_panel_plot.pdf"), width = 8, height = 4)





