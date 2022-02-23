##%######################################################%##
#                                                          #
####        Organising and mapping climate data         ####
#                                                          #
##%######################################################%##

# This script explores the climate data, produces anomalies and associated maps.
# This script can also be used to generate datasets of the anomaly using different
# temperature thresholds. The threshold used in the paper is 10 degrees C. 


# directories
dataDir <- "0_data/"
outDir <- "3_PrepareClimateIndexMaps/"
#outDir <- "10_SCA_Baseline_testing/"
#outDir <- "14_Additional_Tests/"

if(!dir.exists(outDir)) dir.create(outDir)

sink(paste0(outDir,"log.txt"))

t.start <- Sys.time()

print(t.start)


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
library(snow)


# load in the mean temperature data from CRU
tmp <- stack(paste0(dataDir,"cru_ts4.03.1901.2018.tmp.dat.nc"),varname="tmp")

# take names of values for 1901 to 1930
tmp1901_1930 <- tmp[[names(tmp)[1:360]]]

# for testing baseline lengths:
#tmp1901_1905 <- tmp[[names(tmp)[1:60]]]
#tmp1901_1910 <- tmp[[names(tmp)[1:120]]]
#tmp1901_1920 <- tmp[[names(tmp)[1:240]]]


# present day temperature range (mid year 2005)
tmp2004_6 <- tmp[[names(tmp)[1237:1272]]]

# more recent data, used for maps in Extended Data Figure 7.
#tmp2016_18<- tmp[[names(tmp)[1381:1416]]]

# what is the active months threshold
thresh <- 10 # 6, 8, 10 degrees C

# which raster to use as present
pre_ras <- tmp2004_6 
#pre_ras <- tmp2016_18

# get vector of positions in values of raster that are not NA
ras <- pre_ras[[1]]
vals <- values(ras)

# create a set of points based on the non NA cells

wgs84 <- crs(tmp)

pnts <- rasterToPoints(ras, spatial = T)

SP <- SpatialPoints(pnts, proj4string=wgs84)  



nCores <- parallel::detectCores()

st1 <- Sys.time()

cl <- snow::makeCluster(nCores-1)

# export to clusters
snow::clusterExport(
  cl = cl,
  list = c('pre_ras', 'values', 'names', 'length', 'mean', 'sd',
           'tmp', 'SP','rasterize','crop','trim', 'grep', 'sapply', 'strsplit',
           'cellStats', 'thresh', 'tmp1901_1930', 'tmp1901_1905', 'tmp1901_1910', 'tmp1901_1920'),envir = environment())

temperatureVars <- data.frame(t(parSapply(
  cl = cl,X = (1:length(SP)),FUN = function(i){
    #cl = cl,X = (20208:20500),FUN = function(i){
      
      # #for testing
      # #anomVars <- NULL
      #  temperatureVars <-NULL
      # #for(i in 1:length(SP)){
      # for(i in 20208:20500){
      # print(i)
        
   # focus on this cell only, mask to improve speed
        
        ## Mask to improve speed
        mask <- trim(rasterize(SP[i, ], pre_ras[[1]]))
        mapCrop <- crop(pre_ras, mask)
        
  if(!length(names(mapCrop)[values(mapCrop) >= thresh]) == 0 & length(values(mapCrop)[!is.na(values(mapCrop))]) > 0 ){
          
# first identify insect active months for that cell. 
        
        # Get the average temperature for each month across 5 years
        vals <- NULL
        
        # for each month, get the average temp over the 5 years
        for(j in 1:12){
          
          if(j < 10){ mon <- paste0(0, j) }else {mon <- j}
          
          monthmean <- values(mean(mapCrop[[grep(mon, sapply(strsplit(names(mapCrop), "[.]"), "[[", 2))  ]]))
          
          vals <- rbind(vals, c(mon, monthmean))
          
        }
        
        vals <- as.data.frame(vals)
        vals$V2 <- as.numeric(as.character(vals$V2))
        
        # which months are the 5 year average >= the threshold
        vals <- vals[vals$V2 >= thresh, ]
        
        
        # sometimes the mean values don't make it over the threshold even if odd months are above
        if(nrow(vals) == 0){
          
          avg_temp = NA
          n_months = NA
          Anom <- NA
          StdAnom <- NA
          
          #temperatureVars <- rbind(temperatureVars, c(n_months = n_months, Anom = Anom, StdAnom = StdAnom))
          
          return(c(avg_temp = avg_temp, n_months = n_months, Anom = Anom, StdAnom = StdAnom))
          
          
        }else{
        
        
        
        # which are the good months
        months <- vals$V1
        
        # how many months are at or above the threshold?
        n_months <- length(months)
        
        # calculate the "present day" mean and sd
        avg_temp <- mean(vals$V2)

        ### now work out the baseline mean and sd for the active months ###
        
        # get the values for that grid cell across all years
        baseline <- crop(tmp1901_1930, mask)
        #baseline <- crop(tmp1901_1905, mask)
        #baseline <- crop(tmp1901_1910, mask)
        #baseline <- crop(tmp1901_1920, mask)
        
        # subset the baseline to just the required months
        baseline <-  baseline[[names(baseline)[sapply(strsplit(names(baseline), "[.]"), "[[", 2) %in% months]]]
        
        # get the mean and sd
        mean_baseline <- mean(values(baseline))
        sd_mean_baseline <- sd(values(baseline))
        
        
        # now calc the anomaly and standardised anomaly
        Anom <- avg_temp - mean_baseline
        StdAnom <-  Anom/sd_mean_baseline
        
        
        #temperatureVars <- rbind(temperatureVars, c(n_months = n_months, Anom = Anom, StdAnom = StdAnom))
        
        return(c(avg_temp = avg_temp, n_months = n_months, Anom = Anom, StdAnom = StdAnom))
        
  }}else{ # after 0/NA check
  
    avg_temp = NA
    n_months = NA
    Anom <- NA
    StdAnom <- NA
    
    #temperatureVars <- rbind(temperatureVars, c(n_months = n_months, Anom = Anom, StdAnom = StdAnom))
    
    return(c(avg_temp = avg_temp, n_months = n_months, Anom = Anom, StdAnom = StdAnom))
    
    
}

        
} # end of function




)))



snow::stopCluster(cl)

st2 <- Sys.time()

print(st2 - st1) # Time difference of 1.013413 hours

# save
#save(temperatureVars, file = paste0(outDir, "Map_data_tempvars_2004_06_thresh_ALLMonths.rdata"))
save(temperatureVars, file = paste0(outDir, "Map_data_tempvars_2004_06_thresh_", thresh, ".rdata"))
#save(temperatureVars, file = paste0(outDir, "Map_data_tempvars_2016_18_thresh_", thresh, ".rdata"))
#save(temperatureVars, file = paste0(outDir, "Map_data_tempvars_2004_06_thresh_", thresh, "_baseline5.rdata"))
#save(temperatureVars, file = paste0(outDir, "Map_data_tempvars_2004_06_thresh_", thresh, "_baseline10.rdata"))
#save(temperatureVars, file = paste0(outDir, "Map_data_tempvars_2004_06_thresh_", thresh, "_baseline20.rdata"))


#### take a look at correlations between the different metrics ####

# load(file = paste0(outDir, "Map_data_tempvars_2004_06.rdata"))
# load(file = paste0(outDir, "Map_data_tempvars_2016_18.rdata"))

# Using the 2004_06 data here

temperatureVars <- as.data.frame(temperatureVars) # 67420 rows

# remove the NAs
temp_data <- temperatureVars[!is.na(temperatureVars$avg_temp), ] # 58319 rows

cor(temp_data$avg_temp, temp_data$Anom) # -0.21, 2004-6 version, thresh 10

ggplot(data = temp_data, aes(x = avg_temp, y = Anom)) + 
  geom_point( size = 0.5) + 
  geom_smooth(method = "lm", size = 2) +
  theme_bw() + 
  xlab("Average temperature of the location\n (2004-2006, active months)") +
  ylab("Anomaly (difference between \npresent and baseline temperatures)") +
  ggtitle(paste0("Global values:\nPearson's correlation coefficient = ", round(cor(temp_data$avg_temp, temp_data$Anom), digits = 2)))

ggsave(filename = paste0(outDir, "Correlation_global_avgtemp_Anom.pdf"))

# remove the few outliers, some Infs
nrow(temp_data[temp_data$StdAnom >12, ]) # 6 rows
temp_data <- temp_data[temp_data$StdAnom < 12, ] # 58313 rows

cor(temp_data$avg_temp, temp_data$StdAnom) # -0.15, 2004-06, thresh 10

ggplot(data = temp_data, aes(x = avg_temp, y = StdAnom)) + 
  geom_point( size = 0.5) + 
  geom_smooth(method = "lm", size = 2) +
  theme_bw() + 
  xlab("Average temperature of the location\n (2004-2006, active months)") +
  ylab("Standardised climate anomaly") +
  ggtitle(paste0("Global values:\nPearson's correlation coefficient = ", round(cor(temp_data$avg_temp, temp_data$StdAnom), digits = 2)))

ggsave(filename = paste0(outDir, "Correlation_global_avgtemp_StdAnom.pdf"))

# remove STA values above 3
ggplot(data = temp_data[temp_data$StdAnom <=3, ], aes(x = avg_temp, y = StdAnom)) + 
  geom_point( size = 0.5) + 
  geom_smooth(method = "lm", size = 2) +
  theme_bw() + 
  xlab("Average temperature of the location\n (2004-2006, active months)") +
  ylab("Standardised climate anomaly") +
  ggtitle(paste0("Global values, outliers removed:\nPearson's correlation coefficient = ", round(cor(temp_data[temp_data$StdAnom <=3, 'avg_temp'], temp_data[temp_data$StdAnom <=3, 'StdAnom']), digits = 2)))

ggsave(filename = paste0(outDir, "Correlation_global_avgtemp_StdAnom_outliersrem.pdf"))


cor(temp_data$Anom, temp_data$StdAnom) # 0.36, 2004-06, thresh 10

# now take a look at the realms
# add in the lat/lons
# add data to the point lat/lons
SP_df <- as.data.frame(SP)

SP_df <- cbind(SP_df, temperatureVars)

# split by realm
SP_df$Tropical <- NA

SP_df[SP_df$y > -23.44 & SP_df$y < 23.44, 'Tropical'] <- "Tropical"
SP_df[is.na(SP_df$Tropical), 'Tropical'] <- "Temperate"

SP_df <- SP_df[!is.na(SP_df$avg_temp), ]

table(SP_df$Tropical)
# Temperate  Tropical 
# 40042     18277

cor_temp <- round(cor(SP_df[SP_df$Tropical == "Temperate", "avg_temp"], SP_df[SP_df$Tropical == "Temperate", "Anom"]), digits = 2) # 0.02
cor_trop <- round(cor(SP_df[SP_df$Tropical == "Tropical", "avg_temp"], SP_df[SP_df$Tropical == "Tropical", "Anom"]), digits = 2) # 0.14

SP_df$Tropical <- sub("Temperate", paste0("Temperate, cor = ", cor_temp), SP_df$Tropical)
SP_df$Tropical <- sub("Tropical", paste0("Tropical, cor = ", cor_trop), SP_df$Tropical)

ggplot(data = SP_df, aes(x = avg_temp, y = Anom)) + 
  geom_point( size = 0.5) + 
  geom_smooth(method = "lm", size = 2) +
  facet_wrap(~ Tropical) +
  theme_bw() + 
  xlab("Average temperature of the location\n (2016-2018, active months)") +
  ylab("Anomaly (difference between \npresent and baseline temperatures)") 
  
ggsave(filename = paste0(outDir, "Correlation_global_avgtemp_Anom_REALM.pdf"))


# split by realm
SP_df$Tropical <- NA

SP_df[SP_df$y > -23.44 & SP_df$y < 23.44, 'Tropical'] <- "Tropical"
SP_df[is.na(SP_df$Tropical), 'Tropical'] <- "Temperate"

# remove the few outliers
SP_df <- SP_df[SP_df$StdAnom < 12, ] # 58313 rows

cor_temp <- round(cor(SP_df[SP_df$Tropical == "Temperate", "avg_temp"], SP_df[SP_df$Tropical == "Temperate", "StdAnom"]), digits = 2) # -0.55
cor_trop <- round(cor(SP_df[SP_df$Tropical == "Tropical", "avg_temp"], SP_df[SP_df$Tropical == "Tropical", "StdAnom"]), digits = 2) # 0.08

SP_df$Tropical2 <- sub("Temperate", paste0("Temperate, cor = ", cor_temp), SP_df$Tropical)
SP_df$Tropical2 <- sub("Tropical", paste0("Tropical, cor = ", cor_trop), SP_df$Tropical)


ggplot(data = SP_df, aes(x = avg_temp, y = StdAnom)) + 
  geom_point( size = 0.5) + 
  geom_smooth(method = "lm", size = 2) +
  facet_wrap(~ Tropical2, scales = "free") +
  theme_bw() + 
  xlab("Average temperature of the location\n (2016-2018, active months)") +
  ylab("Standardised climate anomaly") 

ggsave(filename = paste0(outDir, "Correlation_global_avgtemp_StdAnom_REALM_1618.pdf"))


cor_temp <- round(cor(SP_df[SP_df$Tropical == "Temperate", "Anom"], SP_df[SP_df$Tropical == "Temperate", "StdAnom"]), digits = 2) # 0.37, 2004-6


##%######################################################%##
#                                                          #
####                 Manuscript Figures                 ####
#                                                          #
##%######################################################%##

# load the data for threshold and year range required
# load(file = paste0(outDir, "Map_data_tempvars_2004_06.rdata"))

# convert to dataframe
temperatureVars2 <- as.data.frame(temperatureVars)

# add data to the point lat/lons
SP_df <- as.data.frame(SP)

SP_df <- cbind(SP_df, temperatureVars2)

# quick look at plots of the global data
avg_temp_ras <- rasterFromXYZ(SP_df[ , 1:3])
Anom_ras <- rasterFromXYZ(SP_df[ , c(1,2,5)])
StdAnom_ras <- rasterFromXYZ(SP_df[ , c(1,2,6)])
n_months <- rasterFromXYZ(SP_df[ , c(1,2,4)])

plot(avg_temp_ras)
plot(Anom_ras)
plot(StdAnom_ras)
plot(n_months)



### first , the absolute change ###

# convert raster to dataframe
plot_data <- SP_df[, c(1,2,5)]

# organise breaks, colours and labels
brks <- c(-0.6,-0.2,-0.1,0,0.1,0.5,0.75,1,1.5,3,5)
cols <- c(rev(brewer.pal(n = 8,name = "Greens"))[5:8],
          (brewer.pal(n = 8,name = "Purples"))[4:6],
          (brewer.pal(n = 8,name = "Oranges"))[3:5])
labs <- c("-0.6 : -0.2","-0.2 : -0.1","-0.1 : 0",
          "0 : 0.1","0.1 : 0.5","0.5 : 0.75","0.75 : 1","1 : 1.5","1.5 : 3", ">3")

# assign values into bins
plot_data$bins <- cut(plot_data$Anom, 
                      breaks = brks, 
                      labels = labs,
                      include.lowest = TRUE)


# get worldmap for outline
world_map <- map_data("world")


# trying alternative legend
p1 <- ggplot(plot_data[!is.na(plot_data$Anom),]) + 
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), colour = "lightgrey", fill = "white",  size = 0.1) +
  geom_tile(aes(x = x, y = y, fill = bins), na.rm = TRUE) +
  scale_fill_manual(values = cols) + 
  xlab("") +
  ylab("") +
  labs(fill = "Absolute\nTemperature\nChange") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.3), 
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_blank(),
        #legend.key.width = unit(3, "cm"),
        axis.ticks = element_blank(), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 7), 
        legend.key.size = unit(0.2,"cm"), legend.direction = "vertical",
        title = element_text(size = 8, face = "bold")) +
  guides(fill = guide_legend(title.position = "top") ) + 
  ggtitle("a")




# get the mean climate value for each row of the dataset
rows <- init(Anom_ras, v='row')
ravg <- zonal(Anom_ras, rows, fun = 'mean', na.rm = T)
ravg[is.nan(ravg)] <- NA
ravg <- as.data.frame(ravg)

# plot the marginal plot
p2 <- ggplot(data = ravg) +
  geom_line( aes(x = zone, y = mean), col = c("#8379BD")) +
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
plot_data2 <- SP_df[, c(1,2,6)]

# organise breaks, colours and labels
brks2 <- c(-0.65,-0.2,-0.1,0,0.1,0.5,0.75,1,1.5,3, 5)
cols2 <- c(rev(brewer.pal(n = 8,name = "Greens"))[5:8],
           (brewer.pal(n = 8,name = "Purples"))[4:6],
           (brewer.pal(n = 8,name = "Oranges"))[3:5])
labs2 <- c("-0.6 : -0.2","-0.2 : -0.1","-0.1 : 0",
           "0 : 0.1","0.1 : 0.5","0.5 : 0.75","0.75 : 1","1 : 1.5","1.5 : 3", "> 3")

# assign values into bins
plot_data2$bins <- cut(plot_data2$StdAnom, 
                       breaks = brks2, 
                       labels = labs2,
                       include.lowest = TRUE)

plot_data2 <- plot_data2[!is.na(plot_data2$bins), ]


p3 <- ggplot(plot_data2[!is.na(plot_data2$StdAnom),]) + 
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), colour = "lightgrey", fill = "white", size = 0.1) +
  geom_raster(aes(x = x, y = y, fill = bins), na.rm = TRUE) +
  scale_fill_manual(values = cols2) + 
  xlab("") +
  ylab("") +
  labs(fill = "Standardised\nTemperature\nAnomaly") +
  theme_bw() +
  theme(legend.position = c(0.15, 0.3), 
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_blank(),
        #legend.key.width = unit(3, "cm"),
        axis.ticks = element_blank(), 
        legend.text = element_text(size = 6), 
        legend.title = element_text(size = 7), 
        legend.key.size = unit(0.2,"cm"), legend.direction = "vertical",
        title = element_text(size = 8, face = "bold")) +
  guides(fill = guide_legend(title.position = "top") ) + 
  ggtitle("b")

# get the mean climate value for each row of the dataset
rows2 <- init(StdAnom_ras, v='row')
ravg2 <- zonal(StdAnom_ras, rows2, fun = 'mean', na.rm = T)
ravg2[is.nan(ravg2)] <- NA
ravg2 <- as.data.frame(ravg2)

# remove the 3 infinite values
ravg2 <- ravg2[!ravg2$mean == Inf, ]

# plot the marginal plot
p4 <- ggplot(data = ravg2) +
  geom_line( aes(x = zone, y = mean), col = c("#8379BD")) +
  geom_ribbon(aes(ymin = min(ravg2$mean, na.rm = T), ymax = mean, x = zone), fill = c("#473C8B"), alpha = 0.7) +
  theme_bw() + 
  scale_x_reverse(limits = c(300, 1), expand = c(0,0)) +
  scale_y_continuous(limits = c(min(ravg2$mean, na.rm = T), max(ravg2$mean, na.rm = T)), expand = c(0,0)) +
  theme(panel.border = element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks = element_blank()
  ) + 
  coord_flip()


# organise the plots and legends into one object


final_plot <- cowplot::plot_grid(
  #cowplot::plot_grid(
    cowplot::plot_grid(
      p1 
      , p2
      , nrow = 1
      , align = "hv"
      , rel_widths = c(3,1)),

  cowplot::plot_grid(
      p3 
      , p4
      , nrow = 1
      , align = "hv"
      , rel_widths = c(3,1)),

  nrow = 2
)

# save as a pdf
ggsave(filename = paste0(outDir, "Extended_Data1_maps_thresh_", thresh, ".pdf"), plot = last_plot(), width = 183, height = 200, units = "mm", dpi = 300)
ggsave(filename = paste0(outDir, "Extended_Data1_maps_thresh_", thresh, ".jpeg"), plot = last_plot(), width = 183, height = 200, units = "mm", dpi = 300)



#### Figure - map of number of active months ####

map_data <- SP_df[, c(1,2,4)]

ggplot(map_data[!is.na(map_data$n_months),]) + 
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), colour = "lightgrey", fill = "white", size = 0.1) +
  geom_raster(aes(x = x, y = y, fill = n_months), na.rm = TRUE) +
  scale_fill_gradient(low = c("#79CDCD"), high = c("#00688B"), na.value = NA, limits = c(1, 12)) +
  xlab("") +
  ylab("") +
  labs(fill = "Number of months\nabove 10 degrees C") +
  theme_bw() +
  theme(legend.position = 'bottom', 
        panel.border = element_blank(), 
        panel.grid = element_blank(),
        axis.text = element_blank(),
        axis.ticks = element_blank(), 
        legend.text = element_text(size = 8), 
        legend.title = element_text(size = 10)) +
  guides(colour = guide_colourbar(show.limits = TRUE)) 


ggsave(filename = paste0(outDir, "Nmonths_plot_thresh_", thresh, ".pdf"), width = 4, height = 3)


t.end <- Sys.time()

print(round(t.end - t.start,0))

sink()




