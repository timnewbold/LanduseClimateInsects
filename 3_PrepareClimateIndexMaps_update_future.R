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
library(snow)


# directories
dataDir <- "0_data/"
outDir <- "3_PrepareClimateIndexMaps/"

# load in the mean temperature data from CRU
tmp <- stack(paste0(dataDir,"cru_ts4.03.1901.2018.tmp.dat.nc"),varname="tmp")

# take names of values for 1901 to 1930
tmp1901_1930 <- tmp[[names(tmp)[1:360]]]

# take the current predicts time data
tmp2004_6 <- tmp[[names(tmp)[1237:1272]]]

# more recent data
tmp2016_18<- tmp[[names(tmp)[1381:1416]]]

# what is the active months threshold
thresh <- 10 # 6, 8, 10

# create a function that can be run in parallel to produce maps

### getting the future temperature anomalies to be added to the present day values ###

# future temperature estimates from ISIMIP
futureClimateDir <- "0_data/ISIMIPAnomalies"

# selection of years
years <- 2069:2071

# file path for ISIMIP data
all.files <- dir(path = futureClimateDir,recursive = TRUE,full.names = TRUE)

# Behrman projection
behrCRS <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')

# using RCP 8.5
mean.anom.2069.2071 <- stack(lapply(X = years,FUN = function(yr){
  
  print(yr)
  
  all.model.files <- all.files[grepl("rcp85",all.files) & grepl(yr,all.files)]
  
  # Check that there are the same files for each scenario-year combination
  stopifnot(all(sapply(
    X = gsub("0_data/ISIMIPAnomalies/","",all.model.files),function(f) return(strsplit(x = f,split = "[-_]",fixed = FALSE)[[1]][1]))==
      c("GFDL","HadGEM2","IPSL","MIROC5")))
  
  # what are each of these files? One for each model, taking an average across the 4.
  
  #
  meant.anom <- mean(stack(lapply(X = all.model.files,function(f){
    
    ras <- stack(f)$"X0.1"
    
  })),na.rm=TRUE)/10
  
  #meant <- hist.mean.temp.1979.2013 + (meant.anom/10)
  
  return(meant.anom)
  
}))

# mean across years, these are the values to add to each month of the present day (2016-18 temps)
mean.anom.2069.2071 <- stackApply(x = mean.anom.2069.2071,indices = rep(1,3),fun = mean)

plot(mean.anom.2069.2071)


# which raster to use as present
# pre_ras <- tmp2004_6 
pre_ras <- tmp2016_18

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
           'cellStats', 'thresh', 'tmp1901_1930', 'mean.anom.2069.2071'),envir = environment())

temperatureVars <- data.frame(t(parSapply(
  cl = cl,X = (1:length(SP)),FUN = function(i){
    #cl = cl,X = (20208:20500),FUN = function(i){
    
    # #for testing
    #anomVars <- NULL
    #temperatureVars <-NULL
    #for(i in 1:length(SP)){
    #for(i in 20208:20500){
    #print(i)
    
    # focus on this cell only, mask to improve speed
    
    ## Mask to improve speed
    mask <- trim(rasterize(SP[i, ], pre_ras[[1]]))
    mapCrop <- crop(pre_ras, mask)
    
    
    ### add the future anomaly values to the present day months to get the future temperatures ###
    mean.anom.2069.2071_crop <- crop(mean.anom.2069.2071, mask)
    
    mapCrop_fut <- mapCrop + mean.anom.2069.2071_crop
    names(mapCrop_fut) <- names(mapCrop)
    
    ### then select those that meet the threshold and calculate the anomaly values ###

    
    if(!length(names(mapCrop_fut)[values(mapCrop_fut) >= thresh]) == 0 & length(values(mapCrop_fut)[!is.na(values(mapCrop_fut))]) > 0 ){
      
      # first identify insect active months for that cell. 
      
      # Get the average temperature for each month across 5 years
      vals <- NULL
      
      # for each month, get the average temp over the 5 years
      for(j in 1:12){
        
        if(j < 10){ mon <- paste0(0, j) }else {mon <- j}
        
        monthmean <- values(mean(mapCrop_fut[[grep(mon, sapply(strsplit(names(mapCrop_fut), "[.]"), "[[", 2))  ]]))
        
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
save(temperatureVars, file = paste0(outDir, "Map_data_tempvars_2016_18_thresh_", thresh, "_FUTURE.rdata"))


# load(file = paste0(outDir, "Map_data_tempvars_2016_18_thresh_", thresh, "_FUTURE.rdata"))

# 25 rows have infs from where the sd from baseline with 1 month of data is 0.




##%################################################################%##
#                                                                    #
####                 Manuscript Figure Future Map                 ####
#                                                                    #
##%################################################################%##



# convert to dataframe
temperatureVars2 <- as.data.frame(temperatureVars)

# add data to the point lat/lons
SP_df <- as.data.frame(SP)

SP_df <- cbind(SP_df, temperatureVars2)


# load in the 2017 data
load(file = paste0(outDir, "Map_data_tempvars_2016_18.rdata"))

# add data to the point lat/lons
SP_df_2018 <- as.data.frame(SP)

SP_df_2018 <- cbind(SP_df_2018, temperatureVars)

# add year labels
SP_df$year <- 2070
SP_df_2018$year <- 2018

plot_data <- rbind(SP_df, SP_df_2018)
#### now the standardised anomaly ####

# convert raster to dataframe
plot_data2 <- plot_data[, c(1,2,6,7)]

# organise breaks, colours and labels
brks <- c(-0.5,-0.2,-0.1,0,0.1,0.5,0.75,1,1.5,2,5,10,50)
cols <- c(rev(brewer.pal(n = 8,name = "Greens"))[5:8],
          (brewer.pal(n = 8,name = "Purples"))[4:6],
          (brewer.pal(n = 9,name = "Oranges"))[5:9])
labs <- c("-0.5 : -0.2","-0.2 : -0.1","-0.1 : 0",
          "0 : 0.1","0.1 : 0.5","0.5 : 0.75","0.75 : 1","1 : 1.5","1.5 : 2","2 : 5","5 : 10", "> 10")

# assign values into bins
plot_data2$bins <- cut(plot_data2$StdAnom, 
                     breaks = brks, 
                     labels = labs,
                     include.lowest = TRUE)


plot_data2 <- plot_data2[!is.na(plot_data2$bins), ]

# plot the raster
p3 <- ggplot(plot_data2[!is.na(plot_data2$StdAnom),]) + 
  geom_polygon(data = world_map, aes(x = long, y = lat, group = group), colour = "lightgrey", fill = "white", size = 0.1) +
  geom_raster(aes(x = x, y = y, fill = bins), na.rm = TRUE) +
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



# save as a pdf
ggsave(filename = paste0(outDir, "/Figure4_Current_Future_panel_plot.pdf"), width = 8, height = 4)

