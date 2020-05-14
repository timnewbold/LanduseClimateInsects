###Map of standardized anomaly

##Packages
library(raster)
library(sp)
library(dismo)
library(rgdal)
library(RColorBrewer)
library(ncdf4)

dataDir <- "0_data/"

outDir <- "3_PrepareClimateIndexMaps/"

tmp <- stack(paste0(dataDir,"cru_ts4.03.1901.2018.tmp.dat.nc"),varname="tmp")

tmp1901_1905 <- tmp[[names(tmp)[1:60]]]

tmp1901_1905mean <- calc(tmp1901_1905, base::mean)
tmp1901_1905sd <- calc(tmp1901_1905, stats::sd)

###2005 is mean year for insect data

par(mfrow=c(1,1))
par(oma=c(5,5,5,5))

###2004 to 2006 maps
tmp2004_6 <- tmp[[names(tmp)[1237:1272]]]

#Standardised climate anomaly
tmp2004_6mean <-  calc(tmp[[names(tmp)[1237:1272]]], base::mean)
tmp2004_6_climate_anomaly <- (calc(tmp2004_6, base::mean)-tmp1901_1905mean)
tmp2004_6std_climate_anomaly <- (calc(tmp2004_6, base::mean)-tmp1901_1905mean)  / tmp1901_1905sd
breaks <- c(-2,0,0.25, 0.5,0.75, 1,1.5, 2,2.5, 3,4,5, 10)
pallete <- colorRampPalette(c("lightblue","red", "black"))
plot(tmp2004_6std_climate_anomaly, breaks=breaks, col=pallete(12), main="Mean Standardized Climate Anomaly 2004 to 2006")

#Climate anomaly (non standardised, for comparison)
warming2004_6 <- calc(tmp2004_6, base::mean)-tmp1901_1905mean
breaks2 <- c(0,0.25,0.5,0.75, 1 ,1.5,2,2.5,3, 4)
length(breaks2)
pallete2 <- colorRampPalette(c("lightblue","red", "black"))
plot(warming2004_6, breaks=breaks2, col=pallete2(11), main="Mean warming 2004 to 2006 Since 1901")

##1970

names(tmp)[829:840] ##1970
tmp1970sd <-  calc(tmp[[names(tmp)[829:840]]], stats::sd)
tmp1970mean <-  calc(tmp[[names(tmp)[829:840]]], base::mean)
tmp1970std_climate_anomaly <- (tmp1970mean-tmp1901_1905mean)  / tmp1901_1905sd
breaks <- c(-2,0,0.25, 0.5,0.75, 1,1.5, 2,2.5, 3,4,5)
pallete <- colorRampPalette(c("lightblue","red", "black"))
plot(tmp1970std_climate_anomaly, breaks=breaks, col=pallete(12), main="1970")



##2016 to 2018

names(tmp)[1381:1416] ## 2016 to 2018 
tmp2016_18sd <-  calc(tmp[[names(tmp)[1381:1416]]], stats::sd)
tmp2016_18mean <-  calc(tmp[[names(tmp)[1381:1416]]], base::mean)
tmp2016_18std_climate_anomaly_2 <-  (tmp2016_18mean - tmp1901_1905mean) / tmp2016_18sd
tmp2016_18std_climate_anomaly <- (tmp2016_18mean-tmp1901_1905mean)  / tmp1901_1905sd
breaks <- c(-2,0,0.25, 0.5,0.75, 1,1.5, 2,2.5, 3,4,5, 15)
pallete <- colorRampPalette(c("lightblue","red", "black"))
plot(tmp2016_18std_climate_anomaly, breaks=breaks, col=pallete(12), main="2016_18")

months.1979.2013 <- 937:1356

hist.mean.temp.1979.2013 <- stack(stackApply(x = tmp[[months.1979.2013]],
                                       indices = (rep(1:35,each=12)),fun = mean))
hist.mean.temp.1979.2013 <- stackApply(x = hist.mean.temp.1979.2013,indices = rep(1,35),
                                       fun = mean)

futureClimateDir <- "0_data/ISIMIPAnomalies/"

years <- 2069:2071

all.files <- dir(path = futureClimateDir,recursive = TRUE,full.names = TRUE)

behrCRS <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')

mean.temp.2069.2071 <- stack(lapply(X = years,FUN = function(yr){
  
  print(yr)
  
  all.model.files <- all.files[grepl("rcp85",all.files) & grepl(yr,all.files)]
  
  # Check that there are the same files for each scenario-year combination
  stopifnot(all(sapply(
    X = gsub("0_data/ISIMIPAnomalies/","",all.model.files),function(f) return(strsplit(x = f,split = "[-_]",fixed = FALSE)[[1]][1]))==
      c("GFDL","HadGEM2","IPSL","MIROC5")))
  
  meant.anom <- mean(stack(lapply(X = all.model.files,function(f){
    
    ras <- stack(f)$"X0.1"
    
  })),na.rm=TRUE)
  
  meant <- hist.mean.temp.1979.2013 + (meant.anom/10)
  
  return(meant)
  
}))

mean.temp.2069.2071 <- stackApply(x = mean.temp.2069.2071,indices = rep(1,3),fun = mean)

tmp2069_71_climate_anomaly <- (mean.temp.2069.2071-tmp1901_1905mean)
tmp2069_71std_climate_anomaly <- (mean.temp.2069.2071-tmp1901_1905mean)  / tmp1901_1905sd

brks <- c(-50,-5,-2,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,2,5,50)
brks2 <- c(-1,-0.75,-0.5,-0.25,-0.1,0,0.1,0.25,0.5,0.75,1,1.5,3.1)
cols <- c(rev(brewer.pal(n = 9,name = "Greens"))[3:9],
          (brewer.pal(n = 9,name = "Purples"))[3:6],
          (brewer.pal(n = 9,name = "Oranges"))[7:9])
cols2 <- c(rev(brewer.pal(n = 9,name = "Blues"))[5:9],
          (brewer.pal(n = 9,name = "Reds"))[3:9])

pdf(file = paste0(outDir,"ClimateIndexMaps.pdf"),width = 12.5/2.54,height = 19.84/2.54)

layout.mat <- matrix(data = 1:5,nrow = 5,ncol = 1,byrow = TRUE)
layout(mat = layout.mat,widths = 17.5,heights = c(2,5.28,5.28,5.28,2))

# par(mfrow=c(3,1))

par(mar=c(0,0,0,0))

# image(tmp1970std_climate_anomaly,breaks=brks,col=cols,xaxt="n",yaxt="n",bty="n")
plot.new()
legend(0.1,1,c("< -0.75","-0.75 : -0.5","-0.5 : -0.25","-0.25 : -0.1","-0.1 : 0",
       "0 : 0.1","0.1 : 0.25","0.25 : 0.5","0.5 : 0.75","0.75 : 1","1 : 1.5","> 1.5"),
       ncol=5,fill=cols2,bty="n")
image(tmp2004_6_climate_anomaly,breaks=brks2,col=cols2,xaxt="n",yaxt="n",bty="n",
      xlim=c(-180,180),ylim=c(-68,84))
text(-175,-20,"2005 (Absolute)",pos=4)
image(tmp2004_6std_climate_anomaly,breaks=brks,col=cols,xaxt="n",yaxt="n",bty="n")
text(-175,-20,"2005 (Standardised)",pos=4)
image(tmp2069_71std_climate_anomaly,breaks=brks,col=cols,xaxt="n",yaxt="n",bty="n")
text(-175,-20,"2070 - RCP 8.5\n (Standardised)",pos=4)
plot.new()
legend(0.1,1,c("< -5","-5 : -2","-2 : -1","-1 : -0.5","-0.5 : -0.2","-0.2 : -0.1","-0.1 : 0",
             "0 : 0.1","0.1 : 0.2","0.2 : 0.5","0.5 : 1","1 : 2","2 : 5","> 5"),
       ncol=5,fill=cols,bty="n")
# image(tmp2016_18std_climate_anomaly,breaks=brks,col=cols,xaxt="n",yaxt="n",bty="n")

invisible(dev.off())
