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


brks <- c(-50,-5,-2,-1,-0.5,-0.2,-0.1,0,0.1,0.2,0.5,1,2,5,50)
cols <- c(rev(brewer.pal(n = 9,name = "Blues"))[3:9],(brewer.pal(n = 9,name = "Reds"))[3:9])

pdf(file = paste0(outDir,"ClimateIndexMaps.pdf"),width = 17.5,height = 26.25)

par(mfrow=c(3,1))

par(mar=c(0,0,0,0))

image(tmp1970std_climate_anomaly,breaks=brks,col=cols,xaxt="n",yaxt="n",bty="n")
image(tmp2004_6std_climate_anomaly,breaks=brks,col=cols,xaxt="n",yaxt="n",bty="n")
image(tmp2016_18std_climate_anomaly,breaks=brks,col=cols,xaxt="n",yaxt="n",bty="n")

invisible(dev.off())
