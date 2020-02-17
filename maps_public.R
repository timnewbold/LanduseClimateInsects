###Map of standardized anomaly

##Packages
library(raster)
library(sp)
library(dismo)
library(rgdal)
library(RColorBrewer)

tmp <- stack("Your CRU v4.03 tmp file name here")

tmp1901_1905 <- tmp[[names(tmp)[1:60]]]
names(tmp1901_1905)

tmp1901_1905mean <- calc(tmp1901_1905, base::mean)
tmp1901_1905sd <- calc(tmp1901_1905, function(x) stats::sd(x, na.rm = T))
tmp1901_1905sd

###2005 is mean year for insect data
###2004 to 2006 maps
names(tmp)[1237:1272] ##2004_6
tmp2004_6 <- tmp[[names(tmp)[1237:1272]]]
tmp2004_6sd <-  calc(tmp[[names(tmp)[1237:1272]]], stats::sd)
tmp2004_6mean <-  calc(tmp[[names(tmp)[1237:1272]]], base::mean)
tmp2004_6std_climate_anomaly <- (calc(tmp2004_6, base::mean)-tmp1901_1905mean)  / tmp1901_1905sd
plot(tmp2004_6std_climate_anomaly)
breaks <- c(-2,0,0.25, 0.5,0.75, 1,1.5, 2,2.5, 3,4,5, 10)
pallete <- colorRampPalette(c("lightblue","red", "black"))
plot(tmp2004_6std_climate_anomaly, breaks=breaks, col=pallete(12), main="2004_6 historic")



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

