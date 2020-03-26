###Calculate proportion NH surrounding each site in PREDICTS at multiple spatial scales, 1,3,5 and 10km.


library(raster)
library(sp)
library(predictsFunctions)

wgs84 <- "+proj=longlat +datum=WGS84 +ellps=WGS84 +towgs84=0,0,0"
  
predicts_sites <- create_sites(predicts)
predicts_sites2 <- predicts_sites[!is.na(predicts_sites$Latitude), ]
predicts_sp <- SpatialPointsDataFrame(coords = cbind(predicts_sites2$Longitude, predicts_sites2$Latitude), data = predicts_sites2, proj4string = wgs84)


PV <- raster("Data/Global_30s_resolution_land_use_for_2005/PRI_2005/PRI_1km_2005_0ice.bil", crs=wgs84)
SV <- raster("Data/Global_30s_resolution_land_use_for_2005/SEC_2005/SEC_1km_2005_0ice.bil", crs=wgs84)


nrow(predicts_sp)

##Check number of rows needed
for(i in 1:22661) {
  buff1 <- raster::buffer(predicts_sp[i, ], width=1000)
  PV_crop1 <- raster::crop(PV, buff1)
  PV_mean1 <- cellStats(PV_crop1, stat="mean")
  SV_crop1 <- raster::crop(SV, buff1)
  SV_mean1 <- cellStats(SV_crop1, stat="mean")
  predicts_sp$PV_1000[i] <- PV_mean1
  predicts_sp$SV_1000[i] <-  SV_mean1
  
  
  buff3 <- raster::buffer(predicts_sp[i, ], width=3000)
  PV_crop3 <- raster::crop(PV, buff3)
  PV_mean3 <- cellStats(PV_crop3, stat="mean")
  SV_crop3 <- raster::crop(SV, buff3)
  SV_mean3 <- cellStats(SV_crop3, stat="mean")
  predicts_sp$PV_3000[i] <- PV_mean3
  predicts_sp$SV_3000[i] <-  SV_mean3
  
  
  buff5 <- raster::buffer(predicts_sp[i, ], width=5000)
  PV_crop5 <- raster::crop(PV, buff5)
  PV_mean5 <- cellStats(PV_crop5, stat="mean")
  SV_crop5 <- raster::crop(SV, buff5)
  SV_mean5 <- cellStats(SV_crop5, stat="mean")
  predicts_sp$PV_5000[i] <- PV_mean5
  predicts_sp$SV_5000[i] <-  SV_mean5
  
  
  buff10 <- raster::buffer(predicts_sp[i, ], width=10000)
  PV_crop10 <- raster::crop(PV, buff10)
  PV_mean10 <- cellStats(PV_crop10, stat="mean")
  SV_crop10 <- raster::crop(SV, buff1)
  SV_mean10 <- cellStats(SV_crop10, stat="mean")
  predicts_sp$PV_10000[i] <- PV_mean10
  predicts_sp$SV_10000[i] <-  SV_mean10
  print(i)
  
}
saveRDS(predicts_sp, "Outputs/predicts_NH_multiple_scales.rds")




