library(dplyr)

nat_habitat <- readRDS("Outputs/predicts_NH_multiple_scales.rds")
nat_habitat
nat_habitat <- dplyr::select(as.data.frame(nat_habitat), SSBS, PV_1000:SV_10000)
nat_habitat
predicts_climate_info <- readRDS("Outputs/predicts_climate_info.rds")
predicts_climate_info <- as.data.frame(predicts_climate_info)
head(predicts_climate_info)
predicts_climate_info <- dplyr::select(predicts_climate_info, Source_ID:historic_sd_tmax)
predicts_sites_NH_combined <- dplyr::left_join(predicts_climate_info, nat_habitat, by="SSBS")

predicts_sites_NH_combined
saveRDS(predicts_sites_NH_combined, "Outputs/predicts_sites_info.rds")


