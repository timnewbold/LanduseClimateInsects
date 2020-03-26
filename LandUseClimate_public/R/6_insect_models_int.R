library(StatisticalModels)
library(sjPlot)

##insect models interactions

climate_int <- GLMER(insect_sites, 
                     responseVar = "log_scaled_abundance", fitFamily = "gaussian",
                     fixedStruct = "Land*std_climate_anomaly.s*NH_5000.s",
                     randomStruct = "(1|SS) + (1|SSB)", saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_5000"))
summary(climate_int$model)
vif.mer(climate_int$model)
par(mfrow=c(3,1))
for(UI in levels(climate_int$data$Land)[2:4]) {
  
  #Mean and SD NAT HABITAT
  mean.ag2NH <- mean(na.omit(insect_sites$NH_5000))
  sd.ag2NH <- sd(na.omit(insect_sites$NH_5000))
  
  ##MEAN AND SD for centerd anomaly
  mean.anom <- mean(na.omit(insect_sites$std_climate_anomaly))
  sd.anom <- sd(na.omit(insect_sites$std_climate_anomaly))
  
  std_climate_anomaly <- seq(from=0, to=2, by=0.1)
  std_climate_anomaly_t <- (std_climate_anomaly -mean.anom) / sd.anom
  
  
  
  newdat10 <- data.frame(std_climate_anomaly.s=std_climate_anomaly_t, NH_5000.s=rep((0.1-mean.ag2NH) /sd.ag2NH , 1),
                         log_scaled_abundance = rep(0, 1), 
                         Land=factor(rep((UI), 1), levels= levels(climate_int$data$Land)))
  
  newdat10 <- cbind(newdat10, exp(PredictGLMER(climate_int$model, data = newdat10, se.fit = T, seMultiplier = 1)))
  
  newdat50 <- data.frame(std_climate_anomaly.s=std_climate_anomaly_t, NH_5000.s=rep((0.5-mean.ag2NH) /sd.ag2NH , 1),
                         log_scaled_abundance = rep(0, 1), 
                         Land=factor(rep((UI), 1), levels= levels(climate_int$data$Land)))
  
  newdat50 <- cbind(newdat50, exp(PredictGLMER(climate_int$model, data = newdat50, se.fit = T, seMultiplier = 1)))
  
  
  newdat70 <- data.frame(std_climate_anomaly.s=std_climate_anomaly_t, NH_5000.s=rep((0.7-mean.ag2NH) /sd.ag2NH , 1),
                         log_scaled_abundance = rep(0, 1), 
                         Land=factor(rep((UI), 1), levels= levels(climate_int$data$Land)))
  
  newdat70 <- cbind(newdat70, exp(PredictGLMER(climate_int$model, data = newdat70, se.fit = T, seMultiplier = 1)))
  
  
  plot(newdat10$std_climate_anomaly, newdat10$y, type="l", 
       col="red", main=paste(UI) , ylim=c(0,0.5),
       xlab="Standardised climate anomaly", ylab="Scaled Abundance", xaxt="n")
  axis(1, at=std_climate_anomaly_t[c(1, 6, 11, 16, 21)], labels=seq(from=0, to=2, by=0.5))
  lines(newdat50$std_climate_anomaly, newdat50$y, type="l", col="tan2")
  lines(newdat70$std_climate_anomaly, newdat70$y, type="l", col="green")
  
  
  rug(jitter(climate_int$data$std_climate_anomaly.s[climate_int$data$Land==UI], amount = 0.01),
      side = 3,  ticksize = 0.04, lwd = 0.1)
  #rug(climate_int$data$std_climate_anomaly.s[climate_int$data$Land==UI],
  #   side=3, ticksize = 0.04, lwd = 0.01)
  #text(std_climate_anomaly_t[11], 0.5, paste(UI, "NH_5000"))
  polygon(c(newdat10$std_climate_anomaly, rev(newdat10$std_climate_anomaly)), c(newdat10$yplus, rev(newdat10$yminus)), 
          border = "tomato3", col=adjustcolor("tomato3", alpha.f = 0.4))
  polygon(c(newdat50$std_climate_anomaly, rev(newdat50$std_climate_anomaly)), c(newdat50$yplus, rev(newdat50$yminus)), 
          border = "tan2", col=adjustcolor("tan2", alpha.f = 0.4))
  polygon(c(newdat70$std_climate_anomaly, rev(newdat70$std_climate_anomaly)), c(newdat70$yplus, rev(newdat70$yminus)), 
          border = "forestgreen", col=adjustcolor("forestgreen", alpha.f = 0.4))
  legend("bottomleft", col =c("forestgreen", "tan2", "tomato3"), 
         pch=19, legend = c("70%", "50%", "10%"), title="Percentage NH", cex=1.5)
  
}
tab_model(climate_int$model, show.stat = T, collapse.ci = T)


tmax_int <- GLMER(insect_sites, 
                  responseVar = "log_scaled_abundance", fitFamily = "gaussian",
                  fixedStruct = "Land*std_tmax_anomaly.s*NH_5000.s",
                  randomStruct = "(1|SS) + (1|SSB)", saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_5000", "std_tmax_anomaly"))
summary(tmax_int$model)
vif.mer(tmax_int$model)

sjPlot::tab_model(climate_int$model, tmax_int$model, show.stat = T, collapse.ci = T)

par(mfrow=c(3,1))
for(UI in levels(tmax_int$data$Land)[2:4]) {
  
  #Mean and SD NAT HABITAT
  mean.ag2NH <- mean(na.omit(insect_sites$NH_5000))
  sd.ag2NH <- sd(na.omit(insect_sites$NH_5000))
  
  ##MEAN AND SD for centerd anomaly
  mean.anom <- mean(na.omit(insect_sites$std_tmax_anomaly))
  sd.anom <- sd(na.omit(insect_sites$std_tmax_anomaly))
  
  std_tmax_anomaly <- seq(from=0, to=2, by=0.1)
  std_tmax_anomaly_t <- (std_tmax_anomaly -mean.anom) / sd.anom
  
  
  
  newdat10 <- data.frame(std_tmax_anomaly.s=std_tmax_anomaly_t, NH_5000.s=rep((0.1-mean.ag2NH) /sd.ag2NH , 1),
                         log_scaled_abundance = rep(0, 1), 
                         Land=factor(rep((UI), 1), levels= levels(tmax_int$data$Land)))
  
  newdat10 <- cbind(newdat10, exp(PredictGLMER(tmax_int$model, data = newdat10, se.fit = T, seMultiplier = 1)))
  
  newdat50 <- data.frame(std_tmax_anomaly.s=std_tmax_anomaly_t, NH_5000.s=rep((0.5-mean.ag2NH) /sd.ag2NH , 1),
                         log_scaled_abundance = rep(0, 1), 
                         Land=factor(rep((UI), 1), levels= levels(tmax_int$data$Land)))
  
  newdat50 <- cbind(newdat50, exp(PredictGLMER(tmax_int$model, data = newdat50, se.fit = T, seMultiplier = 1)))
  
  
  newdat70 <- data.frame(std_tmax_anomaly.s=std_tmax_anomaly_t, NH_5000.s=rep((0.7-mean.ag2NH) /sd.ag2NH , 1),
                         log_scaled_abundance = rep(0, 1), 
                         Land=factor(rep((UI), 1), levels= levels(tmax_int$data$Land)))
  
  newdat70 <- cbind(newdat70, exp(PredictGLMER(tmax_int$model, data = newdat70, se.fit = T, seMultiplier = 1)))
  
  
  plot(newdat10$std_tmax_anomaly, newdat10$y, type="l", 
       col="red", main=paste(UI) , ylim=c(0,0.5),
       xlab="Standardised tmax anomaly", ylab="Scaled Abundance", xaxt="n")
  axis(1, at=std_tmax_anomaly_t[c(1, 6, 11, 16, 21)], labels=seq(from=0, to=2, by=0.5))
  lines(newdat50$std_tmax_anomaly, newdat50$y, type="l", col="tan2")
  lines(newdat70$std_tmax_anomaly, newdat70$y, type="l", col="green")
  
  
  rug(jitter(tmax_int$data$std_tmax_anomaly.s[tmax_int$data$Land==UI], amount = 0.01),
      side = 3,  ticksize = 0.04, lwd = 0.1)
  #rug(tmax_int$data$std_tmax_anomaly.s[tmax_int$data$Land==UI],
  #   side=3, ticksize = 0.04, lwd = 0.01)
  #text(std_tmax_anomaly_t[11], 0.5, paste(UI, "NH_5000"))
  polygon(c(newdat10$std_tmax_anomaly, rev(newdat10$std_tmax_anomaly)), c(newdat10$yplus, rev(newdat10$yminus)), 
          border = "tomato3", col=adjustcolor("tomato3", alpha.f = 0.4))
  polygon(c(newdat50$std_tmax_anomaly, rev(newdat50$std_tmax_anomaly)), c(newdat50$yplus, rev(newdat50$yminus)), 
          border = "tan2", col=adjustcolor("tan2", alpha.f = 0.4))
  polygon(c(newdat70$std_tmax_anomaly, rev(newdat70$std_tmax_anomaly)), c(newdat70$yplus, rev(newdat70$yminus)), 
          border = "forestgreen", col=adjustcolor("forestgreen", alpha.f = 0.4))
  legend("bottomleft", col =c("forestgreen", "tan2", "tomato3"), 
         pch=19, legend = c("70%", "50%", "10%"), title="Percentage NH", cex=1.5)
  
}


climate_rich_add <- GLMER(insect_sites, 
                          responseVar = "Species_richness", fitFamily = "poisson",
                          fixedStruct = "Land*std_climate_anomaly.s+Land*NH_5000.s",
                          randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_5000"))
summary(climate_rich_add$model)

climate_rich_int <- GLMER(insect_sites, 
                          responseVar = "Species_richness", fitFamily = "poisson",
                          fixedStruct = "Land*std_climate_anomaly.s*NH_5000.s",
                          randomStruct = "(1|SS) + (1|SSB) + (1|SSBS)", saveVars = c("Species_richness", "Total_abundance", "SSBS", "NH_5000"))




summary(climate_rich_int$model)
vif.mer(climate_rich_int$model)
