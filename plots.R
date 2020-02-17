##Log scaled abundance vs climate anomaly additive model

PlotGLMERContinuous(model = mod_add$model,data = mod_add$data,
                    effects = "std_climate_anomaly", otherFactors = list(LandUse=c("Primary vegetation")),
                    otherContEffects="percNH.s",
                    xlab = "Standardised Climate anomaly",ylab = "Scaled Abundance",
                    logLink = "e", main ="Additive model all insects standardized for percNH", ylim = c(0,0.6),
                    line.cols ="#003300", plotRug = TRUE,seMultiplier = 1)

PlotGLMERContinuous(model = mod_add$model,data = mod_add$data,
                    effects = "std_climate_anomaly", otherFactors = list(LandUse=c("Cropland")),
                    otherContEffects="percNH.s",
                    logLink = "e",line.cols ="#ff9999", plotRug = TRUE, add=TRUE, seMultiplier = 1)

PlotGLMERContinuous(model = mod_add$model,data = mod_add$data,
                    effects = "std_climate_anomaly", otherFactors = list(LandUse=c("Pasture")), logLink = "e",
                    otherContEffects="percNH.s",
                    line.cols ="#999966", plotRug = TRUE, add=TRUE, seMultiplier = 1)

PlotGLMERContinuous(model = mod_add$model,data = mod_add$data,
                    effects = "std_climate_anomaly", otherFactors = list(LandUse=c("Plantation forest")), logLink = "e",
                    line.cols ="#994d00", plotRug = TRUE, add=TRUE, seMultiplier = 1, otherContEffects="percNH.s")

legend(-0.5, 0.15, legend = c("Primary vegetation", "Cropland", "Pasture", "Plantation forest"),
       col=c("#003300", "#ff9999", "#999966", "#994d00"), pch = 19)


##interactive model plot


mod_int <- GLMER(insect_sites, responseVar = "log_scaled_abundance", fitFamily = "gaussian",
                 fixedStruct = "LandUse*std_climate_anomaly*percNH",
                 randomStruct = "(1|SS) + (1|SSB)")

newdat10 <- data.frame(LandUse=factor(rep("Cropland", 31), levels = levels(mod_int$data$LandUse)),
                       std_climate_anomaly=seq(from=0, to=3, by=0.1), percNH=rep(100 , 31),
                       log_scaled_abundance = rep(0, 31))

NH10 <- cbind(newdat10, exp(PredictGLMER(mod_int$model, data = newdat10, se.fit = T, seMultiplier = 1)))

newdat50 <- data.frame(LandUse=factor(rep("Cropland", 31), levels = levels(mod_int$data$LandUse)),
                       std_climate_anomaly=seq(from=0, to=3, by=0.1), percNH=rep(500, 31),
                       log_scaled_abundance = rep(0, 31))

NH50 <- cbind(newdat50, exp(PredictGLMER(mod_int$model, data = newdat50, se.fit = T, seMultiplier = 1)))


newdat70 <- data.frame(LandUse=factor(rep("Cropland", 31), levels = levels(mod_int$data$LandUse)),
                       std_climate_anomaly=seq(from=0, to=3, by=0.1), percNH=rep(700 , 31),
                       log_scaled_abundance = rep(0, 31))

NH70 <- cbind(newdat70, exp(PredictGLMER(mod_int$model, data = newdat70, se.fit = T, seMultiplier = 1)))


plot(NH70$std_climate_anomaly, NH70$y, type="l", col="forestgreen", ylim=c(0, 0.35),
     main=" percNH surrounding cropland reduces climate change effects on insect abundance",
     xlab="Standardised Climate anomaly", ylab="Scaled Abundance")
lines(NH50$std_climate_anomaly, NH50$y, col="darkorange")
lines(NH10$std_climate_anomaly, NH10$y, col="darkred")
polygon(c(NH10$std_climate_anomaly, rev(NH10$std_climate_anomaly)), c(NH10$yplus, rev(NH10$yminus)), 
        border = "firebrick1", col=adjustcolor("firebrick1", alpha.f = 0.4))
polygon(c(NH10$std_climate_anomaly, rev(NH10$std_climate_anomaly)), c(NH50$yplus, rev(NH50$yminus)), 
        border = "darkorange", col=adjustcolor("darkorange", alpha.f = 0.4))
polygon(c(NH70$std_climate_anomaly, rev(NH70$std_climate_anomaly)), c(NH70$yplus, rev(NH70$yminus)), 
        border = "forestgreen", col=adjustcolor("forestgreen", alpha.f = 0.4))
legend(0.1, 0.1, col =c("forestgreen", "darkorange", "firebrick1"), 
                        pch=19, legend = c("70%", "50%", "10%"), title="Percentage NH")

rug(mod_int$data$std_climate_anomaly[mod_int$data$std_climate_anomaly>=0])       

