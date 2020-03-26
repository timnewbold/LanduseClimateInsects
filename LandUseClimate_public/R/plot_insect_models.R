
cols <- c("#003300", "#00cc99", "#ffcc00", "#800000")

##PLots of insect models

##Axist labels
std_climate_anomaly_t <- c(-0.5, 0, 0.5 , 1 ,1.5, 2)
std_climate_anomaly_t <- (std_climate_anomaly_t - mean(insect_sites$std_climate_anomaly, na.rm=T)) / sd(insect_sites$std_climate_anomaly, na.rm=T)

std_tmax_anomaly_t <- c(-0.5, 0, 0.5 , 1 ,1.5, 2, 2.5)
std_tmax_anomaly_t <- (std_tmax_anomaly_t - mean(insect_sites$std_tmax_anomaly, na.rm=T)) / sd(insect_sites$std_tmax_anomaly, na.rm=T)

std_climate_anomaly_t

par(mfrow=c(2,2))
par(oma=c(1,1,1,1))
PlotGLMERContinuous(model = climate$model,data = climate$data,
                    effects = "std_climate_anomaly.s", byFactor = "Land",
                    xlab = "Standardised Climate Anomaly",ylab = "Scaled Abundance",
                    logLink = "e", ylim = c(0,0.5),
                    line.cols =cols, plotRug = F,seMultiplier = 1, params=list(xaxt="n"))
par(xaxt="s")
legend(-2, 0.2, legend = levels(climate$data$Land), col = cols, pch=19, cex=1)
axis(side = 1, at=std_climate_anomaly_t, labels=c(-0.5, 0,0.5,1,1.5, 2))
#mtext(side=3, text="Standardised climate anomaly", padj=-1.3, cex=1)

rug(climate$data$std_climate_anomaly.s[climate$data$Land=="PV"], col = cols[1])
rug(climate$data$std_climate_anomaly.s[climate$data$Land=="SV"], col = cols[2])
rug(climate$data$std_climate_anomaly.s[climate$data$Land=="Low Agriculture"], col = cols[3])
rug(climate$data$std_climate_anomaly.s[climate$data$Land=="High Agriculture"], col = cols[4])




PlotGLMERContinuous(model = tmax$model,data = tmax$data,
                    effects = "std_tmax_anomaly.s", byFactor = "Land",
                    xlab = "Standardised Tmax Anomaly",ylab = "Scaled Abundance",
                    logLink = "e", ylim = c(0,0.5),
                    line.cols =cols, plotRug = FALSE,seMultiplier = 1, params=list(xaxt="n"))

par(xaxt="s")
legend(-2, 0.2, legend = levels(tmax$data$Land), col = cols, pch=19, cex=1)
axis(side = 1, at=std_tmax_anomaly_t, labels=c(-0.5, 0,0.5,1,1.5, 2, 2.5))
#mtext(side=3, text="Tmax", padj=-1.3, cex=1)

rug(tmax$data$std_tmax_anomaly.s[tmax$data$Land=="PV"], col = cols[1])
rug(tmax$data$std_tmax_anomaly.s[tmax$data$Land=="SV"], col = cols[2])
rug(tmax$data$std_tmax_anomaly.s[tmax$data$Land=="Low Agriculture"], col = cols[3])
rug(tmax$data$std_tmax_anomaly.s[tmax$data$Land=="High Agriculture"], col = cols[4])




PlotGLMERContinuous(model = climate_rich$model,data = climate_rich$data,
                    effects = "std_climate_anomaly.s", byFactor = "Land",
                    xlab = "Standardised Climate Anomaly",ylab = "Species Richness",
                    logLink = "e", ylim = c(0,22),
                    line.cols =cols, plotRug = FALSE,seMultiplier = 1, param=list(xaxt="n"))

par(xaxt="s")
legend(-2, 7, legend = levels(climate$data$Land), col = cols, pch=19, cex=1)
axis(side = 1, at=std_climate_anomaly_t, labels=c(-0.5, 0,0.5,1,1.5, 2))
#mtext(side=3, text="Standardised climate anomaly", padj=-1.3, cex=1)

rug(climate_rich$data$std_climate_anomaly.s[climate_rich$data$Land=="PV"], col = cols[1])
rug(climate_rich$data$std_climate_anomaly.s[climate_rich$data$Land=="SV"], col = cols[2])
rug(climate_rich$data$std_climate_anomaly.s[climate_rich$data$Land=="Low Agriculture"], col = cols[3])
rug(climate_rich$data$std_climate_anomaly.s[climate_richx$data$Land=="High Agriculture"], col = cols[4])



PlotGLMERContinuous(model = tmax_rich$model,data = tmax_rich$data,
                    effects = "std_tmax_anomaly.s", byFactor = "Land",
                    xlab = "Standardised Tmax Anomaly",ylab = "Species richness",
                    logLink = "e", ylim = c(0,20),
                    line.cols =cols, plotRug = FALSE,seMultiplier = 1, param=list(xaxt="n"))


par(xaxt="s")
legend(-2, 6, legend = levels(tmax$data$Land), col = cols, pch=19, cex=1)
axis(side = 1, at=std_tmax_anomaly_t, labels=c(-0.5, 0,0.5,1,1.5, 2, 2.5))
#mtext(side=3, text="Tmax", padj=-1.3, cex=1)

rug(tmax_rich$data$std_tmax_anomaly.s[tmax_rich$data$Land=="PV"], col = cols[1])
rug(tmax_rich$data$std_tmax_anomaly.s[tmax_rich$data$Land=="SV"], col = cols[2])
rug(tmax_rich$data$std_tmax_anomaly.s[tmax_rich$data$Land=="Low Agriculture"], col = cols[3])
rug(tmax_rich$data$std_tmax_anomaly.s[tmax_rich$data$Land=="High Agriculture"], col = cols[4])





