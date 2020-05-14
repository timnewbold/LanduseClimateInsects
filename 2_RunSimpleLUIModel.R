library(StatisticalModels)

inDir <- "1_PreparePREDICTSData/"

outDir <- "2_RunSimpleLUIModel/"

sites <- readRDS(file = paste0(inDir,"PREDICTSSiteData.rds"))

model_data <- na.omit(sites[,c('Species_richness','LandUse','Use_intensity','UI2','SS','SSB','SSBS')])

print(table(model_data$UI2))

sm0 <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "1",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)
sm1 <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "LandUse",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)
sm2 <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)
sm3 <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "UI2",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)
sm4 <- GLMER(modelData = model_data,responseVar = "Species_richness",fitFamily = "poisson",
             fixedStruct = "LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)+(1|SSBS)",REML = FALSE)

print(AIC(sm0$model,sm1$model,sm2$model,sm3$model,sm4$model))

model_data <- na.omit(sites[,c('LogAbund','LandUse','Use_intensity','UI2','SS','SSB','SSBS')])

am0 <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "1",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)
am1 <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "LandUse",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)
am2 <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "LandUse+Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)
am3 <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "UI2",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)
am4 <- GLMER(modelData = model_data,responseVar = "LogAbund",fitFamily = "gaussian",
             fixedStruct = "LandUse*Use_intensity",randomStruct = "(1|SS)+(1|SSB)",REML = FALSE)

print(AIC(am0$model,am1$model,am2$model,am3$model,am4$model))

nd <- data.frame(UI2=factor(c("Primary vegetation","Secondary vegetation",
                              "Agriculture_Low","Agriculture_High","Urban")),
                 Species_richness=0,
                 LogAbund=0)

s.preds <- PredictGLMERRandIter(model = sm3$model,data = nd)

s.preds <- exp(s.preds)

s.preds <- sweep(x = s.preds,MARGIN = 2,STATS = s.preds[1,],FUN = '/')

s.preds.median <- ((apply(X = s.preds,MARGIN = 1,FUN = median))*100)-100
s.preds.upper <- ((apply(X = s.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
s.preds.lower <- ((apply(X = s.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

a.preds <- PredictGLMERRandIter(model = am3$model,data = nd)

a.preds <- exp(a.preds)-1

a.preds <- sweep(x = a.preds,MARGIN = 2,STATS = a.preds[1,],FUN = '/')

a.preds.median <- ((apply(X = a.preds,MARGIN = 1,FUN = median))*100)-100
a.preds.upper <- ((apply(X = a.preds,MARGIN = 1,FUN = quantile,probs = 0.975))*100)-100
a.preds.lower <- ((apply(X = a.preds,MARGIN = 1,FUN = quantile,probs = 0.025))*100)-100

pdf(file = paste0(outDir,"LUI_Plot.pdf"),width = 8.5/2.54,height = 12/2.54)

par(mfrow=c(2,1))
par(cex=1)
par(cex.lab=1)
par(cex.axis=1)
par(cex.main=1)
par(ps=10)
par(las=1)
par(mgp=c(1.6,0.2,0))
par(mar=c(2.6,2.6,0.5,0.5))
par(tck=-0.01)


errbar.cols <- c("#009E73","#0072B2","#E69F00","#D55E00","#CC79A7")

errbar(x = 1:5,y = s.preds.median,yplus = s.preds.upper,yminus = s.preds.lower,
       col=errbar.cols,errbar.col = errbar.cols,ylim=c(min(s.preds.lower),max(s.preds.upper)),xaxt="n",
       ylab="Species richness (%)",xlab="Land use",bty="l")

axis(side = 1,at = 1:5,labels = c("PV","SV","AG.Low","AG.Hi","URB"))

abline(h=0,col="#00000077",lty=2)

errbar(x = 1:5,y = a.preds.median,yplus = a.preds.upper,yminus = a.preds.lower,
       col=errbar.cols,errbar.col = errbar.cols,ylim=c(min(a.preds.lower),max(a.preds.upper)),xaxt="n",
       ylab="Total abundance (%)",xlab="Land use",bty="l")

axis(side = 1,at = 1:5,labels = c("PV","SV","AG.Low","AG.Hi","URB"))

abline(h=0,col="#00000077",lty=2)

invisible(dev.off())
