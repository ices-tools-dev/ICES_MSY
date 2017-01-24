####################################
# Length-based Spawner per recruit #
####################################
rm(list = ls())
source("~/ICES_MSY/R/LB-SPR_functions.R")

# # Test LBSPR Function 
# FM <- 0.9
# SL50 <- 40
# SL95 <- 50
# # 
# testFun <- LBSPRFunc(MK=1.5, Linf=100, CVLinf=0.1, L50=55, L95=60, Beta=3, FM=FM, SL50=SL50, SL95=SL95,
# 	By=2, MaxLMult=1.25, P=0.01, Nage=100)
# testFun$SPR
# plot(testFun$LenMids, testFun$PropLen, type="l")

#########################
# Fit to Empirical Data #
#########################
# Import Length Data
RawDat <- read.csv("~/ICES_MSY/data/nep-2829/LFD/nep2829_LFD_males.csv")
colnames(RawDat)[1] <- "CL"
colnames(RawDat)[-1] <- gsub("X", "", colnames(RawDat)[-1])

# Biological Parameters 
MK <- 1.5
Linf <- 70  
CVLinf <- 0.1  
L50 <- 28.4
L95 <- 28.5
Beta <- 3.2229 

# Other pararameters
Nage <- 100 
P <- 0.01

# To store results
rslts <- data.frame(year=2000:2014,
                    FM=rep(0,length(2000:2014)),
                    FMl=rep(0,length(2000:2014)),
                    FMh=rep(0,length(2000:2014)),
                    SL50=rep(0,length(2000:2014)),
                    SL50l=rep(0,length(2000:2014)),
                    SL50h=rep(0,length(2000:2014)),
                    SL95=rep(0,length(2000:2014)),
                    SL95l=rep(0,length(2000:2014)),
                    SL95h=rep(0,length(2000:2014)),
                    SPR=rep(0,length(2000:2014)),
                    SPRl=rep(0,length(2000:2014)),
                    SPRh=rep(0,length(2000:2014)))

col <- 2

# To plot fits
png(filename="~/ICES_MSY/output/nep2829_fits.png", 
    width = 8,
    height = 6,
    units="in",
    res= 100)

par(mfrow=c(5,3), mar=c(2,2,3,1))
year = 2005
# Run through each year
for (year in 2000:2014){

  col <- grep(year, colnames(RawDat))
  # Binning parameters 
  By <- 1
  MaxLMult <- 1.25
  LenMids <- seq(from=By *0.5, to= Linf * MaxLMult, by=By)
  LenBins <- seq(from=0, by=By, length.out=length(LenMids)+1)
  #LenDat <- as.vector(table(cut(RawDat, LenBins)))
  LenDat <- rep(0,length(LenMids))
  LenDat[is.element(LenMids, RawDat$CL)] <- RawDat[,col]
  N <- sum(LenDat) 
  
  # Starting parameters
  SL50Start <- LenMids[which.max(LenDat)]/Linf 
  DeltaStart <- LenMids[which.max(LenDat)+1]/Linf / SL50Start
  Starts <- log(c(1, SL50Start, DeltaStart))
  
  # Run optimisation
  Opt <- nlminb(Starts, 
                OptRoutine, 
                LenDat = LenDat, 
                MK = MK,
                Linf = Linf, 
                CVLinf = CVLinf, 
                L50 = L50,
                L95 = L95, 
                Beta = Beta, 
                By = By, 
                MaxLMult = MaxLMult, 
                P = P, 
                Nage = Nage)
  
  # Estimated parameters
  estFM <- exp(Opt$par[1])	
  estSL50 <- exp(Opt$par[2]) * Linf
  estSL95 <- estSL50 + (exp(Opt$par[3]) * estSL50	)
  
  # Error estimable parameters
  Err <- sqrt(diag(solve(hessian(OptRoutine,
                                 Opt$par,
                                 LenDat = LenDat, 
                                 MK = MK, 
                                 Linf = Linf, 
                                 CVLinf = CVLinf, 
                                 L50 = L50, 
                                 L95 = L95, 
                                 Beta = Beta, 
                                 By = By, 
                                 MaxLMult = MaxLMult, 
                                 P = P, 
                                 Nage = Nage))))
  
  estFMl <- exp(Opt$par[1]-1.96*Err[1])
  estFMh <- exp(Opt$par[1]+1.96*Err[1])
  estSL50l <- exp(Opt$par[2]-1.96*Err[2]) * Linf
  estSL50h <- exp(Opt$par[2]+1.96*Err[2]) * Linf
  estSL95l <- estSL50l + (exp(Opt$par[3]-1.96*Err[3]) * estSL50l)
  estSL95h <- estSL50h + (exp(Opt$par[3]+1.96*Err[3]) * estSL50h)
  
  # Variance SPR
  SPRvar <- delta.meth.1(Opt$par, hessian(OptRoutine,Opt$par, LenDat=LenDat, MK=MK, Linf=Linf, CVLinf=CVLinf, L50=L50, L95=L95,Beta=Beta, By=By, MaxLMult=MaxLMult, P=P, Nage=Nage), LBSPRFuncwrap, MK=MK, Linf=Linf, CVLinf=CVLinf, L50=L50, L95=L95, 
               Beta=Beta, By=By, MaxLMult=MaxLMult, P=P, Nage=Nage)
  
  # Run model with estimated parameters 
  runMod <- LBSPRFunc(MK=MK, Linf=Linf, CVLinf=CVLinf, L50=L50, L95=L95, 
                      Beta=Beta, FM=estFM, SL50=estSL50, SL95=estSL95, By=By, MaxLMult=MaxLMult, P=P, Nage=Nage)
  
  tt <- barplot(LenDat, names.arg=LenMids, ylim=c(0, max(LenDat,runMod$PropLen * N)), main=year)  # plot binned data 
  lines(tt, runMod$PropLen * N, lwd=4)
  
  # Store results
  rslts$FM[rslts$year==year] <- estFM 
  rslts$FMl[rslts$year==year] <- estFMl
  rslts$FMh[rslts$year==year] <- estFMh
  rslts$SL50[rslts$year==year] <- estSL50
  rslts$SL50l[rslts$year==year] <- estSL50l
  rslts$SL50h[rslts$year==year] <- estSL50h
  rslts$SL95[rslts$year==year] <- estSL95
  rslts$SL95l[rslts$year==year] <- estSL95l
  rslts$SL95h[rslts$year==year] <- estSL95h
  rslts$SPR[rslts$year==year] <- runMod$SPR
  rslts$SPRl[rslts$year==year] <- runMod$SPR - 1.96*sqrt(SPRvar)
  rslts$SPRh[rslts$year==year] <- runMod$SPR + 1.96*sqrt(SPRvar)
  
  # col <- col + 1
}

dev.off()

# Estimates 
# c(FM=estFM, SL50=estSL50, SL95=estSL95, SPR=runMod$SPR)

## PLOTS ###
# Plot SPR
plot(rslts$year, rslts$SPR, type='p', lwd=2, xlab="", ylab="SPR", ylim=c(0,1), main="nep2829_males")
arrows(rslts$year, rslts$SPRl, rslts$year, rslts$SPRh, length=0.1, angle=90, code=3, lwd=2)
abline(h=0.4, lwd=2, col='red')
abline(h=0.35, lwd=2, lty=2, col='red')
abline(h=0.3, lwd=2, col='red')

# Plot F/M
plot(rslts$year, rslts$FM, type='p', lwd=2, xlab="", ylab="F/M", main="nep2829_males", ylim=c(min(rslts$FMl),max(rslts$FMh)))
arrows(rslts$year, rslts$FMl, rslts$year, rslts$FMh, length=0.1, angle=90, code=3, lwd=2)
abline(h=1, lwd=2, col='red')

# Plot selectivity
plot(rslts$year, rslts$SL50, type='p', lwd=2, xlab="", ylab="Selectivity", main="nep2829_males", ylim=c(min(rslts$SL50l),max(rslts$SL95h)))
arrows(rslts$year, rslts$SL50l, rslts$year, rslts$SL50h, length=0.1, angle=90, code=3, lwd=2)
points(rslts$year, rslts$SL95, lwd=2, col='blue')
arrows(rslts$year, rslts$SL95l, rslts$year, rslts$SL95h, length=0.1, angle=90, code=3, lwd=2, col='blue')
abline(h=mean(rslts$SL50), lwd=2)
abline(h=mean(rslts$SL95), lwd=2, col='blue')

# Selectivity at length by year
SelMat <- RawDat
SelMat[,2:length(SelMat)] <- 0
for(col in 1:(ncol(SelMat)-1)){
  SelMat[, col+1] <- SLen(SelMat$CL, rslts$SL50[col], rslts$SL95[col])
}

# Mean length and mean selected length
rslts$Lbar <- apply(RawDat[,2:length(RawDat)], 2, function(x) sum(x*RawDat$CL)/sum(x))
rslts$Lsel <- apply(RawDat[,2:length(RawDat)]*SelMat[,2:length(SelMat)], 2, function(x) sum(x*RawDat$CL)/sum(x))

# Plot mean length
plot(rslts$year, rslts$Lbar, type='p', lwd=2, xlab="", ylab="L & Lsel", main="nep2829_males", ylim=c(min(rslts$Lbar),max(rslts$Lsel)))
points(rslts$year, rslts$Lsel, lwd=2, col='blue')
abline(h=mean(rslts$Lbar), lwd=2)
abline(h=mean(rslts$Lsel), lwd=2, col='blue')

# Plot mean length against SPR
rslts$DSelMat <- rslts$SL50-L50
rslts$DSelMat[rslts$year==2011] <- NA
plot(rslts$DSelMat, rslts$SPR, lwd=2, xlab="LS50-L50", ylab="SPR", main="nep2829_males")
fit <- lm(rslts$SPR ~ rslts$DSelMat)
summary(fit)$r.squared
lines(rslts$DSelMat[-12], predict(fit), lwd=2)
text(-0.5,0.55, paste("r^2=", round(summary(fit)$r.squared, digits=3)), cex=.8)

# Plot change mean selected length against change SPR
SPRy <- rslts$SPR[-c(1,12)] # get rid of first year and 2011
SPRy1 <- rslts$SPR[-c(12,15)] # get rid of 2011 and last year
DSPR <- SPRy-SPRy1
Lsely <- rslts$Lsel[-c(1,12)]
Lsely1 <- rslts$Lsel[-c(12,15)]
DLsel <- Lsely-Lsely1
plot(DLsel, DSPR, lwd=2, xlab="DLsel", ylab="DSPR", main="nep2829_males")
fit <- lm(DSPR ~ DLsel)
summary(fit)$r.squared
lines(DLsel, predict(fit), lwd=2)
text(-2.5,0.15, paste("r^2=", round(summary(fit)$r.squared, digits=2)), cex=.8)
