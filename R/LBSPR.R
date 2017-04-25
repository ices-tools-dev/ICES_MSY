rm(list = ls())

# install.packages("LBSPR")
library(LBSPR)

dat_path <- "data/nep/LFD/nep_lenFreq.csv"

MyPars <- new("LB_pars")
MyPars@Linf <- 70
MyPars@L50 <- 28.4
MyPars@L95 <- 28.5
MyPars@MK <- 1.5 
MyPars@CVLinf <- 0.1
MyPars@BinWidth <- 1
MyPars@Wbeta <- 3.2229

Len1 <- new("LB_lengths", 
            LB_pars = MyPars, 
            file = dat_path, 
            dataType = "freq",
            header = TRUE)

ModelFit <- LBSPRfit(MyPars, Len1)


Results <- matrix(c(ModelFit@SL50, ModelFit@SL95, ModelFit@FM, ModelFit@SPR),
                  ncol=4, byrow=FALSE)

# 95% confidence intervals #
CIlower <- Results[,1:4] - 1.96 * sqrt(ModelFit@Vars)
CIupper <- Results[,1:4] + 1.96 * sqrt(ModelFit@Vars)
CIlower[!is.finite(CIlower)] <- 0
CIupper[!is.finite(CIupper)] <- 0
CIlower[CIlower <0 ] <- 0
CIupper[CIupper <0 ] <- 0

# correct bounded parameters - dodgy I know!
CIlower[CIlower[,3]<0,3] <- 0
CIupper[CIupper[,4]>1,4] <- 1
CIlower[CIlower[,4]<0,4] <- 0

CIlower <- round(CIlower,2)
CIupper <- round(CIupper,2)

DF <- data.frame(Years = ModelFit@Years,
                 SPR = round(ModelFit@SPR, 2),
                 SPR_l =  CIlower[,4],
                 SPR_h = CIupper[,4],
                 SL50 = round(ModelFit@SL50, 2),
                 SL50_l = CIlower[,1],
                 SL50_h = CIupper[,1], 
                 SL95 = round(ModelFit@SL95, 2),
                 SL95_l = CIlower[,2],
                 SL95_h = CIupper[,2],
                 FM = round(ModelFit@FM, 2),
                 FM_l = CIlower[,3],
                 FM_h = CIupper[,4])