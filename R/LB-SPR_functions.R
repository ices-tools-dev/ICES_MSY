## PACKAGES ###
library(numDeriv)

## FUNCTIONS ###
# Length-based SPR Function 
LBSPRFunc <- function(MK = 1.5,
                      Linf = 100, 
                      CVLinf = 0.1, 
                      L50 = 55, 
                      L95 = 60, 
                      Beta = 3, 
                      FM = 1, 
                      SL50 = 30, 
                      SL95 = 35, 
                      By = 2, 
                      MaxLMult = 1.25, 
                      P = 0.01, 
                      Nage = 100) {
  
  #LenMids <- seq(from=By *0.5, to= Linf * MaxLMult, by=By)
  #LenBins <- seq(from=0, by=By, length.out=length(LenMids)+1)	
  x <- seq(from=0, to=1, length.out=Nage) # relative age vector
  EL <- (1-P^(x/MK)) * Linf # length at relative age 
  rLens <- EL/Linf # relative length 
  SDL <- EL * CVLinf # standard deviation of length-at-age
  
  Nlen <- length(LenMids) 
  Prob <- matrix(NA, nrow=Nage, ncol=Nlen)
  Prob[,1] <- pnorm((LenBins[2] - EL)/SDL, 0, 1) # probablility of length-at-age
  for (i in 2:(Nlen-1)) {
    Prob[,i] <- pnorm((LenBins[i+1] - EL)/SDL, 0, 1) - pnorm((LenBins[i] - EL)/SDL, 0, 1)
  }
  Prob[,Nlen] <- 1 - pnorm((LenBins[Nlen] - EL)/SDL, 0, 1)
  
  SL <- 1/(1+exp(-log(19)*(LenMids-SL50)/(SL95-SL50))) # Selectivity at length
  Sx <- apply(t(Prob) * SL, 2, sum) # Selectivity at relative age 
  MSX <- cumsum(Sx) / seq_along(Sx) # Mean cumulative selectivity for each age 
  Ns <- (1-rLens)^(MK+(MK*FM)*MSX) # number at relative age in population
  
  Cx <- t(t(Prob) * SL) # Conditional catch length-at-age probablilities
  # for (X in seq_along(x)) {
  # if (sum(Cx[X,]) > 0 )Cx[X,] <- Cx[X,]/sum(Cx[X,])
  # if (sum(Cx[X,]) == 0 )Cx[X,] <- 0  
  # }
  
  Nc <- apply(Ns * Cx, 2, sum)
  
  Ml <- 1/(1+exp(-log(19)*(LenMids-L50)/(L95-L50))) # Maturity at length
  Ma <-  apply(t(Prob) * Ml, 2, sum) # Maturity at relative age 
  
  N0 <- (1-rLens)^MK # Unfished numbers-at-age 
  SPR <- sum(Ma * Ns * rLens^Beta)/sum(Ma * N0 * rLens^Beta)
  
  Output <- NULL 
  Output$SPR <- SPR 
  Output$LenMids <- LenMids
  Output$PropLen <- Nc/sum(Nc)
  return(Output)
}

# Function wrapper
LBSPRFuncwrap <- function(par, 
                          MK = 1.5, 
                          Linf = 100, 
                          CVLinf = 0.1, 
                          L50 = 55, 
                          L95 = 60, 
                          Beta = 3, 
                          By = 2,
                          MaxLMult = 1.25, 
                          P = 0.01, 
                          Nage = 100) {
  estFM <- exp(par[1])  
  estSL50 <- exp(par[2]) * Linf
  estSL95 <- estSL50 + (exp(par[3]) * estSL50	)
  
  LBSPRFunc(MK=MK, Linf=Linf, CVLinf=CVLinf, L50=L50, L95=L95,Beta=Beta, FM=estFM, SL50=estSL50, SL95=estSL95, By=By, MaxLMult=MaxLMult, P=P, Nage=Nage)$SPR
}

# Optimisation Function   
OptRoutine <- function(Pars, LenDat, MK=1.5, Linf=100, CVLinf=0.1, L50=55, L95=60, Beta=3, By=2, MaxLMult=1.25, P=0.01, Nage=100) {
  FM <- exp(Pars[1])
  SL50 <- exp(Pars[2]) * Linf
  SL95 <- SL50 + (exp(Pars[3]) * SL50)
  runMod <- LBSPRFunc(MK, Linf, CVLinf, L50, L95, Beta, FM, SL50, SL95, By, MaxLMult, P, Nage)
  LenProb <- LenDat/sum(LenDat)
  ind <- (LenProb > 0)
  return(-sum(LenDat[ind] * log(runMod$PropLen[ind]/LenProb[ind])))
} 

# Delta method function
delta.meth.1 <- function(mle, info, f, MK, Linf, CVLinf, L50, L95, 
                         Beta, By, MaxLMult, P, Nage)
{
  theta.hat <- f(mle,  MK=MK, Linf=Linf, CVLinf=CVLinf, L50=L50, L95=L95, 
                 Beta=Beta, By=By, MaxLMult=MaxLMult, P=P, Nage=Nage)
  epsilon <- 1.0*10^(-6)
  alpha1 <- mle
  alpha2 <- mle
  alpha3 <- mle
  alpha1[1] <- alpha1[1] + epsilon
  alpha2[2] <- alpha2[2] + epsilon
  alpha3[3] <- alpha3[3] + epsilon 
  d1 <- ( f(alpha1, MK=MK, Linf=Linf, CVLinf=CVLinf, L50=L50, L95=L95, 
            Beta=Beta, By=By, MaxLMult=MaxLMult, P=P, Nage=Nage) - theta.hat )/epsilon
  d2 <- ( f(alpha2,MK=MK, Linf=Linf, CVLinf=CVLinf, L50=L50, L95=L95, 
            Beta=Beta, By=By, MaxLMult=MaxLMult, P=P, Nage=Nage) - theta.hat )/epsilon
  d3 <- ( f(alpha3, MK=MK, Linf=Linf, CVLinf=CVLinf, L50=L50, L95=L95, 
            Beta=Beta, By=By, MaxLMult=MaxLMult, P=P, Nage=Nage) - theta.hat )/epsilon
  part.der <- c(d1,d2,d3)
  theta.var <- part.der %*% solve(info) %*% part.der
  answer <- theta.var
  answer
}

# Selectivity at length function
SLen <- function(i,SL50,SL95){
  sel <- 1/(1+(exp(-log(19)*(i-SL50)/(SL95-SL50))))
  return(sel)
}