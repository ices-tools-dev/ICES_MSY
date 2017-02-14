# YPR - function to calculate F-01 from life history parameters from natural mortality M, 
#       growth parameters (von Bertalanffy: Linf, K, t0; length-weight: a, b), and Lc.
#       
#       Maxage of YPR is currently set to be age at which 1% of cohort survives given M
#       maxF is the maximum F for searching for F01
#       nsteps is the increment of F in the YPR search

YPR <- function(Linf, K, t0, M, a, b, Lc, maxage = -log(0.01)/M, maxF = 2, nsteps = 0.01, 
                graph = TRUE) {
  
  maxage <- as.integer(maxage)
  frates <- seq(0, maxF, nsteps) # range of F for YPR
  tc <- -log(1 - Lc/Linf)/K + t0 # convert Lc to age, based on von Bertalanffy parameters
  tc <- round(tc, 0)
  tc[tc < 1] <- 1
  tc[tc > maxage] <- maxage
  
  #average weight at age - follow von Bertalanffy growth
  age <- 1:maxage
  la <- Linf * (1 - exp(-K*(age - t0)))
  wa <- a * la^b
  
  #vulnerability schedule - assumes knife-edge vulnerability, where all individuals age tc to maxage are fully vulnerbale
  #all individulas less than age tc are not vulnerable
  vul <- ifelse(age < tc, 0, 1)
  
  lx <- numeric(maxage)
  ypr <- numeric(length(frates))
  
  lx[1]<-1
  for(k in 1:length(frates)) {
    for(i in 2:maxage) {
      lx[i] <- lx[i-1] * exp(-(M + vul[i-1] * frates[k]))
    }
    
    phi_vb <- sum(lx*wa*vul) # Equilibrium biomass per recruit
    
    ypr[k] <- (1 - exp(-frates[k])) * phi_vb
  }
  
  # More code that derived F0.1 in 'per recruit analysis.R' (Meaghan Bryan)
  slope.origin <- (ypr[2] - ypr[1])/(frates[2] - frates[1])
  slope.10 <- 0.1*slope.origin
  
  slope <- numeric(length(ypr))
  slope[1] <- slope.origin
  for(i in 3:length(ypr)) slope[i-1] <- (ypr[i] - ypr[i-1])/(frates[i] - frates[i-1])
  dif <- abs(slope - slope.10)
  dif[is.na(dif)] <- 10e10
  
  F01 <- frates[which.min(dif)]
  
  output <- list(summary = data.frame(frates = frates, YPR = ypr), F01 = F01)
  
  if(graph) {
    plot(YPR ~ frates, output$summary, xlab = "Fishing Mortality (F)", ylab = "YPR", typ = "l")
    points(F01, ypr[which.min(dif)], col = "red", pch = 16)
  }
  
  return(output)
}


# MLZ - wrapper function to use mean length estimator and YPR to obtain reference point and F/Fmsy.
#
#       Maxage of YPR is currently set to be age at which 1% of cohort survives given M
#       maxF is the maximum F for searching for F01
#       nsteps is the increment of F in the YPR search

MLZ <- function(year, mlen, ss, K, Linf, t0, Lc, nbreaks, styrs, stZ, 
                M, a, b, maxage = -log(0.01)/M, maxF = 2, nsteps = 0.01, graph = TRUE) {
  
  # estimate F from mean length estimator
  Z.estimator <- gh(year, mlen, ss, K, Linf, Lc, nbreaks, styrs, stZ, graph)
  Z.recent <- Z.estimator$summary$Estimate[nbreaks + 1]
  F.benchmark <- Z.recent - M
  
  if(F.benchmark <= 0) stop("F = Z - M results in F < 0.")
  
  # obtain F-01 from YPR
  YPR.analysis <- YPR(Linf, K, t0, M, a, b, Lc, maxage, maxF, nsteps, graph)
  F01 <- YPR.analysis$F01
  
  if(F01 == maxF) warning("F-01 was the maximum allowed F. Does F-01 exist?")
  
  # F/Fmsy
  Fstatus <- data.frame(F.benchmark, F01, F.benchmark, F01)
  names(Fstatus) <- c("F", "F-01", "F/Fmsy")
  
  output <- list(Zestimator = Z.estimator, YPR = YPR.analysis, Status = Fstatus)
  return(output)
}


