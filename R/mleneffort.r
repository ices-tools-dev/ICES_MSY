# mlen_effort is the function to implement Amy's model

# Year - vector of years in data
# Mlen - vector of mean lengths (can have NA)
# Effort - vector of effort (must have data for all years)
# Linf, K, t0, Lc - numeric, t0 by default = 0 because it is often not reported
# est.M - boolean as to whether model should estimate M (natural mortality) or fix M
# stq - starting value for catchability coefficient. Scales effort to esimate F (fishing mortality).
# stM - starting value for M, if est.M = F, then this is the value of M assumed in model
# n_age - number of ages above Lc to model
# n_cohort = number of cohorts per year. Set = 1 for annual reproduction, = 12 for monthly 
#            (approximates continuous)
# init_effort - a vector of length n_age*n_cohort of effort prior to first year of time series.
#               For example, for virgin conditions prior to year 1, use: rep(0, n_age)

mlen_effort <- function(year, mlen, effort, ss, Linf, K, t0 = 0, Lc, est.M = TRUE, stq, stM, n_age,
                        n_cohort = 12, init_effort, graph = TRUE) {
    
  if(est.M) {
    results <- optim(c(stq, stM), return_LL_estM, mlen = mlen, effort = effort, ss = ss,
                     Linf = Linf, K = K, t0 = t0, Lc = Lc, n_age = n_age,
                     n_cohort = n_cohort, init_effort = init_effort, method = "BFGS", hessian = TRUE)
  }
  else results <- optim(stq, return_LL_fixM, mlen = mlen, effort = effort, ss = ss,
                        Linf = Linf, K = K, t0 = t0, Lc = Lc, n_age = n_age,
                        n_cohort = n_cohort, init_effort = init_effort, M = stM, method = "BFGS", hessian = TRUE)
  
  # Process results
  var <- solve(results$hessian)
  SE <- sqrt(diag(var))
  
  if(est.M) {
    summary.par <- data.frame(Variable = c("q", "M", "sigma"), 
                              Estimate = c(results$par, sig),
                              Std.dev = c(SE, NA))
  }
  else summary.par <- data.frame(Variable = c("q", "sigma"), 
                                 Estimate = c(results$par, sig),
                                 Std.dev = c(SE, NA))
  
  predicted <- data.frame(year = year, observed = mlen, predicted = Lpred,
                          residual = mlen - Lpred)
  
  if(results$convergence == 0) 
    status <- paste("Optimizer converged after", results$counts[1], 
                    "function calls.")
  else status <- "Optimizer did not converge."
  
  if (all(eigen(results$hessian)$values > 0)) 
    hess.status <- "Hessian matrix is positive definite."
  else hess.status <- "Hessian matrix is not positive definite. Check parameter estimates."
  
  if(est.M) cor.qM <- cov2cor(solve(results$hessian))[1, 2]
  else cor.qM <- "M is fixed."
  
  # Draw figure
  if(graph) {
    par(mfrow = c(3, 1))
    plot(year, mlen, xlab = "Year", ylab = "Mean Length", 
         typ = "o", pch = 16)
    lines(year, Lpred, col = "red")
    plot(year, mlen - Lpred, xlab = "Year", ylab = "Residual", 
         typ = "o", pch = 16)
    abline(h = 0, col = "red")
    plot(year, results$par[1] * effort, xlab = "Year", ylab = "Fishing Mortality (F)",
         typ = "o", pch = 16, col = "red")
  }
  
  list(summary = summary.par, convergence = status, hessian = hess.status, 
       correlation = cor.qM,
       results = data.frame(year = year, observed = mlen, predicted = Lpred, 
                            residual = mlen - Lpred, F = results$par[1] * effort))
}


return_LL_estM <- function(x, mlen, ss, effort, Linf, K, Lc, t0, n_age, n_cohort, init_effort) {
  
  q <- x[1]
  M <- x[2]
  
  age_min <- t0 - log(1 - Lc/Linf)/K
  
  age.vector <- seq(age_min, age_min + n_age, by = 1/n_cohort)[1:(n_age*n_cohort)]
  length.vector <- Linf * (1 - exp(-K*(age.vector-t0)))
  
  n_yr <- length(mlen)
  n_data <- sum(ss != 0)
  
  ####### Initial Z
  Zinit <- numeric(length = n_age * n_cohort)
  for(a in 1:n_age) {
    for(k in 1:n_cohort) Zinit[(a-1)*n_cohort+k] <- (q*init_effort[a] + M)/n_cohort
  }
  
  ####### Z during time series
  Z <- q * effort + M
  
  ####### Create abundance array
  Nt <- array(data = 1, dim = c(n_yr, n_cohort*n_age, n_cohort))
  NtLt <- array(data = 1, dim = c(n_yr, n_cohort*n_age, n_cohort))
  Lpred <<- numeric(length = n_yr)
  
  ######## First year of mortality
  # First cohort of the first year
  Nt[1, 1, 1] <- 1
  for(a in 1:(n_age*n_cohort-1)) Nt[1, a+1, 1] <- Nt[1, a, 1] * exp(-Zinit[a])
  
  # Subsequent cohorts of the first year
  if(n_cohort > 1) {
    for(k in 1:(n_cohort-1)) {
      Nt[1, 1, k+1] <- 1
      for(a in 1:(n_age*n_cohort-1)) Nt[1, a+1, k+1] <- Nt[1, a, k] * exp(-Z[1]/n_cohort)
    }
  }
  
  ######## Mortality schedule for all other years
  for(t in 1:(n_yr-1)) {
    Nt[t+1, 1, 1] <- 1
    
    for(a in 1:(n_age*n_cohort - 1)) Nt[t+1, a+1, 1] <- Nt[t, a, n_cohort] * exp(-Z[t]/n_cohort)
    
    if(n_cohort > 1) {
      for(k in 1:(n_cohort-1)) {
        Nt[t+1, 1, k+1] <- 1
        
        for(a in 1:(n_age*n_cohort - 1)) Nt[t+1, a+1, k+1] <- Nt[t+1, a, k] * exp(-Z[t]/n_cohort)
      }
    }
  }
  
  for(t in 1:n_yr) {
    NtLt[t, , ] <- Nt[t, , ] * length.vector
    Lpred[t] <<- sum(NtLt[t, , ])/sum(Nt[t, , ])
  }
  
  sig <<- sqrt(sum(ss * (mlen - Lpred)^2, na.rm = TRUE)/n_data)
  
  neg.LL <- -n_data * log(sig) - 0.5 * sum(ss * (mlen - Lpred)^2, na.rm = TRUE) / (sig^2)
  
  -1 * neg.LL
  
}


return_LL_fixM <- function(x, mlen, ss, effort, Linf, K, Lc, t0, n_age, n_cohort, init_effort, M) {
  
  q <- x
  
  age_min <- t0 - log(1 - Lc/Linf)/K
  
  age.vector <- seq(age_min, age_min + n_age, by = 1/n_cohort)[1:(n_age*n_cohort)]
  length.vector <- Linf * (1 - exp(-K*(age.vector-t0)))
  
  n_yr <- length(mlen)
  n_data <- sum(ss != 0)
  
  ####### Initial Z
  Zinit <- numeric(length = n_age * n_cohort)
  for(a in 1:n_age) {
    for(k in 1:n_cohort) Zinit[(a-1)*n_cohort+k] <- (q*init_effort[a] + M)/n_cohort
  }
  
  ####### Z during time series
  Z <- q * effort + M
  
  ####### Create abundance array
  Nt <- array(data = 1, dim = c(n_yr, n_cohort*n_age, n_cohort))
  NtLt <- array(data = 1, dim = c(n_yr, n_cohort*n_age, n_cohort))
  Lpred <<- numeric(length = n_yr)
  
  ######## First year of mortality
  # First cohort of the first year
  Nt[1, 1, 1] <- 1
  for(a in 1:(n_age*n_cohort-1)) Nt[1, a+1, 1] <- Nt[1, a, 1] * exp(-Zinit[a])
  
  # Subsequent cohorts of the first year
  if(n_cohort > 1) {
    for(k in 1:(n_cohort-1)) {
      Nt[1, 1, k+1] <- 1
      for(a in 1:(n_age*n_cohort-1)) Nt[1, a+1, k+1] <- Nt[1, a, k] * exp(-Z[1]/n_cohort)
    }
  }
  
  ######## Mortality schedule for all other years
  for(t in 1:(n_yr-1)) {
    Nt[t+1, 1, 1] <- 1
    
    for(a in 1:(n_age*n_cohort - 1)) Nt[t+1, a+1, 1] <- Nt[t, a, n_cohort] * exp(-Z[t]/n_cohort)
    
    if(n_cohort > 1) {
      for(k in 1:(n_cohort-1)) {
        Nt[t+1, 1, k+1] <- 1
        
        for(a in 1:(n_age*n_cohort - 1)) Nt[t+1, a+1, k+1] <- Nt[t+1, a, k] * exp(-Z[t]/n_cohort)
      }
    }
  }
  
  for(t in 1:n_yr) {
    NtLt[t, , ] <- Nt[t, , ] * length.vector
    Lpred[t] <- sum(NtLt[t, , ])/sum(Nt[t, , ])
  }
  
  sig <<- sqrt(sum(ss * (mlen - Lpred)^2, na.rm = TRUE)/n_data)
  
  neg.LL <- -n_data * log(sig) - 0.5 * sum(ss * (mlen - Lpred)^2, na.rm = TRUE) / (sig^2)
  
  return(-1 * neg.LL)
  
}



