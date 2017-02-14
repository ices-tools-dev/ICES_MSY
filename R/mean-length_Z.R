#################################################
# WK Proxy
# Mean Length Estimates of Z
# 11/3/15 - 11/6/15
#################################################

# Overview

# This script estimates mortality using mean-length-based methods:
#     Gedamke Hoenig (2006) - non-equlibrium extension of Beverton-Holt mortality
#     THoG - incorporates effort into the Gedamke-Hoenig formulation
# This is a general script and can be applied to data that follows the attached
# format. Data may need to  be reformatted into the appropriate manner before
# running.
# As an example, data for nep-2829 is attached and can be tried.


#################################################
# Run all of this first
#################################################

# path <- "C:/path"
# setwd(path)

library(fishmethods)
library(dplyr)
library(lattice)
library(data.table)
library(reshape2)

# Modification for discrete (annual) reproduction
source('/R/gedamke.hoenig_BFGS_discrete.r')
# THoG model
source('/R/mleneffort.r')
# YPR estimator
source('/R/ypr_ml.r')

# Function to calculate mean lengths
neph_Lbar <- function(input, Lc){
  sub <- subset(input, Length >= Lc & !(value==0))
  sub <- data.table(sub)
  fin <- sub[, list(
    mlen = sum(value * Length) / sum(value),
    ss_bins = length(value),
    ss = 1
  ),
  by = list(Year)]
  fin$Year <- as.numeric(fin$Year)
  out <- fin
}

#################################################
# Read in the data and format it
#################################################

input_lenfreq <- read.csv('/data/nep/LFD/nep_LFD_males.csv')
input_effort <- read.csv('/data/nep/LFD/nep_effort.csv')
# input_lenfreq <- read.csv('length_freq_file.csv')
# input_effort <- read.csv('effort_file.csv')
colnames(input_lenfreq)[1] <- "Length"
input_effort <- input_effort[input_effort$Year %in% c(2000:2014),]

name <- "nep"
sex <- "M"
K <- 0.2
Linf <- 70
max_age <- 10 # just a guess

input_melt <- melt(input_lenfreq, id.vars="Length", variable.name="Year")
input_melt$Year <- as.numeric(substr(input_melt$Year, 2, 5))


#################################################
# Examine histograms
#################################################

stockname <- paste(name, sex)

par(mfrow=c(1, 1))
input_tot <- rowSums(input_lenfreq[, -1])
input_lengths <- input_lenfreq[,1]

# Length of full selectivity (Lc)
Lc <- input_lengths[which.max(input_tot)] ; Lc
barplot(input_tot, names.arg=input_lengths, 
        main=paste(stockname, "Length Frequency Distribution \nPeak=", Lc))

# If the distribution has more than one 'peak', use your judgement for what the
# Lc should be, and change it accordingly.
# Lc <- NULL

# Length frequency for each year
vert_line <- function(x, y, ...) {
  panel.abline(v = Lc, col = "red")
  panel.xyplot(x, y, ...)
}
xyplot(value ~ Length | as.factor(Year), input_melt, panel = vert_line, typ = "l", 
       main=paste(stockname, "Length Frequency Distribution \nPeak=", Lc), as.table=T)

#################################################
# Plot peaks for each year
#################################################

input_meltdt <- data.table(input_melt)
Lc_year <- input_meltdt[, list(
  peak = Length[which.max(value)]
),
by = list(Year)]

mean(Lc_year$peak)

plot(Lc_year, type='b', pch=19, main=paste(stockname, "Peak Over Time"))

#################################################
# Calculate mean lengths
#################################################

# If you want to run the model for more than one Lc (which is suggested for
# diagnostic purposes), you can input the different Lcs here.

ml_gh <- neph_Lbar(input=input_melt, Lc=Lc)

# ml_gh2 <- neph_Lbar(input=input_melt, Lc=Lc2)



#################################################
# Gedamke-Hoenig model
#################################################

# Annual recruitment: use the function gh()
# Continuous recruitment: use the function bhnoneq() in the fishmethods package.

# You need to put in your own starting guesses for year of change and Zs.
# Error messages mean that you should try a different guess.

# Single Z
gh(year = ml_gh$Year, 
   mlen = ml_gh$mlen, 
   ss = ml_gh$ss, 
   K = K, 
   Linf = Linf, 
   Lc = Lc, 
   nbreaks = 0,
   styrs = c(), 
   stZ = c(.3))

# Two Zs
# One change: Z = .48 -> .44, aic = 37
gh(year=ml_gh$Year,
   mlen = ml_gh$mlen, 
   ss = ml_gh$ss, 
   K = K, 
   Linf = Linf, 
   Lc = Lc, 
   nbreaks = 1,
   styrs = c(2007),
   stZ = c(.9, .5))

# Three Zs
# two change: Z = .41 -> .67 -> .43, aic = 37 
gh(year = ml_gh$Year,
   mlen = ml_gh$mlen, 
   ss = ml_gh$ss,
   K = K, 
   Linf = Linf, 
   Lc = Lc, 
   nbreaks = 2,
   styrs = c(2002, 2006),
   stZ = c(.3, .6, .3))

#################################################
# THoG model
#################################################

# You need to put starting values for q and M.
# The value of q will depend on the effort values, so if you keep getting errors
# then change the magnitude of your starting q (or effort).

# For continuous recruitment, change to n_cohort=12.

# If you think that there was no effort prior to the start of the time series,
# replace init_effort=rep(0, max_age).

ml_thog <- left_join(input_effort, 
                     ml_gh, 
                     by = "Year")
ml_thog$ss_bins[is.na(ml_thog$ss_bins)] <- 0
ml_thog$ss[is.na(ml_thog$ss)] <- 0

# Estimate M
max_age <- 10
mlen_effort(year = ml_thog$Year,
            mlen = ml_thog$mlen,
            effort = ml_thog$effort, 
            ss = ml_thog$ss,
            K = K,
            Linf = Linf, 
            Lc = Lc, 
            n_cohort = 1,
            n_age = max_age, 
            init_effort = rep(ml_thog$effort[1], max_age),
            stq = c(1),
            stM = c(.3))

# Fix M
mlen_effort(year = ml_thog$Year, 
            mlen = ml_thog$mlen,
            effort = ml_thog$effort, 
            ss = ml_thog$ss,
            K = K, 
            Linf = Linf, 
            Lc = Lc, 
            n_cohort = 1,
            n_age = max_age, 
            init_effort = rep(ml_thog$effort[1], max_age),
            stq = c(1), 
            stM = c(.3), est.M = F)

######################################
# YPR estimation for reference points
######################################

# Three Zs
# two change: Z = .41 -> .67 -> .43, aic = 37 
MLZ(year = ml_gh$Year,
   mlen = ml_gh$mlen, 
   ss = ml_gh$ss,
   K = K, 
   Linf = Linf, 
   t0 = -0.15,
   Lc = Lc, 
   nbreaks = 2,
   styrs = c(2002, 2006),
   stZ = c(.3, .6, .3),
   M = 0.3,
   a = 0.00043,
   b = 3.160)
