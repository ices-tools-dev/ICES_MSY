rm(list = ls())
###########################
# Length-based indicators #
###########################
# Adapted from T. Miethe and C. Silva, WKLIFE-V, Oct2015
library(reshape2)
source("https://raw.githubusercontent.com/ices-tools-dev/LBIndicator_shiny/master/utilities.R")

wal <- url("https://raw.githubusercontent.com/ices-tools-dev/LBIndicator_shiny/master/data/walNeph.csv")
freq <- url("https://raw.githubusercontent.com/ices-tools-dev/LBIndicator_shiny/master/data/freqNeph.csv")

wal <- read.csv(wal, 
                stringsAsFactors = FALSE)

freq <- read.csv(freq, 
                 stringsAsFactors = FALSE)


DF <- lb_ind(data =  freq, 
             binwidth = 5,
             linf = 70, 
             lmat = 28,
             mk_ratio = 1.5, # m/k ratio
             wal) 