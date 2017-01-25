rm(list = ls())
###########################
# Length-based indicators #
###########################
# Adapted from T. Miethe and C. Silva, WKLIFE-V, Oct2015

####################
# Install Packages #
####################
library(reshape2)
library(lattice)
source("~/ICES_MSY/R/ICE_LFD.R")

####################
# Input parameters #
####################
# Sexes for which analyses have to be performed on, M males, F females, N unsexed
S <- c("M","F")

# Life history parameters (M, F, N)
Linf <- c(70, 65) 
Lmat <- c(28.4, 30)

# Years of available data
startyear <- 2000
endyear <- 2014
stock = "nep"

##########################################
## Load data: Numbers-at-length per year #
##########################################

# 1st column name MeanLength, 2nd and remaining columns named according to the years.
# Use the same format for the MeanWeight-at-length dataframe

#MeanLength 2002  2003  2004 ...
#20.5        0     0     0   ...
#21.5        1.24  0.2   0.8 ...
#22.5        5.1   3.1   5.2 ...
#...
# define
# m, f, ns - dataframes names for numbers for males, females or unsexed 
# mw, fw, nsw - dataframes names for mean weights for males, females or unsexed

# ICE_LFD() will format InterCatch Exchange (what submitters use) data files.
ic_list <- ICE_LFD(file_path = "~/ICES_MSY/data/nep/",
                   startyear = startyear,
                   endyear = endyear,
                   country = "XYZ",
                   stock = stock,
                   ctype = "L",# for landings L, discard D
                   ftype = ".dat", # file type (should be .dat)
                   S = S)

m <- ic_list$male_length
f <- ic_list$female_length
ns <- ic_list$unsexed_length

mw <- ic_list$male_weight
fw <- ic_list$female_weight
nsw <- ic_list$unsexed_weight

# Choose bin size (ClassInt) width for estimation of Lc

######################################################################################################

length_data <- m
MeanLengthLB = FALSE # TRUE  if the MeanLength is not the class midpoint but the class lower bound and bin size 1
sex = "M"

if(MeanLengthLB){
  length_data$MeanLength <- length_data$MeanLength + 0.5
}

############################################################
# step 1 check length distribution plots to decide whether 
# regrouping is necessary to determine Lc (Length at first catch= 50% of mode)
############################################################

length_plot <- function(length_data, 
                        ClassInt = 2,
                        filename = "~/length_distribution.png",
                        save_plot = FALSE) {


  df0 <- length_data

  df0.long <- melt(df0, id.vars = 'MeanLength')
  df0.long$variable <- gsub("X", "", as.character(df0.long$variable))
  
  minCL <- floor((min(df0$MeanLength) - .5) / ClassInt) * ClassInt  # original data 1mm length class
  maxCL <- ceiling((max(df0$MeanLength) + .5) / ClassInt) * ClassInt
  
  df0$LC <- cut(df0$MeanLength, 
                breaks = seq(minCL,
                             maxCL, 
                             ClassInt),
                include.lowest = T)
  
  df0.gr <- aggregate(df0[, 2:ncol(df0)-1], by=list(df0$LC), sum)
  names(df0.gr)[1] <- 'lclass'
  df0.gr <- cbind(lclass=df0.gr$lclass,
                  lmidp = as.numeric(substr(df0.gr$lclass, 
                                            2, 3)) + ClassInt / 2,
                  df0.gr[, 3:ncol(df0.gr)])
  df0.gr.long <- melt(df0.gr[ ,-1], id.var = 'lmidp')
  df0.gr.long$variable <- gsub("X", "", as.character(df0.gr.long$variable))
  
  names(df0.gr.long)[2:3] <- c('year', 'Number')
  df0.gr.long$year <- as.numeric(as.character(df0.gr.long$year))
  
  
  length_bars <- barchart(Number ~ lmidp|as.factor(year),
                          data = df0.gr.long, 
                          horizontal = F, 
                          as.table = T, 
                          ylim = c(0, NA),
                          xlab = 'Length',
                          ylab = 'Number (10^3)',
                          scales = list(x = list(at = seq(2,
                                                          length(unique(df0.gr.long$lmidp)), 4),
                                                 labels = seq(min(df0.gr.long$lmidp) + ClassInt, 
                                                              max(df0.gr.long$lmidp), 4 * ClassInt))),
                          main = paste0("Length class: ", ClassInt, " mm"),
                          cex.main = 1.2)

  if(save_plot) {
  png(filename = filename, 
      bg = "white",
      pointsize = 5,
      units = "cm",
      width = 35, 
      height = 18, 
      res = 600)
  
  print(length_bars)
  dev.off()
  } else {
    print(length_bars)
  }
}

length_plot(length_data,
            ClassInt = 2,
            save_plot = FALSE)

# step 6 final decision on regrouping fill in! 
  ClassInt <- 2


######################################################################################################################
# step 2 Calculate indicators per sex

  
  Year <- c(startyear:endyear)
  Year1 <- c((startyear+1):endyear)
  
  for(s in 1:length(S)){
    sex <- S[s] 
    if(sex=="M") final <- m  #numbers
    if(sex=="F") final <- f
    if(sex=="N") final <- ns
    
    if(sex=="M") weight <- mw   #mean weights
    if(sex=="F") weight <- fw
    if(sex=="N") weight <- nsw
    
    Ind <- data.frame(matrix(ncol=24, nrow=endyear-startyear+1)) 
    names(Ind) <- c('Year','L75','L25','Lmed', 'L90', 'L95', 'Lmean','Lc','LFeM','Lmaxy' ,'Lmat', 'Lopt','Linf', 'Lmax5',  'Lmean_LFeM','Lc_Lmat','L25_Lmat','Lmean_Lmat','Lmean_Lopt', 'L95_Linf', 'Lmaxy_Lopt','Lmax5_Linf','Pmega','Pmegaref')
    Ind$Year <- startyear:endyear
    
    #  regrouping with selected length class width
    df0 <- final
    minCL <- floor((min(df0$MeanLength)-.5)/ClassInt)*ClassInt  #originaldat 1mm length class
    maxCL <- ceiling((max(df0$MeanLength)+.5)/ClassInt)*ClassInt
    df0$LC <- cut(df0$MeanLength, breaks=seq(minCL,maxCL,ClassInt), include.lowest=T)
    df0.gr <- aggregate(df0[,2:ncol(df0)-1], by=list(df0$LC), sum)
    names(df0.gr)[1] <- 'lclass'
    df0.gr <- cbind(lclass=df0.gr$lclass, lmidp=as.numeric(substr(df0.gr$lclass,2,3))+ClassInt/2, df0.gr[,3:ncol(df0.gr)])
    df0.gr.long <- melt(df0.gr[,-1], id.var='lmidp')
    names(df0.gr.long)[2:3] <- c('year', 'Number')
    df0.gr.long$year <- as.numeric(as.character(df0.gr.long$year))
    
    df0.gr.long <- melt(df0.gr[,-1], id.var='lmidp')
    names(df0.gr.long)[2:3] <- c('year', 'Number')
    df0.gr.long$year <- as.numeric(as.character(df0.gr.long$year))
    res <- data.frame(year=min(as.numeric(df0.gr.long$year)):max(as.numeric(df0.gr.long$year)), lmidp=NA, nmax=NA, lc=NA)
    
    for (j in 3:ncol(df0.gr)) {
      for (i in 2:nrow(df0.gr)) {
        if(df0.gr[i+1,j]-df0.gr[i,j]>=0) {
          next
        } else {
          res$lmidp[j-2] = df0.gr$lmidp[i]
          res$nmax[j-2] = df0.gr[i,j]
          a = res$nmax[j-2]/2
          df1 = df0.gr[,c(2,j)]
          for (k in 1:nrow(df1)) {
            if (df1[k,2] < a) {
              next
            } else {  
              res$lc[j-2] = df1[k,1]
            }
            break
          }
        }
        break
      }
    }
    
    Ind$Lc <- res$lc
    Ind$Lmat <- Lmat[s]
    Ind$Lopt <- 2/3*Linf[s]
    Ind$Linf <- Linf[s]
    
    for(jj in (1:length(Year))+1){
      j <- jj-1 
      
      final2 <- final[,c(1,jj)]
      colnames(final2) <- c("lngth","number")
      
      final2$cumsum <- cumsum(final2[,2])
      final2$cumsum_perc <- final2$cumsum/sum(final2$number)
      
      # find mean top 5%
      numb <- as.data.frame(final2[rev(order(final2$lngth)),"number"])    # from largest starting
      colnames(numb) <- "number"
      numb$cum <- cumsum(numb$number) 
      numb$lngth <- final2[rev(order(final2$lngth)),"lngth"] 
      numb$cumperc <- round(numb$cum/sum(numb$number),5)  
      numb$num5 <- 0
      numb[numb$cumperc<=0.05,"num5"] <- numb[numb$cumperc<=0.05,"number"]
      numb[max(which(numb$cumperc<=0.05))+1,"num5"] <- (0.05-numb[max(which(numb$cumperc<=0.05)),"cumperc"])*sum(numb$number)
      Ind[j,"Lmax5"] <- sum(numb$num5*numb$lngth)/sum(numb$num5)
      
      # indicators
      Ind[j, "L75"] <- min(final2[which(final2$cumsum_perc >= 0.75), "lngth"])
      Ind[j, "L25"] <- min(final2[which(final2$cumsum_perc >= 0.25), "lngth"])
      Ind[j, "Lmed"] <- min(final2[which(final2$cumsum_perc >= 0.5), "lngth"])
      Ind[j, "L95"] <- min(final2[which(final2$cumsum_perc >= 0.95), "lngth"])
      Ind[j, "L90"] <- min(final2[which(final2$cumsum_perc >= 0.90), "lngth"])
      
      final3 <- final2[final2$lngth >= Ind[j, "Lc"], ]    # calculate mean of individuals above Lc
      Ind[j, "Lmean"] <- sum(final3$lngth * final3$number) / sum(final3$number)
      
      final2$biomass <- final2$number * weight[, jj]
      Ind[j, "Lmaxy"] <- final2[final2$biomass == max(final2$biomass), "lngth"]  # length class with max yield
      
      Lopt <- (2 / 3) * Linf[s]
      
      Ind[j, "Pmega"] <- sum(final2[which(final2$lngth >= (Lopt + 0.1 * Lopt)),
                                    "number"]) / sum(final2$number)   # proportion larger Lopt+10%
      Ind[j, "Year"] <- Year[j]
      Ind[j, "Pmegaref"] <- 0.3   # proxy reference point of 30% in catch
      Ind[j, "LFeM"] <- 0.75 * Ind[j, "Lc"] + 0.25 * Ind[j, "Linf"]
    }
    
    #calculate various ratios
    Ind$Lmaxy_Lopt <- Ind$Lmaxy / Ind$Lopt
    Ind$L95_Linf <- Ind$L95 / Ind$Linf
    Ind$Lmean_LFeM <- Ind$Lmean / Ind$LFeM
    Ind$Lmean_Lmat <- Ind$Lmean / Ind$Lmat
    Ind$Lmean_Lopt <- Ind$Lmean / Ind$Lopt
    Ind$Lmax5_Linf <- Ind$Lmax5 / Ind$Linf
    Ind$Lc_Lmat <- Ind$Lc / Ind$Lmat
    Ind$L25_Lmat <- Ind$L25 / Ind$Lmat
    
    if(sex == "M") Males <- Ind
    if(sex == "F") Females <- Ind
    if(sex == "N") Unsexed <- Ind
    
    # write.csv(Ind, file = paste(stock,
    #                             "_",
    #                             sex,
    #                             "_IndicatorRatios_table.csv",
    #                             sep = ""),
    #           row.names = F)
    
  }
  


###################################################################################################
## step 3 plot indicator time series per sex
  file_path <- "~/ICES_MSY/output/"
  
  
for(s in 1:length(S)){
    sex <- S[s] 
    if(sex == "M") Ind <- Males  
    if(sex == "F") Ind <- Females
    if(sex == "N") Ind <- Unsexed
    
png(paste(file_path,
          stock,
          "_",
          sex,
          "_timeseries.png",
          sep = ""),
    bg = "white",
    pointsize = 5,
    units = "cm",
    width = 10,
    height = 18,
    res = 600)

  par(mar = c(5, 4, 3, 4), 
      mfrow = c(3, 1), 
      family = "serif",
      cex = 1.5)
  
    plot(Linf ~ Year,
         data = Ind,
         ylab = "Length",
         col = "transparent",
         main = "(a) Conservation",
         xlab = "Year",
         xlim = c(Year[1],
                  tail(Year, 1) + 1), ylim = c(min(Ind$Lc) * .9,
                                               unique(Ind$Linf) * 1.1),
         bty = "l")
      axis(1, at = Ind$Year, 
           labels = FALSE, 
           cex.axis = 0.1,
           tick = TRUE)
      lines(L95 ~ Year,
            data = Ind, 
            lwd = 2,
            col = "purple")
      text(tail(Year, 1) + 1,
           tail(Ind$L95, 1),
           expression(L["95%"]),
           col = "purple",
           cex = 1.1)
      lines(Lmax5 ~ Year, 
            data = Ind,
            lwd = 2,
            col = "black")
      text(tail(Year, 1) + 1, 
           tail(Ind$Lmax5, 1),
           expression(L["max5%"]),
           col = "black",
           cex = 0.9)
      lines(Lmat ~ Year, 
            data = Ind, 
            lwd = 1, 
            col = "black",
            lty = "dashed")
      text(tail(Year, 1) + 1,
           tail(Ind$Lmat, 1),
           expression(L["mat"]),
           col = "black",
           cex = 1.1)
      lines(Lc~Year, data=Ind, lwd=2, col="blue")
      text(tail(Year,1)+1, tail(Ind$Lc,1), expression(L["c"]), col="blue", cex=1.1)
      lines(Linf~Year, data=Ind, lwd=1, col="black", lty="dashed")
      text(tail(Year,1)+1, tail(Ind$Linf,1), expression(L["inf"]), col="black", cex=1.1)
      lines(L25~Year, data=Ind, lwd=1, col="red")
      text(tail(Year,1)+1, tail(Ind$L25,1), expression(L["25%"]), col="red", cex=1.1)
  
    plot(Linf~Year, data=Ind, ylab="Length", main="(b) Optimal Yield", col="transparent", xlab="Year", xlim=c(Year[1],tail(Year,1)+1), ylim=c(min(Ind$Lc)*.9, unique(Ind$Linf)*1.1), bty="l")
      axis(1, at=Ind$Year, labels=FALSE, cex.axis=0.1, tick=TRUE)
      lines(L75~Year, data=Ind, lwd=1, col="red")
      text(tail(Year,1)+1, tail(Ind$L75,1), expression(L["75%"]), col="red", cex=1.1)
      lines(Lmean~Year, data=Ind, lwd=2, col="darkred")
      text(tail(Year,1)+1, tail(Ind$Lmean,1), expression(L["mean"]), col="darkred", cex=1.1)
      lines(Lopt~Year, data=Ind, lwd=1, col="black", lty="dashed")
      text(tail(Year,1)+1, tail(Ind$Lopt,1), expression(L["opt"]), col="black", cex=1.2)
      lines(Lmaxy~Year, data=Ind, lwd=2, col="green")
      text(tail(Year,1)+1, tail(Ind$Lmaxy,1), expression(L["maxy"]), col="green", cex=1.2)
      lines(Lmat~Year, data=Ind, lwd=1, col="black", lty="dashed")
      text(tail(Year,1)+1, tail(Ind$Lmat,1), expression(L["mat"]), col="black", cex=1.1)
      lines(L25~Year, data=Ind, lwd=1, col="red")
      text(tail(Year,1)+1, tail(Ind$L25,1), expression(L["25%"]), col="red", cex=1.1)
    
    plot(Lmat~Year, data=Ind, type="l", ylab="Length", main="(c) Maximum Sustainable Yield", col="black", lty="dashed", xlab="Year", xlim=c(Year[1],tail(Year,1)+1), ylim=c(min(Ind$Lc)*.9, unique(Ind$Linf)*1.1), bty="l")
      axis(1, at=Ind$Year, labels=FALSE, cex.axis=0.1, tick=TRUE)
      text(tail(Year,1)+1, tail(Ind$Lmat,1), expression(L["mat"]), col="black", cex=1.2)
      lines(Lmean~Year, data=Ind, lwd=2, col="darkred")
      text(tail(Year,1)+1, tail(Ind$Lmean,1), expression(L["mean"]), col="darkred", cex=1.1)
      lines(LFeM~Year, data=Ind, lwd=2, col="blue", lty="dashed")
      text(tail(Year,1)+1, tail(Ind$LFeM,1), expression(L["F=M"]), col="blue", cex=1.2, lty="dashed")

dev.off()

png(paste(file_path,
          stock,
          "_",
          sex,
          "_timeseries_ratios.png",
          sep=""), 
    bg="white", 
    pointsize=5,
    units="cm",
    width=10, 
    height=18, 
    res = 600)

  par( mar = c(5, 4, 3, 4), mfrow=c(3,1), family="serif", cex=1.5)
  
    plot(c(Year[1], tail(Year,1)+3), c(0, 1.5), ylab="Indicator Ratio", col="transparent", main="(a) Conservation", xlab="Year", xlim=c(Year[1], tail(Year,1)+3), ylim=c(0,1.5), bty="l")
      axis(1, at=Ind$Year, labels=FALSE, cex.axis=0.1, tick=TRUE)
      lines(Lmax5_Linf~Year, data=Ind, lwd=2, col="black")
      text(tail(Year,1)+2, tail(Ind$Lmax5_Linf,1), expression(L["max5%"]/L["inf"]), col="black", cex=1.0)
      lines(L95_Linf~Year, data=Ind, lwd=2, col="purple")
      text(tail(Year,1)+2, tail(Ind$L95_Linf,1), expression(L["95%"]/L["inf"]), col="purple", cex=1.1)
      lines(Pmega~Year, data=Ind, lwd=2, col="blue")
      text(tail(Year,1)+2, tail(Ind$Pmega,1), expression(P["mega"]), col="blue", cex=1.2)
      lines(Pmegaref~Year, data=Ind, lwd=1, col="black", lty="dashed")
      text(tail(Year,1)+2, tail(Ind$Pmegaref,1), expression("30%"), col="black", cex=1.0)
      lines(Lc_Lmat~Year, data=Ind, lwd=2, col="red")
      text(tail(Year,1)+2, tail(Ind$Lc_Lmat,1), expression(L["c"]/L["mat"]), col="red", cex=1.1)
      lines(L25_Lmat~Year, data=Ind, lwd=2, col="darkred")
      text(tail(Year,1)+2, tail(Ind$L25_Lmat,1), expression(L["25"]/L["mat"]), col="darkred", cex=1.1)
    
    plot(c(Year[1], tail(Year,1)+3), c(0, 1.6), ylab="Indicator Ratio", col="transparent", main="(b) Optimal yield", xlab="Year" ,xlim=c(Year[1], tail(Year,1)+3), ylim=c(0,1.6), bty="l")
      axis(1, at=Ind$Year, labels=FALSE, cex.axis=0.1, tick=TRUE)
      lines(Lmean_Lopt~Year, data=Ind, lwd=2, col="darkred")
      text(tail(Year,1)+2, tail(Ind$Lmean_Lopt,1), expression(L["mean"]/L["opt"]), col="darkred", cex=1.1)
      lines(Lmaxy_Lopt~Year, data=Ind, lwd=2, col="green")
      text(tail(Year,1)+2, tail(Ind$Lmaxy_Lopt,1), expression(L["maxy"]/L["opt"]), col="green", cex=1.1)
    
    plot(c(Year[1], tail(Year,1)+3), c(0, 1.6), ylab="Indicator Ratio", col="transparent", main="(c) Maximum sustainable yield", xlab="Year", xlim=c(Year[1], tail(Year,1)+3), ylim=c(0,1.6), bty="l")
      axis(1, at=Ind$Year, labels=FALSE, cex.axis=0.1, tick=TRUE)
      lines(Lmean_LFeM~Year, data=Ind, lwd=2, col="blue")
      text(tail(Year,1)+2, tail(Ind$Lmean_LFeM,1), expression(L["mean"]/L["F=M"]), col="blue",cex=1.1)

dev.off()

}


