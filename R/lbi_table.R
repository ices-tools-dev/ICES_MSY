lbi_table <- function(Ind) {

  refList <- c("Year", "Lc_Lmat", "L25_Lmat" ,"Lmax5_Linf","Pmega",
               "Lmean_Lopt", "Lmean_LFeM")
  cols <- c(refList, names(Ind)[-which(names(Ind) %in% refList)])
  Ind <- Ind[cols]
  
  refTable <- Ind[Year %in% seq(from = max(Year) - 2, 
                                to = max(Year),
                                by = 1),
                  colnames(Ind) %in% refList]
  refTable <- round(refTable, 2)
  
  refTable$Lmean_LFeM[refTable$Lmean_LFeM >= 1] <- paste0("**",
                                                          refTable$Lmean_LFeM[refTable$Lmean_LFeM >= 1],
                                                          "**")
  
  refTable$Lc_Lmat[refTable$Lc_Lmat > 1] <- paste0("**",
                                                   refTable$Lc_Lmat[refTable$Lc_Lmat > 1],
                                                   "**")
  
  refTable$L25_Lmat[refTable$L25_Lmat > 1] <- paste0("**",
                                                     refTable$L25_Lmat[refTable$L25_Lmat > 1],
                                                     "**")
  
  refTable$Lmax5_Linf[refTable$Lmax5_Linf > .8] <- paste0("**",
                                                          refTable$Lmax5_Linf[refTable$Lmax5_Linf > .8],
                                                          "**")
  
  refTable$Pmega[refTable$Pmega > 1] <- paste0("**",
                                               refTable$Pmega[refTable$Pmega > .3],
                                               "**")
  
  refTable$Lmean_Lopt[refTable$Lmean_Lopt > .9] <- paste0("**",
                                                          refTable$Lmean_Lopt[refTable$Lmean_Lopt > .9],
                                                          "**")
  
  
  return(kable(refTable))
  
}
