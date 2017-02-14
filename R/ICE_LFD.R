#####################################################################
# Calculate length-frequency distributions from InterCatch Exchange #
#####################################################################
# Adapted from T. Miethe and C. Silva, WKLIFE-V, Oct2015

ICE_LFD <- function(file_path = "",
                   startyear = 2000,
                   endyear = 2014,
                   country = "XYZ",
                   stock = "nep",
                   ctype = "L",  # for landings L, discard D
                   ftype = ".dat",
                   S = c("M","F")) {
  
  if(!file.exists(paste(file_path,
                        country, 
                        "_",
                        startyear,
                        "_DataCall2015_DLS_",
                        stock,
                        ftype,
                        sep = ""))) {
    stop(paste0("Check that the file name and path are correct. \n Current file.path: ", file_path, 
                " \n Required file name structure: /XYZ_2000_DataCall2015_DLS_nep.dat"))
  }
  
  # Load first year data, Intercatch SD file
  icdata <- read.table(paste(file_path,
                             country, 
                             "_",
                             startyear,
                             "_DataCall2015_DLS_",
                             stock,
                             ftype,
                             sep = ""),
                       sep = ",",
                       header = FALSE, 
                       col.names = as.character(1:33),
                       fill = T)
  
  colnames(icdata) <- c("RecordType", "Country", "Year", "SeasonType",
                        "Season", "Fleet", "AreaType", "FishingArea",
                        "DepthRange", "Species", "Stock", "CatchCategory",
                        "ReportingCategory", "Sex", "CANUMtype", 
                        "AgeLength", "PlusGroup", "sampledCatch", 
                        "NumSamplesLngt", "NumLngtMeas", "NumSamplesAge",
                        "NumAgeMeas", "unitMeanWeight","unitCanum",
                        "UnitAgeOrLength","UnitMeanLength", "Maturity",
                        "NumberCaught", "MeanWeight", "MeanLength", 
                        "varNumLanded", "varWgtLanded", "varLgtLanded")
  
  icdata <- subset(icdata, RecordType == 'SD' & CatchCategory == ctype)
  
  if(icdata$unitMeanWeight[1] == "g") {
    multiplier <- 0.001  
  } else {
    multiplier <- 1   # adapt MeanWeight multiplier
  }
  
  m <- icdata[icdata$Sex == "M", c("MeanLength", "NumberCaught")]
  f <- icdata[icdata$Sex == "F", c("MeanLength", "NumberCaught")]
  ns <- icdata[icdata$Sex == "N", c("MeanLength", "NumberCaught")]
  
  mw <- icdata[icdata$Sex == "M", c("MeanLength", "MeanWeight")]
  fw <- icdata[icdata$Sex == "F", c("MeanLength", "MeanWeight")]
  nsw <- icdata[icdata$Sex == "N", c("MeanLength", "MeanWeight")]
  
  colnames(mw) <- c("MeanLength", startyear)
  colnames(m) <- c("MeanLength", startyear)
  colnames(fw) <- c("MeanLength", startyear)
  colnames(f) <- c("MeanLength", startyear)
  colnames(nsw) <- c("MeanLength", startyear)
  colnames(ns) <- c("MeanLength", startyear)
  
  
  Year <- c(startyear:endyear)
  Year1 <- c((startyear+1):endyear)
  
  #load Intercatch format data and merge years
  for(year in Year1){
    for(i in 1:length(S)){
      sex <- S[i]  
      
      # select data file name for Intercatch SD files
      icdata <- read.table(paste(file_path, country, "_",
                                 year, "_DataCall2015_DLS_",
                                 stock, ftype, 
                                 sep = ""),
                           sep = ",",
                           header = FALSE, 
                           col.names = as.character(1:33),
                           fill = T)
      colnames(icdata) <- c("RecordType", "Country", "Year", "SeasonType",
                            "Season", "Fleet", "AreaType", "FishingArea",
                            "DepthRange", "Species", "Stock", "CatchCategory",
                            "ReportingCategory", "Sex", "CANUMtype", "AgeLength", 
                            "PlusGroup", "sampledCatch", "NumSamplesLngt", "NumLngtMeas",
                            "NumSamplesAge", "NumAgeMeas", "unitMeanWeight", "unitCanum",
                            "UnitAgeOrLength", "UnitMeanLength", "Maturity", "NumberCaught",
                            "MeanWeight", "MeanLength", "varNumLanded",
                            "varWgtLanded", "varLgtLanded")
      icdata <- subset(icdata, RecordType == 'SD' & CatchCategory == ctype)
      
      icdatanum <- icdata[icdata$Sex == sex,
                          c("MeanLength", "NumberCaught")]
      colnames(icdatanum) <- c("MeanLength", year)
      
      icdataw <- icdata[icdata$Sex == sex,
                        c("MeanLength", "MeanWeight")]
      colnames(icdataw) <- c("MeanLength", year)
      
      if(sex=="M") m <- merge(m, icdatanum, by="MeanLength", all=TRUE)
      if(sex=="F") f <- merge(f, icdatanum, by="MeanLength", all=TRUE)
      if(sex=="N") ns <- merge(ns, icdatanum, by="MeanLength", all=TRUE)
      
      if(sex=="M") mw <- merge(mw, icdataw, by="MeanLength", all=TRUE)
      if(sex=="F") fw <- merge(fw, icdataw, by="MeanLength", all=TRUE)
      if(sex=="N") nsw <- merge(nsw, icdataw, by="MeanLength", all=TRUE)
      
    }
  }
  
  ns <- ns[ns$MeanLength > 0,]
  m <- m[m$MeanLength > 0,]
  f <- f[f$MeanLength > 0,]
  
  mw <- mw[mw$MeanLength > 0,]
  fw <- fw[fw$MeanLength > 0,]
  nsw <- nsw[nsw$MeanLength > 0,]
  
  # remove NAs
  fw[is.na(fw)] <- 0
  mw[is.na(mw)] <- 0
  nsw[is.na(nsw)] <- 0
  
  f[is.na(f)] <- 0
  m[is.na(m)] <- 0
  ns[is.na(ns)] <- 0
  
  fw <- fw*multiplier
  mw <- mw*multiplier
  nsw <- nsw*multiplier
  
  out_list <- list(male_weight = mw,
                   female_weight = fw,
                   unsexed_weight = nsw,
                   male_length = m,
                   female_length = f,
                   unsexed_length = ns)
  return(out_list)
}

