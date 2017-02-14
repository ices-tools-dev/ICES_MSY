rm(list = ls())
# Download and unzip data from InterCatch into a /data folder. The code below will extract the length data and
# aggregate into a length frequency distribution with corresponding mean weight at length.
# 
library(data.table)

# Unzip each of the InterCatch download files to a folder "/data" in your working directory
wd <- getwd()
file_list <- paste0(wd,
                    "/data/",
                    grep("CatchAndSampleDataTables.txt", 
                         list.files("data/", recursive = TRUE),
                         value = TRUE))

 
dat <- data.frame()
for (fl in file_list){
  get_text <- scan(fl,
               what = 'character',
               sep = '\t')
  
  table2 <- get_text[(which(get_text == "TABLE 2.") + 3) : length(get_text)]
  tmp <- table2[-c(1:56)]			  
  
  dat_fl <- data.frame(matrix(tmp,
                                  ncol = 27,
                                  byrow = T), 
                           stringsAsFactors =F)
  colnames(dat_fl) <- table2[1:27]
  dat_fl <- data.table(dat_fl)
  # dat_fl <- dat_fl[, CATON := as.numeric(as.character(CATON))]
  dat_fl <- dat_fl[, CANUM := as.numeric(as.character(CANUM))]
  dat_fl <- dat_fl[, WECA := as.numeric(as.character(WECA))]
  dat_fl <- dat_fl[, AgeOrLength := as.numeric(as.character(AgeOrLength))]
  dat_fl <- dat_fl[, Year := as.numeric(Year)]
  # dat_fl <- dat_fl[, Area := as.factor(Area)]
  # dat_fl <- dat_fl[, Fleet := factor(Fleet)]
  # dat_fl <- dat_fl[, Season := factor(Season)]
  # dat_fl <- dat_fl[, Country := factor(Country)]

  dat <- rbind(dat, dat_fl)
}


DF <- dat[, list(CANUM = sum(CANUM),
                 WECA = weighted.mean(WECA, CANUM)), 
           by = c('Year','Stock','AgeOrLength')]

DF <- as.data.frame(DF)

# First column as length remaining as numbers at length
df_number <- as.data.frame(DF[,c("Year", "AgeOrLength", "CANUM")])
col_name <- c("Length", paste0("X", sort(unique(df_number$Year))))
df_number <- reshape(df_number, idvar = c('AgeOrLength'), 
                     timevar = "Year", 
                     direction = 'wide')
colnames(df_number) <- col_name
df_number <- df_number[order(df_number[, "Length"]),]


# First column as length remaining are mean weight at length
df_weight <- as.data.frame(DF[,c("Year", "AgeOrLength", "WECA")])
col_name <- c("MeanLength", paste0("X", sort(unique(df_weight$Year))))
df_weight <- reshape(df_weight, idvar = c('AgeOrLength'), 
                     timevar = "Year", 
                     direction = 'wide')
colnames(df_weight) <- col_name
df_weight <- df_weight[order(df_weight[, "MeanLength"]),]
