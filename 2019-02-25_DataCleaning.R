
require(tidyverse)
require(data.table)
library(stringr)



#Read files

filenames <- list.files(recursive = TRUE, pattern=".CSV")

##Function to add path names as a variable
read_csv_filename <- function(filename){
  ret <- read.csv(filename, header = FALSE, col.names=c("V1", "V2", "V3"), colClasses = c("character", "character", "character"))
  ret$Source <- filename
  ret
}



RFIDdata <- as.data.frame(rbindlist(lapply(filenames, read_csv_filename), fill=T))

##Use this when doing multiple folders
RFIDdataClean <-RFIDdata %>%
  select("V1", "V2", "Source", "V3") %>% 
  separate(Source, into = c("Folder", "Logger"), sep = "/") %>% 
  setnames(c("Datetime", "AntennaRFID", "DateCollected", "Logger", "Nutty_error")) 


#making things prettier
RFIDdataClean$Logger = gsub('.CSV', '', RFIDdataClean$Logger) #remove csv from logger column
RFIDdataClean$DateCollected =gsub("_FeederData", "", RFIDdataClean$DateCollected) #shorten folder name to just date collected 

########### 
## Fix Nutty_error
# RFIDdata %>% if(nchar(RFIDdata$AntennaRFID) > 12){
#   separate(RFIDdata$AntennaRFID, into = c("AntennaRFID1", "Datetime2"), sep = 12)
# }else{print("NA")}
### remove blank rows
RFIDdataClean$AntennaRFID <-
  as.character(RFIDdataClean$AntennaRFID)


##SOMETHING IS WRONG HERE AND I CAN'T FIGURE OUT WHAT, WHERE ARE ALL THE NUTTY ERROR RFID TAGS?

### separate datetime and RFID in rows with nutty_error

errors <- RFIDdataClean %>%
  dplyr::filter(nchar(AntennaRFID) > 12) %>%
  tidyr::separate(AntennaRFID, into = c("AntennaRFID2", "Datetime2"), sep = 12)
## pull one column pair into a separate df
fixed <- errors %>%
  select(-Datetime2, -Nutty_error)
colnames(fixed)[colnames(fixed) == "AntennaRFID2"] <- "AntennaRFID"##fix colnames to match
## pull the other column pair into separate df
nutty <- errors %>%
  select(-Datetime, -AntennaRFID2)
colnames(nutty)[colnames(nutty) == "Nutty_error"] <- "AntennaRFID"##fix colnames to match
colnames(nutty)[colnames(nutty) == "Datetime2"] <- "Datetime"
## grab correct values from original df
normal <- RFIDdataClean %>%
  filter(nchar(AntennaRFID) == 12) %>%
  select(-Nutty_error)
## bind them all together
fixedrows <- rbind(fixed, nutty, normal)
##########
## Split 
RFIDdf <- fixedrows %>%
  separate(AntennaRFID, into = c("Antenna", "RFID"), sep = 2)

##trying to figure out weird RFID tags - looks like majority are coming from Logger04 on 2019-02-01 where the antenna code is messed (either 10 or 00)

filter(RFIDdf, Antenna=="10")
filter(RFIDdf, Antenna=="00")


#Data through 3-10
write.csv(RFIDdf, "CleanFeederData3_10.csv")

