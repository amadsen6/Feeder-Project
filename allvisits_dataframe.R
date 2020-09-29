## Data cleaning 01/28/2019
## Annie Madsen
##packages
.libPaths("C:/Users/Shizuka Lab/Documents/R/win-library/3.4")
require(tidyverse)
require(data.table)
setwd("C:/Users/Shizuka Lab/Documents/RFID_data")

#Read files
##Read all csv files in the working directory, list as string in filenames
filenames <- list.files(recursive = TRUE)

##Function to add path names as a variable
read_csv_filename <- function(filename){
  ret <- read.csv(filename, header = FALSE, col.names = c("V1", "V2", "V3"), colClasses = c("character","character","character"))
  ret$Source <- filename
  ret
}

RFIDdata <- as.data.frame(rbindlist(lapply(filenames, read_csv_filename), fill = TRUE))


RFIDdata <- RFIDdata %>%
  separate(Source, into = c("Folder", "Logger"), sep = "/") %>%
  select(-Folder)
RFIDdata$Logger <- gsub('.CSV', '', RFIDdata$Logger)##remove filetype text
RFIDdata$Logger <- gsub('.csv', '', RFIDdata$Logger)
colnames(RFIDdata)[c(1:3)] <- c("Datetime", "AntennaRFID", "Nutty_error")## rename columns
RFIDdata <- RFIDdata[,c(4,1,2,3)]##rearrange columns

## Fix antenna error
RFIDdata$AntennaRFID <- as.character(RFIDdata$AntennaRFID)

RFIDdata$AntennaRFID <- ifelse(nchar(RFIDdata$AntennaRFID) == 11, paste0("0", RFIDdata$AntennaRFID), paste0(RFIDdata$AntennaRFID))

########### 
## Fix Nutty_error
### remove blank rows
RFIDdata <- RFIDdata %>%
  filter(nchar(AntennaRFID) > 11)

### separate datetime and RFID in rows with nutty_error
errors <- RFIDdata %>%
  filter(nchar(AntennaRFID) > 12) %>%
  separate(AntennaRFID, into = c("AntennaRFID2", "Datetime2"), sep = 12)
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
normal <- RFIDdata %>%
  filter(nchar(AntennaRFID) == 12) %>%
  select(-Nutty_error)
## bind them all together
fixedrows <- rbind(fixed, nutty, normal)
## fix antenna error
fixedrows$AntennaRFID <- ifelse(nchar(fixedrows$AntennaRFID) == 11, paste0("0", fixedrows$AntennaRFID), paste0(fixedrows$AntennaRFID))

##########
require(lubridate)
RFIDdf <- fixedrows %>%
  separate(AntennaRFID, into = c("Antenna", "RFID"), sep = 2)

RFIDdf$Datetime <- ifelse(nchar(RFIDdf$Datetime) == 9, paste0("0", RFIDdf$Datetime), paste0(RFIDdf$Datetime))
RFIDdf$Datetime <- as.POSIXct(RFIDdf$Datetime, format = "%m%d%H%M%S")

RFIDdf$RFID <- as.character(RFIDdf$RFID)

demo <- read.csv("C:/Users/Shizuka Lab/Documents/Madsen/Pioneers Park Project/Bird Banding Data/RFID_Records_03052019.csv")

RFIDdem <- RFIDdf %>%
  left_join(demo, by = "RFID") %>%
  select(-Date, -Time) %>%
  filter(Species != "TestTag") %>%
  mutate(Date = date(Datetime))

## Look at average length of time between hits
datasort <- arrange(RFIDdem %>% group_by(Logger, RFID), Datetime, .by_group = TRUE)
datasort$values = 0
datasort <- datasort %>% 
  select(Logger, Datetime, Antenna, RFID, Species, Weight, Date, values)

## calculate time difference between consecutive hits for each individual
for(i in 2:nrow(datasort)){
  
  datasort$values[i] <- difftime(time1 = datasort$Datetime[i], time2 = datasort$Datetime[i-1], units = "secs")
  
}

## clean it up
datasort$values[datasort$values < 0] <- "0"
datasort$values <- as.numeric(datasort$values)

datasort <- datasort %>%
  mutate(minutes = as.numeric(values)/60)

## collapse consecutive hits
require(lubridate)
visits <- datasort %>%
  mutate(consec = ifelse(values <= 2, paste0("0"), paste0("1")))

dif_visits <- visits %>%
  filter(consec == "1")

same_visits <- visits %>%
  filter(consec == "0") %>%
  mutate(seqrle = sequence(rle(consec)$lengths)) %>%
  filter(seqrle == 1)

all_visits <- rbind(dif_visits, same_visits) %>%
  select(-seqrle, -consec) %>%
  mutate(Datehour = floor_date(Datetime, unit = "hours"))