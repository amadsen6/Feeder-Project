.libPaths("C:/Users/Laura/Desktop/RLibrary")
library(tidyverse)
library(stringr)
library(igraph)
library(asnipe)
library(lubridate)
library(assortnet)

#Feeder RFID data for all species unsorted
datnut = read.csv("CleanFeederData3_10.csv")


datnut$Lubridate = parse_date_time(datnut$Datetime, "%m%d%H%M%S") #make date time object
datnut$Timestamp=as.numeric(datnut$Lubridate) #get seconds passed from date time object, neet to do this from start time of project
datnut$Day=day(datnut$Lubridate) #pull day from datetime object
datnut$Month=month(datnut$Lubridate)
datnut = unite(datnut, "LoggerDate", c("Logger", "Month", "Day"), remove=F) #make logger day, month object
datnut=datnut[complete.cases(datnut),]


##filter out logger 04 date collect 02012019

datno4 = datnut[!(datnut$Logger=="LOGGER04"& datnut$DateCollected =="2019-02-01"),]

#Need to sort by species
birddat = read.csv("RFID_Records_fixed.csv", colClasses = c("RFID" ="character"))

dat = merge(datno4, birddat[, c("RFID", "Species")], by="RFID")

datDOWO = filter(dat, Species == "DOWO")
datDOWO$X = NULL
datDOWO = distinct(datDOWO)

datWBNU = filter(dat, Species == "WBNU")
datWBNU$X = NULL
datWBNU = distinct(datWBNU)


#testing issues
WBNUtest = filter(datWBNU, LoggerDate == "LOGGER02_1_29")
test = gmmevents(WBNUtest$Timestamp, WBNUtest$RFID, WBNUtest$LoggerDate)

WBNUtest = filter(datWBNU, Logger != "LOGGER02")



#running events by species

DOWOflocks = gmmevents(datDOWO$Timestamp, datDOWO$RFID, datDOWO$LoggerDate) 
saveRDS(DOWOflocks, "conspecificDOWOflocks.rds")


WBNUflocks = gmmevents(datWBNU$Timestamp, datWBNU$RFID, datWBNU$LoggerDate) 
saveRDS(WBNUflocks, "conspecificWBNUflocks.rds")
