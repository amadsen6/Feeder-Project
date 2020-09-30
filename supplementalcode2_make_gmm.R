library(tidyverse)
library(stringr)
#library(igraph)
library(asnipe)
library(lubridate)
#library(assortnet)

#Feeder RFID data for all species unsorted
datnut = read.csv("CleanFeederData3_10.csv")


datnut$Lubridate = parse_date_time(datnut$Datetime, "%m%d%H%M%S") #make date time object
datnut$Timestamp=as.numeric(datnut$Lubridate) #get seconds passed from date time object, neet to do this from start time of project
datnut$Day=day(datnut$Lubridate) #pull day from datetime object
datnut$Month=month(datnut$Lubridate)
datnut = unite(datnut, "LoggerDate", c("Logger", "Month", "Day"), remove=F) #make logger day, month object
datnut=datnut[complete.cases(datnut),]


##filter out logger 04 date collect 02012019

# datno4 = datnut[!(datnut$Logger=="LOGGER04"& datnut$DateCollected =="2019-02-01"),]
# 
# #Need to sort by species
# birddat = read.csv("RFID_Records_fixed.csv", colClasses = c("RFID" ="character"))
# 
# dat = merge(datno4, birddat[, c("RFID", "Species")], by="RFID")


##filter out logger 04 date collect 02012019

#Need to sort by species
birddat = read.csv("RFID_Records_fixed.csv", colClasses = c("RFID" ="character"))

dat = merge(datnut, birddat[, c("RFID", "Species")], by="RFID")


datDOWO = filter(dat, Species == "DOWO")
datDOWO = distinct(datDOWO)

datWBNU = filter(dat, Species == "WBNU")
datWBNU = datWBNU[,which(colnames(datWBNU)!="X")]
datWBNU = distinct(datWBNU)

#testing issues
#WBNUtest = filter(datWBNU, LoggerDate == "LOGGER02_1_29")
test = gmmevents(WBNUtest$Timestamp, WBNUtest$RFID, WBNUtest$LoggerDate)

#WBNUtest2 = filter(datWBNU, LoggerDate == "LOGGER01_1_29")
test2 = gmmevents(WBNUtest2$Timestamp, WBNUtest2$RFID, WBNUtest2$LoggerDate)

#DOWOtest = filter(datDOWO, LoggerDate == "LOGGER02_1_29")
test.dowo = gmmevents(DOWOtest$Timestamp, DOWOtest$RFID, DOWOtest$LoggerDate)
#running events by species

DOWOflocks_replicate4 = gmmevents(datDOWO$Timestamp, datDOWO$RFID, datDOWO$LoggerDate) 

WBNUflocks = gmmevents(datWBNU$Timestamp, datWBNU$RFID, datWBNU$LoggerDate) #stopeed after logger 2 1-29 replacement has length zero


###
gmmDOWO=readRDS("conspecificDOWOflocks.rds")
#load("DOWOflocks_replicate")
gbi1=gmmDOWO$gbi
gbi2=DOWOflocks$gbi
gbi3=DOWOflocks_replicate$gbi
gbi5=DOWOflocks_replicate2$gbi
gbi6=DOWOflocks_replicate4$gbi
dim(gbi1)
dim(gbi2)
dim(gbi3)
dim(gbi5)
dim(gbi6)
dim(DOWOflocks_replicate2$gbi)
###
load('all_visits.dat')
all_visits <- all_visits %>%
  ungroup() %>%
  mutate(index = as.numeric(Date- min(Date)+1)) %>%
  mutate(Timestamp = as.numeric(Datetime))

dat.dowo <- all_visits %>%
  ungroup() %>%
  filter(Species == "DOWO") %>%
  mutate(LoggerDate = paste0(Logger, Date))
dat.wbnu <- all_visits %>%
  ungroup() %>%
  filter(Species == "WBNU") %>%
  mutate(LoggerDate = paste0(Logger, Date))

dowoflocks_rep2=gmmevents(dat.dowo$Timestamp, dat.dowo$RFID, dat.dowo$LoggerDate)
gbi44=dowoflocks_rep2$gbi
dim(gbi4)

dowo_net=graph_from_adjacency_matrix(get_network(gbi4), mode="undirected", weighted=T)
indivs=read.csv("RFID_Records_fixed.csv")
indivs$RFID=as.character(indivs$RFID)
head(indivs)

V(dowo_net)$sex=indivs[match(V(dowo_net)$name, indivs$RFID),"Sex"]
sex_color=data.frame(sex=c("F", "M", "U"), color=c("yellow", "purple", "white"))
plot(dowo_net, vertex.color=sex_color[match(V(dowo_net)$sex, sex_color$sex), "color"], vertex.label="", edge.width=E(dowo_net)$weight*30)

sexassort_dowo=assortment.discrete(as_adj(dowo_net, sparse=F, attr="weight"), V(dowo_net)$sex, SE=T)
sexassort_dowo$r