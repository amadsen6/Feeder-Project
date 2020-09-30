library(tidyverse)
library(stringr)
#library(igraph)
library(asnipe)
library(lubridate)
#library(assortnet)

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