### Re-doing the gmmevents function...

library(tidyverse)
library(stringr)
library(asnipe)
library(lubridate)


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

#Run Gaussian Mixture Model. Note that this is a stochastic process so end result can be slightly different each time.
dowoflocks=gmmevents(dat.dowo$Timestamp, dat.dowo$RFID, dat.dowo$LoggerDate)
wbnuflocks=gmmevents(dat.wbnu$Timestamp, dat.wbnu$RFID, dat.wbnu$LoggerDate)

saveRDS(dowoflocks, "conspecificDOWOflocks_v2.rds")
saveRDS(wbnuflocks, "conspecificWBNUflocks_v2.rds")