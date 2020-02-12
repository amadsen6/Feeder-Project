### generating exploratory tendency metrics for Kris Hans UCARE project

require(tidyverse)
require(lubridate)
##data
load("all_visits.dat")
head(all_visits)

all_visits$Date
wbnu=all_visits %>% filter(Species=="WBNU") %>% arrange(RFID, Logger) 

wbnu.visits.daily=wbnu%>% group_by(RFID, Date) %>% count()
wbnu.visits.daily.feeders=wbnu%>% group_by(RFID, Date) %>% n_distinct(Logger)
