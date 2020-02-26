### generating exploratory tendency metrics for Kris Hans UCARE project

require(tidyverse)
require(lubridate)
##data
load("all_visits.dat")
head(all_visits)

all_visits$Date

usedat=all_visits %>% filter(Species=="WBNU"|Species=="DOWO") %>% arrange(RFID, Logger) 

visits.daily=usedat%>% group_by(RFID, Date) %>% count()

##visits.daily is the total number of visits to any feeder for a given bird on a given day


visits.daily.feeders=usedat%>% group_by(RFID, Date) %>% select(RFID, Date, Logger) %>% summarise(n_distinct(Logger))
##wbnu.visits.daily.feeders is the number of distinct feeders used by each bird on a given day

logger_array=usedat%>% group_by(RFID, Date) %>% select(Logger,  Date,RFID) %>% table()

apply(logger_array, 3, function(x) {
  if (colSums(x)==0) p=NA else p=x/colSums(x)
  p
  })

  #ln_p=log(p)
  #sum(p*ln_p, na.rm=T)}

indiv_date_sum=apply(logger_array, c(2,3), sum)

indiv_avg_visit=colMeans(indiv_date_sum)

#write.csv(indiv_date_sum, "visits_indiv_x_date.csv")


indiv_date_feeders=pivot_wider(visits.daily.feeders, id_cols=Date, names_from = RFID, values_from = `n_distinct(Logger)`)

#write.csv(indiv_date_feeders, "feeders_indivxdate.csv")

#######
head(usedat)

usedat%>% select(Datetime, RFID, Logger) %>% group_by(RFID)
