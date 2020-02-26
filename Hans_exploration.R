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
#write.csv(visits.daily, "dailyvisits.csv")

visits.daily.feeders=usedat%>% group_by(RFID, Date) %>% select(RFID, Date, Logger) %>% summarise(n_distinct(Logger))
#write.csv(visits.daily.feeders, "dailyfeeders.csv")
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

wide.dat=usedat%>% select(Date, Datetime, RFID, Logger) %>% group_by(RFID, Date, Logger) %>% count() %>% pivot_wider(values_from = n, names_from = Logger)

shannon.wide.dat= wide.dat  %>% mutate(sum=sum(LOGGER04,LOGGER08,LOGGER07,LOGGER05,LOGGER02,LOGGER06,LOGGER01,LOGGER03, na.rm=T)) %>% filter(sum>10) %>% mutate(H=sum(-LOGGER01/sum*log(LOGGER01/sum), -LOGGER02/sum*log(LOGGER02/sum), -LOGGER03/sum*log(LOGGER03/sum), -LOGGER04/sum*log(LOGGER04/sum), -LOGGER05/sum*log(LOGGER05/sum), -LOGGER06/sum*log(LOGGER06/sum), -LOGGER07/sum*log(LOGGER07/sum), -LOGGER08/sum*log(LOGGER08/sum), na.rm=T))                                                                                              
shannon=shannon.wide.dat %>% select(Date, RFID, H)

boxplot(H~RFID, data=shannon)

write.csv(shannon, "dailyshannon.csv")

daily_dat=inner_join(visits.daily, visits.daily.feeders) %>% inner_join(shannon)

write.csv(daily_dat, "daily_data.csv")
