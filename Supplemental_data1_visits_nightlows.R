#### Analysis code for Madsen, Vander Meiden & Shizuka: Social partners and temperature jointly affect morning foraging activity of small birds in winter#####

require(tidyverse)
require(lubridate)

##data
load("all_visits.dat")
load("weather.dat")

weather$Date <- as.Date(weather$Date, format = "%m/%d/%Y")
weather <- weather %>%
  mutate(Hour = hour(datetime)) %>%
  filter(Hour >= 19 | Hour <= 4) %>%
  mutate(Lagdate = ifelse(Hour <= 4, paste0(lag(Date)), paste0(Date)))  %>%
  group_by(Date) %>%
  summarize(nightlows = min(Temp_C)) %>%
  ungroup() %>%
  mutate(index = as.numeric(Date-min(Date)+1))


demo <- read.csv("RFID_Records_fixed.csv")
demo <- demo %>%
  select(-Date) %>%
  filter(Species != "SCJU")


#dates with no more than one feeder has NA
usedates=as.Date(c("2019-01-29", "2019-01-30", "2019-02-04", "2019-02-15", "2019-02-16", "2019-02-17", "2019-02-18", "2019-02-19", "2019-02-20", "2019-02-22", "2019-02-26", "2019-02-27", "2019-02-28", "2019-03-02", "2019-03-04", "2019-03-05", "2019-03-06", "2019-03-07", "2019-03-08"))

### Species Models for DOWO and WBNU
morn_visits_to_publish <- all_visits %>%
  mutate(Hour = hour(Datetime)) %>%
  mutate(Date = date(Datetime)) %>%
  filter(Hour >= 6 | Hour <= 11) %>%
  group_by(RFID, Date) %>%
  summarise(sumvisits = n()) %>%
  ungroup() %>%
  left_join(weather, by = "Date") %>%
  left_join(demo, by = "RFID") %>%
  select(Date, RFID, Species, Weight, Tarsus, Wing, Sex, Age, Culmen, nightlows, sumvisits) %>%
  filter(Date < "2019-03-11" & Date > "2019-01-25") %>%
  filter(Species == "DOWO" | Species == "WBNU") %>%
  select(Date, RFID, Species, Sex, sumvisits,nightlows)

morn_visits_to_publish$RFID <- as.factor(morn_visits$RFID)

write.csv(morn_visits_to_publish, "Supplementaldata_morningvisits.csv", row.names = F)

