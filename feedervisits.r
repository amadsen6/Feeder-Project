library(tidyverse)
library(lme4)
library(lmerTest)

dat=read.csv("morn_visits_to_publish.csv")

head(dat)

dowo=dat %>% filter(Species=="DOWO") 

mod=lm(sumvisits~nightlows, data=dowo)
summary(mod)

mod=lmer(sumvisits~nightlows + (1|RFID), data=dowo)
anova(mod)

dowo_summary=dat %>% filter(Species=="DOWO") %>% group_by(Date) %>% summarise(mean_visit=mean(sumvisits), temp=mean(nightlows))


ggplot(dowo_summary, aes(x=temp, y=mean_visit))+
  geom_point() +
  geom_smooth(method="lm")

