
## This time with single-species networks
library(asnipe)
library(igraph)

#import files
gmmDOWO=readRDS("conspecificDOWOflocks.rds") #import gmm results file. Will need this for permuting group-by-individual matrices.
gmmWBNU=readRDS("conspecificWBNUflocks.rds") #import gmm results file. Will need this for permuting group-by-individual matrices.


gbi_dowo=gmmDOWO$gbi
gbi_wbnu=gmmWBNU$gbi

dowoadj=get_network(gbi_dowo)
wbnuadj=get_network(gbi_wbnu)

diag(dowoadj)=NA
diag(wbnuadj)=NA

###
load("all_visits.dat")
load("weather.dat")
demo=read.csv("RFID_Records_fixed.csv")
all_visits
require(tidyverse)

## daily visits per individual
dv <- all_visits %>%
  filter(Date < "2019-03-11" & Date > "2019-01-25") %>%
  group_by(RFID, Date) %>%
  summarise(dailyvisits = n()) 

## mean daily visits per individual
ref <- dv %>%
  group_by(RFID) %>%
  summarise(mdv = mean(dailyvisits)) %>%
  ungroup()

## sd of daily visits
dvsd <- sd(dv$dailyvisits)

## build the df
require(reshape)
mat_ddm <- cast(dv, Date ~ RFID, value = "dailyvisits")
mat_ddm[is.na(mat_ddm)] <- 0
mat <- mat_ddm[2:46]/dvsd

myvec <- (ref$mdv[match(names(mat_ddm[2:45]), ref$RFID)])/dvsd

mat_final <- mat[1] - myvec[1]
for(i in 2:44){
  mat_temp <- mat[i] - myvec[i]
  mat_final <- cbind(mat_final, mat_temp)
}

mat_final

#weather dat
library(lubridate)
weather$Date <- as.Date(weather$Date, format = "%m/%d/%Y")
weather <- weather %>%
  mutate(Hour = hour(datetime)) %>%
  filter(Hour >= 19 | Hour <= 4) %>%
  mutate(Lagdate = ifelse(Hour <= 4, paste0(lag(Date)), paste0(Date))) %>%
  group_by(Date) %>%
  summarize(nightlows = min(Temp_C)) %>%
  ungroup()
weather

weather.use=weather[which(weather$Date>"2019-01-26"&weather$Date<"2019-03-10"),]
weather.use

##
morn_visits <- all_visits %>%
  mutate(Hour = hour(Datetime)) %>%
  mutate(Date = date(Datetime)) %>%
  filter(Hour >= 6 | Hour <= 11) %>%
  group_by(RFID, Date) %>%
  summarise(sumvisits = n()) %>%
  ungroup() %>%
  left_join(weather, by = "Date") %>%
  filter(Date < "2019-03-11" & Date > "2019-01-25")

morn_visits$RFID <- as.factor(morn_visits$RFID)



##dowo
dowovisits=mat_final[,which(colnames(mat_final)%in%rownames(dowoadj))]
dowovisits.mat=as.matrix(dowovisits)

dowoadj.bin=dowoadj
dowoadj.bin[which(dowoadj.bin>0)]=1
colSums(dowoadj.bin[1,]*t(dowovisits.mat), na.rm=T)

dowo.friend.activity=apply(dowoadj, 1, function(x) colSums(x*t(dowovisits.mat), na.rm=T))
dowo.friend.activity=as.data.frame(dowo.friend.activity)
dowo.friend.activity


library(lme4)
library(lmerTest)
library(MuMIn)

dowovisits.dat=as.data.frame(dowovisits)
dowovisits.dat$Date=weather.use$Date
dowovisits.dat$nightlows=weather.use$nightlows

dowo.a= dowovisits.dat %>% gather("ID", "z_score", -Date, -nightlows)

dowo.friend.activity$Date=weather.use$Date
dowo.friend.activity$nightlows=weather.use$nightlows
dowo.b = dowo.friend.activity %>% gather("ID", "z_score_friends", -Date, -nightlows)

dowo.final.dat=merge(dowo.a,dowo.b)

dowomod=lmer(z_score~scale(nightlows)*scale(z_score_friends)+(1|ID), data=dowo.final.dat)
summary(dowomod)

r.squaredGLMM(dowomod)

mod0=lmer(z_score~scale(nightlows) + (1|ID), data=dowo.final.dat)
summary(mod0)
r.squaredGLMM(mod0)

##wbnu
wbnuvisits=mat_final[,which(colnames(mat_final)%in%rownames(wbnuadj))]
wbnuvisits.mat=as.matrix(wbnuvisits)
wbnuadj.trim=wbnuadj[which(rownames(wbnuadj)%in%colnames(mat_final)),which(rownames(wbnuadj)%in%colnames(mat_final))]

colSums(wbnuadj.trim[1,]*t(wbnuvisits.mat), na.rm=T)

wbnu.friend.activity=apply(wbnuadj.trim, 1, function(x) colSums(x*t(wbnuvisits.mat), na.rm=T))
wbnu.friend.activity=as.data.frame(wbnu.friend.activity)
wbnu.friend.activity


library(lme4)
library(lmerTest)

wbnuvisits.dat=as.data.frame(wbnuvisits)
wbnuvisits.dat$Date=weather.use$Date
wbnuvisits.dat$nightlows=weather.use$nightlows

wbnu.a= wbnuvisits.dat %>% gather("ID", "z_score", -Date, -nightlows)

wbnu.friend.activity$Date=weather.use$Date
wbnu.friend.activity$nightlows=weather.use$nightlows
wbnu.b = wbnu.friend.activity %>% gather("ID", "z_score_friends", -Date, -nightlows)

wbnu.final.dat=merge(wbnu.a,wbnu.b)

mod=lmer(z_score~scale(nightlows)*scale(z_score_friends)+(1|ID), data=wbnu.final.dat)
summary(mod)
anova(mod)
r.squaredGLMM(mod)

mod0=lmer(z_score~scale(nightlows) + (1|ID), data=wbnu.final.dat)
summary(mod0)
r.squaredGLMM(mod0)