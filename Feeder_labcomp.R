#### PPNC Network and Temperature Paper
require(tidyverse)
require(lubridate)
##data
load("C:/Users/Annie Madsen/Documents/Madsen/R/all_visits.dat")
load("C:/Users/Annie Madsen/Documents/Madsen/R/weather.dat")
weather$Date <- as.Date(weather$Date, format = "%m/%d/%Y")
weather <- weather %>%
  mutate(Hour = hour(datetime)) %>%
  filter(Hour >= 19 | Hour <= 4) %>%
  mutate(Lagdate = ifelse(Hour <= 4, paste0(lag(Date)), paste0(Date)))  %>%
  group_by(Date) %>%
  summarize(nightlows = min(Temp_C)) %>%
  ungroup() %>%
  mutate(index = as.numeric(Date-min(Date)+1))


demo <- read.csv("C:/Users/ShizukaLab/Downloads/RFID_Records_03052019.csv")
demo <- demo %>%
  select(-Date) %>%
  filter(Species != "SCJU")

## quick stats for results section
ndowo <- all_visits %>%
  filter(Species == "DOWO")
nwbnu <- all_visits %>%
  filter(Species == "WBNU")

### Species Models for DOWO and WBNU
morn_visits <- all_visits %>%
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
  filter(Species == "DOWO" | Species == "WBNU")

morn_visits$RFID <- as.factor(morn_visits$RFID)

# newdata <- morn_visits[-9]
# write.csv(newdata, file = "C:/Users/Annie Madsen/Documents/Madsen/Coursework/R Class/feederdata.csv")

require(mgcv)
md <- gamm(sumvisits ~ s(nightlows, fx = FALSE, bs = "tp") + s(RFID, bs = "re"),
           family = poisson,
           data = morn_visits %>% filter(Species == "DOWO"))
summary(md$gam)
summary(md$lme)
plot(md$gam)

mw <- gamm(sumvisits ~ s(nightlows, fx = FALSE, bs = "tp") + s(RFID, bs = "re"),
           family = poisson,
           data = morn_visits %>% filter(Species == "WBNU"))
summary(mw$gam)
summary(mw$lme)
plot(mw$gam)

### species level LMMs
require(lme4)
test_dowo <- glmer(sumvisits ~ scale(nightlows) + (1|RFID), data = morn_visits %>% filter(Species == "DOWO"), family = "poisson")
plot(resid(test))

test_wbnu <- glmer(sumvisits ~ scale(nightlows) + (1|RFID), data = morn_visits %>% filter(Species == "WBNU"), family = "poisson")
plot(resid(test))

library(rsq)
library(MuMIn)
library(arm) 

r.squaredGLMM(test_dowo)
r.squaredGLMM(test_wbnu)

nsim = 10000
bsim = sim(test_dowo, n.sim = nsim)
apply(bsim@fixef,2,mean)
apply(bsim@fixef,2,quantile, prob=c(0.025,0.975))
apply(bsim@ranef,2,mean)
apply(bsim@ranef,2,quantile, prob=c(0.025,0.975))

bsim = sim(test_wbnu, n.sim = nsim)
apply(bsim@fixef,2,mean)
apply(bsim@fixef,2,quantile, prob=c(0.025,0.975))
apply(bsim@ranef,2,mean)
apply(bsim@ranef,2,quantile, prob=c(0.025,0.975))

## Bayes
require(MCMCglmm)
poismod_dowo <- MCMCglmm(sumvisits ~ nightlows,
                         random = ~RFID,
                         family = "poisson", 
                         data = morn_visits %>% filter(Species == "DOWO"),
                         nitt = 220000, ###nitt = burnin + thin*(number of samples to keep)
                         thin = 20,
                         burnin = 2000)

poismod_wbnu <- MCMCglmm(sumvisits ~ nightlows,
                         random = ~RFID,
                         family = "poisson", 
                         data = morn_visits %>% filter(Species == "WBNU"),
                         nitt = 220000, ###nitt = burnin + thin*(number of samples to keep)
                         thin = 20,
                         burnin = 2000)
#####
## Figure 1
## Panel A
md_temp <- morn_visits %>%
  filter(Species == "DOWO")
md_gp <- as.data.frame(predict(md$gam, re.form = TRUE, se = TRUE, type = "response", exclude = s(RFID)))
md_pred <- cbind(md_temp, md_gp)

pana <- ggplot(md_pred) +
  geom_line(aes(nightlows, fit)) +
  geom_ribbon(data = md_pred, aes(x = nightlows, ymin = (fit - 2*se.fit), ymax = (fit + 2*se.fit)), linetype = 2, alpha = 0.2) +
  theme_bw() +
  labs(x = "", y = "", title = "Downy Woodpeckers") +
  theme(plot.title = element_text(size = 20, hjust = 0.5))

## Panel B
mw_temp <- morn_visits %>%
  filter(Species == "WBNU")
mw_gp <- as.data.frame(predict(mw$gam, re.form = TRUE, se = TRUE, type = "response", exclude = s(RFID)))
mw_pred <- cbind(mw_temp, mw_gp)

panb <- ggplot(mw_pred) +
  geom_line(aes(nightlows, fit)) +
  geom_ribbon(data = mw_pred, aes(x = nightlows, ymin = (fit - 2*se.fit), ymax = (fit + 2*se.fit)), linetype = 2, alpha = 0.2) +
  theme_bw() +
  labs(x = "", y = "", title = "White-breasted Nuthatches") +
  theme(plot.title = element_text(size = 20, hjust = 0.5))

test <- rbind(md_pred, mw_pred)
test <- test %>%
  mutate(SpeciesName = ifelse(Species == "DOWO", paste0("Downy Woodpecker"), paste0("White-breasted Nuthatch")))

ggplot(test) +
  geom_line(aes(nightlows, fit)) +
  geom_ribbon(data = test, aes(x = nightlows, ymin = (fit - 2*se.fit), ymax = (fit + 2*se.fit)), linetype = 2, alpha = 0.2) +
  theme_bw() +
  labs(x = "Lowest Overnight Temperature (\u00B0C)", y = "Morning Foraging Activity (Visits/Individual)") +
  facet_wrap("SpeciesName", ncol = 2) +
  theme(strip.text.x = element_text(size = 16), strip.background = element_blank(), strip.placement = "outside", axis.title.x = element_text(size = 13), axis.title.y = element_text(size = 13))


#####
### individual models
require(rsq)
detach(package:mgcv)
require(gam)
RFIDref <- unique(morn_visits$RFID)
vals <- data.frame(RFID = 0, rs = 1:45, pval = 0)

for(i in 1:31){
  ## model
  temp <- morn_visits %>% filter(RFID == RFIDref[i])
  Mtemp <- gam(sumvisits ~ lo(nightlows, span = 0.6), data = temp)
  Rtemp <- rsq(Mtemp)
  
  ## plot predictions over raw data
  predMtemp <- as.data.frame(predict(Mtemp, re.form = TRUE, se = TRUE))
  dftemp <- cbind(temp, predMtemp)
  
  jpeg(file = paste("GAM Fit for ", RFIDref[i], ".jpeg"),
       width = 800,
       height = 480)
  
  plot <- ggplot() +
    geom_point(data = dftemp, aes(x = nightlows, y = sumvisits)) +
    geom_line(data = dftemp, aes(x = nightlows, y = fit)) +
    geom_ribbon(data = dftemp, aes(x = nightlows, ymin = (fit - 2*se.fit), ymax = (fit + 2*se.fit)), linetype = 2, alpha = 0.2, color = "black") +
    labs(x = "Lowest Overnight Temperature (C)", y = "Morning Visitation Rate", title = paste0("GAM Fit for ", RFIDref[i]))
  print(plot)
  graphics.off()
  
  assign(paste0("M_", RFIDref[i]), Mtemp)
  assign(paste0("Plot_", i), plot)
  assign(paste0("R", i), Rtemp)
  
  vals$RFID[i] <- RFIDref[i]
  vals$rs[i] <- Rtemp
  vals$pval[i] <- (summary(Mtemp)$anova[,3][2])
}


## Figure S1
dowo_df <- morn_visits %>% 
  filter(Species == "DOWO")
ggplot(dowo_df, aes(nightlows, sumvisits)) +
  geom_point() +
  geom_smooth() +
  facet_wrap("RFID") +
  theme_bw() +
  labs(title = "Downy Woodpeckers", x = "Lowest Overnight Temperature (\u00B0C)" , y = "Morning Foraging Activity (total visits)") +
  theme(plot.title = element_text(size = 22, face = "bold"), axis.title.x = element_text(size = 22, face = "bold"), axis.title.y = element_text(size = 22, face = "bold"), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))

wbnu_df <- morn_visits %>% 
  filter(Species == "WBNU")
ggplot(wbnu_df, aes(nightlows, sumvisits)) +
  geom_point() +
  geom_smooth() +
  facet_wrap("RFID") +
  theme_bw() +
  labs(title = "White-breasted Nuthatches", x = "Lowest Overnight Temperature (\u00B0C)" , y = "Morning Foraging Activity (total visits)") +
  theme(plot.title = element_text(size = 22, face = "bold"), axis.title.x = element_text(size = 22, face = "bold"), axis.title.y = element_text(size = 22, face = "bold"), axis.text.x = element_text(size = 13), axis.text.y = element_text(size = 13))


#####
## Figure 4? 
## Start by making daily networks
## index data by day
all_visits <- all_visits %>%
  ungroup() %>%
  mutate(LoggerDate = paste0(Logger, Date)) %>%
  mutate(index = as.numeric(Date- min(Date)+1)) %>%
  mutate(Timestamp = as.numeric(Datetime))

## reformat visitation data for gmmevents
require(asnipe)
dat.dowo <- all_visits %>%
  ungroup() %>%
  filter(Species == "DOWO")
dat.wbnu <- all_visits %>%
  ungroup() %>%
  filter(Species == "WBNU")
DOWOflocks = list()
for(i in 1:max(dat.dowo$index)){
  dowo_temp = dat.dowo %>% filter(index == i)
  if (nrow(dowo_temp)<15 | !is.null(dat.dowo$index[i])) next else{
    DOWOflocks[[i]] = gmmevents(dowo_temp$Timestamp, dowo_temp$RFID, dowo_temp$LoggerDate) 
  }
}
save(DOWOflocks, file = "C:/Users/ShizukaLab/Downloads/DOWOflocks.r")

WBNUflocks = list()
for(i in 1:max(dat.wbnu$index)){
  wbnu_temp = dat.wbnu %>% filter(index == i)
  if (nrow(wbnu_temp)<15) next else{
    WBNUflocks[[i]] = gmmevents(wbnu_temp$Timestamp, wbnu_temp$RFID, wbnu_temp$LoggerDate) 
  }
}
save(WBNUflocks, file = "C:/Users/ShizukaLab/Downloads/WBNUflocks.r")

## Daily networks
library(igraph)
library(ndtv)
library(asnipe)
library(visNetwork)
library(ggmap)



birddat = demo[which(demo$Species!=""),] #getting rid of false RFID tags

## DOWO
dowo.gbi = list()
dowo.net = list()
dowo.g = list()
for(i in 1:length(DOWOflocks)){
  if(!is.null(DOWOflocks[[i]])){
    dowo.gbi[[i]]=DOWOflocks[[i]]$gbi[,which(colnames(DOWOflocks[[i]]$gbi)%in%birddat$RFID)]
    dowo.net[[i]] = get_network(dowo.gbi[[i]], association_index = "SRI")
    dowo.g[[i]] = graph_from_adjacency_matrix(dowo.net[[i]], "undirected", weighted=T)
    plot(dowo.g[[i]], vertex.label="", edge.width=E(dowo.g[[i]])$weight*20) ## all daily networks
  }
}

### calculate edge density
dowo.n = list()
dowo.m = list()
dowo.dyads = list()
dowo.density = list()
for(i in 1:length(dowo.g)){
  if(!is.null(dowo.g[[i]])){
    dowo.n[[i]]=vcount(dowo.g[[i]]) ## total number of individuals for the day
    dowo.m[[i]]=ecount(dowo.g[[i]])
    dowo.dyads[[i]]=dowo.n[[i]]*(dowo.n[[i]]-1)/2 ## number of possible edges based on number of vertices
    dowo.density[[i]]=m[[i]]/dowo.dyads[[i]] ## ratio of the number of edges and the number of possible edges 
  }
}

## WBNU
wbnu.gbi = list()
wbnu.net = list()
wbnu.g = list()
for(i in 1:length(WBNUflocks)){
  if(!is.null(WBNUflocks[[i]])){
    wbnu.gbi[[i]]=WBNUflocks[[i]]$gbi[,which(colnames(WBNUflocks[[i]]$gbi)%in%birddat$RFID)]
    wbnu.net[[i]] = get_network(wbnu.gbi[[i]], association_index = "SRI")
    wbnu.g[[i]] = graph_from_adjacency_matrix(wbnu.net[[i]], "undirected", weighted=T)
    plot(wbnu.g[[i]], vertex.label="", edge.width=E(wbnu.g[[i]])$weight*20) ## all daily networks
  }
}

### calculate edge density
wbnu.n = list()
wbnu.m = list()
wbnu.dyads = list()
wbnu.density = list()
for(i in 1:length(wbnu.g)){
  if(!is.null(wbnu.g[[i]])){
    wbnu.n[[i]]=vcount(wbnu.g[[i]]) ## total number of individuals for the day
    wbnu.m[[i]]=ecount(wbnu.g[[i]])
    wbnu.dyads[[i]]=wbnu.n[[i]]*(wbnu.n[[i]]-1)/2 ## number of possible edges based on number of vertices
    wbnu.density[[i]]=m[[i]]/wbnu.dyads[[i]] ## ratio of the number of edges and the number of possible edges 
  }
}

## filter weather days we need
weath_net <- weather[which(weather$Date %in% all_visits$Date),]
weath_net <- weath_net %>%
  mutate(index = as.numeric(Date-min(Date)+1)) %>%
  filter(index != 48)

##plot density vs. temperature
plot(weath_net$nightlows, unlist(dowo.density))
plot(weath_net$nightlows, unlist(dowo.n))
plot(unlist(dowo.n), unlist(dowo.density)) ## should be correlated, not cancelling out and causing a null relationship

weath_net.wbnu <- weath_net %>%
  filter(index != 45)
plot(weath_net.wbnu$nightlows, unlist(wbnu.density))
plot(weath_net.wbnu$nightlows, unlist(wbnu.n))
plot(unlist(wbnu.n), unlist(wbnu.density))


######
## 01/07/2020 Need to keep same order of individuals in the circles between warm & cold plots
## 5 warmest days
## filter data 
warm <- all_visits %>%
  filter(index == 7 | index == 40 | index == 33 | index == 2 | index == 39)
## flocking events
## DOWO
dowo_warm <- warm %>%
  filter(Species == "DOWO")
dowo_wgmm = gmmevents(dowo_warm$Timestamp, dowo_warm$RFID, dowo_warm$LoggerDate)
## WBNU
wbnu_warm <- warm %>%
  filter(Species == "WBNU")
wbnu_wgmm = gmmevents(wbnu_warm$Timestamp, wbnu_warm$RFID, wbnu_warm$LoggerDate)

## networks
## DOWO
dowo_wgbi = dowo_wgmm$gbi[,which(colnames(dowo_wgmm$gbi)%in%birddat$RFID)]
dowo_wnet = get_network(dowo_wgbi, association_index = "SRI")
dowo_wg = graph_from_adjacency_matrix(dowo_wnet, "undirected", weighted=T)
set.seed(5)
plot(dowo_wg, layout = layout_in_circle, vertex.label="", vertex.color = "red", edge.width=E(dowo_wg)$weight*20)
## WBNU
wbnu_wgbi = wbnu_wgmm$gbi[,which(colnames(wbnu_wgmm$gbi)%in%birddat$RFID)]
wbnu_wnet = get_network(wbnu_wgbi, association_index = "SRI")
wbnu_wg = graph_from_adjacency_matrix(wbnu_wnet, "undirected", weighted=T)
set.seed(5)
plot(wbnu_wg, layout = layout_in_circle, vertex.label="", vertex.color = "red", edge.width=E(wbnu_wg)$weight*20)

## 5 coldest days
## filter data 
cold <- all_visits %>%
  filter(index == 10 | index == 11 | index == 45 | index == 18 | index == 19)
## flocking events
## DOWO
dowo_cold <- cold %>%
  filter(Species == "DOWO")
dowo_cgmm = gmmevents(dowo_cold$Timestamp, dowo_cold$RFID, dowo_cold$LoggerDate)
## WBNU
wbnu_cold <- cold %>%
  filter(Species == "WBNU")
wbnu_cgmm = gmmevents(wbnu_cold$Timestamp, wbnu_cold$RFID, wbnu_cold$LoggerDate)

## networks
## DOWO
dowo_cgbi = dowo_cgmm$gbi[,which(colnames(dowo_cgmm$gbi)%in%birddat$RFID)]
## deal with missing individuals
new_dowo_cgi = matrix(data = rep(0, 335), nrow = 335, ncol = 3)
row_names = c(1:335)
col_names = c("011017396E","0700E0FFDB", "0700EE0805")
dimnames(new_dowo_cgi) <- list(row_names, col_names)
dowo_cgbi_final <- cbind(dowo_cgbi, new_dowo_cgi)
## moving on
dowo_cnet = get_network(dowo_cgbi_final, association_index = "SRI")
dowo_cg = graph_from_adjacency_matrix(dowo_cnet, "undirected", weighted=T)
set.seed(5)
plot(dowo_cg, layout = layout_in_circle, vertex.label = "", vertex.color = "blue", edge.width=E(dowo_cg)$weight*20)
## WBNU
wbnu_cgbi = wbnu_cgmm$gbi[,which(colnames(wbnu_cgmm$gbi)%in%birddat$RFID)]
## deal with missing individuals
new_wbnu_cgi = matrix(data = rep(0, 204), nrow = 204, ncol = 4)
row_names = c(1:204)
col_names = c("01101706AD","0110173C26", "0110175F3E","0700EE3187")
dimnames(new_wbnu_cgi) <- list(row_names, col_names)
wbnu_cgbi_final <- cbind(wbnu_cgbi, new_wbnu_cgi)
## moving on
wbnu_cnet = get_network(wbnu_cgbi_final, association_index = "SRI")
wbnu_cg = graph_from_adjacency_matrix(wbnu_cnet, "undirected", weighted=T)
set.seed(5)
plot(wbnu_cg, layout = layout_in_circle, vertex.label = "", vertex.color ="blue", edge.width=E(wbnu_cg)$weight*20)



##daily average flock size; row sums/# cols for each day
avg_flocks <- list()
for(i in 1:length(DOWOflocks)){
  if(!is.null(DOWOflocks[[i]]$gbi)){
    for(j in DOWOflocks[[i]]$gbi[j,]){ ### messed up something here
      avg_flocks[[i]] = sum(DOWOflocks[[i]]$gbi[j,])/ncol(DOWOflocks[[i]]$gbi)
    }
  }
}


##number flocks/feeder/day; # rows = #gmm events for that day, do this per logger



## temp responses vs. assortativity analysis:positive if similar vertices tend to connect to each other, negative if not




