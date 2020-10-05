### NUll model testing of MRQAP results

## This time with single-species networks
library(asnipe)
library(igraph)
library(iterators)
library(foreach)
library(parallel)
library(doParallel)
library(lubridate)
library(tidyverse)
library(reshape)
load("all_visits.dat")
#import files
gmmDOWO=readRDS("conspecificDOWOflocks.rds") #import gmm results file. Will need this for permuting group-by-individual matrices.
gmmWBNU=readRDS("conspecificWBNUflocks.rds") #import gmm results file. Will need this for permuting group-by-individual matrices.

gbi_dowo=gmmDOWO$gbi
gbi_wbnu=gmmWBNU$gbi

dowoadj=get_network(gbi_dowo)
wbnuadj=get_network(gbi_wbnu)


diag(dowoadj)=NA #make diagonal of adjacency matrices NA so we don't count these when normalizing values later.
diag(wbnuadj)=NA #make diagonal of adjacency matrices NA so we don't count these when normalizing values later.

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

mat_ddm <- cast(dv, Date ~ RFID, value = "dailyvisits")
mat_ddm[is.na(mat_ddm)] <- 0
mat <- mat_ddm[2:46]/dvsd

myvec <- (ref$mdv[match(names(mat_ddm[2:45]), ref$RFID)])/dvsd

mat_final <- mat[1] - myvec[1]
for(i in 2:44){
  mat_temp <- mat[i] - myvec[i]
  mat_final <- cbind(mat_final, mat_temp)
}
#mat_final <- cbind(mat_ddm$Date, mat_final)

### Similarity Matrix
## how similar are individuals' daily visitation z-scores?

simmat <- as.matrix(simil(mat_final, by_rows = FALSE))

#########Spatial Overlap Matrix
### summarise number of visits at each feeder for each bird
vis <- all_visits %>%
  group_by(Logger, RFID) %>%
  summarise(logvis = n()) %>%
  ungroup()

### reshape data frame and calculate proportion of visits at each feeder
logsums <- cast(vis, Logger ~ RFID, value = "logvis")
logsums[is.na(logsums)] <- 0
y = colSums(logsums)
fin <- as.data.frame(mapply("/", logsums[-1], y))

### make correlation/similarity matrix
require(proxy)
logmat <- as.matrix(simil(fin, by_rows = FALSE))



dowosim=simmat[match(rownames(dowoadj), rownames(simmat)), match(rownames(dowoadj), rownames(simmat))] #sort the activity correlation matrix so rows/columns match adjacency matrix
dowospat=logmat[match(rownames(dowoadj), rownames(logmat)), match(rownames(dowoadj), rownames(logmat))] #sort spatial correlation matrix so rows/columns match adjacency matrix. 


#same for WBNU
wbnusim=simmat[match(rownames(wbnuadj), rownames(simmat)), match(rownames(wbnuadj), rownames(simmat))]
wbnuspat=logmat[match(rownames(wbnuadj), rownames(logmat)), match(rownames(wbnuadj), rownames(logmat))]


#MRQAP-without normalizing values.

mrqap.dsp(dowosim~dowoadj+dowospat) #using activity similarity as the response variable and network + logger similarity as a covariates. This is raw numbers, not normalized.

mrqap.dsp(wbnusim~wbnuadj+wbnuspat)

# #now, normalize all matrix values so that minimum number = 0 and maximum number = 1
normalize_matrix=function(m){
  (m-min(m, na.rm=T))/(max(m, na.rm=T)-min(m, na.rm=T))
}

dowosim.norm=normalize_matrix(dowosim)
dowoadj.norm=normalize_matrix(dowoadj)
dowospat.norm=normalize_matrix(dowospat)

dowo.mrqap.norm=mrqap.dsp(dowosim.norm~dowoadj.norm+dowospat.norm) #same test, but now with normalized values. The results are the same but the coefficient is different.
dowo.mrqap.norm
str(dowo.mrqap.norm)
mean(dowo.mrqap.norm$fitted.values)


###dowo permutations

dowo.ids=rownames(dowoadj)

dowogbi=gbi_dowo[,which(colnames(gbi_dowo)%in%dowo.ids)] #get gbi with only DOWOs
dowogbi.filt=dowogbi[which(rowSums(dowogbi)>0),] #remove groups that no DOWOs belong to.
  
dowometadata.filt=gmmDOWO$metadata[which(rowSums(dowogbi)>0),] #remove the same groups in the group metadata

#using the filtered metadata, we can extract the feeder & date of the group. This will be useful when we constrain the permutation by day
dowo.locations=as.numeric(as.factor(substr(dowometadata.filt$Location, start=1, stop=8)))
dowo.days=as.numeric(as.factor(substr(dowometadata.filt$Location, start=10, stop=13)))

#Run 10,000 group membership permutations with swaps constrained by date (but not feeder)
# dowoperm=network_permutation(dowogbi.filt, data_format="GBI", permutation=10000, association_matrix = dowoadj, locations=dowo.locations, days=dowo.days, within_day=TRUE, within_location=FALSE)
# 
# plot(apply(dowoperm, 1, max), type="l", ylab="maximum edge weight", xlab="# swaps", main="DOWO permutations")

#store the results of MRQAP with empirical network
emp.mod.dowo=mrqap.dsp(dowosim.norm~dowoadj.norm+dowospat.norm)
emp.coef.dowo=emp.mod.dowo$coefficients[2]

#now do 1000 sets of swaps and store results
normalize_matrix=function(m){
  (m-min(m, na.rm=T))/(max(m, na.rm=T)-min(m, na.rm=T))
}


### DOWO with parallel processing!
library(foreach)
library(parallel)
library(doParallel)
times=1000
n.cores=detectCores()
system.time({
  registerDoParallel(n.cores)
  dowo.results.parallel=foreach(i = 1:times) %dopar% as.matrix(network_swap(dowogbi.filt, data_format="GBI", swaps=10000, association_matrix = dowoadj, locations=dowo.locations, days=dowo.days, within_day=TRUE, within_location=FALSE)$Association_index)
})
stopImplicitCluster()

dowoperm.adjs=dowo.results.parallel

dowoperm.mrqap.coefs=vector(length=times)
for(j in 1:times){
  adj=dowoperm.adjs[[j]]
  diag(adj)=NA
  adj.norm=normalize_matrix(adj)
  dowoperm.mrqap.coefs[j]=mrqap.dsp(dowosim.norm~adj.norm+dowospat.norm, randomisations=1)$coefficients[2]
}

dowoperm.mrqap.coefs

p.dowo=(length(which(dowoperm.mrqap.coefs>=emp.coef.dowo))+1)/(times+1)
p.dowo #p-value

ci.dowo=quantile(dowoperm.mrqap.coefs, probs=c(0.025, 0.975))
ci.dowo #confidence interval of null model -- this likely will not overlap the empirical coefficient value


save(dowoperm.adjs, dowoperm.mrqap.coefs, p.dowo, ci.dowo, file="dowo_results_20201001_rdat")
###WBNU
wbnu.ids=rownames(wbnuadj)

wbnugbi=gbi_wbnu[,which(colnames(gbi_wbnu)%in%wbnu.ids)] #get gbi with only wbnus
wbnugbi.filt=wbnugbi[which(rowSums(wbnugbi)>0),] #remove groups that no wbnus belong to.

wbnumetadata.filt=gmmWBNU$metadata[which(rowSums(wbnugbi)>0),] #remove the same groups in the group metadata

#using the filtered metadata, we can extract the feeder & date of the group. This will be useful when we constrain the permutation by day
wbnu.locations=as.numeric(as.factor(substr(wbnumetadata.filt$Location, start=1, stop=8)))
wbnu.days=as.numeric(as.factor(substr(wbnumetadata.filt$Location, start=10, stop=13)))


#store the results of MRQAP with empirical network
wbnusim.norm=normalize_matrix(wbnusim)
wbnuadj.norm=normalize_matrix(wbnuadj)
wbnuspat.norm=normalize_matrix(wbnuspat)

emp.mod.wbnu=mrqap.dsp(wbnusim.norm~wbnuadj.norm+wbnuspat.norm)
emp.coef.wbnu=emp.mod.wbnu$coefficients[2]

### WBNU with parallel processing!
library(foreach)
library(parallel)
library(doParallel)
times=1000
n.cores=detectCores()
system.time({
  registerDoParallel(n.cores)
wbnu.results.parallel=foreach(i = 1:times) %dopar% as.matrix(network_swap(wbnugbi.filt, data_format="GBI", swaps=10000, association_matrix = wbnuadj, locations=wbnu.locations, days=wbnu.days, within_day=TRUE, within_location=FALSE)$Association_index)
})
stopImplicitCluster()

wbnuperm.adjs=wbnu.results.parallel

wbnuperm.mrqap.coefs=vector(length=times)
for(j in 1:times){
  adj=wbnuperm.adjs[[j]]
  diag(adj)=NA
  adj.norm=normalize_matrix(adj)
  wbnuperm.mrqap.coefs[j]=mrqap.dsp(wbnusim.norm~adj.norm+wbnuspat.norm, randomisations=1)$coefficients[2]
}

wbnuperm.mrqap.coefs

p.wbnu=(length(which(wbnuperm.mrqap.coefs>=emp.coef.wbnu))+1)/(times+1)
p.wbnu #p-value

ci.wbnu=quantile(wbnuperm.mrqap.coefs, probs=c(0.025, 0.975))
ci.wbnu #confidence interval of null model -- this likely will not overlap the empirical coefficient value


save(wbnuperm.adjs, wbnuperm.mrqap.coefs, p.wbnu, ci.wbnu, file="wbnu_results_20201001_rdat")
