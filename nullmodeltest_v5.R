### NUll model testing of MRQAP results

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


diag(dowoadj)=NA #make diagonal of adjacency matrices NA so we don't count these when normalizing values later.
diag(wbnuadj)=NA #make diagonal of adjacency matrices NA so we don't count these when normalizing values later.
simmat=as.matrix(read.csv("simmat.csv", row.names = 1, check.names = F)) #import correlation matrix of daily activity (Z-scores)
logmat=as.matrix(read.csv("logmat.csv", row.names = 1, check.names = F)) #import correlation matrix of proportional feeder use

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

times=1000
dowoperm.adjs=list()
for(i in 1:times){
dowoperm_repeat=network_swap(dowogbi.filt, data_format="GBI", swaps=10000, association_matrix = dowoadj, locations=dowo.locations, days=dowo.days, within_day=TRUE, within_location=FALSE)
dowoperm.adjs[[i]]=dowoperm_repeat$Association_index
}

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

save(dowoperm.adjs, dowoperm.mrqap.coefs, p.dowo, ci.dowo, file="dowo_results_20190814_rdat")


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

#now do 1000 sets of swaps and store results
# times=1000
# wbnuperm.adjs=list()
# for(i in 1:times){
#   wbnuperm_repeat=network_swap(wbnugbi.filt, data_format="GBI", swaps=10000, association_matrix = wbnuadj, locations=wbnu.locations, days=wbnu.days, within_day=TRUE, within_location=FALSE)
#   wbnuperm.adjs[[i]]=wbnuperm_repeat$Association_index
# }
# 
# wbnuperm.mrqap.coefs=vector(length=times)
# for(j in 1:times){
#   adj=wbnuperm.adjs[[j]]
#   diag(adj)=NA
#   adj.norm=normalize_matrix(adj)
#   wbnuperm.mrqap.coefs[j]=mrqap.dsp(wbnusim.norm~adj.norm+wbnuspat.norm, randomisations=1)$coefficients[2]
# }
# 
# wbnuperm.mrqap.coefs
# 
# p.wbnu=(length(which(wbnuperm.mrqap.coefs>=emp.coef.wbnu))+1)/(times+1)
# p.wbnu #p-value
# 
# ci.wbnu=quantile(wbnuperm.mrqap.coefs, probs=c(0.025, 0.975))
# ci.wbnu #confidence interval of null model -- this likely will not overlap the empirical coefficient value
# 
# 
# save(wbnuperm.adjs, wbnuperm.mrqap.coefs, p.wbnu, ci.wbnu, file="dowo_results_20190814_rdat")


### WBNU with parallel processing!
library(foreach)
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


save(wbnuperm.adjs, wbnuperm.mrqap.coefs, p.wbnu, ci.wbnu, file="wbnu_results_20190814_rdat")
