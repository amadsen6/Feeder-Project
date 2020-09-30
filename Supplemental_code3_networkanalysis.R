#### Analysis code, Part 1 for Madsen, Vander Meiden & Shizuka: Social partners and temperature jointly affect morning foraging activity of small birds in winter#####

#generates network plots for Figure 1
#measures assortment and social differentiation

library(asnipe)
library(igraph)
library(assortnet)

#import Gaussian Mixture Model Results
dowo_gmm=readRDS("conspecificDOWOflocks.rds")
wbnu_gmm=readRDS("conspecificWBNUflocks.rds")

#extract group-by-individual matrices
dowo_gbi=dowo_gmm$gbi
wbnu_gbi=wbnu_gmm$gbi

#make networks
dowo_net=graph_from_adjacency_matrix(get_network(dowo_gbi), mode="undirected", weighted=T)
wbnu_net=graph_from_adjacency_matrix(get_network(wbnu_gbi), mode="undirected", weighted=T)

#import individual attributes data & match up with network vertices
indivs=read.csv("RFID_Records_fixed.csv")

V(dowo_net)$sex=indivs[match(V(dowo_net)$name, indivs$RFID),"Sex"]
V(wbnu_net)$sex=indivs[match(V(wbnu_net)$name, indivs$RFID),"Sex"]
V(dowo_net)$age=indivs[match(V(dowo_net)$name, indivs$RFID),"Age"]

#make network plots for Figure 1
sex_color=data.frame(sex=c("F", "M", "U"), color=c("yellow", "purple", "white"))
plot(dowo_net, vertex.color=sex_color[match(V(dowo_net)$sex, sex_color$sex), "color"], vertex.label="", edge.width=E(dowo_net)$weight*30)

plot(wbnu_net, vertex.color=sex_color[match(V(wbnu_net)$sex, sex_color$sex), "color"], vertex.label="", edge.width=E(wbnu_net)$weight*30)

#assortment by sex for DOWO
sexassort_dowo=assortment.discrete(as_adj(dowo_net, sparse=F, attr="weight"), V(dowo_net)$sex, SE=T)
sexassort_dowo$r


#do node permutations and generate p-value and confidence interval for DOWO
random_sex_dowo=lapply(1:1000, function(x) sample(V(dowo_net)$sex, length(V(dowo_net)$sex), replace=F))
assort_rand_dowo=sapply(random_sex_dowo, function(x) assortment.discrete(as_adj(dowo_net, sparse=F, attr="weight"), x, SE=F)$r)
p_assort_dowo=length(which(assort_rand_dowo<sexassort_dowo$r))/1001
ci_assort_rand_dowo=quantile(assort_rand_dowo, probs = c(0.025, 0.925))
p_assort_dowo
ci_assort_rand_dowo

#assortment by sex for WBNU
sexassort_wbnu=assortment.discrete(as_adj(wbnu_net, attr="weight", sparse=F), V(wbnu_net)$sex, SE=T)
sexassort_wbnu$r

#do node permutations and generate p-value and confidence interval for WBNU
random_sex_wbnu=lapply(1:1000, function(x) sample(V(wbnu_net)$sex, length(V(wbnu_net)$sex), replace=F))
assort_rand_wbnu=sapply(random_sex_wbnu, function(x) assortment.discrete(as_adj(wbnu_net, attr="weight", sparse=F), x, SE=F)$r)
p_assort_wbnu=length(which(assort_rand_wbnu<sexassort_wbnu$r))/1001
ci_assort_rand_wbnu=quantile(assort_rand_wbnu, probs = c(0.025, 0.925))
mean(assort_rand_wbnu)
p_assort_wbnu
ci_assort_rand_wbnu


##assortment from group permutations.
load("dowo_results_20190814.rdat")
load("wbnu_results_20190814.rdat")
# 

cv=function(x) sd(x)/mean(x)
cv_dowo_emp=cv(E(dowo_net)$weight)
cv_dowo_rand=sapply(dowoperm.adjs, function(y) cv(y))
hist(cv_dowo_rand, xlim=c(0.65, 1.3))
lines(c(cv_dowo_emp, cv_dowo_emp), c(0,200), lty=2, col="red", lwd=2)
quantile(cv_dowo_rand, probs=c(0.025, 0.975))

cv_wbnu_emp=cv(E(wbnu_net)$weight)
cv_wbnu_rand=sapply(wbnuperm.adjs, function(y) cv(y))
hist(cv_wbnu_rand, xlim=c(0.8, 2.2))
lines(c(cv_wbnu_emp, cv_wbnu_emp), c(0,200), lty=2, col="red", lwd=2)
quantile(cv_wbnu_rand, probs=c(0.025, 0.975))
