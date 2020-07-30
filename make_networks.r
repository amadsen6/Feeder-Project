#networks

library(asnipe)
library(igraph)

dowo_gmm=readRDS("conspecificDOWOflocks.rds")
str(dowo_gmm)

dowo_gbi=dowo_gmm$gbi

wbnu_gmm=readRDS("conspecificWBNUflocks.rds")

wbnu_gbi=wbnu_gmm$gbi

dowo_net=get_network(dowo_gbi)

