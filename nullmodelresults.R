#extract results from null model analysis

load("dowo_results_20190814.rdat")


mean(dowoperm.mrqap.coefs)
ci.dowo

load("wbnu_results_20190814.rdat")
mean(wbnuperm.mrqap.coefs)
ci.wbnu
