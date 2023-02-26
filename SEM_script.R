library(semPlot)
library(nlme)
library(lavaan)
library(lavaanPlot)
SEMM<-read.table("clipboard", header = T, sep = "\t", dec = ",")

lav_ECMR1<- ' EcMFsppTot ~Genus+EcM_Pielou+Age+ECM_BA+EcM_perc+ph+ph2+delta15N+EcMPspp+YEAR
delta15N ~Genus+ECM_BA+EcM_perc+ph+ph2
ph~ Genus+Age+ECM_BA+EcM_perc
ph2~ Genus+Age+ECM_BA+EcM_perc
EcM_Pielou~ Genus+EcMPspp+ EcMFsppTot+Age+ECM_BA+EcM_perc+ph+ph2+delta15N+YEAR
ECM_BA~ Genus+EcMPspp
EcM_perc~ Genus+EcMPspp+Age+ECM_BA+ ph+ph2+delta15N+YEAR
EcMPspp ~Age'

summary(fit.lav_ECMR1,fit.measures=T,standardized=T,rsq=T)

semPaths(object = fit.lav_ECMR1,
         layout = "spring",
         rotation = 1,
         whatLabels = "std",
         edge.label.cex = 1,
         what = "std",
         edge.color = "navy",
         sizeMan = 10, nCharNodes=8,title=T)
semPaths(fit.lav_ECMR1, whatLabels="std", style="lisrel", exoCov = T, curvePivot = TRUE, sizeMan = 8, sizeInt = 12, 
         residuals=F,nCharNodes=8) 

modificationindices(fit.lav_mod1, sort = T)

lavaanPlot(model = fit.lav_ECMR6, node_options = list(shape = "box", fontname = "Helvetica"), edge_options = list(color = "grey"), coefs = T,covs = TRUE)
		 
		 