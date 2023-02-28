#Biparite Network analysis

library(igraph)
library(ggplot2)
library(reshape2)

# Data
taxon.data=read.csv("EcMFall.csv", header = T)
OTU.ee1=taxon.data
rownames(OTU.ee1)=OTU.ee1$OTU
OTU.ee1=as.data.frame(t(OTU.ee1))
taxon2=OTU.ee1
taxon2$OTU=NULL
taxon2$lineage=as.factor(taxon2$lineage)
taxon3 <- data.frame(sapply(taxon2, function(x) as.numeric(as.character(x))))
taxon3$lineage=NULL
taxon2=taxon3
row.names(taxon2)=OTU.ee1$OTU
table(taxon2)
#keep only common OTUs
taxon3=as.data.frame(t(taxon3))
colnames(taxon3)=OTU.ee1$OTU
d1=sapply(taxon3, function(x) sum(x>0) )
d1=as.data.frame(d1)
sum(d1$d1>5)
tax5 =d1$d1>5
#more abundant NOG data
d2=taxon3[,tax5]
taxon4=taxon2[colnames(d2),]

# ABUNDANT OTUs
OTU.ee=taxon4
OTU.ee$lineage=NULL
d3=aggregate(OTU.ee, by=list(taxon4$lineage), FUN=sum)
row.names(d3)=d3$Group.1
d3$Group.1=NULL
Lineage.d=as.data.frame(t(d3))

### Network
otu.t=as.data.frame(t(OTU.ee))
OTU.net=aggregate(otu.t, by=list(data.ee$genus), FUN=sum)
class(OTU.net)
row.names(OTU.net)=OTU.net$Group.1
OTU.net$Group.1=NULL

otu.net2=as.matrix(t(OTU.net))
# make the network
two.mod.network2= graph.incidence(otu.net2)
summary(two.mod.network)
plot.igraph(two.mod.network2, vertex.color=c("green","pink")[V(two.mod.network)$type+1],
            layout = layout_with_mds, vertex.label.cex=0.5,edge.arrow.size=0.3,vertex.label=NA)
# export network graph
write_graph(two.mod.network2,
            "OTU two mode network.graphml",
            format = "graphml")
dev.off()

#add variables to the graph in gephi
otu.names=colnames(OTU.net)
otu.tax=taxon2[,otu.names]
otu.tax=t(otu.tax)
otu.tax=otu.tax[,c(1:3)]
otu.tax=as.data.frame(otu.tax)
otu.tax$lineage=as.factor(otu.tax$lineage)
otu.tax$lineage2=otu.tax$lineage
levels(otu.tax$lineage2)=list("tuber-helvella"="tuber-helvella","sebacina"="sebacina",
                              "amphinema-tylospora"="amphinema-tylospora","hebeloma-alnicola"="hebeloma-alnicola",
                              "russula-lactarius"="russula-lactarius","cortinarius"="cortinarius","tomentella-thelephora"="tomentella-thelephora",
                              "inocybe"="inocybe","Others"=c("tulasnella1",
                                                             "ceratobasidium2",
                                                             "serendipita2",
                                                             "helotiales1",
                                                             "tremellodendropsis",
                                                             "byssocorticium",
                                                             "leotia",
                                                             "rhodoscypha",
                                                             "hydnellum-sarcodon",
                                                             "hysterangium",
                                                             "tomentellopsis",
                                                             "marcelleina-peziza_gerardii",
                                                             "elaphomyces",
                                                             "pachyphloeus-amylascus",
                                                             "pustularia",
                                                             "hygrophorus",
                                                             "entoloma",
                                                             "sphaerosporella",
                                                             "serendipita1",
                                                             "boletus",
                                                             "cantharellus",
                                                             "otidea",
                                                             "tarzetta",
                                                             "terfezia-peziza_depressa",
                                                             "pisolithus-scleroderma",
                                                             "amanita",
                                                             "galactinia",
                                                             "paxillus-gyrodon",
                                                             "sordariales2",
                                                             "pulvinula",
                                                             "geopora",
                                                             "cenococcum",
                                                             "clavulina",
                                                             "suillus-rhizopogon",
                                                             "laccaria",
                                                             "tricholoma",
                                                             "genea-humaria",
                                                             "pseudotomentella",
                                                             "piloderma",
                                                             "meliniomyces",
                                                             "wilcoxina"))
#export for Gephi
write.csv(otu.tax,"OTUtax network data.csv")
write_graph(two.mod.network2,
            "OTU two mode network with lineage.graphml",
            format = "graphml")
dev.off()


