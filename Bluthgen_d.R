## Script to estimate Blüthgen d' index

library(data.table)
library(plyr)
library(bipartite) 

## Load data with genus-level host trees
# datt
# str(datt)
#    'data.frame':
#     $ SampleID  : chr [1:427]
#     $ no_reads  : num [1:427]
#     $ genus     : chr [1:427]
#     $ ...OTUs

# setDT(datt)
# datt[, no_reads := NULL ]


## Function to prepare data
prep_data <- function(x){

  host_var_type <- colnames(x)[2]
  colnames(x)[2] <- "Host"

  ## Reshape data
  x[, SampleID := NULL ]
  dl <- melt(data = x, id.vars = "Host", variable.name = "OTU", value.name = "Abundance")
  dl <- dl[ Abundance > 0, ]
  dl[ , Occurrence := 1 ]

  ## Reshape to adjacency matrix (rows = host plants, columns = OTUs)
  adjm <- dcast(data = dl,
    formula = Host ~ OTU,
    value.var = "Occurrence", fun.aggregate = sum, fill = 0)

  adjd <- copy(adjm)
  setDF(adjd)
  rownames(adjd) <- adjd$Host
  adjd$Host <- NULL
  adjd <- as.matrix(adjd)

  attr(adjd, which = "host_var_type") <- host_var_type
  return(adjd)
}


DATT <- prep_data(x = datt, with_test_otus = FALSE),

## Estimate d for all datasets
est_d <- function(x){ specieslevel(x, level = "both", index = "d") }
D <- est_d(DATT)

## Extract table with d' values
# metagMisc::dfRowName(x = D$`higher level`, name = "OTU")
