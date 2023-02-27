
## Partner specificity and association networks in ectomycorrhizal symbiosis


##############################################
############################################## Indicator values for all species, genus and lineage/family levels
##############################################

library(indicspecies)
library(data.table)
library(plyr)
library(ggplot2)
library(openxlsx)

set.seed(14789)
theme_set(theme_classic(base_size = 14))

### Load the data

## Species-level host trees
datt_sp <- readxl::read_excel(
  path = "EstLatPhi_v2.xlsx",
  sheet = "R_species17",
  na = c("", " ", "NA", "#N/A"))

## Genus-level host trees
datt_gn <- readxl::read_excel(
  path = "EstLatPhi_v2.xlsx",
  sheet = "R_genus9",
  na = c("", " ", "NA", "#N/A"))

## Linage/Family-level host trees
datt_li <- readxl::read_excel(
  path = "EstLatPhi_v2.xlsx",
  sheet = "R_lineage4",
  na = c("", " ", "NA", "#N/A"))

## Angiosperms vs Gymnosperms
datt_dc <- readxl::read_excel(
  path = "EstLatPhi_v2.xlsx",
  sheet = "R_phyl2",
  na = c("", " ", "NA", "#N/A"))

setDT(datt_sp)
setDT(datt_gn)
setDT(datt_li)
setDT(datt_dc)

datt_sp[, no_reads := NULL ]
datt_gn[, no_reads := NULL ]
datt_li[, no_reads := NULL ]
datt_dc[, no_reads := NULL ]

## Convert to presence-absences
colz <- colnames(datt_sp)[ -c(1,2) ]
datt_sp[, c(colz) := lapply(.SD, function(x){ ifelse(x > 0, 1, 0) }), .SDcols = colz]

colz <- colnames(datt_gn)[ -c(1,2) ]
datt_gn[, c(colz) := lapply(.SD, function(x){ ifelse(x > 0, 1, 0) }), .SDcols = colz]

colz <- colnames(datt_li)[ -c(1,2) ]
datt_li[, c(colz) := lapply(.SD, function(x){ ifelse(x > 0, 1, 0) }), .SDcols = colz]

colz <- colnames(datt_dc)[ -c(1,2) ]
datt_dc[, c(colz) := lapply(.SD, function(x){ ifelse(x > 0, 1, 0) }), .SDcols = colz]

any(is.na(datt_sp))
any(is.na(datt_gn))
any(is.na(datt_li))
any(is.na(datt_dc))


DATT <- list(
  "Species" = datt_sp,
  "Genus" = datt_gn,
  "Lineage" = datt_li,
  "DecidConif" = datt_dc
  )

## Export data
saveRDS(object = DATT,
  file = "Data_EstLatPhi.RData",
  compress = "xz")

## Start up a local cluster (to run estimation in parallel)
library(doFuture)
registerDoFuture()
# plan(multisession, workers = 4)      # for RStudio
plan(multicore, workers = 4)           # will crash RStudio
options(future.globals.maxSize = 5e9)  # 5GB; default = 500 * 1024 ^ 2 = 500 MiB


## Main workflow function
calc_phi <- function(x, perm = 10000){
  # the second column in x should contain groupping variable

  setDF(x)

  res <- list()

  res$PHI <- try( multipatt(
    x = x[, -c(1,2)],
    cluster = x[,2],
    func = "r.g",
    control = how(nperm = perm),
    duleg = TRUE) )

  res$ASC <- try( strassoc(
    X = x[, -c(1,2)],
    cluster = x[,2],
    func = "r.g",
    group = NULL,
    nboot.ci = perm,
    alpha.ci = 0.05) )

  ## Calculate complementary fidelity values for the tree species (UbinB metric)
  res$III <- try( multipatt(
    x = x[, -c(1,2)],
    cluster = x[,2],
    func = "IndVal.g",
    control = how(nperm = perm),
    duleg = TRUE) )

  return(res)
}
# tst <- calc_phi(DATT$DecidConif)


RES <- llply(.data = DATT, .fun = function(x){
  calc_phi(x)
  }, .parallel = TRUE)


saveRDS(object = RES,
  file = "Data_EstLatPhi_RESULTS.RData",
  compress = "xz")



############### Parse the results


library(indicspecies)
library(data.table)
library(plyr)
library(ggplot2)
library(openxlsx)

set.seed(14789)
theme_set(theme_classic(base_size = 14))

## Load data
DATT <- readRDS("Data_EstLatPhi.RData")

## Load results
RES <- readRDS("Data_EstLatPhi_RESULTS.RData")

## Extract Phi
PHI <- llply(.data = RES, .fun = function(x){
  phi <- melt(x$ASC$stat, value.name = "Phi")
  setDT(phi)
  colnames(phi)[1:2] <- c("OTU", "TreeSpecies")

  ## Reshape long to wide (by phi columns)
  res <- dcast(data = phi, formula = OTU ~ TreeSpecies, fun.aggregate = sum, fill = NA, value.var = "Phi")

  # return(phi)
  return(res)
  })

write.xlsx(PHI,
  file = "Indicspecies_Results__PhiValues_AllHosts.xlsx", colNames = TRUE)


## Extract significant Phi
PHIs <- llply(.data = RES, .fun = function(x){
  
  phi <- metagMisc::dfRowName(x = x$PHI$sign, name = "OTU")
  setDT(phi)

  treesp <- colnames(phi)
  treesp <- treesp[ ! treesp %in% c("OTU", "index", "stat", "p.value") ]

  phi[, BestHost := treesp[ index ] ]
  setcolorder(x = phi, neworder = c("OTU", "BestHost", "index", "stat", "p.value"))

  return(phi)
  })

write.xlsx(PHIs,
  file = "Indicspecies_Results__PhiValues_BestHost.xlsx", colNames = TRUE)



## Extract indvals
IIIs <- llply(.data = RES, .fun = function(x){
  
  iii <- metagMisc::dfRowName(x = x$III$sign, name = "OTU")
  setDT(iii)

  treesp <- colnames(iii)
  treesp <- treesp[ ! treesp %in% c("OTU", "index", "stat", "p.value") ]

  iii[, BestHost := treesp[ index ] ]
  setcolorder(x = iii, neworder = c("OTU", "BestHost", "index", "stat", "p.value"))

  return(iii)
  })

write.xlsx(IIIs,
  file = "Indicspecies_Results__IndValValues_BestHost.xlsx", colNames = TRUE)



## summary(III, alpha = 1, indvalcomp = TRUE)

## Extract summary in tabular form (based on `summary.multipatt`)
get_stats_by_cluster <- function(x,
  alpha = 0.05, minstat = NULL, At = NULL, Bt = NULL, indvalcomp = FALSE, ...) {
  
    ncomb = ncol(x$str)
    ncolsign = ncol(x$sign)
    nsps = nrow(x$sign)

    cat("\n Multilevel pattern analysis")
    cat("\n ---------------------------\n")
    cat("\n Association function:", x$func)
    cat("\n Significance level (alpha):", alpha)
    if (x$func == "IndVal" || x$func == "IndVal.g") {
        if (!is.null(At)) { cat("\n Minimum positive predictive value (At):", At) }
        if (!is.null(Bt)) { cat("\n Minimum sensitivity (Bt):", Bt) }
    }
    cat("\n\n Total number of species:", nsps)

    sel = !is.na(x$sign$p.value) & x$sign$p.value <= alpha
    if (!is.null(minstat)) { sel = sel & (x$sign$stat >= minstat) }
    if (!is.null(Bt) && !is.null(x$B)) {
        for (i in 1:nrow(x$sign)){ sel[i] = sel[i] && (x$B[i, x$sign$index[i]] >= Bt) }
    }
    if (!is.null(At) && !is.null(x$A)) {
        for (i in 1:nrow(x$sign)) { sel[i] = sel[i] && (x$A[i, x$sign$index[i]] >= At) }
    }

    a = x$sign[sel, ]

    cat("\n Selected number of species:", nrow(a), "\n")
    
    cols = (ncolsign - 1):ncolsign
    
    if (indvalcomp && !is.null(x$B) && !is.null(x$A)) {
        As = numeric(nrow(x$sign))
        Bs = numeric(nrow(x$sign))
        for (i in 1:nrow(x$sign)) {
            As[i] = x$A[i, x$sign$index[i]]
            Bs[i] = x$B[i, x$sign$index[i]]
        }
        y = cbind(x$sign, As, Bs)
        cols = c(ncol(y) - 1, ncol(y), cols)
        names(y) = c(names(x$sign), "A", "B")
    }
    else y = x$sign
    
    for (k in 1:(ncolsign - 4)) {
        cat(" Number of species associated to", k, if (k == 1) 
            "group:"
        else "groups:", sum(rowSums(a[, 1:(ncolsign - 3)]) == k), "\n")
    }
    

    res <- list()

    for (i in 1:ncomb) {
        sel = x$sign$index == i & !is.na(x$sign$p.value) & x$sign$p.value <= alpha
        if (!is.null(minstat)) { sel = sel & (x$sign$stat >= minstat) }
        if (!is.null(Bt) && !is.null(x$B)) {
            for (j in 1:nrow(x$sign)) { sel[j] = sel[j] && (x$B[j, x$sign$index[j]] >= Bt) }
        }
        if (!is.null(At) && !is.null(x$A)) {
            for (j in 1:nrow(x$sign)) { sel[j] = sel[j] && (x$A[j, x$sign$index[j]] >= At) }
        }
        m = y[sel, ]
        if (nrow(m) > 0) {
          m <- m[ order(m$stat, decreasing = TRUE), cols ]
          m <- metagMisc::dfRowName(x = m, name = "OTU")
          setDT(m)
          res[[ colnames(x$comb)[i] ]] <- m
        }
    }
  return(res)
}


IIS <- llply(.data = RES,
  .fun = function(x){ 
    get_stats_by_cluster(x$III, alpha = 1, indvalcomp = TRUE)
  })


names(IIS)

write.xlsx(IIS$Species,
  file = "Indicspecies_Results__Fidelity_Species.xlsx", colNames = TRUE)

write.xlsx(IIS$Genus,
  file = "Indicspecies_Results__Fidelity_Genus.xlsx", colNames = TRUE)

write.xlsx(IIS$Lineage,
  file = "Indicspecies_Results__Fidelity_Lineage.xlsx", colNames = TRUE)

write.xlsx(IIS$DecidConif,
  file = "Indicspecies_Results__Fidelity_DecidConif.xlsx", colNames = TRUE)



ggplot(data = phi, aes(x = OTU, y = Phi)) + 
  geom_hline(yintercept=0, color="darkgrey", linetype = "longdash") +
  geom_point(aes(color = TreeSpecies)) +
  facet_wrap(~ TreeSpecies, scales = "free") + 
  coord_flip()





x <- DATT$Species
occ <- data.table(Sp = colnames(x)[-c(1,2)], Occ = colSums(x[, -c(1,2)]))

r <- PHIs$Species

r <- merge(r, occ, by.x = "OTU", by.y = "Sp")
r[, invp := 1/p.value]
r[, sign := p.value < 0.01 ]
r[  grepl(pattern = "TEST", x = OTU), test := "test" ]
r[ !grepl(pattern = "TEST", x = OTU), test := "real" ]
r$test <- factor(r$test, levels = c("real", "test"))

ggplot(data = r, aes(x = Occ, y = stat)) + 
  geom_point(aes(size = invp, color = sign, shape = test)) + 
  facet_wrap(~ BestHost)


r[, BestHost := gsub(pattern = "^s.", replacement = "", x = BestHost)]
r[, BestHost := gsub(pattern = "_", replacement = " ", x = BestHost)]


pp <- ggplot(data = r[ Occ < 8 & p.value < 0.5 ], aes(x = p.value, y = stat)) + 
  geom_vline(xintercept=0.05, color="darkgrey", linetype = "longdash") + 
  geom_point(aes(size = Occ, color = Occ, shape = test), alpha = 0.5) + 
  scale_size(range = c(0.2, 3)) +
  scale_color_distiller(direction = 1) + 
  scale_shape_manual(values = c(16, 18)) +
  facet_wrap(~ BestHost, nrow = 2) + 
  scale_y_continuous(breaks=c(0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.6)) +
  annotate("segment", x=-Inf, xend=Inf, y=-Inf, yend=-Inf) +
  annotate("segment", x=-Inf, xend=-Inf, y=-Inf, yend=Inf) + 
  labs(x = "p-value", y = "Phi", color = "Number of\noccurrences", size = "Number of\noccurrences") +
  theme(strip.background = element_blank())

ggsave(filename = "Phi_threshold_plot.png", 
  plot = pp, width = 24, height = 10)



##############################################
############################################## PhiNM - LT, 12 Jan 2023
##############################################
# PhiNM sheet for generating genus-level Phi-values for each 9 host genera for each of the fungal OTU.
# In a later step, we also need to generate average PHI - both weighted and unweighted values - 
# based on the community sheet. 


library(indicspecies)
library(data.table)
library(plyr)
library(ggplot2)
library(openxlsx)

set.seed(14789)
theme_set(theme_classic(base_size = 14))


### Load the data

datt_nm <- readxl::read_excel(
  path = "Phi_NM.xlsx",
  sheet = "PhiNM",
  na = c("", " ", "NA", "#N/A"))

tmp <- unlist(as.vector(datt_nm[1,]))[-1]
datt_nm_meta <- data.table(SampleID = names(tmp), Genus = tmp)
rm(tmp)

datt_nm <- datt_nm[-1, ]
colnames(datt_nm)[1] <- "OTU"
setDT(datt_nm)
colz_to_fix <- colnames(datt_nm)[-1]
datt_nm[, (colz_to_fix) := lapply(.SD, as.numeric), .SDcols = colz_to_fix]
rm(colz_to_fix)


## Prepare data for indicspecies
NM <- as.data.frame( t(datt_nm[,-1]) )
setDT(NM)
colnames(NM) <- datt_nm[,1]$OTU
colz <- colnames(NM)
NM[, c(colz) := lapply(.SD, function(x){ ifelse(x > 0, 1, 0) }), .SDcols = colz]
NM <- data.table(SampleID = colnames(datt_nm)[-1], NM)
NM$genus <- datt_nm_meta$Genus
setcolorder(x = NM, neworder = c("SampleID", "genus"))


## Main workflow function
calc_phi <- function(x, perm = 10000){
  # the second column in x should contain groupping variable

  setDF(x)

  res <- list()

  res$PHI <- try( multipatt(
    x = x[, -c(1,2)],
    cluster = x[,2],
    func = "r.g",
    control = how(nperm = perm),
    duleg = TRUE) )

  res$ASC <- try( strassoc(
    X = x[, -c(1,2)],
    cluster = x[,2],
    func = "r.g",
    group = NULL,
    nboot.ci = perm,
    alpha.ci = 0.05) )

  ## Calculate complementary fidelity values for the tree species (UbinB metric)
  res$III <- try( multipatt(
    x = x[, -c(1,2)],
    cluster = x[,2],
    func = "IndVal.g",
    control = how(nperm = perm),
    duleg = TRUE) )

  return(res)
}


NM_RES <- calc_phi(NM)

saveRDS(object = NM_RES, file = "PhiNM_results.RData", compress = "xz")



## Extract Phi
ext_phi <- function(x){
  phi <- melt(x$ASC$stat, value.name = "Phi")
  setDT(phi)
  colnames(phi)[1:2] <- c("OTU", "Genus")

  ## Reshape long to wide (by phi columns)
  res <- dcast(data = phi, formula = OTU ~ Genus, fun.aggregate = sum, fill = NA, value.var = "Phi")

  # return(phi)
  return(res)
  }

PHI <- ext_phi(NM_RES)

## Extract significant Phi
ext_phi_sig <- function(x){
  
  phi <- metagMisc::dfRowName(x = x$PHI$sign, name = "OTU")
  setDT(phi)

  treesp <- colnames(phi)
  treesp <- treesp[ ! treesp %in% c("OTU", "index", "stat", "p.value") ]

  phi[, BestHost := treesp[ index ] ]
  setcolorder(x = phi, neworder = c("OTU", "BestHost", "index", "stat", "p.value"))

  return(phi)
  }


PHIs <- ext_phi_sig(NM_RES)


write.xlsx(list(
  "PhiNM_genus" = PHI,
  "PhiNM_BestHost" = PHIs
  ),
  file = "PhiNM_results_Phi-by-Genus.xlsx", colNames = TRUE)



######### average Phi per community

## Load community data
datt_sp <- readxl::read_excel(
  path = "Phi_NM.xlsx",
  sheet = "NM_matrix",
  na = c("", " ", "NA", "#N/A"))

setDT(datt_sp)
datt_sp <- datt_sp[ , -"LTID"]

## Sample metadata (genus level)
meta <- readxl::read_excel(
  path = "EstLatPhi_v2.xlsx",
  sheet = "R_genus9",
  na = c("", " ", "NA", "#N/A"))

## Sample metadata tree hosts proportions
host <- readxl::read_excel(
  path = "host.xlsx",
  na = c("", " ", "NA", "#N/A"))

setDT(meta)
setDT(host)
meta <- meta[, .(SampleID, genus)]

## Find the hosts in mixed-host communities
colnames(host)[1] <- "genus"
smps <- colnames(host)[ ! colnames(host) %in% meta$SampleID ]
mixhosts <- host[ , ..smps ]
mixhosts <- melt(mixhosts, id.vars = "genus", variable.name = "SampleID", value.name = "Proportion")

# mixhosts <- mixhosts[ Proportion >= 0.1 ]   # exclude hosts with <0.1% abundance
mixhosts <- mixhosts[ Proportion > 0 ]
mixhosts <- mixhosts[ , .(genus = paste(genus, collapse = ";")), by = "SampleID" ]

meta <- rbind(meta, mixhosts)
any( duplicated(meta$SampleID) )
meta <- meta[, .(SampleID, genus)]


## Convert to long format, add sample metadata
datl <- melt(data = datt_sp, id.vars = "OTU", variable.name = "SampleID", value.name = "Abundance")
datl <- merge(x = datl, y = meta, by = "SampleID", all.x = TRUE)
datl <- merge(x = datl, y = PHI, by = "OTU", all.x = TRUE)
datl <- merge(x = datl, y = PHIs[, .(OTU, stat)], by = "OTU", all.x = TRUE)
setnames(x = datl, old = "stat", new = "PhiMax")

## Remove zero-abundancence
datl <- datl[ Abundance > 0, ]

## Estimate relative abundances within a sample
datl[ , RelAbundance := Abundance / sum(Abundance), by = "SampleID"]

## Estimate wigthed phi values
datl[ , WPhiMax      := PhiMax  * RelAbundance ]
datl[ , WPhi.Alnus   := Alnus   * RelAbundance ]
datl[ , WPhi.Betula  := Betula  * RelAbundance ]
datl[ , WPhi.Corylus := Corylus * RelAbundance ]
datl[ , WPhi.Picea   := Picea   * RelAbundance ]
datl[ , WPhi.Pinus   := Pinus   * RelAbundance ]
datl[ , WPhi.Populus := Populus * RelAbundance ]
datl[ , WPhi.Quercus := Quercus * RelAbundance ]
datl[ , WPhi.Salix   := Salix   * RelAbundance ]
datl[ , WPhi.Tilia   := Tilia   * RelAbundance ]


## Average max-Phi for a genus (+ Weighted average)
avg <- datl[ , .(
  N = .N,
  MaxPhi_Mean = mean(PhiMax, na.rm = TRUE),
  MaxPhi_MeanW = weighted.mean(x = PhiMax, w = Abundance, na.rm = TRUE),
  MaxPhi_Median = median(PhiMax, na.rm = TRUE),
  MaxPhi_LCI = Hmisc::smean.cl.boot(PhiMax, na.rm = TRUE)[[ 1 ]],
  MaxPhi_UCI = Hmisc::smean.cl.boot(PhiMax, na.rm = TRUE)[[ 3 ]]
  ), by = "genus" ]

write.xlsx(list(
  "maxphi_summary" = avg
  ),
  file = "PhiNM_results_Weighted_PHI_summary_byGenus.xlsx", colNames = TRUE)


## Average phi per community
datll <- melt(data = datl,
  id.vars = c("OTU", "SampleID", "Abundance", "RelAbundance", "genus"),
  measure.vars = c("PhiMax", "Alnus", "Betula", "Corylus", "Picea", "Pinus", "Populus", "Quercus", "Salix", "Tilia"), 
    # "WPhi.Alnus", "WPhi.Betula", "WPhi.Corylus", "WPhi.Picea", "WPhi.Pinus", "WPhi.Populus", "WPhi.Quercus", "WPhi.Salix", "WPhi.Tilia", , "WPhiMax"),
  variable.name = "Variable", value.name = "Value")

datll <- datll[ ! is.na(Value)]

avgc <- datll[ , .(
  N = .N,
  Mean = mean(Value, na.rm = TRUE),
  MeanW = weighted.mean(x = Value, w = Abundance, na.rm = TRUE),
  Median = median(Value, na.rm = TRUE),
  LCI = Hmisc::smean.cl.boot(Value, na.rm = TRUE)[[ 1 ]],
  UCI = Hmisc::smean.cl.boot(Value, na.rm = TRUE)[[ 3 ]]
  ), by = c("SampleID", "genus", "Variable") ]

avgcw <- dcast(data = avgc, )

write.xlsx(list(
  "per_sample" =
  "long_tab" = avgc
  ),
  file = "PhiNM_results_Phi-per-community.xlsx", colNames = TRUE)





## Un-weighted table
write.xlsx(list(
  "PhiMax"      = dcast(data = datl, formula = OTU ~ SampleID, value.var = "PhiMax"),
  "Phi.Alnus"   = dcast(data = datl, formula = OTU ~ SampleID, value.var = "Alnus"),
  "Phi.Betula"  = dcast(data = datl, formula = OTU ~ SampleID, value.var = "Betula"),
  "Phi.Corylus" = dcast(data = datl, formula = OTU ~ SampleID, value.var = "Corylus"),
  "Phi.Picea"   = dcast(data = datl, formula = OTU ~ SampleID, value.var = "Picea"),
  "Phi.Pinus"   = dcast(data = datl, formula = OTU ~ SampleID, value.var = "Pinus"),
  "Phi.Populus" = dcast(data = datl, formula = OTU ~ SampleID, value.var = "Populus"),
  "Phi.Quercus" = dcast(data = datl, formula = OTU ~ SampleID, value.var = "Quercus"),
  "Phi.Salix"   = dcast(data = datl, formula = OTU ~ SampleID, value.var = "Salix"),
  "Phi.Tilia"   = dcast(data = datl, formula = OTU ~ SampleID, value.var = "Tilia")
  ),
  file = "PhiNM_results_PHI_matrix_NoWeigth.xlsx", colNames = TRUE)


## Relative-abundance weighted table
write.xlsx(list(
  "WPhiMax"      = dcast(data = datl, formula = OTU ~ SampleID, value.var = "WPhiMax"),
  "WPhi.Alnus"   = dcast(data = datl, formula = OTU ~ SampleID, value.var = "WPhi.Alnus"),
  "WPhi.Betula"  = dcast(data = datl, formula = OTU ~ SampleID, value.var = "WPhi.Betula"),
  "WPhi.Corylus" = dcast(data = datl, formula = OTU ~ SampleID, value.var = "WPhi.Corylus"),
  "WPhi.Picea"   = dcast(data = datl, formula = OTU ~ SampleID, value.var = "WPhi.Picea"),
  "WPhi.Pinus"   = dcast(data = datl, formula = OTU ~ SampleID, value.var = "WPhi.Pinus"),
  "WPhi.Populus" = dcast(data = datl, formula = OTU ~ SampleID, value.var = "WPhi.Populus"),
  "WPhi.Quercus" = dcast(data = datl, formula = OTU ~ SampleID, value.var = "WPhi.Quercus"),
  "WPhi.Salix"   = dcast(data = datl, formula = OTU ~ SampleID, value.var = "WPhi.Salix"),
  "WPhi.Tilia"   = dcast(data = datl, formula = OTU ~ SampleID, value.var = "WPhi.Tilia")
  ),
  file = "PhiNM_results_PHI_matrix_WeigthByRelAbund.xlsx", colNames = TRUE)

# There could be blank cells in the table - it was not possible to estimate Phi for some OTUs in several cases.


write.xlsx(list(
  "phi_long" = datl
  ),
  file = "PhiNM_results_PHI_LongTable.xlsx", colNames = TRUE)

