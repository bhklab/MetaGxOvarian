#' Function to load ovarian cancer expression sets from the Experiment Hub
#'
#' This function returns ovarian cancer datasets from the hub and a vector of patients from the datasets that are most likely duplicates
#' @param remove.duplicates remove patients with a Spearman correlation greater than or equal to 0.98 with other patient expression profiles (default TRUE)
#' @param quantile.cutoff A nueric between 0 and 1 specifying to remove genes with standard deviation below the required quantile (default 0)
#' @param rescale apply centering and scaling to the expression sets (default FALSE)
#' @param min.number.of.genes an integer specifying to remove expression sets with less genes than this number (default 0)
#' @param min.number.of.events an integer specifying how man survival events must be in the dataset to keep the dataset (default 0)
#' @param min.sample.size an integer specifying the minimum number of patients required in an eset (default 0)
#' @param remove.retracted remove datasets from retracted papers (default TRUE, currently just PMID17290060 dataset)
#' @param remove.subsets remove datasets that are a subset of other datasets (defeault TRUE, currently just PMID19318476)
#' @param keep.common.only remove probes not common to all datasets (default FALSE)
#' @param impute.missing remove patients from datasets with missing expression values
#' @return a list with 2 elements. The First element named esets contains the datasets. The second element named duplicates contains
#' a vector with patient IDs for the duplicate patients (those with  Spearman correlation greater than or equal to 0.98 with other patient expression profiles).
#' @export
#' @examples
#'
#' esetsAndDups = loadOvarianEsets()


loadOvarianEsets = function(remove.duplicates = TRUE, quantile.cutoff = 0, rescale = FALSE, min.number.of.genes = 0,
                            min.number.of.events = 0, min.sample.size = 0, remove.retracted = TRUE, remove.subsets = TRUE,
                            keep.common.only = FALSE, impute.missing = FALSE)
{
  ## -----------------------------------------------------------------------------
  ## needed functions
  ## -----------------------------------------------------------------------------
  filterQuantile <- function(object, q){
    if (!identical(q >=0 && q < 1, TRUE))
      stop("require 0 <= q < 1")
    if (!identical(class(object) == "ExpressionSet", TRUE))
      stop("object must be an ExpressionSet")
    gene.sd <- Biobase::esApply(object,1,sd, na.rm=TRUE)
    gene.quantile <- stats::quantile(gene.sd, probs=q)
    actual.makescutoff <- sum(gene.sd < gene.quantile) / length(gene.sd)
    ##make sure the correct number of genes are getting filtered:
    if (abs(q - actual.makescutoff) > 0.01){
      stop("Not scaling this object, likely pre-scaled.")
    }else{
      object <- object[gene.sd > gene.quantile, ]
    }
    return(object)
  }
  ##recursive intersect function
  intersectMany <- function(lst){
    ## Find the intersection of multiple vectors stored as elements of a
    ## list, through a tail-recursive function.
    if (length(lst)==2){
      return(intersect(lst[[1]],lst[[2]]))
    }else{
      return(intersectMany(c(list(intersect(lst[[1]],lst[[2]])),lst[-1:-2])))
    }
  }

  ##Split out non-specific probe sets
  expandProbesets <- function (eset, sep = "///"){
    x <- lapply(Biobase::featureNames(eset), function(x) strsplit(x, sep)[[1]])
    eset <- eset[order(sapply(x, length)), ]
    x <- lapply(Biobase::featureNames(eset), function(x) strsplit(x, sep)[[1]])
    idx <- unlist(sapply(1:length(x), function(i) rep(i, length(x[[i]]))))
    xx <- !duplicated(unlist(x))
    idx <- idx[xx]
    x <- unlist(x)[xx]
    eset <- eset[idx, ]
    Biobase::featureNames(eset) <- x
    eset
  }

  ## -----------------------------------------------------------------------------
  ##load the esets
  ## -----------------------------------------------------------------------------

  hub = ExperimentHub()
  AnnotationHub::possibleDates(hub)
  ovarianData = query(hub, "MetaGxOvarian")
  esets <- list()
  for(i in 1:length(ovarianData))
  {
    esets[[i]] = ovarianData[[names(ovarianData)[i]]]
    names(esets)[i] = ovarianData[i]$title
  }

  ## -----------------------------------------------------------------------------
  ##Explicit removal of samples from specified datasets:
  ## -----------------------------------------------------------------------------
  delim <- ":"   ##This is the delimiter used to specify dataset:sample,

  ## same as used in metagx getbrcadata
  #load("inst\\extdata\\BenDuplicate.rda")
  #source(system.file("extdata", "patientselection.config", package="MetaGxOvarian"))
  load(system.file("extdata", "BenDuplicate.rda", package="MetaGxOvarian"))

  rmix <- duplicates
  ii <- 1
  while (length(rmix) > ii){
    rmix <- rmix [!is.element(names(rmix), rmix[[ii]])]
    ii <- ii+1
  }
  rmix <- unique(unlist(rmix))

  message("Clean up the esets.")
  for (i in 1:length(esets)){
    eset <- esets[[i]]

    ##filter genes with standard deviation below the required quantile
    if(quantile.cutoff > 0 && quantile.cutoff < 1){
      eset <- filterQuantile(eset, q=quantile.cutoff)
    }
    ##rescale to z-scores
    if(rescale == TRUE){
      exprs(eset) <- t(scale(t(exprs(eset))))
    }

    if(remove.duplicates == TRUE){
      keepix <- setdiff(Biobase::sampleNames(eset), rmix)
      Biobase::exprs(eset) <- Biobase::exprs(eset)[, keepix, drop=FALSE]
      Biobase::pData(eset) <- Biobase::pData(eset)[keepix, , drop=FALSE]

    }

    ##include study if it has enough samples and events:
    if (!is.na(min.number.of.events)
        && exists("min.sample.size") && !is.na(min.sample.size)
        && min.number.of.events > 0
        && sum(eset$vital_status == "deceased") < min.number.of.events
        || ncol(eset) < min.sample.size)
    {
      message(paste("excluding",
                    "(min.number.of.events or min.sample.size)"))
      next
    }
    if(nrow(eset) < min.number.of.genes) {
      message(paste("excluding experiment hub dataset",ovarianData[i]$title,"(min.number.of.genes)"))
      next
    }
    if(remove.retracted && length(grep("retracted", experimentData(eset)@other$warnings$warnings)) > 0){
      message(paste("excluding experiment hub dataset",ovarianData[i]$title,"(remove.retracted)"))
      next
    }
    if(remove.subsets && length(grep("subset", experimentData(eset)@other$warnings$warnings)) > 0){
      message(paste("excluding experiment hub dataset",ovarianData[i]$title,"(remove.subsets)"))
      next
    }
    message(paste("including experiment hub dataset",ovarianData[i]$title))
    ##    featureNames(eset) <- make.names(featureNames(eset))  ##should not do this, it is irreversible.
    esets[[i]] <- eset
    rm(eset)
  }

  ##optionally take the intersection of genes common to all platforms:
  if(keep.common.only){
    features.per.dataset <- lapply(esets, Biobase::featureNames)
    intersect.genes <- intersectMany(features.per.dataset)
    esets <- lapply(esets, function(eset){
      eset <- eset[intersect.genes, ]
      return(eset)
    })
  }

  ids.with.missing.data <- which(sapply(esets, function(X)
    sum(!complete.cases(exprs(X))) > 0))
  message(paste("Ids with missing data:", paste(names(ids.with.missing.data),
                                                collapse=", ")))

  if (length(ids.with.missing.data) > 0 && impute.missing) {
    for (i in ids.with.missing.data) {
      exprs(esets[[i]]) = impute.knn(exprs(esets[[i]]))$data
    }
  }

  retList = list(esets, duplicates)
  names(retList) = c("esets", "duplicates")
  return(retList)
}
