library(wCI)
library(PharmacoGx)
library(ggplot2)
library(data.table)
library(dplyr)
library(rmutil)
library(R.utils)
library(Biobase)

get_sim_mat <- function(ds1, ds2, metric, naive=1){
  print(sprintf("metric = %s", metric))
  common_drugs <- intersect(rownames(ds1), rownames(ds2))
  common_cells <- intersect(colnames(ds1), colnames(ds2))
  
  aac1 <- ds1[common_drugs, common_cells]
  aac2 <- ds2[, common_cells]
  
  print(sprintf("AAC1: %dx%d, AAC2: %dx%d", dim(aac1)[1], dim(aac1)[2], dim(aac2)[1], dim(aac2)[2]))
  
  cellcounts <- get_cellline_isects(aac1, aac2)
  
  if (naive){
    simmat <- matrix(numeric(dim(aac1)[1] * dim(aac2)[1]), nrow=dim(aac1)[1], dimnames=list(rownames(aac1), rownames(aac2)))

    if (metric == "pearson"){
      simmat <- cor(t(aac1), t(aac2), method="pearson", use = "p")
    } else if (metric == "spearman"){
      simmat <- cor(t(aac1), t(aac2), method="spearman", use="p")
    } else if (metric == "CI"){
      for (ii in seq(dim(aac1)[1])){
        if (ii %% 100 == 0){
          print(sprintf("ii = %d/%d", ii, dim(aac1)[1]))
        }
        for (jj in seq(dim(aac2)[1])){
          ix <- complete.cases(aac1[ii,], aac2[jj,])
          if (sum(ix) > 1){
            ct <- paired.concordance.index.new(aac1[ii,ix], aac2[jj,ix], delta.pred=0, delta.obs=0, p_method="Asymptotic", conf_int_method="Asymptotic", outx=TRUE, outy=TRUE)
            simmat[ii,jj] <- ct$cindex
          }
        }
      }
    } else if (metric == "rCI"){
      for (ii in seq(dim(aac1)[1])){
        if (ii %% 100 == 0){
          print(sprintf("ii = %d/%d", ii, dim(aac1)[1]))
        }
        for (jj in seq(dim(aac2)[1])){
          ix <- complete.cases(aac1[ii,], aac2[jj,])
          if (sum(ix) > 1){
            ct <- paired.concordance.index.new(aac1[ii,ix], aac2[jj,ix], delta.pred=0.12, delta.obs=0.12, p_method="Asymptotic", conf_int_method="Asymptotic", outx=TRUE, outy=TRUE)
            simmat[ii,jj] <- ct$cindex
          }
        }
      }
    } else if (metric == "kCI"){
      for (ii in seq(dim(aac1)[1])){
        for (jj in seq(dim(aac2)[1])){
          ix <- complete.cases(aac1[ii,], aac2[jj,])
          if (sum(ix) > 1){
            ct <- wCI:::KCI(as.numeric(aac1[ii,ix]), as.numeric(aac2[jj,ix]), as.numeric(length(aac1[ii,ix])), 
                            psymmetric=1L, pkern1=-27.52480305, pkern2=0.06463092)
            simmat[ii,jj] <- ct[1]/(ct[1] + ct[2])
          }
        }
      }
    }
    return(list(simmat=simmat, cellcounts=cellcounts, metric=metric))
  }
  else { # Subsample for different values of common cell lines
  }
}


get_cellline_isects <- function(ds1, ds2){
  common_drugs <- intersect(rownames(ds1), rownames(ds2))
  common_cells <- intersect(colnames(ds1), colnames(ds2))
  
  aac1 <- ds1[, common_cells]
  aac2 <- ds2[, common_cells]
  
  cellcounts <- matrix(numeric(dim(aac1)[1] * dim(aac2)[1]), nrow=dim(aac1)[1], dimnames=list(rownames(aac1), rownames(aac2)))
  # for loop is gross
  for (ii in seq(dim(aac1)[1])){
    for (jj in seq(dim(aac2)[1])){
      cellcounts[ii,jj] <- sum(complete.cases(aac1[ii,], aac2[jj,]))
    }
  }
  return(cellcounts)
}


drug_recall_dir <- function(psetdir = ".", outdir = "."){
  mypsets <- list.files(psetdir, pattern="rds")
  cnames <- strtrim(mypsets, 5)
  metrics <- c("pearson", "spearman", "CI", "rCI", "kCI")
  
  for (ii in seq_along(mypsets)){
    for (jj in setdiff(seq_along(mypsets), ii)){
      print(sprintf("pset1: %s, pset2: %s", mypsets[ii], mypsets[jj]))
      ds1 <- readRDS(paste(psetdir, mypsets[ii], sep="/"))
      ds2 <- readRDS(paste(psetdir, mypsets[jj], sep="/"))
      
      aac1 <- summarizeSensitivityProfiles(ds1, sensitivity.measure="aac_recomputed")
      aac2 <- summarizeSensitivityProfiles(ds2, sensitivity.measure="aac_recomputed")
      
      res <- lapply(metrics, FUN=function(x) get_sim_mat(aac1, aac2, metric=x, naive=1))
      saveRDS(res, file=paste(outdir, sprintf("%s_%s_drugrecallquery_naive=1.rds", cnames[ii], cnames[jj]), sep="/"))
    }
  }
}

drug_recall_psetdir <- function(pset1, psetdir=".", outdir="."){
  mypsets <- setdiff(list.files(psetdir, pattern="rds"), pset1)
  cname1 <- strtrim(pset1, 5)
  cnames <- strtrim(mypsets, 5)
  metrics <- c("pearson", "spearman", "CI", "rCI", "kCI")
  
  ds1 <- readRDS(paste(psetdir, pset1, sep="/"))

  for (ii in seq_along(mypsets)){
    print(sprintf("pset1: %s, pset2: %s", pset1, mypsets[ii]))
    ds2 <- readRDS(paste(psetdir, mypsets[ii], sep="/"))
      
    aac1 <- summarizeSensitivityProfiles(ds1, sensitivity.measure="aac_recomputed")
    aac2 <- summarizeSensitivityProfiles(ds2, sensitivity.measure="aac_recomputed")
      
    res <- lapply(metrics, FUN=function(x) get_sim_mat(aac1, aac2, metric=x, naive=1))
    saveRDS(res, file=paste(outdir, sprintf("%s_%s_drugrecallquery_naive=1.rds", cname1, cnames[ii]), sep="/"))
  }
}


drug_recall_calc <- function(pset1="", pset2="", dataname="mydata", outdir = "."){
  metrics <- c("pearson", "spearman", "CI", "rCI", "kCI")
  
  ds1 <- readRDS(pset1)
  ds2 <- readRDS(pset2)
  
  aac1 <- summarizeSensitivityProfiles(ds1, sensitivity.measure="aac_recomputed")
  aac2 <- summarizeSensitivityProfiles(ds2, sensitivity.measure="aac_recomputed")
  
  res <- lapply(metrics, FUN=function(x) get_sim_mat(aac1, aac2, metric=x, naive=1))
  saveRDS(res, file=paste(outdir, sprintf("%s_drugrecallquery_naive=1.rds", dataname), sep="/"))
}
