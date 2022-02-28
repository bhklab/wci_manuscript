library(parallel)

library(fastCI)
library(wCI)
library(data.table)

run_fpr_tests <- function(npts, iter, myfunc, globargs, params, label, outdir){
  printf("npts = %d \n", npts)
  myargs <- c(globargs, n=npts)
  
  pearson_vals <- t(sapply(seq(iter), FUN=function(x) {
    x <- do.call(myfunc, myargs)
    y <- do.call(myfunc, myargs)
    cor.test(x, y, method="pearson")[c("estimate", "p.value")]
  }))
  spearman_vals <- t(sapply(seq(iter), FUN=function(x) {
    x <- do.call(myfunc, myargs)
    y <- do.call(myfunc, myargs)
    cor.test(x, y, method="spearman")[c("estimate", "p.value")]
  }))
  kendall_vals <- t(sapply(seq(iter), FUN=function(x) {
    x <- do.call(myfunc, myargs)
    y <- do.call(myfunc, myargs)
    cor.test(x, y, method="kendall")[c("estimate", "p.value")]
  }))
  ci_vals <- t(sapply(seq(iter), FUN=function(x) {
    x <- do.call(myfunc, myargs)
    y <- do.call(myfunc, myargs)
    paired.concordance.index.new(predictions = x, observations = y, delta.pred=0, delta.obs=0, outx=TRUE, logic.operator="and", 
                                 CPP=TRUE, p_method="Asymptotic", conf_int_method="Asymptotic")[c("cindex", "p.value", "relevant.pairs.no")]
  }))
  rci_vals_20 <- t(sapply(seq(iter), FUN=function(x) {
    x <- do.call(myfunc, myargs)
    y <- do.call(myfunc, myargs)
    paired.concordance.index.new(predictions = x, observations = y, delta.pred=0.2, delta.obs=0.2, outx=TRUE, logic.operator="and", 
                                 CPP=TRUE, p_method="Asymptotic", conf_int_method="Asymptotic")[c("cindex", "p.value", "relevant.pairs.no")]
  }))
  rci_vals_10 <- t(sapply(seq(iter), FUN=function(x) {
    x <- do.call(myfunc, myargs)
    y <- do.call(myfunc, myargs)
    paired.concordance.index.new(predictions = x, observations = y, delta.pred=0.1, delta.obs=0.1, outx=TRUE, logic.operator="and", 
                                 CPP=TRUE, p_method="Asymptotic", conf_int_method="Asymptotic")[c("cindex", "p.value", "relevant.pairs.no")]
  }))
  rci_vals_05 <- t(sapply(seq(iter), FUN=function(x) {
    x <- do.call(myfunc, myargs)
    y <- do.call(myfunc, myargs)
    paired.concordance.index.new(predictions = x, observations = y, delta.pred=0.05, delta.obs=0.05, outx=TRUE, logic.operator="and", 
                                 CPP=TRUE, p_method="Asymptotic", conf_int_method="Asymptotic")[c("cindex", "p.value", "relevant.pairs.no")]
  }))
  
  # The output of the above is a matrix of lists.  This solves that, and converts them to numeric:
  pearson_vals <- as.data.frame(matrix(sapply(pearson_vals, unlist), ncol=2, dimnames=list(rownames(pearson_vals), colnames(pearson_vals))))
  spearman_vals <- as.data.frame(matrix(sapply(spearman_vals, unlist), ncol=2, dimnames=list(rownames(spearman_vals), colnames(spearman_vals))))
  kendall_vals <- as.data.frame(matrix(sapply(kendall_vals, unlist), ncol=2, dimnames=list(rownames(kendall_vals), colnames(kendall_vals))))
  ci_vals <- as.data.frame(matrix(sapply(ci_vals, unlist), ncol=3, dimnames=list(rownames(ci_vals), colnames(ci_vals))))
  rci_vals_20 <- as.data.frame(matrix(sapply(rci_vals_20, unlist), ncol=3, dimnames=list(rownames(rci_vals_20), colnames(rci_vals_20))))
  rci_vals_10 <- as.data.frame(matrix(sapply(rci_vals_10, unlist), ncol=3, dimnames=list(rownames(rci_vals_10), colnames(rci_vals_10))))
  rci_vals_05 <- as.data.frame(matrix(sapply(rci_vals_05, unlist), ncol=3, dimnames=list(rownames(rci_vals_05), colnames(rci_vals_05))))
  
  fpr_vals <- list(pearson_vals=pearson_vals, 
                   spearman_vals=spearman_vals, 
                   kendall_vals=kendall_vals,
                   ci_vals=ci_vals, 
                   rci_vals_20=rci_vals_20, 
                   rci_vals_10=rci_vals_10, 
                   rci_vals_05=rci_vals_05, 
                   params=params)
  
  saveRDS(fpr_vals, file=paste(outdir, sprintf("%s_n=%dk_allmetrics_rcivariants_npts=%d.RDS", label, iter/1000, npts), sep="/"))
  return(npts)
}


kci_tests <- function(npts, iter, myfunc, globargs, params, label, outdir, parallel=0){
  if (parallel == 0){
    res <- list(npts=npts)
    for (ii in npts){
      printf("npts = %d \n", ii)
      myargs <- c(globargs, n=ii)
      kci_vals <- t(sapply(seq(iter), FUN=function(x) {
        x <- do.call(myfunc, myargs)
        y <- do.call(myfunc, myargs)
        naiveKCI(x, y, p.method="Asymptotic")[c("cindex", "p.value", "relevant.pairs.no")]
      }))
    
      kci_vals <- as.data.frame(matrix(sapply(kci_vals, unlist), ncol=3, dimnames=list(rownames(kci_vals), colnames(kci_vals))))
      
      res <- c(res, list(kci_vals))
    }
    names(res) <- c("npts", sprintf("kci_npts_%d", npts))
    saveRDS(res, file=paste(outdir, sprintf("%s_n=%dk_kci_n=%d_%d.RDS", label, iter/1000, min(npts), max(npts)), sep="/"))
  }
  if (parallel == 1){
    res <- list(npts=npts)
    for (ii in npts){
      printf("npts = %d \n", ii)
      myargs <- c(globargs, n=ii)
      kci_vals <- mclapply(seq(iter), FUN=function(z) {
        x <- do.call(myfunc, myargs)
        y <- do.call(myfunc, myargs)
        as.data.frame(naiveKCI(x,y, p.method="Asymptotic")[c("cindex", "p.value", "relevant.pairs.no")])
        }, 
        mc.cores=detectCores())
      
      kci_vals <- as.data.frame(rbindlist(kci_vals))
      #kci_vals <- as.data.frame(matrix(sapply(kci_vals, unlist), ncol=3, dimnames=list(rownames(kci_vals), colnames(kci_vals))))
      
      res <- c(res, list(kci_vals))
    }
    names(res) <- c("npts", sprintf("kci_npts_%d", npts))
    saveRDS(res, file=paste(outdir, sprintf("%s_n=%dk_kci_n=%d_%d.RDS", label, iter/1000, min(npts), max(npts)), sep="/"))
  }
}