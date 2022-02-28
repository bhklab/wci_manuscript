library(R.utils)
library(wCI)
library(tictoc)
library(ggplot2)
library(parallel)

source("run_fpr_tests.R")

#npts = 100000
covar <- 0
iter <- 200000
perms <- 1
outdir <- "/data/results"
setwd(outdir)

distribution <- "norm"
use_parallel <- TRUE

if (distribution == "beta") {
  label <- "beta_1-2_4-5"
  dsname <- "Beta(1.2, 4.5)"
  myfunc <- rbeta
  globargs <- list(shape1=1.2, shape2=4.5)
  params <- c(globargs, list(iter=iter, perms=perms, xdist="Beta, shape1=1.2, shape2=4.5", ydist="Beta, shape1=1.2, shape2=4.5"))
  npts_vec <- c(25, 40, 50, 100, 200, 300, 400, 600, 800, 1000)
} else if (distribution == "beta_xtreme"){
  label <- "beta_1_10"
  dsname <- "Beta(1, 10)"
  myfunc <- rbeta
  globargs <- list(shape1=1, shape2=10)
  params <- c(globargs, list(iter=iter, perms=perms, xdist="Beta, shape1=1, shape2=10", ydist="Beta, shape1=1, shape2=10"))
  npts_vec <- c(25, 40, 50, 100, 200, 300, 400, 600, 800, 1000)
} else {
  label <- "norm_sd=1"
  dsname <- "Norm(0,1)"
  myfunc <- rnorm
  globargs <- list(mean=0, sd=1)
  params <- c(globargs, list(iter=iter, perms=perms, xdist="Norm, mean=0, sd=1", ydist="Norm, mean=0, sd=1"))
  npts_vec <- c(25, 40, 50, 100, 200, 300, 400, 600, 800, 1000)
}


tic()
if (use_parallel){
  ### Parallel implementation
  res <- mclapply(npts_vec, FUN=function(x) run_fpr_tests(npts=x, iter=iter, myfunc=myfunc, 
                                                          globargs=globargs, params=params, 
                                                          label=label, outdir=outdir), 
                  mc.cores=detectCores())
} else {
  ### Serial implementation
  for (npts in npts_vec) {
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
      paired.concordance.index.new(predictions = x, observations = y, delta.pred=0.1, delta.obs=0.1, outx=TRUE, logic.operator="and", 
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
  }
}
toc()

### Run KCI separately, as it's quite slow:
kci_range <- c(25, 40, 50, 100, 200, 300, 400, 600, 800)
kci_iter <- 100000
kci_tests(npts=kci_range, iter=kci_iter, myfunc=myfunc, globargs=globargs, params=params, 
          label=label, outdir=outdir, parallel=use_parallel)

