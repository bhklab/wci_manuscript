#####################################################################
### Code to investigate whether rCI is more robust against loosing power 
### as noise in the measurements increases - bivariate normal case
###
### Petr Smirnov, 2020
#####################################################################

library(wCI)
library(MASS)
library(tictoc)
library(iterators)
library(mvtnorm)
library(abind)
library(doParallel)
source("quickstop.R")
source("perm_test_funs.R")




#' Simulate Power of CI at differentiating bivariate normal nulls and alternates on noised data
#' 
#' This code runs with a doParallel backend. Please set one up prior 
#' to running code. 
#' 
#' @param noise the noise added to X^3
#' @param N Sample size for samples from null or alt
#' @param sampleN How many times to resample/simulate
#' @param delta If missing, then equal to noise added. Otherwise, a numeric value should be provided
#' @param propTrue If only interested in power, can set to 1 for all to be Alternates. Otherwise, can investigate FDR as well. 
runPowerNormalCubic <- function(noise,#exp_CI = 0.6, 
                           Ns = c(100), 
                           sampleN = 1000, 
                           delta, 
                           propTrue=1, req_alpha=0.001){
  # tau = 2*exp_CI - 1
  # rho <- sin((pi*tau)/(2))
  
  if(missing(delta)) {
    delta <- "noise"
  }
  combineFun <- function(...) return(abind(..., along = -1))
  res <- foreach(noise = noise) %:% foreach(N = Ns, .final = combineFun) %:% foreach(i = seq_len(sampleN), .final = combineFun,
                 .export=c("delta_vector", "nullSigma", "altSigma", "propTrue", "N", "req_alpha")) %dopar% {
    # if(i %% 10 == 0){print(i)}
    suppressMessages(require(MASS, quietly = TRUE, warn.conflicts = FALSE))
    require(wCI, quietly = TRUE, warn.conflicts = FALSE)

    if(delta == "noise"){
      delta_vector <- noise
    } else {
      delta_vector <- delta
    }



    ci_p <- numeric(length(delta_vector))
    mci_p <- numeric(length(delta_vector))
    pearson_p <- numeric(length(delta_vector))
    
    
    spearman_p <- numeric(length(delta_vector))
    kendall_p <- numeric(length(delta_vector))

    if(N >=500){
      rCI.perm.test <- rCI.perm.test.large
    }

    if(runif(1) < propTrue){
      truth <- rep(1, times=length(delta_vector))

      x <- rnorm(N)
      y <- x^3+rnorm(N,0,noise)
    } else {
      truth <- rep(0, times=length(delta_vector))
      
      x <- rnorm(N)
      y <- rnorm(N)^3+rnorm(N,0,noise)
    }
    pearson_p_1 <- pearson.perm.test(x,y,req_alpha = req_alpha)$p.value
    spearman_p_1 <- cor.test(x, y, method="spearman")$p.value
    kendall_p_1 <- cor.test(x,y, method="kendall")$p.value
    ci_p_1 <- rCI.perm.test(x, y, 0, req_alpha)$p.value
    for(j in seq_along(delta_vector)){
      sgm <- delta_vector[j]
      if(diff(range(x)) <= sgm || diff(range(y)) <= sgm){
        mci_p[j] <- NA_real_
      } else {
        mci <- rCI.perm.test(x, y, sgm, req_alpha)
        mci_p[j] <- mci$p.value
      }
      
      
      ci_p[j] <- ci_p_1

      pearson_p[j] <- pearson_p_1
      spearman_p[j] <- spearman_p_1
      kendall_p[j] <- kendall_p_1

      
    }


    res <- cbind(truth, ci_p, mci_p, pearson_p, spearman_p, kendall_p)
    colnames(res) <- c("Alternative", "CI_p", "rCI_p", "Pearson_p", "Spearman_p", "Kendall_p")
    rownames(res) <- delta_vector
    res
  }



  return(res)
}
