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
source("bivariateBeta.R")
source("perm_test_funs.R")

alpha <- 0.001


#' Simulate Power of CI at differentiating bivariate beta nulls and alternates 
#' 
#' This code runs with a doParallel backend. Please set one up prior 
#' to running code. 
#' 
#' @param rho the expected alternate pearson
#' @param shape the shape parameters for the marginals (currently symetric)
#' @param N Sample size for samples from null or alt
#' @param sampleN How many times to resample/simulate
#' @param delta_vector Vector of delta values for rCI 
#' @param propTrue If only interested in power, can set to 1 for all to be Alternates. Otherwise, can investigate FDR as well. 
runPowerBetaNull <- function(rho = 0.6, shape = c(1,10),
                           N = 50, 
                           sampleN = 1000, 
                           delta_vector = seq(0, 1, by = .05), 
                           propTrue=1, req_alpha=0.001,
                           withKCI = FALSE){
  
  samplePars <- optimizeShapeForPearson(rho, shape, shape)

  

  res <- foreach(i = seq_len(sampleN), 
                 .export=c("delta_vector", "nullSigma", "altSigma", "propTrue", "N", "req_alpha", "samplePars")) %dopar% {
    # if(i %% 10 == 0){print(i)}
    suppressMessages(require(MASS, quietly = TRUE, warn.conflicts = FALSE))
    require(wCI, quietly = TRUE, warn.conflicts = FALSE)
    require(rmutil, quietly = TRUE, warn.conflicts = FALSE)
    ci_p <- numeric(length(delta_vector))
    mci_p <- numeric(length(delta_vector))
    pearson_p <- numeric(length(delta_vector))
    kci_p <- numeric(length(delta_vector))
    
    
    spearman_p <- numeric(length(delta_vector))
    
    if(runif(1) < propTrue){
      truth <- rep(1, times=length(delta_vector))
      sample_mat <- rbivariateBeta(n = N, samplePars)
      
      x <- sample_mat[,1]
      y <- sample_mat[,2]
    } else {
      truth <- rep(0, times=length(delta_vector))
      sample_mat <- mvrnorm(n = N, mu = c(0,0), Sigma = nullSigma)
      
      x <- sample_mat[,1]
      y <- sample_mat[,2]
    }
    pearson_p_1 <- pearson.perm.test(x,y,req_alpha = req_alpha)$p.value
    spearman_p_1 <- cor.test(x, y, method="spearman")$p.value
    ci_p_1 <- rCI.perm.test(x, y, 0, req_alpha)$p.value
    if(withKCI){
      kci_p_1 <- kCI.perm.test(x,y, req_alpha)$p.value

    } else {
      kci_p_1 <- NA_real_ 
    }

    for(j in seq_along(delta_vector)){
      sgm <- delta_vector[j]
      if(diff(range(x)) <= sgm || diff(range(y)) <= sgm){
        mci_p[j] <- NA_real_
      } else {
        mci <- rCI.perm.test(x, y, sgm, req_alpha)
        mci_p[j] <- mci$p.value
      }
      
      
      ci_p[j] <- ci_p_1
      kci_p[j] <- kci_p_1
      pearson_p[j] <- pearson_p_1
      spearman_p[j] <- spearman_p_1
      
    }


    res <- cbind(truth, ci_p, mci_p, kci_p, pearson_p, spearman_p)
    colnames(res) <- c("Alternative", "CI_p", "rCI_p", "KCI_p", "Pearson_p", "Spearman_p")
    rownames(res) <- delta_vector
    res
  }
  return(res)
}
