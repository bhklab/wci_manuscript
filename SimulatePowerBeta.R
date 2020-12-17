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
                           propTrue=1, req_alpha=0.001){
  
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

    kci_p_1 <- kCI.perm.test(x,y, req_alpha)$p.value

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

registerDoParallel(40)

exprho <- seq(0.32, 0.5, .01)
nsamples_loop <- c(100)
list_mat <- matrix(list(), nrow = length(exprho), ncol=length(nsamples_loop), dimnames = list(exprho, nsamples_loop))
if(file.exists("beta_1_10_dist_power_analysis_withkci.RData")) load("beta_1_10_dist_power_analysis_withkci.RData")

for(nsamples in nsamples_loop){
  for(erho in exprho){
    print(c(erho, nsamples))
    test <- runPowerBetaNull(rho = erho, N = nsamples, sampleN=1000, delta_vector = c(0, 0.13), req_alpha = alpha)


    list_mat[as.character(erho), as.character(nsamples)] <- list(test)
    save(list_mat, file="beta_1_10_dist_power_analysis_withkci.RData")

  }
}


exprho <- seq(0, 0.5, .01)
nsamples_loop <- c(100)
list_mat <- matrix(list(), nrow = length(exprho), ncol=length(nsamples_loop), dimnames = list(exprho, nsamples_loop))

if(file.exists("beta_1_2__4_5_dist_power_analysis_withkci.RData")) load("beta_1_2__4_5_dist_power_analysis_withkcismirnov96
  .RData")


for(nsamples in nsamples_loop){
  for(erho in exprho){
    print(c(erho, nsamples))
    test <- runPowerBetaNull(rho = erho, N = nsamples, shape = c(1.2, 4.5) ,sampleN=1000, delta_vector = c(0, 0.13), req_alpha = alpha)


    list_mat[as.character(erho), as.character(nsamples)] <- list(test)
    save(list_mat, file="beta_1_2__4_5_dist_power_analysis_withkci.RData")

  }
}




load("beta_1_10_dist_power_analysis_withkci.RData")

n100 <- list_mat[,"100"]
# n100 <- lapply(n100, unlist, recursive=FALSE)
n100 <- lapply(n100, abind, along=-1)

contingency_helper <- function(x1, x2){
  
  resTbl <- matrix(0, nrow=2, ncol=2)
  dimnames(resTbl) <- list(c("0", "1"), c("0","1"))
  resTbl["0","0"] <- sum(x1==0 & x2==0)
  resTbl["0","1"] <- sum(x1==0 & x2==1)
  resTbl["1","0"] <- sum(x1==1 & x2==0)
  resTbl["1","1"] <- sum(x1==1 & x2==1)
  return(resTbl)
}

n100_power <- lapply(n100, function(x) {
  x[,,2:6] <- apply(x[,,2:6], c(1,2,3), function(x) (as.numeric(x < alpha)))
  powers <- sapply(dimnames(x)[[2]], function(xx){

      conf_mat_CI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"CI_p"])
      conf_mat_rCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"rCI_p"])
      conf_mat_KCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"KCI_p"])

      conf_mat_pearson <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Pearson_p"])
      conf_mat_spearman <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Spearman_p"])
      # browser()
      pow_CI <- conf_mat_CI["1","1"]/sum(conf_mat_CI["1",])
      pow_rCI <- conf_mat_rCI["1","1"]/sum(conf_mat_rCI["1",])
      pow_KCI <- conf_mat_KCI["1","1"]/sum(conf_mat_KCI["1",])

      pow_pearson <- conf_mat_pearson["1","1"]/sum(conf_mat_pearson["1",])
      pow_spearman <- conf_mat_spearman["1","1"]/sum(conf_mat_spearman["1",])
      
      return(c("CI Power" = pow_CI, "rCI Power" = pow_rCI, "KCI Power" = pow_KCI, "Pearson Power" = pow_pearson, "Spearman Power" = pow_spearman))
    })
  return(powers)
})

n100_power <- abind(n100_power, along = -1)
library(reshape2)

toPlot <- melt(n100_power[,,"0.13"])
colnames(toPlot) <- c("Effect Size", "Method", "Power")
library(ggplot2)

pres_ready <- theme_bw() + 
  theme(axis.title = element_text(size=20), axis.text = element_text(size=20), legend.text = element_text(size=20), title = element_text(size=24),legend.key.height = unit(1.0, 'cm'))

pdf("power_beta_1_10_with_pearson_delta_0_13_withkci.pdf", height = 6, width=9)
ggplot(toPlot, aes(x=`Effect Size`, y=Power, col=Method)) + geom_line(size=1.5) + ggtitle(paste0("Beta(1,10), N=100,\nDelta = 0.13, alpha = ", alpha)) + pres_ready
dev.off()



load("beta_1_2__4_5_dist_power_analysis.RData")


n100 <- list_mat[,"100"]
# n100 <- lapply(n100, unlist, recursive=FALSE)
n100 <- lapply(n100, abind, along=-1)

contingency_helper <- function(x1, x2){
  
  resTbl <- matrix(0, nrow=2, ncol=2)
  dimnames(resTbl) <- list(c("0", "1"), c("0","1"))
  resTbl["0","0"] <- sum(x1==0 & x2==0)
  resTbl["0","1"] <- sum(x1==0 & x2==1)
  resTbl["1","0"] <- sum(x1==1 & x2==0)
  resTbl["1","1"] <- sum(x1==1 & x2==1)
  return(resTbl)
}

n100_power <- lapply(n100, function(x) {
  x[,,2:5] <- apply(x[,,2:5], c(1,2,3), function(x) (as.numeric(x < alpha)))
  powers <- sapply(dimnames(x)[[2]], function(xx){

      conf_mat_CI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"CI_p"])
      conf_mat_rCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"rCI_p"])
      conf_mat_pearson <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Pearson_p"])
      conf_mat_spearman <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Spearman_p"])
      # browser()
      pow_CI <- conf_mat_CI["1","1"]/sum(conf_mat_CI["1",])
      pow_rCI <- conf_mat_rCI["1","1"]/sum(conf_mat_rCI["1",])
      pow_pearson <- conf_mat_pearson["1","1"]/sum(conf_mat_pearson["1",])
      pow_spearman <- conf_mat_spearman["1","1"]/sum(conf_mat_spearman["1",])
      
      return(c("CI_power" = pow_CI, "rCI_power" = pow_rCI, "pearson_power" = pow_pearson, "spearman_power" = pow_spearman))
    })
  return(powers)
})

n100_power <- abind(n100_power, along = -1)
library(reshape2)

toPlot <- melt(n100_power[,,"0.13"])
colnames(toPlot) <- c("Effect Size", "Method", "Power")
library(ggplot2)

pdf("power_beta_1_2__4_5_with_pearson_delta_0_13.pdf", height = 6, width=9)
ggplot(toPlot, aes(x=`Effect Size`, y=Power, col=Method)) + geom_line(size=1.5) + ggtitle(paste0("Beta(1.2, 4.5), N=100,\nDelta = 0.13, alpha = ", alpha)) + pres_ready
dev.off()
 