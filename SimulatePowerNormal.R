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

alpha <- 0.001



#' Simulate Power of CI at differentiating bivariate normal nulls and alternates on noised data
#' 
#' This code runs with a doParallel backend. Please set one up prior 
#' to running code. 
#' 
#' @param exp_CI the expected alternate CI
#' @param N Sample size for samples from null or alt
#' @param sampleN How many times to resample/simulate
#' @param delta_vector Vector of delta values for rCI 
#' @param propTrue If only interested in power, can set to 1 for all to be Alternates. Otherwise, can investigate FDR as well. 
runPowerNormalNull <- function(rho,#exp_CI = 0.6, 
                           N = 50, 
                           sampleN = 1000, 
                           delta_vector = seq(0, 1, by = .05), 
                           propTrue=1, req_alpha=0.001){
  # tau = 2*exp_CI - 1
  # rho <- sin((pi*tau)/(2))
  
  altSigma <- c(1,rho,rho,1)
  dim(altSigma) <- c(2,2)
  
  nullSigma <- c(1,0,0,1)
  dim(nullSigma) <- c(2,2)
  

  res <- foreach(i = seq_len(sampleN), 
                 .export=c("delta_vector", "nullSigma", "altSigma", "propTrue", "N", "req_alpha")) %dopar% {
    # if(i %% 10 == 0){print(i)}
    suppressMessages(require(MASS, quietly = TRUE, warn.conflicts = FALSE))
    require(wCI, quietly = TRUE, warn.conflicts = FALSE)
    require(rmutil, quietly = TRUE, warn.conflicts = FALSE)
    ci_p <- numeric(length(delta_vector))
    mci_p <- numeric(length(delta_vector))
    pearson_p <- numeric(length(delta_vector))
    
    
    spearman_p <- numeric(length(delta_vector))
    
    if(N >=500){
      rCI.perm.test <- rCI.perm.test.large
    }

    if(runif(1) < propTrue){
      truth <- rep(1, times=length(delta_vector))
      sample_mat <- mvrnorm(n = N, mu = c(0,0), Sigma = altSigma)
      
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
      
    }


    res <- cbind(truth, ci_p, mci_p, pearson_p, spearman_p)
    colnames(res) <- c("Alternative", "CI_p", "rCI_p", "Pearson_p", "Spearman_p")
    rownames(res) <- delta_vector
    res
  }
  return(res)
}

registerDoParallel(40)

exprhos <- seq(0.0, 0.5, .01)
nsamples_loop <- c(50, 100)
list_mat <- matrix(list(), nrow = length(exprhos), ncol=length(nsamples_loop), dimnames = list(exprhos, nsamples_loop))


if(file.exists("normal_dist_power_analysis_rho_n_50_deltasearch.RData")) load("normal_dist_power_analysis_rho_n_50_deltasearch.RData")

for(nsamples in nsamples_loop[1]){
  for(rho in exprhos){
    print(c(rho, nsamples))
    test <- runPowerNormalNull(rho = rho, N = nsamples,sampleN=1000, delta_vector = c(0, 0.5, 1, 1.5, 2, 2.5), req_alpha = alpha)
    

    list_mat[as.character(rho), as.character(nsamples)] <- list(test)
    save(list_mat, file="normal_dist_power_analysis_rho_n_50_deltasearch.RData")

  }
}

# exprhos <- seq(0.7, 0.9, .01)
# nsamples_loop <- c(20, 50)
# list_mat_extension <- matrix(list(), nrow = length(exprhos), ncol=length(nsamples_loop), dimnames = list(exprhos, nsamples_loop))
# for(nsamples in nsamples_loop){
#   for(rho in exprhos){
#     print(c(rho, nsamples))
#     test <- runPowerNormalNull(rho = rho, N = nsamples,sampleN=1000, delta_vector = c(0, 0.5, 1), req_alpha = alpha)
    

#     list_mat_extension[as.character(rho), as.character(nsamples)] <- list(test)
#     save(list_mat_extension, file="normal_dist_power_analysis_extension_rho.RData")

#   }
# }

# exprhos <- seq(0.5, 0.7, .01)
# nsamples_loop <- c(100, 200)
# list_mat <- matrix(list(), nrow = length(exprhos), ncol=length(nsamples_loop), dimnames = list(exprhos, nsamples_loop))
# for(nsamples in nsamples_loop){
#   for(rho in exprhos){
#     print(c(rho, nsamples))
#     test <- runPowerNormalNull(rho = rho, N = nsamples,sampleN=1000, delta_vector = c(0, 0.5, 1), req_alpha = alpha)
    

#     list_mat[as.character(rho), as.character(nsamples)] <- list(test)
#    save(list_mat, file="normal_dist_power_analysis_larger_rho.RData")

#   }
# }
# load("normal_dist_power_analysis.RData")
# load("normal_dist_power_analysis_extension.RData")




# n50 <- c(list_mat[,"50"], list_mat_extension[,"50"])
# # n50 <- lapply(n50, unlist, recursive=FALSE)
# n50 <- lapply(n50, abind, along=-1)

# contingency_helper <- function(x1, x2){
  
#   resTbl <- matrix(0, nrow=2, ncol=2)
#   dimnames(resTbl) <- list(c("0", "1"), c("0","1"))
#   resTbl["0","0"] <- sum(x1==0 & x2==0)
#   resTbl["0","1"] <- sum(x1==0 & x2==1)
#   resTbl["1","0"] <- sum(x1==1 & x2==0)
#   resTbl["1","1"] <- sum(x1==1 & x2==1)
#   return(resTbl)
# }

# n50_power <- lapply(n50, function(x) {
#   x[,,2:5] <- apply(x[,,2:5], c(1,2,3), function(x) (as.numeric(x < alpha)))
#   powers <- sapply(dimnames(x)[[2]], function(xx){

#       conf_mat_CI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"CI_p"])
#       conf_mat_rCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"rCI_p"])
#       conf_mat_pearson <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Pearson_p"])
#       conf_mat_spearman <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Spearman_p"])
#       # browser()
#       pow_CI <- conf_mat_CI["1","1"]/sum(conf_mat_CI["1",])
#       pow_rCI <- conf_mat_rCI["1","1"]/sum(conf_mat_rCI["1",])
#       pow_pearson <- conf_mat_pearson["1","1"]/sum(conf_mat_pearson["1",])
#       pow_spearman <- conf_mat_spearman["1","1"]/sum(conf_mat_spearman["1",])
      
#       return(c("CI_power" = pow_CI, "rCI_power" = pow_rCI, "pearson_power" = pow_pearson, "spearman_power" = pow_spearman))
#     })
#   return(powers)
# })

# n50_power <- abind(n50_power, along = -1)
# library(reshape2)

# toPlot <- melt(n50_power[,,"1"])
# colnames(toPlot) <- c("Effect Size", "Method", "Power")
# library(ggplot2)

# pres_ready <- theme_bw() + 
#   theme(axis.title = element_text(size=24), axis.text = element_text(size=24), legend.text = element_text(size=24), title = element_text(size=28),legend.key.height = unit(1.0, 'cm'))

# pdf("power_gaussian_with_pearson_delta_1.pdf", height = 6, width=9)
# ggplot(toPlot, aes(x=`Effect Size`, y=Power, col=Method)) + geom_line(size=1) + ggtitle(paste0("N=50, Delta = 1, alpha = ", alpha)) + pres_ready
# dev.off()



# toPlot <- melt(n50_power[,,"0"])
# colnames(toPlot) <- c("Effect Size", "Method", "Power")

# library(ggplot2)
# pdf("power_gaussian_with_pearson_delta_0.pdf", height = 6, width=9)
# ggplot(toPlot, aes(x=`Effect Size`, y=Power, col=Method)) + geom_line(size=1.5) + ggtitle(paste0("N=50, Delta = 0, alpha = ", alpha)) + pres_ready
# dev.off()


# n20 <- list_mat[,"20"]
# # n20 <- lapply(n20, unlist, recursive=FALSE)
# n20 <- lapply(n20, abind, along=-1)

# contingency_helper <- function(x1, x2){
  
#   resTbl <- matrix(0, nrow=2, ncol=2)
#   dimnames(resTbl) <- list(c("0", "1"), c("0","1"))
#   resTbl["0","0"] <- sum(x1==0 & x2==0)
#   resTbl["0","1"] <- sum(x1==0 & x2==1)
#   resTbl["1","0"] <- sum(x1==1 & x2==0)
#   resTbl["1","1"] <- sum(x1==1 & x2==1)
#   return(resTbl)
# }

# n20_power <- lapply(n20, function(x) {
#   x[,,2:5] <- apply(x[,,2:5], c(1,2,3), function(x) (as.numeric(x < alpha)))
#   powers <- sapply(dimnames(x)[[2]], function(xx){

#       conf_mat_CI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"CI_p"])
#       conf_mat_rCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"rCI_p"])
#       conf_mat_pearson <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Pearson_p"])
#       conf_mat_spearman <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Spearman_p"])
#       # browser()
#       pow_CI <- conf_mat_CI["1","1"]/sum(conf_mat_CI["1",])
#       pow_rCI <- conf_mat_rCI["1","1"]/sum(conf_mat_rCI["1",])
#       pow_pearson <- conf_mat_pearson["1","1"]/sum(conf_mat_pearson["1",])
#       pow_spearman <- conf_mat_spearman["1","1"]/sum(conf_mat_spearman["1",])
      
#       return(c("CI_power" = pow_CI, "rCI_power" = pow_rCI, "pearson_power" = pow_pearson, "spearman_power" = pow_spearman))
#     })
#   return(powers)
# })

# n20_power <- abind(n20_power, along = -1)
# library(reshape2)

# toPlot <- melt(n20_power[,,"1"])
# colnames(toPlot) <- c("Effect Size", "Method", "Power")
# library(ggplot2)

# pres_ready <- theme_classic() + 
#   theme(axis.title = element_text(size=24), axis.text = element_text(size=24), legend.text = element_text(size=24), title = element_text(size=28),legend.key.height = unit(1.0, 'cm'))

# pdf("power_gaussian_with_pearson_delta_1_n20.pdf", height = 6, width=9)
# ggplot(toPlot, aes(x=`Effect Size`, y=Power, col=Method)) + geom_line(size=1.5) + ggtitle(paste0("N=20, Delta = 1, alpha = ", alpha)) + pres_ready
# dev.off()



# toPlot <- melt(n20_power[,,"0"])
# colnames(toPlot) <- c("Effect Size", "Method", "Power")

# library(ggplot2)
# pdf("power_gaussian_with_pearson_delta_0.pdf", height = 6, width=9)
# ggplot(toPlot, aes(x=`Effect Size`, y=Power, col=Method)) + geom_line(size=1.5) + ggtitle(paste0("N=200, Delta = 0, alpha = ", alpha)) + pres_ready
# dev.off()

#                                                                                   