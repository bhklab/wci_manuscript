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

pearson.perm.test <- function(x, y, req_alpha){

    t0 <- cor(x,y, method="pearson")


    x <- scale(x)
    y <- scale(y)
    x <- x/sum(x^2)
    sampleFun <- function(){
      y <- sample(y)
      return(abs(crossprod(x,y)) > abs(t0))
    }

    p.value <- quickstop(sampleFun, req_alpha = req_alpha)

    return(c(list(t0 = t0), p.value))
}


rCI.perm.test <- function(x, y, delta, req_alpha, outties = 1L){
    outties <- as.integer(!outties)
    t0res <- wCI:::newPCI(as.numeric(x), as.numeric(y), 
                          as.numeric(length(x)), 
                          pdeltaX = as.numeric(delta), 
                          pdeltaY = as.numeric(delta), 
                          pxties = as.integer(outties), 
                          pyties = as.integer(outties), plogic = 1L)
    t0 <- t0res[1]/(t0res[1] + t0res[2])

    sampleFun <- function(){
      y <- sample(y)
      tres <- wCI:::newPCI(as.numeric(x), as.numeric(y), 
                          as.numeric(length(x)), 
                          pdeltaX = as.numeric(delta), 
                          pdeltaY = as.numeric(delta), 
                          pxties = as.integer(outties), 
                          pyties = as.integer(outties), plogic = 1L)
      t <- tres[1]/(tres[1] + tres[2])
      return(abs(t-0.5) > abs(t0 - 0.5))
    }

    p.value <- quickstop(sampleFun, req_alpha = req_alpha)

    return(c(list(t0 = t0), p.value))
}

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
runPowerNormalNull <- function(exp_CI = 0.6, 
                           N = 50, 
                           sampleN = 1000, 
                           delta_vector = seq(0, 1, by = .05), 
                           propTrue=1, req_alpha=0.001){
  tau = 2*exp_CI - 1
  rho <- sin((pi*tau)/(2))
  
  altSigma <- c(1,rho,rho,1)
  dim(altSigma) <- c(2,2)
  
  nullSigma <- c(1,0,0,1)
  dim(nullSigma) <- c(2,2)
  

  res <- foreach(i = seq_len(sampleN), 
                 .export=c("delta_vector", "nullSigma", "altSigma", "propTrue", "N", "req_alpha")) %dopar% {
    # if(i %% 10 == 0){print(i)}
    require(MASS)
    require(wCI)
    require(rmutil)
    ci_p <- numeric(length(delta_vector))
    mci_p <- numeric(length(delta_vector))
    pearson_p <- numeric(length(delta_vector))
    
    
    spearman_p <- numeric(length(delta_vector))
    
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
    ci_p_1 <- cor.test(x, y, method="kendal")$p.value
    for(j in seq_along(delta_vector)){
      sgm <- delta_vector[j]

      mci <- paired.concordance.index(predictions = x,  observations = y, delta.pred = sgm, delta.obs = sgm)
      mci_p[j] <- mci$p.value
      
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

expcis <- seq(0.5, 0.7, .01)
nsamples_loop <- c(20, 50)
list_mat <- matrix(list(), nrow = length(expcis), ncol=length(nsamples_loop), dimnames = list(expcis, nsamples_loop))
for(nsamples in nsamples_loop){
  for(eci in expcis){
    test <- runPowerNormalNull(exp_CI = eci, N = nsamples,sampleN=1000, delta_vector = c(0, 0.5, 1), req_alpha = 0.001)
    

    list_mat[as.character(eci), as.character(nsamples)] <- test
 
  }
}
save(list_mat, file="normal_dist_laplace_noise_JobIds_finemap.RData")
res_list_mat <- apply(list_mat, c(1,2), function(x){
  return(list(getJobResult(x)))
})

save(res_list_mat, file="saveres/powers_res_list_pearson_finemap_laplace.RData")

alpha <- 0.0001

n200 <- res_list_mat[,"200"]
n200 <- lapply(n200, unlist, recursive=FALSE)
n200 <- lapply(n200, abind, along=-1)

contingency_helper <- function(x1, x2){
  
  resTbl <- matrix(0, nrow=2, ncol=2)
  dimnames(resTbl) <- list(c("0", "1"), c("0","1"))
  resTbl["0","0"] <- sum(x1==0 & x2==0)
  resTbl["0","1"] <- sum(x1==0 & x2==1)
  resTbl["1","0"] <- sum(x1==1 & x2==0)
  resTbl["1","1"] <- sum(x1==1 & x2==1)
  return(resTbl)
}


n200_power <- lapply(n200, function(x) {
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

n200_power <- abind(n200_power, along = -1)
library(reshape2)

toPlot <- melt(n200_power[,,"1"])
colnames(toPlot) <- c("Effect Size", "Method", "Power")

pres_ready <- theme_classic() + 
  theme(axis.title = element_text(size=24), axis.text = element_text(size=24), legend.text = element_text(size=24), title = element_text(size=28),legend.key.height = unit(1.0, 'cm'))

library(ggplot2)
pdf("power_gaussian_with_pearson.pdf", height = 6, width=9)
ggplot(toPlot, aes(x=`Effect Size`, y=Power, col=Method)) + geom_line(size=1.5) + ggtitle(paste0("N=200, Delta = 1, alpha = ", alpha)) + pres_ready
dev.off()



toPlot <- melt(n200_power[,,"0"])
colnames(toPlot) <- c("Effect Size", "Method", "Power")

library(ggplot2)
pdf("power_gaussian_with_pearson_nonoise.pdf", height = 6, width=9)
ggplot(toPlot, aes(x=`Effect Size`, y=Power, col=Method)) + geom_line(size=1.5) + ggtitle(paste0("N=200, Delta = 0, alpha = ", alpha)) + pres_ready
dev.off()


### Lets look at n=1000

alpha <- 1e-10

n1000 <- res_list_mat[,"1000"]
n1000 <- lapply(n1000, unlist, recursive=FALSE)
n1000 <- lapply(n1000, abind, along=-1)



n1000_power <- lapply(n1000, function(x) {
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

n1000_power <- abind(n1000_power, along = -1)

toPlot <- melt(n1000_power[,,"1"])
colnames(toPlot) <- c("Effect Size", "Method", "Power")

library(ggplot2)

ggplot(toPlot, aes(x=`Effect Size`, y=Power, col=Method)) + geom_line(size=1.5) + ggtitle(paste0("N=1000, Delta = 1, alpha = ", alpha)) + pres_ready

                                                                                  