library(wCI)
library(MASS)
library(tictoc)
library(iterators)
library(mvtnorm)
library(abind)
library(doParallel)

source("SimulatePowerBetaWithNoise.R")

alpha <- 0.001


registerDoParallel(40)

exprho <- seq(0.0, 0.75, .01)
nsamples_loop <- c(100)
list_mat <- matrix(list(), nrow = length(exprho), ncol=length(nsamples_loop), dimnames = list(exprho, nsamples_loop))
#if(file.exists("beta_1_10_dist_power_analysis_withkci_withnoise.RData")) load("beta_1_10_dist_power_analysis_withkci_withnoise.RData")

for(nsamples in nsamples_loop){
  for(erho in exprho){
    print(c(erho, nsamples))
    test <- runPowerBetaNullNoise(rho = erho, N = nsamples, sampleN=1000, delta_vector = c(0, 0.1), req_alpha = alpha)


    list_mat[as.character(erho), as.character(nsamples)] <- list(test)
    save(list_mat, file="beta_1_10_dist_power_analysis_withkci_withnoise.RData")

  }
}


exprho <- seq(0.0, 0.6, .01)
nsamples_loop <- c(100)
list_mat <- matrix(list(), nrow = length(exprho), ncol=length(nsamples_loop), dimnames = list(exprho, nsamples_loop))

#if(file.exists("beta_1_2__4_5_dist_power_analysis_withkci_withnoise.RData")) load("beta_1_2__4_5_dist_power_analysis_withkci_withnoise.RData")


for(nsamples in nsamples_loop){
  for(erho in exprho){
    print(c(erho, nsamples))
    test <- runPowerBetaNullNoise(rho = erho, N = nsamples, shape = c(1.2, 4.5) ,sampleN=1000, delta_vector = c(0, 0.1), req_alpha = alpha)


    list_mat[as.character(erho), as.character(nsamples)] <- list(test)
    save(list_mat, file="beta_1_2__4_5_dist_power_analysis_withkci_withnoise.RData")

  }
}

