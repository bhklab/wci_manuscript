
library(wCI)
library(MASS)
library(tictoc)
library(iterators)
library(mvtnorm)
library(abind)
library(doParallel)

source("SimulatePowerBeta.R")

alpha <- 0.001


registerDoParallel(40)



# exprho <- c(0.2, 0.3, 0.4)
# nsamples_loop <- c(100)
# list_mat <- matrix(list(), nrow = length(exprho), ncol=length(nsamples_loop), dimnames = list(exprho, nsamples_loop))
# if(file.exists("beta_1_2__4_5__dist_power_analysis_deltasearch.RData")) load("beta_1_2__4_5__dist_power_analysis_deltasearch.RData")

# for(nsamples in nsamples_loop){
#   for(erho in exprho){
#     print(c(erho, nsamples))
#     test <- runPowerBetaNull(rho = erho, N = nsamples, shape = c(1.2,4.5), sampleN=1000, delta_vector = seq(0,0.3, 0.05), req_alpha = alpha)

#     list_mat[as.character(erho), as.character(nsamples)] <- list(test)
#     save(list_mat, file="beta_1_2__4_5__dist_power_analysis_deltasearch.RData")

#   }
# }



exprho <- c(0.2, 0.3, 0.4)
nsamples_loop <- c(100)
list_mat <- matrix(list(), nrow = length(exprho), ncol=length(nsamples_loop), dimnames = list(exprho, nsamples_loop))
if(file.exists("beta_1_2__4_5__dist_power_analysis_deltafinemap.RData")) load("beta_1_2__4_5__dist_power_analysis_deltafinemap.RData")

for(nsamples in nsamples_loop){
  for(erho in exprho){
    print(c(erho, nsamples))
    test <- runPowerBetaNull(rho = erho, N = nsamples, shape = c(1.2,4.5), sampleN=10000, delta_vector = c(seq(0,0.05, 0.005),  0.1, 0.15, 0.2, 0.25, 0.3), req_alpha = alpha)

    list_mat[as.character(erho), as.character(nsamples)] <- list(test)
    save(list_mat, file="beta_1_2__4_5__dist_power_analysis_deltafinemap.RData")

  }
}





# exprho <- seq(0.32, 0.5, .01)
# nsamples_loop <- c(100)
# list_mat <- matrix(list(), nrow = length(exprho), ncol=length(nsamples_loop), dimnames = list(exprho, nsamples_loop))
# if(file.exists("beta_1_10_dist_power_analysis_withkci.RData")) load("beta_1_10_dist_power_analysis_withkci.RData")

# for(nsamples in nsamples_loop){
#   for(erho in exprho){
#     print(c(erho, nsamples))
#     test <- runPowerBetaNull(rho = erho, N = nsamples, sampleN=1000, delta_vector = c(0, 0.13), req_alpha = alpha)


#     list_mat[as.character(erho), as.character(nsamples)] <- list(test)
#     save(list_mat, file="beta_1_10_dist_power_analysis_withkci.RData")

#   }
# }


# exprho <- seq(0, 0.5, .01)
# nsamples_loop <- c(100)
# list_mat <- matrix(list(), nrow = length(exprho), ncol=length(nsamples_loop), dimnames = list(exprho, nsamples_loop))

# if(file.exists("beta_1_2__4_5_dist_power_analysis_withkci.RData")) load("beta_1_2__4_5_dist_power_analysis_withkcismirnov96
#   .RData")


# for(nsamples in nsamples_loop){
#   for(erho in exprho){
#     print(c(erho, nsamples))
#     test <- runPowerBetaNull(rho = erho, N = nsamples, shape = c(1.2, 4.5) ,sampleN=1000, delta_vector = c(0, 0.13), req_alpha = alpha)


#     list_mat[as.character(erho), as.character(nsamples)] <- list(test)
#     save(list_mat, file="beta_1_2__4_5_dist_power_analysis_withkci.RData")

#   }
# }




# load("beta_1_10_dist_power_analysis_withkci.RData")

# n100 <- list_mat[,"100"]
# # n100 <- lapply(n100, unlist, recursive=FALSE)
# n100 <- lapply(n100, abind, along=-1)

# contingency_helper <- function(x1, x2){
  
#   resTbl <- matrix(0, nrow=2, ncol=2)
#   dimnames(resTbl) <- list(c("0", "1"), c("0","1"))
#   resTbl["0","0"] <- sum(x1==0 & x2==0)
#   resTbl["0","1"] <- sum(x1==0 & x2==1)
#   resTbl["1","0"] <- sum(x1==1 & x2==0)
#   resTbl["1","1"] <- sum(x1==1 & x2==1)
#   return(resTbl)
# }

# n100_power <- lapply(n100, function(x) {
#   x[,,2:6] <- apply(x[,,2:6], c(1,2,3), function(x) (as.numeric(x < alpha)))
#   powers <- sapply(dimnames(x)[[2]], function(xx){

#       conf_mat_CI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"CI_p"])
#       conf_mat_rCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"rCI_p"])
#       conf_mat_KCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"KCI_p"])

#       conf_mat_pearson <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Pearson_p"])
#       conf_mat_spearman <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Spearman_p"])
#       # browser()
#       pow_CI <- conf_mat_CI["1","1"]/sum(conf_mat_CI["1",])
#       pow_rCI <- conf_mat_rCI["1","1"]/sum(conf_mat_rCI["1",])
#       pow_KCI <- conf_mat_KCI["1","1"]/sum(conf_mat_KCI["1",])

#       pow_pearson <- conf_mat_pearson["1","1"]/sum(conf_mat_pearson["1",])
#       pow_spearman <- conf_mat_spearman["1","1"]/sum(conf_mat_spearman["1",])
      
#       return(c("CI Power" = pow_CI, "rCI Power" = pow_rCI, "KCI Power" = pow_KCI, "Pearson Power" = pow_pearson, "Spearman Power" = pow_spearman))
#     })
#   return(powers)
# })

# n100_power <- abind(n100_power, along = -1)
# library(reshape2)

# toPlot <- melt(n100_power[,,"0.13"])
# colnames(toPlot) <- c("Effect Size", "Method", "Power")
# library(ggplot2)

# pres_ready <- theme_bw() + 
#   theme(axis.title = element_text(size=20), axis.text = element_text(size=20), legend.text = element_text(size=20), title = element_text(size=24),legend.key.height = unit(1.0, 'cm'))

# pdf("power_beta_1_10_with_pearson_delta_0_13_withkci.pdf", height = 6, width=9)
# ggplot(toPlot, aes(x=`Effect Size`, y=Power, col=Method)) + geom_line(size=1.5) + ggtitle(paste0("Beta(1,10), N=100,\nDelta = 0.13, alpha = ", alpha)) + pres_ready
# dev.off()



# load("beta_1_2__4_5_dist_power_analysis.RData")


# n100 <- list_mat[,"100"]
# # n100 <- lapply(n100, unlist, recursive=FALSE)
# n100 <- lapply(n100, abind, along=-1)

# contingency_helper <- function(x1, x2){
  
#   resTbl <- matrix(0, nrow=2, ncol=2)
#   dimnames(resTbl) <- list(c("0", "1"), c("0","1"))
#   resTbl["0","0"] <- sum(x1==0 & x2==0)
#   resTbl["0","1"] <- sum(x1==0 & x2==1)
#   resTbl["1","0"] <- sum(x1==1 & x2==0)
#   resTbl["1","1"] <- sum(x1==1 & x2==1)
#   return(resTbl)
# }

# n100_power <- lapply(n100, function(x) {
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

# n100_power <- abind(n100_power, along = -1)
# library(reshape2)

# toPlot <- melt(n100_power[,,"0.13"])
# colnames(toPlot) <- c("Effect Size", "Method", "Power")
# library(ggplot2)

# pdf("power_beta_1_2__4_5_with_pearson_delta_0_13.pdf", height = 6, width=9)
# ggplot(toPlot, aes(x=`Effect Size`, y=Power, col=Method)) + geom_line(size=1.5) + ggtitle(paste0("Beta(1.2, 4.5), N=100,\nDelta = 0.13, alpha = ", alpha)) + pres_ready
# dev.off()
#  
