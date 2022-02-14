library(wCI)
library(MASS)
library(tictoc)
library(iterators)
library(mvtnorm)
library(abind)
library(doParallel)
source("quickstop.R")
source("perm_test_funs.R")
source("SimulatePowerNormal.R4")
alpha <- 0.001

registerDoParallel(20)

exprhos <- seq(0,0.5,0.01)
nsamples_loop <- c(100)
list_mat <- matrix(list(), nrow = length(exprhos), ncol=length(nsamples_loop), dimnames = list(exprhos, nsamples_loop))


if(file.exists("normal_dist_power_analysis_n_100_withkendall.RData")) load("normal_dist_power_analysis_n_100_withkendall.RData")

for(nsamples in nsamples_loop){
  for(rho in exprhos){
    print(c(rho, nsamples))

    test <- runPowerNormalNull(rho = rho, N = nsamples, sampleN=1000, delta_vector = c(1), req_alpha = alpha)
    

    list_mat[as.character(rho), as.character(nsamples)] <- list(test)
    save(list_mat, file="normal_dist_power_analysis_n_100_withkendall.RData")
  }
}


# exprhos <- seq(0.0, 0.5, .01)

#nsamples_loop <- c(100)
#list_mat <- matrix(list(), nrow = length(exprhos), ncol=length(nsamples_loop), dimnames = list(exprhos, nsamples_loop))


#if(file.exists("normal_dist_power_analysis_rho_n_100_deltasearch.RData")) load("normal_dist_power_analysis_rho_n_100_deltasearch.RData")

#for(nsamples in nsamples_loop[1]){
#  for(rho in exprhos){
#    print(c(rho, nsamples))
#    test <- runPowerNormalNull(rho = rho, N = nsamples,sampleN=1000, delta_vector = c(0, 0.5, 1, 1.5, 2), req_alpha = alpha)
    

#    list_mat[as.character(rho), as.character(nsamples)] <- list(test)
#    save(list_mat, file="normal_dist_power_analysis_rho_n_100_deltasearch.RData")

#  }
#}

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