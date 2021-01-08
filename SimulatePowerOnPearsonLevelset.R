library(pwr)
alpha <- 0.001

source("SimulatePowerNormal.R")


nsamples_loop <- c(50,100,150,200,250,300)

exprhos <- sapply(nsamples_loop, function(n) {


	pwr.r.test(n, sig.level=alpha, power=0.5)$r


})

list_mat <- matrix(list(), nrow = length(exprhos), ncol=length(nsamples_loop), dimnames = list(exprhos, nsamples_loop))


if(file.exists("normal_dist_power_analysis_delta_1_0_fixed_pearson_power_0_5.RData")) load("normal_dist_power_analysis_delta_1_0_fixed_pearson_power_0_5.RData")

for(nsamples in nsamples_loop[1]){
  for(rho in exprhos){
    print(c(rho, nsamples))

    test <- runPowerNormalNull(rho = rho, N = nsamples,sampleN=10000, delta_vector = c(1), req_alpha = alpha)
    

    list_mat[as.character(rho), as.character(nsamples)] <- list(test)
    save(list_mat, file="normal_dist_power_analysis_delta_1_0_fixed_pearson_power_0_5.RData")
  }
}
