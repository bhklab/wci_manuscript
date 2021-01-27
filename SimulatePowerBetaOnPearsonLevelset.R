library(pwr)
alpha <- 0.001

source("SimulatePowerBeta.R")

registerDoParallel(40)

nsamples_loop <- c(100,150,200,250,300)

exprhos <- sapply(nsamples_loop, function(n) {


	pwr.r.test(n, sig.level=alpha, power=0.5)$r


})

list_mat <- list()


if(file.exists("beta_dist_power_analysis_delta_0_1_fixed_pearson_power_0_5_large.RData")) load("beta_dist_power_analysis_delta_0_1_fixed_pearson_power_0_5_large.RData")

for(ii in seq_along(nsamples_loop)){

	rho <- exprhos[ii]
	nsamples <- nsamples_loop[ii]

    print(c(rho, nsamples))



    test <- runPowerBetaNull(rho = rho, N = nsamples,shape = c(1.2,4.5),sampleN=1000, delta_vector = c(0.1), req_alpha = alpha)
    

    list_mat[ii] <- list(test)
    save(list_mat, nsamples_loop, exprhos, file="beta_dist_power_analysis_delta_0_1_fixed_pearson_power_0_5_large.RData")
}
