library(pwr)
alpha <- 0.001

source("SimulatePowerNormal.R")

registerDoParallel(80)

nsamples_loop <- c(50,100,150,200,250,300)

exprhos <- sapply(nsamples_loop, function(n) {


	pwr.r.test(n, sig.level=alpha, power=0.5)$r


})

list_mat <- list()


if(file.exists("normal_dist_power_analysis_delta_1_0_fixed_pearson_power_0_5_large.RData")) load("normal_dist_power_analysis_delta_1_0_fixed_pearson_power_0_5_large.RData")

for(ii in seq_along(nsamples_loop)){

	rho <- exprhos[ii]
	nsamples <- nsamples_loop[ii]

    print(c(rho, nsamples))



    test <- runPowerNormalNull(rho = rho, N = nsamples,sampleN=10000, delta_vector = c(1), req_alpha = alpha)
    

    list_mat[ii] <- list(test)
    save(list_mat, nsamples_loop, exprhos, file="normal_dist_power_analysis_delta_1_0_fixed_pearson_power_0_5_large.RData")
}
