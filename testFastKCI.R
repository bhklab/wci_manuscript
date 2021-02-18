library(wCI)
source("../../fastCI/R/naiveKCI.R")

library(foreach)

foreach(i = 1:100, .combine = rbind) %do% {

	x <- runif(100)
	y <- runif(100)
	test1 <- wCI:::KCI(as.numeric(x), as.numeric(y), 
                          as.numeric(length(x)), 
                          psymmetric = 1L, pkern1 = -27.5512, pkern2 = 0.08000) 
	test1 <- test1[[1]]/(test1[[1]] + test1[[2]])
	test2 <- naiveKCI(x, y, compute.p = FALSE)[[1]]

	c(test1, test2)


}

