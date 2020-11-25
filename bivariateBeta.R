library(MCMCpack)


rbivariateBeta <- function(n, shape){
  
  if(length(shape) != 4){
    stop("Bivariate beta requires 4 shape parameters")
  }
  
  smpls <- rdirichlet(n, shape)

  return(cbind(x=smpls[,1] + smpls[,2], y=smpls[,1] + smpls[,3]))
} 


expectedPearsonBiBeta <- function(shape){
  
  a_1_p = shape[1] + shape[2]
  a_p_1 = shape[1] + shape[3]
  a_0_p = shape[4] + shape[3]
  a_p_0 = shape[4] + shape[2]
  
  return((shape[1]*shape[4] - shape[2]*shape[3])/sqrt(a_1_p*a_p_1*a_0_p*a_p_0))
}

optimizeShapeForPearson <- function(target_pearson, xshape, yshape){
  
  if(sum(xshape) != sum(yshape)){
    stop("Sum of marginal shape parameters should be the same for this distribution to work. ")
  }
  b1 <- xshape[1]
  b2 <- xshape[2]
  b3 <- yshape[1]
  rho_optim <- function(a1){
    return(abs(target_pearson - (a1*(b2-b3+a1) - (b1-a1)*(b3-a1))/sqrt(b1*b2*b3*(b1+b2-b3))))
  }
  a1init <- b1/2
  
  res <- optim(par = a1init,fn =  rho_optim, lower = max(0, b3 - b2), upper = min(c(b1,b3)), method="Brent")
  # browser()
  ## constraints are that a1 > 0,  a2 = b1 - a1 >0 a3 = b3 - a1 > 0, a4 = b2-b3 - a1
  
  a1 <- res$par
  a2 <- b1 - a1 
  a3 <- b3 - a1
  a4 <- b2 - a3
  stopifnot(all(c(a1,a2,a3,a4)>0))
  return(c(a1,a2,a3,a4))
}


mixBiVarBeta <- function(N, shape1, shape2, pi1){
  
  samples1 <- rbivariateBeta(N, shape1)
  samples2 <- rbivariateBeta(N, shape2)
  
  pivec <- runif(N) < pi1
  pivec2 <- runif(N) < pi1
  
  samples <- samples1*cbind(pivec,pivec2) + samples2*(!cbind(pivec,pivec2))
  return(samples)
}

