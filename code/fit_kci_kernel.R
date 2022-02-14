
fit_kci_kernel <- function(all_deltas, rep_deltas, make_plot=FALSE){
  max_delta <- max(abs(all_deltas), abs(rep_deltas))*1.001
  delta_step <- max_delta/1000
  
  params <- list()
  RMSE <- list()

  # This approximates the CDF discretely, even for very large vectors:
  tcdf <- table(c(seq(0, max_delta, delta_step), ceiling((1/delta_step)*abs(all_deltas)) * delta_step)) - 1
  allcdf_int <- 1 - cumsum(tcdf)/sum(tcdf)
  rcdf <- table(c(seq(0, max_delta, delta_step), ceiling((1/delta_step)*abs(rep_deltas)) * delta_step)) - 1
  repcdf_int <- 1 - cumsum(rcdf)/sum(rcdf)
  
  kervals <- data.frame(x=seq(0, max_delta, delta_step), y=1-repcdf_int/(allcdf_int+1e-8))
  
  ## Fit Kernel
  #sigmodel <- nls(y~1/(1 + exp(-b * (x-c))), data=kervals[1:600,], start=list(b=20,c=0.2))
  sigFunOptim <- function(par, data) {
    x <- data[,"x"]
    y <- data[,"y"]
    return(sum((y-par[[1]]/(1 + exp(-par[[2]] * (x-par[[3]]))))^2))
  }
  sigFunPred <- function(par, x) {
    return(par[[1]]/(1 + exp(-par[[2]] * (x-par[[3]]))))
  }


  sigmodel <- optim(list(a=1, b=1, c=0),sigFunOptim,data=kervals)


  params[["sigmodel 3-param"]] <- sigmodel$par
  RMSE[["sigmodel 3-param"]] <- sqrt(sigmodel$value)
  

  sig3y <- sigFunPred(sigmodel$par, kervals$x)
  
  sigFun2Optim <- function(par, data) {
    x <- data[,"x"]
    y <- data[,"y"]
    return(sum((y-1/(1 + exp(-par[[1]] * (x-par[[2]]))))^2))
  }
  sigFun2Pred <- function(par, x) {
    return(1/(1 + exp(-par[[1]] * (x-par[[2]]))))
  }

  sigmodel2 <- optim(list(b=1, c=0),sigFun2Optim,data=kervals)

  params[["sigmodel 2-param"]] <- sigmodel2$par
  sig2y <- sigFun2Pred(sigmodel2$par, kervals$x)
  RMSE[["sigmodel 2-param"]] <- sqrt(sigmodel2$value)


  expFunOptim <- function(par, data) {
    x <- data[,"x"]
    y <- data[,"y"]
    return(sum((y-1 + exp(-par[[1]] * x))^2))
  }
  expFunPred <- function(par, x) {
    return(1 - exp(-par[[1]] * x))
  }


  expmodel <- optim(list(b=1),expFunOptim,data=kervals)
  params[["exp 1-param"]] <- expmodel$par
  exp1 <- expFunPred(expmodel$par, kervals$x)
  RMSE[["exp 1-param"]] <- sqrt(expmodel$value)


  expFun2Optim <- function(par, data) {
    x <- data[,"x"]
    y <- data[,"y"]
    return(sum((y-1 + exp(-par[[1]] * (x - par[[2]])))^2))
  }
  expFun2Pred <- function(par, x) {
    return(1 - exp(-par[[1]] * (x - par[[2]])))
  }

  expmodel2 <- optim(list(b=1, c=0),expFun2Optim,data=kervals)
  params[["exp 2-param"]] <- expmodel2$par
  exp2 <- expFun2Pred(expmodel2$par, kervals$x)
  RMSE[["exp 2-param"]] <- sqrt(expmodel2$value)

  if (make_plot){
    plot(kervals$x, kervals$y, col="black", type="l", lwd=3, xlab="Delta", ylab="Weight", main=" Fitted Kernel comparison")
    lines(kervals$x, sig2y, col="red", lty=2, lwd=2)
    lines(kervals$x, sig3y, col="purple", lty=2, lwd=2)
    lines(kervals$x, exp1, col="blue", lty=2, lwd=2)
    lines(kervals$x, exp2, col="green", lty=2, lwd=2)
    legend(x="bottomright", legend=c("Empirical Weight", "2-param sigmoid", "3-sigmoid", "1-exponential", "2-exponential"), 
           col=c("black", "red", "purple", "blue", "green"), lwd=c(2,2,2,2,2))
  }
  
  return(list(params=params, RMSE=RMSE))
}