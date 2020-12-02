
fit_kci_kernel <- function(all_deltas, rep_deltas, make_plot=FALSE){
  max_delta <- max(abs(all_deltas), abs(rep_deltas))*1.001
  delta_step <- max_delta/1000
  
  params <- list()
  
  # This approximates the CDF discretely, even for very large vectors:
  tcdf <- table(c(seq(0, max_delta, delta_step), ceiling((1/delta_step)*abs(all_deltas)) * delta_step)) - 1
  allcdf_int <- 1 - cumsum(tcdf)/sum(tcdf)
  rcdf <- table(c(seq(0, max_delta, delta_step), ceiling((1/delta_step)*abs(rep_deltas)) * delta_step)) - 1
  repcdf_int <- 1 - cumsum(rcdf)/sum(rcdf)
  
  kervals <- data.frame(x=seq(0, max_delta, delta_step), y=1-repcdf_int/(allcdf_int+1e-8))
  
  ## Fit Kernel
  #sigmodel <- nls(y~1/(1 + exp(-b * (x-c))), data=kervals[1:600,], start=list(b=20,c=0.2))
  sigmodel <- nls(y~a/(1 + exp(-b * (x-c))), data=kervals, start=list(a=1, b=1, c=1), trace=FALSE, control=list(maxiter=500))
  params[["sigmodel 3-param"]] <- coef(sigmodel)
  sig3y <- predict(sigmodel, kervals$x)
  
  sigmodel2 <- nls(y~1/(1 + exp(-b * (x-c))), data=kervals, start=list(b=1, c=1), trace=FALSE, control=list(maxiter=500))
  params[["sigmodel 2-param"]] <- coef(sigmodel2)
  sig2y <- predict(sigmodel2, kervals$x)
  
  expmodel <- nls(y ~ 1 - exp(-b*(x)), data=kervals, start=list(b=20), trace=FALSE, control=list(maxiter=100))
  params[["exp 1-param"]] <- coef(expmodel)
  exp1 <- predict(expmodel, kervals$x)
  
  expmodel2 <- nls(y ~ 1 - exp(-b*(x-c)), data=kervals, start=list(b=20,c=0), trace=FALSE, control=list(maxiter=100))
  params[["exp 2-param"]] <- coef(expmodel2)
  exp2 <- predict(expmodel2, kervals$x)
  
  if (make_plot){
    plot(kervals$x, kervals$y, col="black", type="l", lwd=3, xlab="Delta", ylab="Weight", main=" Fitted Kernel comparison")
    lines(kervals$x, sig2y, col="red", lty=2, lwd=2)
    lines(kervals$x, sig3y, col="purple", lty=2, lwd=2)
    lines(kervals$x, exp1, col="blue", lty=2, lwd=2)
    lines(kervals$x, exp2, col="green", lty=2, lwd=2)
    legend(x="bottomright", legend=c("Empirical Weight", "2-param sigmoid", "3-sigmoid", "1-exponential", "2-exponential"), 
           col=c("black", "red", "purple", "blue", "green"), lwd=c(2,2,2,2,2))
  }
  
  params
}