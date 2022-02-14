library(R.utils)

optimize_delta <- function(rep_delta, all_delta){
  # Note: this function makes no assumptions about the distribution of the values of deltas.  
  # Inputs:
  #   replicate delta - a vector of all delta measurements between pairs of replicates
  #   all delta - a vector of delta measurements between all pairs (or at least a random sampling) of data
  # Returns:
  #   
  
  if (sum(is.na(rep_delta)) + sum(is.na(all_delta)) > 0){
    printf("Removing %d replicate and %d all delta NAs or infinities", sum(is.na(rep_delta)) + sum(is.infinite(rep_delta)), sum(is.na(all_delta)) + sum(is.infinite(all_delta)))
  }
  rep_delta <- rep_delta[!is.na(rep_delta) & !is.infinite(rep_delta)]
  all_delta <- all_delta[!is.na(all_delta) & !is.infinite(all_delta)]
  
  k <- 1000
  yvals <- seq(1/(2*k), 1-1/(2*k), 1/k)
  
  allcdf <- sort(abs(all_delta))[seq(length(all_delta)/(2*k), length(all_delta)*(1- 1/(2*k)), length(all_delta)/k)]
  repcdf <- sort(abs(rep_delta))[seq(length(rep_delta)/(2*k), length(rep_delta)*(1- 1/(2*k)), length(rep_delta)/k)]
  
  # Myquants chooses a reasonable range for computing fixed values of the CDFs based on the inputs:
  myquants <- sort(c(allcdf, repcdf))
  
  # This computes a (fast) approximation of the CDF.
  #   allfixed$x contains values of x
  #   allfixed$y contains the fraction of abs(all_delta) exceeding x
  allfixed <- approx(x=allcdf, y=1-yvals, xout=myquants, n=10, yright=0)
  repfixed <- approx(x=repcdf, y=1-yvals, xout=myquants, n=10, yright=0)
  
  # Prop null computes in the neighborhood of 0 (at the k/5 points closest to 0) the ratio of 
  # the fraction of all_deltas in that neighborhood to the fraction of rep_deltas in that neighborhood.  
  # Under the hypothesis that deltas from the alternate distribution are greater than those from the null
  # (i.e. replicates), this is exactly equivalent to Storey's pi_{0}, i.e. the fraction of samples from 
  # the "all delta" distribution that are actually drawn from the replicate distribution in the limit. 
  prop_null <- t(sapply(seq(k/5), function(x) approxExtrap(allfixed$x[x:k], (1-allfixed$y[x:k])/(1-repfixed$y[x:k]), xout=0, n=100)))
  prop_null <- as.numeric(prop_null[,2])
  
  # This is ghetto; use loess
  # p_naught - is the estimate of the total proportion of samples in all_delta that are drawn from the rep_delta distribution
  # psame - the (fuzzy) fraction around each value of allfixed$x that are drawn from the rep_delta distribution
  # prop_valid - 1 - psame, i.e. the fraction around each value of allfixed$x that are drawn from the (unknown) alternate distribution
  # pfalsenegs - a vector representing for each bin of allfixed$x 
  
  p_naught <- median(prop_null[20:70])
  psame <- repfixed$y / allfixed$y
  
  #pfalsenegs <- 1 - (1 - repfixed$y)/(1-allfixed$y)
  pfalsenegs <- 1/(1 - p_naught) * ((1-allfixed$y) - p_naught *(1-repfixed$y))
  prop_valid <- 1 - psame
  
  myf1 <- data.frame(x=repfixed$x, f1=2*(1-pfalsenegs)*prop_valid/((1-pfalsenegs) + prop_valid))
  mymcc <- data.frame(x=repfixed$x, mcc=get_mcc(tp=allfixed$y - p_naught*repfixed$y, tn=p_naught*(1-repfixed$y), 
                                                fp=p_naught*repfixed$y, fn=(1-p_naught)-(allfixed$y - p_naught*repfixed$y)))
  
  return(list(f1=myf1, mcc=mymcc, kernel=data.frame(x=allfixed$x, kernel=1-repfixed$y/allfixed$y), repinvcdf=repfixed, allinvcdf=allfixed))
}


get_mcc <- function(tp, tn, fp, fn){
  return((tp * tn - fp * fn)/sqrt((tp+fp)*(tp+fn)*(tn+fp)*(tn+fn)))
}