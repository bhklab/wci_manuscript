### Code to calculate the expected number of inversions after applying K adjecent swaps
# library(fastCI)

#' Compute alpha variable for each summand
#' 
#' @param k index of summand
#' @param m total size of list
alpha <- function(k, m){
  return((2*k+1)*pi/(2*m+2))
}

#' Compute the expected number of inversions for N swaps
#' 
#' See paper: http://arxiv.org/abs/0909.0103
#' 
#' @param n Number of swaps
#' @param m Sample/list size
exp_inv <- function(n, m){
  m <- m - 1
  k <- seq(from=0, to=m)
  c <- cos(alpha(k,m))
  s <- sin(alpha(k,m))
  
  run_sum <- 0
  
  for(k in seq_along(c)){
    for(j in seq_along(s)){
      run_sum <- run_sum + 1/(8*(m+1)^2)*(c[k] + c[j])^2/(s[k]^2*s[j]^2)*(1-4/m*(1-c[k]*c[j]))^n
    }
  }
  
  invs <- m*(m+1)/4 - run_sum
  return(invs)
}


#' Compute required K
#' 
#' Computing the required number of swaps to achieve desired
#' number of inversions by binary searching the space between 0 and M^3. 
#' M^3 is theoretically the mixing time of this Markov chain. 
#' See paper: http://arxiv.org/abs/0909.0103
#' 
#' @param ninv Desired number of inversions 
#' @param M Sample size

computeRequiredK <- function(ninv, M){
  
  max_k <- M^3
  min_k <- 0
  if(exp_inv(max_k, M) < ninv){
    return(max_k)
  }
  # zero does not work as condition, may never be achieved
  while(max_k - min_k > 1){ 
    midpoint <- ceiling((max_k + min_k )/2)
    E_inv <- exp_inv(midpoint, M)
    if(E_inv < ninv){
      min_k = midpoint
    } else {
      max_k = midpoint
    }
  }
  E_max <- exp_inv(max_k, M)
  E_min <- exp_inv(min_k, M)
  if(abs(E_max - ninv) < abs(E_min - ninv)){
    return(max_k)
  } else {
    return(min_k)
  }
}

# Code below is some interactive tests, not run
if(FALSE){
  # N is number of swaps
  # M is length of list
  
  M <- 100
  N <- 20434
  
  nperm <- 100
  invs <- numeric(nperm)
  
  for(jjj in seq_len(nperm)){
    cat("\r", jjj)
    x <- as.numeric(seq_len(M))
    
    for(ii in seq_len(N)){
      
      i <- sample.int(n = M - 1, size=1)
      
      temp <- x[i]
      x[i] <- x[i+1]
      x[i+1] <- temp
      
    }
    invs[jjj] <- (1-fastCI(as.numeric(seq_len(M)), x, noise.ties = FALSE)[[1]])*choose(M, 2)
    
  }
  
  exp_inv(N,M)
  hist(invs)
  
  t2 <- seq(0,500)
  t1 <- sapply(t2, function(x) return(exp_inv(x,5)))
  plot(t2,t1)
}
