library(wCI)


log_denom <- function(suc, total, p){
    tmp <- 0;
    # if(log10(p)>-32)
    # {
     tmp <- (suc)*log(p)
    # } else {
    #  warning("p reached precision threshold")
    # }
    # if(log10(1.0-p)>-32)
    # {
     tmp <- tmp + (total - suc)*log(1 - p)
    # } else {
    #    warning("1-p reached precision threshold")
    # }
    return(tmp)

}
log_denom <- function(suc, total, p){
    return((suc)*log(p) + (total - suc)*log(1 - p))

}

### Implementing algorithm from quick stop paper

## decision boundary is inverse of error prob
## do everything in log scale because of numerical precision.


#' Quickstop An algorithm for early stopping during permutation testing
#' 
#' This is implemented from the paper: A flexible and nearly optimal sequential testing approach to randomized testing: QUICK-STOP
#' by Hecker et al.
#' 
#' @param sample_function A function which does a permutation sample, and returns whether or not that sample was more extreme that the observed
#'                        statistic. This simply is a pass-fail for the random data to be "more interesting" than the null, and allows for the 
#'                        permutation test logic to be abstracted away from the stopping procedure
#' @param req_alpha The threshold below which we will deem a p-value to be significant/reject the null
#' @param tolerance_par The algorithm requires a tolerance region where the hypothesis remains undecided. Between req_alpha and req_alpha*(1+tolerance_par)
#'                      The algorithm cannot decide on the significance of the results. For this case, we implement a maximum number of interations equal to
#'                      (1/req_alpha)*100. 
#' @param log_decision_boundary Can be interpereted as the log10 of the odds ratio between the two pairs of hypotheses 
#'                              (p less than req_alpha or not, and p greater than  req_alpha*(1+tolerance_par)) at which the algorithm stops. 
#'                              In effect, if the algorithm return a significant result, we have a 10^(-lDB) probability of that being a false positive. 
#' @param max_iter A maximum number of iterations after which this algorithm quits. Here to take into account the very rare case that the real p value falls
#'                 into the in-decision region. Defaults to 1/(req_alpha)*100 
quickstop <- function(sample_function, req_alpha=0.05, tolerance_par = req_alpha*0.001, log_decision_boundary = 10, max_iter = 1/req_alpha*100){

  num.larger <- 0

  cur_success <- 0
  cur_iter <- 1

  log_cur_PiN <- log(1) # Initialization so the logic can stay the same through all loops

  p1 <- req_alpha
  p2 <- p1 + tolerance_par

  pr_min_1 <- 1/2

  while(cur_iter < max_iter){
    # vec1 <- sample(vec1)
    # perm.cor <- coop::pcor(vec1, vec2, use="complete.obs")

    # success <- abs(perm.cor) > abs(obs.cor)

    success <- sample_function()

    if(success){
      cur_success <- cur_success + 1
      log_cur_PiN <- log_cur_PiN + log_denom(1,1,pr_min_1)
    } else {
      log_cur_PiN <- log_cur_PiN + log_denom(0,1,pr_min_1)
    }
    # if(pr_min_1 >= p2){
    #   log_cur_suph1 <- log_denom(cur_success, cur_iter, p1)
    #   log_cur_suph2 <- log_denom(cur_success, cur_iter, pr_min_1)
    # } else if(pr_min_1 <= p1){
    #   log_cur_suph1 <- log_denom(cur_success, cur_iter, pr_min_1)
    #   log_cur_suph2 <- log_denom(cur_success, cur_iter, p2)
    # } else {
    #   log_cur_suph1 <- log_denom(cur_success, cur_iter, p1)
    #   log_cur_suph2 <- log_denom(cur_success, cur_iter, p2)
    # }
    if(pr_min_1<p1) {
      log_cur_suph1 <- log_denom(cur_success, cur_iter, pr_min_1)
    } else {
      log_cur_suph1 <- log_denom(cur_success, cur_iter, p1)
    }
    if(pr_min_1>p2) {
      log_cur_suph2 <- log_denom(cur_success, cur_iter, pr_min_1)
    } else {
      log_cur_suph2 <- log_denom(cur_success, cur_iter, p2)
    }

    cur_iter <- cur_iter + 1
    pr_min_1 <- (cur_success + 1/2)/cur_iter

    # if(cur_success == 0){
    #   next
    # }
    
    if(log_cur_PiN - log_cur_suph2 > log_decision_boundary){
      return(list(significant = TRUE, "p.value" = pr_min_1, num_iter=cur_iter, num_larger=cur_success))
    }
    if(log_cur_PiN - log_cur_suph1 > log_decision_boundary){
      return(list(significant = FALSE, "p.value" = pr_min_1, num_iter=cur_iter, num_larger=cur_success))
    } 
  }
  return(list(significant = NA, "p.value" = pr_min_1, num_iter=cur_iter, num_larger=cur_success))
}
