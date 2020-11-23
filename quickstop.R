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

# partialPermute <- function(obs.cor, data, formula1, formula2, req_alpha = 0.05, tolerance_par = req_alpha*0.1, error_prob = 1e-10){
  
#   num.larger <- 0


#   for(i in seq_len(req_alpha)){
#     data[,1] <- sample(data[,1], nrow(data))
#     data[,2] <- sample(data[,2], nrow(data))

#     partial.dp <- residuals(lm(formula(formula1), data))
#     partial.x <- residuals(lm(formula(formula2), data))

#     perm.cor <- coop::pcor(partial.dp, partial.x, use="complete.obs")

#     if(abs(perm.cor) >= abs(obs.cor)){
#       num.larger <- num.larger + 1
#     }
#   }

#   return(num.larger/req_alpha)
# }

