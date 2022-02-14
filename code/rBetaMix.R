library(iterators)
library(foreach)

library(profmem)

# Looping is very slow and also memory intensive??? 
# RmixtureBeta <- function(N, pi_0=0.8, shape.vec.1=c(1.2,4), shape.vec.2=c(8,2)){
#  
#   i <- 1
#   res <- numeric(N)
#   while(i <= N){
#     if(runif(1) < pi_0){
#       res[i] <- rbeta(1, shape.vec.1[1], shape.vec.1[2])
#     } else {
#       res[i] <- rbeta(1, shape.vec.2[1], shape.vec.2[2])
#     }
#     i <- i + 1
#   } 
# 
#   return(res)
# }
# 
# RmixtureBeta.2 <- function(N, pi_0=0.8, shape.vec.1=c(1.2,4), shape.vec.2=c(8,2)){
#   
#   rb1 <- rbeta(N, shape.vec.1[1], shape.vec.1[2])
#   rb2 <- rbeta(N, shape.vec.2[1], shape.vec.2[2])
#   res <- ifelse(runif(N) < pi_0, rb1, rb2)
#   
#   return(res)
# }
# 

RmixtureBeta <- function(N, pi_0=0.8, shape.vec.1=c(1.2,4), shape.vec.2=c(8,2)){
  
  rb1 <- rbeta(N, shape.vec.1[1], shape.vec.1[2])
  rb2 <- rbeta(N, shape.vec.2[1], shape.vec.2[2])
  pivec <- (runif(N) < pi_0)
  res <-  pivec* rb1 + (!pivec)*rb2
  
  return(res)
}

# 
# test2 <- sapply(1:10, function(x) {system.time(RmixtureBeta.2(5000000, pi_0=0.9))[3]})
# test <- sapply(1:10, function(x) {system.time(RmixtureBeta(5000000, pi_0=0.9))[3]})
# boxplot(test, test2)
test <- numeric(10000)
for(i in seq_along(test)){
  test[i] <- cor.test(RmixtureBeta(1000), RmixtureBeta(1000))$p.value
}
