pearson.perm.test <- function(x, y, req_alpha){

    t0 <- cor(x,y, method="pearson")


    x <- scale(x)
    y <- scale(y)
    x <- x/sum(x^2)
    sampleFun <- function(){
      y <- sample(y)
      return(abs(crossprod(x,y)) > abs(t0))
    }

    p.value <- quickstop(sampleFun, req_alpha = req_alpha)

    return(c(list(t0 = t0), p.value))
}

rCI.perm.test <- function(x, y, delta, req_alpha, outties = 1L){
    outties <- as.integer(!outties)
    t0res <- wCI:::newPCI(as.numeric(x), as.numeric(y), 
                          as.numeric(length(x)), 
                          pdeltaX = as.numeric(delta), 
                          pdeltaY = as.numeric(delta), 
                          pxties = as.integer(outties), 
                          pyties = as.integer(outties), plogic = 1L)
    t0 <- t0res[1]/(t0res[1] + t0res[2])


    ## This function returns T or F depending on whether the permutation was larger than the observed statistic
    ## It saves the x and y variables in its closure, and does the resampling using "sample"
    sampleFun <- function(){
      y <- sample(y)
      t <- NaN
      while(is.nan(t)){
        tres <- wCI:::newPCI(as.numeric(x), as.numeric(y), 
                    as.numeric(length(x)), 
                    pdeltaX = as.numeric(delta), 
                    pdeltaY = as.numeric(delta), 
                    pxties = as.integer(outties), 
                    pyties = as.integer(outties), plogic = 1L)
        t <- tres[1]/(tres[1] + tres[2])
        if(is.nan(t)) stop("T is nan")
      }
      ## sampling from only the permutations where I can calculate rCI
      return(abs(t-0.5) > abs(t0 - 0.5))
    }

    if(is.nan(t0)) {return(t=t0, p.value = 1)}
    p.value <- quickstop(sampleFun, req_alpha = req_alpha)

    return(c(list(t0 = t0), p.value))
}


kCI.perm.test <- function(x, y, req_alpha) {

  t0res <- wCI:::KCI(as.numeric(x), as.numeric(y), 
                          as.numeric(length(x)), 
                          psymmetric = 1L, pkern1 = -27.5512, pkern2 = 0.8000)
    t0 <- t0res[1]/(t0res[1] + t0res[2])

    sampleFun <- function(){
          y <- sample(y)
          t <- NaN
 
          tres <- wCI:::KCI(as.numeric(x), as.numeric(y), 
                          as.numeric(length(x)), 
                          psymmetric = 1L, pkern1 = -27.5512, pkern2 = 0.8000)
          t <- tres[1]/(tres[1] + tres[2])

          return(abs(t-0.5) > abs(t0 - 0.5))
        }
    p.value <- quickstop(sampleFun, req_alpha = req_alpha)

    return(c(list(t0 = t0), p.value))
}

