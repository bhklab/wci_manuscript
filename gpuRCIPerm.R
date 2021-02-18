



# .Call("permCUDA", as.numeric(1:100), runif(100), 1e7, 100, 0, 0, runif(1))



rCI.GPU.perm.test <- function(x, y, delta, req_alpha, outties = 1L){
    tryCatch(dyn.load("naivePermGPU.so"), error=function(e) warning("GPU library load failed"))
    outties <- as.integer(!outties)
    t0res <- wCI:::newPCI(as.numeric(x), as.numeric(y), 
                          as.numeric(length(x)), 
                          pdeltaX = as.numeric(delta), 
                          pdeltaY = as.numeric(delta), 
                          pxties = as.integer(outties), 
                          pyties = as.integer(outties), plogic = 1L)
    t0 <- t0res[1]/(t0res[1] + t0res[2])

    R <- ceiling(1/req_alpha)*100.0
	  t <- .Call("permCUDA", as.numeric(x), as.numeric(y), R, as.numeric(length(x)), delta, delta, runif(1))
    browser()
    p.value <- mean(abs(t-0.5)>abs(t0-0.5))

    return(c(list(t0 = t0,p.value= p.value)))
}


