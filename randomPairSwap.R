#####################################################################
### Code to investigate the behaviours of *CIs with increasing random 
### pair swaps.
###
### Petr Smirnov, 2019
#####################################################################


### Just a dirty script to plot *CI values over random pair swaps, checking that they
### Evolve down to 0.5 as expected. Dont see much reuse out of this. 

library(wCI)
library(MASS)
library(tictoc)
library(doParallel)
registerDoParallel(cores=20)
source("adjecentSwaps.R")
source("~/Code/mCI_paper_simulations/kmci.R")
library(abind)


## the dyn.load command should be replacable by loading the latest fastCI package
## dyn.load("src/swapRandomPairs.so")
# library(fastCI)

corFuns <- list('CI' = function(x,y) return(paired.concordance.index.new(x,y,delta.pred = 0, delta.obs = 0, p_method="As", conf_int="As", outx=T, outy=T)[[1]]),
  'rCI' = function(x,y) return(paired.concordance.index.new(x,y,delta.pred = 0.1, delta.obs = 0.1, p_method="As", conf_int="As", outx=T, outy=T)[[1]]),
  'kCI' = function(x,y) {
    res <- wCI:::KCI(as.numeric(x), as.numeric(y), 
                          as.numeric(length(x)), 
                          psymmetric = 1L, pkern1 = -27.52480305, pkern2 = 0.06463092) 
    return(res[1]/(res[1] + res[2]))},
  'pearson' = function(x,y) cor(x,y, use="pairwise"),
  'spearman' = function(x,y)cor(x,y, use="pairwise", method="spearman"))

bigSwappy <- function(x,k) return(wCI:::do_random_swaps_c(x, as.numeric(length(x)), as.numeric(k), as.numeric(runif(2))))
littleSwappy <- function(x,k) return(wCI:::do_swaps_c(x, as.numeric(length(x)), as.numeric(k), as.numeric(runif(2))))

runEvolutionOfStats <- function(sampleFun = runif, N = 1000, sampleN = 100, 
                                noisePar = seq(1, 5, by = 0.2), corFunctions = corFuns, noiseFun = bigSwappy){


  res <- foreach(k = 10^noisePar) %:% foreach(j = seq_len(sampleN), .combine=rbind) %dopar% {


    x <- sort(sampleFun(N))
    y <- x
    x <- noiseFun(x,k)
    sapply(corFuns, function(f) return(f(x,y)))
  }
  res <- abind(res, along=-1)
  dimnames(res)[[1]] <- 10^noisePar
  return(res)
}


plotSwapRes <- function(resArray, title="Completely Random Swaps") {

  library(reshape2)
  library(data.table)
  plt <- data.table(melt(resArray))

  colnames(plt) <- c("Swaps", "perm","Metric", "Value")

  plt[,Stat := mean(Value),.(Swaps, Metric)]
  plt[,pos := Stat + (1.96*sd(Value)/sqrt(.N)),.(Swaps, Metric)]
  plt[,neg := Stat - (1.96*sd(Value)/sqrt(.N)),.(Swaps, Metric)]


  toPlot <- unique(plt[Metric %in% c("CI", "kCI", "rCI"),.(Swaps, Stat, Metric, pos, neg)])

  library(ggplot2)


  pres_ready <- theme_bw() + 
    theme(axis.title = element_text(size=16), axis.text = element_text(size=16), legend.text = element_text(size=16), title = element_text(size=20),legend.key.height = unit(1.0, 'cm')) 
  ggplot(toPlot, aes(x=Swaps, y=Stat, color=Metric)) + geom_errorbar(aes(ymin=neg, ymax=pos), width=.1) + geom_line() + scale_x_log10() + ggtitle(title) + pres_ready

}

mySeq <- seq(1, 5, by = 0.2)
bigSwappyUniformRes <- runEvolutionOfStats(noisePar=mySeq)

save(bigSwappyUniformRes, file="big_swappy_CIs_up_to_10e5_and.RData")
# load("big_swappy_CIs_up_to_10e5.RData")
mySeq <- seq(1, 5, by = 0.2)




pdf("figures/CI_evolution_over_noise_uniform_random_swaps_and.pdf", height=5, width=10)
print(plotSwapRes(bigSwappyUniformRes))
dev.off()


beta110 <- function(N) return(rbeta(N, 1, 10))

beta110res <- runEvolutionOfStats(sampleFun = beta110, noisePar=mySeq)


pdf("figures/CI_evolution_over_noise_beta_1_10_random_swaps_and.pdf", height=5, width=10)
print(plotSwapRes(beta110res))
dev.off()



beta1245 <- function(N) return(rbeta(N, 1.2, 4.5))

beta1245res <- runEvolutionOfStats(sampleFun = beta1245, noisePar=mySeq)
plotSwapRes(beta1245res)
####################
### Now we look at neighbour swaps
####################


mySeq <- seq(1, 9, by = 0.5)
littleSwappyUniformRes <- runEvolutionOfStats(noisePar=mySeq, noiseFun=littleSwappy)
save(littleSwappyUniformRes, file="pairwise_swappy_CIs_up_to_10e9_and.RData")


pdf("figures/CI_evolution_over_noise_uniform_adjecent_swaps_and.pdf", height=5, width=10)
print(plotSwapRes(littleSwappyUniformRes[-(1:8),,]), title="Adjacent Swaps")
dev.off()




littleSwappybeta110Res <- runEvolutionOfStats(sampleFun = beta110, noisePar=mySeq, noiseFun=littleSwappy)



pdf("figures/CI_evolution_over_noise_beta_1_10_adjecent_swaps_and.pdf", height=5, width=10)
print(plotSwapRes(littleSwappybeta110Res[-(1:8),,]), title="Adjacent Swaps")
dev.off()

####################
### Now lets look at adding noise
####################
# 
# pltCIs <- melt(CIs)
# colnames(pltCIs) <- c("Swaps", "perm","Value")
# pltCIs$Metric <- "CI"
# 
# pltrCIs <- melt(rCIs)
# colnames(pltrCIs) <- c("Swaps","perm","Value")
# pltrCIs$Metric <- "rCI"
# # 
# # pltPearsons <- melt(pearsons)
# # colnames(pltPearsons) <- c("Swaps","perm","Value")
# # pltPearsons$Metric <- "Pearson Rescaled"
# # pltPearsons[,3] <- pltPearsons[,3]/2 + 0.5
# # 
# # 
# # pltSpearmans <- melt(spearmans)
# # colnames(pltSpearmans) <- c("Swaps","perm","Value")
# # pltSpearmans$Metric <- "Spearman Rescaled"
# # pltSpearmans[,3] <- pltSpearmans[,3]/2 + 0.5
# 
# 
# toPlot <- rbind(pltCIs, pltrCIs)#, pltPearsons, pltSpearmans)
# 
# 
# pdf("figures/AllMetrics_evolution_over_noise_uniform_random_swaps_and.pdf", height=5, width=10)
# ggplot(toPlot, aes(x=Swaps, y=Value, color=Metric)) + geom_point() + scale_x_log10() + geom_smooth() + ggtitle("Completely Random Swaps") + pres_ready
# dev.off()


# dyn.load("src/swapPairs.so")



# N <- 1000
# sampleN <- 100

# mySeq <- seq(1, 9, by = 0.5)


# CIs <- matrix(NA_real_, nrow=length(mySeq), ncol = sampleN)
# rCIs <- matrix(NA_real_, nrow=length(mySeq), ncol = sampleN)
# pearsons <- matrix(NA_real_, nrow=length(mySeq), ncol = sampleN)
# spearmans <- matrix(NA_real_, nrow=length(mySeq), ncol = sampleN)
# i <- 1

# # mySeq <- seq(1, 9, by = 0.5)

# for(k in 10^mySeq){
#   for(j in seq_len(100)){
#     x <- sort(runif(N))
#     y <- x
#     x <- wCI:::do_swaps_c(x, as.numeric(length(x)), as.numeric(k), as.numeric(runif(2)))
    
#     CIs[i, j] <- paired.concordance.index.new(x,y,delta.pred = 0, delta.obs = 0, p_method="As", conf_int="As", outx=T, outy=T)[[1]]
#     rCIs[i, j] <- paired.concordance.index.new(x,y,delta.pred = 0.2, delta.obs = 0.2, p_method="As", conf_int="As", outx=T, outy=T)[[1]]
#     pearsons[i,j] <- cor(x,y, use="pairwise")
#     spearmans[i,j] <- cor(x,y, use="pairwise", method="spearman")
    
#   }
#   i <- i + 1
#   print(i)
# }

# rownames(CIs) <- rownames(rCIs) <- rownames(pearsons) <- rownames(spearmans) <- 10^mySeq


# save(CIs, rCIs, pearsons, spearmans, file="pairwise_swappy_CIs_up_to_10e9_and.RData")
# # load("pairwise_swappy_CIs_up_to_10e9.RData")
# mySeq <- seq(1, 9, by = 0.5)

# library(reshape2)
# library(data.table)
# pltCIs <- data.table(melt(CIs))
# colnames(pltCIs) <- c("Swaps", "perm","Value")
# pltCIs$Metric <- "CI"
# pltCIs[,Stat := mean(Value),Swaps]
# pltCIs[,pos := Stat + (1.96*sd(Value)/sqrt(.N)),Swaps]
# pltCIs[,neg := Stat - (1.96*sd(Value)/sqrt(.N)),Swaps]



# pltrCIs <- melt(rCIs)
# pltrCIs <- data.table(pltrCIs)
# colnames(pltrCIs) <- c("Swaps","perm","Value")
# pltrCIs$Metric <- "rCI"
# pltrCIs[,Stat := mean(Value),Swaps]
# pltrCIs[,pos := Stat + (1.96*sd(Value)/sqrt(.N)),Swaps]
# pltrCIs[,neg := Stat - (1.96*sd(Value)/sqrt(.N)),Swaps]

# toPlot <- rbind(pltCIs, pltrCIs)

# library(ggplot2)


# # pltCIs <- melt(CIs)
# # colnames(pltCIs) <- c("Swaps", "perm","Value")
# # pltCIs$Metric <- "CI"
# # 
# # pltrCIs <- melt(rCIs)
# # colnames(pltrCIs) <- c("Swaps","perm","Value")
# # pltrCIs$Metric <- "rCI"
# # 
# # pltPearsons <- melt(pearsons)
# # colnames(pltPearsons) <- c("Swaps","perm","Value")
# # pltPearsons$Metric <- "Pearson Rescaled"
# # pltPearsons[,3] <- pltPearsons[,3]/2 + 0.5
# # 
# # 
# # pltSpearmans <- melt(spearmans)
# # colnames(pltSpearmans) <- c("Swaps","perm","Value")
# # pltSpearmans$Metric <- "Spearman Rescaled"
# # pltSpearmans[,3] <- pltSpearmans[,3]/2 + 0.5
# # 
# # 
# # toPlot <- rbind(pltCIs, pltrCIs, pltPearsons, pltSpearmans)
# # 
# # 
# # pdf("figures/AllMetrics_evolution_over_noise_uniform_adjecent_swaps.pdf", height=5, width=10)
# # ggplot(toPlot, aes(x=Swaps, y=Value, color=Metric)) + geom_point() + scale_x_log10() + geom_smooth() + ggtitle("Adjacent Random Swaps") + pres_ready
# # dev.off()
# # 
# # 



# N <- 1000
# sampleN <- 100
# noise_vector <- seq(0,3, by=0.05)
# CIs <- matrix(NA_real_, nrow=length(noise_vector), ncol = sampleN)
# rCIs <- matrix(NA_real_, nrow=length(noise_vector), ncol = sampleN)

# pearsons <- matrix(NA_real_, nrow=length(noise_vector), ncol = sampleN)
# spearmans <- matrix(NA_real_, nrow=length(noise_vector), ncol = sampleN)

# i <- 1
# for(sgm in noise_vector){
#   for(j in seq_len(100)){
#     x <- sort(runif(N))
#     y <- x
#     x <- x + rnorm(length(x), 0, sgm/1.959964)
    
#     CIs[i, j] <- paired.concordance.index.new(x,y,delta.pred = 0, delta.obs = 0, p_method="As", conf_int="As", outx=T, outy=T)[[1]]
#     rCIs[i, j] <- paired.concordance.index.new(x,y,delta.pred = 0.2, delta.obs = 0.2, p_method="As", conf_int="As", outx=T, outy=T)[[1]]
#     pearsons[i,j] <- cor(x,y, use="pairwise")
#     spearmans[i,j] <- cor(x,y, use="pairwise", method="spearman")
#   }
#   i <- i + 1
#   print(i)
# }

# rownames(CIs) <- rownames(rCIs) <- rownames(pearsons) <- rownames(spearmans) <- noise_vector
# save(CIs, rCIs, pearsons, spearmans, file="noising_up_CIs_and.RData")

# # load("noising_up_CIs.RData")

# library(reshape2)
# library(data.table)
# pltCIs <- data.table(melt(CIs))
# colnames(pltCIs) <- c("Noise", "perm","Value")
# pltCIs$Metric <- "CI"
# pltCIs[,Stat := mean(Value),Noise]
# pltCIs[,pos := Stat + (1.96*sd(Value)/sqrt(.N)),Noise]
# pltCIs[,neg := Stat - (1.96*sd(Value)/sqrt(.N)),Noise]



# pltrCIs <- melt(rCIs)
# pltrCIs <- data.table(pltrCIs)
# colnames(pltrCIs) <- c("Noise","perm","Value")
# pltrCIs$Metric <- "rCI"
# pltrCIs[,Stat := mean(Value),Noise]
# pltrCIs[,pos := Stat + (1.96*sd(Value)/sqrt(.N)),Noise]
# pltrCIs[,neg := Stat - (1.96*sd(Value)/sqrt(.N)),Noise]
# toPlot <- rbind(pltCIs, pltrCIs)


# library(ggplot2)
# pdf("figures/CI_evolution_over_noise_adding_gaussian_and.pdf", height=5, width=10)
# ggplot(toPlot, aes(x=Noise, y=Stat, color=Metric)) + geom_errorbar(aes(ymin=neg, ymax=pos), width=.1) + geom_line() +  ggtitle("Adding Gaussian Noise") + pres_ready + scale_x_log10()
# dev.off()

# library(ggplot2)
# pdf("figures/CI_evolution_over_noise_adding_gaussian_nolog_and.pdf", height=5, width=10)
# ggplot(toPlot, aes(x=Noise, y=Stat, color=Metric)) + geom_errorbar(aes(ymin=neg, ymax=pos), width=.1) + geom_line() +  ggtitle("Adding Gaussian Noise") + pres_ready 
# dev.off()

# library(reshape2)
# pltCIs <- melt(CIs)
# colnames(pltCIs) <- c("Noise", "perm","Value")
# pltCIs$Metric <- "CI"

# pltrCIs <- melt(rCIs)
# colnames(pltrCIs) <- c("Noise","perm","Value")
# pltrCIs$Metric <- "rCI"

# pltPearsons <- melt(pearsons)
# colnames(pltPearsons) <- c("Noise","perm","Value")
# pltPearsons$Metric <- "Pearson Rescaled"
# pltPearsons[,3] <- pltPearsons[,3]/2 + 0.5


# pltSpearmans <- melt(spearmans)
# colnames(pltSpearmans) <- c("Noise","perm","Value")
# pltSpearmans$Metric <- "Spearman Rescaled"
# pltSpearmans[,3] <- pltSpearmans[,3]/2 + 0.5


# toPlot <- rbind(pltCIs, pltrCIs, pltPearsons, pltSpearmans)


# library(ggplot2)
# pdf("figures/AllMetrics_evolution_over_noise_adding_gaussian_and.pdf", height=5, width=10)
# ggplot(toPlot, aes(x=Noise, y=Value, color=Metric)) + geom_point() + scale_x_log10() + geom_smooth() + ggtitle("Adding Gaussian Noise") + pres_ready
# dev.off()
