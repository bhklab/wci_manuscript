
library(rmutil)

library(data.table)
pset.dir <- "PSets"
library(PharmacoGx)

# ## Loading in drug sensitivity data to serve as test cases
# CCLE <- readRDS(file.path(pset.dir, "CCLE.rds"))
GDSCv2 <- readRDS(file.path(pset.dir, "GDSC2.rds"))
# rownames(GDSCv2@molecularProfiles$rna) <- gsub(rownames(GDSCv2@molecularProfiles$rna), pat="_at$", rep="")

# AZ <- readRDS(file.path(pset.dir, "AZ.rds"))
CTRPv2 <- readRDS(file.path(pset.dir, "CTRPv2.rds"))
# rownames(CCLE.CTRPv2@molecularProfiles$rnaseq) <- gsub(rownames(CCLE.CTRPv2@molecularProfiles$rnaseq), pat="\\.[0-9]*$", rep="")

gCSI <- readRDS(file.path(pset.dir, "gCSI2.rds"))
# rownames(gCSI@molecularProfiles$rnaseq) <- gsub(rownames(gCSI@molecularProfiles$rnaseq), pat="\\.[0-9]*$", rep="")

GRAY2017 <- readRDS(file.path(pset.dir, "GRAY2017.rds"))
# AZ <- readRDS("~/Data/PSets/AZ.rds")



allAAC <- c(sensitivityProfiles(GDSCv2)[,"aac_recomputed"],
			 sensitivityProfiles(CTRPv2)[,"aac_recomputed"],
			 sensitivityProfiles(gCSI)[,"aac_recomputed"],
			 sensitivityProfiles(GRAY2017)[,"aac_recomputed"])
allAAC <- allAAC[!is.na(allAAC)]
write.csv(allAAC, file="allAACs.csv")





exps <- paste0(GRAY2017@sensitivity$info$cellid, "~", GRAY2017@sensitivity$info$drugid)


unique_dups <- unique(exps[duplicated(exps)])
tbl_exps <- table(exps)

nAACs <- sum(sapply(tbl_exps[tbl_exps>1], function(x) return(choose(x,2))))
delta_AACs <- numeric(nAACs)
i <- 1
paired_aacs <- matrix(NA_real_, nrow=0, ncol=2)
drug_cell <- matrix(NA_character_, nrow=0, ncol=3)
for (exp in unique_dups){
	myx <- exps == exp
	aacs <- GRAY2017@sensitivity$profiles[myx,"aac_recomputed"]
	paired_aacs <- rbind(paired_aacs, t(combn(aacs, 2)))
	deltas <- combn(rev(aacs), 2, FUN = function(x) {return(x[1] - x[2])})
	delta_AACs[seq(i,i+length(deltas)-1)] <- deltas
	drug_cell <- rbind(drug_cell, matrix(c("GRAY2017", unlist(strsplit(exp, split="~"))), nrow=length(deltas), ncol=3, byrow=TRUE))

	i <- i + length(deltas)
}

GRAY2017.replicates <- data.table(data.frame(drug_cell, paired_aacs))
# rm(GRAY2017)

# hist(apply(paired_aacs, 1, function(x) return(x[1]-x[2])), breaks=200)

#GDSCv2

exps <- paste0(GDSCv2@sensitivity$info$cellid, "~", GDSCv2@sensitivity$info$drugid)


unique_dups <- unique(exps[duplicated(exps)])
tbl_exps <- table(exps)

nAACs <- sum(sapply(tbl_exps[tbl_exps>1], function(x) return(choose(x,2))))
delta_AACs <- numeric(nAACs)
i <- 1
paired_aacs <- matrix(NA_real_, nrow=nAACs, ncol=2)
drug_cell <- matrix(NA_character_, nrow=nAACs, ncol=3)
pb <- txtProgressBar(0, nAACs, style=3)

for (exp in unique_dups){
	myx <- exps == exp
	aacs <- GDSCv2@sensitivity$profiles[myx,"aac_recomputed"]
	deltas <- combn(rev(aacs), 2, FUN = function(x) {return(x[1] - x[2])})
	delta_AACs[seq(i,i+length(deltas)-1)] <- deltas
	paired_aacs[seq(i,i+length(deltas)-1),] <- t(combn(aacs, 2))
	drug_cell[seq(i,i+length(deltas)-1),] <- matrix(c("GDSC2", unlist(strsplit(exp, split="~"))), nrow=length(deltas), ncol=3, byrow=TRUE)

	i <- i + length(deltas)
	setTxtProgressBar(pb, i)	
}

GDSCv2.replicates <- data.table(data.frame(drug_cell, paired_aacs))
# rm(GDSCv2)
# hist(apply(paired_aacs, 1, function(x) return(x[1]-x[2])), breaks=200)


#CTRPv2

exps <- paste0(CTRPv2@sensitivity$info$cellid, "~", CTRPv2@sensitivity$info$drugid)


unique_dups <- unique(exps[duplicated(exps)])
tbl_exps <- table(exps)

nAACs <- sum(sapply(tbl_exps[tbl_exps>1], function(x) return(choose(x,2))))
delta_AACs <- numeric(nAACs)
i <- 1
paired_aacs <- matrix(NA_real_, nrow=nAACs, ncol=2)
drug_cell <- matrix(NA_character_, nrow=nAACs, ncol=3)
pb <- txtProgressBar(0, nAACs, style=3)

for (exp in unique_dups){
	myx <- exps == exp
	aacs <- CTRPv2@sensitivity$profiles[myx,"aac_recomputed"]
	deltas <- combn(rev(aacs), 2, FUN = function(x) {return(x[1] - x[2])})
	delta_AACs[seq(i,i+length(deltas)-1)] <- deltas
	paired_aacs[seq(i,i+length(deltas)-1),] <- t(combn(aacs, 2))
	drug_cell[seq(i,i+length(deltas)-1),] <- matrix(c("CTRPv2", unlist(strsplit(exp, split="~"))), nrow=length(deltas), ncol=3, byrow=TRUE)

	i <- i + length(deltas)
	setTxtProgressBar(pb, i)	
}

CTRPv2.replicates <- data.table(data.frame(drug_cell, paired_aacs))
# rm(CTRPv2)

## gCSI

exps <- paste0(gCSI@sensitivity$info$cellid, "~", gCSI@sensitivity$info$drugid)


unique_dups <- unique(exps[duplicated(exps)])
tbl_exps <- table(exps)

nAACs <- sum(sapply(tbl_exps[tbl_exps>1], function(x) return(choose(x,2))))
delta_AACs <- numeric(nAACs)
i <- 1
paired_aacs <- matrix(NA_real_, nrow=nAACs, ncol=2)
drug_cell <- matrix(NA_character_, nrow=nAACs, ncol=3)
pb <- txtProgressBar(0, nAACs, style=3)

for (exp in unique_dups){
	myx <- exps == exp
	aacs <- gCSI@sensitivity$profiles[myx,"aac_recomputed"]
	deltas <- combn(rev(aacs), 2, FUN = function(x) {return(x[1] - x[2])})
	delta_AACs[seq(i,i+length(deltas)-1)] <- deltas
	paired_aacs[seq(i,i+length(deltas)-1),] <- t(combn(aacs, 2))
	drug_cell[seq(i,i+length(deltas)-1),] <- matrix(c("gCSI2", unlist(strsplit(exp, split="~"))), nrow=length(deltas), ncol=3, byrow=TRUE)

	i <- i + length(deltas)
	setTxtProgressBar(pb, i)	
}

gCSI.replicates <- data.table(data.frame(drug_cell, paired_aacs))
# rm(gCSI)

## AZ

# exps <- paste0(AZ@sensitivity$info$cellid, "~", AZ@sensitivity$info$drugid)


# unique_dups <- unique(exps[duplicated(exps)])
# tbl_exps <- table(exps)

# nAACs <- sum(sapply(tbl_exps[tbl_exps>1], function(x) return(choose(x,2))))
# delta_AACs <- numeric(nAACs)
# i <- 1
# paired_aacs <- matrix(NA_real_, nrow=nAACs, ncol=2)
# drug_cell <- matrix(NA_character_, nrow=nAACs, ncol=3)
# pb <- txtProgressBar(0, nAACs, style=3)

# for (exp in unique_dups){
# 	myx <- exps == exp
# 	aacs <- AZ@sensitivity$profiles[myx,"aac_recomputed"]
# 	deltas <- combn(rev(aacs), 2, FUN = function(x) {return(x[1] - x[2])})
# 	delta_AACs[seq(i,i+length(deltas)-1)] <- deltas
# 	paired_aacs[seq(i,i+length(deltas)-1),] <- t(combn(aacs, 2))
# 	drug_cell[seq(i,i+length(deltas)-1),] <- matrix(c("AZ", unlist(strsplit(exp, split="~"))), nrow=length(deltas), ncol=3, byrow=TRUE)

# 	i <- i + length(deltas)
# 	setTxtProgressBar(pb, i)	
# }

# AZ.replicates <- data.table(data.frame(drug_cell, paired_aacs))
# # rm(AZ)



### Make plot per drug first

all.replicates <- rbindlist(list(GRAY2017.replicates, GDSCv2.replicates, CTRPv2.replicates, gCSI.replicates))


colnames(all.replicates) <- c("Dataset", "Cell Line", "Drug", "AAC1", "AAC2")
all.replicates[,deltaAAC := AAC1 - AAC2]


fwrite(all.replicates, file="~/Code/mCI_paper_simulations/all.replicates.csv")


## now we sample 700000 random pairs


i = 0
random_pairs <- matrix(NA_real_, nrow=7e5, ncol=2)

while(i < 7e5){
	i <- i + 1
	i1 <- sample(length(allAAC), 1)
	i2 <- sample(length(allAAC), 1)
	while(i2 == i1) {i2 <- sample(length(allAAC), 1)}
	random_pairs[i, ] <- c(allAAC[i1], allAAC[i2])
	if(i %% 10000) print(i)
}
write.csv(random_pairs, "random_pairs_aac.csv")

source("optimize_delta.R")
library(Hmisc)
optimal_delta <- optimize_delta(all.replicates$deltaAAC, random_pairs[,1] - random_pairs[,2])


optimal_delta$mcc[which.max(optimal_delta$mcc[,2]),]

optimal_delta$f1[which.max(optimal_delta$f1[,2]),]


source("fit_kci_kernel.R")

fit_kci_kernel(na.omit(random_pairs[,1] - random_pairs[,2]), na.omit(all.replicates$deltaAAC), make_plot=TRUE)

ecdfRep <- ecdf(na.omit(all.replicates$deltaAAC))

optfun <- function(s) {

	logprobs <- log(1/2/s) - abs(na.omit(all.replicates$deltaAAC))/s

	maxVal <- max(logprobs)

	return(maxVal+ log(sum(exp(logprobs - maxVal))))
}


optfun <- function(s) {
	
	logprobs <- log(1/2/s) - abs(na.omit(all.replicates$deltaAAC))/s

	maxVal <- max(logprobs)

	return(sum(logprobs))
}


optfunG <- function(s) {

	logprobs <- log(1/sqrt(2*pi)/s) - 0.5*(na.omit(all.replicates$deltaAAC)/s)^2

	maxVal <- max(logprobs)

	return(sum(logprobs))
}


optLaplacian <- optimise( optfun, c(0,100),maximum=TRUE)$maximum

optGaussian <- optimise( optfunG, c(0,100),maximum=TRUE)$maximum

plotX <- seq(-1,1,0.01)
plot(ecdfRep)
lines(plotX, pnorm(plotX, sd=optGaussian), col="red")
lines(plotX, plaplace(plotX, s=optLaplacian), col="blue")



plotX <- seq(-1,1,0.01)
plot(density(na.omit(all.replicates$deltaAAC)))
lines(plotX, dnorm(plotX, sd=optGaussian), col="red")
lines(plotX, dlaplace(plotX, s=optLaplacian), col="blue")
lines(plotX, dlaplace(plotX, s=mad(na.omit(all.replicates$deltaAAC))), col="green")


print(optLaplacian)




# optBeta <- function(par){
# 	return(sum(-dbeta(aac[aac>0], par[[1]], par[[2]], log=TRUE), na.rm=TRUE))
# }


# optim(list(1,2), optBeta)


# optTN <- function(par){
# 	return(sum(-log(dtruncnorm(aac[aac>0], 0,1, par[[1]], par[[2]])+.Machine$double.eps), na.rm=TRUE))
# }


# optim(list(0.1,3), optTN)



