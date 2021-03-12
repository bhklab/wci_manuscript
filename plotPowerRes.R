library(abind)
library(data.table)
library(ggplot2)
library(reshape2)

MakePowerOverEffectPlot <- function(power_res, sampleSize, alpha = 0.001, delta=1, fileName){

	if(missing(fileName)){
		fileName = paste("power_plot_n", sampleSize, "delta", delta, ".pdf", sep="_")
	}
	n50 <- power_res[,as.character(sampleSize)]

	n50 <- lapply(n50, unlist, recursive=FALSE)
	n50 <- lapply(n50, unlist, recursive=FALSE)
	n50 <- lapply(n50, abind, along=-1)


	contingency_helper <- function(x1, x2){
	  
	  resTbl <- matrix(0, nrow=2, ncol=2)
	  dimnames(resTbl) <- list(c("0", "1"), c("0","1"))
	  resTbl["0","0"] <- sum(x1==0 & x2==0)
	  resTbl["0","1"] <- sum(x1==0 & x2==1)
	  resTbl["1","0"] <- sum(x1==1 & x2==0)
	  resTbl["1","1"] <- sum(x1==1 & x2==1)
	  return(resTbl)
	}

	n50_power <- lapply(n50, function(x) {
	  x[,,2:5] <- apply(x[,,2:5, drop=FALSE], c(1,2,3), function(x) (as.numeric(x < alpha)))
	  powers <- sapply(dimnames(x)[[2]], function(xx){

	      conf_mat_CI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"CI_p"])
	      conf_mat_rCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"rCI_p"])
	      conf_mat_pearson <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Pearson_p"])
	      conf_mat_spearman <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Spearman_p"])
	      # browser()
	      pow_CI <- conf_mat_CI["1","1"]/sum(conf_mat_CI["1",])
	      pow_rCI <- conf_mat_rCI["1","1"]/sum(conf_mat_rCI["1",])
	      pow_pearson <- conf_mat_pearson["1","1"]/sum(conf_mat_pearson["1",])
	      pow_spearman <- conf_mat_spearman["1","1"]/sum(conf_mat_spearman["1",])
	      
	      return(c("CI" = pow_CI, "rCI" = pow_rCI, "Pearson" = pow_pearson, "Spearman" = pow_spearman))
	    })
	  return(powers)
	})

	n50_power <- abind(n50_power, along = -1)

	toPlot <- melt(n50_power[,,as.character(delta)])
	colnames(toPlot) <- c("Effect Size", "Method", "Power")

	pres_ready <- theme_bw() + 
	  theme(axis.title = element_text(size=24), axis.text = element_text(size=24), legend.text = element_text(size=24), title = element_text(size=28),legend.key.height = unit(1.0, 'cm'))

	pdf(fileName, height = 6, width=9)
	p <- ggplot(toPlot, aes(x=`Effect Size`, y=Power, col=Method)) + geom_line(size=1) + ggtitle(paste0("N=", sampleSize, ", Delta = ",delta,", alpha = ", alpha)) + 
				pres_ready + theme(
								    legend.position = c(0.05, 0.95),
								    legend.justification = c("left", "top"),
								    legend.box.just = "left"#,
								    # legend.margin = margin(6, 6, 6, 6)
									)
	print(p)	
	dev.off()


}

MakePowerOverEffectPlotBeta <- function(power_res, sampleSize, alpha = 0.001, delta=0.1,  fileNamePrefix="", fileName){

	if(missing(fileName)){
		fileName = paste("power_plot",fileNamePrefix,"beta_n", sampleSize, "delta", delta, "_withKCI.pdf", sep="_")
	}
	n50 <- power_res[,as.character(sampleSize)]

	# n50 <- lapply(n50, unlist, recursive=FALSE)
	# n50 <- lapply(n50, unlist, recursive=FALSE)
	n50 <- lapply(n50, abind, along=-1)


	contingency_helper <- function(x1, x2){
	  
	  resTbl <- matrix(0, nrow=2, ncol=2)
	  dimnames(resTbl) <- list(c("0", "1"), c("0","1"))
	  resTbl["0","0"] <- sum(x1==0 & x2==0)
	  resTbl["0","1"] <- sum(x1==0 & x2==1)
	  resTbl["1","0"] <- sum(x1==1 & x2==0)
	  resTbl["1","1"] <- sum(x1==1 & x2==1)
	  return(resTbl)
	}

	n50_power <- lapply(n50, function(x) {
	  x[,,2:6] <- apply(x[,,2:6, drop=FALSE], c(1,2,3), function(x) (as.numeric(x < alpha)))
	  powers <- sapply(dimnames(x)[[2]], function(xx){

	      conf_mat_CI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"CI_p"])
	      conf_mat_rCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"rCI_p"])
	      conf_mat_pearson <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Pearson_p"])
	      conf_mat_spearman <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Spearman_p"])
  	      conf_mat_KCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"KCI_p"])

	      # browser()
	      pow_CI <- conf_mat_CI["1","1"]/sum(conf_mat_CI["1",])
	      pow_rCI <- conf_mat_rCI["1","1"]/sum(conf_mat_rCI["1",])
	      pow_pearson <- conf_mat_pearson["1","1"]/sum(conf_mat_pearson["1",])
	      pow_spearman <- conf_mat_spearman["1","1"]/sum(conf_mat_spearman["1",])
	      pow_KCI <- conf_mat_KCI["1","1"]/sum(conf_mat_KCI["1",])

	      return(c("CI" = pow_CI, "rCI" = pow_rCI, "Pearson" = pow_pearson, "Spearman" = pow_spearman, KCI = pow_KCI))
	    })
	  return(powers)
	})

	n50_power <- abind(n50_power, along = -1)

	toPlot <- melt(n50_power[,,as.character(delta)])
	colnames(toPlot) <- c("Effect Size", "Method", "Power")

	pres_ready <- theme_bw() + 
	  theme(axis.title = element_text(size=24), axis.text = element_text(size=24), legend.text = element_text(size=24), title = element_text(size=28),legend.key.height = unit(1.0, 'cm'))

	pdf(fileName, height = 6, width=9)
	p <- ggplot(toPlot, aes(x=`Effect Size`, y=Power, col=Method)) + geom_line(size=1) + ggtitle(paste0("N=", sampleSize, ", Delta = ",delta,", alpha = ", alpha)) + 
				pres_ready + theme(
								    legend.position = c(0.05, 0.95),
								    legend.justification = c("left", "top"),
								    legend.box.just = "left"#,
								    # legend.margin = margin(6, 6, 6, 6)
									)
	print(p)	
	dev.off()


}



MakePowerOverEffectPlotWithKendall <- function(power_res, sampleSize, alpha = 0.001, delta=1, fileName){

	if(missing(fileName)){
		fileName = paste("power_plot_kendall_n", sampleSize, "delta", delta, ".pdf", sep="_")
	}
	n50 <- power_res[,as.character(sampleSize)]

	n50 <- lapply(n50, unlist, recursive=FALSE)
	n50 <- lapply(n50, unlist, recursive=FALSE)
	n50 <- lapply(n50, abind, along=-1)

	contingency_helper <- function(x1, x2){
	  
	  resTbl <- matrix(0, nrow=2, ncol=2)
	  dimnames(resTbl) <- list(c("0", "1"), c("0","1"))
	  resTbl["0","0"] <- sum(x1==0 & x2==0)
	  resTbl["0","1"] <- sum(x1==0 & x2==1)
	  resTbl["1","0"] <- sum(x1==1 & x2==0)
	  resTbl["1","1"] <- sum(x1==1 & x2==1)
	  return(resTbl)
	}

	
	n50_power <- lapply(n50, function(x) {
	  x[,,2:6] <- apply(x[,,2:6, drop=FALSE], c(1,2,3), function(x) (as.numeric(x < alpha)))
	  powers <- sapply(dimnames(x)[[2]], function(xx){

	      conf_mat_CI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"CI_p"])
	      conf_mat_rCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"rCI_p"])
	      conf_mat_pearson <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Pearson_p"])
	      conf_mat_spearman <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Spearman_p"])
  	      conf_mat_kendall <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Kendall_p"])

	      # browser()
	      pow_CI <- conf_mat_CI["1","1"]/sum(conf_mat_CI["1",])
	      pow_rCI <- conf_mat_rCI["1","1"]/sum(conf_mat_rCI["1",])
	      pow_pearson <- conf_mat_pearson["1","1"]/sum(conf_mat_pearson["1",])
	      pow_spearman <- conf_mat_spearman["1","1"]/sum(conf_mat_spearman["1",])
  	      pow_kendall <- conf_mat_kendall["1","1"]/sum(conf_mat_kendall["1",])

	      return(c("CI_power" = pow_CI, "rCI_power" = pow_rCI, "pearson_power" = pow_pearson, "spearman_power" = pow_spearman, "kendall_power" = pow_kendall))
	    })
	  return(powers)
	})

	n50_power <- abind(n50_power, along = -1)

	toPlot <- melt(n50_power[,,as.character(delta)])
	colnames(toPlot) <- c("Effect Size", "Method", "Power")

	pres_ready <- theme_bw() + 
	  theme(axis.title = element_text(size=24), axis.text = element_text(size=24), legend.text = element_text(size=24), title = element_text(size=28),legend.key.height = unit(1.0, 'cm'))

	pdf(fileName, height = 6, width=9)
	p <- ggplot(toPlot, aes(x=`Effect Size`, y=Power, col=Method)) + geom_line(size=1) + ggtitle(paste0("N=", sampleSize, " Delta = ",delta,", alpha = ", alpha)) + pres_ready
	print(p)	
	dev.off()



	# toPlot <- melt(n50_power[,,"0"])
	# colnames(toPlot) <- c("Effect Size", "Method", "Power")

	# library(ggplot2)
	# pdf("power_gaussian_with_pearson_delta_0.pdf", height = 6, width=9)
	# ggplot(toPlot, aes(x=`Effect Size`, y=Power, col=Method)) + geom_line(size=1.5) + ggtitle(paste0("N=50", sampleSize, " Delta = 0, alpha = ", alpha)) + pres_ready
	# dev.off()

}


#TODO:: actually implement this
MakePowerOverDeltaPlot <- function(power_res, sampleSize, effectSize, alpha = 0.001,  fileName){

	if(missing(fileName)){
		fileName = paste("power_plot_n", sampleSize, "effectSize", effectSize, ".pdf", sep="_")
	}
	# browser()
	n50 <- power_res[,as.character(sampleSize)]
	n50 <- lapply(n50, unlist, recursive=FALSE)
	n50 <- lapply(n50, unlist, recursive=FALSE)

	n50 <- lapply(n50, abind, along=-1)

	contingency_helper <- function(x1, x2){
	  
	  resTbl <- matrix(0, nrow=2, ncol=2)
	  dimnames(resTbl) <- list(c("0", "1"), c("0","1"))
	  resTbl["0","0"] <- sum(x1==0 & x2==0)
	  resTbl["0","1"] <- sum(x1==0 & x2==1)
	  resTbl["1","0"] <- sum(x1==1 & x2==0)
	  resTbl["1","1"] <- sum(x1==1 & x2==1)
	  return(resTbl)
	}

	n50_power <- lapply(n50, function(x) {
	  x[,,2:5] <- apply(x[,,2:5], c(1,2,3), function(x) (as.numeric(x < alpha)))
	  powers <- sapply(dimnames(x)[[2]], function(xx){

	      conf_mat_CI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"CI_p"])
	      conf_mat_rCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"rCI_p"])
	      conf_mat_pearson <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Pearson_p"])
	      conf_mat_spearman <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Spearman_p"])
	      # browser()
	      pow_CI <- conf_mat_CI["1","1"]/sum(conf_mat_CI["1",])
	      pow_rCI <- conf_mat_rCI["1","1"]/sum(conf_mat_rCI["1",])
	      pow_pearson <- conf_mat_pearson["1","1"]/sum(conf_mat_pearson["1",])
	      pow_spearman <- conf_mat_spearman["1","1"]/sum(conf_mat_spearman["1",])
	      
	      return(c("CI_power" = pow_CI, "rCI_power" = pow_rCI, "pearson_power" = pow_pearson, "spearman_power" = pow_spearman))
	    })
	  return(powers)
	})

	n50_power <- abind(n50_power, along = -1)
	# browser()
	toPlot <- melt(n50_power[as.character(effectSize),,])
	colnames(toPlot) <- c("Method", "Delta", "Power")

	pres_ready <- theme_bw() + 
	  theme(axis.title = element_text(size=24), axis.text = element_text(size=24), legend.text = element_text(size=24), title = element_text(size=28),legend.key.height = unit(1.0, 'cm'))

	pdf(fileName, height = 6, width=9)
	p <- ggplot(toPlot, aes(x=`Delta`, y=Power, col=Method)) + geom_line(size=1) + ggtitle(paste0("N=", sampleSize, " Effect Size = ",effectSize,", alpha = ", alpha)) + pres_ready
	print(p)	
	dev.off()

}


MakePowerOverDeltaPlotBeta <- function(power_res, sampleSize, effectSize, alpha = 0.001,  fileName){

	if(missing(fileName)){
		fileName = paste("power_plot_beta_n", sampleSize, "effectSize", effectSize, ".pdf", sep="_")
	}
	# browser()
	n50 <- power_res[,as.character(sampleSize)]
	# n50 <- lapply(n50, unlist, recursive=FALSE)
	# n50 <- lapply(n50, unlist, recursive=FALSE)

	n50 <- lapply(n50, abind, along=-1)

	contingency_helper <- function(x1, x2){
	  
	  resTbl <- matrix(0, nrow=2, ncol=2)
	  dimnames(resTbl) <- list(c("0", "1"), c("0","1"))
	  resTbl["0","0"] <- sum(x1==0 & x2==0)
	  resTbl["0","1"] <- sum(x1==0 & x2==1)
	  resTbl["1","0"] <- sum(x1==1 & x2==0)
	  resTbl["1","1"] <- sum(x1==1 & x2==1)
	  return(resTbl)
	}

	n50_power <- lapply(n50, function(x) {
	  x[,,2:6] <- apply(x[,,2:6], c(1,2,3), function(x) (as.numeric(x < alpha)))
	  powers <- sapply(dimnames(x)[[2]], function(xx){

	      conf_mat_CI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"CI_p"])
	      conf_mat_rCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"rCI_p"])
	      conf_mat_pearson <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Pearson_p"])
	      conf_mat_spearman <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Spearman_p"])
	      conf_mat_KCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"KCI_p"])

	      # browser()
	      pow_CI <- conf_mat_CI["1","1"]/sum(conf_mat_CI["1",])
	      pow_rCI <- conf_mat_rCI["1","1"]/sum(conf_mat_rCI["1",])
	      pow_pearson <- conf_mat_pearson["1","1"]/sum(conf_mat_pearson["1",])
	      pow_spearman <- conf_mat_spearman["1","1"]/sum(conf_mat_spearman["1",])
	      pow_KCI <- conf_mat_KCI["1","1"]/sum(conf_mat_KCI["1",])

	      return(c("CI" = pow_CI, "rCI" = pow_rCI, "Pearson" = pow_pearson, "Spearman" = pow_spearman, KCI = pow_KCI))
	    })
	  return(powers)
	})

	n50_power <- abind(n50_power, along = -1)
	# browser()
	toPlot <- melt(n50_power[as.character(effectSize),,])
	colnames(toPlot) <- c("Method", "Delta", "Power")

	pres_ready <- theme_bw() + 
	  theme(axis.title = element_text(size=24), axis.text = element_text(size=24), legend.text = element_text(size=24), title = element_text(size=28),legend.key.height = unit(1.0, 'cm'))

	pdf(fileName, height = 6, width=9)
	p <- ggplot(toPlot, aes(x=`Delta`, y=Power, col=Method)) + geom_line(size=1) + ggtitle(paste0("N=", sampleSize, " Effect Size = ",effectSize,", alpha = ", alpha)) + pres_ready
	print(p)	
	dev.off()

}


#TODO:: actually implement this
MakePowerOverDeltaPlotWithKendall <- function(power_res, sampleSize, effectSize, alpha = 0.001,  fileName){

	if(missing(fileName)){
		fileName = paste("power_plot_kendall_n", sampleSize, "effectSize", effectSize, ".pdf", sep="_")
	}
	# browser()
	n50 <- power_res[,as.character(sampleSize)]
	n50 <- lapply(n50, unlist, recursive=FALSE)
	n50 <- lapply(n50, unlist, recursive=FALSE)

	n50 <- lapply(n50, abind, along=-1)

	contingency_helper <- function(x1, x2){
	  
	  resTbl <- matrix(0, nrow=2, ncol=2)
	  dimnames(resTbl) <- list(c("0", "1"), c("0","1"))
	  resTbl["0","0"] <- sum(x1==0 & x2==0)
	  resTbl["0","1"] <- sum(x1==0 & x2==1)
	  resTbl["1","0"] <- sum(x1==1 & x2==0)
	  resTbl["1","1"] <- sum(x1==1 & x2==1)
	  return(resTbl)
	}

	n50_power <- lapply(n50, function(x) {
	  x[,,2:6] <- apply(x[,,2:6], c(1,2,3), function(x) (as.numeric(x < alpha)))
	  powers <- sapply(dimnames(x)[[2]], function(xx){

	      conf_mat_CI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"CI_p"])
	      conf_mat_rCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"rCI_p"])
	      conf_mat_pearson <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Pearson_p"])
	      conf_mat_spearman <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Spearman_p"])
  	      conf_mat_kendall <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Kendall_p"])

	      # browser()
	      pow_CI <- conf_mat_CI["1","1"]/sum(conf_mat_CI["1",])
	      pow_rCI <- conf_mat_rCI["1","1"]/sum(conf_mat_rCI["1",])
	      pow_pearson <- conf_mat_pearson["1","1"]/sum(conf_mat_pearson["1",])
	      pow_spearman <- conf_mat_spearman["1","1"]/sum(conf_mat_spearman["1",])
  	      pow_kendall <- conf_mat_kendall["1","1"]/sum(conf_mat_kendall["1",])

	      return(c("CI_power" = pow_CI, "rCI_power" = pow_rCI, "pearson_power" = pow_pearson, "spearman_power" = pow_spearman, "kendall_power" = pow_kendall))
	    })
	  return(powers)
	})

	n50_power <- abind(n50_power, along = -1)
	# browser()
	toPlot <- melt(n50_power[as.character(effectSize),,])
	colnames(toPlot) <- c("Method", "Delta", "Power")

	pres_ready <- theme_bw() + 
	  theme(axis.title = element_text(size=24), axis.text = element_text(size=24), legend.text = element_text(size=24), title = element_text(size=28),legend.key.height = unit(1.0, 'cm'))

	pdf(fileName, height = 6, width=9)
	p <- ggplot(toPlot, aes(x=`Delta`, y=Power, col=Method)) + geom_line(size=1) + ggtitle(paste0("N=", sampleSize, " Effect Size = ",effectSize,", alpha = ", alpha)) + pres_ready
	print(p)	
	dev.off()

}

#TODO:: actually implement this
MakePowerOverDeltaColorEffectPlot <- function(power_res, sampleSize, effectSize, alpha = 0.001,  fileName){

	if(missing(fileName)){
		fileName = paste("power_plot_n", sampleSize, "colorEffect", ".pdf", sep="_")
	}
	# browser()
	n50 <- power_res[,as.character(sampleSize)]
	n50 <- lapply(n50, unlist, recursive=FALSE)
	n50 <- lapply(n50, unlist, recursive=FALSE)

	n50 <- lapply(n50, abind, along=-1)

	contingency_helper <- function(x1, x2){
	  
	  resTbl <- matrix(0, nrow=2, ncol=2)
	  dimnames(resTbl) <- list(c("0", "1"), c("0","1"))
	  resTbl["0","0"] <- sum(x1==0 & x2==0)
	  resTbl["0","1"] <- sum(x1==0 & x2==1)
	  resTbl["1","0"] <- sum(x1==1 & x2==0)
	  resTbl["1","1"] <- sum(x1==1 & x2==1)
	  return(resTbl)
	}

	n50_power <- lapply(n50, function(x) {
	  x[,,2:5] <- apply(x[,,2:5], c(1,2,3), function(x) (as.numeric(x < alpha)))
	  powers <- sapply(dimnames(x)[[2]], function(xx){

	      conf_mat_CI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"CI_p"])
	      conf_mat_rCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"rCI_p"])
	      conf_mat_pearson <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Pearson_p"])
	      conf_mat_spearman <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Spearman_p"])
	      # browser()
	      pow_CI <- conf_mat_CI["1","1"]/sum(conf_mat_CI["1",])
	      pow_rCI <- conf_mat_rCI["1","1"]/sum(conf_mat_rCI["1",])
	      pow_pearson <- conf_mat_pearson["1","1"]/sum(conf_mat_pearson["1",])
	      pow_spearman <- conf_mat_spearman["1","1"]/sum(conf_mat_spearman["1",])
	      
	      return(c("CI_power" = pow_CI, "rCI_power" = pow_rCI, "pearson_power" = pow_pearson, "spearman_power" = pow_spearman))
	    })
	  return(powers)
	})

	n50_power <- abind(n50_power, along = -1)
	# browser()
	toPlot <- melt(n50_power[,"rCI_power",])
	# browser()
	colnames(toPlot) <- c("Effect Size", "Delta", "Power")
	toPlot[,"Effect Size"] <- as.character(toPlot[,"Effect Size"])
	toPlot <- data.table(toPlot)
	toPlot[,Power := Power/max(Power), .(`Effect Size`)]
	pres_ready <- theme_bw() + 
	  theme(axis.title = element_text(size=24), axis.text = element_text(size=24), legend.text = element_text(size=24), title = element_text(size=28),legend.key.height = unit(1.0, 'cm'))

	pdf(fileName, height = 6, width=9)
	p <- ggplot(toPlot, aes(x=`Delta`, y=Power, col=`Effect Size`)) + geom_line(size=1) + ggtitle(paste0("N=", sampleSize, ", alpha = ", alpha)) + pres_ready + ylab("Power as % of max power")
	print(p)	
	dev.off()


}


#TODO:: actually implement this
MakePowerOverDeltaBetaColorEffectPlot <- function(power_res, sampleSize, effectSize, alpha = 0.001,  fileName){

	if(missing(fileName)){
		fileName = paste("power_plot_beta_n", sampleSize, "colorEffect", ".pdf", sep="_")
	}
	# browser()
	n50 <- power_res[,as.character(sampleSize)]
	# n50 <- lapply(n50, unlist, recursive=FALSE)
	# n50 <- lapply(n50, unlist, recursive=FALSE)

	n50 <- lapply(n50, abind, along=-1)

	contingency_helper <- function(x1, x2){
	  
	  resTbl <- matrix(0, nrow=2, ncol=2)
	  dimnames(resTbl) <- list(c("0", "1"), c("0","1"))
	  resTbl["0","0"] <- sum(x1==0 & x2==0)
	  resTbl["0","1"] <- sum(x1==0 & x2==1)
	  resTbl["1","0"] <- sum(x1==1 & x2==0)
	  resTbl["1","1"] <- sum(x1==1 & x2==1)
	  return(resTbl)
	}

	n50_power <- lapply(n50, function(x) {
	  x[,,2:6] <- apply(x[,,2:6,drop=FALSE], c(1,2,3), function(x) (as.numeric(x < alpha)))
	  powers <- sapply(dimnames(x)[[2]], function(xx){

	      conf_mat_CI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"CI_p"])
	      conf_mat_rCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"rCI_p"])
	      conf_mat_pearson <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Pearson_p"])
	      conf_mat_spearman <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Spearman_p"])
	      conf_mat_KCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"KCI_p"])

	      # browser()
	      pow_CI <- conf_mat_CI["1","1"]/sum(conf_mat_CI["1",])
	      pow_rCI <- conf_mat_rCI["1","1"]/sum(conf_mat_rCI["1",])
	      pow_pearson <- conf_mat_pearson["1","1"]/sum(conf_mat_pearson["1",])
	      pow_spearman <- conf_mat_spearman["1","1"]/sum(conf_mat_spearman["1",])
	      pow_KCI <- conf_mat_KCI["1","1"]/sum(conf_mat_KCI["1",])

	      return(c("CI" = pow_CI, "rCI" = pow_rCI, "Pearson" = pow_pearson, "Spearman" = pow_spearman, KCI = pow_KCI))
	    })
	  return(powers)
	})

	n50_power <- abind(n50_power, along = -1)
	# browser()
	toPlot <- melt(n50_power[,"rCI",])
	# browser()
	colnames(toPlot) <- c("Effect Size", "Delta", "Power")
	toPlot[,"Effect Size"] <- as.character(toPlot[,"Effect Size"])
	toPlot <- data.table(toPlot)
	toPlot[,Power := Power/max(Power), .(`Effect Size`)]
	pres_ready <- theme_bw() + 
	  theme(axis.title = element_text(size=24), axis.text = element_text(size=24), legend.text = element_text(size=24), title = element_text(size=28),legend.key.height = unit(1.0, 'cm'))

	pdf(fileName, height = 6, width=9)
	p <- ggplot(toPlot, aes(x=`Delta`, y=Power, col=`Effect Size`)) + geom_line(size=1) + ggtitle(paste0("N=", sampleSize, ", alpha = ", alpha)) + pres_ready + ylab("Power as % of max power")
	print(p)	
	dev.off()


}



#TODO:: actually implement this
MakePowerLevelSet <- function(power_res, sampleSizes = c(50,100,150,200,250,300), 
							  powerLevel = 0.5, alpha = 0.001, plotType= c("percent", "nominal"), fileName){
	require(pwr)
	plotType <- match.arg(plotType)
	if(missing(fileName)){
		fileName = paste("power_plot_levelset", plotType, ".pdf", sep="_")
	}


	effect_sizes <-  sapply(sampleSizes, function(n) {
		pwr.r.test(n, sig.level=alpha, power=powerLevel)$r
	})

	n50 <- lapply(power_res, unlist, recursive=FALSE)
	n50 <- lapply(n50, unlist, recursive=FALSE)

	n50 <- lapply(n50, abind, along=-1)

	contingency_helper <- function(x1, x2){
	  
	  resTbl <- matrix(0, nrow=2, ncol=2)
	  dimnames(resTbl) <- list(c("0", "1"), c("0","1"))
	  resTbl["0","0"] <- sum(x1==0 & x2==0)
	  resTbl["0","1"] <- sum(x1==0 & x2==1)
	  resTbl["1","0"] <- sum(x1==1 & x2==0)
	  resTbl["1","1"] <- sum(x1==1 & x2==1)
	  return(resTbl)
	}

	n50_power <- lapply(n50, function(x) {
	  x[,,2:5] <- apply(x[,,2:5, drop=FALSE], c(1,2,3), function(x) (as.numeric(x < alpha)))
	  powers <- sapply(dimnames(x)[[2]], function(xx){

	      conf_mat_CI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"CI_p"])
	      conf_mat_rCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"rCI_p"])
	      conf_mat_pearson <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Pearson_p"])
	      conf_mat_spearman <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Spearman_p"])
	      # browser()
	      pow_CI <- conf_mat_CI["1","1"]/sum(conf_mat_CI["1",])
	      pow_rCI <- conf_mat_rCI["1","1"]/sum(conf_mat_rCI["1",])
	      pow_pearson <- conf_mat_pearson["1","1"]/sum(conf_mat_pearson["1",])
	      pow_spearman <- conf_mat_spearman["1","1"]/sum(conf_mat_spearman["1",])
	      
	      return(c("CI_power" = pow_CI, "rCI_power" = pow_rCI, "pearson_power" = pow_pearson, "spearman_power" = pow_spearman))
	    })
	  return(powers)
	})

	n50_power <- abind(n50_power, along = -1)
	if(plotType == "percent"){
		n50_power <- n50_power/n50_power[,"pearson_power",]
	}
	# browser()
	toPlot <- melt(n50_power)
	# browser()
	colnames(toPlot) <- c("Effect Size", "Method", "Sample Size", "Power")
	toPlot[,"Effect Size"] <- rep(effect_sizes, ncol(n50_power))
	toPlot[,"Sample Size"] <- rep(sampleSizes, ncol(n50_power))


	pres_ready <- theme_bw() + 
	  theme(axis.title = element_text(size=24), axis.text = element_text(size=24), legend.text = element_text(size=24), title = element_text(size=28),legend.key.height = unit(1.0, 'cm'))
	
	pdf(fileName, height = 6, width=9)
	p <- ggplot(toPlot, aes(x=`Sample Size`, y=Power, col=`Method`)) + geom_line(size=1) + ggtitle(paste0("Power = ", powerLevel, ", Delta = 1, alpha = ", alpha)) + pres_ready 
	if(plotType == "percent"){
		p <- p + ylab("% of Pearson Power") 
	}
	print(p)	
	dev.off()


}





MakePowerLevelSetBeta <- function(power_res, sampleSizes = c(100,150,200,250,300), 
							  powerLevel = 0.5, alpha = 0.001, plotType= c("percent", "nominal"), fileName){
	require(pwr)
	plotType <- match.arg(plotType)
	if(missing(fileName)){
		fileName = paste("power_plot_beta_levelset", plotType, ".pdf", sep="_")
	}


	effect_sizes <-  sapply(sampleSizes, function(n) {
		pwr.r.test(n, sig.level=alpha, power=powerLevel)$r
	})
	# browser()
	# n50 <- lapply(power_res, unlist, recursive=FALSE)
	# n50 <- lapply(n50, unlist, recursive=FALSE)

	n50 <- lapply(power_res, abind, along=-1)
	# browser()
	contingency_helper <- function(x1, x2){
	  
	  resTbl <- matrix(0, nrow=2, ncol=2)
	  dimnames(resTbl) <- list(c("0", "1"), c("0","1"))
	  resTbl["0","0"] <- sum(x1==0 & x2==0)
	  resTbl["0","1"] <- sum(x1==0 & x2==1)
	  resTbl["1","0"] <- sum(x1==1 & x2==0)
	  resTbl["1","1"] <- sum(x1==1 & x2==1)
	  return(resTbl)
	}

	n50_power <- lapply(n50, function(x) {
	  x[,,2:6] <- apply(x[,,2:6, drop=FALSE], c(1,2,3), function(x) (as.numeric(x < alpha)))
	  powers <- sapply(dimnames(x)[[2]], function(xx){

	      conf_mat_CI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"CI_p"])
	      conf_mat_rCI <- contingency_helper(x[,xx,"Alternative"],x[,xx,"rCI_p"])
	      conf_mat_pearson <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Pearson_p"])
	      conf_mat_spearman <- contingency_helper(x[,xx,"Alternative"],x[,xx,"Spearman_p"])
	      # browser()
	      pow_CI <- conf_mat_CI["1","1"]/sum(conf_mat_CI["1",])
	      pow_rCI <- conf_mat_rCI["1","1"]/sum(conf_mat_rCI["1",])
	      pow_pearson <- conf_mat_pearson["1","1"]/sum(conf_mat_pearson["1",])
	      pow_spearman <- conf_mat_spearman["1","1"]/sum(conf_mat_spearman["1",])
	      
	      return(c("CI_power" = pow_CI, "rCI_power" = pow_rCI, "pearson_power" = pow_pearson, "spearman_power" = pow_spearman))
	    })
	  return(powers)
	})

	n50_power <- abind(n50_power, along = -1)
	if(plotType == "percent"){
		n50_power <- n50_power/n50_power[,"pearson_power",]
	}
	# browser()
	toPlot <- melt(n50_power)
	# browser()
	colnames(toPlot) <- c("Effect Size", "Method", "Sample Size", "Power")
	toPlot[,"Effect Size"] <- rep(effect_sizes, ncol(n50_power))
	toPlot[,"Sample Size"] <- rep(sampleSizes, ncol(n50_power))


	pres_ready <- theme_bw() + 
	  theme(axis.title = element_text(size=24), axis.text = element_text(size=24), legend.text = element_text(size=24), title = element_text(size=28),legend.key.height = unit(1.0, 'cm'))
	
	pdf(fileName, height = 6, width=9)
	p <- ggplot(toPlot, aes(x=`Sample Size`, y=Power, col=`Method`)) + geom_line(size=1) + ggtitle(paste0("Power = ", powerLevel, ", Delta = 1, alpha = ", alpha)) + pres_ready 
	if(plotType == "percent"){
		p <- p + ylab("% of Pearson Power") 
	}
	print(p)	
	dev.off()


}




