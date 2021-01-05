library(abind)
library(data.table)
library(ggplot2)
library(reshape2)
library(ggplot2)

MakePowerOverEffectPlot <- function(power_res, sampleSize, alpha = 0.001, delta=1, fileName){

	if(missing(fileName)){
		fileName = paste("power_plot_n", sampleSize, "delta", delta, ".pdf", sep="_")
	}

	n50 <- c(list_mat[,as.character(sampleSize)], power_res[,as.character(sampleSize)])
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

	n50 <- c(list_mat[,as.character(sampleSize)], power_res[,as.character(sampleSize)])
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



	# toPlot <- melt(n50_power[,,"0"])
	# colnames(toPlot) <- c("Effect Size", "Method", "Power")

	# library(ggplot2)
	# pdf("power_gaussian_with_pearson_delta_0.pdf", height = 6, width=9)
	# ggplot(toPlot, aes(x=`Effect Size`, y=Power, col=Method)) + geom_line(size=1.5) + ggtitle(paste0("N=50", sampleSize, " Delta = 0, alpha = ", alpha)) + pres_ready
	# dev.off()

}
