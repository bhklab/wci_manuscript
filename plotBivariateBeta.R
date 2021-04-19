library(ggExtra)

source("bivariateBeta.R")


samplePars <- optimizeShapeForPearson(0.7, c(1.2,4.5), c(1.2,4.5))


test <- data.frame(rbivariateBeta(100000, samplePars))


test$col = densCols(test$x, test$y, colramp = colorRampPalette(c("#2b8cbe", "#00FEFF", "#45FE4F",
                            "#FCFF00", "#FF9400", "#FF3100")))

p <- ggplot(as.data.frame(test), aes(x,y, color=col)) + geom_point(alpha=0.5) + theme_bw() + 
theme(axis.text = element_text(size=14), axis.title=element_text(size=18), legend.position="none")

pdf("bivariateBeta_1_2__4_5.pdf", width=5, height=5)
ggMarginal(p)
dev.off()

