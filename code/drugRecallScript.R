
# Drug recall analysis for RCI manuscript
# This function takes as input rds outputs produced by get_sim_mat (in drug_recall.R) and summarizes
# the results. Generates Figures 4a and 4b. 

drug_recall_analysis <- function(resdir=".", outdir="."){
  results <- list(pearson_rk=list(), spearman_rk=list(), ci_rk=list(), rci_rk=list(), kci_rk=list())
  
  f <- list.files(resdir, pattern=".rds")
  for (myfile in f){
    print(myfile)
    z <- readRDS(file.path(resdir, myfile))
    
    for (metric in seq(5)){
      mymat <- z[[metric]]$simmat
      trank <- numeric(dim(mymat)[1])
      
      for (jj in seq(dim(mymat)[1])){
        trank[jj] <- (sum(!is.na(mymat[jj,])) - rank(mymat[jj,], na.last="keep")[which(colnames(mymat) == rownames(mymat)[jj])])/(sum(!is.na(mymat[jj,])) - 1)
      }
      
      results[[metric]] <- c(results[[metric]], list(trank))
    }
  }
  
  cols <- rainbow(5)
  
  pdf(file=file.path(outdir, "fig4a_drugrecall_ecdf.pdf"), width=8, height=6.5)
  plot(ecdf(unlist(results$pearson_rk)), xlim=c(0,1), ylim=c(0,1), xlab="Rank", ylab="Cumulative Fraction", col=cols[1], main="Drug Recall Analysis, all drugs")
  lines(ecdf(unlist(results$spearman_rk)), col=cols[2], type="b")
  lines(ecdf(unlist(results$ci_rk)), col=cols[3], type="b")
  lines(ecdf(unlist(results$rci_rk)), col=cols[4], type="b")
  lines(ecdf(unlist(results$kci_rk)), col=cols[5], type="b")
  legend(x=0.8, y=0.25, legend=c("Pearson", "Spearman", "CI", "rCI", "kCI"), lwd=c(2,2,2,2,2), col=cols)
  dev.off()
  
  # Slice with mean AAC and N cell lines:
  cellcounts <- list()
  f <- list.files(resdir, pattern = ".rds")
  for (myfile in f){
    print(myfile)
    z <- readRDS(file.path(resdir, myfile))
    
    mycts <- z[[1]]$cellcounts
    trank <- numeric(dim(mycts)[1])
      
    for (jj in seq(dim(mycts)[1])){
      trank[jj] <- mycts[jj, which(colnames(mycts) == rownames(mycts)[jj])]
    }
      
    cellcounts <- c(cellcounts, list(trank))
  }
  
  cellcountmat <- list()
  f <- list.files(resdir, pattern = ".rds")
  for (myfile in f){
    print(myfile)
    z <- readRDS(file.path(resdir, myfile))
    
    mycts <- z[[1]]$cellcounts
    cellcountmat <- c(cellcountmat, list(mycts))
  }
  
  cellcountnorm <- list()
  for (mymat in cellcountmat){
    
  }
  
  pdf(file=file.path(outdir, "fig4b_drugrecall_ecdf.pdf"), width=8, height=6.5)
  plot(ecdf(unlist(results$pearson_rk)[unlist(cellcounts) > 50]), xlim=c(0,1), ylim=c(0,1), xlab="Rank", ylab="Cumulative Fraction", col=cols[1], main="Drug Recall Analysis, N > 50 cell lines")
  lines(ecdf(unlist(results$spearman_rk)[unlist(cellcounts) > 50]), col=cols[2], type="b")
  lines(ecdf(unlist(results$ci_rk)[unlist(cellcounts) > 50]), col=cols[3], type="b")
  lines(ecdf(unlist(results$rci_rk)[unlist(cellcounts) > 50]), col=cols[4], type="b")
  lines(ecdf(unlist(results$kci_rk)[unlist(cellcounts) > 50]), col=cols[5], type="b")
  legend(x=0.8, y=0.25, legend=c("Pearson", "Spearman", "CI", "rCI", "kCI"), lwd=c(2,2,2,2,2), col=cols)
  dev.off()

  
  return(results, cellcounts)
}
