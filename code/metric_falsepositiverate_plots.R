library(ggplot2)

labels <- c("beta_1-2_4-5", "norm_sd=1")
iter <- 200000
npts_vec <- c(25, 40, 50, 100, 200, 300, 400, 600, 800, 1000)
kci_range <- c(25, 40, 50, 100, 200, 300, 400, 600, 800)
kci_iter <- 100000

outdir <- "/results"

for (label in labels){
  # Load results computed by metric_falsepositiverate.R
  
  # Read saved KCI results
  xk <- readRDS(file=file.path("data/results", sprintf("%s_n=%dk_kci_n=%d_%d.RDS", label, kci_iter/1000, min(kci_range), max(kci_range))))
  
  
  summary_rates <- data.frame(npts=numeric(0), metric=character(0), p_lt_4=numeric(0), p_lt_3=numeric(0), p_lt_2=numeric(0))
  for (npts in npts_vec){
    x <- readRDS(file.path("data/results", sprintf("%s_n=%dk_allmetrics_rcivariants_npts=%d.RDS", label, iter/1000, npts)))
    
    #### Figures ####
    pearsp <- unlist(x$pearson_vals$p.value)
    spearp <- unlist(x$spearman_vals$p.value)
    kendallp <- unlist(x$kendall_vals$p.value)
    cip <- unlist(x$ci_vals$p.value)
    
    rcip <- unlist(x$rci_vals_20$p.value)
    rcipairs <- unlist(x$rci_vals_20$relevant.pairs.no)
    rcip_10 <- unlist(x$rci_vals_10$p.value)
    rcipairs_10 <- unlist(x$rci_vals_10$relevant.pairs.no)
    rcip_05 <- unlist(x$rci_vals_05$p.value)
    rcipairs_05 <- unlist(x$rci_vals_05$relevant.pairs.no)
    
    maxy <- min(12, max(-log10(pearsp), -log10(spearp), -log10(cip), -log10(rcip)), na.rm = TRUE)
    
    summary_rates <- rbind(summary_rates, 
                           data.frame(npts=npts, metric="pearson", p_lt_4=mean(pearsp < 1e-4), p_lt_3=mean(pearsp < 1e-3), p_lt_2=mean(pearsp < 1e-2)))
    summary_rates <- rbind(summary_rates, 
                           data.frame(npts=npts, metric="spearman", p_lt_4=mean(spearp < 1e-4), p_lt_3=mean(spearp < 1e-3), p_lt_2=mean(spearp < 1e-2)))
    summary_rates <- rbind(summary_rates, 
                           data.frame(npts=npts, metric="kendall", p_lt_4=mean(spearp < 1e-4), p_lt_3=mean(spearp < 1e-3), p_lt_2=mean(spearp < 1e-2)))
    summary_rates <- rbind(summary_rates, 
                           data.frame(npts=npts, metric="ci", p_lt_4=mean(cip < 1e-4), p_lt_3=mean(cip < 1e-3), p_lt_2=mean(cip < 1e-2)))
    summary_rates <- rbind(summary_rates, 
                           data.frame(npts=npts, metric="rci_d=0.2", p_lt_4=mean(rcip < 1e-4), p_lt_3=mean(rcip < 1e-3), p_lt_2=mean(rcip < 1e-2)))
    summary_rates <- rbind(summary_rates, 
                           data.frame(npts=npts, metric="rci_d=0.1", p_lt_4=mean(rcip_10 < 1e-4), p_lt_3=mean(rcip_10 < 1e-3), p_lt_2=mean(rcip_10 < 1e-2)))
    summary_rates <- rbind(summary_rates, 
                           data.frame(npts=npts, metric="rci_d=0.05", p_lt_4=mean(rcip_05 < 1e-4), p_lt_3=mean(rcip_05 < 1e-3), p_lt_2=mean(rcip_05 < 1e-2)))
    summary_rates <- rbind(summary_rates, 
                           data.frame(npts=npts, metric="rci_d=0.2_pairs>3k", p_lt_4=mean(rcip[rcipairs > 3000] < 1e-4), 
                                      p_lt_3=mean(rcip[rcipairs > 3000] < 1e-3), p_lt_2=mean(rcip[rcipairs > 3000] < 1e-2)))
    summary_rates <- rbind(summary_rates,
                           data.frame(npts=npts, metric="rci_d=0.1_pairs>3k", p_lt_4=mean(rcip_10[rcipairs_10 > 3000] < 1e-4), 
                                      p_lt_3=mean(rcip_10[rcipairs_10 > 3000] < 1e-3), p_lt_2=mean(rcip_10[rcipairs_10 > 3000] < 1e-2)))
    summary_rates <- rbind(summary_rates, 
                           data.frame(npts=npts, metric="rci_d=0.05_pairs>3k", p_lt_4=mean(rcip_05[rcipairs_05 > 3000] < 1e-4), 
                                      p_lt_3=mean(rcip_05[rcipairs_05 > 3000] < 1e-3), p_lt_2=mean(rcip_05[rcipairs_05 > 3000] < 1e-2)))
    
    pdf(sprintf("false_positives_%s_rcivariants_n=%dk_npts=%d.pdf", label, iter/1000, npts), width=10, height=8)
    plot(-log10(rank(pearsp)[rank(pearsp) < 0.1*length(pearsp)]/length(pearsp)), -log10(pearsp)[rank(pearsp) < 0.1*length(pearsp)], xlab="-log10 Empirical Quantile (Rank)",
         ylab="-log10 Theoretical Quantile (computed p-value)", main=sprintf("QQ Plot of p-values of permutations of %s, n=%dk; length=%d", dsname, iter/1000, npts),
         xlim=c(0, maxy), ylim=c(0, maxy), col="black", pch=16)
    grid()
    points(-log10(rank(spearp)[rank(spearp) < 0.1*length(spearp)]/length(spearp)), -log10(spearp)[rank(spearp) < 0.1*length(spearp)], pch=16, col="green")
    points(-log10(rank(kendallp)[rank(kendallp) < 0.1*length(kendallp)]/length(kendallp)), -log10(kendallp)[rank(kendallp) < 0.1*length(kendallp)], pch=16, col="purple")
    points(-log10(rank(cip)[rank(cip) < 0.1*length(cip)]/length(cip)), -log10(cip)[rank(cip) < 0.1*length(cip)], pch=16, col="blue")
    points(-log10(rank(rcip)[rank(rcip) < 0.1*length(rcip)]/length(rcip)), -log10(rcip)[rank(rcip) < 0.1*length(rcip)], pch=16, col="red")
    points(-log10(rank(rcip_10)[rank(rcip_10) < 0.1*length(rcip_10)]/length(rcip_10)), -log10(rcip_10)[rank(rcip_10) < 0.1*length(rcip_10)], pch=16, col="orangered")
    points(-log10(rank(rcip_05)[rank(rcip_05) < 0.1*length(rcip_05)]/length(rcip_05)), -log10(rcip_05)[rank(rcip_05) < 0.1*length(rcip_05)], pch=16, col="orange")
    lines(c(0,maxy), c(0,maxy), lty=2, col="grey")
    legend(x="bottomright", legend=c("Pearson", "Spearman", "Kendall", "CI", "RCI d=0.2", "RCI d=0.1", "RCI d=0.05"),
           col=c("black", "green", "purple", "blue", "red", "orangered", "orange"), pch=16)
    dev.off()
  }
  
  if (file.exists (sprintf("%s_n=%dk_kci_n=%d_%d.RDS", label, kci_iter/1000, min(kci_range), max(kci_range)), sep="/")){
    for (ii in kci_range){
      kcip <- xk[[sprintf("kci_npts_%d", ii)]]$p.value
      kcipairs <- xk[[sprintf("kci_npts_%d", ii)]]
      summary_rates <- rbind(summary_rates, 
                             data.frame(npts=ii, metric="kci", p_lt_4=mean(kcip < 1e-4), p_lt_3=mean(kcip < 1e-3), p_lt_2=mean(kcip < 1e-2)))
    }
  }
  
  ggplot(summary_rates[!(summary_rates$metric %in% c("rci_d=0.2_pairs>3k", "rci_d=0.1_pairs>3k", "rci_d=0.05_pairs>3k")), ], 
         aes(x=npts, y=p_lt_3, group=metric)) + geom_line(aes(color=metric), size=2, alpha=0.6) + geom_point(aes(fill=metric), pch=21, color="black", size=3) + 
    theme_minimal() + ggtitle(sprintf("False Positive Rates for p < 1e-3, %d random draws from %s", iter, dsname)) + xlab("N points") +
    geom_hline(yintercept=1e-3, linetype="dashed", color="black", size=1) + xlim(c(0, 1000)) + 
    ylim(c(0, max(summary_rates$p_lt_3[summary_rates$p_lt_3 < 0.2], na.rm = 1))) + ylab("Fraction with p < 1e-3)")
  dev.print(pdf, paste(outdir, sprintf("summary_falsepositives_%s_rcivariants_n=%dk_p=3.pdf", label, iter/1000), sep="/"), width=10, height=8)
  
  ggplot(summary_rates[!(summary_rates$metric %in% c("rci_d=0.2_pairs>3k", "rci_d=0.1_pairs>3k", "rci_d=0.05_pairs>3k")), ], 
         aes(x=npts, y=p_lt_4, group=metric)) + geom_line(aes(color=metric), size=2, alpha=0.6) + geom_point(aes(fill=metric), pch=21, color="black", size=3) + 
    theme_minimal() + ggtitle(sprintf("False Positive Rates at p < 1e-4, %d random draws from %s", iter, dsname)) + xlab("N points") +
    geom_hline(yintercept=1e-4, linetype="dashed", color="black", size=1) + xlim(c(0, 1000)) + 
    ylim(c(0, max(summary_rates$p_lt_4[summary_rates$p_lt_3 < 0.2], na.rm = 1))) + ylab("Fraction with p < 1e-4)")
  dev.print(pdf, paste(outdir, sprintf("summary_falsepositives_%s_rcivariants_n=%dk_p=4.pdf", label, iter/1000), sep="/"), width=10, height=8)
  
  ggplot(summary_rates[!(summary_rates$metric %in% c("rci_d=0.2_pairs>3k", "rci_d=0.1_pairs>3k", "rci_d=0.05_pairs>3k")), ], 
         aes(x=npts, y=p_lt_2, group=metric)) + geom_line(aes(color=metric), size=2, alpha=0.6) + geom_point(aes(fill=metric), pch=21, color="black", size=3) + 
    theme_minimal() + ggtitle(sprintf("False Positive Rates at p < 1e-2, %d random draws from %s", iter, dsname)) + xlab("N points") + 
    geom_hline(yintercept=1e-2, linetype="dashed", color="black", size=1) + xlim(0, 1000) + ylab("Fraction with p < 1e-2")
  dev.print(pdf, paste(outdir, sprintf("summary_falsepositives_%s_rcivariants_n=%dk_p=2.pdf", label, iter/1000), sep="/"), width=10, height=8)
}
