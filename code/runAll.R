library(ggplot2)
source("/code/plotPowerRes.R")


## Output figure 3

load("/data/results/normal_dist_power_analysis_rho_n_100_deltafinemap.RData")
MakePowerOverDeltaPlot(list_mat, effectSize = 0.3, sampleSize = 100, fileName = "/results/figure_3a.pdf")
MakePowerOverDeltaColorEffectPlot(list_mat, sampleSize = 100,fileName = "/results/figure_3b.pdf")



load("/data/results/normal_dist_power_analysis_n_100_withkendall.RData")
MakePowerOverEffectPlot(list_mat, 100, fileName = "/results/figure_3c.pdf")



load("/data/results/normal_dist_power_analysis_delta_1_0_fixed_pearson_power_0_5_large.RData")
MakePowerLevelSet(list_mat, fileName = "/results/figure_3d.pdf")



## Output supplementary figure 3 
source("/code/plotBivariateBeta.R")

load("/data/results/beta_1_2__4_5_dist_power_analysis_withkci.RData")

MakePowerOverEffectPlotBeta(list_mat, 100, delta=0.1,fileName = "/results/figure_s3b.pdf")


load("/data/results/beta_dist_power_analysis_delta_0_1_fixed_pearson_power_0_5_large.RData")
MakePowerLevelSetBeta(list_mat, fileName = "/results/figure_s3c.pdf")

  load("/data/results/beta_dist_power_analysis_delta_0_1_fixed_pearson_power_0_25_large.RData")
MakePowerLevelSetBeta(list_mat, powerLevel = 0.25, fileName = "/results/figure_s3d.pdf")

load("/data/results/beta_1_2__4_5_dist_power_analysis_withkci_withnoise.RData")

MakePowerOverEffectPlotBeta(list_mat, 100, delta=0.1,fileName = "/results/figure_s4b.pdf")
