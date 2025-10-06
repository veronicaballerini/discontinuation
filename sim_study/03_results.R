################################################################################ 
# Evaluating causal effects on time-to-event outcomes in an RCT in oncology
# with treatment discontinuation 
#
# Authors: V. Ballerini, B. Bornkamp, F .Mealli, C. Wang, Y. Zhang, A. Mattei
#
# Replication Package for results in Section 5 
# In particular: Figures 7, 8, 9, 10, 11, 12
#
# Code authors: Veronica Ballerini, Alessandra Mattei 
# Last modified: October 6, 2025 
################################################################################ 

# Optional: clean the envrionment 
# rm(list=ls())

# set the working directory: change here you wd
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

# process the results

# choose whether you want to reproduce the full simulation (full <- TRUE) study 
# or the intermediate one (full <- FALSE)
full <- FALSE

source("tables.R")  # this script contains tables in the Appendix

# upload libraries and other settings for plots
library(ggplot2)
library(RColorBrewer)

okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", 
               "#D55E00", "#CC79A7", "#000000")

############ DISTRIBUTIONAL CAUSAL EFFECTS ############################

if(full == TRUE){
  
####### COVERAGE: Figure 7 of the manuscript
coverage <- data.frame(values = rbind(as.matrix(coverage_dceitt1, ncol = 1), 
                                      as.matrix(coverage_dcend1, ncol = 1), 
                                      as.matrix(coverage_dceD1, ncol = 1), 
                                      as.matrix(coverage_dced1[1, ], ncol = 1), 
                                      as.matrix(coverage_dced1[2, ], ncol = 1), 
                                      as.matrix(coverage_dced1[3, ], ncol = 1), 
                                      as.matrix(coverage_dced1[4, ], ncol = 1), 
                                      as.matrix(coverage_dceitt1_nocov, ncol = 1), 
                                      as.matrix(coverage_dcend1_nocov, ncol = 1), 
                                      as.matrix(coverage_dceD1_nocov, ncol = 1), 
                                      as.matrix(coverage_dced1_nocov[1, ], ncol = 1), 
                                      as.matrix(coverage_dced1_nocov[2, ], ncol = 1), 
                                      as.matrix(coverage_dced1_nocov[3, ], ncol = 1), 
                                      as.matrix(coverage_dced1_nocov[4, ], ncol = 1), 
                                      as.matrix(coverage_dceitt2, ncol = 1), 
                                      as.matrix(coverage_dcend2, ncol = 1), 
                                      as.matrix(coverage_dceD2, ncol = 1), 
                                      as.matrix(coverage_dced2[1, ], ncol = 1), 
                                      as.matrix(coverage_dced2[2, ], ncol = 1), 
                                      as.matrix(coverage_dced2[3, ], ncol = 1), 
                                      as.matrix(coverage_dced2[4, ], ncol = 1), 
                                      as.matrix(coverage_dceitt2_nocov, ncol = 1), 
                                      as.matrix(coverage_dcend2_nocov, ncol = 1), 
                                      as.matrix(coverage_dceD2_nocov, ncol = 1), 
                                      as.matrix(coverage_dced2_nocov[1, ], ncol = 1), 
                                      as.matrix(coverage_dced2_nocov[2, ], ncol = 1), 
                                      as.matrix(coverage_dced2_nocov[3, ], ncol = 1), 
                                      as.matrix(coverage_dced2_nocov[4, ], ncol = 1)), 
                       scenario = rep(c("1", "1 w/o covariates", "2", "2 w/o covariates"), 
                                      each = 20 * 7), 
                       effect = rep(c("ITT", "ND", "D", "d=1", "d=2", "d=3", "d=4"),
    each = 20))

coverage$effect <- factor(coverage$effect, 
                          levels = c("ITT", "ND", "D", "d=1", "d=2", "d=3", "d=4", "d=5"))

  jpeg("plots/FULL/Figure7.jpeg", 
       width = 7000, height = 5000, units = "px", res = 1200)
  print(ggplot(coverage, aes(x = scenario, y = values, fill = effect)) + 
    labs(x = "Scenario", y = "Coverage", fill = NULL) + geom_boxplot() + 
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") + 
    scale_fill_manual(values = okabe_ito[1:length(unique(coverage$effect))]) + 
    theme_minimal())
  dev.off()
  
  ####### WIDTH: Figure 8 of the manuscript
  avg_width <- data.frame(values = rbind(as.matrix(avg_width_dceitt1, ncol = 1), 
                                         as.matrix(avg_width_dcend1, ncol = 1), 
                                         as.matrix(avg_width_dceD1, ncol = 1), 
                                         as.matrix(avg_width_dced1[1, ], ncol = 1), 
                                         as.matrix(avg_width_dced1[2, ], ncol = 1), 
                                         as.matrix(avg_width_dced1[3, ], ncol = 1), 
                                         as.matrix(avg_width_dced1[4, ], ncol = 1), 
                                         as.matrix(avg_width_dceitt1_nocov, ncol = 1), 
                                         as.matrix(avg_width_dcend1_nocov, ncol = 1), 
                                         as.matrix(avg_width_dceD1_nocov, ncol = 1), 
                                         as.matrix(avg_width_dced1_nocov[1, ], ncol = 1), 
                                         as.matrix(avg_width_dced1_nocov[2, ], ncol = 1), 
                                         as.matrix(avg_width_dced1_nocov[3, ], ncol = 1), 
                                         as.matrix(avg_width_dced1_nocov[4, ], ncol = 1), 
                                         as.matrix(avg_width_dceitt2, ncol = 1), 
                                         as.matrix(avg_width_dcend2, ncol = 1), 
                                         as.matrix(avg_width_dceD2, ncol = 1), 
                                         as.matrix(avg_width_dced2[1, ], ncol = 1), 
                                         as.matrix(avg_width_dced2[2, ], ncol = 1), 
                                         as.matrix(avg_width_dced2[3, ], ncol = 1), 
                                         as.matrix(avg_width_dced2[4, ], ncol = 1), 
                                         as.matrix(avg_width_dceitt2_nocov, ncol = 1), 
                                         as.matrix(avg_width_dcend2_nocov, ncol = 1), 
                                         as.matrix(avg_width_dceD2_nocov, ncol = 1), 
                                         as.matrix(avg_width_dced2_nocov[1, ], ncol = 1), 
                                         as.matrix(avg_width_dced2_nocov[2, ], ncol = 1), 
                                         as.matrix(avg_width_dced2_nocov[3, ], ncol = 1), 
                                         as.matrix(avg_width_dced2_nocov[4, ], ncol = 1)), 
                          scenario = rep(c("1", "1 w/o covariates", "2", "2 w/o covariates"), 
                                         each = 20 * 7), 
                          effect = rep(c("ITT", "ND", "D", "d=1", "d=2", "d=3", "d=4"), 
                                       each = 20))
  
  avg_width$effect <- factor(avg_width$effect, 
                             levels = c("ITT", "ND", "D", "d=1", "d=2", "d=3", "d=4", "d=5"))
  
  jpeg("plots/FULL/Figure8.jpeg", 
       width = 7000, height = 5000, units = "px", res = 1200)
  print(ggplot(avg_width, aes(x = scenario, y = values, fill = effect)) + 
          labs(x = "Scenario", y = "HPD width", fill = NULL) + 
          geom_boxplot() + 
          scale_fill_manual(values = okabe_ito[1:length(unique(coverage$effect))]) +
          theme_minimal() + 
          coord_cartesian(ylim = c(0, 2)))
  dev.off()
  
  ####### BIAS: Figure 11 of the manuscript
  bias <- data.frame(values = rbind(as.matrix(dceitt_mean_bias1, ncol = 1), 
                                    as.matrix(dce.nd_mean_bias1, ncol = 1), 
                                    as.matrix(dce.D_mean_bias1, ncol = 1), 
                                    as.matrix(dce.d_mean_bias1[1, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias1[2, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias1[3, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias1[4, ], ncol = 1), 
                                    as.matrix(dceitt_mean_bias1_nocov, ncol = 1), 
                                    as.matrix(dce.nd_mean_bias1_nocov, ncol = 1), 
                                    as.matrix(dce.D_mean_bias1_nocov, ncol = 1), 
                                    as.matrix(dce.d_mean_bias1_nocov[1, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias1_nocov[2, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias1_nocov[3, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias1_nocov[4, ], ncol = 1), 
                                    as.matrix(dceitt_mean_bias1_nocov, ncol = 1), 
                                    as.matrix(dceitt_mean_bias2, ncol = 1), 
                                    as.matrix(dce.nd_mean_bias2, ncol = 1), 
                                    as.matrix(dce.D_mean_bias2, ncol = 1), 
                                    as.matrix(dce.d_mean_bias2[1, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias2[2, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias2[3, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias2[4, ], ncol = 1), 
                                    as.matrix(dce.nd_mean_bias2_nocov, ncol = 1), 
                                    as.matrix(dce.D_mean_bias2_nocov, ncol = 1), 
                                    as.matrix(dce.d_mean_bias2_nocov[1, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias2_nocov[2, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias2_nocov[3, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias2_nocov[4, ], ncol = 1)), 
                     scenario = rep(c("1", "1 w/o covariates", "2", "2 w/o covariates"), 
                                    each = 20 * 7), 
                     effect = rep(c("ITT", "ND", "D", "d=1", "d=2", "d=3", "d=4"), 
                                  each = 20))
  
  bias$effect <- factor(bias$effect, 
                        levels = c("ITT", "ND", "D", "d=1", "d=2", "d=3", "d=4", "d=5"))
  
  jpeg("plots/FULL/Figure11.jpeg", 
       width = 7000, height = 5000, units = "px", res = 1200)
  print(ggplot(bias, aes(x = scenario, y = values, fill = effect)) + 
          labs(x = "Scenario", y = "Bias", fill = NULL) + 
          geom_boxplot() + 
          scale_fill_manual(values = okabe_ito[1:length(unique(coverage$effect))]) +
          theme_minimal() + coord_cartesian(ylim = c(-0.2, 1.6)))
  dev.off()
  
  ##### RESTRICTED MEAN SURVIVAL TIME ############################
  
  ####### COVERAGE: Figure 9 of the manuscript
  coverage_rmste <- data.frame(values = rbind(as.matrix(coverage_rmsteitt1[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(coverage_rmstend1[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(coverage_rmsteD1[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(coverage_rmsteitt1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(coverage_rmstend1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(coverage_rmsteD1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(coverage_rmsteitt2[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(coverage_rmstend2[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(coverage_rmsteD2[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(coverage_rmsteitt2_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(coverage_rmstend2_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(coverage_rmsteD2_nocov[c(4, 8, 12, 18)], ncol = 1)), 
                               scenario = rep(c("1", "1 w/o covariates", "2", "2 w/o covariates"), 
                                              each = 4 * 3), 
                               effect = rep(c("ITT", "ND", "D"), each = 4), 
                               tau = rep(c("4", "8", "12", "18")))
  
  coverage_rmste$effect <- factor(coverage_rmste$effect, levels = c("ITT", "ND", "D"))
  coverage_rmste$tau <- factor(coverage_rmste$tau, levels = c("4", "8", "12", "18"))
  
    jpeg("plots/FULL/Figure9.jpeg",
         width = 7000, height = 5000, units = "px", res = 1200)
    print(ggplot(coverage_rmste, 
                 aes(x = scenario, y = values, shape = tau, color = effect, group = effect)) + 
            geom_point(size = 3, position = position_dodge(width = 0.6)) +
            scale_shape_manual(values = c(21, 22, 23, 24)) + 
            scale_color_manual(values = okabe_ito[c(1:3)]) + 
            labs(x = "Scenario", y = "Coverage", color = NULL, shape = expression(tau)) +
            geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") + 
            theme_minimal() + theme(legend.position = "right") + 
            coord_cartesian(ylim = c(0, 1)))
    dev.off()
    
    ####### WIDTH: Figure 10 of the manuscript
    rmste_avg_width <- data.frame(values = rbind(as.matrix(avg_width_rmsteitt1[c(4, 8, 12, 18)], ncol = 1), 
                                                 as.matrix(avg_width_rmstend1[c(4, 8, 12, 18)], ncol = 1), 
                                                 as.matrix(avg_width_rmsteD1[c(4, 8, 12, 18)], ncol = 1), 
                                                 as.matrix(avg_width_rmsteitt1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                                 as.matrix(avg_width_rmstend1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                                 as.matrix(avg_width_rmsteD1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                                 as.matrix(avg_width_rmsteitt2[c(4, 8, 12, 18)], ncol = 1), 
                                                 as.matrix(avg_width_rmstend2[c(4, 8, 12, 18)], ncol = 1), 
                                                 as.matrix(avg_width_rmsteD2[c(4, 8, 12, 18)], ncol = 1), 
                                                 as.matrix(avg_width_rmsteitt2_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                                 as.matrix(avg_width_rmstend2_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                                 as.matrix(avg_width_rmsteD2_nocov[c(4, 8, 12, 18)], ncol = 1)), 
                                  scenario = rep(c("1", "1 w/o covariates", "2", "2 w/o covariates"), 
                                                 each = 4 * 3), 
                                  effect = rep(c("ITT", "ND", "D"), 
                                               each = 4), 
                                  tau = rep(c("4", "8", "12", "18")))
    
    rmste_avg_width$effect <- factor(rmste_avg_width$effect, levels = c("ITT", "ND", "D"))
    rmste_avg_width$tau <- factor(rmste_avg_width$tau, levels = c("4", "8", "12", "18"))
    
      jpeg("plots/FULL/Figure10.jpeg",
           width = 7000, height = 5000, units = "px", res = 1200)
      print(ggplot(rmste_avg_width, 
                   aes(x = scenario, y = values, shape = tau, color = effect, group = effect)) + 
              geom_point(size = 3, position = position_dodge(width = 0.6)) +
              scale_shape_manual(values = c(21, 22, 23, 24)) + 
              scale_color_manual(values = okabe_ito) + 
              labs(x = "Scenario", y = "HPD width", color = NULL, shape = expression(tau)) +
              theme_minimal() + theme(legend.position = "right") + 
              coord_cartesian(ylim = c(0, 10)))
      dev.off()
      
      ####### BIAS: Figure 12 of the manuscript
      rmste_bias <- data.frame(values = rbind(as.matrix(rmsteitt_mean_bias1[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(rmste.nd_mean_bias1[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(rmste.D_mean_bias1[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(rmsteitt_mean_bias1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(rmste.nd_mean_bias1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(rmste.D_mean_bias1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(rmsteitt_mean_bias2[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(rmste.nd_mean_bias2[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(rmste.D_mean_bias2[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(rmsteitt_mean_bias2_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(rmste.nd_mean_bias2_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(rmste.D_mean_bias2_nocov[c(4, 8, 12, 18)], ncol = 1)), 
                               scenario = rep(c("1", "1 w/o covariates", "2", "2 w/o covariates"), 
                                              each = 4 * 3), effect = rep(c("ITT", "ND", "D"), 
                                                                          each = 4), 
                               tau = rep(c("4", "8", "12", "18")))
      
      rmste_bias$effect <- factor(rmste_bias$effect, levels = c("ITT", "ND", "D"))
      rmste_bias$tau <- factor(rmste_bias$tau, levels = c("4", "8", "12", "18"))
      
        jpeg("plots/FULL/Figure12.jpeg",
             width = 7000, height = 5000, units = "px", res = 1200)
        print(ggplot(rmste_bias, 
                     aes(x = scenario, y = values, shape = tau, color = effect, group = effect)) + 
                geom_point(size = 3, position = position_dodge(width = 0.6)) + 
                scale_shape_manual(values = c(21, 22, 23, 24)) + 
                scale_color_manual(values = okabe_ito) + 
                labs(x = "Scenario", y = "Bias", color = NULL, shape = expression(tau)) + 
                theme_minimal() + 
                theme(legend.position = "right"))
        dev.off()
}else{
  coverage <- data.frame(values = rbind(as.matrix(coverage_dceitt1_nocov, ncol = 1), 
                                        as.matrix(coverage_dcend1_nocov, ncol = 1), 
                                        as.matrix(coverage_dceD1_nocov, ncol = 1), 
                                        as.matrix(coverage_dced1_nocov[1, ], ncol = 1), 
                                        as.matrix(coverage_dced1_nocov[2, ], ncol = 1), 
                                        as.matrix(coverage_dced1_nocov[3, ], ncol = 1), 
                                        as.matrix(coverage_dced1_nocov[4, ], ncol = 1), 
                                        as.matrix(coverage_dceitt2_nocov, ncol = 1), 
                                        as.matrix(coverage_dcend2_nocov, ncol = 1), 
                                        as.matrix(coverage_dceD2_nocov, ncol = 1), 
                                        as.matrix(coverage_dced2_nocov[1, ], ncol = 1), 
                                        as.matrix(coverage_dced2_nocov[2, ], ncol = 1), 
                                        as.matrix(coverage_dced2_nocov[3, ], ncol = 1), 
                                        as.matrix(coverage_dced2_nocov[4, ], ncol = 1)), 
                         scenario = rep(c("1 w/o covariates", "2 w/o covariates"), 
                                        each = 20 * 7), 
                         effect = rep(c("ITT", "ND", "D", "d=1", "d=2", "d=3", "d=4"),
                                      each = 20))
  
  coverage$effect <- factor(coverage$effect, 
                            levels = c("ITT", "ND", "D", "d=1", "d=2", "d=3", "d=4", "d=5"))
  
  jpeg("plots/INTERMEDIATE/Figure7.jpeg",
       width = 7000, height = 5000, units = "px", res = 1200)
  print(ggplot(coverage, aes(x = scenario, y = values, fill = effect)) + 
    labs(x = "Scenario", y = "Coverage", fill = NULL) + geom_boxplot() + 
    geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") + 
    scale_fill_manual(values = okabe_ito[1:length(unique(coverage$effect))]) + 
    theme_minimal())
  dev.off()
  
  ####### WIDTH: Figure 8 of the manuscript
  avg_width <- data.frame(values = rbind(as.matrix(avg_width_dceitt1_nocov, ncol = 1), 
                                         as.matrix(avg_width_dcend1_nocov, ncol = 1), 
                                         as.matrix(avg_width_dceD1_nocov, ncol = 1), 
                                         as.matrix(avg_width_dced1_nocov[1, ], ncol = 1), 
                                         as.matrix(avg_width_dced1_nocov[2, ], ncol = 1), 
                                         as.matrix(avg_width_dced1_nocov[3, ], ncol = 1), 
                                         as.matrix(avg_width_dced1_nocov[4, ], ncol = 1), 
                                         as.matrix(avg_width_dceitt2_nocov, ncol = 1), 
                                         as.matrix(avg_width_dcend2_nocov, ncol = 1), 
                                         as.matrix(avg_width_dceD2_nocov, ncol = 1), 
                                         as.matrix(avg_width_dced2_nocov[1, ], ncol = 1), 
                                         as.matrix(avg_width_dced2_nocov[2, ], ncol = 1), 
                                         as.matrix(avg_width_dced2_nocov[3, ], ncol = 1), 
                                         as.matrix(avg_width_dced2_nocov[4, ], ncol = 1)), 
                          scenario = rep(c("w/o covariates", "2 w/o covariates"), 
                                         each = 20 * 7), 
                          effect = rep(c("ITT", "ND", "D", "d=1", "d=2", "d=3", "d=4"), 
                                       each = 20))
  
  avg_width$effect <- factor(avg_width$effect, 
                             levels = c("ITT", "ND", "D", "d=1", "d=2", "d=3", "d=4", "d=5"))
  
  jpeg("plots/INTERMEDIATE/Figure8.jpeg", 
       width = 7000, height = 5000, units = "px", res = 1200)
  print(ggplot(avg_width, aes(x = scenario, y = values, fill = effect)) + 
          labs(x = "Scenario", y = "HPD width", fill = NULL) + 
          geom_boxplot() + 
          scale_fill_manual(values = okabe_ito[1:length(unique(coverage$effect))]) +
          theme_minimal() + 
          coord_cartesian(ylim = c(0, 2)))
  dev.off()

  ####### BIAS: Figure 11 of the manuscript
  bias <- data.frame(values = rbind(as.matrix(dceitt_mean_bias1_nocov, ncol = 1), 
                                    as.matrix(dce.nd_mean_bias1_nocov, ncol = 1), 
                                    as.matrix(dce.D_mean_bias1_nocov, ncol = 1), 
                                    as.matrix(dce.d_mean_bias1_nocov[1, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias1_nocov[2, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias1_nocov[3, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias1_nocov[4, ], ncol = 1), 
                                    as.matrix(dceitt_mean_bias1_nocov, ncol = 1), 
                                    as.matrix(dce.nd_mean_bias2_nocov, ncol = 1), 
                                    as.matrix(dce.D_mean_bias2_nocov, ncol = 1), 
                                    as.matrix(dce.d_mean_bias2_nocov[1, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias2_nocov[2, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias2_nocov[3, ], ncol = 1), 
                                    as.matrix(dce.d_mean_bias2_nocov[4, ], ncol = 1)), 
                     scenario = rep(c("1 w/o covariates", "2 w/o covariates"), 
                                    each = 20 * 7), 
                     effect = rep(c("ITT", "ND", "D", "d=1", "d=2", "d=3", "d=4"), 
                                  each = 20))
  
  bias$effect <- factor(bias$effect, 
                        levels = c("ITT", "ND", "D", "d=1", "d=2", "d=3", "d=4", "d=5"))
  
  jpeg("plots/INTERMEDIATE/Figure11.jpeg",
       width = 7000, height = 5000, units = "px", res = 1200)
  print(ggplot(bias, aes(x = scenario, y = values, fill = effect)) + 
    labs(x = "Scenario", y = "Bias", fill = NULL) + 
    geom_boxplot() + 
    scale_fill_manual(values = okabe_ito[1:length(unique(coverage$effect))]) +
    theme_minimal() + coord_cartesian(ylim = c(-0.2, 1.6)))
  dev.off()
  
  ##### RESTRICTED MEAN SURVIVAL TIME ############################
  
  ####### COVERAGE: Figure 9 of the manuscript
  coverage_rmste <- data.frame(values = rbind(as.matrix(coverage_rmsteitt1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(coverage_rmstend1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(coverage_rmsteD1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(coverage_rmsteitt2_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(coverage_rmstend2_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                              as.matrix(coverage_rmsteD2_nocov[c(4, 8, 12, 18)], ncol = 1)), 
                               scenario = rep(c("1 w/o covariates", "2 w/o covariates"), 
                                              each = 4 * 3), 
                               effect = rep(c("ITT", "ND", "D"), each = 4), 
                               tau = rep(c("4", "8", "12", "18")))
  
  coverage_rmste$effect <- factor(coverage_rmste$effect, levels = c("ITT", "ND", "D"))
  coverage_rmste$tau <- factor(coverage_rmste$tau, levels = c("4", "8", "12", "18"))
  
  jpeg("plots/INTERMEDIATE/Figure9.jpeg",
       width = 7000, height = 5000, units = "px", res = 1200)
  print(ggplot(coverage_rmste, 
               aes(x = scenario, y = values, shape = tau, color = effect, group = effect)) + 
          geom_point(size = 3, position = position_dodge(width = 0.6)) +
          scale_shape_manual(values = c(21, 22, 23, 24)) + 
          scale_color_manual(values = okabe_ito[c(1:3)]) + 
          labs(x = "Scenario", y = "Coverage", color = NULL, shape = expression(tau)) +
          geom_hline(yintercept = 0.95, linetype = "dashed", color = "black") + 
          theme_minimal() + 
          theme(legend.position = "right") + 
          coord_cartesian(ylim = c(0, 1)))
  dev.off()
  
  ####### WIDTH: Figure 10 of the manuscript
  rmste_avg_width <- data.frame(values = rbind(as.matrix(avg_width_rmsteitt1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                               as.matrix(avg_width_rmstend1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                               as.matrix(avg_width_rmsteD1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                               as.matrix(avg_width_rmsteitt2_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                               as.matrix(avg_width_rmstend2_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                               as.matrix(avg_width_rmsteD2_nocov[c(4, 8, 12, 18)], ncol = 1)), 
                                scenario = rep(c("1 w/o covariates", "2 w/o covariates"), 
                                               each = 4 * 3), 
                                effect = rep(c("ITT", "ND", "D"), 
                                             each = 4), 
                                tau = rep(c("4", "8", "12", "18")))
  
  rmste_avg_width$effect <- factor(rmste_avg_width$effect, levels = c("ITT", "ND", "D"))
  rmste_avg_width$tau <- factor(rmste_avg_width$tau, levels = c("4", "8", "12", "18"))
  
  jpeg("plots/INTERMEDIATE/Figure10.jpeg",
       width = 7000, height = 5000, units = "px", res = 1200)
  print(ggplot(rmste_avg_width, 
               aes(x = scenario, y = values, shape = tau, color = effect, group = effect)) + 
          geom_point(size = 3, position = position_dodge(width = 0.6)) +
          scale_shape_manual(values = c(21, 22, 23, 24)) + 
          scale_color_manual(values = okabe_ito) + 
          labs(x = "Scenario", y = "HPD width", color = NULL, shape = expression(tau)) +
          theme_minimal() + theme(legend.position = "right") + 
          coord_cartesian(ylim = c(0, 10)))
  dev.off()
  
  ####### BIAS: Figure 12 of the manuscript
  rmste_bias <- data.frame(values = rbind(as.matrix(rmsteitt_mean_bias1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                          as.matrix(rmste.nd_mean_bias1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                          as.matrix(rmste.D_mean_bias1_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                          as.matrix(rmsteitt_mean_bias2_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                          as.matrix(rmste.nd_mean_bias2_nocov[c(4, 8, 12, 18)], ncol = 1), 
                                          as.matrix(rmste.D_mean_bias2_nocov[c(4, 8, 12, 18)], ncol = 1)), 
                           scenario = rep(c("1 w/o covariates", "2 w/o covariates"), 
                                          each = 4 * 3), effect = rep(c("ITT", "ND", "D"), 
                                                                      each = 4), 
                           tau = rep(c("4", "8", "12", "18")))
  
  rmste_bias$effect <- factor(rmste_bias$effect, levels = c("ITT", "ND", "D"))
  rmste_bias$tau <- factor(rmste_bias$tau, levels = c("4", "8", "12", "18"))
  
  jpeg("plots/INTERMEDIATE/Figure12.jpeg", 
       width = 7000, height = 5000, units = "px", res = 1200)
  print(ggplot(rmste_bias, 
         aes(x = scenario, y = values, shape = tau, color = effect, group = effect)) + 
    geom_point(size = 3, position = position_dodge(width = 0.6)) + 
    scale_shape_manual(values = c(21, 22, 23, 24)) + 
    scale_color_manual(values = okabe_ito) + 
    labs(x = "Scenario", y = "Bias", color = NULL, shape = expression(tau)) + 
    theme_minimal() + theme(legend.position = "right"))
  dev.off()
}
