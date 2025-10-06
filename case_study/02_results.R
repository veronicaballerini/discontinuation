################################################################################ 
# Evaluating causal effects on time-to-event outcomes in an RCT in Oncology 
# with treatment discontinuation
#
# Authors: V. Ballerini, B., Bornkamp, F. Mealli, C. Wang, Y. Zhang, A. Mattei
#
# Replication Package for results in Section 6 
# In particular: a) from source "results_processing," delta interval mentioned 
# in Section 6, and PPPV for KMdm described in text in Section 6.1 
# [.rds outputs in wd]; b) Figures 13, 14, 15.
#
# Code authors: Veronica Ballerini
# Last modified: October 6, 2025 
################################################################################ 

wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

library(dplyr)
library(RColorBrewer)

source("results_processing.R") # for processing results. It takes some minutes.

okabe_ito <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00",
               "#CC79A7", "#000000")

# Figure 13
dce <- data.frame(DCE_d = rbind(as.matrix(colMeans(dce.itt)), 
                                as.matrix(colMeans(dce.nd)), 
                                as.matrix(colMeans(dce.D))), 
                  CI_lower = rbind(as.matrix(apply(dce.itt, 2, hdi)[1, ]), 
                                   as.matrix(apply(dce.nd, 2, hdi)[1, ]), 
                                   as.matrix(apply(dce.D, 2, hdi)[1, ])), 
                  CI_upper = rbind(as.matrix(apply(dce.itt, 2, hdi)[2, ]), 
                                   as.matrix(apply(dce.nd, 2, hdi)[2, ]), 
                                   as.matrix(apply(dce.D, 2, hdi)[2, ])), 
                  Discontinuation = rep(c("ITT", "ND", "D"), each = 20), 
                  Time = seq(1, 20))

dce$Discontinuation <- factor(dce$Discontinuation, levels = c("ITT", "ND", "D"))

jpeg("figures/Figure13.jpeg", width = 7000, height = 5000, units = "px", res = 1200)
ggplot(dce, aes(x = Time, y = DCE_d, color = Discontinuation, group = Discontinuation)) + 
  geom_line(size = 1) +
    geom_point(size = 2) + 
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.3, alpha = 0.7) +
    facet_wrap(~Discontinuation, ncol = 1) + 
  scale_color_manual(values = okabe_ito[1:3]) + 
  labs(x = "y", y = "DCE", color = NULL) + 
  theme_minimal() + 
  coord_cartesian(xlim = c(1, 20), ylim = c(-0.2, 0.7)) + 
  theme(legend.position = "none")
dev.off()

# Figure 14
dce_d <- data.frame(DCE_d = rbind(as.matrix(colMeans(dce.d_1)), 
                                  as.matrix(colMeans(dce.d_2)), 
                                  as.matrix(colMeans(dce.d_3)), 
                                  as.matrix(colMeans(dce.d_4))), 
                    CI_lower = rbind(as.matrix(apply(dce.d_1, 2, hdi)[1, ]), 
                                     as.matrix(apply(dce.d_2, 2, hdi)[1, ]), 
                                     as.matrix(apply(dce.d_3, 2, hdi)[1, ]), 
                                     as.matrix(apply(dce.d_4, 2, hdi)[1, ])), 
                    CI_upper = rbind(as.matrix(apply(dce.d_1, 2, hdi)[2, ]), 
                                     as.matrix(apply(dce.d_2, 2, hdi)[2, ]), 
                                     as.matrix(apply(dce.d_3, 2, hdi)[2, ]), 
                                     as.matrix(apply(dce.d_4, 2, hdi)[2, ])), 
                    Discontinuation = rep(c("d=1", "d=2", "d=3", "d=4"), each = 20), 
                    Time = seq(1, 20))

jpeg("figures/Figure14.jpeg", width = 7000, height = 5000, units = "px", res = 1200)
ggplot(dce_d, aes(x = Time, y = DCE_d, color = Discontinuation, group = Discontinuation)) + 
  geom_line(size = 1) +
    geom_point(size = 2) + 
  geom_errorbar(aes(ymin = CI_lower, ymax = CI_upper), width = 0.3, alpha = 0.7) +
    facet_wrap(~Discontinuation, ncol = 1) + 
  scale_fill_manual(values = okabe_ito[4:7]) + 
  scale_color_manual(values = okabe_ito[4:7]) + 
  labs(x = "y", y = expression(DCE[D] * "(" * y * "|" * d * ")"), fill = NULL, color = NULL) +
    coord_cartesian(xlim = c(1, 20), ylim = c(-0.2, 0.7)) + 
  theme_minimal()
dev.off()

# Figure 15
d <- rep(c("ITT", "ND", "D", "d=1", "d=2", "d=3", "d=4"), each = 2001 * 4)
time <- rep(c("4", "8", "12", "18"), each = 2001)

rmste <- rbind(as.matrix(rmste.itt[, 4], ncol = 1), 
               as.matrix(rmste.itt[, 8], ncol = 1), 
               as.matrix(rmste.itt[, 12], ncol = 1), 
               as.matrix(rmste.itt[, 18], ncol = 1), 
               as.matrix(rmste.nd[, 4], ncol = 1), 
               as.matrix(rmste.nd[, 8], ncol = 1), 
               as.matrix(rmste.nd[, 12], ncol = 1), 
               as.matrix(rmste.nd[, 18], ncol = 1), 
               as.matrix(rmste.D[, 4], ncol = 1), 
               as.matrix(rmste.D[, 8], ncol = 1), 
               as.matrix(rmste.D[, 12], ncol = 1), 
               as.matrix(rmste.D[, 18], ncol = 1), 
               as.matrix(rmste.d_1[, 4], ncol = 1), 
               as.matrix(rmste.d_1[, 8], ncol = 1), 
               as.matrix(rmste.d_1[, 12], ncol = 1), 
               as.matrix(rmste.d_1[, 18], ncol = 1), 
               as.matrix(rmste.d_2[, 4], ncol = 1), 
               as.matrix(rmste.d_2[, 8], ncol = 1), 
               as.matrix(rmste.d_2[, 12], ncol = 1), 
               as.matrix(rmste.d_2[, 18], ncol = 1), 
               as.matrix(rmste.d_3[, 4], ncol = 1), 
               as.matrix(rmste.d_3[, 8], ncol = 1), 
               as.matrix(rmste.d_3[, 12], ncol = 1), 
               as.matrix(rmste.d_3[, 18], ncol = 1), 
               as.matrix(rmste.d_4[, 4], ncol = 1), 
               as.matrix(rmste.d_4[, 8], ncol = 1), 
               as.matrix(rmste.d_4[, 12], ncol = 1), 
               as.matrix(rmste.d_4[, 18], ncol = 1))
df <- data.frame(time, d, rmste)
df$d <- factor(df$d, levels = c("ITT", "ND", "D", "d=1", "d=2", "d=3", "d=4"))
df$time <- factor(df$time, levels = c("4", "8", "12", "18"))

summary_df <- df %>%
    group_by(d, time) %>%
    summarize(mean = mean(rmste), 
              lower = hdi(rmste, credMass = 0.95)[1], 
              upper = hdi(rmste, credMass = 0.95)[2])

jpeg("figures/Figure15.jpeg", width = 7000, height = 5000, units = "px", res = 1200)
ggplot(summary_df, aes(x = time, y = mean, color = d, group = d)) + geom_point(position = position_dodge(width = 0.5),
    size = 1.8) + geom_errorbar(aes(ymin = lower, ymax = upper), width = 0.2, position = position_dodge(width = 0.5)) +
    labs(x = expression(tau), y = expression(RMSTE[d](tau)), color = NULL) + scale_color_manual(values = okabe_ito[1:length(unique(summary_df$d))]) +
    theme_minimal()
dev.off()
