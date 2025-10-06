################################################################################ 
# Evaluating causal effects on time-to-event outcomes in an RCT in oncology
# with treatment discontinuation 
# Authors: V. Ballerini, B. Bornkamp, A. Mattei, F .Mealli, C. Wang, Y. Zhang 
# Replication Package for results in Section 5 
# In this script: Figures 3, 4, 5, 6
# Code authors: Veronica Ballerini, Alessandra Mattei 
# Last modified: September 18, 2025 
################################################################################ 

# load useful libraries
library(ggfortify)
library(survival)

############################## SCENARIO 1 ######################################

# load one dataset from SCENARIO 1
realdata1 <- read.table("sim/realdata1_4798.txt", header = TRUE)
observed.data1 <- read.table("sim/simdata1_4798.txt", header = TRUE)

n <- nrow(realdata1)
y1 <- c(realdata1$Y0, realdata1$Y1)
z <- c(rep(0, n), rep(1, n))
id1 <- c(realdata1$ID, realdata1$ID)

# Figure 3(a) of the manuscript:
jpeg(paste(wd, "/plots/Figure3a.jpeg", sep = ""), width = 5000, height = 5000, 
     units = "px", res = 1200)
autoplot(survfit(Surv(y1[id1 == 1]) ~ z[id1 == 1]), xlim = c(0, 33)) + 
  ggtitle("") + theme_light() + guides(fill = FALSE) + 
  labs(x = "Months", y = "Survival Probability") +
    theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) + 
  labs(colour = "Z")
dev.off()

# Figure 3(b) of the manuscript:
jpeg(paste(wd, "/plots/Figure3b.jpeg", sep = ""), width = 5000, height = 5000, 
     units = "px", res = 1200)
autoplot(survfit(Surv(y1[id1 == 0]) ~ z[id1 == 0]), xlim = c(0, 33)) + 
  ggtitle("") + theme_light() + guides(fill = FALSE) + 
  labs(x = "Months", y = "Survival Probability") +
    theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) + 
  labs(colour = "Z")
dev.off()

status <- ifelse(observed.data1$f.up == observed.data1$Y, 0, 1)

# Figure 4 of the manuscript:
jpeg(paste(wd, "/plots/Figure4.jpeg", sep = ""), width = 5000, height = 5000, 
     units = "px", res = 1200)
autoplot(survfit(Surv(observed.data1$Y, status) ~ observed.data1$Z), xlim = c(0, 33)) + 
  theme_light() + guides(fill = FALSE) + labs(x = "Months", y = "Survival Probability") +
    labs(colour = "Z")
dev.off()

############################## SCENARIO 2 ######################################

# load one dataset from SCENARIO 2
realdata2 <- read.table("sim/realdata2_4800.txt", header = TRUE)
observed.data2 <- read.table("sim/simdata2_4800.txt", header = TRUE)

y2 <- c(realdata2$Y0, realdata2$Y1)
id2 <- c(realdata2$ID, realdata2$ID)

# Figure 5(a) of the manuscript:
jpeg(paste(wd, "/plots/Figure5a.jpeg", sep = ""), width = 5000, height = 5000, 
     units = "px", res = 1200)
autoplot(survfit(Surv(y2[id2 == 1]) ~ z[id2 == 1]), xlim = c(0, 33)) + 
  ggtitle("") + theme_light() + guides(fill = FALSE) + 
  labs(x = "Months", y = "Survival Probability") +
    theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) + 
  labs(colour = "Z")
dev.off()

# Figure 5(b) of the manuscript:
jpeg(paste(wd, "/plots/Figure5b.jpeg", sep = ""), width = 5000, height = 5000, 
     units = "px", res = 1200)
autoplot(survfit(Surv(y2[id2 == 0]) ~ z[id2 == 0]), xlim = c(0, 30)) + 
  ggtitle("") + theme_light() + guides(fill = FALSE) + 
  labs(x = "Months", y = "Survival Probability") +
    theme(legend.title = element_text(size = 15), legend.text = element_text(size = 15)) + 
  labs(colour = "Z")
dev.off()

status <- ifelse(observed.data2$f.up == observed.data2$Y, 0, 1)

# Figure 6 of the manuscript:
jpeg(paste(wd, "/plots/Figure6.jpeg", sep = ""), width = 5000, height = 5000, 
     units = "px", res = 1200)
autoplot(survfit(Surv(observed.data2$Y, status) ~ observed.data2$Z), xlim = c(0, 30)) + 
  theme_light() + guides(fill = FALSE) + labs(x = "Months", y = "Survival Probability") +
    labs(colour = "Z")
dev.off()
