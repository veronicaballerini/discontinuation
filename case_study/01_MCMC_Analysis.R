################################################################################ 
# Evaluating causal effects on time-to-event outcomes in an RCT in Oncology 
# with treatment discontinuation
#
# Authors: V. Ballerini, B., Bornkamp, F. Mealli, C. Wang, Y. Zhang, A. Mattei
#
# Replication Package for results in Sections 2 and 6
# In particular: This script produces Tables 1 and 2 (left and right), and Figure 
# 2 in Section 6, and implement the model for the Exp-Exp specification.
#
# Code author: Veronica Ballerini
# Last modified: October 6, 2025 
################################################################################

# Optional: clean the environment 
# rm(list = ls())
wd <- "/home/BiomJourn_RepPack/case_study"  # change here your wd
setwd(wd)

# Load libraries and functions
source("Functions.R")

library(HDInterval)
library(survival)
library(mvtnorm)
library(parallel)
library(ggplot2)
library(ggfortify)
library(xtable)
library(dplyr)
library(RColorBrewer)

#### Load your data
observed.data <- read.csv("synthetic_data.csv", sep = ",")

source("descriptive.R")

## Initialization
niter <- 50000
burn <- 30000
thin <- 10

# Exp-Exp model: Exponential distributions for the potential time to
# discontinuation and the potential PFS
source("MCMC_expexp.R")
source("CompleteLogPost_expexp.R")
source("DataAugmentation_expexp.R")
source("MCMC_initialization_expexp.R")

set.seed(474756)
chain_expexp <- mcmc.tdiscontinue_expexp(niter = niter, burn = burn, thin = thin,
    par.start = par.start_expexp, proposal = proposal_expexp, theta.prior = theta.prior_expexp,
    observed.data = observed.data)

#### Save results
save(chain_expexp, file = "final_casestudy.RData")
