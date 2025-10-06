################################################################################ 
# Evaluating causal effects on time-to-event outcomes in an RCT in oncology
# with treatment discontinuation 
#
# Authors: V. Ballerini, B. Bornkamp, F .Mealli, C. Wang, Y. Zhang, A. Mattei
#
# Replication Package for results in Section 5 
# Use this script to run the simulation study (full or intermediate).
#
# Code authors: Veronica Ballerini, Alessandra Mattei 
# Last modified: October 6, 2025 
################################################################################ 

# Optional: clean the environment 
# rm(list=ls())

# set the working directory: change here you wd
wd <- "/home/as218929/BiomJourn_RepPack/sim_study"
setwd(wd)

# recall functions
source("Functions.R")
source("estimands_apply.R")
source("estimands_apply_nocov.R")
source("Estimands.R")
source("applyfun_1.R")
source("applyfun_1nocov.R")
source("applyfun_2.R")
source("applyfun_2nocov.R")
source("MCMC.R")
source("MCMC_nocov.R")
source("MCMC_initialization.R")

# upload libraries
library(HDInterval)
library(survival)
library(parallel)

# Choose whether you want to reproduce the full simulation (full <- TRUE) study,
# with 150 datasets, or an "intermediate" one, with 2 datasets (full <- FALSE). 
full <- FALSE

# set the number of iterations, burnin and thinning
niter <- 50000
burn <- 30000
thin <- 10

if (full == TRUE) {
  # upload vectors of seeds for reproducibility
  seeds1 <- read.table("seeds1.txt")
  seeds1_list <- list()
  for (i in 1:150) {
    seeds1_list[[i]] <- seeds1[i, 1]
  }
  
  seeds1_nocov <- read.table("seeds1_nocov.txt", header = T)
  seeds1_nocov_list <- list()
  for (i in 1:150) {
    seeds1_nocov_list[[i]] <- seeds1_nocov[i, 1]
  }
  
  seeds2 <- read.table("seeds2.txt")
  seeds2_list <- list()
  for (i in 1:150) {
    seeds2_list[[i]] <- seeds2[i, 1]
  }
  
  seeds2_nocov <- read.table("seeds2_nocov.txt")
  seeds2_nocov_list <- list()
  for (i in 1:150) {
    seeds2_nocov_list[[i]] <- seeds2_nocov[i, 1]
  }
  
  # MCMC initialization
  source("MCMC_initialization.R")
  
  start <- proc.time()
  
  # Run the simulation (using 150 cores) - Scenario 1, with and without covariates
  sim_SCENARIO1 <- mclapply(seeds1_list, applyfun_1, mc.cores = 150)
  sim_SCENARIO1_nocov <- mclapply(seeds1_nocov_list, applyfun_1nocov, mc.cores = 150)
  
  #                                      - Scenario 2, with and without covariates
  sim_SCENARIO2 <- mclapply(seeds2_list, applyfun_2, mc.cores = 150)
  sim_SCENARIO2_nocov <- mclapply(seeds2_nocov_list, applyfun_2nocov, mc.cores = 150)
  
  end <- proc.time()
  
  howlong <- end - start
  
  # save simulation data
  rm(list = setdiff(ls(), c("sim_SCENARIO1", "sim_SCENARIO1_nocov", "sim_SCENARIO2", 
                          "sim_SCENARIO2_nocov", "howlong")))
  save.image("final_simdata.RData")
} else {
    # upload vectors of seeds for reproducibility
    seeds1_nocov <- read.table("seeds1_nocov.txt", header = T)
    seeds1_nocov_list <- list()
    for (i in 1:2) {
      seeds1_nocov_list[[i]] <- seeds1_nocov[i, 1]
    }
    
    seeds2_nocov <- read.table("seeds2_nocov.txt")
    seeds2_nocov_list <- list()
    for (i in 1:2) {
      seeds2_nocov_list[[i]] <- seeds2_nocov[i, 1]
    }
    
    start <- proc.time()
    
    # Run the simulation (using 1 core) - Scenario 1, with and without covariates
    sim_SCENARIO1_nocov <- lapply(seeds1_nocov_list, applyfun_1nocov)
    
    #                                   - Scenario 2, with and without covariates
    sim_SCENARIO2_nocov <- lapply(seeds2_nocov_list, applyfun_2nocov)
    
    end <- proc.time()
    
    howlong <- end - start
    
    # save simulation data
    rm(list = setdiff(ls(), c("sim_SCENARIO1_nocov", "sim_SCENARIO2_nocov", "howlong")))
    save.image("final_simdata_intermediate.RData")
}
