################################################################################
####    Evaluating causal effects on time-to-event outcomes in an RCT in   #####
####                Oncology with treatment discontinuation                #####
####                                                                       #####
####             Authors: V. Ballerini, B. Bornkamp, A. Mattei,            #####
####                      F. Mealli, C. Wang, Y. Zhang                     #####
################################################################################

rm(list=ls())
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

#### Load your data. Read the README file for information about the required structure.
observed.data<- read.table("example_data.txt",header=T)
x <- cbind(1,as.matrix(observed.data[,grepl("x",names(observed.data))]))
n_cov <- dim(x)[2] - 1

#### Load libraries and functions
source("Functions.R")
source("ce_apply.R")
source("MCMC.R")
source("CompleteLogPost.R")
source("DataAugmentation.R")
source("MCMC_initialization.R")
source("apply_cov.R")

library(HDInterval)
library(survival)
library(mvtnorm)
library(parallel)
library(ggplot2)
  
## Initialization
niter<- 50000
burn<- 30000
thin<-2

## set seed for chain 1
set.seed(474747)

## Run chain 1
chain1<- mcmc.tdiscontinue(niter, burn, thin, par.start, theta.prior, 
                           observed.data, n_cov)

## Check the acceptance rates: they must be between 0.25 and 0.5. If too low or
## too high, modify values of the proposal, e.g.:
# proposal$sd.0 <- 0.2
colMeans(chain1$jump[burn:niter,])

# ## Chain 2
# set.seed(474748)
# chain2<- mcmc.tdiscontinue(niter, burn, thin, par.start, theta.prior, 
#                            observed.data, n_cov)
# colMeans(chain1$jump[burn:niter,])
# 
# ## Chain 3
# set.seed(474749)
# chain3<- mcmc.tdiscontinue(niter, burn, thin, par.start, theta.prior, 
#                            observed.data, n_cov)
# colMeans(chain3$jump[burn:niter,])

#### Save results
string1<-substr(date(),5,7)
string2<-substr(date(),9,10)
string3<-substr(date(),12,13)
string4<-substr(date(),15,16)

file.name<-paste(wd,"/",string1,string2,"_",string3,"-",
                 string4,".RData",sep="")

save.image(file.name)
