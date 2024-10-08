################################################################################
####    Evaluating causal effects on time-to-event outcomes in an RCT in   #####
####                Oncology with treatment discontinuation                #####
####                                                                       #####
####             Authors: V. Ballerini, B. Bornkamp, A. Mattei,            #####
####                      F. Mealli, C. Wang, Y. Zhang                     #####
####                                                                       #####
####   This code reproduces results of scenario 1 included in the paper    #####
################################################################################

rm(list=ls())
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

#### Load libraries and functions
source("Functions.R")
source("ce_apply.R")
source("Estimands.R")
source("MCMC.R")
source("CompleteLogPost.R")
source("DataAugmentation.R")

#install.packages("TruncatedDistributions", repos="http://R-Forge.R-project.org")
library(HDInterval)
library(TruncatedDistributions)
library(survival)
library(mvtnorm)
library(parallel)
library(ggplot2)

#### Load data
observed.data<- read.table("simdata1_4816.txt",header=T)

#### MCMC preparation

## Set hyperparameters
a0<- .5#1
b0<- .5#1

aa0<-a0
bb0<-b0

m0<-0
s20 <- 10^2#10
s20.du <- 5^2#5

a.start <- 1.5
scale.start <- 1
m.start<-0
m.start_y<--3
sd.start<- 1

theta.prior <- list(mu.0 = .6, sigma2.0 = s20,
                    mu.1=m0, sigma2.1=s20.du,
                    mu.2=m0, sigma2.2=s20.du,
                    mu.3=m0, sigma2.3=s20.du,
                    a.d=a0, b.d=b0, 
                    mu.d=m0, sigma2.d= s20,
                    mu.d1=m0, sigma2.d1=s20.du,
                    mu.d2=m0, sigma2.d2=s20.du,
                    mu.d5=m0, sigma2.d5=s20.du,
                    a.y1nd=a0, b.y1nd=b0, 
                    mu.y1nd=m0, sigma2.y1nd=s20,
                    a.y1d=a0,  b.y1d=b0, 
                    mu.y1d=m0, sigma2.y1d=s20,
                    a.y0nd=aa0, b.y0nd=bb0, 
                    mu.y0nd=m0, sigma2.y0nd=s20,
                    a.y0d=aa0, b.y0d=bb0, 
                    mu.y0d=m0, sigma2.y0d=s20,
                    mu.y1=m0, sigma2.y1=s20.du,
                    mu.y2=m0, sigma2.y2=s20.du,
                    mu.y5=m0, sigma2.y5=s20.du,
                    mu.lambda=m.start,  sigma2.lambda=s20)

#### Set starting values
par.start <- list(mu.0 = m.start, sigma.0 = sd.start,
                  mu.1 = m.start, sigma.1 = sd.start,
                  mu.2 = m.start, sigma.2 = sd.start,
                  mu.3 = m.start, sigma.3 = sd.start,
                  alpha.d=1.2,  beta.d  = -2,
                  mu.d1 = m.start, sigma.d1 = sd.start,
                  mu.d2 = m.start, sigma.d2 = sd.start,
                  mu.d5 = m.start, sigma.d5 = sd.start,
                  shape.y1nd = a.start, scale.y1nd = scale.start, 
                  mu.y1nd = m.start_y, sigma.y1nd = sd.start,
                  shape.y1d  = a.start, scale.y1d = scale.start, 
                  mu.y1d = m.start_y,  sigma.y1d = sd.start,
                  shape.y0nd = a.start, scale.y0nd = scale.start, 
                  mu.y0nd = m.start_y, sigma.y0nd = sd.start, 
                  shape.y0d  = a.start,  scale.y0d = scale.start, 
                  mu.y0d  = m.start_y,   sigma.y0d  = sd.start,
                  mu.y1 = m.start, sigma.y1 = sd.start,
                  mu.y2 = m.start, sigma.y2 = sd.start,
                  mu.y5 = m.start, sigma.y5 = sd.start,
                  mu.lambda = m.start,  sigma.lambda = sd.start)

#### Set parameters of the proposal distributions
proposal <- list(sd.0 = 0.3, sd.1 = 0.3, sd.2 = 0.4, sd.3 = 0.7,
                 scale.d = 0.02, sd.d = 0.3,
                 sd.d1 = 0.3, sd.d2 = 0.5, sd.d5 = 0.6,
                 scale.y1nd=0.01, sd.y1nd =0.2,
                 scale.y1d=0.02,  sd.y1d =0.4,
                 # scale.y0nd=4.7, sd.y0nd = 2.3,
                 scale.y0nd=0.015, sd.y0nd = 0.2,
                 scale.y0d=0.03,  sd.y0d = 0.3,
                 sd.y1 = 0.15, sd.y2 = 0.2, sd.y5 = 0.3,
                 sd.lambda = 0.2)

## Initialization
niter<- 50000
burn<- 30000
thin<-5

## set seed for the chain
set.seed(474747)

## Run the chain
chain1<- mcmc.tdiscontinue(niter, burn, thin, par.start, theta.prior, 
                           observed.data, thetatrue=thetatrue)

#### Save results
string1<-substr(date(),5,7)
string2<-substr(date(),9,10)
string3<-substr(date(),12,13)
string4<-substr(date(),15,16)

save.image(paste(wd,"/",string1,string2,"_",string3,"-",
                 string4,".RData",sep=""))
