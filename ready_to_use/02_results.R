################################################################################
####    Evaluating causal effects on time-to-event outcomes in an RCT in   #####
####                Oncology with treatment discontinuation                #####
####                                                                       #####
####             Authors: V. Ballerini, B. Bornkamp, A. Mattei,            #####
####                      F. Mealli, C. Wang, Y. Zhang                     #####
####                                                                       #####
#### Use this file to obtain a table of results and as a necessary step before #
#### running script 03_graphs.R. If you have more than 3 covariates, you have ##
#### to slightly modify the object at line 56 of this script.              #####
################################################################################

rm(list=setdiff(ls(), "file.name"))
wd <- dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(wd)

load(file.name)

library(HDInterval)
library(parallel)
library(ggplot2)
library(xtable)

THETA<-rbind(chain1$Theta[seq(burn, niter, by=thin),])
# if more chains:
# THETA<-rbind(chain1$Theta[seq(burn, niter, by=thin),],
#              chain2$Theta[seq(burn, niter, by=thin),],
#              chain3$Theta[seq(burn, niter, by=thin),])

# Parameters' mean
apply(THETA,2,mean)

# Trace plots of the parameters
for(i in 1:dim(THETA)[2]){
  plot(THETA[,i],type="l",xlab="iter",
       ylab=paste(colnames(THETA)[i]))
}

DID<-array(NA,dim=c(dim(chain1$DID[[1]])[1],
                    dim(chain1$DID[[1]])[2],
                    dim(THETA)[1]))
j<-NULL
for(i in seq({burn}, niter, by=thin)){
  j<-sum(j,1)
  DID[,,j]<-chain1$DID[[i]]
}

## Modify the following if you have more than 3 covariates, simply adding 
# eta.4 = THETA[i, "eta.4"]
# ...
# eta.d4 = THETA[i, "eta.d4"]
# ...
#
THETA_list<-list()
for(i in 1:dim(THETA)[1]){
  THETA_list[[i]]<-list(eta.0=THETA[i,"eta.0"],
                        eta.1=THETA[i,"eta.1"],
                        eta.2=THETA[i,"eta.2"],
                        eta.3=THETA[i,"eta.3"], ## add parameters here if you have n_cov > 3
                        alpha.d=THETA[i,"alpha.d"],
                        beta.d=THETA[i,"beta.d"],
                        eta.d1=THETA[i,"eta.d1"],
                        eta.d2=THETA[i,"eta.d2"],
                        eta.d3=THETA[i,"eta.d3"], ## add parameters here if you have n_cov > 3
                        alpha.y1nd=THETA[i,"alpha.y1nd"],
                        beta.y1nd=THETA[i,"beta.y1nd"],
                        alpha.y1d=THETA[i,"alpha.y1d"],
                        beta.y1d=THETA[i,"beta.y1d"],
                        alpha.y0nd=THETA[i,"alpha.y0nd"],
                        beta.y0nd=THETA[i,"beta.y0nd"],
                        alpha.y0d=THETA[i,"alpha.y0d"],
                        beta.y0d=THETA[i,"beta.y0d"],
                        eta.y1=THETA[i,"eta.y1"],
                        eta.y2=THETA[i,"eta.y2"],
                        eta.y3=THETA[i,"eta.y3"], ## add parameters here if you have n_cov > 3
                        lambda=THETA[i,"lambda"],
                        D=DID[,,i][,1],
                        ID=DID[,,i][,2])
}

prova<-lapply(THETA_list, ce_apply)
cov_distr<-lapply(THETA_list, apply_cov)

ace.d<-matrix(NA,dim(THETA)[1],10)
for(i in 1:dim(THETA)[1]){ace.d[i,]<-prova[[i]]$ace.d}

ace.fd<-rep(NA,dim(THETA)[1])
for(i in 1:dim(THETA)[1]){ace.fd[i]<-sum(prova[[i]]$ace.fd)}

ace.nd<-rep(NA,dim(THETA)[1])
for(i in 1:dim(THETA)[1]){ace.nd[i]<-prova[[i]]$ace.nd}

ace.overall<-rep(NA,dim(THETA)[1])
for(i in 1:dim(THETA)[1]){ace.overall[i]<-prova[[i]]$aceoverall}

ND_mat<-DID[,2,]
pc_ND<-mean(apply(ND_mat,2,mean))
hdi_ND<-hdi(apply(ND_mat,2,mean))

cfr1<-rbind(c(pc_ND,hdi_ND[1],hdi_ND[2]),
            c(mean(ace.overall),hdi(ace.overall,.95)[1], hdi(ace.overall,.95)[2]),
            c(mean(ace.nd),hdi(ace.nd,.95)[1],hdi(ace.nd,.95)[2]),
            c(mean(ace.fd),hdi(ace.fd,.95)[1],hdi(ace.fd,.95)[2]))

cfr2<-cbind(colMeans(ace.d)[1:5],apply(ace.d,2,hdi)[1,1:5],apply(ace.d,2,hdi)[2,1:5])
cfr<-rbind(cfr1,cfr2)

colnames(cfr)<-c("post mean", "lower HDI", "upper HDI")
rownames(cfr)<-c("%ND","ACE","ACE_ND","ACE_D","ACE(D(1)=1)","ACE(D(1)=2)",
                 "ACE(D(1)=3)","ACE(D(1)=4)","ACE(D(1)=5)")

xtable(cfr)

prova_dce<-lapply(THETA_list, dce_apply)

dce.nd<-matrix(NA,dim(THETA)[1],30)
for(i in 1:dim(THETA)[1]){dce.nd[i,]<-prova_dce[[i]]$dce.nd}

dce.d_1<-matrix(NA,dim(THETA)[1],30)
for(i in 1:dim(THETA)[1]){dce.d_1[i,]<-prova_dce[[i]]$dce.d[1,]}

dce.d_2<-matrix(NA,dim(THETA)[1],30)
for(i in 1:dim(THETA)[1]){dce.d_2[i,]<-prova_dce[[i]]$dce.d[2,]}

dce.d_3<-matrix(NA,dim(THETA)[1],30)
for(i in 1:dim(THETA)[1]){dce.d_3[i,]<-prova_dce[[i]]$dce.d[3,]}

dce.d_4<-matrix(NA,dim(THETA)[1],30)
for(i in 1:dim(THETA)[1]){dce.d_4[i,]<-prova_dce[[i]]$dce.d[4,]}

dce.d_5<-matrix(NA,dim(THETA)[1],30)
for(i in 1:dim(THETA)[1]){dce.d_5[i,]<-prova_dce[[i]]$dce.d[5,]}

