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

source("Functions.R")
source("ce_apply.R")
source("Estimands.R")

library(HDInterval)
library(TruncatedDistributions)
library(parallel)
library(ggplot2)
library(xtable)

#### Load data
# change the name of the .RData file at the following line
load(paste(wd,"/Oct 8_18-15.RData",sep="")) ## change with the name of your file
realdata<-read.table("realdata1_4816.txt",header=T)
thetatrue<-myload("thetatrue1_4816.RData")

THETA<-rbind(chain1$Theta[seq(burn, niter, by=thin),])

apply(THETA,2,mean)

int95<-t(apply(THETA,2,hdi))
int99<-t(apply(THETA,2,hdi99))
truev<-unlist(thetatrue[-c(22)])
check<-cbind(truev>=int95[,1]&truev<=int95[,2],
             truev>=int99[,1]&truev<=int99[,2],int95,truev)
colnames(check)<-c("in HPD95%", "in HPD99%","lower 95","upper 95","true value")
check

for(i in 1:21){
  plot(THETA[,i],type="l",xlab="iter",
       ylab=paste(colnames(THETA)[i]))
  abline(h=truev[i],col="red")
}

dim(THETA)

DID<-array(NA,dim=c(dim(chain1$DID[[1]])[1],
                    dim(chain1$DID[[1]])[2],dim(THETA)[1]))
j<-NULL
for(i in seq({burn}, niter, by=thin)){
  j<-sum(j,1)
  DID[,,j]<-chain1$DID[[i]]
}

THETA_list<-list()
for(i in 1:dim(THETA)[1]){
  THETA_list[[i]]<-list(eta.0=THETA[i,1],
                        eta.1=THETA[i,2],
                        eta.2=THETA[i,3],
                        eta.3=THETA[i,4],
                        alpha.d=THETA[i,5],
                        beta.d=THETA[i,6],
                        eta.d1=THETA[i,7],
                        eta.d2=THETA[i,8],
                        eta.d3=THETA[i,9],
                        alpha.y1nd=THETA[i,10],
                        beta.y1nd=THETA[i,11],
                        alpha.y1d=THETA[i,12],
                        beta.y1d=THETA[i,13],
                        alpha.y0nd=THETA[i,14],
                        beta.y0nd=THETA[i,15],
                        alpha.y0d=THETA[i,16],
                        beta.y0d=THETA[i,17],
                        eta.y1=THETA[i,18],
                        eta.y2=THETA[i,19],
                        eta.y3=THETA[i,20],
                        lambda=THETA[i,21],
                        D=DID[,,i][,1],
                        ID=DID[,,i][,2])
}

prova<-lapply(THETA_list, ce_apply)
rce<-ce(thetatrue,seq(1,10),realdata,scenario=1)

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
prop_real<-sum(realdata$ID)/335

cfr1<-rbind(c(prop_real,pc_ND,hdi_ND[1],hdi_ND[2]),
            c(rce$aceoverall,mean(ace.overall),hdi(ace.overall,.95)[1],
              hdi(ace.overall,.95)[2]),
            c(rce$ace.nd,mean(ace.nd),hdi(ace.nd,.95)[1],hdi(ace.nd,.95)[2]),
            c(sum(rce$ace.fd),
              mean(ace.fd),hdi(ace.fd,.95)[1],hdi(ace.fd,.95)[2]) )

cfr2<-cbind(rce$ace.d[1:5],colMeans(ace.d)[1:5],apply(ace.d,2,hdi)[1,1:5],
            apply(ace.d,2,hdi)[2,1:5])
cfr<-rbind(cfr1,cfr2)
colnames(cfr)<-c("true value", "post mean", "lower HDI", "upper HDI")
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

