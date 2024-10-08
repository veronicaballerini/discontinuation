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

library(ggfortify)

##############################################################################
############################## DATA SIMULATION ###############################
##############################################################################

#Pre-treatment variables:
#X1 = Continuous variable; the bigger the number, the higher the risk (Age)
#X2 = Presence of specific metastatic status before enrolling the study
#X5 = Disease burden indicator.
#Treatment variable:
#Z=1 investigational trt + Standard of Care (SOC)
#Z=0 SOC 
#Discontinuation indicator:
#ID = 1 if the discontinuation cannot be observed (due to death or progression)
#ID = 0 if the discontinuation is observed, even after the follow-up period
#Discontinuation variable:
#D1 = Continuous variable; time until the discontinuation under active treatment
#Outcome variable:
#Y = survival

seed<-4816
set.seed(seed) # Set the seed for reproducibility

n<-335 #total sample size
nT<-181 #number of treated individuals
nC<-n-nT #number of controls

x1<-rnorm(n,63.27,10.50) #age
x1st<-(x1-mean(x1))/sd(x1) #age - standardised
x2<-rbinom(n,1,.44) #1 if metastatic status prior to the enrolment, 0 otherwise
x5<-rbinom(n,1,1-.76) #1 if the desease burden is >=3, 0 if <3

Z<-matrix(c(rep(1,nT),rep(0,nC)),n,1) #treatment

############################### Discontinuation ################################
#Indicator I{D_i(1)=D} = 1 if the potential discontinuation under treatment does
#                        not apply [Never Discontinued]
#                      = 0 otherwise
#I{D_i(1)=D} ~ Bernoulli(pix_i)

eta<-list(int=.6,x1st=-.5,x2=.45,x5=.55) #intercept:  #previously x1st=.3
#we assumed the proportion of ND to be approximately equal to .65: lower bound 
#is 90/181 
#eta$x1st: we assumed the older the patient the lower the probability to be an 
#ND patient, since discontinuation is competing with progression;
#eta$x2,x5: we assumed the more severe the condition, the higher the probability 
#of being an ND patient (as above)

pix<-exp(eta$int+
           eta$x1*x1st+
           eta$x2*x2+
           eta$x5*x5)/(1+exp(eta$int+
                               eta$x1st*x1st+
                               eta$x2*x2+
                               eta$x5*x5))
ID <- rbinom(n,1,pix)

#Discontinuation variable:
#D1 ~ Weibull(shape=alpha.d, scale=beta.d+eta.d1*x1st+eta.d2*x2+eta.d5*x5)

#parameters of the Weibulls:
shape.d <- alpha.d <- 1.2

mean.d1 <- 4 #time to adverse event, average
par.scale.d <- list(beta.d = -alpha.d*log(mean.d1/gamma(1+1/alpha.d)),
                    eta.d1 = .45, eta.d2 = .25, eta.d5 = .35)
#eta.d1,d2,d5 = we assumed that individuals with higher risk
#discontinue sooner, given that they will discontinue
scale.d <- exp(-(par.scale.d$beta.d+
                   par.scale.d$eta.d1*x1st+
                   par.scale.d$eta.d2*x2+
                   par.scale.d$eta.d5*x5)/alpha.d)

D1 <- rep(NA,n) #potential discontinuation under treatment
D1[ID==0] <- rweibull(sum(1-ID),shape=shape.d,
                      scale=scale.d[ID==0]) #potential discontinuation under
#treatment for discontinuers
tc<-sample(0:1,n,replace = TRUE)
entry<-tc*rbinom(n,23,0.9)+(1-tc)*sample(0:23,n,replace=TRUE)
f.up<-33-entry #follow up period, from entry to study end (33 months)

#Outcome variable:
#Y1|D1=NA ~ Weibull(shape = alpha.y1nd,
#                   scale = beta.y1nd + eta.y1*x1st + eta.y2*x2 + eta.y5*x5)
#Y1|D1=d ~ Weibull(shape = alpha.y1d,
#                  scale = beta.y1d + lambda*log(D1) + eta.y1*x1st + eta.y2*x2 +
#                   eta.y5*x5)
#Y0|D1=NA ~ Weibull(shape = alpha.y0nd,
#                   scale = beta.y0nd + eta.y1*x1st + eta.y2*x2 + eta.y5*x5)
#Y0|D1=d ~ Weibull(shape=alpha.y0d,
#                  scale = beta.y0d + lambda*log(D1) + eta.y1*x1st + eta.y2*x2 +
#                   eta.y5*x5)

################### parameters under H0 of positive effect #####################
#parameters of the Weibulls:
shape.y1nd <- alpha.y1nd <- 1.7
shape.y0nd <- alpha.y0nd <- 1.6
shape.y1d <- alpha.y1d <- 1.6
shape.y0d <- alpha.y0d <- 1.5

mean.y1nd <- 11 #time to primary event, average
mean.y0nd <- 5 #time to primary event, average

mean.y1d <- 7 #time to primary event, average
mean.y0d <- 4  #time to primary event, average

beta.y1nd <- -alpha.y1nd*log(mean.y1nd/gamma(1+1/alpha.y1nd))
beta.y0nd <- -alpha.y0nd*log(mean.y0nd/gamma(1+1/alpha.y0nd))
beta.y1d <- -alpha.y1d*log(mean.y1d/gamma(1+1/alpha.y1d))
beta.y0d <- -alpha.y0d*log(mean.y0d/gamma(1+1/alpha.y0d))

#eta.y1,y2,y5 = we assumed that individuals with higher risk
#discontinue sooner, given that they will discontinue
par.scale.y1nd<-list(beta.y=beta.y1nd,
                     eta.y1=.25,
                     eta.y2=.7,
                     eta.y5=.45)
par.scale.y0nd<-list(beta.y=beta.y0nd,
                     eta.y1=.25,
                     eta.y2=.7,
                     eta.y5=.45)
par.scale.y1d<-list(beta.y=beta.y1d,
                    eta.y1=.25,
                    eta.y2=.7,
                    eta.y5=.45)
par.scale.y0d<-list(beta.y=beta.y0d,
                    eta.y1=.25,
                    eta.y2=.7,
                    eta.y5=.45)

theta<-list(alpha.y1nd = alpha.y1nd, beta.y1nd = beta.y1nd,
            eta.y1 = par.scale.y1nd$eta.y1,
            eta.y2 = par.scale.y1nd$eta.y2,
            eta.y5 = par.scale.y1nd$eta.y5,
            alpha.y1d = alpha.y1d, beta.y1d = beta.y1d,
            alpha.y0nd = alpha.y0nd, beta.y0nd = beta.y0nd,
            alpha.y0d = alpha.y0d, beta.y0d = beta.y0d,
            lambda = 0.3)

scale.y1nd <- exp(-(theta$beta.y1nd+theta$eta.y1*x1st+theta$eta.y2*x2+
                      theta$eta.y5*x5)/theta$alpha.y1nd)
scale.y1d <- exp(-(theta$beta.y1d+theta$lambda*log(D1)+theta$eta.y1*x1st+
                     theta$eta.y2*x2+theta$eta.y5*x5)/theta$alpha.y1d)
scale.y0nd <- exp(-(theta$beta.y0nd+theta$eta.y1*x1st+theta$eta.y2*x2+
                      theta$eta.y5*x5)/theta$alpha.y0nd)
scale.y0d <- exp(-(theta$beta.y0d+theta$lambda*log(D1)+theta$eta.y1*x1st+
                     theta$eta.y2*x2+theta$eta.y5*x5)/theta$alpha.y0d)

# Hazard function for NDs and Ds

Y1<-NULL #potential outcome under treatment
Y0<-NULL #potential outcome under control

#for ND:
Y1[ID==1]<-rweibull(sum(ID),shape=theta$alpha.y1nd,
                    scale=scale.y1nd[ID==1])
Y0[ID==1]<-rweibull(sum(ID),shape=theta$alpha.y0nd,
                    scale=scale.y0nd[ID==1])

Y1[ID==0]<-rtweibull(sum(1-ID),
                     shape=theta$alpha.y1d,
                     scale=scale.y1d[ID==0],
                     a=D1[ID==0])

if(sum(is.infinite(Y1)==TRUE)>0){
  Y1[which(is.infinite(Y1)==TRUE)]<-D1[which(is.infinite(Y1)==TRUE)]
  D1[which(is.infinite(Y1)==TRUE)]<-NA
  ID[which(is.infinite(Y1)==TRUE)]<-1
}

Y0[ID==0]<-rweibull(sum(1-ID),
                    shape=theta$alpha.y0d,
                    scale=scale.y0d[ID==0])

if(sum(is.infinite(Y0)==TRUE)>0){
  Y0[which(is.infinite(Y0)==TRUE)]<-D1[which(is.infinite(Y0)==TRUE)]
  D1[which(is.infinite(Y0)==TRUE)]<-NA
  ID[which(is.infinite(Y0)==TRUE)]<-1
}

y<-c(Y0,Y1)
d<-D1
z<-c(rep(0,n),rep(1,n))
id<-c(ID,ID)

jpeg("KM_ND_SC1.jpeg",width = 5000, height = 5000, units="px", res=1200)
autoplot(survfit(Surv(y[id==1])~z[id==1]),xlim=c(0,30))+
  ggtitle("")+
  theme_light()+
  guides(fill=FALSE) +
  labs(x = "Months", y = "Survival Probability") +
  theme(legend.title = element_text(size=15), legend.text = element_text(size=15))+
  labs(colour = "Z")
survfit(Surv(y[id==1])~z[id==1])
dev.off()

jpeg("KM_D_SC1.jpeg",width = 5000, height = 5000, units="px", res=1200)
autoplot(survfit(Surv(y[id==0])~z[id==0]),xlim=c(0,30))+
  ggtitle("") +
  theme_light()+
  guides(fill=FALSE) +
  labs(x = "Months", y = "Survival Probability") +
  theme(legend.title = element_text(size=15), legend.text = element_text(size=15))+
  labs(colour = "Z")
dev.off()

######## observed data ########
Dobs<-RD<-Dobs<-RY<-Yobs<-rep(0, n)

Dobs[Z==1 & ID==0]<- D1[Z==1 & ID==0]*(D1[Z==1 & ID==0]<f.up[Z==1 & ID==0]) +
  f.up[Z==1 & ID==0]*(D1[Z==1 & ID==0]>=f.up[Z==1 & ID==0])
Dobs[Z==1 & ID==1]<- f.up[Z==1 & ID==1]
Yobs[Z==0] <- Y0[Z==0]*(Y0[Z==0]<f.up[Z==0]) + f.up[Z==0]*(Y0[Z==0]>=f.up[Z==0])
Yobs[Z==1] <- Y1[Z==1]*(Y1[Z==1]<f.up[Z==1]) + f.up[Z==1]*(Y1[Z==1]>=f.up[Z==1])

RD[Z==1 & ID==0] <- 1*(D1[Z==1 & ID==0]<f.up[Z==1 & ID==0])
# indicator variable: 1 if the
# discontinuation is observed
RY[Z==0] <- 1*(Y0[Z==0]<=f.up[Z==0])        # indicator variables: 1 if the
RY[Z==1] <- 1*(Y1[Z==1]<=f.up[Z==1])        # death is observed

ID.obs <- rep(0,n) #0= Discontinued/Not available, 1=Never discontinued
ID.obs[Z==1 & RD==0 & RY==1] <- 1
ID.obs[Z==1 & RD==1] <- 0

RD.obs<-RD
Y<-Yobs
cens<-rep(33,n)

data<-data.frame(entry,Z,f.up,cens,RD.obs,Dobs,RY,Y,ID.obs,x1,x1st,x2,x5)

status<-ifelse(f.up==Yobs,0,1)

jpeg("KM_obs_SC1.jpeg",width = 5000, height = 5000, units="px", res=1200)
autoplot(survfit(Surv(Yobs,status)~Z),xlim=c(0,30))+
  theme_light()+
  guides(fill=FALSE) +
  labs(x = "Months", y = "Survival Probability") +
  labs(colour = "Z")
survfit(Surv(Yobs,status)~Z)
dev.off()

thetatrue<-list(eta.0 = eta$int,
                eta.1 = eta$x1st,
                eta.2 = eta$x2,
                eta.5 = eta$x5,
                alpha.d = alpha.d, beta.d = par.scale.d$beta.d,
                eta.d1 = par.scale.d$eta.d1,
                eta.d2 = par.scale.d$eta.d2,
                eta.d5 = par.scale.d$eta.d5,
                alpha.y1nd = theta$alpha.y1nd, beta.y1nd = theta$beta.y1nd,
                alpha.y1d = theta$alpha.y1d, beta.y1d = theta$beta.y1d,
                alpha.y0nd = theta$alpha.y0nd, beta.y0nd = theta$beta.y0nd,
                alpha.y0d = theta$alpha.y0d, beta.y0d = theta$beta.y0d,
                eta.y1 = theta$eta.y1,
                eta.y2 = theta$eta.y2,
                eta.y5 = theta$eta.y5,
                lambda = theta$lambda,
                pi=pix)

save(thetatrue, file = paste("thetatrue1_",seed,".RData",sep=""))
rm(theta)

write.table(data,paste("simdata1_",seed,".txt",sep=""))

realdata<-data.frame(Z=Z,Y1=Y1,Y0=Y0,x1=x1st,x2=x2,x5=x5,ID=ID,D=D1)
write.table(realdata,paste("realdata1_",seed,".txt",sep=""))