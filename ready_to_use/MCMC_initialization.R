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

theta.prior <- list()

theta.prior <- c(theta.prior, 
                 list(mu.0 = .6, sigma2.0 = s20,
                    a.d=a0, b.d=b0, 
                    mu.d=m0, sigma2.d= s20,
                    a.y1nd=a0, b.y1nd=b0, 
                    mu.y1nd=m0, sigma2.y1nd=s20,
                    a.y1d=a0,  b.y1d=b0, 
                    mu.y1d=m0, sigma2.y1d=s20,
                    a.y0nd=aa0, b.y0nd=bb0, 
                    mu.y0nd=m0, sigma2.y0nd=s20,
                    a.y0d=aa0, b.y0d=bb0, 
                    mu.y0d=m0, sigma2.y0d=s20,
                    mu.delta=m.start,  sigma2.delta=s20))

mu_covariates<-c()
mu_d_covariates<-c()
mu_y_covariates<-c()
sigma_covariates<-c()
sigma_d_covariates<-c()
sigma_y_covariates<-c()

for(cov in 1:n_cov){
  mu_covariates[cov]<-paste("mu.",cov,sep="")
  mu_d_covariates[cov]<-paste("mu.d",cov,sep="")
  mu_y_covariates[cov]<-paste("mu.y",cov,sep="")
  sigma_covariates[cov]<-paste("sigma2.",cov,sep="")
  sigma_d_covariates[cov]<-paste("sigma2.d",cov,sep="")
  sigma_y_covariates[cov]<-paste("sigma2.y",cov,sep="")
}

for(cov in mu_covariates){
  theta.prior[[cov]]<-m0
}
for(cov in mu_d_covariates){
  theta.prior[[cov]]<-m0
}
for(cov in mu_y_covariates){
  theta.prior[[cov]]<-m0
}
for(cov in sigma_covariates){
  theta.prior[[cov]]<-s20.du
}
for(cov in sigma_d_covariates){
  theta.prior[[cov]]<-s20.du
}
for(cov in sigma_y_covariates){
  theta.prior[[cov]]<-s20.du
}

#### Set starting values
par.start <- list()
par.start <- c(par.start, list(mu.0 = m.start, sigma.0 = sd.start,
                  alpha.d=1.2,  beta.d  = -2,
                  shape.y1nd = a.start, scale.y1nd = scale.start, 
                  mu.y1nd = m.start_y, sigma.y1nd = sd.start,
                  shape.y1d  = a.start, scale.y1d = scale.start, 
                  mu.y1d = m.start_y,  sigma.y1d = sd.start,
                  shape.y0nd = a.start, scale.y0nd = scale.start, 
                  mu.y0nd = m.start_y, sigma.y0nd = sd.start, 
                  shape.y0d  = a.start,  scale.y0d = scale.start, 
                  mu.y0d  = m.start_y,   sigma.y0d  = sd.start,
                  mu.delta = m.start,  sigma.delta = sd.start))

for(cov in mu_covariates){
  par.start[[cov]]<-m.start
}
for(cov in mu_d_covariates){
  par.start[[cov]]<-m.start
}
for(cov in mu_y_covariates){
  par.start[[cov]]<-m.start
}
for(cov in sigma_covariates){
  par.start[[cov]]<-sd.start
}
for(cov in sigma_d_covariates){
  par.start[[cov]]<-sd.start
}
for(cov in sigma_y_covariates){
  par.start[[cov]]<-sd.start
}

#### Set parameters of the proposal distributions
proposal <- list()
proposal <- c(proposal, list(sd.0 = 0.15, 
                             scale.d = 0.02, sd.d = 0.3,
                             scale.y1nd=0.01, sd.y1nd =0.2,
                             scale.y1d=0.02,  sd.y1d =0.4,
                             scale.y0nd=0.015, sd.y0nd = 0.2,
                             scale.y0d=0.03,  sd.y0d = 0.3,
                             sd.delta = 0.2))

sd_covariates<-c()
sd_d_covariates<-c()
sd_y_covariates<-c()

for(cov in 1:n_cov){
  sd_covariates[cov]<-paste("sd.",cov,sep="")
  sd_d_covariates[cov]<-paste("sd.d",cov,sep="")
  sd_y_covariates[cov]<-paste("sd.y",cov,sep="")
}

for(cov in sigma_covariates){
  proposal[[cov]]<-0.15
}

for(cov in sigma_d_covariates){
  proposal[[cov]]<-0.2
}
for(cov in sigma_y_covariates){
  proposal[[cov]]<-0.1
}

rm(cov,mu_covariates,mu_d_covariates,mu_y_covariates,
   sd_covariates,sd_d_covariates,sd_y_covariates,
   sigma_covariates,sigma_d_covariates,sigma_y_covariates)

