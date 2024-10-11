mcmc.tdiscontinue <- function(niter, burn, thin, par.start, 
                              theta.prior, observed.data, n_cov = 3){
  
  prop_etas_pos<-grepl("sigma2.",names(proposal))+
    grepl("sigma2.d",names(proposal))+
    grepl("sigma2.y",names(proposal))
  proposal_eta<-NULL
  proposal_eta<-cbind(proposal_eta,proposal$sd.0)
  for(i in which(prop_etas_pos==1)){
    proposal_eta<-rbind(proposal_eta,as.numeric(proposal[i]))
  }
  
  prop_etas_d_pos<-grepl("sigma2.d",names(proposal))
  proposal_eta_d<-NULL
  for(i in which(prop_etas_d_pos==1)){
    proposal_eta_d<-rbind(proposal_eta_d,as.numeric(proposal[i]))
  }
  
  prop_etas_y_pos<-grepl("sigma2.y",names(proposal))
  proposal_eta_y<-NULL
  for(i in which(prop_etas_y_pos==1)){
    proposal_eta_y<-rbind(proposal_eta_y,as.numeric(proposal[i]))
  }

  theta.start <- list()
  theta.start <- c(theta.start,
                   list(
    eta.0      = rnorm(1, mean = par.start$mu.0, sd = par.start$sigma.0),
    alpha.d    = par.start$alpha.d + abs(rnorm(1,0,0.1)),
    beta.d     = par.start$beta.d + rnorm(1,0,0.1),
    alpha.y1nd = par.start$shape.y1nd*par.start$scale.y1nd,
    beta.y1nd  = rnorm(1, mean  = par.start$mu.y1nd, sd = par.start$sigma.y1nd),
    alpha.y1d  = par.start$shape.y1d*par.start$scale.y1d,
    beta.y1d   = rnorm(1, mean  = par.start$mu.y1d, sd = par.start$sigma.y1d),
    alpha.y0nd = par.start$shape.y0nd*par.start$scale.y0nd,
    beta.y0nd  = rnorm(1, mean  = par.start$mu.y0nd, sd = par.start$sigma.y0nd),
    alpha.y0d  = par.start$shape.y0d*par.start$scale.y0d,
    beta.y0d   = rnorm(1, mean  = par.start$mu.y0d, sd = par.start$sigma.y0d),
    delta     = rnorm(1, mean  = par.start$mu.delta, sd = par.start$sigma.delta)
  ))
  
  covariates<-c()
  d_covariates<-c()
  y_covariates<-c()
  
  for(cov in 1:n_cov){
    covariates[cov]<-paste("eta.",cov,sep="")
    d_covariates[cov]<-paste("eta.d",cov,sep="")
    y_covariates[cov]<-paste("eta.y",cov,sep="")
  }
  
  c<-NULL
  for(cov in covariates){
    c<-sum(c,1)
    theta.start[[cov]]<-rnorm(1, mean  = unlist(par.start[paste("mu.",c,sep="")]), 
                              sd = unlist(par.start[paste("sigma2.",c,sep="")]))
  }
  
  c<-NULL
  for(cov in d_covariates){
    c<-sum(c,1)
    theta.start[[cov]]<-rnorm(1, mean  = unlist(par.start[paste("mu.d",c,sep="")]), 
                              sd = unlist(par.start[paste("sigma2.d",c,sep="")]))
  }
  
  c<-NULL
  for(cov in y_covariates){
    c<-sum(c,1)
    theta.start[[cov]]<-rnorm(1, mean  = unlist(par.start[paste("mu.y",c,sep="")]), 
                              sd = unlist(par.start[paste("sigma2.y",c,sep="")]))
  }
  
  rm(covariates,d_covariates,y_covariates)
  
  etas_pos<-grepl("eta.",names(theta.start))+
    grepl("beta",names(theta.start))+
    grepl("eta.d",names(theta.start))+
    grepl("eta.y",names(theta.start))
  
  eta.start<-NULL

  for(i in which(etas_pos==1)){
    eta.start<-rbind(eta.start,as.numeric(theta.start[i]))
  }
  
  eta.start<-as.matrix(eta.start)
  
  # x<-cbind(1,as.matrix(observed.data[,grepl("x",names(observed.data))]))
  
  pi <- exp(x%*%eta.start)/(1+exp(x%*%eta.start))
  
  ID.start<- D.start <- NULL
  ID.start[observed.data$Z==1 & 
             observed.data$RD.obs==0 & 
             observed.data$RY==1]<-1
  ID.start[observed.data$Z==1 & 
             observed.data$RD.obs==1]<-0
  ID.start[observed.data$Z==1 & 
             observed.data$RD.obs==0 & 
             observed.data$RY==0]<-
    rbinom(sum((observed.data$Z)*(1-observed.data$RD.obs)*(1-observed.data$RY)),
           1,
           pi[observed.data$Z==1 & observed.data$RD.obs==0 & observed.data$RY==0])
  ID.start[observed.data$Z==0]<-rbinom(sum(observed.data$Z==0),
                                       1,
                                       pi[observed.data$Z==0])

  D.start[ID.start==1]<-0
  D.start[ID.start==0 & observed.data$Z==1]  <- 
    observed.data$Dobs[ID.start==0 & observed.data$Z==1]
  
  etas.d_pos<-grepl("eta.d",names(theta.start))+grepl("beta.d",names(theta.start))
  
  eta.dstart<-NULL
  
  for(i in which(etas.d_pos==1)){
    eta.dstart<-rbind(eta.dstart,as.numeric(theta.start[i]))
  }
  
  eta.dstart<-as.matrix(eta.dstart)

  D.start[ID.start==0 & observed.data$Z==0] <-
    rweibull(sum(ID.start==0 & observed.data$Z==0), 
             shape=theta.start$alpha.d,
             scale=exp(-(x[ID.start==0 & observed.data$Z==0,]%*%matrix(c(theta.start$beta.d,eta.dstart)))/theta.start$alpha.d))

  complete.data <- data.frame(cens=observed.data$f.up, Z=observed.data$Z,
                              ID=ID.start, RD=observed.data$RD.obs, D=D.start, 
                              RY=observed.data$RY, Y=observed.data$Y)
  
  # pos<-1
  # for(i in names(observed.data[,x_names==1])){
  #   pos<-pos+1
  #   complete.data[i]<-x[,pos]
  # }
  
  theta <- list()

  pos<-NULL
  for(cov in names(theta.start[which(etas_pos==1)])){
    pos<-sum(pos,1)
    theta[[cov]]<-eta.start[pos,1]
  }
  
  theta <- c(theta,list(alpha.d=theta.start$alpha.d,beta.d=theta.start$beta.d))
  
  pos<-NULL
  for(cov in names(theta.start[which(etas.d_pos==1)])){
    pos<-sum(pos,1)
    theta[[cov]]<-eta.dstart[pos,1]
  }
  
  theta<-c(theta,list(alpha.y1nd=theta.start$alpha.y1nd,beta.y1nd=theta.start$beta.y1nd,
                alpha.y1d=theta.start$alpha.y1d,beta.y1d=theta.start$beta.y1d,
                alpha.y0nd=theta.start$alpha.y0nd,beta.y0nd=theta.start$beta.y0nd,
                alpha.y0d=theta.start$alpha.y0d,beta.y0d=theta.start$beta.y0d))
  
  etas.y_pos<-grepl("eta.y",names(theta.start))+grepl("beta.y",names(theta.start))
  
  eta.ystart<-NULL
  
  for(i in which(etas.y_pos==1)){
    eta.ystart<-rbind(eta.ystart,as.numeric(theta.start[i]))
  }
  
  eta.ystart<-as.matrix(eta.ystart)
  
  pos<-NULL
  for(cov in names(theta.start[which(etas.y_pos==1)])){
    pos<-sum(pos,1)
    theta[[cov]]<-eta.ystart[pos,1]
  }
  
  theta <- c(theta, list(pi=pi, delta=theta.start$delta))

  thetanop<-theta
  thetanop[["pi"]]<-NULL
  npar  <- length(unlist(thetanop))
  jump  <- matrix(NA, niter, npar)
  Theta <- matrix(NA, niter, npar)
  DID <- list()
  
  colnames(jump)  <- names(unlist(thetanop))
  colnames(Theta) <- names(unlist(thetanop))
  
  nojump<-NA
  
  rm(eta.dstart,eta.start,eta.ystart,c,cov,etas_pos,etas.d_pos,etas.y_pos)
  
  etas_pos<-grepl("eta.",names(theta))+
    grepl("beta",names(theta))+
    grepl("eta.d",names(theta))+
    grepl("eta.y",names(theta))
  
  etas_d_pos<-grepl("beta.d",names(theta))+
    grepl("eta.d",names(theta))
  
  etas_y_pos<-grepl("beta.y",names(theta))+
    grepl("eta.y",names(theta))
  
  for(j in 1:niter){
    
    ID.old <- complete.data$ID
    D.old <- complete.data$D
    if(j%%1000==0){
      print(j)
    }
    
      discontinuedata <- da.discontinuation(theta, theta.prior, observed.data, 
                                            ID.old, D.old)
      complete.data$ID <- discontinuedata$ID
      complete.data$D <- discontinuedata$D
      
      rm(ID.old, D.old)
    
    ###eta0
      etas<-NULL
      for(i in which(etas_pos==1)){
        etas<-rbind(etas,as.numeric(theta[i]))
      }
      etas<-as.matrix(etas)
      
    eta <- rmvnorm(1, mean=etas, sigma=diag(c(proposal_eta^2),nrow=n_cov+1))
    eta<-matrix(eta,ncol=1)
    
    thetaprop <- theta
    
    i<-NULL
    for(pos in which(etas_pos==1)){
      i<-sum(i,1)
      thetaprop[[pos]]<-eta[i]
    }
    
    thetaprop$pi <- as.vector(exp(x%*%eta)/(1+exp(x%*%eta)))
    
    log.post.eta.c   <- logpost(thetaprop,
                                theta.prior,
                                complete.data)
    log.post.eta.old <- logpost(theta,
                                theta.prior,
                                complete.data)
    
    pr.num <-  log.post.eta.c    
    pr.den <-  log.post.eta.old   
    pr <- exp(pr.num - pr.den)
    
    u<- runif(1)
    if(u<=pr  & is.nan(pr)==FALSE){
      theta <- thetaprop
    }
    
    jump[j,c(1:(n_cov+1))]<- as.numeric(u<=pr & is.nan(pr)==FALSE)
    rm(thetaprop, u, pr.num, pr.den, pr, log.post.eta.c, log.post.eta.old)
    
    ###alpha.d
    alpha.d <-  rgamma(1, shape=theta$alpha.d/proposal$scale.d, scale=proposal$scale.d)
    
    thetaprop <- theta
    thetaprop$alpha.d <- alpha.d
    
    log.post.alpha.d.c   <-	logpost(thetaprop,theta.prior, complete.data)
    log.post.alpha.d.old <-	logpost(theta,theta.prior, complete.data)
    
    pr.num <- log.post.alpha.d.c    + log(dgamma(theta$alpha.d, shape=alpha.d/proposal$scale.d,       scale=proposal$scale.d))
    pr.den <- log.post.alpha.d.old  + log(dgamma(alpha.d,       shape=theta$alpha.d/proposal$scale.d, scale=proposal$scale.d))
    pr <- exp(pr.num - pr.den)
    
    u<- runif(1)
    theta$alpha.d <- ifelse((u<=pr & is.nan(pr)==FALSE), alpha.d, theta$alpha.d)
    jump[j,n_cov+2]<- as.numeric(u<=pr & is.nan(pr)==FALSE)
    rm(alpha.d, thetaprop, u, pr.num, pr.den, pr, log.post.alpha.d.c, log.post.alpha.d.old)
    
    ###beta.d
    beta.d <-  rnorm(1,theta$beta.d, proposal$sd.d)
    
    thetaprop <- theta
    thetaprop$beta.d <- beta.d
    
    log.post.beta.d.c   <- logpost(thetaprop, theta.prior, complete.data)
    log.post.beta.d.old <- logpost(theta, theta.prior, complete.data)
    
    pr.num <-  log.post.beta.d.c 
    pr.den <-  log.post.beta.d.old  
    
    pr <- exp(pr.num - pr.den)
    
    u<- runif(1)
    theta$beta.d <- ifelse((u<=pr & is.nan(pr)==FALSE), beta.d, theta$beta.d)
    jump[j,n_cov+3]<- as.numeric(u<=pr & is.nan(pr)==FALSE)
    rm(beta.d, thetaprop, u, pr.num, pr.den, pr, log.post.beta.d.c, log.post.beta.d.old)

    ###eta.d
    etas_d<-NULL
    for(i in which(etas_d_pos==1)){
      etas_d<-rbind(etas_d,as.numeric(theta[i]))
    }
    etas_d<-as.matrix(etas_d)
    
    eta_d <- rmvnorm(1, mean=etas_d, sigma=diag(c(proposal_eta_d^2), nrow=n_cov))
    eta_d<-matrix(eta_d,ncol=1)
    
    thetaprop <- theta
    
    i<-NULL
    for(pos in which(etas_d_pos==1)){
      i<-sum(i,1)
      thetaprop[[pos]]<-eta_d[i]
    }
    
    log.post.eta.d.c   <- logpost(thetaprop, theta.prior, complete.data)
    log.post.eta.d.old <- logpost(theta, theta.prior, complete.data)
    
    pr.num <-  log.post.eta.d.c   
    pr.den <-  log.post.eta.d.old 
    pr <- exp(pr.num - pr.den)
    
    u<- runif(1)
    if(u<=pr  & is.nan(pr)==FALSE){
      theta <- thetaprop
    }
    
    jump[j,c((n_cov+4):(n_cov+4+n_cov-1))]<- as.numeric(u<=pr & is.nan(pr)==FALSE)
    rm(eta_d, thetaprop, u, pr.num, pr.den, pr, log.post.eta.d.c, log.post.eta.d.old)
 
    ###alpha.y1nd
    alpha.y1nd <-  rgamma(1, shape=theta$alpha.y1nd/proposal$scale.y1nd, scale=proposal$scale.y1nd)
    
    thetaprop <- theta
    thetaprop$alpha.y1nd <- alpha.y1nd
    
    log.post.alpha.y1nd.c   <- logpost(thetaprop, theta.prior, complete.data)
    log.post.alpha.y1nd.old <- logpost(theta, theta.prior, complete.data)
    
    pr.num <- log.post.alpha.y1nd.c    + log(dgamma(theta$alpha.y1nd, 
                                                     shape=alpha.y1nd/proposal$scale.y1nd,       
                                                     scale=proposal$scale.y1nd))
    pr.den <- log.post.alpha.y1nd.old  + log(dgamma(alpha.y1nd,       
                                                     shape=theta$alpha.y1nd/proposal$scale.y1nd, 
                                                     scale=proposal$scale.y1nd))
    
    pr <- exp(pr.num - pr.den)
    
    u<- runif(1)
    theta$alpha.y1nd <- ifelse((u<=pr & is.nan(pr)==FALSE), alpha.y1nd, theta$alpha.y1nd)
    jump[j,n_cov+4+n_cov]<- as.numeric(u<=pr & is.nan(pr)==FALSE)
    rm(alpha.y1nd, thetaprop, u, pr.num, pr.den, pr, log.post.alpha.y1nd.c, log.post.alpha.y1nd.old)
    
    ###beta.y1nd
    beta.y1nd <-  rnorm(1,theta$beta.y1nd, proposal$sd.y1nd)
    
    thetaprop <- theta
    thetaprop$beta.y1nd <- beta.y1nd
    
    log.post.beta.y1nd.c   <- logpost(thetaprop, theta.prior, complete.data)
    log.post.beta.y1nd.old <- logpost(theta, theta.prior, complete.data)
    
    pr.num <-  log.post.beta.y1nd.c   
    pr.den <-  log.post.beta.y1nd.old 
    pr <- exp(pr.num - pr.den)
    
    u<- runif(1)
    theta$beta.y1nd <- ifelse((u<=pr & is.nan(pr)==FALSE), beta.y1nd, theta$beta.y1nd)
    jump[j,n_cov+4+n_cov+1]<- as.numeric(u<=pr & is.nan(pr)==FALSE)
    rm(beta.y1nd, thetaprop, u, pr.num, pr.den, pr, log.post.beta.y1nd.c, log.post.beta.y1nd.old)
    
    ###alpha.y1d
    alpha.y1d <-  rgamma(1, shape=theta$alpha.y1d/proposal$scale.y1d, scale=proposal$scale.y1d)
    
    thetaprop <- theta
    thetaprop$alpha.y1d <- alpha.y1d
    
    log.post.alpha.y1d.c   <- logpost(thetaprop,theta.prior, complete.data)
    log.post.alpha.y1d.old <- logpost(theta, theta.prior, complete.data)
    
    pr.num <- log.post.alpha.y1d.c    + log(dgamma(theta$alpha.y1d, 
                                                    shape=alpha.y1d/proposal$scale.y1d,        
                                                    scale=proposal$scale.y1d))
    pr.den <- log.post.alpha.y1d.old  + log(dgamma(alpha.y1d,       
                                                    shape=theta$alpha.y1d/proposal$scale.y1d,  
                                                    scale=proposal$scale.y1d))
    pr <- exp(pr.num - pr.den)
    
    u<- runif(1)
    theta$alpha.y1d <- ifelse((u<=pr & is.nan(pr)==FALSE), alpha.y1d, theta$alpha.y1d)
    jump[j,n_cov+4+n_cov+2]<- as.numeric(u<=pr & is.nan(pr)==FALSE)
    rm(alpha.y1d, thetaprop, u, pr.num, pr.den, pr, log.post.alpha.y1d.c, log.post.alpha.y1d.old)
    
    ###beta.y1d
    beta.y1d <-  rnorm(1,theta$beta.y1d, proposal$sd.y1d)

    thetaprop <- theta
    thetaprop$beta.y1d <- beta.y1d
    
    log.post.beta.y1d.c   <- logpost(thetaprop, theta.prior, complete.data)
    log.post.beta.y1d.old <- logpost(theta, theta.prior, complete.data)
    
    pr.num <-  log.post.beta.y1d.c  
    pr.den <-  log.post.beta.y1d.old 
    
    pr <- exp(pr.num - pr.den)
    u<- runif(1)
    theta$beta.y1d <- ifelse((u<=pr & is.nan(pr)==FALSE), beta.y1d, theta$beta.y1d)
    jump[j,n_cov+4+n_cov+3]<- as.numeric(u<=pr & is.nan(pr)==FALSE)
    rm(beta.y1d, thetaprop, u, pr.num, pr.den, pr, log.post.beta.y1d.c, log.post.beta.y1d.old)
    
    ###alpha.y0nd
    alpha.y0nd <-  rgamma(1, shape=theta$alpha.y0nd/proposal$scale.y0nd, scale=proposal$scale.y0nd)
    

    thetaprop <- theta
    thetaprop$alpha.y0nd <- alpha.y0nd

    log.post.alpha.y0nd.c   <- logpost(thetaprop,theta.prior, complete.data)
    log.post.alpha.y0nd.old <- logpost(theta, theta.prior, complete.data)

    pr.num <- log.post.alpha.y0nd.c   + log(dgamma(theta$alpha.y0nd,
                                                    shape=alpha.y0nd/proposal$scale.y0nd,
                                                    scale=proposal$scale.y0nd))
    pr.den <- log.post.alpha.y0nd.old + log(dgamma(alpha.y0nd,
                                                    shape=theta$alpha.y0nd/proposal$scale.y0nd,
                                                    scale=proposal$scale.y0nd))
    pr <- exp(pr.num - pr.den)

    u<- runif(1)
    theta$alpha.y0nd <- ifelse((u<=pr & is.nan(pr)==FALSE), alpha.y0nd, theta$alpha.y0nd)
    jump[j,n_cov+4+n_cov+4]<- as.numeric(u<=pr & is.nan(pr)==FALSE)
    rm(alpha.y0nd, thetaprop, u, pr.num, pr.den, pr, log.post.alpha.y0nd.c, log.post.alpha.y0nd.old)

    ###beta.y0nd
    beta.y0nd <-  rnorm(1,theta$beta.y0nd, proposal$sd.y0nd)

    thetaprop <- theta
    thetaprop$beta.y0nd <- beta.y0nd

    log.post.beta.y0nd.c   <- logpost(thetaprop, theta.prior, complete.data)
    log.post.beta.y0nd.old <- logpost(theta, theta.prior, complete.data)

    pr.num <-  log.post.beta.y0nd.c
    pr.den <-  log.post.beta.y0nd.old

    pr <- exp(pr.num - pr.den)
    u<- runif(1)
    theta$beta.y0nd <- ifelse((u<=pr& is.nan(pr)==FALSE), beta.y0nd, theta$beta.y0nd)
    jump[j,n_cov+4+n_cov+5]<- as.numeric(u<=pr & is.nan(pr)==FALSE)
    rm(beta.y0nd, thetaprop, u, pr.num, pr.den, pr, log.post.beta.y0nd.c, log.post.beta.y0nd.old)

    ### alpha.y0d
    alpha.y0d <-  rgamma(1, shape=theta$alpha.y0d/proposal$scale.y0d, scale=proposal$scale.y0d)
    
    thetaprop <- theta
    thetaprop$alpha.y0d <- alpha.y0d

    log.post.alpha.y0d.c   <- logpost(thetaprop,  theta.prior, complete.data)
    log.post.alpha.y0d.old <- logpost(theta,   theta.prior, complete.data)

    pr.num <- log.post.alpha.y0d.c    + log(dgamma(theta$alpha.y0d, 
                                                    shape=alpha.y0d/proposal$scale.y0d,       
                                                    scale=proposal$scale.y0d))
    pr.den <- log.post.alpha.y0d.old  + log(dgamma(alpha.y0d,       
                                                    shape=theta$alpha.y0d/proposal$scale.y0d, 
                                                    scale=proposal$scale.y0d))
    pr <- exp(pr.num - pr.den)

    u<- runif(1)

    theta$alpha.y0d <- ifelse((u<=pr & is.nan(pr)==FALSE), alpha.y0d, theta$alpha.y0d)
    jump[j,n_cov+4+n_cov+6]<- as.numeric(u<=pr  & is.nan(pr)==FALSE)
    rm(alpha.y0d, thetaprop, u, pr.num, pr.den, pr, log.post.alpha.y0d.c, log.post.alpha.y0d.old)

    ###beta.y0d
    beta.y0d <-  rnorm(1,theta$beta.y0d, proposal$sd.y0d)

    thetaprop <- theta
    thetaprop$beta.y0d <- beta.y0d

    log.post.beta.y0d.c   <- logpost(thetaprop, theta.prior, complete.data)
    log.post.beta.y0d.old <- logpost(theta,    theta.prior, complete.data)

    pr.num <-  log.post.beta.y0d.c 
    pr.den <-  log.post.beta.y0d.old 
    pr <- exp(pr.num - pr.den)

    u<- runif(1)
    theta$beta.y0d <- ifelse((u<=pr  & is.nan(pr)==FALSE), beta.y0d, theta$beta.y0d)
    jump[j,n_cov+4+n_cov+7]<- as.numeric(u<=pr & is.nan(pr)==FALSE)
    rm(beta.y0d, thetaprop, u, pr.num, pr.den, pr, log.post.beta.y0d.c, log.post.beta.y0d.old)

    #y1
    etas_y<-NULL
    for(i in which(etas_y_pos==1)){
      etas_y<-rbind(etas_y,as.numeric(theta[i]))
    }
    etas_y<-as.matrix(etas_y)
    
    eta_y <- rmvnorm(1, mean=etas_y, sigma=diag(c(proposal_eta_y^2), nrow=n_cov))
    eta_y<-matrix(eta_y,ncol=1)
    
    thetaprop <- theta
    
    i<-NULL
    for(pos in which(etas_y_pos==1)){
      i<-sum(i,1)
      thetaprop[[pos]]<-eta_y[i]
    }
    
    log.post.eta.y.c   <- logpost(thetaprop, theta.prior, complete.data)
    log.post.eta.y.old <- logpost(theta, theta.prior, complete.data)
    
    pr.num <-  log.post.eta.y.c 
    pr.den <-  log.post.eta.y.old 
    pr <- exp(pr.num - pr.den)
    
    u<- runif(1)
    if(u<=pr  & is.nan(pr)==FALSE){
      theta <- thetaprop
    }
    
    jump[j,c((n_cov+4+n_cov+8):(n_cov+4+n_cov+8+n_cov-1))]<- as.numeric(u<=pr & is.nan(pr)==FALSE)
    rm(thetaprop, eta_y, u, pr.num, pr.den, pr, log.post.eta.y.c, log.post.eta.y.old)
    
    ###delta
    delta <-  rnorm(1,theta$delta, proposal$sd.delta)

    thetaprop <- theta
    thetaprop$delta <- delta

    log.post.delta.c   <- logpost(thetaprop,  theta.prior, complete.data)
    log.post.delta.old <- logpost(theta, theta.prior, complete.data)

    pr.num <-  log.post.delta.c  
    pr.den <-  log.post.delta.old 
    pr <- exp(pr.num - pr.den)

    u<- runif(1)
    theta$delta <- ifelse((u<=pr & is.nan(pr)==FALSE), delta, theta$delta)
    jump[j,n_cov+4+n_cov+8+n_cov]<- as.numeric(u<=pr & is.nan(pr)==FALSE)
    rm(delta, thetaprop, u, pr.num, pr.den, pr, log.post.delta.c, log.post.delta.old)
    # theta$delta <- 0
    #print(unlist(theta))
    
    thetanop<-theta
    thetanop$pi<-NULL
    Theta[j,]<- unlist(thetanop)
    
    DID[[j]]<-cbind(complete.data$D, complete.data$ID)
    
    # jump[j,]<-0
    if(j==niter){
      if(sum(colMeans(jump[c(burn:niter),]))==0){
        nojump<-1
        # break
      }
    }
    
  } #End loop over i =1 ... niter
  
  list(nojump=nojump, jump=jump, Theta=Theta, DID=DID)
}#End function mcmc.tdiscontinue
