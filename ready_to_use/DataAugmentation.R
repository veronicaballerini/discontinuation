##################
##DA for D(1)
##################

da.discontinuation<-function(theta, 
                             theta.prior, 
                             observed.data, 
                             ID.old, 
                             D.old){
  
  n      <- nrow(observed.data)
  Z  	   <- observed.data$Z
  RY     <- observed.data$RY
  Y      <- observed.data$Y
  ID.obs <- observed.data$ID.obs
  RD.obs <- observed.data$RD.obs
  Dobs   <- observed.data$Dobs
  
  etas_pos<-grepl("eta.",names(theta))+
    grepl("beta",names(theta))+
    grepl("eta.d",names(theta))+
    grepl("eta.y",names(theta))
  
  etas_d_pos<-grepl("beta.d",names(theta))+
    grepl("eta.d",names(theta))
  
  etas_y_pos<-grepl("beta.y",names(theta))+
    grepl("eta.y",names(theta))
  
  eta.d<-NULL
  for(i in which(etas_d_pos==1)){
    eta.d<-rbind(eta.d,as.numeric(theta[i]))
  }
  eta.d<-as.matrix(eta.d)
  
  etas.y<-NULL
  for(i in which(etas_y_pos==1)){
    etas.y<-rbind(etas.y,as.numeric(theta[i]))
  }
  etas.y<-as.matrix(etas.y)
  
  ID <- D <- rep(NA,n)
  
  #Under treatment
  ID[Z==1 & RD.obs==0 & RY==1]<- 1
  D[Z==1 & RD.obs==0 & RY==1] <- 0
  
  ID[Z==1 & RD.obs==1]        <- 0
  D[Z==1 & RD.obs==1]         <- Dobs[Z==1 & RD.obs ==1]

  id1 <- (Z==1 & RD.obs==0 & RY==0)
  
  num <- theta$pi[id1==1]*Sweib(Y[id1==1], 
                                theta$alpha.y1nd, 
                                theta$beta.y1nd+
                                  x[id1==1,-1]%*%etas.y)
  den <- num + (1-theta$pi[id1==1])*Sweib(Dobs[id1==1],
                                          theta$alpha.d,
                                          theta$beta.d+
                                            x[id1==1,-1]%*%etas.y)
  pr  <- num/den
  ID[id1==1] <- rbinom(sum(id1),1,pr)
  
  # D[Z==1 & ID==0]  <- Dobs[Z==1 & ID==0] ### e di quelli censurati?? ecco sotto

  ii<-as.numeric(id1==1&ID==0)
  dii<-rweibull(sum(ii==1),
                  shape=theta$alpha.d,
                  scale=exp(-(theta$beta.d+
                                x[ii==1,-1]%*%eta.d)/theta$alpha.d))
  D[ii==1]<-dii*(dii>=Dobs[ii==1])+Dobs[ii==1]*(dii<Dobs[ii==1])
  
  rm(ii)
  ii<-as.numeric(id1==1&ID==1)
  D[ii==1]<-0
  
  rm(id1, num, den, pr, ii)
  
  #Under Control
  n0 <- sum(Z==0)
  ID.c <- NULL
  ID.c[Z==1] <- ID[Z==1]
  ID.c[Z==0] <- rbinom(n0,1,theta$pi[Z==0])
  
  D.c <- rep(0, n)
  D.c[Z==1] <- D[Z==1]
  D.c[Z==0 & ID.c==1]<-0
  
  ii <- as.numeric(Z==0 & ID.c==0)
  D.c[ii==1] <- rweibull(sum(ii),
                         shape=theta$alpha.d,
                         scale=exp(-(theta$beta.d+
                                       x[ii==1,-1]%*%eta.d)/theta$alpha.d))
  rm(ii)
  
  complete.data.old <- data.frame(cens=observed.data$f.up, Z=observed.data$Z,
                                  ID=ID.old, RD=observed.data$RD.obs, D=D.old, 
                                  RY=observed.data$RY, Y=observed.data$Y)
  
  complete.data.c   <- data.frame(cens=observed.data$f.up, Z=observed.data$Z,
                                  ID=ID.c, RD=observed.data$RD.obs, D=D.c, 
                                  RY=observed.data$RY, Y=observed.data$Y)
  
  # pos<-1
  # for(i in names(observed.data[,x_names==1])){
  #   pos<-pos+1
  #   complete.data.old[i]<-x[,pos]
  #   complete.data.c[i]<-x[,pos]
  # }
  
  rr <- exp(logpost.i(theta, theta.prior, complete.data.c) - 
              logpost.i(theta, theta.prior, complete.data.old))
  ratio<-rep(0, n)
  
  ii<-as.numeric(Z==0 & ID.old==1  & ID.c ==1)
  ratio[ii==1] <- rr[ii==1]
  
  ii<-as.numeric(Z==0 & ID.old==1 & ID.c==0)
  ratio[ii==1] <- rr[ii==1]*{theta$pi[ii==1]/
    {(1-theta$pi[ii==1])*dweib(D.c[ii==1],
                                 theta$alpha.d,
                                 theta$beta.d+
                                 x[ii==1,-1]%*%eta.d)}}
  
  ii<-as.numeric(Z==0 & ID.old==0 & ID.c==1)
  ratio[ii==1] <- rr[ii==1]*
    {{(1-theta$pi[ii==1])*
        dweib(D.old[ii==1],
              theta$alpha.d,
              theta$beta.d+
                x[ii==1,-1]%*%eta.d)}/theta$pi[ii==1]}
  
  ii<-as.numeric(Z==0 & ID.old==0 & ID.c==0)
  ratio[ii==1] <- rr[ii==1]*
    {dweib(D.old[ii==1],
              theta$alpha.d,
              theta$beta.d+
             x[ii==1,-1]%*%eta.d)/dweib(D.c[ii==1],
                                              theta$alpha.d,
                                              theta$beta.d+
                                          x[ii==1,-1]%*%eta.d)}
  
  ru<-runif(n0)
  ID[Z==0] <- ID.c[Z==0]*(ru <= ratio[Z==0] & is.nan(ratio[Z==0])==FALSE) + 
    ID.old[Z==0]*(ru > ratio[Z==0] | is.nan(ratio[Z==0])==TRUE)  
  D[Z==0] <- D.c[Z==0]*(ru <= ratio[Z==0] & is.nan(ratio[Z==0])==FALSE) + 
    D.old[Z==0]*(ru > ratio[Z==0] | is.nan(ratio[Z==0])==TRUE)   
  
  rm(n0, ID.c, D.c, rr, ratio, ru, complete.data.old, complete.data.c)
  list(ID=ID, D=D)
  
}

