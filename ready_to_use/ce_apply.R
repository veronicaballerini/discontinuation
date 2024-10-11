ce_apply <- function(x){
  
  d<-seq(1,10) 
  ID<-x$ID
  D<-x$D
  xs<-cbind(1,as.matrix(observed.data[,grepl("x",names(observed.data))]))
  
  etas_pos<-grepl("eta.",names(x))+
    grepl("beta",names(x))+
    grepl("eta.d",names(x))+
    grepl("eta.y",names(x))
  
  etas<-NULL
  for(i in which(etas_pos==1)){
    etas<-rbind(etas,as.numeric(x[i]))
  }
  etas<-as.matrix(etas)
  
  pi<-exp(xs%*%etas)/(1+exp(xs%*%etas))
  
  etas_y_pos<-grepl("beta.y",names(x))+grepl("eta.y",names(x))
  
  etas_y<-NULL
  for(i in which(etas_y_pos==1)){
    etas_y<-rbind(etas_y,as.numeric(x[i]))
  }
  etas_y<-as.matrix(etas_y)
  
  etas_d_pos<-grepl("beta.d",names(x))+grepl("eta.d",names(x))
  
  etas_d<-NULL
  for(i in which(etas_d_pos==1)){
    etas_d<-rbind(etas_d,as.numeric(x[i]))
  }
  etas_d<-as.matrix(etas_d)
  
  EY1<-NULL
  EY0<-NULL
  
  for(i in which(ID==1)){
    EY1 <- sum(EY1,exp(-(x$beta.y1nd+xs[i,-1]%*%etas_y)/x$alpha.y1nd)*gamma(1+{1/x$alpha.y1nd})*pi[i])
    EY0 <- sum(EY0,exp(-(x$beta.y0nd+xs[i,-1]%*%etas_y)/x$alpha.y0nd)*gamma(1+{1/x$alpha.y0nd})*pi[i])
  }
  
  pG<-sum(pi[ID==1])
  
  EY1.nd <- EY1/pG
  EY0.nd <- EY0/pG
  
  rm(EY1,EY0,pG)
  
  ace.nd<- EY1.nd-EY0.nd
  
  pG<-sum(1-pi[ID==0])
  
  ace.d <- ace.fd <- NULL
  
  for(h in 1:length(d)){
    
    EY0 <- EY1 <- NULL
    
    for(i in which(ID==0)){
      EY1 <- sum(EY1,Weib.Trunc.Moments(order=1,
                                        shape=x$alpha.y1d,
                                        scale=exp(-(x$beta.y1d+xs[i,-1]%*%etas_y+
                                                      x$lambda*log(d[h]))/x$alpha.y1d),
                                        left=d[h])*(1-pi[i]))
      EY0 <- sum(EY0,(exp(-(x$beta.y0d+
                              x$lambda*log(d[h])+
                              xs[i,-1]%*%etas_y)/x$alpha.y0d)*gamma(1+{1/x$alpha.y0d}))*(1-pi[i]))
      
    }
    
    EY1.d <- EY1/pG
    EY0.d <- EY0/pG
    
    rm(EY0,EY1)
    
    ace.d[h] <- EY1.d-EY0.d
    
    ace.fd[h] <- ace.d[h]*dweib(h,
                                a=x$alpha.d,
                                b=x$beta.d+xs[i,-1]%*%etas_d)
    
    rm(EY1.d,EY0.d)
  }
  
  pND<-sum(ID==1)/nrow(observed.data)
  pD<-1-pND
  
  aceoverall<-ace.nd*pND + sum(ace.fd)*pD
  
  list(ace.nd=ace.nd, ace.d=ace.d, ace.fd = ace.fd,
       aceoverall=aceoverall)
}

dce_apply <- function(x){
  
  d<-seq(1,6)
  y<-seq(1,30)
  ID<-x$ID
  D<-x$D
  xs<-cbind(1,as.matrix(observed.data[,grepl("x",names(observed.data))]))
  
  etas_pos<-grepl("eta.",names(x))+
    grepl("beta",names(x))+
    grepl("eta.d",names(x))+
    grepl("eta.y",names(x))
  
  etas<-NULL
  for(i in which(etas_pos==1)){
    etas<-rbind(etas,as.numeric(x[i]))
  }
  etas<-as.matrix(etas)
  
  pi<-exp(xs%*%etas)/(1+exp(xs%*%etas))
  
  etas_y<-NULL
  for(i in which(etas_y_pos==1)){
    etas_y<-rbind(etas_y,as.numeric(x[i]))
  }
  etas_y<-as.matrix(etas_y)
  
  etas_d_pos<-grepl("beta.d",names(x))+grepl("eta.d",names(x))
  
  etas_d<-NULL
  for(i in which(etas_d_pos==1)){
    etas_d<-rbind(etas_d,as.numeric(x[i]))
  }
  etas_d<-as.matrix(etas_d)
  
  dce.nd1<-dce.nd0<-NULL
  
  for(i in which(ID==1)){
    dce.nd1 <- rbind(dce.nd1,
                     Sweib(y, x$alpha.y1nd, 
                           x$beta.y1nd+xs[i,-1]%*%etas_y)*pi[i] ) 
    dce.nd0 <- rbind(dce.nd0,
                     Sweib(y, x$alpha.y0nd,
                           x$beta.y0nd+xs[i,-1]%*%etas_y)*pi[i])
  }
  
  dce.nd1s<-colSums(dce.nd1)
  dce.nd0s<-colSums(dce.nd0)
  
  pG<-sum(pi[ID==1])
  
  dcend1 <- dce.nd1s/pG
  dcend0 <- dce.nd0s/pG
  
  rm(dce.nd1,dce.nd0,pG)
  
  dce.nd<- dcend1-dcend0
  dce.nd<-t(as.vector(dce.nd))
  
  dce.d <- NULL
  
  for(h in 1:length(d)){
    
    dce.d0 <- dce.d1 <- NULL
    
    for(i in which(ID==0)){
      # dce.d1 <- rbind(dce.d1,
      #                 Stweib(y,
      #                         sh=x$alpha.y1d,
      #                         sc=exp(-(x$beta.y1d+xs[i,-1]%*%etas_y+
      #                                    x$delta*log(d[h]))/x$alpha.y1d),
      #                         a=d[h])*(1-pi[i]))
      dce.d1 <- rbind(dce.d1,
                      Stweib(y,
                             a=x$alpha.y1d,
                             b=x$beta.y1d+xs[i,-1]%*%etas_y+
                               x$lambda*log(d[h]),
                             l=d[h])*(1-pi[i]))
      
      dce.d0 <- rbind(dce.d0,
                      (Sweib(y,
                             x$alpha.y0d,
                             x$beta.y0d+
                               x$lambda*log(d[h])+xs[i,-1]%*%etas_y))*(1-pi[i]))
      
    }
    
    pG<-sum(1-pi[ID==0])
    
    dced0 <- colSums(dce.d0)/pG
    dced1 <- colSums(dce.d1)/pG
    
    dce.d <- rbind(dce.d,dced1-dced0)
    
    rm(dced0,dced1)
    
  }
  
  list(dce.nd=dce.nd, dce.d=dce.d)
  
}
