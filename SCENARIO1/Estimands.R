ce <- function(theta, d,complete.data, scenario){
  
  x1<-complete.data$x1
  x2<-complete.data$x2
  x3<-complete.data$x3
  ID<-complete.data$ID
  D<-complete.data$D

  EY1<-NULL
  EY0<-NULL
  
  for(i in nrow(complete.data)){
    EY1 <- sum(EY1,exp(-(theta$beta.y1nd+
                           theta$eta.y1*x1[i]+
                           theta$eta.y2*x2[i]+
                           theta$eta.y3*x3[i])/theta$alpha.y1nd)*gamma(1+{1/theta$alpha.y1nd})*theta$pi[i])
    EY0 <- sum(EY0,exp(-(theta$beta.y0nd+
                           theta$eta.y1*x1[i]+
                           theta$eta.y2*x2[i]+
                           theta$eta.y3*x3[i])/theta$alpha.y0nd)*gamma(1+{1/theta$alpha.y0nd})*theta$pi[i])
  }
  
  pG<-sum(theta$pi)
  
  EY1.nd <- EY1/pG
  EY0.nd <- EY0/pG
  
  EY0.nd_save <- EY0.nd
  EY0.nd_save2 <- EY0
  
  rm(EY1,EY0,pG)
  
  ace.nd<- EY1.nd-EY0.nd
    
  ace.d <- ace.fd <- NULL
  
  for(h in 1:length(d)){
    
    EY0 <- EY1 <- ph <- NULL
    
    for(i in nrow(complete.data)){
      ph <- rbind(ph, dweib(d[h],
                            a=theta$alpha.d,
                            b=theta$beta.d+
                              theta$eta.y1*x1[i]+
                              theta$eta.y2*x2[i]+
                              theta$eta.y3*x5[i]))
      
      EY1 <- sum(EY1,Weib.Trunc.Moments(order=1,
                                        shape=theta$alpha.y1d,
                                        scale=exp(-(theta$beta.y1d+
                                                      theta$eta.y1*x1[i]+
                                                      theta$eta.y2*x2[i]+
                                                      theta$eta.y3*x3[i]+
                                                      theta$lambda*log(d[h]))/theta$alpha.y1d),
                                        left=d[h]))
      
      if(scenario==1){
        EY0 <- sum(EY0,(exp(-(theta$beta.y0d+
                                theta$lambda*log(d[h])+
                                theta$eta.y1*x1[i]+
                                theta$eta.y2*x2[i]+
                                theta$eta.y3*x3[i])/theta$alpha.y0d)*gamma(1+{1/theta$alpha.y0d})))
      }else{
        EY0 <- sum(EY0,Weib.Trunc.Moments(order=1,
                                          shape=theta$alpha.y0d,
                                          scale=exp(-(theta$beta.y0d+
                                                        theta$eta.y1*x1[i]+
                                                        theta$eta.y2*x2[i]+
                                                        theta$eta.y3*x3[i]+
                                                        theta$lambda*log(d[h]))/theta$alpha.y0d),
                                          left=d[h]))
      }

    }

    pG<-sum(ph)
        
    ace.d[h] <- sum((EY1-EY0)*(ph/pG))
    
    ace.fd[h] <- ace.d[h]*pG*(1/nrow(complete.data))
}
  
  pND<-sum(ID==1)/nrow(complete.data)
  pD<-1-pND
  
  aceoverall<-ace.nd*pND + sum(ace.fd)*pD
  
  list(ace.nd=ace.nd, ace.d=ace.d, ace.fd=ace.fd,
       aceoverall=aceoverall)
}

dce <- function(theta, d, y, complete.data){

  x1<-complete.data$x1
  x2<-complete.data$x2
  x3<-complete.data$x3
  ID<-complete.data$ID
  D<-complete.data$D

  dce.nd1<-dce.nd0<-NULL
  
  for(i in which(ID==1)){
    dce.nd1 <- rbind(dce.nd1,
                     Sweib(y, theta$alpha.y1nd, 
                           theta$beta.y1nd+
                             theta$eta.y1*x1[i]+
                             theta$eta.y2*x2[i]+
                             theta$eta.y3*x3[i])*theta$pi[i] ) 
    dce.nd0 <- rbind(dce.nd0,
                     Sweib(y, theta$alpha.y0nd, 
                           theta$beta.y0nd+
                             theta$eta.y1*x1[i]+
                             theta$eta.y2*x2[i]+
                             theta$eta.y3*x3[i])*theta$pi[i])
  }
  
  dce.nd1s<-colSums(dce.nd1)
  dce.nd0s<-colSums(dce.nd0)
  
  pG<-sum(theta$pi[ID==1])
  
  dcend1 <- dce.nd1s/pG
  dcend0 <- dce.nd0s/pG
  
  rm(dce.nd1,dce.nd0,pG)
  
  dce.nd<- dcend1-dcend0
  dce.nd<-t(as.vector(dce.nd))
  
  dce.d <- NULL
  
  for(h in 1:length(d)){
    
    dce.d0 <- dce.d1 <- NULL
    
    for(i in which(ID==0)){
      dce.d1 <- rbind(dce.d1,
                      Stweib(y,
                            sh=theta$alpha.y1d,
                            sc=exp(-(theta$beta.y1d+
                                       theta$eta.y1*x1[i]+
                                       theta$eta.y2*x2[i]+
                                       theta$eta.y3*x3[i]+
                                       theta$lambda*log(d[h]))/theta$alpha.y1d),
                            a=d[h])*(1-theta$pi[i]))
      dce.d0 <- rbind(dce.d0,
                      Sweib(y,
                            theta$alpha.y0d,
                            theta$beta.y0d+
                              theta$lambda*log(d[h])+
                              theta$eta.y1*x1[i]+
                              theta$eta.y2*x2[i]+
                              theta$eta.y3*x3[i])*(1-theta$pi[i]))
      
    }
    
    pG<-sum(1-theta$pi[ID==0])
    
    dced0 <- colSums(dce.d0)/pG
    dced1 <- colSums(dce.d1)/pG

    dce.d <- rbind(dce.d,dced1-dced0)
    
    rm(dced0,dced1)

  }
  
  # list(dce.nd=dce.nd)
  list(dce.nd=dce.nd, dce.d=dce.d)

}
