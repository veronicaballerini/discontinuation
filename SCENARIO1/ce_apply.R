ce_apply <- function(x){
  
  d<-seq(1,10) 
  x1<-observed.data$x1st
  x2<-observed.data$x2
  x3<-observed.data$x3
  ID<-x$ID
  D<-x$D
  
  pi<-exp(x$eta.0+
            x$eta.1*x1 + 
            x$eta.2*x2 + 
            x$eta.3*x3)/(1+
                           exp(x$eta.0+
                                 x$eta.1*x1 + 
                                 x$eta.2*x2 + 
                                 x$eta.3*x3))
  
  m1_d<-mean(x1[ID==0])
  m2_d<-mean(x2[ID==0])
  m3_d<-mean(x3[ID==0])
  
  m1_nd<-mean(x1[ID==1])
  m2_nd<-mean(x2[ID==1])
  m3_nd<-mean(x3[ID==1])
  
  m1_ed<-mean(x1[ID==0&D<median(D[ID==0])])
  m2_ed<-mean(x2[ID==0&D<median(D[ID==0])])
  m3_ed<-mean(x3[ID==0&D<median(D[ID==0])])
  
  m1_ld<-mean(x1[ID==0&D>=median(D[ID==0])])
  m2_ld<-mean(x2[ID==0&D>=median(D[ID==0])])
  m3_ld<-mean(x3[ID==0&D>=median(D[ID==0])])
  
  EY1<-NULL
  EY0<-NULL
  
  for(i in which(ID==1)){
    EY1 <- sum(EY1,exp(-(x$beta.y1nd+
                           x$eta.y1*x1[i]+
                           x$eta.y2*x2[i]+
                           x$eta.y3*x3[i])/x$alpha.y1nd)*gamma(1+{1/x$alpha.y1nd})*pi[i])
    EY0 <- sum(EY0,exp(-(x$beta.y0nd+
                           x$eta.y1*x1[i]+
                           x$eta.y2*x2[i]+
                           x$eta.y3*x3[i])/x$alpha.y0nd)*gamma(1+{1/x$alpha.y0nd})*pi[i])
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
                                        scale=exp(-(x$beta.y1d+
                                                      x$eta.y1*x1[i]+
                                                      x$eta.y2*x2[i]+
                                                      x$eta.y3*x3[i]+
                                                      x$delta*log(d[h]))/x$alpha.y1d),
                                        left=d[h])*(1-pi[i]))
      EY0 <- sum(EY0,(exp(-(x$beta.y0d+
                              x$delta*log(d[h])+
                              x$eta.y1*x1[i]+
                              x$eta.y2*x2[i]+
                              x$eta.y3*x3[i])/x$alpha.y0d)*gamma(1+{1/x$alpha.y0d}))*(1-pi[i]))
      
    }
    
    EY1.d <- EY1/pG
    EY0.d <- EY0/pG
    
    rm(EY0,EY1)
    
    ace.d[h] <- EY1.d-EY0.d
    
    ace.fd[h] <- ace.d[h]*dweib(h,
                                a=x$alpha.d,
                                b=x$beta.d+
                                  x$eta.d1*x1[i]+
                                  x$eta.d2*x2[i]+
                                  x$eta.d3*x3[i])
    
    rm(EY1.d,EY0.d)
  }
  
  pND<-sum(ID==1)/nrow(observed.data)
  pD<-1-pND
  
  aceoverall<-ace.nd*pND + sum(ace.fd)*pD
  
  list(ace.nd=ace.nd, ace.d=ace.d, ace.fd = ace.fd,
       aceoverall=aceoverall, 
       m1_d=m1_d, m2_d=m2_d, m3_d=m3_d, 
       m1_nd=m1_nd, m2_nd=m2_nd, m3_nd=m3_nd, 
       m1_ed=m1_ed, m2_ed=m2_ed, m3_ed=m3_ed,
       m1_ld=m1_ld, m2_ld=m2_ld, m3_ld=m3_ld)
}

dce_apply <- function(x){
  
  d<-seq(1,6)
  y<-seq(1,30)
  
  x1<-observed.data$x1st
  x2<-observed.data$x2
  x3<-observed.data$x3
  ID<-x$ID
  D<-x$D
  
  pi<-exp(x$eta.0+
            x$eta.1*x1 + 
            x$eta.2*x2 + 
            x$eta.3*x3)/(1+
                           exp(x$eta.0+
                                 x$eta.1*x1 + 
                                 x$eta.2*x2 + 
                                 x$eta.3*x3))
  
  
  dce.nd1<-dce.nd0<-NULL
  
  for(i in which(ID==1)){
    dce.nd1 <- rbind(dce.nd1,
                     Sweib(y, x$alpha.y1nd, 
                           x$beta.y1nd+
                             x$eta.y1*x1[i]+
                             x$eta.y2*x2[i]+
                             x$eta.y3*x3[i])*pi[i] ) 
    dce.nd0 <- rbind(dce.nd0,
                     Sweib(y, x$alpha.y0nd,
                           x$beta.y0nd+
                             x$eta.y1*x1[i]+
                             x$eta.y2*x2[i]+
                             x$eta.y3*x3[i])*pi[i])
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
      dce.d1 <- rbind(dce.d1,
                      Stweib(y,
                              sh=x$alpha.y1d,
                              sc=exp(-(x$beta.y1d+
                                         x$eta.y1*x1[i]+
                                         x$eta.y2*x2[i]+
                                         x$eta.y3*x3[i]+
                                         x$delta*log(d[h]))/x$alpha.y1d),
                              a=d[h])*(1-pi[i]))
      # dce.d1 <- rbind(dce.d1,
      #                 Stweib(y,
      #                        a=x$alpha.y1d,
      #                        b=x$beta.y1d+
      #                          x$eta.y1*x1[i]+
      #                          x$eta.y2*x2[i]+
      #                          x$eta.y3*x3[i]+
      #                          x$delta*log(d[h]),
      #                        l=d[h])*(1-pi[i]))
      
      dce.d0 <- rbind(dce.d0,
                      (Sweib(y,
                             x$alpha.y0d,
                             x$beta.y0d+
                               x$delta*log(d[h])+
                               x$eta.y1*x1[i]+
                               x$eta.y2*x2[i]+
                               x$eta.y3*x3[i]))*(1-pi[i]))
      
    }
    
    pG<-sum(1-pi[ID==0])
    
    dced0 <- colSums(dce.d0)/pG
    dced1 <- colSums(dce.d1)/pG
    
    dce.d <- rbind(dce.d,dced1-dced0)
    
    rm(dced0,dced1)
    
  }
  
  list(dce.nd=dce.nd, dce.d=dce.d)
  
}
