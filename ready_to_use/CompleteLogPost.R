logpost<-function(theta, theta.prior, complete.data){
  
  Z  <- complete.data$Z   #Treatment indicator
  ID <- complete.data$ID  #Discontinuation indicator: 
                          #ID=1 for never discontinued
  RD <- complete.data$RD  #RD=1 for treated units who discontinue
  D  <- complete.data$D   #Time to discontinuation for those who discontinue
  RY <- complete.data$RY  #RY=1 for units who die
  Y  <- complete.data$Y   #survival time

  etas_pos<-grepl("eta.",names(theta))+
    grepl("beta",names(theta))+
    grepl("eta.d",names(theta))+
    grepl("eta.y",names(theta))
  
  etas_d_pos<-grepl("beta.d",names(theta))+
    grepl("eta.d",names(theta))
  
  etas_y_pos<-grepl("beta.y",names(theta))+
    grepl("eta.y",names(theta))
  
  theta_eta<-NULL
  for(i in which(etas_pos==1)){
    theta_eta<-rbind(theta_eta,as.numeric(theta[i]))
  }
  theta_eta<-as.matrix(theta_eta)
  
  theta_deta<-NULL
  for(i in which(etas_d_pos==1)){
    theta_deta<-rbind(theta_deta,as.numeric(theta[i]))
  }
  theta_deta<-as.matrix(theta_deta)
  
  theta_yeta<-NULL
  for(i in which(etas_y_pos==1)){
    theta_yeta<-rbind(theta_yeta,as.numeric(theta[i]))
  }
  theta_yeta<-as.matrix(theta_yeta)
  
  loglik<-0
  #Z=1 ID=1 RY=1 ---> Never discontinued who die under treatment
  ii <- Z*ID*RY==1
  loglik <- loglik + 
    sum({log(theta$pi[ii]) + 
        log(dweib(Y[ii], 
                  theta$alpha.y1nd, 
                  theta$beta.y1nd+
                    x[ii,-1]%*%theta_yeta))},
        na.rm = TRUE)
  
  #Z=1 ID=1 RY=0 ---> Never discontinued who do not die under treatment
  ii <- Z*ID*(1-RY)==1
  loglik <- loglik + 
    sum({log(theta$pi[ii]) + 
        log(Sweib(Y[ii], 
                  theta$alpha.y1nd, 
                  theta$beta.y1nd+
                    x[ii,-1]%*%theta_yeta))},
        na.rm = TRUE)
  
  #Z=1 ID=0 RD=1 RY=1 ---> Discontinued who die under treatment for whom we 
  #observe the discontinuation
  ii <- Z*(1-ID)*RD*RY==1 
  loglik <- loglik + 
    sum({log(1-theta$pi[ii]) + 
        log(dweib(D[ii], 
                  theta$alpha.d, 
                  theta$beta.d+
                    x[ii,-1]%*%theta_deta))+
        # log(dtweibull(Y[ii],
        #           shape=theta$alpha.y1d,
        #           scale=exp(-(theta$beta.y1d+
        #                   theta$eta.y1*x1[ii]+
        #                   theta$eta.y2*x2[ii]+
        #                   theta$eta.y3*x3[ii]+
        #                     theta$delta*log(D[ii]))/theta$alpha.y1d),
        #           a=D[ii]))},
        log(dtweib(Y[ii],
                   sh=theta$alpha.y1d,
                   sc=exp(-(theta$beta.y1d+
                              x[ii,-1]%*%theta_yeta+
                              theta$delta*log(D[ii]))/theta$alpha.y1d),
                   a=D[ii]))},
        na.rm = TRUE)
  
  #Z=1 ID=0 RD=1 RY=0 ---> Discontinued who do not die under treatment for whom 
  #we observe the discontinuation
  ii <- Z*(1-ID)*RD*(1-RY)==1
  loglik <- loglik + 
    sum({log(1-theta$pi[ii]) + 
        log(dweib(D[ii], 
                  theta$alpha.d, 
                  theta$beta.d+x[ii,-1]%*%theta_deta))+
        # log(Stweib(Y[ii],
        #           sh=theta$alpha.y1d,
        #           sc=exp(-(theta$beta.y1d+
        #                               theta$eta.y1*x1[ii]+
        #                               theta$eta.y2*x2[ii]+
        #                               theta$eta.y3*x3[ii]+
        #                               theta$delta*log(D[ii]))/theta$alpha.y1d),
        #           a=D[ii]))
        log(Stweib(Y[ii],
                   a=theta$alpha.y1d,
                   b=theta$beta.y1d+x[ii,-1]%*%theta_yeta+
                     theta$delta*log(D[ii]),
                   l=D[ii]))
      },
        na.rm = TRUE)
  
  #Z=1 ID=0 RD=0 RY=0 ---> Discontinued who discontinue after the end of the 
  #study and do not die under treatment
  ii <- Z*(1-ID)*(1-RD)*(1-RY)==1
  loglik <- loglik + 
    sum({log(1-theta$pi[ii]) + 
        log(Sweib(D[ii], 
                  theta$alpha.d, 
                  theta$beta.d+x[ii,-1]%*%theta_deta))},
        na.rm = TRUE)
  
  
  #Z=0 ID=1 RY=1 ---> Never discontinued who die under control
  ii <- (1-Z)*ID*RY==1
  loglik <- loglik + 
    sum({log(theta$pi[ii]) + 
        log(dweib(Y[ii], 
                    theta$alpha.y0nd, 
                    theta$beta.y0nd+x[ii,-1]%*%theta_yeta))
      },na.rm = TRUE)

  #Z=0 ID=1 RY=0 ---> Never discontinued who do not die under control
  ii <- (1-Z)*ID*(1-RY)==1
  loglik <- loglik + 
    sum({log(theta$pi[ii]) + 
        log(Sweib(Y[ii], 
                    theta$alpha.y0nd, 
                    theta$beta.y0nd+x[ii,-1]%*%theta_yeta))},
        na.rm = TRUE)
  
  #Z=0 ID=0 RY=1 ---> Discontinued who die under control
  ii <- (1-Z)*(1-ID)*RY==1
  loglik <- loglik + 
    sum({log(1-theta$pi[ii]) + 
        log(dweib(D[ii],
                  theta$alpha.d,
                  theta$beta.d+x[ii,-1]%*%theta_deta))+
        log(dweib(Y[ii], 
                    theta$alpha.y0d, 
                    {theta$beta.y0d+
                        theta$delta*log(D[ii])+x[ii,-1]%*%theta_yeta}))
      },na.rm = TRUE)
  
  #Z=0 ID=0 RY=0 ---> Discontinued who do not die under control
  ii <- (1-Z)*(1-ID)*(1-RY)==1
  loglik <- loglik + 
    sum({log(1-theta$pi[ii]) + 
        log(dweib(D[ii],
                  theta$alpha.d,
                  theta$beta.d+x[ii,-1]%*%theta_deta))+
        log(Sweib(Y[ii], 
                    theta$alpha.y0d, 
                    {theta$beta.y0d+
                        theta$delta*log(D[ii])+x[ii,-1]%*%theta_yeta}))},
        na.rm = TRUE)
  
  logprior<-0
  
  sigmaprior_pos<-grepl("sigma2.",names(theta.prior))+
    grepl("sigma2.d",names(theta.prior))+
    grepl("sigma2.y",names(theta.prior))+
    grepl("sigma2.delta",names(theta.prior))
  sigmaprior<-NULL
  for(i in which(sigmaprior_pos==1)){
    sigmaprior<-rbind(sigmaprior,as.numeric(theta.prior[i]))
  }
  sigmaprior<-as.matrix(sigmaprior)
  sigmaprior <- diag(c(sigmaprior)) 
  
  muprior_pos<-grepl("mu.",names(theta.prior))+
    grepl("mu.d",names(theta.prior))+
    grepl("mu.y",names(theta.prior))+
    grepl("mu.delta",names(theta.prior))
  muprior<-NULL
  for(i in which(muprior_pos==1)){
    muprior<-rbind(muprior,as.numeric(theta.prior[i]))
  }

  logprior <- logprior + 
    dmvnorm(c(theta_eta), 
            c(muprior), 
            sigmaprior, log=TRUE)
  
  sigmaprior.d_pos<-grepl("sigma2.d",names(theta.prior))
  sigmaprior.d<-NULL
  for(i in which(sigmaprior.d_pos==1)){
    sigmaprior.d<-rbind(sigmaprior.d,as.numeric(theta.prior[i]))
  }
  sigmaprior.d<-as.matrix(sigmaprior.d)
  sigmaprior.d <- diag(c(sigmaprior.d)) 
  
  mupriord_pos<-grepl("mu.d",names(theta.prior))
  mupriord<-NULL
  for(i in which(mupriord_pos==1)){
    mupriord<-rbind(mupriord,as.numeric(theta.prior[i]))
  }
  
  logprior <- logprior + 
    dgamma(theta$alpha.d, 
           shape=theta.prior$a.d, 
           scale=theta.prior$b.d, 
           log=TRUE) +
    dmvnorm(c(theta$beta.d, theta_deta), 
            c(mupriord), 
            sigmaprior.d, log=TRUE)
  
  logprior <- logprior + 
    dgamma(theta$alpha.y1nd, 
           shape=theta.prior$a.y1nd, 
           scale=theta.prior$b.y1nd, 
           log=TRUE) +
    dnorm(theta$beta.y1nd, 
          theta.prior$mu.y1nd, 
          sqrt(theta.prior$sigma2.y1nd), 
          log=TRUE)
  
  logprior <- logprior + 
    dgamma(theta$alpha.y1d, 
           shape=theta.prior$a.y1d, 
           scale=theta.prior$b.y1d, 
           log=TRUE) +
    dnorm(theta$beta.y1d, 
          theta.prior$mu.y1d, 
          sqrt(theta.prior$sigma2.y1d), 
          log=TRUE)
  
  logprior <- logprior + 
    dgamma(theta$alpha.y0nd, 
           shape=theta.prior$a.y0nd, 
           scale=theta.prior$b.y0nd, 
           log=TRUE) +
    dnorm(theta$beta.y0nd, 
          theta.prior$mu.y0nd, 
          sqrt(theta.prior$sigma2.y0nd), 
          log=TRUE)
  
  logprior <- logprior + 
    dgamma(theta$alpha.y0d, 
           shape=theta.prior$a.y0d, 
           scale=theta.prior$b.y0d, 
           log=TRUE) +
    dnorm(theta$beta.y0d, 
          theta.prior$mu.y0d, 
          sqrt(theta.prior$sigma2.y0d), 
          log=TRUE)
  
  sigmaprior.y_pos<-grepl("sigma2.y",names(theta.prior))+
    grepl("sigma2.y1nd",names(theta.prior))+
    grepl("sigma2.y0nd",names(theta.prior))+
    grepl("sigma2.y1d",names(theta.prior))+
    grepl("sigma2.y0d",names(theta.prior))
  
  sigmaprior.y<-NULL
  for(i in which(sigmaprior.y_pos==1)){
    sigmaprior.y<-rbind(sigmaprior.y,as.numeric(theta.prior[i]))
  }
  sigmaprior.y<-as.matrix(sigmaprior.y)
  sigmaprior.y <- diag(c(sigmaprior.y)) 
  
  mupriory_pos<-grepl("mu.y",names(theta.prior))+
    grepl("mu.y1nd",names(theta.prior))+
    grepl("mu.y0nd",names(theta.prior))+
    grepl("mu.y1d",names(theta.prior))+
    grepl("mu.y0d",names(theta.prior))
  
  mupriory<-NULL
  for(i in which(mupriory_pos==1)){
    mupriory<-rbind(mupriory,as.numeric(theta.prior[i]))
  }
  
  logprior <- logprior + 
    dmvnorm(c(theta_yeta), 
            c(mupriory), 
            sigmaprior.y, log=TRUE)
  
  logprior <- logprior +
    dnorm(theta$delta,
          theta.prior$mu.delta,
          sqrt(theta.prior$sigma2.delta),
          log=TRUE)
  
  logprior + loglik
}


logpost.i <-function(theta, theta.prior, complete.data){
  
  n  <- nrow(complete.data)
  Z  <- complete.data$Z
  ID  <- complete.data$ID
  RD <- complete.data$RD
  D  <- complete.data$D
  RY <- complete.data$RY
  Y  <- complete.data$Y

  etas_pos<-grepl("eta.",names(theta))+
    grepl("beta",names(theta))+
    grepl("eta.d",names(theta))+
    grepl("eta.y",names(theta))
  
  etas_d_pos<-grepl("beta.d",names(theta))+
    grepl("eta.d",names(theta))
  
  etas_y_pos<-grepl("beta.y",names(theta))+
    grepl("eta.y",names(theta))
  
  theta_eta<-NULL
  for(i in which(etas_pos==1)){
    theta_eta<-rbind(theta_eta,as.numeric(theta[i]))
  }
  theta_eta<-as.matrix(theta_eta)
  
  theta_deta<-NULL
  for(i in which(etas_d_pos==1)){
    theta_deta<-rbind(theta_deta,as.numeric(theta[i]))
  }
  theta_deta<-as.matrix(theta_deta)
  
  theta_yeta<-NULL
  for(i in which(etas_y_pos==1)){
    theta_yeta<-rbind(theta_yeta,as.numeric(theta[i]))
  }
  theta_yeta<-as.matrix(theta_yeta)
  
  loglik<-rep(0, n)
  #Z=1 ID=1 RY=1 ---> Never discontinued who die under treatment
  ii<-Z*ID*RY==1 
  loglik[ii] <- log(theta$pi[ii]) + 
    log(dweib(Y[ii], 
              theta$alpha.y1nd, 
              theta$beta.y1nd+x[ii,-1]%*%theta_yeta))
  
  #Z=1 ID=1 RY=0 ---> Never discontinued who do not die under treatment
  ii<-Z*ID*(1-RY)==1
  loglik[ii] <- log(theta$pi[ii]) + 
    log(Sweib(Y[ii], 
              theta$alpha.y1nd, 
              theta$beta.y1nd+x[ii,-1]%*%theta_yeta))
  
  #Z=1 ID=0 RD=1 RY=1 ---> Discontinued who die under treatment for whom we 
  #observe the discontinuation
  ii <- Z*(1-ID)*RD*RY==1
  loglik[ii] <- log(1-theta$pi[ii]) + 
    log(dweib(D[ii], 
              theta$alpha.d, 
              theta$beta.d+x[ii,-1]%*%theta_deta))+
    # log(dtweibull(Y[ii],
    #           shape=theta$alpha.y1d,
    #           scale=exp(-(theta$beta.y1d+
    #                               theta$eta.y1*x1[ii]+
    #                               theta$eta.y2*x2[ii]+
    #                               theta$eta.y3*x3[ii]+
    #                               theta$delta*log(D[ii]))/theta$alpha.y1d),
    #           a=D[ii]))
    log(dtweib(Y[ii],
               sh=theta$alpha.y1d,
               sc=exp(-(theta$beta.y1d+x[ii,-1]%*%theta_yeta+
                          theta$delta*log(D[ii]))/theta$alpha.y1d),
               a=D[ii]))
  
  
  #Z=1 ID=0 RD=1 RY=0 ---> Discontinued who do not die under treatment for whom 
  #we observe the discontinuation
  ii<-Z*(1-ID)*RD*(1-RY)==1
  loglik[ii] <- log(1-theta$pi[ii]) + 
    log(dweib(D[ii], 
              theta$alpha.d, 
              theta$beta.d+x[ii,-1]%*%theta_deta))+
    # log(Stweib(Y[ii],
    #           sh=theta$alpha.y1d,
    #           sc=exp(-(theta$beta.y1d+
    #                         theta$eta.y1*x1[ii]+
    #                         theta$eta.y2*x2[ii]+
    #                         theta$eta.y3*x3[ii]+
    #                         theta$delta*log(D[ii]))/theta$alpha.y1d),
    #           a=D[ii]))
    log(Stweib(Y[ii],
               a=theta$alpha.y1d,
               b=theta$beta.y1d+x[ii,-1]%*%theta_yeta+
                 theta$delta*log(D[ii]),
               l=D[ii]))
  
  
  #Z=1 ID=0 RD=0 RY=0 ---> Discontinued who discontinue after the end of the 
  #study and do not die under treatment
  ii<-Z*(1-ID)*(1-RD)*(1-RY)==1
  loglik[ii] <- log(1-theta$pi[ii]) + 
    log(Sweib(D[ii], 
              theta$alpha.d, 
              theta$beta.d+x[ii,-1]%*%theta_deta))
  
  
  #Z=0 ID=1 RY=1 ---> Never discontinued who die under control
  ii <- (1-Z)*ID*RY==1
  loglik[ii] <- log(theta$pi[ii]) + 
    log(dweib(Y[ii], 
                theta$alpha.y0nd, 
                theta$beta.y0nd+x[ii,-1]%*%theta_yeta))
  
  #Z=0 ID=1 RY=0 Y>Y1*k ---> Never discontinued who do not die under control
  ii <- (1-Z)*ID*(1-RY)==1
  loglik[ii] <- log(theta$pi[ii]) + 
    log(Sweib(Y[ii], 
                theta$alpha.y0nd, 
                theta$beta.y0nd+x[ii,-1]%*%theta_yeta))
  
  #Z=0 ID=0 RY=1 ---> Discontinued who die under control
  ii <- (1-Z)*(1-ID)*RY==1
  loglik[ii] <- log(1-theta$pi[ii]) + 
    log(dweib(D[ii],
              theta$alpha.d,
              theta$beta.d+x[ii,-1]%*%theta_deta))+
    log(dweib(Y[ii], 
                theta$alpha.y0d, 
                {theta$beta.y0d+
                    theta$delta*log(D[ii])+x[ii,-1]%*%theta_yeta}))
  
  #Z=0 ID=0  RY=0 ---> Discontinued who do not die under control
  ii <- (1-Z)*(1-ID)*(1-RY)==1
  loglik[ii] <- log(1-theta$pi[ii]) + 
    log(dweib(D[ii],
              theta$alpha.d,
              theta$beta.d+x[ii,-1]%*%theta_deta))+
    log(Sweib(Y[ii], 
                theta$alpha.y0d, 
                {theta$beta.y0d+
                    theta$delta*log(D[ii])+x[ii,-1]%*%theta_yeta}))
  
  rm(ii)
  
  logprior<-0
  
  sigmaprior_pos<-grepl("sigma2.",names(theta.prior))+
    grepl("sigma2.d",names(theta.prior))+
    grepl("sigma2.y",names(theta.prior))+
    grepl("sigma2.delta",names(theta.prior))
  sigmaprior<-NULL
  for(i in which(sigmaprior_pos==1)){
    sigmaprior<-rbind(sigmaprior,as.numeric(theta.prior[i]))
  }
  sigmaprior<-as.matrix(sigmaprior)
  sigmaprior <- diag(c(sigmaprior)) 
  
  muprior_pos<-grepl("mu.",names(theta.prior))+
    grepl("mu.d",names(theta.prior))+
    grepl("mu.y",names(theta.prior))+
    grepl("mu.delta",names(theta.prior))
  muprior<-NULL
  for(i in which(muprior_pos==1)){
    muprior<-rbind(muprior,as.numeric(theta.prior[i]))
  }
  
  logprior <- logprior + 
    dmvnorm(c(theta_eta), 
            c(muprior), 
            sigmaprior, log=TRUE)
  
  sigmaprior.d_pos<-grepl("sigma2.d",names(theta.prior))
  sigmaprior.d<-NULL
  for(i in which(sigmaprior.d_pos==1)){
    sigmaprior.d<-rbind(sigmaprior.d,as.numeric(theta.prior[i]))
  }
  sigmaprior.d<-as.matrix(sigmaprior.d)
  sigmaprior.d <- diag(c(sigmaprior.d)) 
  
  mupriord_pos<-grepl("mu.d",names(theta.prior))
  mupriord<-NULL
  for(i in which(mupriord_pos==1)){
    mupriord<-rbind(mupriord,as.numeric(theta.prior[i]))
  }
  
  logprior <- logprior + 
    dgamma(theta$alpha.d, 
           shape=theta.prior$a.d, 
           scale=theta.prior$b.d, 
           log=TRUE) +
    dmvnorm(c(theta$beta.d, theta_deta), 
            c(mupriord), 
            sigmaprior.d, log=TRUE)
  
  logprior <- logprior + 
    dgamma(theta$alpha.y1nd, 
           shape=theta.prior$a.y1nd, 
           scale=theta.prior$b.y1nd, 
           log=TRUE) +
    dnorm(theta$beta.y1nd, 
          theta.prior$mu.y1nd, 
          sqrt(theta.prior$sigma2.y1nd), 
          log=TRUE)
  
  logprior <- logprior + 
    dgamma(theta$alpha.y1d, 
           shape=theta.prior$a.y1d, 
           scale=theta.prior$b.y1d, 
           log=TRUE) +
    dnorm(theta$beta.y1d, 
          theta.prior$mu.y1d, 
          sqrt(theta.prior$sigma2.y1d), 
          log=TRUE)
  
  logprior <- logprior + 
    dgamma(theta$alpha.y0nd, 
           shape=theta.prior$a.y0nd, 
           scale=theta.prior$b.y0nd, 
           log=TRUE) +
    dnorm(theta$beta.y0nd, 
          theta.prior$mu.y0nd, 
          sqrt(theta.prior$sigma2.y0nd), 
          log=TRUE)
  
  logprior <- logprior + 
    dgamma(theta$alpha.y0d, 
           shape=theta.prior$a.y0d, 
           scale=theta.prior$b.y0d, 
           log=TRUE) +
    dnorm(theta$beta.y0d, 
          theta.prior$mu.y0d, 
          sqrt(theta.prior$sigma2.y0d), 
          log=TRUE)
  
  sigmaprior.y_pos<-grepl("sigma2.y",names(theta.prior))+
    grepl("sigma2.y1nd",names(theta.prior))+
    grepl("sigma2.y0nd",names(theta.prior))+
    grepl("sigma2.y1d",names(theta.prior))+
    grepl("sigma2.y0d",names(theta.prior))
  
  sigmaprior.y<-NULL
  for(i in which(sigmaprior.y_pos==1)){
    sigmaprior.y<-rbind(sigmaprior.y,as.numeric(theta.prior[i]))
  }
  sigmaprior.y<-as.matrix(sigmaprior.y)
  sigmaprior.y <- diag(c(sigmaprior.y)) 
  
  mupriory_pos<-grepl("mu.y",names(theta.prior))+
    grepl("mu.y1nd",names(theta.prior))+
    grepl("mu.y0nd",names(theta.prior))+
    grepl("mu.y1d",names(theta.prior))+
    grepl("mu.y0d",names(theta.prior))
  
  mupriory<-NULL
  for(i in which(mupriory_pos==1)){
    mupriory<-rbind(mupriory,as.numeric(theta.prior[i]))
  }
  
  logprior <- logprior + 
    dmvnorm(c(theta_yeta), 
            c(mupriory), 
            sigmaprior.y, log=TRUE)
  
  logprior <- logprior +
    dnorm(theta$delta,
          theta.prior$mu.delta,
          sqrt(theta.prior$sigma2.delta),
          log=TRUE)
  
  logprior + loglik
}
