logpost_expweib<-function(theta, theta.prior, complete.data){
  
  Z  <- complete.data$Z   #Treatment indicator
  ID <- complete.data$ID  #Discontinuation indicator: ID=1 for never discontinued
  RD <- complete.data$RD  #RD=1 for treated units who discontinue
  D  <- complete.data$D   #Time to discontinuation for those who discontinue
  RY <- complete.data$RY  #RY=1 for units who die
  Y  <- complete.data$Y   #survival time
  x1 <- complete.data$x1
  x2 <- complete.data$x2
  x3 <- complete.data$x3
  
  
  loglik<-0
  #Z=1 ID=1 RY=1 ---> Never discontinued who die under treatment
  ii <- Z*ID*RY==1
  loglik <- loglik + 
    sum({log(theta$pi[ii]) + 
        log(dweib(Y[ii], 
                  theta$alpha.y1nd, 
                  theta$beta.y1nd+
                    theta$eta.y1*x1[ii]+
                    theta$eta.y2*x2[ii]+
                    theta$eta.y3*x3[ii]))},
        na.rm = TRUE)
  
  #Z=1 ID=1 RY=0 ---> Never discontinued who do not die under treatment
  ii <- Z*ID*(1-RY)==1
  loglik <- loglik + 
    sum({log(theta$pi[ii]) + 
        log(Sweib(Y[ii], 
                  theta$alpha.y1nd, 
                  theta$beta.y1nd+
                    theta$eta.y1*x1[ii]+
                    theta$eta.y2*x2[ii]+
                    theta$eta.y3*x3[ii]))},
        na.rm = TRUE)
  
  #Z=1 ID=0 RD=1 RY=1 ---> Discontinued who die under treatment for whom we 
  #observe the discontinuation
  ii <- Z*(1-ID)*RD*RY==1 
  loglik <- loglik + 
    sum({log(1-theta$pi[ii]) + 
        log(dexp(D[ii], 
                  exp(theta$beta.d+
                    theta$eta.d1*x1[ii]+
                    theta$eta.d2*x2[ii]+
                    theta$eta.d3*x3[ii])))+
        log(dtweib(Y[ii],
                   sh=theta$alpha.y1d,
                   sc=exp(-(theta$beta.y1d+
                              theta$eta.y1*x1[ii]+
                              theta$eta.y2*x2[ii]+
                              theta$eta.y3*x3[ii]+
                              theta$delta*log(D[ii]))/theta$alpha.y1d),
                   a=D[ii]))},
        na.rm = TRUE)
  
  #Z=1 ID=0 RD=1 RY=0 ---> Discontinued who do not die under treatment for whom 
  #we observe the discontinuation
  ii <- Z*(1-ID)*RD*(1-RY)==1
  loglik <- loglik + 
    sum({log(1-theta$pi[ii]) + 
        log(dexp(D[ii], 
                  exp(theta$beta.d+
                    theta$eta.d1*x1[ii]+
                    theta$eta.d2*x2[ii]+
                    theta$eta.d3*x3[ii])))+
        log(Stweib(Y[ii],
                   a=theta$alpha.y1d,
                   b=theta$beta.y1d+
                     theta$eta.y1*x1[ii]+
                     theta$eta.y2*x2[ii]+
                     theta$eta.y3*x3[ii]+
                     theta$delta*log(D[ii]),
                   l=D[ii]))
    },
    na.rm = TRUE)
  
  #Z=1 ID=0 RD=0 RY=0 ---> Discontinued who discontinue after the end of the 
  #study and do not die under treatment
  ii <- Z*(1-ID)*(1-RD)*(1-RY)==1
  loglik <- loglik + 
    sum({log(1-theta$pi[ii]) + 
        log(Sexp(D[ii], 
                  exp(theta$beta.d+
                    theta$eta.d1*x1[ii]+
                    theta$eta.d2*x2[ii]+
                    theta$eta.d3*x3[ii])))},
        na.rm = TRUE)
  
  
  #Z=0 ID=1 RY=1 ---> Never discontinued who die under control
  ii <- (1-Z)*ID*RY==1
  loglik <- loglik + 
    sum({log(theta$pi[ii]) + 
        log(dweib(Y[ii], 
                  theta$alpha.y0nd, 
                  theta$beta.y0nd+
                    theta$eta.y1*x1[ii]+
                    theta$eta.y2*x2[ii]+
                    theta$eta.y3*x3[ii]))
    },na.rm = TRUE)
  
  #Z=0 ID=1 RY=0 ---> Never discontinued who do not die under control
  ii <- (1-Z)*ID*(1-RY)==1
  loglik <- loglik + 
    sum({log(theta$pi[ii]) + 
        log(Sweib(Y[ii], 
                  theta$alpha.y0nd, 
                  theta$beta.y0nd+
                    theta$eta.y1*x1[ii]+
                    theta$eta.y2*x2[ii]+
                    theta$eta.y3*x3[ii]))},
        na.rm = TRUE)
  
  #Z=0 ID=0 RY=1 ---> Discontinued who die under control
  ii <- (1-Z)*(1-ID)*RY==1
  loglik <- loglik + 
    sum({log(1-theta$pi[ii]) + 
        log(dexp(D[ii],
                  exp(theta$beta.d+
                    theta$eta.d1*x1[ii]+
                    theta$eta.d2*x2[ii]+
                    theta$eta.d3*x3[ii])))+
        log(dweib(Y[ii], 
                  theta$alpha.y0d, 
                  {theta$beta.y0d+
                      theta$delta*log(D[ii])+
                      theta$eta.y1*x1[ii]+
                      theta$eta.y2*x2[ii]+
                      theta$eta.y3*x3[ii]}))
    },na.rm = TRUE)
  
  #Z=0 ID=0 RY=0 ---> Discontinued who do not die under control
  ii <- (1-Z)*(1-ID)*(1-RY)==1
  loglik <- loglik + 
    sum({log(1-theta$pi[ii]) + 
        log(dexp(D[ii],
                  exp(theta$beta.d+
                    theta$eta.d1*x1[ii]+
                    theta$eta.d2*x2[ii]+
                    theta$eta.d3*x3[ii])))+
        log(Sweib(Y[ii], 
                  theta$alpha.y0d, 
                  {theta$beta.y0d+
                      theta$delta*log(D[ii])+
                      theta$eta.y1*x1[ii]+
                      theta$eta.y2*x2[ii]+
                      theta$eta.y3*x3[ii]}))},
        na.rm = TRUE)
  
  logprior<-0
  
  sigmaprior <- diag(c(theta.prior$sigma2.0,
                       theta.prior$sigma2.1,
                       theta.prior$sigma2.2,
                       theta.prior$sigma2.3))  
  
  logprior <- logprior + 
    dmvnorm(c(theta$eta.0,
              theta$eta.1,
              theta$eta.2,
              theta$eta.3), 
            c(theta.prior$mu.0,
              theta.prior$mu.1,
              theta.prior$mu.2,
              theta.prior$mu.3), 
            sigmaprior, log=TRUE)
  
  sigmaprior.d <- diag(c(theta.prior$sigma2.d1,
                         theta.prior$sigma2.d2,
                         theta.prior$sigma2.d3))  
  
  logprior <- logprior + 
    dnorm(theta$beta.d, 
          theta.prior$mu.d, 
          sqrt(theta.prior$sigma2.d), 
          log=TRUE) + 
    dmvnorm(c(theta$eta.d1,
              theta$eta.d2,
              theta$eta.d3), 
            c(theta.prior$mu.d1,
              theta.prior$mu.d2,
              theta.prior$mu.d3), 
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
  
  sigmaprior.y <- diag(c(theta.prior$sigma2.y1,
                         theta.prior$sigma2.y2,
                         theta.prior$sigma2.y3))  
  
  logprior <- logprior + 
    dmvnorm(c(theta$eta.y1,
              theta$eta.y2,
              theta$eta.y3), 
            c(theta.prior$mu.y1,
              theta.prior$mu.y2,
              theta.prior$mu.y3), 
            sigmaprior.y, log=TRUE)
  
  logprior <- logprior +
    dnorm(theta$delta,
          theta.prior$mu.delta,
          sqrt(theta.prior$sigma2.delta),
          log=TRUE)
  
  logprior + loglik
}


logpost.i_expweib <-function(theta, theta.prior, complete.data){
  
  n  <- nrow(complete.data)
  Z  <- complete.data$Z
  ID  <- complete.data$ID
  RD <- complete.data$RD
  D  <- complete.data$D
  RY <- complete.data$RY
  Y  <- complete.data$Y
  x1 <- complete.data$x1
  x2 <- complete.data$x2
  x3 <- complete.data$x3
  
  loglik<-rep(0, n)
  #Z=1 ID=1 RY=1 ---> Never discontinued who die under treatment
  ii<-Z*ID*RY==1 
  loglik[ii] <- log(theta$pi[ii]) + 
    log(dweib(Y[ii], 
              theta$alpha.y1nd, 
              theta$beta.y1nd+
                theta$eta.y1*x1[ii]+
                theta$eta.y2*x2[ii]+
                theta$eta.y3*x3[ii]))
  
  #Z=1 ID=1 RY=0 ---> Never discontinued who do not die under treatment
  ii<-Z*ID*(1-RY)==1
  loglik[ii] <- log(theta$pi[ii]) + 
    log(Sweib(Y[ii], 
              theta$alpha.y1nd, 
              theta$beta.y1nd+
                theta$eta.y1*x1[ii]+
                theta$eta.y2*x2[ii]+
                theta$eta.y3*x3[ii]))
  
  #Z=1 ID=0 RD=1 RY=1 ---> Discontinued who die under treatment for whom we 
  #observe the discontinuation
  ii <- Z*(1-ID)*RD*RY==1
  loglik[ii] <- log(1-theta$pi[ii]) + 
    log(dexp(D[ii], 
              exp(theta$beta.d+
                theta$eta.d1*x1[ii]+
                theta$eta.d2*x2[ii]+
                theta$eta.d3*x3[ii])))+
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
               sc=exp(-(theta$beta.y1d+
                          theta$eta.y1*x1[ii]+
                          theta$eta.y2*x2[ii]+
                          theta$eta.y3*x3[ii]+
                          theta$delta*log(D[ii]))/theta$alpha.y1d),
               a=D[ii]))
  
  
  #Z=1 ID=0 RD=1 RY=0 ---> Discontinued who do not die under treatment for whom 
  #we observe the discontinuation
  ii<-Z*(1-ID)*RD*(1-RY)==1
  loglik[ii] <- log(1-theta$pi[ii]) + 
    log(dexp(D[ii], 
              exp(theta$beta.d+
                theta$eta.d1*x1[ii]+
                theta$eta.d2*x2[ii]+
                theta$eta.d3*x3[ii])))+
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
               b=theta$beta.y1d+
                 theta$eta.y1*x1[ii]+
                 theta$eta.y2*x2[ii]+
                 theta$eta.y3*x3[ii]+
                 theta$delta*log(D[ii]),
               l=D[ii]))
  
  
  #Z=1 ID=0 RD=0 RY=0 ---> Discontinued who discontinue after the end of the 
  #study and do not die under treatment
  ii<-Z*(1-ID)*(1-RD)*(1-RY)==1
  loglik[ii] <- log(1-theta$pi[ii]) + 
    log(Sexp(D[ii], 
              exp(theta$beta.d+
                theta$eta.d1*x1[ii]+
                theta$eta.d2*x2[ii]+
                theta$eta.d3*x3[ii])))
  
  
  #Z=0 ID=1 RY=1 ---> Never discontinued who die under control
  ii <- (1-Z)*ID*RY==1
  loglik[ii] <- log(theta$pi[ii]) + 
    log(dweib(Y[ii], 
              theta$alpha.y0nd, 
              theta$beta.y0nd+
                theta$eta.y1*x1[ii]+
                theta$eta.y2*x2[ii]+
                theta$eta.y3*x3[ii]))
  
  #Z=0 ID=1 RY=0 Y>Y1*k ---> Never discontinued who do not die under control
  ii <- (1-Z)*ID*(1-RY)==1
  loglik[ii] <- log(theta$pi[ii]) + 
    log(Sweib(Y[ii], 
              theta$alpha.y0nd, 
              theta$beta.y0nd+
                theta$eta.y1*x1[ii]+
                theta$eta.y2*x2[ii]+
                theta$eta.y3*x3[ii]))
  
  #Z=0 ID=0 RY=1 ---> Discontinued who die under control
  ii <- (1-Z)*(1-ID)*RY==1
  loglik[ii] <- log(1-theta$pi[ii]) + 
    log(dexp(D[ii],
              exp(theta$beta.d+
                theta$eta.d1*x1[ii]+
                theta$eta.d2*x2[ii]+
                theta$eta.d3*x3[ii])))+
    log(dweib(Y[ii], 
              theta$alpha.y0d, 
              {theta$beta.y0d+
                  theta$delta*log(D[ii])+
                  theta$eta.y1*x1[ii]+
                  theta$eta.y2*x2[ii]+
                  theta$eta.y3*x3[ii]}))
  
  #Z=0 ID=0  RY=0 ---> Discontinued who do not die under control
  ii <- (1-Z)*(1-ID)*(1-RY)==1
  loglik[ii] <- log(1-theta$pi[ii]) + 
    log(dexp(D[ii],
              exp(theta$beta.d+
                theta$eta.d1*x1[ii]+
                theta$eta.d2*x2[ii]+
                theta$eta.d3*x3[ii])))+
    log(Sweib(Y[ii], 
              theta$alpha.y0d, 
              {theta$beta.y0d+
                  theta$delta*log(D[ii])+
                  theta$eta.y1*x1[ii]+
                  theta$eta.y2*x2[ii]+
                  theta$eta.y3*x3[ii]}))
  
  rm(ii)
  
  logprior <- 0
  
  sigmaprior <- diag(c(theta.prior$sigma2.0,
                       theta.prior$sigma2.1,
                       theta.prior$sigma2.2,
                       theta.prior$sigma2.3))  
  
  logprior <- logprior + 
    dmvnorm(c(theta$eta.0,
              theta$eta.1,
              theta$eta.2,
              theta$eta.3), 
            c(theta.prior$mu.0,
              theta.prior$mu.1,
              theta.prior$mu.2,
              theta.prior$mu.3), 
            sigmaprior, log=TRUE)
  
  sigmaprior.d <- diag(c(theta.prior$sigma2.d1,
                         theta.prior$sigma2.d2,
                         theta.prior$sigma2.d3))  
  
  logprior <- logprior + 
    dnorm(theta$beta.d, 
          theta.prior$mu.d, 
          sqrt(theta.prior$sigma2.d), 
          log=TRUE) + 
    dmvnorm(c(theta$eta.d1,
              theta$eta.d2,
              theta$eta.d3), 
            c(theta.prior$mu.d1,
              theta.prior$mu.d2,
              theta.prior$mu.d3), 
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
           scale=theta.prior$b.y0nd, log=TRUE) +
    dnorm(theta$beta.y0nd, 
          theta.prior$mu.y0nd, 
          sqrt(theta.prior$sigma2.y0nd), log=TRUE)
  
  logprior <- logprior + 
    dgamma(theta$alpha.y0d, 
           shape=theta.prior$a.y0d, 
           scale=theta.prior$b.y0d, log=TRUE) +
    dnorm(theta$beta.y0d, 
          theta.prior$mu.y0d, 
          sqrt(theta.prior$sigma2.y0d), log=TRUE)
  
  sigmaprior.y <- diag(c(theta.prior$sigma2.y1,
                         theta.prior$sigma2.y2,
                         theta.prior$sigma2.y3))  
  
  logprior <- logprior + 
    dmvnorm(c(theta$eta.y1,theta$eta.y2,theta$eta.y3), 
            c(theta.prior$mu.y1,theta.prior$mu.y2,theta.prior$mu.y3), 
            sigmaprior.y, log=TRUE)
  
  logprior <- logprior +
    dnorm(theta$delta,
          theta.prior$mu.delta,
          sqrt(theta.prior$sigma2.delta), log=TRUE)
  
  logprior + loglik
}
