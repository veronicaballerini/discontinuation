rmste_apply_expweib <- function(x){
  
  d<-seq(1,5) 
  y<-seq(1,20)
  n_obs <- nrow(observed.data)
  
  # Build covariate grid
  cov <- create_covariate_grid_general(
    continuous_vars = list(x1 = c(-2.9, 2.9, 0.1)),
    binary_count = 2
  )
  
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
  
  pG<-sum(pi)
  
  ###################
  ### DCE for ND ####
  ###################
  
  linpred_y <- xs[, -1, drop = FALSE] %*% etas_y
  
  SY1_ND <- sapply(y, function(h) {
    surv_vals <- Sweib(h, x$alpha.y1nd, x$beta.y1nd + linpred_y)
    sum(surv_vals * pi / pG)
  })
  
  SY0_ND <- sapply(y, function(h) {
    surv_vals <- Sweib(h, x$alpha.y0nd, x$beta.y0nd + linpred_y)
    sum(surv_vals * pi / pG)
  })
  
  DCE_ND <- SY1_ND - SY0_ND
  RMSTE_ND <- cumsum(DCE_ND)
  
  ###################
  #### DCE for D ####
  ###################
  
  linpred_d <- xs[, -1, drop = FALSE] %*% etas_d
  scale_fd <- exp(x$beta.d + linpred_d)
  
  fd <- matrix(0, nrow = n_obs, ncol = length(d))
  for (j in seq_along(d)) {
    fd[, j] <- dexp(d[j],  
                        rate = scale_fd)
  }
  
  # fD from cov grid
  linpred_d_cov <- as.matrix(cov) %*% etas_d
  scale_fD <- exp(x$beta.d + linpred_d_cov)
  
  fD <- matrix(0, nrow = nrow(cov), ncol = length(d))
  for (j in seq_along(d)) {
    fD[, j] <- dexp(d[j], rate = scale_fD)
  }
  
  fD_av <- colMeans(fD)
  
  ##########################
  ### DCE_Y_D  #######
  ##########################
  
  DCE_Y_D <- matrix(0, nrow = length(d), ncol = length(y))
  
  for (j in seq_along(d)) {
    log_j <- log(d[j])
    
    denom_fdj <- sum(fd[, j])  # sum over i for time d[j]
    
    for (h in seq_along(y)) {
      
      # Treatment
      surv1 <- Stweib(
        h,
        a = x$alpha.y1d,
        b = x$beta.y1d + linpred_y + x$delta * log_j,
        l = j
      )
      SY1_D <- sum(surv1 * fd[, j]) / denom_fdj
      
      # Control
      surv0 <- Sweib(
        h,
        x$alpha.y0d,
        x$beta.y0d + linpred_y + x$delta * log_j
      )
      SY0_D <- sum(surv0 * fd[, j]) / denom_fdj
      
      DCE_Y_D[j, h] <- SY1_D - SY0_D
    }
  }
  
  DCE_Y_D[is.nan(DCE_Y_D)] <- 0
  
  # Average across distributions fD_av
  DCE_D <- colSums(DCE_Y_D * matrix(fD_av, nrow = length(d), ncol = length(y), byrow = FALSE))
  
  # Cumulative estimates
  RMSTE_d <- t(apply(DCE_Y_D, 1, cumsum))
  RMSTE_D <- cumsum(DCE_D)
  
  return(list(
    DCE_ND   = DCE_ND,
    DCE_Y_D  = DCE_Y_D,
    DCE_D = DCE_D,
    RMSTE_ND = RMSTE_ND,
    RMSTE_D  = RMSTE_D,
    RMSTE_d  = RMSTE_d
  ))
}

bicfun_expweib <- function(x){
  
  data<-x$data
  n <- nrow(data)
  k<-length(x)-3
  
  Z  <- data$Z   #Treatment indicator
  ID <- x$ID  #Discontinuation indicator: #ID=1 for never discontinued
  RD <- data$RD  #RD=1 for treated units who discontinue
  D  <- x$D   #Time to discontinuation for those who discontinue
  RY <- data$RY  #RY=1 for units who die
  Y  <- data$Y   #survival time
  x1 <- data$x1
  x2 <- data$x2
  x3 <- data$x3
  
  xs<-cbind(1,as.matrix(data[,grepl("x",names(data))]))
  
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
  
  loglik<-0
  #Z=1 ID=1 RY=1 ---> Never discontinued who die under treatment
  ii <- Z*ID*RY==1
  loglik <- loglik + 
    sum({log(pi[ii]) + 
        log(dexp(Y[ii], 
                 rate=exp(x$beta.y1nd+
                            x$eta.y1*x1[ii]+
                            x$eta.y2*x2[ii]+
                            x$eta.y3*x3[ii])))},
        na.rm = TRUE)
  
  #Z=1 ID=1 RY=0 ---> Never discontinued who do not die under treatment
  ii <- Z*ID*(1-RY)==1
  loglik <- loglik + 
    sum({log(pi[ii]) + 
        log(Sexp(Y[ii], 
                 exp(x$beta.y1nd+
                       x$eta.y1*x1[ii]+
                       x$eta.y2*x2[ii]+
                       x$eta.y3*x3[ii])))},
        na.rm = TRUE)
  
  #Z=1 ID=0 RD=1 RY=1 ---> Discontinued who die under treatment for whom we 
  #observe the discontinuation
  ii <- Z*(1-ID)*RD*RY==1 
  loglik <- loglik + 
    sum({log(1-pi[ii]) + 
        log(dexp(D[ii], 
                 exp(x$beta.d+
                       x$eta.d1*x1[ii]+
                       x$eta.d2*x2[ii]+
                       x$eta.d3*x3[ii])))+
        log(dtexp(Y[ii],
                  rate=exp(x$beta.y1d+
                             x$eta.y1*x1[ii]+
                             x$eta.y2*x2[ii]+
                             x$eta.y3*x3[ii]+
                             x$delta*log(D[ii])),
                  a=D[ii]))},
        na.rm = TRUE)
  
  #Z=1 ID=0 RD=1 RY=0 ---> Discontinued who do not die under treatment for whom 
  #we observe the discontinuation
  ii <- Z*(1-ID)*RD*(1-RY)==1
  loglik <- loglik + 
    sum({log(1-pi[ii]) + 
        log(dexp(D[ii], 
                 exp(x$beta.d+
                       x$eta.d1*x1[ii]+
                       x$eta.d2*x2[ii]+
                       x$eta.d3*x3[ii])))+
        log(Stexp(Y[ii],
                  rate=exp(x$beta.y1d+
                             x$eta.y1*x1[ii]+
                             x$eta.y2*x2[ii]+
                             x$eta.y3*x3[ii]+
                             x$delta*log(D[ii])),
                  a=D[ii]))
    },
    na.rm = TRUE)
  
  #Z=1 ID=0 RD=0 RY=0 ---> Discontinued who discontinue after the end of the 
  #study and do not die under treatment
  ii <- Z*(1-ID)*(1-RD)*(1-RY)==1
  loglik <- loglik + 
    sum({log(1-pi[ii]) + 
        log(Sexp(D[ii], 
                 exp(x$beta.d+
                       x$eta.d1*x1[ii]+
                       x$eta.d2*x2[ii]+
                       x$eta.d3*x3[ii])))},
        na.rm = TRUE)
  
  
  #Z=0 ID=1 RY=1 ---> Never discontinued who die under control
  ii <- (1-Z)*ID*RY==1
  loglik <- loglik + 
    sum({log(pi[ii]) + 
        log(dexp(Y[ii], 
                 rate=exp(x$beta.y0nd+
                            x$eta.y1*x1[ii]+
                            x$eta.y2*x2[ii]+
                            x$eta.y3*x3[ii])))
    },na.rm = TRUE)
  
  #Z=0 ID=1 RY=0 ---> Never discontinued who do not die under control
  ii <- (1-Z)*ID*(1-RY)==1
  loglik <- loglik + 
    sum({log(pi[ii]) + 
        log(Sexp(Y[ii], 
                 exp(x$beta.y0nd+
                       x$eta.y1*x1[ii]+
                       x$eta.y2*x2[ii]+
                       x$eta.y3*x3[ii])))},
        na.rm = TRUE)
  
  #Z=0 ID=0 RY=1 ---> Discontinued who die under control
  ii <- (1-Z)*(1-ID)*RY==1
  loglik <- loglik + 
    sum({log(1-pi[ii]) + 
        log(dexp(D[ii],
                 exp(x$beta.d+
                       x$eta.d1*x1[ii]+
                       x$eta.d2*x2[ii]+
                       x$eta.d3*x3[ii])))+
        log(dexp(Y[ii], 
                 exp(x$beta.y0d+
                       x$delta*log(D[ii])+
                       x$eta.y1*x1[ii]+
                       x$eta.y2*x2[ii]+
                       x$eta.y3*x3[ii])))
    },na.rm = TRUE)
  
  #Z=0 ID=0 RY=0 ---> Discontinued who do not die under control
  ii <- (1-Z)*(1-ID)*(1-RY)==1
  loglik <- loglik + 
    sum({log(1-pi[ii]) + 
        log(dexp(D[ii],
                 exp(x$beta.d+
                       x$eta.d1*x1[ii]+
                       x$eta.d2*x2[ii]+
                       x$eta.d3*x3[ii])))+
        log(Sexp(Y[ii], 
                 exp(x$beta.y0d+
                       x$delta*log(D[ii])+
                       x$eta.y1*x1[ii]+
                       x$eta.y2*x2[ii]+
                       x$eta.y3*x3[ii])))},
        na.rm = TRUE)
  
  bic<- -2*loglik + k*log(n)
  return(bic=bic)
}

waicfun_expweib <- function(x){
  
  data<-x$data
  n <- nrow(data)
  
  Z  <- data$Z   #Treatment indicator
  ID <- x$ID  #Discontinuation indicator: #ID=1 for never discontinued
  RD <- data$RD  #RD=1 for treated units who discontinue
  D  <- x$D   #Time to discontinuation for those who discontinue
  RY <- data$RY  #RY=1 for units who die
  Y  <- data$Y   #survival time
  x1 <- data$x1
  x2 <- data$x2
  x3 <- data$x3
  
  xs<-cbind(1,as.matrix(data[,grepl("x",names(data))]))
  
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
  
  loglik<-rep(0, n)
  #Z=1 ID=1 RY=1 ---> Never discontinued who die under treatment
  ii<-Z*ID*RY==1 
  loglik[ii] <- log(pi[ii]) + 
    log(dweib(Y[ii], 
              x$alpha.y1nd, 
              x$beta.y1nd+
                x$eta.y1*x1[ii]+
                x$eta.y2*x2[ii]+
                x$eta.y3*x3[ii]))
  
  #Z=1 ID=1 RY=0 ---> Never discontinued who do not die under treatment
  ii<-Z*ID*(1-RY)==1
  loglik[ii] <- log(pi[ii]) + 
    log(Sweib(Y[ii], 
              x$alpha.y1nd, 
              x$beta.y1nd+
                x$eta.y1*x1[ii]+
                x$eta.y2*x2[ii]+
                x$eta.y3*x3[ii]))
  
  #Z=1 ID=0 RD=1 RY=1 ---> Discontinued who die under treatment for whom we 
  #observe the discontinuation
  ii <- Z*(1-ID)*RD*RY==1
  loglik[ii] <- log(1-pi[ii]) + 
    log(dexp(D[ii], 
             exp(x$beta.d+
                   x$eta.d1*x1[ii]+
                   x$eta.d2*x2[ii]+
                   x$eta.d3*x3[ii])))+
    # log(dtweibull(Y[ii],
    #           shape=x$alpha.y1d,
    #           scale=exp(-(x$beta.y1d+
    #                               x$eta.y1*x1[ii]+
    #                               x$eta.y2*x2[ii]+
    #                               x$eta.y3*x3[ii]+
    #                               x$delta*log(D[ii]))/x$alpha.y1d),
    #           a=D[ii]))
    log(dtweib(Y[ii],
               sh=x$alpha.y1d,
               sc=exp(-(x$beta.y1d+
                          x$eta.y1*x1[ii]+
                          x$eta.y2*x2[ii]+
                          x$eta.y3*x3[ii]+
                          x$delta*log(D[ii]))/x$alpha.y1d),
               a=D[ii]))
  
  
  #Z=1 ID=0 RD=1 RY=0 ---> Discontinued who do not die under treatment for whom 
  #we observe the discontinuation
  ii<-Z*(1-ID)*RD*(1-RY)==1
  loglik[ii] <- log(1-pi[ii]) + 
    log(dexp(D[ii], 
             exp(x$beta.d+
                   x$eta.d1*x1[ii]+
                   x$eta.d2*x2[ii]+
                   x$eta.d3*x3[ii])))+
    # log(Stweib(Y[ii],
    #           sh=x$alpha.y1d,
    #           sc=exp(-(x$beta.y1d+
    #                         x$eta.y1*x1[ii]+
    #                         x$eta.y2*x2[ii]+
    #                         x$eta.y3*x3[ii]+
    #                         x$delta*log(D[ii]))/x$alpha.y1d),
    #           a=D[ii]))
    log(Stweib(Y[ii],
               a=x$alpha.y1d,
               b=x$beta.y1d+
                 x$eta.y1*x1[ii]+
                 x$eta.y2*x2[ii]+
                 x$eta.y3*x3[ii]+
                 x$delta*log(D[ii]),
               l=D[ii]))
  
  
  #Z=1 ID=0 RD=0 RY=0 ---> Discontinued who discontinue after the end of the 
  #study and do not die under treatment
  ii<-Z*(1-ID)*(1-RD)*(1-RY)==1
  loglik[ii] <- log(1-pi[ii]) + 
    log(Sexp(D[ii], 
             exp(x$beta.d+
                   x$eta.d1*x1[ii]+
                   x$eta.d2*x2[ii]+
                   x$eta.d3*x3[ii])))
  
  
  #Z=0 ID=1 RY=1 ---> Never discontinued who die under control
  ii <- (1-Z)*ID*RY==1
  loglik[ii] <- log(pi[ii]) + 
    log(dweib(Y[ii], 
              x$alpha.y0nd, 
              x$beta.y0nd+
                x$eta.y1*x1[ii]+
                x$eta.y2*x2[ii]+
                x$eta.y3*x3[ii]))
  
  #Z=0 ID=1 RY=0 Y>Y1*k ---> Never discontinued who do not die under control
  ii <- (1-Z)*ID*(1-RY)==1
  loglik[ii] <- log(pi[ii]) + 
    log(Sweib(Y[ii], 
              x$alpha.y0nd, 
              x$beta.y0nd+
                x$eta.y1*x1[ii]+
                x$eta.y2*x2[ii]+
                x$eta.y3*x3[ii]))
  
  #Z=0 ID=0 RY=1 ---> Discontinued who die under control
  ii <- (1-Z)*(1-ID)*RY==1
  loglik[ii] <- log(1-pi[ii]) + 
    log(dexp(D[ii],
             exp(x$beta.d+
                   x$eta.d1*x1[ii]+
                   x$eta.d2*x2[ii]+
                   x$eta.d3*x3[ii])))+
    log(dweib(Y[ii], 
              x$alpha.y0d, 
              {x$beta.y0d+
                  x$delta*log(D[ii])+
                  x$eta.y1*x1[ii]+
                  x$eta.y2*x2[ii]+
                  x$eta.y3*x3[ii]}))
  
  #Z=0 ID=0  RY=0 ---> Discontinued who do not die under control
  ii <- (1-Z)*(1-ID)*(1-RY)==1
  loglik[ii] <- log(1-pi[ii]) + 
    log(dexp(D[ii],
             exp(x$beta.d+
                   x$eta.d1*x1[ii]+
                   x$eta.d2*x2[ii]+
                   x$eta.d3*x3[ii])))+
    log(Sweib(Y[ii], 
              x$alpha.y0d, 
              {x$beta.y0d+
                  x$delta*log(D[ii])+
                  x$eta.y1*x1[ii]+
                  x$eta.y2*x2[ii]+
                  x$eta.y3*x3[ii]}))
  
  return(loglik=loglik)
}

