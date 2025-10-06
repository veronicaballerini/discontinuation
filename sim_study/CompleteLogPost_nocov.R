nc_logpost <- function(theta, theta.prior, complete.data) {

    Z <- complete.data$Z  #Treatment indicator
    ID <- complete.data$ID  #Discontinuation indicator: #ID=1 for ND patients
    RD <- complete.data$RD  #RD=1 for treated units who discontinue
    D <- complete.data$D  #Time to discontinuation for D patients
    RY <- complete.data$RY  #RY=1 for units who die
    Y <- complete.data$Y  #survival time

    ############################ LOG LIKELIHOOD ##################################

    loglik <- 0
    # Z=1 ID=1 RY=1 ---> Never discontinued who die under treatment
    ii <- Z * ID * RY == 1
    loglik <- loglik + sum({
        log(theta$pi) + log(dweib(Y[ii], theta$alpha.y1nd, theta$beta.y1nd))
    }, na.rm = TRUE)

    # Z=1 ID=1 RY=0 ---> Never discontinued who do not die under treatment
    ii <- Z * ID * (1 - RY) == 1
    loglik <- loglik + sum({
        log(theta$pi) + log(Sweib(Y[ii], theta$alpha.y1nd, theta$beta.y1nd))
    }, na.rm = TRUE)

    # Z=1 ID=0 RD=1 RY=1 ---> Discontinued who die under treatment for whom we observe the discontinuation
    ii <- Z * (1 - ID) * RD * RY == 1
    loglik <- loglik + sum({
        log(1 - theta$pi) + log(dweib(D[ii], theta$alpha.d, theta$beta.d)) + log(dtweib(Y[ii], sh = theta$alpha.y1d, sc = exp(-(theta$beta.y1d + theta$delta * log(D[ii]))/theta$alpha.y1d),
            a = D[ii]))
    }, na.rm = TRUE)

    # Z=1 ID=0 RD=1 RY=0 ---> Discontinued who do not die under treatment for whom we observe the discontinuation
    ii <- Z * (1 - ID) * RD * (1 - RY) == 1
    loglik <- loglik + sum({
        log(1 - theta$pi) + log(dweib(D[ii], theta$alpha.d, theta$beta.d)) + log(Stweib(Y[ii], a = theta$alpha.y1d, b = theta$beta.y1d + theta$delta * log(D[ii]),
            l = D[ii]))
    }, na.rm = TRUE)

    # Z=1 ID=0 RD=0 RY=0 ---> Discontinued who discontinue after the end of the study and do not die under treatment
    ii <- Z * (1 - ID) * (1 - RD) * (1 - RY) == 1
    loglik <- loglik + sum({
        log(1 - theta$pi) + log(Sweib(D[ii], theta$alpha.d, theta$beta.d))
    }, na.rm = TRUE)


    # Z=0 ID=1 RY=1 ---> Never discontinued who die under control
    ii <- (1 - Z) * ID * RY == 1
    loglik <- loglik + sum({
        log(theta$pi) + log(dweib(Y[ii], theta$alpha.y0nd, theta$beta.y0nd))
    }, na.rm = TRUE)

    # Z=0 ID=1 RY=0 ---> Never discontinued who do not die under control
    ii <- (1 - Z) * ID * (1 - RY) == 1
    loglik <- loglik + sum({
        log(theta$pi) + log(Sweib(Y[ii], theta$alpha.y0nd, theta$beta.y0nd))
    }, na.rm = TRUE)

    # Z=0 ID=0 RY=1 ---> Discontinued who die under control
    ii <- (1 - Z) * (1 - ID) * RY == 1
    loglik <- loglik + sum({
        log(1 - theta$pi) + log(dweib(D[ii], theta$alpha.d, theta$beta.d)) + log(dweib(Y[ii], theta$alpha.y0d, {
            theta$beta.y0d + theta$delta * log(D[ii])
        }))
    }, na.rm = TRUE)

    # Z=0 ID=0 RY=0 ---> Discontinued who do not die under control
    ii <- (1 - Z) * (1 - ID) * (1 - RY) == 1
    loglik <- loglik + sum({
        log(1 - theta$pi) + log(dweib(D[ii], theta$alpha.d, theta$beta.d)) + log(Sweib(Y[ii], theta$alpha.y0d, {
            theta$beta.y0d + theta$delta * log(D[ii])
        }))
    }, na.rm = TRUE)

    ############################### LOG PRIOR ####################################

    logprior <- 0

    logprior <- logprior + dnorm(theta$eta.0, theta.prior$mu.0, theta.prior$sigma2.0, log = TRUE)

    logprior <- logprior + dgamma(theta$alpha.d, shape = theta.prior$a.d, scale = theta.prior$b.d, log = TRUE) + dnorm(theta$beta.d, theta.prior$mu.d, sqrt(theta.prior$sigma2.d),
        log = TRUE)

    logprior <- logprior + dgamma(theta$alpha.y1nd, shape = theta.prior$a.y1nd, scale = theta.prior$b.y1nd, log = TRUE) + dnorm(theta$beta.y1nd, theta.prior$mu.y1nd,
        sqrt(theta.prior$sigma2.y1nd), log = TRUE)

    logprior <- logprior + dgamma(theta$alpha.y1d, shape = theta.prior$a.y1d, scale = theta.prior$b.y1d, log = TRUE) + dnorm(theta$beta.y1d, theta.prior$mu.y1d,
        sqrt(theta.prior$sigma2.y1d), log = TRUE)

    logprior <- logprior + dgamma(theta$alpha.y0nd, shape = theta.prior$a.y0nd, scale = theta.prior$b.y0nd, log = TRUE) + dnorm(theta$beta.y0nd, theta.prior$mu.y0nd,
        sqrt(theta.prior$sigma2.y0nd), log = TRUE)

    logprior <- logprior + dgamma(theta$alpha.y0d, shape = theta.prior$a.y0d, scale = theta.prior$b.y0d, log = TRUE) + dnorm(theta$beta.y0d, theta.prior$mu.y0d,
        sqrt(theta.prior$sigma2.y0d), log = TRUE)

    logprior <- logprior + dnorm(theta$delta, theta.prior$mu.delta, sqrt(theta.prior$sigma2.delta), log = TRUE)

    logprior + loglik
}


nc_logpost.i <- function(theta, theta.prior, complete.data) {

    n <- nrow(complete.data)
    Z <- complete.data$Z
    ID <- complete.data$ID
    RD <- complete.data$RD
    D <- complete.data$D
    RY <- complete.data$RY
    Y <- complete.data$Y

    ############################ LOG LIKELIHOOD ##################################

    loglik <- rep(0, n)
    # Z=1 ID=1 RY=1 ---> Never discontinued who die under treatment
    ii <- Z * ID * RY == 1
    loglik[ii] <- log(theta$pi) + log(dweib(Y[ii], theta$alpha.y1nd, theta$beta.y1nd))

    # Z=1 ID=1 RY=0 ---> Never discontinued who do not die under treatment
    ii <- Z * ID * (1 - RY) == 1
    loglik[ii] <- log(theta$pi) + log(Sweib(Y[ii], theta$alpha.y1nd, theta$beta.y1nd))

    # Z=1 ID=0 RD=1 RY=1 ---> Discontinued who die under treatment for whom we observe the discontinuation
    ii <- Z * (1 - ID) * RD * RY == 1
    loglik[ii] <- log(1 - theta$pi) + log(dweib(D[ii], theta$alpha.d, theta$beta.d)) + log(dtweib(Y[ii], sh = theta$alpha.y1d, sc = exp(-(theta$beta.y1d + theta$delta *
        log(D[ii]))/theta$alpha.y1d), a = D[ii]))


    # Z=1 ID=0 RD=1 RY=0 ---> Discontinued who do not die under treatment for whom we observe the discontinuation
    ii <- Z * (1 - ID) * RD * (1 - RY) == 1
    loglik[ii] <- log(1 - theta$pi) + log(dweib(D[ii], theta$alpha.d, theta$beta.d)) + log(Stweib(Y[ii], a = theta$alpha.y1d, b = theta$beta.y1d + theta$delta *
        log(D[ii]), l = D[ii]))

    # Z=1 ID=0 RD=0 RY=0 ---> Discontinued who discontinue after the end of the study and do not die under treatment
    ii <- Z * (1 - ID) * (1 - RD) * (1 - RY) == 1
    loglik[ii] <- log(1 - theta$pi) + log(Sweib(D[ii], theta$alpha.d, theta$beta.d))


    # Z=0 ID=1 RY=1 ---> Never discontinued who die under control
    ii <- (1 - Z) * ID * RY == 1
    loglik[ii] <- log(theta$pi) + log(dweib(Y[ii], theta$alpha.y0nd, theta$beta.y0nd))

    # Z=0 ID=1 RY=0 Y>Y1*k ---> Never discontinued who do not die under control
    ii <- (1 - Z) * ID * (1 - RY) == 1
    loglik[ii] <- log(theta$pi) + log(Sweib(Y[ii], theta$alpha.y0nd, theta$beta.y0nd))

    # Z=0 ID=0 RY=1 ---> Discontinued who die under control
    ii <- (1 - Z) * (1 - ID) * RY == 1
    loglik[ii] <- log(1 - theta$pi) + log(dweib(D[ii], theta$alpha.d, theta$beta.d)) + log(dweib(Y[ii], theta$alpha.y0d, {
        theta$beta.y0d + theta$delta * log(D[ii])
    }))

    # Z=0 ID=0 RY=0 ---> Discontinued who do not die under control
    ii <- (1 - Z) * (1 - ID) * (1 - RY) == 1
    loglik[ii] <- log(1 - theta$pi) + log(dweib(D[ii], theta$alpha.d, theta$beta.d)) + log(Sweib(Y[ii], theta$alpha.y0d, {
        theta$beta.y0d + theta$delta * log(D[ii])
    }))

    rm(ii)

    ############################### LOG PRIOR ####################################

    logprior <- 0

    logprior <- logprior + dnorm(theta$eta.0, theta.prior$mu.0, theta.prior$sigma2.0, log = TRUE)

    logprior <- logprior + dgamma(theta$alpha.d, shape = theta.prior$a.d, scale = theta.prior$b.d, log = TRUE) + dnorm(theta$beta.d, theta.prior$mu.d, sqrt(theta.prior$sigma2.d),
        log = TRUE)

    logprior <- logprior + dgamma(theta$alpha.y1nd, shape = theta.prior$a.y1nd, scale = theta.prior$b.y1nd, log = TRUE) + dnorm(theta$beta.y1nd, theta.prior$mu.y1nd,
        sqrt(theta.prior$sigma2.y1nd), log = TRUE)

    logprior <- logprior + dgamma(theta$alpha.y1d, shape = theta.prior$a.y1d, scale = theta.prior$b.y1d, log = TRUE) + dnorm(theta$beta.y1d, theta.prior$mu.y1d,
        sqrt(theta.prior$sigma2.y1d), log = TRUE)

    logprior <- logprior + dgamma(theta$alpha.y0nd, shape = theta.prior$a.y0nd, scale = theta.prior$b.y0nd, log = TRUE) + dnorm(theta$beta.y0nd, theta.prior$mu.y0nd,
        sqrt(theta.prior$sigma2.y0nd), log = TRUE)

    logprior <- logprior + dgamma(theta$alpha.y0d, shape = theta.prior$a.y0d, scale = theta.prior$b.y0d, log = TRUE) + dnorm(theta$beta.y0d, theta.prior$mu.y0d,
        sqrt(theta.prior$sigma2.y0d), log = TRUE)

    logprior <- logprior + dnorm(theta$delta, theta.prior$mu.delta, sqrt(theta.prior$sigma2.delta), log = TRUE)

    logprior + loglik
}
