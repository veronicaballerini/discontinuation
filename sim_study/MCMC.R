source("CompleteLogPost.R")
source("DataAugmentation.R")
library(mvtnorm)
library(HDInterval)

mcmc.tdiscontinue <- function(niter, burn, thin, par.start, theta.prior, observed.data, scenario, seed) {

    theta.start <- list(eta.0 = rnorm(1, mean = par.start$mu.0, sd = par.start$sigma.0), eta.1 = rnorm(1, mean = par.start$mu.1, sd = par.start$sigma.1), eta.2 = rnorm(1,
        mean = par.start$mu.2, sd = par.start$sigma.2), eta.3 = rnorm(1, mean = par.start$mu.3, sd = par.start$sigma.3), alpha.d = par.start$alpha.d + abs(rnorm(1,
        0, 0.1)), beta.d = par.start$beta.d + rnorm(1, 0, 0.1), eta.d1 = rnorm(1, mean = par.start$mu.d1, sd = par.start$sigma.d1), eta.d2 = rnorm(1, mean = par.start$mu.d2,
        sd = par.start$sigma.d2), eta.d3 = rnorm(1, mean = par.start$mu.d3, sd = par.start$sigma.d3), alpha.y1nd = par.start$shape.y1nd * par.start$scale.y1nd, beta.y1nd = rnorm(1,
        mean = par.start$mu.y1nd, sd = par.start$sigma.y1nd), alpha.y1d = par.start$shape.y1d * par.start$scale.y1d, beta.y1d = rnorm(1, mean = par.start$mu.y1d,
        sd = par.start$sigma.y1d), alpha.y0nd = par.start$shape.y0nd * par.start$scale.y0nd, beta.y0nd = rnorm(1, mean = par.start$mu.y0nd, sd = par.start$sigma.y0nd),
        alpha.y0d = par.start$shape.y0d * par.start$scale.y0d, beta.y0d = rnorm(1, mean = par.start$mu.y0d, sd = par.start$sigma.y0d), eta.y1 = rnorm(1, mean = par.start$mu.y1,
            sd = par.start$sigma.y1), eta.y2 = rnorm(1, mean = par.start$mu.y2, sd = par.start$sigma.y2), eta.y3 = rnorm(1, mean = par.start$mu.y3, sd = par.start$sigma.y3),
        delta = rnorm(1, mean = par.start$mu.delta, sd = par.start$sigma.delta))

    pi = exp(theta.start$eta.0 + theta.start$eta.1 * observed.data$x1st + theta.start$eta.2 * observed.data$x2 + theta.start$eta.3 * observed.data$x3)/(1 + exp(theta.start$eta.0 +
        theta.start$eta.1 * observed.data$x1st + theta.start$eta.2 * observed.data$x2 + theta.start$eta.3 * observed.data$x3))


    ID.start <- D.start <- NULL
    ID.start[observed.data$Z == 1 & observed.data$RD.obs == 0 & observed.data$RY == 1] <- 1
    ID.start[observed.data$Z == 1 & observed.data$RD.obs == 1] <- 0
    ID.start[observed.data$Z == 1 & observed.data$RD.obs == 0 & observed.data$RY == 0] <- rbinom(sum((observed.data$Z) * (1 - observed.data$RD.obs) * (1 - observed.data$RY)),
        1, pi[observed.data$Z == 1 & observed.data$RD.obs == 0 & observed.data$RY == 0])
    ID.start[observed.data$Z == 0] <- rbinom(sum(observed.data$Z == 0), 1, pi[observed.data$Z == 0])

    D.start[ID.start == 1] <- 0
    D.start[ID.start == 0 & observed.data$Z == 1] <- observed.data$Dobs[ID.start == 0 & observed.data$Z == 1]

    D.start[ID.start == 0 & observed.data$Z == 0] <- rweibull(sum(ID.start == 0 & observed.data$Z == 0), shape = theta.start$alpha.d, scale = exp(-(theta.start$beta.d +
        theta.start$eta.d1 * observed.data$x1st[ID.start == 0 & observed.data$Z == 0] + theta.start$eta.d2 * observed.data$x2[ID.start == 0 & observed.data$Z ==
        0] + theta.start$eta.d3 * observed.data$x3[ID.start == 0 & observed.data$Z == 0])/theta.start$alpha.d))

    complete.data <- data.frame(cens = observed.data$f.up, Z = observed.data$Z, ID = ID.start, RD = observed.data$RD.obs, D = D.start, RY = observed.data$RY, Y = observed.data$Y,
        x1 = observed.data$x1st, x2 = observed.data$x2, x3 = observed.data$x3)

    theta <- list(eta.0 = theta.start$eta.0, eta.1 = theta.start$eta.1, eta.2 = theta.start$eta.2, eta.3 = theta.start$eta.3, alpha.d = theta.start$alpha.d, beta.d = theta.start$beta.d,
        eta.d1 = theta.start$eta.d1, eta.d2 = theta.start$eta.d2, eta.d3 = theta.start$eta.d3, alpha.y1nd = theta.start$alpha.y1nd, beta.y1nd = theta.start$beta.y1nd,
        alpha.y1d = theta.start$alpha.y1d, beta.y1d = theta.start$beta.y1d, alpha.y0nd = theta.start$alpha.y0nd, beta.y0nd = theta.start$beta.y0nd, alpha.y0d = theta.start$alpha.y0d,
        beta.y0d = theta.start$beta.y0d, eta.y1 = theta.start$eta.y1, eta.y2 = theta.start$eta.y2, eta.y3 = theta.start$eta.y3, pi = pi, delta = theta.start$delta)

    thetanop <- theta
    thetanop[[21]] <- NULL
    npar <- length(unlist(thetanop))
    jump <- matrix(NA, niter, npar)
    Theta <- matrix(NA, niter, npar)
    DID <- list()

    colnames(jump) <- names(unlist(thetanop))
    colnames(Theta) <- names(unlist(thetanop))

    nojump <- NA

    for (j in 1:niter) {


        ID.old <- complete.data$ID
        D.old <- complete.data$D
        if (j%%1000 == 0) {
            print(j)
        }
        discontinuedata <- da.discontinuation(theta, theta.prior, observed.data, ID.old, D.old)
        complete.data$ID <- discontinuedata$ID
        complete.data$D <- discontinuedata$D

        rm(ID.old, D.old)

        x <- as.matrix(cbind(1, observed.data$x1st, observed.data$x2, observed.data$x3))

        ### eta0
        if (j == 1) {
            lp_sd0 <- log(proposal$sd.0)
        }
        eta.0 <- rnorm(1, theta$eta.0, exp(lp_sd0))
        eta <- as.matrix(c(eta.0, theta$eta.1, theta$eta.2, theta$eta.3))

        thetaprop <- theta
        thetaprop$eta.0 <- eta.0
        thetaprop$pi <- as.vector(exp(x %*% eta)/(1 + exp(x %*% eta)))

        log.post.eta.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.eta.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.eta.c
        pr.den <- log.post.eta.old
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        if (u <= pr & is.nan(pr) == FALSE) {
            theta$eta.0 <- eta.0
            theta$pi <- thetaprop$pi
        }

        jump[j, 1] <- as.numeric(u <= pr & is.nan(pr) == FALSE)
        rm(eta.0, thetaprop, u, pr.num, pr.den, pr, log.post.eta.c, log.post.eta.old)

        ### eta1
        if (j == 1) {
            lp_sd1 <- log(proposal$sd.1)
        }
        eta.1 <- rnorm(1, theta$eta.1, exp(lp_sd1))
        eta <- as.matrix(c(theta$eta.0, eta.1, theta$eta.2, theta$eta.3))

        thetaprop <- theta
        thetaprop$eta.1 <- eta.1
        thetaprop$pi <- as.vector(exp(x %*% eta)/(1 + exp(x %*% eta)))

        log.post.eta.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.eta.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.eta.c
        pr.den <- log.post.eta.old
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        if (u <= pr & is.nan(pr) == FALSE) {
            theta$eta.1 <- eta.1
            theta$pi <- thetaprop$pi
        }

        jump[j, 2] <- as.numeric(u <= pr & is.nan(pr) == FALSE)
        rm(eta.1, thetaprop, u, pr.num, pr.den, pr, log.post.eta.c, log.post.eta.old)

        ### eta2
        if (j == 1) {
            lp_sd2 <- log(proposal$sd.2)
        }
        eta.2 <- rnorm(1, theta$eta.2, exp(lp_sd2))
        eta <- as.matrix(c(theta$eta.0, theta$eta.1, eta.2, theta$eta.3))

        thetaprop <- theta
        thetaprop$eta.2 <- eta.2
        thetaprop$pi <- as.vector(exp(x %*% eta)/(1 + exp(x %*% eta)))

        log.post.eta.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.eta.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.eta.c
        pr.den <- log.post.eta.old
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        if (u <= pr & is.nan(pr) == FALSE) {
            theta$eta.2 <- eta.2
            theta$pi <- thetaprop$pi
        }

        jump[j, 3] <- as.numeric(u <= pr & is.nan(pr) == FALSE)
        rm(eta.2, thetaprop, u, pr.num, pr.den, pr, log.post.eta.c, log.post.eta.old)

        ### eta3
        if (j == 1) {
            lp_sd3 <- log(proposal$sd.3)
        }

        eta.3 <- rnorm(1, theta$eta.3, exp(lp_sd3))
        eta <- as.matrix(c(theta$eta.0, theta$eta.1, theta$eta.2, eta.3))

        thetaprop <- theta
        thetaprop$eta.3 <- eta.3
        thetaprop$pi <- as.vector(exp(x %*% eta)/(1 + exp(x %*% eta)))

        log.post.eta.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.eta.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.eta.c
        pr.den <- log.post.eta.old
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        if (u <= pr & is.nan(pr) == FALSE) {
            theta$eta.3 <- eta.3
            theta$pi <- thetaprop$pi
        }

        jump[j, 4] <- as.numeric(u <= pr & is.nan(pr) == FALSE)
        rm(eta.3, x, thetaprop, u, pr.num, pr.den, pr, log.post.eta.c, log.post.eta.old)

        ### alpha.d
        alpha.d <- rgamma(1, shape = theta$alpha.d/proposal$scale.d, scale = proposal$scale.d)

        thetaprop <- theta
        thetaprop$alpha.d <- alpha.d

        log.post.alpha.d.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.alpha.d.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.alpha.d.c + log(dgamma(theta$alpha.d, shape = alpha.d/proposal$scale.d, scale = proposal$scale.d))
        pr.den <- log.post.alpha.d.old + log(dgamma(alpha.d, shape = theta$alpha.d/proposal$scale.d, scale = proposal$scale.d))
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        theta$alpha.d <- ifelse((u <= pr & is.nan(pr) == FALSE), alpha.d, theta$alpha.d)
        jump[j, 5] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### beta.d
        if (j == 1) {
            lp_sd.d <- log(proposal$sd.d)
        }
        beta.d <- rnorm(1, theta$beta.d, exp(lp_sd.d))

        thetaprop <- theta
        thetaprop$beta.d <- beta.d

        log.post.beta.d.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.beta.d.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.beta.d.c
        pr.den <- log.post.beta.d.old

        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        theta$beta.d <- ifelse((u <= pr & is.nan(pr) == FALSE), beta.d, theta$beta.d)
        jump[j, 6] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### eta.d1
        if (j == 1) {
            lp_sd.d1 <- log(proposal$sd.d1)
        }

        eta.d1 <- rnorm(1, theta$eta.d1, exp(lp_sd.d1))

        thetaprop <- theta
        thetaprop$eta.d1 <- eta.d1

        log.post.eta.d.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.eta.d.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.eta.d.c
        pr.den <- log.post.eta.d.old
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        if (u <= pr & is.nan(pr) == FALSE) {
            theta$eta.d1 <- eta.d1
        }

        jump[j, 7] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### eta.d2
        if (j == 1) {
            lp_sd.d2 <- log(proposal$sd.d2)
        }

        eta.d2 <- rnorm(1, theta$eta.d2, exp(lp_sd.d2))

        thetaprop <- theta
        thetaprop$eta.d2 <- eta.d2

        log.post.eta.d.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.eta.d.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.eta.d.c
        pr.den <- log.post.eta.d.old
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        if (u <= pr & is.nan(pr) == FALSE) {
            theta$eta.d2 <- eta.d2
        }

        jump[j, 8] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### eta.d3
        if (j == 1) {
            lp_sd.d3 <- log(proposal$sd.d3)
        }

        eta.d3 <- rnorm(1, theta$eta.d3, exp(lp_sd.d3))

        thetaprop <- theta
        thetaprop$eta.d3 <- eta.d3

        log.post.eta.d.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.eta.d.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.eta.d.c
        pr.den <- log.post.eta.d.old
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        if (u <= pr & is.nan(pr) == FALSE) {
            theta$eta.d3 <- eta.d3
        }

        jump[j, 9] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### alpha.y1nd
        alpha.y1nd <- rgamma(1, shape = theta$alpha.y1nd/proposal$scale.y1nd, scale = proposal$scale.y1nd)

        thetaprop <- theta
        thetaprop$alpha.y1nd <- alpha.y1nd

        log.post.alpha.y1nd.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.alpha.y1nd.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.alpha.y1nd.c + log(dgamma(theta$alpha.y1nd, shape = alpha.y1nd/proposal$scale.y1nd, scale = proposal$scale.y1nd))
        pr.den <- log.post.alpha.y1nd.old + log(dgamma(alpha.y1nd, shape = theta$alpha.y1nd/proposal$scale.y1nd, scale = proposal$scale.y1nd))

        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        theta$alpha.y1nd <- ifelse((u <= pr & is.nan(pr) == FALSE), alpha.y1nd, theta$alpha.y1nd)
        jump[j, 10] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### beta.y1nd
        if (j == 1) {
            lp_sd.y1nd <- log(proposal$sd.y1nd)
        }

        beta.y1nd <- rnorm(1, theta$beta.y1nd, exp(lp_sd.y1nd))

        thetaprop <- theta
        thetaprop$beta.y1nd <- beta.y1nd

        log.post.beta.y1nd.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.beta.y1nd.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.beta.y1nd.c
        pr.den <- log.post.beta.y1nd.old
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        theta$beta.y1nd <- ifelse((u <= pr & is.nan(pr) == FALSE), beta.y1nd, theta$beta.y1nd)
        jump[j, 11] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### alpha.y1d
        alpha.y1d <- rgamma(1, shape = theta$alpha.y1d/proposal$scale.y1d, scale = proposal$scale.y1d)

        thetaprop <- theta
        thetaprop$alpha.y1d <- alpha.y1d

        log.post.alpha.y1d.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.alpha.y1d.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.alpha.y1d.c + log(dgamma(theta$alpha.y1d, shape = alpha.y1d/proposal$scale.y1d, scale = proposal$scale.y1d))
        pr.den <- log.post.alpha.y1d.old + log(dgamma(alpha.y1d, shape = theta$alpha.y1d/proposal$scale.y1d, scale = proposal$scale.y1d))
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        theta$alpha.y1d <- ifelse((u <= pr & is.nan(pr) == FALSE), alpha.y1d, theta$alpha.y1d)
        jump[j, 12] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### beta.y1d
        if (j == 1) {
            lp_sd.y1d <- log(proposal$sd.y1d)
        }

        beta.y1d <- rnorm(1, theta$beta.y1d, exp(lp_sd.y1d))

        thetaprop <- theta
        thetaprop$beta.y1d <- beta.y1d

        log.post.beta.y1d.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.beta.y1d.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.beta.y1d.c
        pr.den <- log.post.beta.y1d.old

        pr <- exp(pr.num - pr.den)
        u <- runif(1)
        theta$beta.y1d <- ifelse((u <= pr & is.nan(pr) == FALSE), beta.y1d, theta$beta.y1d)
        jump[j, 13] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### alpha.y0nd
        alpha.y0nd <- rgamma(1, shape = theta$alpha.y0nd/proposal$scale.y0nd, scale = proposal$scale.y0nd)


        thetaprop <- theta
        thetaprop$alpha.y0nd <- alpha.y0nd

        log.post.alpha.y0nd.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.alpha.y0nd.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.alpha.y0nd.c + log(dgamma(theta$alpha.y0nd, shape = alpha.y0nd/proposal$scale.y0nd, scale = proposal$scale.y0nd))
        pr.den <- log.post.alpha.y0nd.old + log(dgamma(alpha.y0nd, shape = theta$alpha.y0nd/proposal$scale.y0nd, scale = proposal$scale.y0nd))
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        theta$alpha.y0nd <- ifelse((u <= pr & is.nan(pr) == FALSE), alpha.y0nd, theta$alpha.y0nd)
        jump[j, 14] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### beta.y0nd
        if (j == 1) {
            lp_sd.y0nd <- log(proposal$sd.y0nd)
        }

        beta.y0nd <- rnorm(1, theta$beta.y0nd, exp(lp_sd.y0nd))

        thetaprop <- theta
        thetaprop$beta.y0nd <- beta.y0nd

        log.post.beta.y0nd.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.beta.y0nd.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.beta.y0nd.c
        pr.den <- log.post.beta.y0nd.old

        pr <- exp(pr.num - pr.den)
        u <- runif(1)
        theta$beta.y0nd <- ifelse((u <= pr & is.nan(pr) == FALSE), beta.y0nd, theta$beta.y0nd)
        jump[j, 15] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### alpha.y0d
        alpha.y0d <- rgamma(1, shape = theta$alpha.y0d/proposal$scale.y0d, scale = proposal$scale.y0d)

        thetaprop <- theta
        thetaprop$alpha.y0d <- alpha.y0d

        log.post.alpha.y0d.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.alpha.y0d.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.alpha.y0d.c + log(dgamma(theta$alpha.y0d, shape = alpha.y0d/proposal$scale.y0d, scale = proposal$scale.y0d))
        pr.den <- log.post.alpha.y0d.old + log(dgamma(alpha.y0d, shape = theta$alpha.y0d/proposal$scale.y0d, scale = proposal$scale.y0d))
        pr <- exp(pr.num - pr.den)

        u <- runif(1)

        theta$alpha.y0d <- ifelse((u <= pr & is.nan(pr) == FALSE), alpha.y0d, theta$alpha.y0d)
        jump[j, 16] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### beta.y0d
        if (j == 1) {
            lp_sd.y0d <- log(proposal$sd.y0d)
        }

        beta.y0d <- rnorm(1, theta$beta.y0d, exp(lp_sd.y0d))

        thetaprop <- theta
        thetaprop$beta.y0d <- beta.y0d

        log.post.beta.y0d.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.beta.y0d.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.beta.y0d.c
        pr.den <- log.post.beta.y0d.old
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        theta$beta.y0d <- ifelse((u <= pr & is.nan(pr) == FALSE), beta.y0d, theta$beta.y0d)
        jump[j, 17] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        # y1
        if (j == 1) {
            lp_sd.y1 <- log(proposal$sd.y1)
        }

        eta1 <- rnorm(1, theta$eta.y1, exp(lp_sd.y1))

        thetaprop <- theta
        thetaprop$eta.y1 <- eta1

        log.post.eta.y.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.eta.y.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.eta.y.c
        pr.den <- log.post.eta.y.old
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        if (u <= pr & is.nan(pr) == FALSE) {
            theta$eta.y1 <- eta1
        }

        jump[j, 18] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        # y2
        if (j == 1) {
            lp_sd.y2 <- log(proposal$sd.y2)
        }

        eta2 <- rnorm(1, theta$eta.y2, exp(lp_sd.y2))

        thetaprop <- theta
        thetaprop$eta.y2 <- eta2

        log.post.eta.y.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.eta.y.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.eta.y.c
        pr.den <- log.post.eta.y.old
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        if (u <= pr & is.nan(pr) == FALSE) {
            theta$eta.y2 <- eta2
        }

        jump[j, 19] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        # y3
        if (j == 1) {
            lp_sd.y3 <- log(proposal$sd.y3)
        }

        eta3 <- rnorm(1, theta$eta.y3, exp(lp_sd.y3))

        thetaprop <- theta
        thetaprop$eta.y3 <- eta3

        log.post.eta.y.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.eta.y.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.eta.y.c
        pr.den <- log.post.eta.y.old
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        if (u <= pr & is.nan(pr) == FALSE) {
            theta$eta.y3 <- eta3
        }

        jump[j, 20] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### delta
        if (j == 1) {
            lp_sd.delta <- log(proposal$sd.delta)
        }

        delta <- rnorm(1, theta$delta, exp(lp_sd.delta))

        thetaprop <- theta
        thetaprop$delta <- delta

        log.post.delta.c <- logpost(thetaprop, theta.prior, complete.data)
        log.post.delta.old <- logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.delta.c
        pr.den <- log.post.delta.old
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        theta$delta <- ifelse((u <= pr & is.nan(pr) == FALSE), delta, theta$delta)
        jump[j, 21] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        thetanop <- theta
        thetanop[[21]] <- NULL
        Theta[j, ] <- unlist(thetanop)

        DID[[j]] <- cbind(complete.data$D, complete.data$ID)

    }  # end for j

    sj <- apply(jump[seq(burn, niter, by = thin), ], 2, mean)
    write.table(sj, paste(wd, "/jumps/jump", scenario, "_", seed, ".txt", sep = ""))
    list(nojump = nojump, jump = jump, Theta = Theta, DID = DID)
}  #End function mcmc.tdiscontinue
