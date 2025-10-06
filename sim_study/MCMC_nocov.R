source("CompleteLogPost_nocov.R")
source("DataAugmentation_nocov.R")

library(mvtnorm)
library(HDInterval)

mcmc.tdiscontinue_nocov <- function(niter, burn, thin, par.start, theta.prior, observed.data, scenario, seed) {

    # restart_counter <- 0 mixing_test_period <- 1000 max_restarts <- 20 # Avoid infinite loop repeat{
    theta.start <- list(eta.0 = rnorm(1, mean = par.start$mu.0, sd = par.start$sigma.0), alpha.d = par.start$alpha.d + abs(rnorm(1, 0, 0.1)), beta.d = par.start$beta.d +
        rnorm(1, 0, 0.1), alpha.y1nd = par.start$shape.y1nd * par.start$scale.y1nd, beta.y1nd = rnorm(1, mean = par.start$mu.y1nd, sd = par.start$sigma.y1nd), alpha.y1d = par.start$shape.y1d *
        par.start$scale.y1d, beta.y1d = rnorm(1, mean = par.start$mu.y1d, sd = par.start$sigma.y1d), alpha.y0nd = par.start$shape.y0nd * par.start$scale.y0nd, beta.y0nd = rnorm(1,
        mean = par.start$mu.y0nd, sd = par.start$sigma.y0nd), alpha.y0d = par.start$shape.y0d * par.start$scale.y0d, beta.y0d = rnorm(1, mean = par.start$mu.y0d,
        sd = par.start$sigma.y0d), delta = rnorm(1, mean = par.start$mu.delta, sd = par.start$sigma.delta))

    pi <- exp(theta.start$eta.0)/(1 + exp(theta.start$eta.0))


    ID.start <- D.start <- NULL
    ID.start[observed.data$Z == 1 & observed.data$RD.obs == 0 & observed.data$RY == 1] <- 1
    ID.start[observed.data$Z == 1 & observed.data$RD.obs == 1] <- 0
    ID.start[observed.data$Z == 1 & observed.data$RD.obs == 0 & observed.data$RY == 0] <- rbinom(sum((observed.data$Z) * (1 - observed.data$RD.obs) * (1 - observed.data$RY)),
        1, pi)
    ID.start[observed.data$Z == 0] <- rbinom(sum(observed.data$Z == 0), 1, pi)

    D.start[ID.start == 1] <- 0
    D.start[ID.start == 0 & observed.data$Z == 1] <- observed.data$Dobs[ID.start == 0 & observed.data$Z == 1]

    D.start[ID.start == 0 & observed.data$Z == 0] <- rweibull(sum(ID.start == 0 & observed.data$Z == 0), shape = theta.start$alpha.d, scale = exp(-(theta.start$beta.d)/theta.start$alpha.d))

    complete.data <- data.frame(cens = observed.data$f.up, Z = observed.data$Z, ID = ID.start, RD = observed.data$RD.obs, D = D.start, RY = observed.data$RY, Y = observed.data$Y)

    theta <- list(eta.0 = theta.start$eta.0, alpha.d = theta.start$alpha.d, beta.d = theta.start$beta.d, alpha.y1nd = theta.start$alpha.y1nd, beta.y1nd = theta.start$beta.y1nd,
        alpha.y1d = theta.start$alpha.y1d, beta.y1d = theta.start$beta.y1d, alpha.y0nd = theta.start$alpha.y0nd, beta.y0nd = theta.start$beta.y0nd, alpha.y0d = theta.start$alpha.y0d,
        beta.y0d = theta.start$beta.y0d, pi = pi, delta = theta.start$delta)

    thetanop <- theta
    thetanop[[12]] <- NULL
    npar <- length(unlist(thetanop))
    jump <- matrix(NA, niter, npar)
    Theta <- matrix(NA, niter, npar)
    DID <- list()

    colnames(jump) <- names(unlist(thetanop))
    colnames(Theta) <- names(unlist(thetanop))

    nojump <- NA

    for (j in 1:niter) {

        # print(j)
        ID.old <- complete.data$ID
        D.old <- complete.data$D
        if (j%%1000 == 0) {
            print(j)
            # print(plot(Theta[1:j,1],type='l'))
        }
        discontinuedata <- da.discontinuation_nocov(theta, theta.prior, observed.data, ID.old, D.old)
        complete.data$ID <- discontinuedata$ID
        complete.data$D <- discontinuedata$D

        # ###eta
        eta.0 <- rnorm(1, theta$eta.0, proposal$sd.0)

        thetaprop <- theta
        thetaprop$eta.0 <- eta.0
        thetaprop$pi <- as.vector(exp(eta.0)/(1 + exp(eta.0)))

        log.post.eta.c <- nc_logpost(thetaprop, theta.prior, complete.data)
        log.post.eta.old <- nc_logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.eta.c
        pr.den <- log.post.eta.old
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        if (u <= pr & is.nan(pr) == FALSE) {
            theta$eta.0 <- eta.0
            theta$pi <- thetaprop$pi
        }

        jump[j, 1] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### alpha.d
        alpha.d <- rgamma(1, shape = theta$alpha.d/proposal$scale.d, scale = proposal$scale.d)

        thetaprop <- theta
        thetaprop$alpha.d <- alpha.d

        log.post.alpha.d.c <- nc_logpost(thetaprop, theta.prior, complete.data)
        log.post.alpha.d.old <- nc_logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.alpha.d.c + log(dgamma(theta$alpha.d, shape = alpha.d/proposal$scale.d, scale = proposal$scale.d))
        pr.den <- log.post.alpha.d.old + log(dgamma(alpha.d, shape = theta$alpha.d/proposal$scale.d, scale = proposal$scale.d))
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        theta$alpha.d <- ifelse((u <= pr & is.nan(pr) == FALSE), alpha.d, theta$alpha.d)
        jump[j, 2] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### beta.d
        beta.d <- rnorm(1, theta$beta.d, proposal$sd.d)

        thetaprop <- theta
        thetaprop$beta.d <- beta.d

        log.post.beta.d.c <- nc_logpost(thetaprop, theta.prior, complete.data)
        log.post.beta.d.old <- nc_logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.beta.d.c
        pr.den <- log.post.beta.d.old

        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        theta$beta.d <- ifelse((u <= pr & is.nan(pr) == FALSE), beta.d, theta$beta.d)
        jump[j, 3] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### alpha.y1nd
        alpha.y1nd <- rgamma(1, shape = theta$alpha.y1nd/proposal$scale.y1nd, scale = proposal$scale.y1nd)

        thetaprop <- theta
        thetaprop$alpha.y1nd <- alpha.y1nd

        log.post.alpha.y1nd.c <- nc_logpost(thetaprop, theta.prior, complete.data)
        log.post.alpha.y1nd.old <- nc_logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.alpha.y1nd.c + log(dgamma(theta$alpha.y1nd, shape = alpha.y1nd/proposal$scale.y1nd, scale = proposal$scale.y1nd))
        pr.den <- log.post.alpha.y1nd.old + log(dgamma(alpha.y1nd, shape = theta$alpha.y1nd/proposal$scale.y1nd, scale = proposal$scale.y1nd))
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        theta$alpha.y1nd <- ifelse((u <= pr & is.nan(pr) == FALSE), alpha.y1nd, theta$alpha.y1nd)
        jump[j, 4] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### beta.y1nd
        beta.y1nd <- rnorm(1, theta$beta.y1nd, proposal$sd.y1nd)

        thetaprop <- theta
        thetaprop$beta.y1nd <- beta.y1nd

        log.post.beta.y1nd.c <- nc_logpost(thetaprop, theta.prior, complete.data)
        log.post.beta.y1nd.old <- nc_logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.beta.y1nd.c
        pr.den <- log.post.beta.y1nd.old
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        theta$beta.y1nd <- ifelse((u <= pr & is.nan(pr) == FALSE), beta.y1nd, theta$beta.y1nd)
        jump[j, 5] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### alpha.y1d
        alpha.y1d <- rgamma(1, shape = theta$alpha.y1d/proposal$scale.y1d, scale = proposal$scale.y1d)

        thetaprop <- theta
        thetaprop$alpha.y1d <- alpha.y1d

        log.post.alpha.y1d.c <- nc_logpost(thetaprop, theta.prior, complete.data)
        log.post.alpha.y1d.old <- nc_logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.alpha.y1d.c + log(dgamma(theta$alpha.y1d, shape = alpha.y1d/proposal$scale.y1d, scale = proposal$scale.y1d))
        pr.den <- log.post.alpha.y1d.old + log(dgamma(alpha.y1d, shape = theta$alpha.y1d/proposal$scale.y1d, scale = proposal$scale.y1d))
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        theta$alpha.y1d <- ifelse((u <= pr & is.nan(pr) == FALSE), alpha.y1d, theta$alpha.y1d)
        jump[j, 6] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### beta.y1d
        beta.y1d <- rnorm(1, theta$beta.y1d, proposal$sd.y1d)

        thetaprop <- theta
        thetaprop$beta.y1d <- beta.y1d

        log.post.beta.y1d.c <- nc_logpost(thetaprop, theta.prior, complete.data)
        log.post.beta.y1d.old <- nc_logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.beta.y1d.c
        pr.den <- log.post.beta.y1d.old

        pr <- exp(pr.num - pr.den)
        u <- runif(1)
        theta$beta.y1d <- ifelse((u <= pr & is.nan(pr) == FALSE), beta.y1d, theta$beta.y1d)
        jump[j, 7] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### alpha.y0nd
        alpha.y0nd <- rgamma(1, shape = theta$alpha.y0nd/proposal$scale.y0nd, scale = proposal$scale.y0nd)

        thetaprop <- theta
        thetaprop$alpha.y0nd <- alpha.y0nd

        log.post.alpha.y0nd.c <- nc_logpost(thetaprop, theta.prior, complete.data)
        log.post.alpha.y0nd.old <- nc_logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.alpha.y0nd.c + log(dgamma(theta$alpha.y0nd, shape = alpha.y0nd/proposal$scale.y0nd, scale = proposal$scale.y0nd))
        pr.den <- log.post.alpha.y0nd.old + log(dgamma(alpha.y0nd, shape = theta$alpha.y0nd/proposal$scale.y0nd, scale = proposal$scale.y0nd))
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        theta$alpha.y0nd <- ifelse((u <= pr & is.nan(pr) == FALSE), alpha.y0nd, theta$alpha.y0nd)
        jump[j, 8] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### beta.y0nd
        beta.y0nd <- rnorm(1, theta$beta.y0nd, proposal$sd.y0nd)

        thetaprop <- theta
        thetaprop$beta.y0nd <- beta.y0nd

        log.post.beta.y0nd.c <- nc_logpost(thetaprop, theta.prior, complete.data)
        log.post.beta.y0nd.old <- nc_logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.beta.y0nd.c
        pr.den <- log.post.beta.y0nd.old

        pr <- exp(pr.num - pr.den)
        u <- runif(1)
        theta$beta.y0nd <- ifelse((u <= pr & is.nan(pr) == FALSE), beta.y0nd, theta$beta.y0nd)
        jump[j, 9] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### alpha.y0d
        alpha.y0d <- rgamma(1, shape = theta$alpha.y0d/proposal$scale.y0d, scale = proposal$scale.y0d)

        thetaprop <- theta
        thetaprop$alpha.y0d <- alpha.y0d

        log.post.alpha.y0d.c <- nc_logpost(thetaprop, theta.prior, complete.data)
        log.post.alpha.y0d.old <- nc_logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.alpha.y0d.c + log(dgamma(theta$alpha.y0d, shape = alpha.y0d/proposal$scale.y0d, scale = proposal$scale.y0d))
        pr.den <- log.post.alpha.y0d.old + log(dgamma(alpha.y0d, shape = theta$alpha.y0d/proposal$scale.y0d, scale = proposal$scale.y0d))
        pr <- exp(pr.num - pr.den)

        u <- runif(1)

        theta$alpha.y0d <- ifelse((u <= pr & is.nan(pr) == FALSE), alpha.y0d, theta$alpha.y0d)
        jump[j, 10] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### beta.y0d
        beta.y0d <- rnorm(1, theta$beta.y0d, proposal$sd.y0d)

        thetaprop <- theta
        thetaprop$beta.y0d <- beta.y0d

        log.post.beta.y0d.c <- nc_logpost(thetaprop, theta.prior, complete.data)
        log.post.beta.y0d.old <- nc_logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.beta.y0d.c
        pr.den <- log.post.beta.y0d.old
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        theta$beta.y0d <- ifelse((u <= pr & is.nan(pr) == FALSE), beta.y0d, theta$beta.y0d)
        jump[j, 11] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        ### delta

        delta <- rnorm(1, theta$delta, proposal$sd.delta)

        thetaprop <- theta
        thetaprop$delta <- delta

        log.post.delta.c <- nc_logpost(thetaprop, theta.prior, complete.data)
        log.post.delta.old <- nc_logpost(theta, theta.prior, complete.data)

        pr.num <- log.post.delta.c
        pr.den <- log.post.delta.old
        pr <- exp(pr.num - pr.den)

        u <- runif(1)
        theta$delta <- ifelse((u <= pr & is.nan(pr) == FALSE), delta, theta$delta)
        jump[j, 12] <- as.numeric(u <= pr & is.nan(pr) == FALSE)

        thetanop <- theta
        thetanop[[12]] <- NULL
        Theta[j, ] <- unlist(thetanop)

        DID[[j]] <- cbind(complete.data$D, complete.data$ID)

        ### [THE KEY MIXING CHECK] if(j %% mixing_test_period == 0 && j > 1){ # check every period jumpblock <- jump[seq_len(j), , drop=FALSE] if(j >= (burn +
        ### mixing_test_period) && j %% mixing_test_period == 0){ # Check acceptance in the *last* mixing_test_period, all after burn-in low <- max(burn+1, j -
        ### mixing_test_period + 1) jumpblock <- jump[low:j, , drop=FALSE] if(all(!is.na(jumpblock))) { accprob <- colMeans(jumpblock)
        ### if(mean(accprob)<0.2|mean(accprob)>0.65) { # If NO parameters have jumped at all, stuck chain: restart restart_counter <- restart_counter + 1
        ### cat('No mixing at iter',j,'- Restarting with new seed:', 474747+restart_counter, '\n') set.seed(474747+restart_counter) if(restart_counter >
        ### max_restarts){ stop('Chain failed to mix after max restarts.') } # break for loop, will restart repeat{} break } } } If last iteration, check for
        ### non-mixing as well: if(j == niter){ jumpblock <- jump[burn:niter, , drop=FALSE] if(all(!is.na(jumpblock))) { accprob <- colMeans(jumpblock)
        ### if(sum(accprob) == 0) { nojump <- 1 cat('No mixing detected in final run (after burn-in).\n') } } }
    }  # end for j

    sj <- apply(jump[seq(burn, niter, by = thin), ], 2, mean)
    write.table(sj, paste(wd, "/jumps/jump_nocov", scenario, "_", seed, ".txt", sep = ""))
    # if(j == niter){ # finished full chain, break repeat{} break } else, jump back to repeat{} and rerun } # end repeat

    list(nojump = nojump, jump = jump, Theta = Theta, DID = DID)
}  #End function mcmc.tdiscontinue
