rmste_apply_expexp <- function(x) {

    data <- x$data
    d <- seq(1, 5)
    y <- seq(1, 20)
    n_obs <- nrow(data)
    ID <- x$ID

    xs <- cbind(1, as.matrix(data[, grepl("x", names(data))]))

    etas_pos <- grepl("eta.", names(x)) + grepl("beta", names(x)) + grepl("eta.d",
        names(x)) + grepl("eta.y", names(x))

    etas <- NULL
    for (i in which(etas_pos == 1)) {
        etas <- rbind(etas, as.numeric(x[i]))
    }
    etas <- as.matrix(etas)

    pi <- exp(xs %*% etas)/(1 + exp(xs %*% etas))

    etas_y_pos <- grepl("beta.y", names(x)) + grepl("eta.y", names(x))

    etas_y <- NULL
    for (i in which(etas_y_pos == 1)) {
        etas_y <- rbind(etas_y, as.numeric(x[i]))
    }
    etas_y <- as.matrix(etas_y)

    etas_d_pos <- grepl("beta.d", names(x)) + grepl("eta.d", names(x))

    etas_d <- NULL
    for (i in which(etas_d_pos == 1)) {
        etas_d <- rbind(etas_d, as.numeric(x[i]))
    }
    etas_d <- as.matrix(etas_d)

    pG <- sum(pi)

    # DCE for ND

    linpred_y <- xs[, -1, drop = FALSE] %*% etas_y

    SY1_ND <- sapply(y, function(h) {
        surv_vals <- Sexp(h, exp(x$beta.y1nd + linpred_y))
        sum(surv_vals * pi/pG)
    })

    SY0_ND <- sapply(y, function(h) {
        surv_vals <- Sexp(h, exp(x$beta.y0nd + linpred_y))
        sum(surv_vals * pi/pG)
    })

    DCE_ND <- SY1_ND - SY0_ND
    RMSTE_ND <- cumsum(DCE_ND)

    # DCE for D

    linpred_d <- xs[, -1, drop = FALSE] %*% etas_d
    rate_fd <- exp(x$beta.d + linpred_d)

    fd <- matrix(0, nrow = n_obs, ncol = length(d))
    for (j in seq_along(d)) {
        fd[, j] <- dexp(d[j], rate = rate_fd)
    }

    fD_av <- colMeans(fd)

    # DCE_Y_D

    DCE_Y_D <- matrix(0, nrow = length(d), ncol = length(y))

    for (j in seq_along(d)) {
        log_j <- log(d[j])

        denom_fdj <- sum(fd[, j])

        for (h in seq_along(y)) {
            surv1 <- Stexp(h, rate = exp(x$beta.y1d + linpred_y + x$delta *
                log_j), a = j)
            SY1_D <- sum(surv1 * fd[, j])/denom_fdj

            surv0 <- Sexp(h, exp(x$beta.y0d + linpred_y + x$delta * log_j))
            SY0_D <- sum(surv0 * fd[, j])/denom_fdj

            DCE_Y_D[j, h] <- SY1_D - SY0_D
        }
    }

    DCE_Y_D[is.nan(DCE_Y_D)] <- 0

    # Average across distributions fD_av
    DCE_D <- colSums(DCE_Y_D * matrix(fD_av, nrow = length(d), ncol = length(y),
        byrow = FALSE))

    # Final weighting
    pND <- mean(ID == 1)
    pD <- 1 - pND

    DCE <- DCE_ND * pND + DCE_D * pD

    # Cumulative estimates
    RMSTE_d <- t(apply(DCE_Y_D, 1, cumsum))
    RMSTE_D <- cumsum(DCE_D)

    RMSTE <- pND * RMSTE_ND + pD * RMSTE_D

    return(list(DCE_ND = DCE_ND, DCE_Y_D = DCE_Y_D, DCE_D = DCE_D, DCE = DCE,
        RMSTE_ND = RMSTE_ND, RMSTE_D = RMSTE_D, RMSTE_d = RMSTE_d, RMSTE = RMSTE))
}


bicfun_expexp <- function(x) {

    data <- x$data
    n <- nrow(data)
    k <- length(x) - 3

    Z <- data$Z
    ID <- x$ID
    RD <- data$RD
    D <- x$D
    RY <- data$RY
    Y <- data$Y
    x1 <- data$x1
    x2 <- data$x2
    x3 <- data$x3

    xs <- cbind(1, as.matrix(data[, grepl("x", names(data))]))

    etas_pos <- grepl("eta.", names(x)) + grepl("beta", names(x)) + grepl("eta.d",
        names(x)) + grepl("eta.y", names(x))

    etas <- NULL
    for (i in which(etas_pos == 1)) {
        etas <- rbind(etas, as.numeric(x[i]))
    }
    etas <- as.matrix(etas)

    pi <- exp(xs %*% etas)/(1 + exp(xs %*% etas))

    loglik <- 0
    # Z=1 ID=1 RY=1 ---> Never discontinued who die under treatment
    ii <- Z * ID * RY == 1
    loglik <- loglik + sum({
        log(pi[ii]) + log(dexp(Y[ii], rate = exp(x$beta.y1nd + x$eta.y1 *
            x1[ii] + x$eta.y2 * x2[ii] + x$eta.y3 * x3[ii])))
    }, na.rm = TRUE)

    # Z=1 ID=1 RY=0 ---> Never discontinued who do not die under
    # treatment
    ii <- Z * ID * (1 - RY) == 1
    loglik <- loglik + sum({
        log(pi[ii]) + log(Sexp(Y[ii], exp(x$beta.y1nd + x$eta.y1 * x1[ii] +
            x$eta.y2 * x2[ii] + x$eta.y3 * x3[ii])))
    }, na.rm = TRUE)

    # Z=1 ID=0 RD=1 RY=1 ---> Discontinued who die under treatment for
    # whom we observe the discontinuation
    ii <- Z * (1 - ID) * RD * RY == 1
    loglik <- loglik + sum({
        log(1 - pi[ii]) + log(dexp(D[ii], exp(x$beta.d + x$eta.d1 * x1[ii] +
            x$eta.d2 * x2[ii] + x$eta.d3 * x3[ii]))) + log(dtexp(Y[ii], rate = exp(x$beta.y1d +
            x$eta.y1 * x1[ii] + x$eta.y2 * x2[ii] + x$eta.y3 * x3[ii] + x$delta *
            log(D[ii])), a = D[ii]))
    }, na.rm = TRUE)

    # Z=1 ID=0 RD=1 RY=0 ---> Discontinued who do not die under
    # treatment for whom we observe the discontinuation
    ii <- Z * (1 - ID) * RD * (1 - RY) == 1
    loglik <- loglik + sum({
        log(1 - pi[ii]) + log(dexp(D[ii], exp(x$beta.d + x$eta.d1 * x1[ii] +
            x$eta.d2 * x2[ii] + x$eta.d3 * x3[ii]))) + log(Stexp(Y[ii], rate = exp(x$beta.y1d +
            x$eta.y1 * x1[ii] + x$eta.y2 * x2[ii] + x$eta.y3 * x3[ii] + x$delta *
            log(D[ii])), a = D[ii]))
    }, na.rm = TRUE)

    # Z=1 ID=0 RD=0 RY=0 ---> Discontinued who discontinue after the
    # end of the study and do not die under treatment
    ii <- Z * (1 - ID) * (1 - RD) * (1 - RY) == 1
    loglik <- loglik + sum({
        log(1 - pi[ii]) + log(Sexp(D[ii], exp(x$beta.d + x$eta.d1 * x1[ii] +
            x$eta.d2 * x2[ii] + x$eta.d3 * x3[ii])))
    }, na.rm = TRUE)


    # Z=0 ID=1 RY=1 ---> Never discontinued who die under control
    ii <- (1 - Z) * ID * RY == 1
    loglik <- loglik + sum({
        log(pi[ii]) + log(dexp(Y[ii], rate = exp(x$beta.y0nd + x$eta.y1 *
            x1[ii] + x$eta.y2 * x2[ii] + x$eta.y3 * x3[ii])))
    }, na.rm = TRUE)

    # Z=0 ID=1 RY=0 ---> Never discontinued who do not die under
    # control
    ii <- (1 - Z) * ID * (1 - RY) == 1
    loglik <- loglik + sum({
        log(pi[ii]) + log(Sexp(Y[ii], exp(x$beta.y0nd + x$eta.y1 * x1[ii] +
            x$eta.y2 * x2[ii] + x$eta.y3 * x3[ii])))
    }, na.rm = TRUE)

    # Z=0 ID=0 RY=1 ---> Discontinued who die under control
    ii <- (1 - Z) * (1 - ID) * RY == 1
    loglik <- loglik + sum({
        log(1 - pi[ii]) + log(dexp(D[ii], exp(x$beta.d + x$eta.d1 * x1[ii] +
            x$eta.d2 * x2[ii] + x$eta.d3 * x3[ii]))) + log(dexp(Y[ii], exp(x$beta.y0d +
            x$delta * log(D[ii]) + x$eta.y1 * x1[ii] + x$eta.y2 * x2[ii] +
            x$eta.y3 * x3[ii])))
    }, na.rm = TRUE)

    # Z=0 ID=0 RY=0 ---> Discontinued who do not die under control
    ii <- (1 - Z) * (1 - ID) * (1 - RY) == 1
    loglik <- loglik + sum({
        log(1 - pi[ii]) + log(dexp(D[ii], exp(x$beta.d + x$eta.d1 * x1[ii] +
            x$eta.d2 * x2[ii] + x$eta.d3 * x3[ii]))) + log(Sexp(Y[ii], exp(x$beta.y0d +
            x$delta * log(D[ii]) + x$eta.y1 * x1[ii] + x$eta.y2 * x2[ii] +
            x$eta.y3 * x3[ii])))
    }, na.rm = TRUE)

    bic <- -2 * loglik + k * log(n)
    return(bic = bic)
}

waicfun_expexp <- function(x) {

    data <- x$data
    n <- nrow(data)

    Z <- data$Z
    ID <- x$ID
    RD <- data$RD
    D <- x$D
    RY <- data$RY
    Y <- data$Y
    x1 <- data$x1
    x2 <- data$x2
    x3 <- data$x3

    xs <- cbind(1, as.matrix(data[, grepl("x", names(data))]))

    etas_pos <- grepl("eta.", names(x)) + grepl("beta", names(x)) + grepl("eta.d",
        names(x)) + grepl("eta.y", names(x))

    etas <- NULL
    for (i in which(etas_pos == 1)) {
        etas <- rbind(etas, as.numeric(x[i]))
    }
    etas <- as.matrix(etas)

    pi <- exp(xs %*% etas)/(1 + exp(xs %*% etas))

    loglik <- rep(0, n)
    # Z=1 ID=1 RY=1 ---> Never discontinued who die under treatment
    ii <- Z * ID * RY == 1
    loglik[ii] <- log(pi[ii]) + log(dexp(Y[ii], exp(x$beta.y1nd + x$eta.y1 *
        x1[ii] + x$eta.y2 * x2[ii] + x$eta.y3 * x3[ii])))

    # Z=1 ID=1 RY=0 ---> Never discontinued who do not die under
    # treatment
    ii <- Z * ID * (1 - RY) == 1
    loglik[ii] <- log(pi[ii]) + log(Sexp(Y[ii], exp(x$beta.y1nd + x$eta.y1 *
        x1[ii] + x$eta.y2 * x2[ii] + x$eta.y3 * x3[ii])))

    # Z=1 ID=0 RD=1 RY=1 ---> Discontinued who die under treatment for
    # whom we observe the discontinuation
    ii <- Z * (1 - ID) * RD * RY == 1
    loglik[ii] <- log(1 - pi[ii]) + log(dexp(D[ii], exp(x$beta.d + x$eta.d1 *
        x1[ii] + x$eta.d2 * x2[ii] + x$eta.d3 * x3[ii]))) + log(dtexp(Y[ii],
        exp(x$beta.y1d + x$eta.y1 * x1[ii] + x$eta.y2 * x2[ii] + x$eta.y3 *
            x3[ii] + x$delta * log(D[ii])), a = D[ii]))


    # Z=1 ID=0 RD=1 RY=0 ---> Discontinued who do not die under
    # treatment for whom we observe the discontinuation
    ii <- Z * (1 - ID) * RD * (1 - RY) == 1
    loglik[ii] <- log(1 - pi[ii]) + log(dexp(D[ii], exp(x$beta.d + x$eta.d1 *
        x1[ii] + x$eta.d2 * x2[ii] + x$eta.d3 * x3[ii]))) + log(Stexp(Y[ii],
        exp(x$beta.y1d + x$eta.y1 * x1[ii] + x$eta.y2 * x2[ii] + x$eta.y3 *
            x3[ii] + x$delta * log(D[ii])), a = D[ii]))


    # Z=1 ID=0 RD=0 RY=0 ---> Discontinued who discontinue after the
    # end of the study and do not die under treatment
    ii <- Z * (1 - ID) * (1 - RD) * (1 - RY) == 1
    loglik[ii] <- log(1 - pi[ii]) + log(Sexp(D[ii], exp(x$beta.d + x$eta.d1 *
        x1[ii] + x$eta.d2 * x2[ii] + x$eta.d3 * x3[ii])))


    # Z=0 ID=1 RY=1 ---> Never discontinued who die under control
    ii <- (1 - Z) * ID * RY == 1
    loglik[ii] <- log(pi[ii]) + log(dexp(Y[ii], exp(x$beta.y0nd + x$eta.y1 *
        x1[ii] + x$eta.y2 * x2[ii] + x$eta.y3 * x3[ii])))

    # Z=0 ID=1 RY=0 Y>Y1*k ---> Never discontinued who do not die under
    # control
    ii <- (1 - Z) * ID * (1 - RY) == 1
    loglik[ii] <- log(pi[ii]) + log(Sexp(Y[ii], exp(x$beta.y0nd + x$eta.y1 *
        x1[ii] + x$eta.y2 * x2[ii] + x$eta.y3 * x3[ii])))

    # Z=0 ID=0 RY=1 ---> Discontinued who die under control
    ii <- (1 - Z) * (1 - ID) * RY == 1
    loglik[ii] <- log(1 - pi[ii]) + log(dexp(D[ii], exp(x$beta.d + x$eta.d1 *
        x1[ii] + x$eta.d2 * x2[ii] + x$eta.d3 * x3[ii]))) + log(dexp(Y[ii],
        exp(x$beta.y0d + x$delta * log(D[ii]) + x$eta.y1 * x1[ii] + x$eta.y2 *
            x2[ii] + x$eta.y3 * x3[ii])))

    # Z=0 ID=0 RY=0 ---> Discontinued who do not die under control
    ii <- (1 - Z) * (1 - ID) * (1 - RY) == 1
    loglik[ii] <- log(1 - pi[ii]) + log(dexp(D[ii], exp(x$beta.d + x$eta.d1 *
        x1[ii] + x$eta.d2 * x2[ii] + x$eta.d3 * x3[ii]))) + log(Sexp(Y[ii],
        exp(x$beta.y0d + x$delta * log(D[ii]) + x$eta.y1 * x1[ii] + x$eta.y2 *
            x2[ii] + x$eta.y3 * x3[ii])))

    return(loglik = loglik)
}


kaplanMeier_expexp <- function(x) {

    data <- x$data
    y <- seq(1, 20)
    d <- seq(1, 20)

    nd <- x$ID == 1
    fit.y.nd <- survfit(Surv(data$Y[nd], data$RY[nd]) ~ data$Z[nd])

    d <- x$ID == 0
    fit.y.d <- survfit(Surv(data$Y[d], data$RY[d]) ~ data$Z[d])

    d1 <- data$Z == 1 & x$ID == 0
    fit.d <- survfit(Surv(data$Dobs[d1], data$RD.obs[d1]) ~ 1)

    G0.nd <- G1.nd <- G0.d <- G1.d <- rep(0, length(y))
    GD.d1 <- rep(0, length(y))

    for (j in 1:length(y)) {
        if (all(fit.y.nd[1]$time > y[j])) {
            G0.nd[j] <- 1
        } else {
            i0 <- max(which(fit.y.nd[1]$time <= y[j]))
            G0.nd[j] <- fit.y.nd[1]$surv[i0]
        }
        if (all(fit.y.nd[2]$time > y[j])) {
            G1.nd[j] <- 1
        } else {
            i1 <- max(which(fit.y.nd[2]$time <= y[j]))
            G1.nd[j] <- fit.y.nd[2]$surv[i1]
        }
        if (all(fit.y.d[1]$time > y[j])) {
            G1.d[j] <- 1
        } else {
            i0 <- max(which(fit.y.d[1]$time <= y[j]))
            G0.d[j] <- fit.y.d[1]$surv[i0]
        }
        if (all(fit.y.d[2]$time > y[j])) {
            G1.d[j] <- 1
        } else {
            i1 <- max(which(fit.y.d[2]$time <= y[j]))
            G1.d[j] <- fit.y.d[2]$surv[i1]
        }
        if (all(fit.d$time > d[j])) {
            GD.d1[j] <- 1
        } else {
            i1 <- max(which(fit.d$time <= d[j]))
            GD.d1[j] <- fit.d$surv[i1]
        }
    }

    Diff.G.nd <- abs(G1.nd - G0.nd)
    Diff.G.d <- abs(G1.d - G0.d)

    list(Diff.G.nd = Diff.G.nd, Diff.G.d = Diff.G.d, GD.d1 = GD.d1)
}
