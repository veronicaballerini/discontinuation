rep.data <- function(x) {

    data <- x$data

    th <- 34 - data$entry
    n <- nrow(data)

    Z.rep <- data$Z

    # select covariates
    xs <- cbind(1, as.matrix(data[, grepl("x", names(data))]))

    # identify the regressors coefficients 1)
    # regression for PS membership probability
    etas_pos <- grepl("eta.", names(x)) + grepl("beta",
        names(x)) + grepl("eta.d", names(x)) + grepl("eta.y",
        names(x))

    etas <- NULL
    for (i in which(etas_pos == 1)) {
        etas <- rbind(etas, as.numeric(x[i]))
    }
    etas <- as.matrix(etas)

    pi <- exp(xs %*% etas)/(1 + exp(xs %*% etas))
    pG <- sum(pi)

    # 2) regression for scale parameter of
    # potential outcomes
    etas_y_pos <- grepl("beta.y", names(x)) + grepl("eta.y",
        names(x))

    etas_y <- NULL
    for (i in which(etas_y_pos == 1)) {
        etas_y <- rbind(etas_y, as.numeric(x[i]))
    }
    etas_y <- as.matrix(etas_y)

    linpred_y <- xs[, -1, drop = FALSE] %*% etas_y

    # 3) regression for scale parameter of
    # D(1)
    etas_d_pos <- grepl("beta.d", names(x)) + grepl("eta.d",
        names(x))

    etas_d <- NULL
    for (i in which(etas_d_pos == 1)) {
        etas_d <- rbind(etas_d, as.numeric(x[i]))
    }
    etas_d <- as.matrix(etas_d)

    linpred_d <- xs[, -1, drop = FALSE] %*% etas_d

    ## ND binary indicator
    ID.rep <- rbinom(n, 1, pi)

    ## Time to discontinuation for D patients
    D1.rep <- rep(0, n)
    D1.rep[ID.rep == 0] <- rexp(sum(ID.rep == 0),
        rate = exp(x$beta.d + linpred_d))

    ## PFS
    Y.rep <- NULL
    # - for those assigned to control that are
    # ND patients
    Y.rep[Z.rep == 0 & ID.rep == 1] <- rexp(sum(Z.rep ==
        0 & ID.rep == 1), rate = exp(x$beta.y0nd +
        linpred_y[Z.rep == 0 & ID.rep == 1]))

    # - for those assigned to control that are
    # D patients
    Y.rep[Z.rep == 0 & ID.rep == 0] <- rexp(sum(Z.rep ==
        0 & ID.rep == 0), rate = exp(x$beta.y0d +
        linpred_y[Z.rep == 0 & ID.rep == 0] + x$delta *
        log(D1.rep[Z.rep == 0 & ID.rep == 0])))

    # - for those assigned to treatment that
    # are ND patients
    Y.rep[Z.rep == 1 & ID.rep == 1] <- rexp(sum(Z.rep ==
        1 & ID.rep == 1), rate = exp(x$beta.y1nd +
        linpred_y[Z.rep == 1 & ID.rep == 1]))

    # - for those assigned to treatment that
    # are D patients
    Y.rep[Z.rep == 1 & ID.rep == 0] <- rtexp(sum(Z.rep ==
        1 & ID.rep == 0), rate = exp(x$beta.y1d +
        linpred_y[Z.rep == 1 & ID.rep == 0] + x$delta *
        log(D1.rep[Z.rep == 1 & ID.rep == 0])),
        a = D1.rep[Z.rep == 1 & ID.rep == 0])

    # censoring
    RY.rep <- rep(1, n)
    RY.rep[Y.rep >= th] <- 0

    Y.rep[RY.rep == 0] <- th[RY.rep == 0]

    RD.rep <- rep(0, n)
    RD.rep[Z.rep == 0 & ID.rep == 0 & D1.rep <=
        th] <- 1

    D1.rep[RD.rep == 0] <- th[RD.rep == 0]

    complete.data.rep <- data.frame(entry = data$entry,
        Z = Z.rep, RD.obs = RD.rep, Dobs = D1.rep,
        RY = RY.rep, Y = Y.rep, ID = ID.rep, x1 = data$x1,
        x2 = data$x2, x3 = data$x3)
    return(list(eta.0 = x$eta.0, eta.1 = x$eta.1,
        eta.2 = x$eta.2, eta.3 = x$eta.3, beta.d = x$beta.d,
        eta.d1 = x$eta.d1, eta.d2 = x$eta.d2, eta.d3 = x$eta.d3,
        beta.y1nd = x$beta.y1nd, beta.y1d = x$beta.y1d,
        beta.y0nd = x$beta.y0nd, beta.y0d = x$beta.y0d,
        eta.y1 = x$eta.y1, eta.y2 = x$eta.y2, eta.y3 = x$eta.y3,
        delta = x$delta, D = D1.rep, ID = ID.rep,
        data = complete.data.rep))

}
