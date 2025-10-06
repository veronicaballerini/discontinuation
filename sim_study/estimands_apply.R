ce_apply <- function(x) {

    d <- seq(1, 5)  # consider 5 times for the discontinuation
    y <- seq(1, 20)  # consider 20 times for the PFS

    ID <- x$ID  # we will need the proportion of ND patients

    observed.data <- x$observed.data

    # select covariates
    xs <- cbind(1, as.matrix(observed.data[, grepl("x", names(observed.data))]))
    xs <- xs[, -2]

    # identify the regressors coefficients 1) regression for PS membership probability
    etas_pos <- grepl("eta.", names(x)) + grepl("beta", names(x)) + grepl("eta.d", names(x)) + grepl("eta.y", names(x))

    etas <- NULL
    for (i in which(etas_pos == 1)) {
        etas <- rbind(etas, as.numeric(x[i]))
    }
    etas <- as.matrix(etas)

    pi <- exp(xs %*% etas)/(1 + exp(xs %*% etas))
    pG <- sum(pi)

    # 2) regression for scale parameter of potential outcomes
    etas_y_pos <- grepl("beta.y", names(x)) + grepl("eta.y", names(x))

    etas_y <- NULL
    for (i in which(etas_y_pos == 1)) {
        etas_y <- rbind(etas_y, as.numeric(x[i]))
    }
    etas_y <- as.matrix(etas_y)

    linpred_y <- xs[, -1, drop = FALSE] %*% etas_y

    # 3) regression for scale parameter of D(1)
    etas_d_pos <- grepl("beta.d", names(x)) + grepl("eta.d", names(x))

    etas_d <- NULL
    for (i in which(etas_d_pos == 1)) {
        etas_d <- rbind(etas_d, as.numeric(x[i]))
    }
    etas_d <- as.matrix(etas_d)

    linpred_d <- xs[, -1, drop = FALSE] %*% etas_d
    scale_fd <- exp(-(x$beta.d + linpred_d)/x$alpha.d)

    #### DISTRIBUTIONAL CAUSAL EFFECT AND RESTRICTED MEAN SURVIVAL TIME EFFECT FOR ND PATIENTS

    # Treatment
    SY1_ND <- sapply(y, function(h) {
        surv_vals <- Sweib(h, x$alpha.y1nd, x$beta.y1nd + linpred_y)
        sum(surv_vals * pi/pG)
    })

    # Control
    SY0_ND <- sapply(y, function(h) {
        surv_vals <- Sweib(h, x$alpha.y0nd, x$beta.y0nd + linpred_y)
        sum(surv_vals * pi/pG)
    })

    # Effects
    DCE_ND <- SY1_ND - SY0_ND
    RMSTE_ND <- cumsum(DCE_ND)

    #### DISTRIBUTIONAL CAUSAL EFFECT AND RESTRICTED MEAN SURVIVAL TIME EFFECT FOR D PATIENTS

    # Compute f_d(d|x)
    fd <- sapply(d, function(t) dweibull(t, shape = x$alpha.d, scale = scale_fd))
    fd <- matrix(fd, ncol = length(d), byrow = FALSE)

    # Integrate over the empirical distribution of the covariates
    fd_av <- colMeans(fd)

    # Distributional effects:

    DCE_Y_D <- matrix(0, nrow = length(d), ncol = length(y))

    for (j in seq_along(d)) {
        # for each d
        log_j <- log(d[j])

        denom_fdj <- sum(fd[, j])  # sum over units for time d[j]

        for (h in seq_along(y)) {
            # for each y

            # Treatment
            surv1 <- Stweib(h, a = x$alpha.y1d, b = x$beta.y1d + linpred_y + x$delta * log_j, l = j)
            SY1_D <- sum(surv1 * fd[, j])/denom_fdj

            # Control
            surv0 <- Sweib(h, x$alpha.y0d, x$beta.y0d + linpred_y + x$delta * log_j)
            SY0_D <- sum(surv0 * fd[, j])/denom_fdj

            # Effect: DCE_D(y|d)
            DCE_Y_D[j, h] <- SY1_D - SY0_D
        }
    }

    DCE_Y_D[is.nan(DCE_Y_D)] <- 0

    # Effect: DCE_D(y)
    DCE_D <- colSums(DCE_Y_D * matrix(fd_av, nrow = length(d), ncol = length(y), byrow = FALSE))  #weight for f_D(d)

    # Final weighting for the overall effect
    pND <- mean(ID == 1)
    pD <- 1 - pND

    DCE <- DCE_ND * pND + DCE_D * pD  # DCE (ITT)

    # Effect: RMSTE for D and overall (ITT)
    RMSTE_d <- t(apply(DCE_Y_D, 1, cumsum))
    RMSTE_D <- cumsum(DCE_D)

    RMSTE <- pND * RMSTE_ND + pD * RMSTE_D

    # Results:
    list(DCE_ND = DCE_ND, DCE_Y_D = DCE_Y_D, DCE_D = DCE_D, DCE = DCE, RMSTE_ND = RMSTE_ND, RMSTE_D = RMSTE_D, RMSTE_d = RMSTE_d, RMSTE = RMSTE)

}
