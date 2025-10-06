dce <- function(theta, d, y, complete.data, scenario) {

    ID <- complete.data$ID  # we will need the proportion of ND patients

    # select covariates
    xs <- cbind(1, as.matrix(complete.data[, grepl("x", names(complete.data))]))

    # identify the regressors coefficients 1) regression for PS membership probability
    etas_pos <- grepl("eta.", names(theta)) + grepl("beta", names(theta)) + grepl("eta.d", names(theta)) + grepl("eta.y", names(theta))

    etas <- NULL
    for (i in which(etas_pos == 1)) {
        etas <- rbind(etas, as.numeric(theta[i]))
    }
    etas <- as.matrix(etas)

    pi <- exp(xs %*% etas)/(1 + exp(xs %*% etas))
    pG <- sum(pi)

    # 2) regression for scale parameter of potential outcomes
    etas_y_pos <- grepl("beta.y", names(theta)) + grepl("eta.y", names(theta))

    etas_y <- NULL
    for (i in which(etas_y_pos == 1)) {
        etas_y <- rbind(etas_y, as.numeric(theta[i]))
    }
    etas_y <- as.matrix(etas_y)

    linpred_y <- xs[, -1, drop = FALSE] %*% etas_y

    # 3) regression for scale parameter of D(1)
    etas_d_pos <- grepl("beta.d", names(theta)) + grepl("eta.d", names(theta))

    etas_d <- NULL
    for (i in which(etas_d_pos == 1)) {
        etas_d <- rbind(etas_d, as.numeric(theta[i]))
    }
    etas_d <- as.matrix(etas_d)

    linpred_d <- xs[, -1, drop = FALSE] %*% etas_d
    scale_fd <- exp(-(theta$beta.d + linpred_d)/theta$alpha.d)

    #### DISTRIBUTIONAL CAUSAL EFFECT AND RESTRICTED MEAN SURVIVAL TIME EFFECT FOR ND PATIENTS

    # Treatment
    SY1_ND <- sapply(y, function(h) {
        surv_vals <- Sweib(h, theta$alpha.y1nd, theta$beta.y1nd + linpred_y)
        sum(surv_vals * pi/pG)
    })

    # Control
    SY0_ND <- sapply(y, function(h) {
        surv_vals <- Sweib(h, theta$alpha.y0nd, theta$beta.y0nd + linpred_y)
        sum(surv_vals * pi/pG)
    })

    # Effects
    DCE_ND <- SY1_ND - SY0_ND
    RMSTE_ND <- cumsum(DCE_ND)

    #### DISTRIBUTIONAL CAUSAL EFFECT AND RESTRICTED MEAN SURVIVAL TIME EFFECT FOR D PATIENTS

    # Compute f_d(d|x)
    fd <- sapply(d, function(t) dweibull(t, shape = theta$alpha.d, scale = scale_fd))
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
            surv1 <- Stweib(h, a = theta$alpha.y1d, b = theta$beta.y1d + linpred_y + theta$delta * log_j, l = j)

            SY1_D <- sum(surv1 * fd[, j])/denom_fdj

            # Control use Weibull for scenario I
            if (scenario == 1) {
                surv0 <- Sweib(h, theta$alpha.y0d, theta$beta.y0d + linpred_y + theta$delta * log_j)
            } else {
                surv0 <- Stweib(h, theta$alpha.y0d, theta$beta.y0d + linpred_y + theta$delta * log_j, l = j)
            }

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

