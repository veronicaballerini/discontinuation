da.discontinuation_nocov <- function(theta, theta.prior, observed.data, ID.old, D.old) {

    n <- nrow(observed.data)
    Z <- observed.data$Z
    RY <- observed.data$RY
    Y <- observed.data$Y
    ID.obs <- observed.data$ID.obs
    RD.obs <- observed.data$RD.obs
    Dobs <- observed.data$Dobs

    ID <- D <- rep(NA, n)

    ############# Under treatment

    ## those who die before discontinuing
    ID[Z == 1 & RD.obs == 0 & RY == 1] <- 1  # are ND patients
    D[Z == 1 & RD.obs == 0 & RY == 1] <- 0  # and their discontinuation time is set to 0

    ## those who are observed discontinuing
    ID[Z == 1 & RD.obs == 1] <- 0  # are D patients
    D[Z == 1 & RD.obs == 1] <- Dobs[Z == 1 & RD.obs == 1]  # and their D(1) = Dobs

    ## for those whose discontinuation time and PFS are censored:
    id1 <- (Z == 1 & RD.obs == 0 & RY == 0)

    # 1) compute the probability to be a ND patient as the ratio between num: the survival of a ND patient at censoring weighted by the prob to be ND den: the
    # sum of the num + the survival of a D patient at censoring time weighted by the prob to be D
    num <- theta$pi * Sweib(Y[id1 == 1], theta$alpha.y1nd, theta$beta.y1nd)
    den <- num + (1 - theta$pi) * Sweib(Dobs[id1 == 1], theta$alpha.d, theta$beta.d)

    pr <- num/den

    # 2) draw the stratum membership
    ID[id1 == 1] <- rbinom(sum(id1), 1, pr)

    # 3) for those assigned to the stratum of D patients at step 2), draw D(1)
    ii <- as.numeric(id1 == 1 & ID == 0)
    dii <- rweibull(sum(ii == 1), shape = theta$alpha.d, scale = exp(-(theta$beta.d)/theta$alpha.d))

    # 4) censor D obs if necessary
    D[ii == 1] <- dii * (dii >= Dobs[ii == 1]) + Dobs[ii == 1] * (dii < Dobs[ii == 1])

    rm(id1, num, den, pr)

    ############# Under control
    n0 <- sum(Z == 0)

    ## 1) create a new vector of PS memberships;
    ID.c <- NULL

    # copy the membership for the treated;
    ID.c[Z == 1] <- ID[Z == 1]

    # propose membership for all controls
    ID.c[Z == 0] <- rbinom(n0, 1, theta$pi)

    ## 2) create a new vector of D(1);
    D.c <- rep(0, n)

    # copy D(1) for the treated;
    D.c[Z == 1] <- D[Z == 1]

    # assign D(1) = 0 to the controls that have been proposed to be ND patients
    D.c[Z == 0 & ID.c == 1] <- 0

    # for those controls proposed to be D patients, draw D(1)
    ii <- as.numeric(Z == 0 & ID.c == 0)
    D.c[ii == 1] <- rweibull(sum(ii), shape = theta$alpha.d, scale = exp(-(theta$beta.d)/theta$alpha.d))

    rm(ii)

    ############# Are we going to accept the proposed values?
    complete.data.old <- data.frame(cens = observed.data$f.up, Z = observed.data$Z, ID = ID.old, RD = observed.data$RD.obs, D = D.old, RY = observed.data$RY, Y = observed.data$Y)
    complete.data.c <- data.frame(cens = observed.data$f.up, Z = observed.data$Z, ID = ID.c, RD = observed.data$RD.obs, D = D.c, RY = observed.data$RY, Y = observed.data$Y)

    # 1) compute the ratio between the posterior evaluated at the proposed values and the posterior evaluated at the previous values
    rr <- exp(nc_logpost.i(theta, theta.prior, complete.data.c) - nc_logpost.i(theta, theta.prior, complete.data.old))
    ratio <- rep(0, n)

    # 2) compute rho's as in Algorithm 2, Supplementary material: - for those whose membership is still ND
    ii <- as.numeric(Z == 0 & ID.old == 1 & ID.c == 1)
    ratio[ii == 1] <- rr[ii == 1]

    # - for those who were ND and now their proposed membership is D
    ii <- as.numeric(Z == 0 & ID.old == 1 & ID.c == 0)
    ratio[ii == 1] <- rr[ii == 1] * {
        theta$pi/{
            (1 - theta$pi) * dweib(D.c[ii == 1], theta$alpha.d, theta$beta.d)
        }
    }

    # - for those who were D and now their proposed membership is ND
    ii <- as.numeric(Z == 0 & ID.old == 0 & ID.c == 1)
    ratio[ii == 1] <- rr[ii == 1] * {
        {
            (1 - theta$pi) * dweib(D.old[ii == 1], theta$alpha.d, theta$beta.d)
        }/theta$pi
    }

    # - for those whose membership is still D
    ii <- as.numeric(Z == 0 & ID.old == 0 & ID.c == 0)
    ratio[ii == 1] <- rr[ii == 1] * {
        dweib(D.old[ii == 1], theta$alpha.d, theta$beta.d)/dweib(D.c[ii == 1], theta$alpha.d, theta$beta.d)
    }

    # 3) accept or reject the proposed memberships and D(1)
    ru <- runif(n0)
    ID[Z == 0] <- ID.c[Z == 0] * (ru <= ratio[Z == 0] & is.nan(ratio[Z == 0]) == FALSE) + ID.old[Z == 0] * (ru > ratio[Z == 0] | is.nan(ratio[Z == 0]) == TRUE)
    D[Z == 0] <- D.c[Z == 0] * (ru <= ratio[Z == 0] & is.nan(ratio[Z == 0]) == FALSE) + D.old[Z == 0] * (ru > ratio[Z == 0] | is.nan(ratio[Z == 0]) == TRUE)

    rm(n0, ID.c, D.c, rr, ratio, ru, complete.data.old, complete.data.c)
    list(ID = ID, D = D)

}

