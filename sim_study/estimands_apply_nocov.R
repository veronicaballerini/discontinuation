ce_apply_nocov <- function(x) {

    d <- seq(1, 5)  # consider 5 times for the discontinuation
    y <- seq(1, 20)  # consider 20 times for the PFS

    ID <- x$ID  # we will need the proportion of ND patients

    #### DISTRIBUTIONAL CAUSAL EFFECT AND RESTRICTED MEAN SURVIVAL TIME EFFECT FOR ND PATIENTS

    # Treatment
    dce.nd1 <- Sweib(y, x$alpha.y1nd, x$beta.y1nd)

    # Control
    dce.nd0 <- Sweib(y, x$alpha.y0nd, x$beta.y0nd)

    # Effects
    DCE_ND <- dce.nd1 - dce.nd0
    RMSTE_ND <- cumsum(DCE_ND)

    #### DISTRIBUTIONAL CAUSAL EFFECT AND RESTRICTED MEAN SURVIVAL TIME EFFECT FOR D PATIENTS

    DCE_Y_D <- DCE_d <- matrix(0, nrow = length(d), ncol = length(y))

    for (j in seq_along(d)) {
        # for each d
        log_j <- log(d[j])

        for (h in seq_along(y)) {
            # for each y

            # Treatment
            SY1_D <- Stweib(h, a = x$alpha.y1d, b = x$beta.y1d + x$delta * log_j, l = j)

            # Control
            SY0_D <- Sweib(h, x$alpha.y0d, x$beta.y0d + x$delta * log_j)

            # Effect: DCE_D(y|d)
            DCE_Y_D[j, h] <- ifelse(is.nan(DCE_Y_D[j, h]) == TRUE, 0, SY1_D - SY0_D)

        }

        DCE_d[j, ] <- DCE_Y_D[j, ] * dweib(j, a = x$alpha.d, b = x$beta.d)  #weight for f_D(d)

    }

    # Effect: DCE_D(y)
    DCE_D <- apply(DCE_d, 2, sum)

    # Final weighting for the overall effect
    pND <- mean(ID == 1)
    pD <- 1 - pND

    DCE <- DCE_ND * pND + DCE_D * pD  # DCE (ITT)

    # Effect: RMSTE for D and overall (ITT)
    RMSTE_d <- t(apply(DCE_Y_D, 1, cumsum))
    RMSTE_D <- cumsum(DCE_D)

    RMSTE <- pND * RMSTE_ND + pD * RMSTE_D

    list(DCE_ND = DCE_ND, DCE_Y_D = DCE_Y_D, DCE_D = DCE_D, DCE = DCE, RMSTE_ND = RMSTE_ND, RMSTE_D = RMSTE_D, RMSTE_d = RMSTE_d, RMSTE = RMSTE)
}
