applyfun_1 <- function(x) {

    # set the seed
    seed <- x

    # if this seed has already generated a dataset, use that dataset
    if (file.exists(paste(wd, "/sim/simdata1_", seed, ".txt", sep = ""))) {
        set.seed(seed)
    } else {
        # otherwise, create a new dataset
        set.seed(seed)

        ############################################################################## DATA SIMULATION ###############################

        # Pre-treatment variables: X1 = Continuous variable; the bigger the number, the higher the risk X2 = Presence of specific metastatic status before
        # enrolling the study X3 = Disease burden indicator. Treatment variable: Z=z, z=0,1 Discontinuation indicator: ID = 1 if the discontinuation cannot be
        # observed (due to death or progression) ID = 0 if the discontinuation is observed, even after the follow-up period Discontinuation variable: D1 =
        # Continuous variable; time until the discontinuation under active treatment Outcome variable: Y = survival

        n <- 389  #total sample size
        nT <- 200  #number of treated individuals
        nC <- n - nT  #number of controls

        x1 <- rnorm(n)  #risk
        x1st <- (x1 - mean(x1))/sd(x1)  #risk score
        x2 <- rbinom(n, 1, 0.491)  #1 if metastatic status prior to the enrollment, 0 otherwise
        x3 <- rbinom(n, 1, 0.319)  #1 if the disease burden is >=3, 0 if <3

        Z <- matrix(c(rep(1, nT), rep(0, nC)), n, 1)  #treatment

        ############################### Discontinuation ################################ Indicator I{D_i(1)=D} = 1 if the potential discontinuation under
        ############################### treatment does not apply [Never Discontinued] = 0 otherwise I{D_i(1)=D} ~ Bernoulli(pix_i)

        eta <- list(int = 0.4, x1st = -0.5, x2 = 0.45, x3 = 0.55)  #intercept:  #previously x1st=.3
        # we assumed the proportion of ND to be approximately equal to .6: the lower bound is 93/200 and the upper bound is 136/200 eta$x1st: we assumed the
        # older the patient the lower the probability to never discontinue, since discontinuation is competing with progression; eta$x2,x3: we assumed the more
        # severe the condition, the higher the probability to never discontinue (as above)

        pix <- exp(eta$int + eta$x1 * x1st + eta$x2 * x2 + eta$x3 * x3)/(1 + exp(eta$int + eta$x1st * x1st + eta$x2 * x2 + eta$x3 * x3))
        ID <- rbinom(n, 1, pix)

        # Discontinuation variable: D1 ~ Weibull(shape=alpha.d, scale=beta.d+eta.d1*x1st+eta.d2*x2+eta.d3*x3)

        # parameters of the Weibulls:
        shape.d <- alpha.d <- 1.75

        median.d1 <- 2.2  #time to adverse event, median
        par.scale.d <- list(beta.d = -log((median.d1^alpha.d)/log(2)), eta.d1 = 0.45, eta.d2 = 0.25, eta.d3 = 0.35)
        # eta.d1,d2,d3 = we assumed that individuals with higher risk discontinue sooner, given that they will discontinue
        scale.d <- exp(-(par.scale.d$beta.d + par.scale.d$eta.d1 * x1st + par.scale.d$eta.d2 * x2 + par.scale.d$eta.d3 * x3)/alpha.d)

        D1 <- rep(NA, n)  #potential discontinuation under treatment
        D1[ID == 0] <- rweibull(sum(1 - ID), shape = shape.d, scale = scale.d[ID == 0])  #potential discontinuation under
        # treatment for discontinuers
        tc <- sample(0:1, n, replace = TRUE)
        entry <- tc * rbinom(n, 23, 0.9) + (1 - tc) * sample(0:23, n, replace = TRUE)
        f.up <- 33 - entry  #follow up period, from entry to study end (33 months)

        # Outcome variable: Y1|D1=NA ~ Weibull(shape = alpha.y1nd, scale = beta.y1nd + eta.y1*x1st + eta.y2*x2 + eta.y3*x3) Y1|D1=d ~ Weibull(shape =
        # alpha.y1d, scale = beta.y1d + delta*log(D1) + eta.y1*x1st + eta.y2*x2 + eta.y3*x3) Y0|D1=NA ~ Weibull(shape = alpha.y0nd, scale = beta.y0nd +
        # eta.y1*x1st + eta.y2*x2 + eta.y3*x3) Y0|D1=d ~ Weibull(shape=alpha.y0d, scale = beta.y0d + delta*log(D1) + eta.y1*x1st + eta.y2*x2 + eta.y3*x3)

        ###################### parameters under H0 of no effect ######################## parameters of the Weibull:
        shape.y1nd <- alpha.y1nd <- 1.2
        shape.y0nd <- alpha.y0nd <- 1.1
        shape.y1d <- alpha.y1d <- 1.2
        shape.y0d <- alpha.y0d <- 1.1

        mean.y1nd <- 14  #time to primary event, average
        mean.y0nd <- 11  #time to primary event, average

        mean.y1d <- 15  #time to primary event, average
        mean.y0d <- 10  #time to primary event, average

        beta.y1nd <- -alpha.y1nd * log(mean.y1nd/gamma(1 + 1/alpha.y1nd))
        beta.y0nd <- -alpha.y0nd * log(mean.y0nd/gamma(1 + 1/alpha.y0nd))
        beta.y1d <- -alpha.y1d * log(mean.y1d/gamma(1 + 1/alpha.y1d))
        beta.y0d <- -alpha.y0d * log(mean.y0d/gamma(1 + 1/alpha.y0d))

        # eta.y1,y2,y3 = we assumed that individuals with higher risk die sooner, given that they will discontinue
        par.scale.y1nd <- list(beta.y = beta.y1nd, eta.y1 = 0.25, eta.y2 = 0.7, eta.y3 = 0.45)
        par.scale.y0nd <- list(beta.y = beta.y0nd, eta.y1 = 0.25, eta.y2 = 0.7, eta.y3 = 0.45)
        par.scale.y1d <- list(beta.y = beta.y1d, eta.y1 = 0.25, eta.y2 = 0.7, eta.y3 = 0.45)
        par.scale.y0d <- list(beta.y = beta.y0d, eta.y1 = 0.25, eta.y2 = 0.7, eta.y3 = 0.45)

        theta <- list(alpha.y1nd = alpha.y1nd, beta.y1nd = beta.y1nd, eta.y1 = par.scale.y1nd$eta.y1, eta.y2 = par.scale.y1nd$eta.y2, eta.y3 = par.scale.y1nd$eta.y3,
            alpha.y1d = alpha.y1d, beta.y1d = beta.y1d, alpha.y0nd = alpha.y0nd, beta.y0nd = beta.y0nd, alpha.y0d = alpha.y0d, beta.y0d = beta.y0d, delta = 0.3)

        scale.y1nd <- exp(-(theta$beta.y1nd + theta$eta.y1 * x1st + theta$eta.y2 * x2 + theta$eta.y3 * x3)/theta$alpha.y1nd)
        scale.y1d <- exp(-(theta$beta.y1d + theta$delta * log(D1) + theta$eta.y1 * x1st + theta$eta.y2 * x2 + theta$eta.y3 * x3)/theta$alpha.y1d)
        scale.y0nd <- exp(-(theta$beta.y0nd + theta$eta.y1 * x1st + theta$eta.y2 * x2 + theta$eta.y3 * x3)/theta$alpha.y0nd)
        scale.y0d <- exp(-(theta$beta.y0d + theta$delta * log(D1) + theta$eta.y1 * x1st + theta$eta.y2 * x2 + theta$eta.y3 * x3)/theta$alpha.y0d)

        Y1 <- NULL  #potential outcome under treatment
        Y0 <- NULL  #potential outcome under control

        # for ND:
        Y1[ID == 1] <- rweibull(sum(ID), shape = theta$alpha.y1nd, scale = scale.y1nd[ID == 1])
        Y0[ID == 1] <- rweibull(sum(ID), shape = theta$alpha.y0nd, scale = scale.y0nd[ID == 1])

        Y1[ID == 0] <- rtweibull(sum(1 - ID), shape = theta$alpha.y1d, scale = scale.y1d[ID == 0], a = D1[ID == 0])

        if (sum(is.infinite(Y1) == TRUE) > 0) {
            Y1[which(is.infinite(Y1) == TRUE)] <- D1[which(is.infinite(Y1) == TRUE)]
            D1[which(is.infinite(Y1) == TRUE)] <- NA
            ID[which(is.infinite(Y1) == TRUE)] <- 1
        }

        Y0[ID == 0] <- rweibull(sum(1 - ID), shape = theta$alpha.y0d, scale = scale.y0d[ID == 0])

        if (sum(is.infinite(Y0) == TRUE) > 0) {
            Y0[which(is.infinite(Y0) == TRUE)] <- D1[which(is.infinite(Y0) == TRUE)]
            D1[which(is.infinite(Y0) == TRUE)] <- NA
            ID[which(is.infinite(Y0) == TRUE)] <- 1
        }

        ######## observed data ########
        Dobs <- RD <- Dobs <- RY <- Yobs <- rep(0, n)

        Dobs[Z == 1 & ID == 0] <- D1[Z == 1 & ID == 0] * (D1[Z == 1 & ID == 0] < f.up[Z == 1 & ID == 0]) + f.up[Z == 1 & ID == 0] * (D1[Z == 1 & ID == 0] >= f.up[Z ==
            1 & ID == 0])
        Dobs[Z == 1 & ID == 1] <- f.up[Z == 1 & ID == 1]
        Yobs[Z == 0] <- Y0[Z == 0] * (Y0[Z == 0] < f.up[Z == 0]) + f.up[Z == 0] * (Y0[Z == 0] >= f.up[Z == 0])
        Yobs[Z == 1] <- Y1[Z == 1] * (Y1[Z == 1] < f.up[Z == 1]) + f.up[Z == 1] * (Y1[Z == 1] >= f.up[Z == 1])

        RD[Z == 1 & ID == 0] <- 1 * (D1[Z == 1 & ID == 0] < f.up[Z == 1 & ID == 0])
        # indicator variable: 1 if the discontinuation is observed
        RY[Z == 0] <- 1 * (Y0[Z == 0] <= f.up[Z == 0])  # indicator variables: 1 if the
        RY[Z == 1] <- 1 * (Y1[Z == 1] <= f.up[Z == 1])  # death is observed

        ID.obs <- rep(0, n)  #0= Discontinued/Not available, 1=Never discontinued
        ID.obs[Z == 1 & RD == 0 & RY == 1] <- 1
        ID.obs[Z == 1 & RD == 1] <- 0

        RD.obs <- RD
        Y <- Yobs
        cens <- rep(33, n)

        data <- data.frame(entry, Z, f.up, cens, RD.obs, Dobs, RY, Y, ID.obs, x1, x1st, x2, x3)

        thetatrue <- list(eta.0 = eta$int, eta.1 = eta$x1st, eta.2 = eta$x2, eta.3 = eta$x3, alpha.d = alpha.d, beta.d = par.scale.d$beta.d, eta.d1 = par.scale.d$eta.d1,
            eta.d2 = par.scale.d$eta.d2, eta.d3 = par.scale.d$eta.d3, alpha.y1nd = theta$alpha.y1nd, beta.y1nd = theta$beta.y1nd, alpha.y1d = theta$alpha.y1d, beta.y1d = theta$beta.y1d,
            alpha.y0nd = theta$alpha.y0nd, beta.y0nd = theta$beta.y0nd, alpha.y0d = theta$alpha.y0d, beta.y0d = theta$beta.y0d, eta.y1 = theta$eta.y1, eta.y2 = theta$eta.y2,
            eta.y3 = theta$eta.y3, delta = theta$delta, pi = pix)

        #### save true values, complete data, and observed data

        save(thetatrue, file = paste(wd, "/sim/thetatrue1_", seed, ".RData", sep = ""))
        rm(theta)

        realdata <- data.frame(Z = Z, Y1 = Y1, Y0 = Y0, x1 = x1st, x2 = x2, x3 = x3, ID = ID, D = D1)
        write.table(realdata, paste(wd, "/sim/realdata1_", seed, ".txt", sep = ""))

        write.table(data, paste(wd, "/sim/simdata1_", seed, ".txt", sep = ""))

    }  # end else 

    ####### POSTERIOR ESTIMATION ########

    # upload the data
    observed.data <- read.table(paste(wd, "/sim/simdata1_", seed, ".txt", sep = ""), header = T)

    # set the seed for the chain
    set.seed(474747)

    # run the chain
    chain1 <- mcmc.tdiscontinue(niter, burn, thin, par.start, theta.prior, observed.data, scenario = 1, seed = seed)

    if (is.na(chain1$nojump) == FALSE) {
        # if the chain didn't jump
        list(ID_save = "nojump")
    } else if (mean(colMeans(chain1$jump[seq({
        burn + thin
    }, niter, by = thin), ])) < 0.2 | mean(colMeans(chain1$jump[seq({
        burn + thin
    }, niter, by = thin), ])) > 0.7) {
        list(ID_save = "toolow|toohigh")  # if the acceptance rate were too low or too high
    } else {
        # else, if everything worked well

        ############################################################################## COMPUTING RESULTS ##############################

        THETA <- chain1$Theta[seq({
            burn + thin
        }, niter, by = thin), ]
        DID <- array(NA, dim = c(dim(chain1$DID[[1]])[1], dim(chain1$DID[[1]])[2], dim(THETA)[1]))
        DID_save <- matrix(NA, nrow = dim(chain1$DID[[1]])[1], ncol = dim(THETA)[1])

        j <- NULL
        for (i in seq({
            burn + thin
        }, niter, by = thin)) {
            j <- sum(j, 1)
            DID[, , j] <- chain1$DID[[i]]
            DID_save[, j] <- chain1$DID[[i]][, 2]
        }

        THETA_list <- list()
        for (i in 1:dim(THETA)[1]) {
            THETA_list[[i]] <- list(eta.0 = THETA[i, 1], eta.1 = THETA[i, 2], eta.2 = THETA[i, 3], eta.3 = THETA[i, 4], alpha.d = THETA[i, 5], beta.d = THETA[i,
                6], eta.d1 = THETA[i, 7], eta.d2 = THETA[i, 8], eta.d3 = THETA[i, 9], alpha.y1nd = THETA[i, 10], beta.y1nd = THETA[i, 11], alpha.y1d = THETA[i, 12],
                beta.y1d = THETA[i, 13], alpha.y0nd = THETA[i, 14], beta.y0nd = THETA[i, 15], alpha.y0d = THETA[i, 16], beta.y0d = THETA[i, 17], eta.y1 = THETA[i,
                  18], eta.y2 = THETA[i, 19], eta.y3 = THETA[i, 20], delta = THETA[i, 21], D = DID[, , i][, 1], ID = DID[, , i][, 2], observed.data = observed.data)
        }

        # apply the estimands' function to the posterior values
        results <- lapply(THETA_list, ce_apply)

        # upload the complete data and the true theta values to compute the 'true' effects
        realdata <- read.table(paste(wd, "/sim/realdata1_", seed, ".txt", sep = ""))
        thetatrue <- myload(paste(wd, "/sim/thetatrue1_", seed, ".RData", sep = ""))
        rdce <- dce(thetatrue, d = seq(1, 5), y = seq(1, 20), realdata, scenario = 1)

        ### PROPORTION OF ND patients
        ID_save <- mean(apply(DID_save, 1, mean, na.rm = TRUE), na.rm = TRUE)

        ### DISTRIBUTIONAL CAUSAL EFFECTS: CI width, coverage, and bias

        # D for each d
        dce.d <- array(NA, dim = c(dim(THETA)[1], 5, 20))
        for (i in 1:dim(THETA)[1]) {
            dce.d[i, , ] <- results[[i]]$DCE_Y_D
        }

        lower_mat_dce.d <- matrix(NA, nrow = 5, ncol = 20)
        upper_mat_dce.d <- matrix(NA, nrow = 5, ncol = 20)
        mean_d_ij <- matrix(NA, 5, 20)

        for (i in 1:5) {
            for (j in 1:20) {
                samples_ij <- dce.d[, i, j]
                hdi_ij <- hdi(samples_ij, credMass = 0.95)
                mean_d_ij[i, j] <- mean(samples_ij, rm.na = TRUE)
                lower_mat_dce.d[i, j] <- hdi_ij[1]
                upper_mat_dce.d[i, j] <- hdi_ij[2]
            }
        }

        width_dce.d <- abs(upper_mat_dce.d - lower_mat_dce.d)  # width

        dced_HPD_include <- matrix(NA, nrow = 5, ncol = 20)

        for (i in 1:5) {
            for (j in 1:20) {
                dced_HPD_include[i, j] <- rdce$DCE_Y_D[i, j] >= lower_mat_dce.d[i, j] & rdce$DCE_Y_D[i, j] <= upper_mat_dce.d[i, j]
            }
        }  # coverage

        dced_bias <- (mean_d_ij - rdce$DCE_Y_D)  # bias

        # average D
        dce.D <- matrix(NA, dim(THETA)[1], 20)
        for (i in 1:dim(THETA)[1]) {
            dce.D[i, ] <- results[[i]]$DCE_D
        }

        hdi_dce.D <- apply(dce.D, 2, hdi)
        mean_dceD <- apply(dce.D, 2, function(x) mean(x, rm.na = TRUE))

        width_dce.D <- abs(hdi_dce.D[2, ] - hdi_dce.D[1, ])  # width

        dceD_HPD_include <- matrix(NA, nrow = 20, ncol = 1)

        for (j in 1:20) {
            dceD_HPD_include[j] <- rdce$DCE_D[j] >= hdi_dce.D[1, j] & rdce$DCE_D[j] <= hdi_dce.D[2, j]
        }  # coverage

        dceD_bias <- (mean_dceD - rdce$DCE_D)  #bias

        # ND
        dce.nd <- matrix(NA, dim(THETA)[1], 20)
        for (i in 1:dim(THETA)[1]) {
            dce.nd[i, ] <- results[[i]]$DCE_ND
        }

        hdi_dce.nd <- hdi(dce.nd)
        mean_dcend <- apply(dce.nd, 2, function(x) mean(x, rm.na = TRUE))

        width_dce.nd <- abs(hdi_dce.nd[2, ] - hdi_dce.nd[1, ])  # width

        dcend_HPD_include <- matrix(NA, nrow = 20, ncol = 1)

        for (j in 1:20) {
            dcend_HPD_include[j] <- rdce$DCE_ND[j] >= hdi_dce.nd[1, j] & rdce$DCE_ND[j] <= hdi_dce.nd[2, j]
        }  # coverage

        dcend_bias <- (mean_dcend - rdce$DCE_ND)  # bias

        # ITT
        itt_dce <- matrix(NA, dim(THETA)[1], 20)
        for (i in 1:dim(THETA)[1]) {
            itt_dce[i, ] <- results[[i]]$DCE
        }

        hdi_itt_dce <- hdi(itt_dce)
        mean_dceitt <- apply(itt_dce, 2, function(x) mean(x, rm.na = TRUE))

        width_itt_dce <- abs(hdi_itt_dce[2, ] - hdi_itt_dce[1, ])  # width

        dceitt_HPD_include <- matrix(NA, nrow = 20, ncol = 1)

        for (j in 1:20) {
            dceitt_HPD_include[j] <- rdce$DCE[j] >= hdi_itt_dce[1, j] & rdce$DCE[j] <= hdi_itt_dce[2, j]
        }  # coverage

        dceitt_bias <- (mean_dceitt - rdce$DCE)  # bias

        ### RESTRICTED MEAN SURVIVAL TIME: CI width, coverage, and bias

        # D for each d
        rmste.d <- array(NA, dim = c(dim(THETA)[1], 5, 20))
        for (i in 1:dim(THETA)[1]) {
            rmste.d[i, , ] <- results[[i]]$RMSTE_d
        }

        lower_mat_rmste.d <- matrix(NA, nrow = 5, ncol = 20)
        upper_mat_rmste.d <- matrix(NA, nrow = 5, ncol = 20)

        mean_rmste.d_ij <- matrix(NA, 5, 20)

        for (i in 1:5) {
            for (j in 1:20) {
                samples_rmste_ij <- rmste.d[, i, j]
                hdi_ij <- hdi(samples_rmste_ij, credMass = 0.95)
                mean_rmste.d_ij[i, j] <- mean(samples_rmste_ij, rm.na = TRUE)
                lower_mat_rmste.d[i, j] <- hdi_ij[1]
                upper_mat_rmste.d[i, j] <- hdi_ij[2]
            }
        }

        width_rmste.d <- abs(upper_mat_rmste.d - lower_mat_rmste.d)  # width

        rmsted_HPD_include <- matrix(NA, nrow = 5, ncol = 20)

        for (i in 1:5) {
            for (j in 1:20) {
                rmsted_HPD_include[i, j] <- rdce$RMSTE_d[i, j] >= lower_mat_rmste.d[i, j] & rdce$RMSTE_d[i, j] <= upper_mat_rmste.d[i, j]
            }
        }  # coverage

        rmsted_bias <- (mean_rmste.d_ij - rdce$RMSTE_d)  # bias

        # average D
        rmste.D <- matrix(NA, dim(THETA)[1], 20)
        for (i in 1:dim(THETA)[1]) {
            rmste.D[i, ] <- results[[i]]$RMSTE_D
        }

        hdi_rmste.D <- apply(rmste.D, 2, hdi)
        mean_rmsteD <- apply(rmste.D, 2, function(x) mean(x, rm.na = TRUE))

        width_rmste.D <- abs(hdi_rmste.D[2, ] - hdi_rmste.D[1, ])  # width

        rmsteD_HPD_include <- matrix(NA, nrow = 20, ncol = 1)

        for (j in 1:20) {
            rmsteD_HPD_include[j] <- rdce$RMSTE_D[j] >= hdi_rmste.D[1, j] & rdce$RMSTE_D[j] <= hdi_rmste.D[2, j]
        }  # coverage

        rmsteD_bias <- (mean_rmsteD - rdce$RMSTE_D)  # bias

        # ND
        rmste.nd <- matrix(NA, dim(THETA)[1], 20)
        for (i in 1:dim(THETA)[1]) {
            rmste.nd[i, ] <- results[[i]]$RMSTE_ND
        }

        hdi_rmste.nd <- hdi(rmste.nd)
        mean_rmstend <- apply(rmste.nd, 2, function(x) mean(x, rm.na = TRUE))

        width_rmste.nd <- abs(hdi_rmste.nd[2, ] - hdi_rmste.nd[1, ])  # width

        rmstend_HPD_include <- matrix(NA, nrow = 20, ncol = 1)

        for (j in 1:20) {
            rmstend_HPD_include[j] <- rdce$RMSTE_ND[j] >= hdi_rmste.nd[1, j] & rdce$RMSTE_ND[j] <= hdi_rmste.nd[2, j]
        }  # coverage

        rmstend_bias <- (mean_rmstend - rdce$RMSTE_ND)  # bias

        # ITT
        itt_rmste <- matrix(NA, dim(THETA)[1], 20)
        for (i in 1:dim(THETA)[1]) {
            itt_rmste[i, ] <- results[[i]]$RMSTE
        }

        hdi_itt_rmste <- hdi(itt_rmste)
        mean_rmsteitt <- apply(itt_rmste, 2, function(x) mean(x, rm.na = TRUE))

        width_itt_rmste <- abs(hdi_itt_rmste[2, ] - hdi_itt_rmste[1, ])  # width

        rmsteitt_HPD_include <- matrix(NA, nrow = 20, ncol = 1)

        for (j in 1:20) {
            rmsteitt_HPD_include[j] <- rdce$RMSTE[j] >= hdi_itt_rmste[1, j] & rdce$RMSTE[j] <= hdi_itt_rmste[2, j]
        }  # coverage

        rmsteitt_bias <- (mean_rmsteitt - rdce$RMSTE)  # bias

        ######## save results

        list(dced_HPD_include = dced_HPD_include, dceD_HPD_include = dceD_HPD_include, dcend_HPD_include = dcend_HPD_include, dceitt_HPD_include = dceitt_HPD_include,
            dced_bias = dced_bias, dceD_bias = dceD_bias, dcend_bias = dcend_bias, dceitt_bias = dceitt_bias, width_hdi_dced = width_dce.d, width_hdi_dceD = width_dce.D,
            width_hdi_dceitt = width_itt_dce, width_hdi_dcend = width_dce.nd, rmsted_HPD_include = rmsted_HPD_include, rmsteD_HPD_include = rmsteD_HPD_include, rmstend_HPD_include = rmstend_HPD_include,
            rmsteitt_HPD_include = rmsteitt_HPD_include, rmsted_bias = rmsted_bias, rmsteD_bias = rmsteD_bias, rmstend_bias = rmstend_bias, rmsteitt_bias = rmsteitt_bias,
            width_hdi_rmsted = width_rmste.d, width_hdi_rmsteD = width_rmste.D, width_hdi_rmsteitt = width_itt_rmste, width_hdi_rmstend = width_rmste.nd, ID_save = ID_save)
    }  # end else 
}  # end function
