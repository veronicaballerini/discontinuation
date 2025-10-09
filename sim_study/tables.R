if(full == TRUE){
  load("final_simdata.RData")
  nsim <- 150

jump1 <- NULL
jump2 <- NULL
jump1_nocov <- NULL
jump2_nocov <- NULL

for (i in 1:nsim) {
    jump1 <- rbind(jump1, is.numeric(sim_SCENARIO1[[i]]$ID_save) == TRUE)
    jump2 <- rbind(jump2, is.numeric(sim_SCENARIO2[[i]]$ID_save) == TRUE)
    jump1_nocov <- rbind(jump1_nocov, is.numeric(sim_SCENARIO1_nocov[[i]]$ID_save) == TRUE)
    jump2_nocov <- rbind(jump2_nocov, is.numeric(sim_SCENARIO2_nocov[[i]]$ID_save) == TRUE)
}

sum(jump1)
sum(jump2)
sum(jump1_nocov)
sum(jump2_nocov)

s1 <- which(jump1 == 1)
s2 <- which(jump2 == 1)
s1_nocov <- which(jump1_nocov == 1)
s2_nocov <- which(jump2_nocov == 1)

########## SCENARIO 1 w covariates

### PROPORTION OF ND patients
ID_save1 <- matrix(NA, nsim, 1)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    ID_save1[j] <- sim_SCENARIO1[[i]]$ID_save
}
pc_ID1 <- mean(ID_save1)
pc_ID1

### DISTRIBUTIONAL CAUSAL EFFECTS: CI width, coverage, and bias

# D for each d
coverage_dced1 <- matrix(NA, 5, 20)

for (d in 1:5) {
    dceHPD_include1_d <- matrix(NA, nsim, 20)
    j <- NULL
    for (i in s1) {
        j <- sum(j, 1)
        dceHPD_include1_d[j, ] <- sim_SCENARIO1[[i]]$dced_HPD_include[d, ]
    }
    coverage_dced1[d, ] <- colMeans(dceHPD_include1_d)  # coverage
}

avg_width_dced1 <- matrix(NA, 5, 20)

for (d in 1:5) {
    dceHPD_width1 <- matrix(NA, nsim, 20)
    j <- NULL
    for (i in s1) {
        j <- sum(j, 1)
        dceHPD_width1[j, ] <- sim_SCENARIO1[[i]]$width_hdi_dced[d, ]
    }
    avg_width_dced1[d, ] <- colMeans(dceHPD_width1)  # width
}

dce.d_mean_bias1 <- matrix(NA, 5, 20)
for (d in 1:5) {
    dce.d_bias1 <- matrix(NA, nrow = nsim, ncol = 20)
    j <- NULL
    for (i in s1) {
        j <- sum(j, 1)
        dce.d_bias1[j, ] <- as.vector(sim_SCENARIO1[[i]]$dced_bias[d, ])
    }
    dce.d_mean_bias1[d, ] <- colMeans(dce.d_bias1)  # bias
}

# average D
dceDHPD_include1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    dceDHPD_include1[j, ] <- sim_SCENARIO1[[i]]$dceD_HPD_include
}
coverage_dceD1 <- colMeans(dceDHPD_include1)  # coverage

dceDHPD_width1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    dceDHPD_width1[j, ] <- sim_SCENARIO1[[i]]$width_hdi_dceD
}

avg_width_dceD1 <- colMeans(dceDHPD_width1)  # width

dce.D_bias1 <- matrix(NA, nrow = nsim, ncol = 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    dce.D_bias1[j, ] <- as.vector(sim_SCENARIO1[[i]]$dceD_bias)
}

dce.D_mean_bias1 <- colMeans(dce.D_bias1)  # bias

# ND
dcendHPD_include1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    dcendHPD_include1[j, ] <- sim_SCENARIO1[[i]]$dcend_HPD_include
}
coverage_dcend1 <- colMeans(dcendHPD_include1)  # coverage

dcendHPD_width1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    dcendHPD_width1[j, ] <- sim_SCENARIO1[[i]]$width_hdi_dcend
}
avg_width_dcend1 <- colMeans(dcendHPD_width1)  # width

dce.nd_bias1 <- matrix(NA, nrow = nsim, ncol = 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    dce.nd_bias1[j, ] <- as.vector(sim_SCENARIO1[[i]]$dcend_bias)
}
dce.nd_mean_bias1 <- colMeans(dce.nd_bias1)  # bias

# ITT
dceittHPD_include1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    dceittHPD_include1[j, ] <- sim_SCENARIO1[[i]]$dceitt_HPD_include
}
coverage_dceitt1 <- colMeans(dceittHPD_include1)  # coverage

dceittHPD_width1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    dceittHPD_width1[j, ] <- sim_SCENARIO1[[i]]$width_hdi_dceitt
}
avg_width_dceitt1 <- colMeans(dceittHPD_width1)  # width

dceitt_bias1 <- matrix(NA, nrow = nsim, ncol = 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    dceitt_bias1[j, ] <- as.vector(sim_SCENARIO1[[i]]$dceitt_bias)
}
dceitt_mean_bias1 <- colMeans(dceitt_bias1)  # bias

### RESTRICTED MEAN SURVIVAL TIME: CI width, coverage, and bias

# D for each d
coverage_rmsted1 <- matrix(NA, 5, 20)

for (d in 1:5) {
    rmsteHPD_include1_d <- matrix(NA, nsim, 20)
    j <- NULL
    for (i in s1) {
        j <- sum(j, 1)
        rmsteHPD_include1_d[j, ] <- sim_SCENARIO1[[i]]$rmsted_HPD_include[d, ]
    }
    coverage_rmsted1[d, ] <- colMeans(rmsteHPD_include1_d)  # coverage
}

avg_width_rmsted1 <- matrix(NA, 5, 20)

for (d in 1:5) {
    rmsteHPD_width1 <- matrix(NA, nsim, 20)
    j <- NULL
    for (i in s1) {
        j <- sum(j, 1)
        rmsteHPD_width1[j, ] <- sim_SCENARIO1[[i]]$width_hdi_rmsted[d, ]
    }
    avg_width_rmsted1[d, ] <- colMeans(rmsteHPD_width1)  # width
}

rmste.d_mean_bias1 <- matrix(NA, 5, 20)
for (d in 1:5) {
    rmste.d_bias1 <- matrix(NA, nrow = nsim, ncol = 20)
    j <- NULL
    for (i in s1) {
        j <- sum(j, 1)
        rmste.d_bias1[j, ] <- as.vector(sim_SCENARIO1[[i]]$rmsted_bias[d, ])
    }
    rmste.d_mean_bias1[d, ] <- colMeans(rmste.d_bias1)  # bias
}

# average D
rmsteDHPD_include1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    rmsteDHPD_include1[j, ] <- sim_SCENARIO1[[i]]$rmsteD_HPD_include
}
coverage_rmsteD1 <- colMeans(rmsteDHPD_include1)  # coverage

rmsteDHPD_width1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    rmsteDHPD_width1[j, ] <- sim_SCENARIO1[[i]]$width_hdi_rmsteD
}

avg_width_rmsteD1 <- colMeans(rmsteDHPD_width1)  # width

rmste.D_bias1 <- matrix(NA, nrow = nsim, ncol = 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    rmste.D_bias1[j, ] <- as.vector(sim_SCENARIO1[[i]]$rmsteD_bias)
}

rmste.D_mean_bias1 <- colMeans(rmste.D_bias1)  # bias

# ND
rmstendHPD_include1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    rmstendHPD_include1[j, ] <- sim_SCENARIO1[[i]]$rmstend_HPD_include
}
coverage_rmstend1 <- colMeans(rmstendHPD_include1)  # coverage

rmstendHPD_width1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    rmstendHPD_width1[j, ] <- sim_SCENARIO1[[i]]$width_hdi_rmstend
}
avg_width_rmstend1 <- colMeans(rmstendHPD_width1)  # width

rmste.nd_bias1 <- matrix(NA, nrow = nsim, ncol = 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    rmste.nd_bias1[j, ] <- as.vector(sim_SCENARIO1[[i]]$rmstend_bias)
}
rmste.nd_mean_bias1 <- colMeans(rmste.nd_bias1)  # bias

# ITT
rmsteittHPD_include1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    rmsteittHPD_include1[j, ] <- sim_SCENARIO1[[i]]$rmsteitt_HPD_include
}
coverage_rmsteitt1 <- colMeans(rmsteittHPD_include1)  # coverage

rmsteittHPD_width1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    rmsteittHPD_width1[j, ] <- sim_SCENARIO1[[i]]$width_hdi_rmsteitt
}
avg_width_rmsteitt1 <- colMeans(rmsteittHPD_width1)  # width

rmsteitt_bias1 <- matrix(NA, nrow = nsim, ncol = 20)
j <- NULL
for (i in s1) {
    j <- sum(j, 1)
    rmsteitt_bias1[j, ] <- as.vector(sim_SCENARIO1[[i]]$rmsteitt_bias)
}
rmsteitt_mean_bias1 <- colMeans(rmsteitt_bias1)  # bias

################################################################################ SCENARIO 1 ################################### w/o covariates
################################################################################ #################################

### PROPORTION OF ND patients
ID_save1_nocov <- matrix(NA, nsim, 1)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    ID_save1_nocov[j] <- sim_SCENARIO1_nocov[[i]]$ID_save_nocov
}
pc_ID1_nocov <- mean(ID_save1_nocov)
pc_ID1_nocov

### DISTRIBUTIONAL CAUSAL EFFECTS: CI width, coverage, and bias

# D for each d
coverage_dced1_nocov <- matrix(NA, 5, 20)
for (d in 1:5) {
    dceHPD_nocov_include1 <- matrix(NA, nsim, 20)
    j <- NULL
    for (i in s1_nocov) {
        j <- sum(j, 1)
        dceHPD_nocov_include1[j, ] <- sim_SCENARIO1_nocov[[i]]$dced_HPD_nocov_include[d, ]
    }
    coverage_dced1_nocov[d, ] <- colMeans(dceHPD_nocov_include1)  # coverage
}

avg_width_dced1_nocov <- matrix(NA, 5, 20)
for (d in 1:5) {
    dceHPD_width1_nocov <- matrix(NA, nsim, 20)
    j <- NULL
    for (i in s1_nocov) {
        j <- sum(j, 1)
        dceHPD_width1_nocov[j, ] <- sim_SCENARIO1_nocov[[i]]$width_hdi_dced_nocov[d, ]
    }
    avg_width_dced1_nocov[d, ] <- colMeans(dceHPD_width1_nocov)  # width
}

dce.d_mean_bias1_nocov <- matrix(NA, 5, 20)
for (d in 1:5) {
    dce.d_bias1_nocov <- matrix(NA, nrow = nsim, ncol = 20)
    j <- NULL
    for (i in s1_nocov) {
        j <- sum(j, 1)
        dce.d_bias1_nocov[j, ] <- as.vector(sim_SCENARIO1_nocov[[i]]$dced_nocov_bias[d, ])
    }
    dce.d_mean_bias1_nocov[d, ] <- colMeans(dce.d_bias1_nocov)  # bias
}

# average D
dceDHPD_nocov_include1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    dceDHPD_nocov_include1[j, ] <- sim_SCENARIO1_nocov[[i]]$dceD_HPD_nocov_include
}
coverage_dceD1_nocov <- colMeans(dceDHPD_nocov_include1)  # coverage

dceDHPD_width1_nocov <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    dceDHPD_width1_nocov[j, ] <- sim_SCENARIO1_nocov[[i]]$width_hdi_dceD_nocov
}
avg_width_dceD1_nocov <- colMeans(dceDHPD_width1_nocov)  # width

dce.D_bias1_nocov <- matrix(NA, nrow = nsim, ncol = 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    dce.D_bias1_nocov[j, ] <- as.vector(sim_SCENARIO1_nocov[[i]]$dceD_nocov_bias)
}
dce.D_mean_bias1_nocov <- colMeans(dce.D_bias1_nocov)  # bias

# ND
dcendHPD_nocov_include1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    dcendHPD_nocov_include1[j, ] <- sim_SCENARIO1_nocov[[i]]$dcend_HPD_nocov_include
}
coverage_dcend1_nocov <- colMeans(dcendHPD_nocov_include1)  # coverage

dcendHPD_width1_nocov <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    dcendHPD_width1_nocov[j, ] <- sim_SCENARIO1_nocov[[i]]$width_hdi_dcend_nocov
}
avg_width_dcend1_nocov <- colMeans(dcendHPD_width1_nocov)  # width

dce.nd_bias1_nocov <- matrix(NA, nrow = nsim, ncol = 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    dce.nd_bias1_nocov[j, ] <- as.vector(sim_SCENARIO1_nocov[[i]]$dcend_nocov_bias)
}
dce.nd_mean_bias1_nocov <- colMeans(dce.nd_bias1_nocov)  # bias

# ITT
dceittHPD_nocov_include1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    dceittHPD_nocov_include1[j, ] <- sim_SCENARIO1_nocov[[i]]$dceitt_HPD_nocov_include
}
coverage_dceitt1_nocov <- colMeans(dceittHPD_nocov_include1)  # coverage

dceittHPD_width1_nocov <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    dceittHPD_width1_nocov[j, ] <- sim_SCENARIO1_nocov[[i]]$width_hdi_dceitt_nocov
}
avg_width_dceitt1_nocov <- colMeans(dceittHPD_width1_nocov)  # width

dceitt_bias1_nocov <- matrix(NA, nrow = nsim, ncol = 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    dceitt_bias1_nocov[j, ] <- as.vector(sim_SCENARIO1_nocov[[i]]$dceitt_nocov_bias)
}
dceitt_mean_bias1_nocov <- colMeans(dceitt_bias1_nocov)  # bias

### RESTRICTED MEAN SURVIVAL TIME: CI width, coverage, and bias

# D for each d
coverage_rmsted1_nocov <- matrix(NA, 5, 20)
for (d in 1:5) {
    rmsteHPD_nocov_include1 <- matrix(NA, nsim, 20)
    j <- NULL
    for (i in s1_nocov) {
        j <- sum(j, 1)
        rmsteHPD_nocov_include1[j, ] <- sim_SCENARIO1_nocov[[i]]$rmsted_HPD_nocov_include[d, ]
    }
    coverage_rmsted1_nocov[d, ] <- colMeans(rmsteHPD_nocov_include1)  # coverage
}

avg_width_rmsted1_nocov <- matrix(NA, 5, 20)
for (d in 1:5) {
    rmsteHPD_width1_nocov <- matrix(NA, nsim, 20)
    j <- NULL
    for (i in s1_nocov) {
        j <- sum(j, 1)
        rmsteHPD_width1_nocov[j, ] <- sim_SCENARIO1_nocov[[i]]$width_hdi_rmsted_nocov[d, ]
    }
    avg_width_rmsted1_nocov[d, ] <- colMeans(rmsteHPD_width1_nocov)  # width
}

rmste.d_mean_bias1_nocov <- matrix(NA, 5, 20)
for (d in 1:5) {
    rmste.d_bias1_nocov <- matrix(NA, nrow = nsim, ncol = 20)
    j <- NULL
    for (i in s1_nocov) {
        j <- sum(j, 1)
        rmste.d_bias1_nocov[j, ] <- as.vector(sim_SCENARIO1_nocov[[i]]$rmsted_nocov_bias[d, ])
    }
    rmste.d_mean_bias1_nocov[d, ] <- colMeans(rmste.d_bias1_nocov)  # bias
}

# average D
rmsteDHPD_nocov_include1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    rmsteDHPD_nocov_include1[j, ] <- sim_SCENARIO1_nocov[[i]]$rmsteD_HPD_nocov_include
}
coverage_rmsteD1_nocov <- colMeans(rmsteDHPD_nocov_include1)  # coverage

rmsteDHPD_width1_nocov <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    rmsteDHPD_width1_nocov[j, ] <- sim_SCENARIO1_nocov[[i]]$width_hdi_rmsteD_nocov
}
avg_width_rmsteD1_nocov <- colMeans(rmsteDHPD_width1_nocov)  # width

rmste.D_bias1_nocov <- matrix(NA, nrow = nsim, ncol = 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    rmste.D_bias1_nocov[j, ] <- as.vector(sim_SCENARIO1_nocov[[i]]$rmsteD_nocov_bias)
}
rmste.D_mean_bias1_nocov <- colMeans(rmste.D_bias1_nocov)  # bias

# ND
rmstendHPD_nocov_include1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    rmstendHPD_nocov_include1[j, ] <- sim_SCENARIO1_nocov[[i]]$rmstend_HPD_nocov_include
}
coverage_rmstend1_nocov <- colMeans(rmstendHPD_nocov_include1)  # coverage

rmstendHPD_width1_nocov <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    rmstendHPD_width1_nocov[j, ] <- sim_SCENARIO1_nocov[[i]]$width_hdi_rmstend_nocov
}
avg_width_rmstend1_nocov <- colMeans(rmstendHPD_width1_nocov)  # width

rmste.nd_bias1_nocov <- matrix(NA, nrow = nsim, ncol = 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    rmste.nd_bias1_nocov[j, ] <- as.vector(sim_SCENARIO1_nocov[[i]]$rmstend_nocov_bias)
}
rmste.nd_mean_bias1_nocov <- colMeans(rmste.nd_bias1_nocov)  # bias

# ITT
rmsteittHPD_nocov_include1 <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    rmsteittHPD_nocov_include1[j, ] <- sim_SCENARIO1_nocov[[i]]$rmsteitt_HPD_nocov_include
}
coverage_rmsteitt1_nocov <- colMeans(rmsteittHPD_nocov_include1)  # coverage

rmsteittHPD_width1_nocov <- matrix(NA, nsim, 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    rmsteittHPD_width1_nocov[j, ] <- sim_SCENARIO1_nocov[[i]]$width_hdi_rmsteitt_nocov
}
avg_width_rmsteitt1_nocov <- colMeans(rmsteittHPD_width1_nocov)  # width

rmsteitt_bias1_nocov <- matrix(NA, nrow = nsim, ncol = 20)
j <- NULL
for (i in s1_nocov) {
    j <- sum(j, 1)
    rmsteitt_bias1_nocov[j, ] <- as.vector(sim_SCENARIO1_nocov[[i]]$rmsteitt_nocov_bias)
}
rmsteitt_mean_bias1_nocov <- colMeans(rmsteitt_bias1_nocov)  # bias

################################################################################ SCENARIO 2 ################################### w covariates
################################################################################ ##################################

n2 <- length(s2)

### PROPORTION OF ND patients
ID_save2 <- matrix(NA, n2, 1)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    ID_save2[j] <- sim_SCENARIO2[[i]]$ID_save
}
pc_ID2 <- mean(ID_save2)

### DISTRIBUTIONAL CAUSAL EFFECTS: CI width, coverage, and bias

# D for each d
coverage_dced2 <- matrix(NA, 5, 20)

for (d in 1:5) {
    dceHPD_include2_d <- matrix(NA, n2, 20)
    j <- NULL
    for (i in s2) {
        j <- sum(j, 1)
        dceHPD_include2_d[j, ] <- sim_SCENARIO2[[i]]$dced_HPD_include[d, ]
    }
    coverage_dced2[d, ] <- colMeans(dceHPD_include2_d)  # coverage
}

avg_width_dced2 <- matrix(NA, 5, 20)

for (d in 1:5) {
    dceHPD_width2 <- matrix(NA, n2, 20)
    j <- NULL
    for (i in s2) {
        j <- sum(j, 1)
        dceHPD_width2[j, ] <- sim_SCENARIO2[[i]]$width_hdi_dced[d, ]
    }
    avg_width_dced2[d, ] <- colMeans(dceHPD_width2)  # width
}

avg_width_dced2  # width

dce.d_mean_bias2 <- matrix(NA, 5, 20)
for (d in 1:5) {
    dce.d_bias2 <- matrix(NA, nrow = n2, ncol = 20)
    j <- NULL
    for (i in s2) {
        j <- sum(j, 1)
        dce.d_bias2[j, ] <- as.vector(sim_SCENARIO2[[i]]$dced_bias[d, ])
    }
    dce.d_mean_bias2[d, ] <- colMeans(dce.d_bias2)  # bias
}

# average D
dceDHPD_include2 <- matrix(NA, n2, 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    dceDHPD_include2[j, ] <- sim_SCENARIO2[[i]]$dceD_HPD_include
}
coverage_dceD2 <- colMeans(dceDHPD_include2)  # coverage

dceDHPD_width2 <- matrix(NA, n2, 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    dceDHPD_width2[j, ] <- sim_SCENARIO2[[i]]$width_hdi_dceD
}
avg_width_dceD2 <- colMeans(dceDHPD_width2)  # width

dce.D_bias2 <- matrix(NA, nrow = n2, ncol = 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    dce.D_bias2[j, ] <- as.vector(sim_SCENARIO2[[i]]$dceD_bias)
}
dce.D_mean_bias2 <- colMeans(dce.D_bias2)  # bias

# ND
dcendHPD_include2 <- matrix(NA, n2, 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    dcendHPD_include2[j, ] <- sim_SCENARIO2[[i]]$dcend_HPD_include
}
coverage_dcend2 <- colMeans(dcendHPD_include2)  # coverage

dcendHPD_width2 <- matrix(NA, n2, 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    dcendHPD_width2[j, ] <- sim_SCENARIO2[[i]]$width_hdi_dcend
}
avg_width_dcend2 <- colMeans(dcendHPD_width2)  # width

dce.nd_bias2 <- matrix(NA, nrow = n2, ncol = 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    dce.nd_bias2[j, ] <- as.vector(sim_SCENARIO2[[i]]$dcend_bias)
}
dce.nd_mean_bias2 <- colMeans(dce.nd_bias2)  # bias

# ITT
dceittHPD_include2 <- matrix(NA, n2, 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    dceittHPD_include2[j, ] <- sim_SCENARIO2[[i]]$dceitt_HPD_include
}
coverage_dceitt2 <- colMeans(dceittHPD_include2)  # coverage

dceittHPD_width2 <- matrix(NA, n2, 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    dceittHPD_width2[j, ] <- sim_SCENARIO2[[i]]$width_hdi_dceitt
}
avg_width_dceitt2 <- colMeans(dceittHPD_width2)  # width

dceitt_bias2 <- matrix(NA, nrow = n2, ncol = 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    dceitt_bias2[j, ] <- as.vector(sim_SCENARIO2[[i]]$dceitt_bias)
}
dceitt_mean_bias2 <- colMeans(dceitt_bias2)  # bias

### RESTRICTED MEAN SURVIVAL TIME: CI width, coverage, and bias

# D for each d
coverage_rmsted2 <- matrix(NA, 5, 20)

for (d in 1:5) {
    rmsteHPD_include2_d <- matrix(NA, n2, 20)
    j <- NULL
    for (i in s2) {
        j <- sum(j, 1)
        rmsteHPD_include2_d[j, ] <- sim_SCENARIO2[[i]]$rmsted_HPD_include[d, ]
    }
    coverage_rmsted2[d, ] <- colMeans(rmsteHPD_include2_d)  # coverage
}

avg_width_rmsted2 <- matrix(NA, 5, 20)

for (d in 1:5) {
    rmsteHPD_width2 <- matrix(NA, n2, 20)
    j <- NULL
    for (i in s2) {
        j <- sum(j, 1)
        rmsteHPD_width2[j, ] <- sim_SCENARIO2[[i]]$width_hdi_rmsted[d, ]
    }
    avg_width_rmsted2[d, ] <- colMeans(rmsteHPD_width2)  # width
}

avg_width_rmsted2  # width

rmste.d_mean_bias2 <- matrix(NA, 5, 20)
for (d in 1:5) {
    rmste.d_bias2 <- matrix(NA, nrow = n2, ncol = 20)
    j <- NULL
    for (i in s2) {
        j <- sum(j, 1)
        rmste.d_bias2[j, ] <- as.vector(sim_SCENARIO2[[i]]$rmsted_bias[d, ])
    }
    rmste.d_mean_bias2[d, ] <- colMeans(rmste.d_bias2)  # bias
}

# average D
rmsteDHPD_include2 <- matrix(NA, n2, 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    rmsteDHPD_include2[j, ] <- sim_SCENARIO2[[i]]$rmsteD_HPD_include
}
coverage_rmsteD2 <- colMeans(rmsteDHPD_include2)  # coverage

rmsteDHPD_width2 <- matrix(NA, n2, 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    rmsteDHPD_width2[j, ] <- sim_SCENARIO2[[i]]$width_hdi_rmsteD
}
avg_width_rmsteD2 <- colMeans(rmsteDHPD_width2)  # width

rmste.D_bias2 <- matrix(NA, nrow = n2, ncol = 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    rmste.D_bias2[j, ] <- as.vector(sim_SCENARIO2[[i]]$rmsteD_bias)
}
rmste.D_mean_bias2 <- colMeans(rmste.D_bias2)  # bias

# ND
rmstendHPD_include2 <- matrix(NA, n2, 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    rmstendHPD_include2[j, ] <- sim_SCENARIO2[[i]]$rmstend_HPD_include
}
coverage_rmstend2 <- colMeans(rmstendHPD_include2)  # coverage

rmstendHPD_width2 <- matrix(NA, n2, 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    rmstendHPD_width2[j, ] <- sim_SCENARIO2[[i]]$width_hdi_rmstend
}
avg_width_rmstend2 <- colMeans(rmstendHPD_width2)  # width

rmste.nd_bias2 <- matrix(NA, nrow = n2, ncol = 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    rmste.nd_bias2[j, ] <- as.vector(sim_SCENARIO2[[i]]$rmstend_bias)
}
rmste.nd_mean_bias2 <- colMeans(rmste.nd_bias2)  # bias

# ITT
rmsteittHPD_include2 <- matrix(NA, n2, 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    rmsteittHPD_include2[j, ] <- sim_SCENARIO2[[i]]$rmsteitt_HPD_include
}
coverage_rmsteitt2 <- colMeans(rmsteittHPD_include2)  # coverage

rmsteittHPD_width2 <- matrix(NA, n2, 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    rmsteittHPD_width2[j, ] <- sim_SCENARIO2[[i]]$width_hdi_rmsteitt
}
avg_width_rmsteitt2 <- colMeans(rmsteittHPD_width2)  # width

rmsteitt_bias2 <- matrix(NA, nrow = n2, ncol = 20)
j <- NULL
for (i in s2) {
    j <- sum(j, 1)
    rmsteitt_bias2[j, ] <- as.vector(sim_SCENARIO2[[i]]$rmsteitt_bias)
}
rmsteitt_mean_bias2 <- colMeans(rmsteitt_bias2)  # bias

################################################################################ SCENARIO 2 ################################### w/o covariates
################################################################################ #################################

n2_nocov <- length(s2_nocov)

### PROPORTION OF ND patients
ID_save2_nocov <- matrix(NA, n2_nocov, 1)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    ID_save2_nocov[j] <- sim_SCENARIO2_nocov[[i]]$ID_save_nocov
}
pc_ID2_nocov <- mean(ID_save2_nocov)
pc_ID2_nocov

### DISTRIBUTIONAL CAUSAL EFFECTS: CI width, coverage, and bias

# D for each d
coverage_dced2_nocov <- matrix(NA, 5, 20)
for (d in 1:5) {
    dceHPD_nocov_include2 <- matrix(NA, n2_nocov, 20)
    j <- NULL
    for (i in s2_nocov) {
        j <- sum(j, 1)
        dceHPD_nocov_include2[j, ] <- sim_SCENARIO2_nocov[[i]]$dced_HPD_nocov_include[d, ]
    }
    coverage_dced2_nocov[d, ] <- colMeans(dceHPD_nocov_include2)  # coverage
}

avg_width_dced2_nocov <- matrix(NA, 5, 20)
for (d in 1:5) {
    dceHPD_width2_nocov <- matrix(NA, n2_nocov, 20)
    j <- NULL
    for (i in s2_nocov) {
        j <- sum(j, 1)
        dceHPD_width2_nocov[j, ] <- sim_SCENARIO2_nocov[[i]]$width_hdi_dced_nocov[d, ]
    }
    avg_width_dced2_nocov[d, ] <- colMeans(dceHPD_width2_nocov)  # width
}

dce.d_mean_bias2_nocov <- matrix(NA, 5, 20)
for (d in 1:5) {
    dce.d_bias2_nocov <- matrix(NA, nrow = n2_nocov, ncol = 20)
    j <- NULL
    for (i in s2_nocov) {
        j <- sum(j, 1)
        dce.d_bias2_nocov[j, ] <- as.vector(sim_SCENARIO2_nocov[[i]]$dced_nocov_bias[d, ])
    }
    dce.d_mean_bias2_nocov[d, ] <- colMeans(dce.d_bias2_nocov)  # bias
}

# average D
dceDHPD_nocov_include2 <- matrix(NA, n2_nocov, 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    dceDHPD_nocov_include2[j, ] <- sim_SCENARIO2_nocov[[i]]$dceD_HPD_nocov_include
}
coverage_dceD2_nocov <- colMeans(dceDHPD_nocov_include2)  # coverage

dceDHPD_width2_nocov <- matrix(NA, n2_nocov, 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    dceDHPD_width2_nocov[j, ] <- sim_SCENARIO2_nocov[[i]]$width_hdi_dceD_nocov
}
avg_width_dceD2_nocov <- colMeans(dceDHPD_width2_nocov)  # width

dce.D_bias2_nocov <- matrix(NA, nrow = n2_nocov, ncol = 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    dce.D_bias2_nocov[j, ] <- as.vector(sim_SCENARIO2_nocov[[i]]$dceD_nocov_bias)
}
dce.D_mean_bias2_nocov <- colMeans(dce.D_bias2_nocov)  # bias 

# ND
dcendHPD_nocov_include2 <- matrix(NA, n2_nocov, 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    dcendHPD_nocov_include2[j, ] <- sim_SCENARIO2_nocov[[i]]$dcend_HPD_nocov_include
}
coverage_dcend2_nocov <- colMeans(dcendHPD_nocov_include2)  # coverage

dcendHPD_width2_nocov <- matrix(NA, n2_nocov, 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    dcendHPD_width2_nocov[j, ] <- sim_SCENARIO2_nocov[[i]]$width_hdi_dcend_nocov
}
avg_width_dcend2_nocov <- colMeans(dcendHPD_width2_nocov)  # width

dce.nd_bias2_nocov <- matrix(NA, nrow = n2_nocov, ncol = 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    dce.nd_bias2_nocov[j, ] <- as.vector(sim_SCENARIO2_nocov[[i]]$dcend_nocov_bias)
}
dce.nd_mean_bias2_nocov <- colMeans(dce.nd_bias2_nocov)  # bias

# ITT
dceittHPD_nocov_include2 <- matrix(NA, n2_nocov, 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    dceittHPD_nocov_include2[j, ] <- sim_SCENARIO2_nocov[[i]]$dceitt_HPD_nocov_include
}
coverage_dceitt2_nocov <- colMeans(dceittHPD_nocov_include2)  # coverage

dceittHPD_width2_nocov <- matrix(NA, n2_nocov, 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    dceittHPD_width2_nocov[j, ] <- sim_SCENARIO2_nocov[[i]]$width_hdi_dceitt_nocov
}
avg_width_dceitt2_nocov <- colMeans(dceittHPD_width2_nocov)  # width 

dceitt_bias2_nocov <- matrix(NA, nrow = n2_nocov, ncol = 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    dceitt_bias2_nocov[j, ] <- as.vector(sim_SCENARIO2_nocov[[i]]$dceitt_nocov_bias)
}
dceitt_mean_bias2_nocov <- colMeans(dceitt_bias2_nocov)  # bias

### RESTRICTED MEAN SURVIVAL TIME: CI width, coverage, and bias

# D for each d
coverage_rmsted2_nocov <- matrix(NA, 5, 20)
for (d in 1:5) {
    rmsteHPD_nocov_include2 <- matrix(NA, n2_nocov, 20)
    j <- NULL
    for (i in s2_nocov) {
        j <- sum(j, 1)
        rmsteHPD_nocov_include2[j, ] <- sim_SCENARIO2_nocov[[i]]$rmsted_HPD_nocov_include[d, ]
    }
    coverage_rmsted2_nocov[d, ] <- colMeans(rmsteHPD_nocov_include2)  # coverage
}

avg_width_rmsted2_nocov <- matrix(NA, 5, 20)
for (d in 1:5) {
    rmsteHPD_width2_nocov <- matrix(NA, n2_nocov, 20)
    j <- NULL
    for (i in s2_nocov) {
        j <- sum(j, 1)
        rmsteHPD_width2_nocov[j, ] <- sim_SCENARIO2_nocov[[i]]$width_hdi_rmsted_nocov[d, ]
    }
    avg_width_rmsted2_nocov[d, ] <- colMeans(rmsteHPD_width2_nocov)  # width
}

rmste.d_mean_bias2_nocov <- matrix(NA, 5, 20)
for (d in 1:5) {
    rmste.d_bias2_nocov <- matrix(NA, nrow = n2_nocov, ncol = 20)
    j <- NULL
    for (i in s2_nocov) {
        j <- sum(j, 1)
        rmste.d_bias2_nocov[j, ] <- as.vector(sim_SCENARIO2_nocov[[i]]$rmsted_nocov_bias[d, ])
    }
    rmste.d_mean_bias2_nocov[d, ] <- colMeans(rmste.d_bias2_nocov)  # bias
}

# average D
rmsteDHPD_nocov_include2 <- matrix(NA, n2_nocov, 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    rmsteDHPD_nocov_include2[j, ] <- sim_SCENARIO2_nocov[[i]]$rmsteD_HPD_nocov_include
}
coverage_rmsteD2_nocov <- colMeans(rmsteDHPD_nocov_include2)  # coverage

rmsteDHPD_width2_nocov <- matrix(NA, n2_nocov, 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    rmsteDHPD_width2_nocov[j, ] <- sim_SCENARIO2_nocov[[i]]$width_hdi_rmsteD_nocov
}
avg_width_rmsteD2_nocov <- colMeans(rmsteDHPD_width2_nocov)  # width

rmste.D_bias2_nocov <- matrix(NA, nrow = n2_nocov, ncol = 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    rmste.D_bias2_nocov[j, ] <- as.vector(sim_SCENARIO2_nocov[[i]]$rmsteD_nocov_bias)
}
rmste.D_mean_bias2_nocov <- colMeans(rmste.D_bias2_nocov)  # bias

# ND
rmstendHPD_nocov_include2 <- matrix(NA, n2_nocov, 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    rmstendHPD_nocov_include2[j, ] <- sim_SCENARIO2_nocov[[i]]$rmstend_HPD_nocov_include
}
coverage_rmstend2_nocov <- colMeans(rmstendHPD_nocov_include2)  # coverage

rmstendHPD_width2_nocov <- matrix(NA, n2_nocov, 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    rmstendHPD_width2_nocov[j, ] <- sim_SCENARIO2_nocov[[i]]$width_hdi_rmstend_nocov
}
avg_width_rmstend2_nocov <- colMeans(rmstendHPD_width2_nocov)  # width

rmste.nd_bias2_nocov <- matrix(NA, nrow = n2_nocov, ncol = 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    rmste.nd_bias2_nocov[j, ] <- as.vector(sim_SCENARIO2_nocov[[i]]$rmstend_nocov_bias)
}
rmste.nd_mean_bias2_nocov <- colMeans(rmste.nd_bias2_nocov)  # bias

# ITT
rmsteittHPD_nocov_include2 <- matrix(NA, n2_nocov, 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    rmsteittHPD_nocov_include2[j, ] <- sim_SCENARIO2_nocov[[i]]$rmsteitt_HPD_nocov_include
}
coverage_rmsteitt2_nocov <- colMeans(rmsteittHPD_nocov_include2)  # coverage

rmsteittHPD_width2_nocov <- matrix(NA, n2_nocov, 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    rmsteittHPD_width2_nocov[j, ] <- sim_SCENARIO2_nocov[[i]]$width_hdi_rmsteitt_nocov
}
avg_width_rmsteitt2_nocov <- colMeans(rmsteittHPD_width2_nocov)  # width

rmsteitt_bias2_nocov <- matrix(NA, nrow = n2_nocov, ncol = 20)
j <- NULL
for (i in s2_nocov) {
    j <- sum(j, 1)
    rmsteitt_bias2_nocov[j, ] <- as.vector(sim_SCENARIO2_nocov[[i]]$rmsteitt_nocov_bias)
}
rmsteitt_mean_bias2_nocov <- colMeans(rmsteitt_bias2_nocov)  # bias

########################################################################## TABLES DCE - coverage
DCE_coverage_scenario1 <- data.frame(DCE_ITT = as.matrix(round(coverage_dceitt1, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(coverage_dcend1, 2), nrow = 20,
    ncol = 1), DCE_D = as.matrix(round(coverage_dceD1, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(coverage_dced1[1, ], 2), nrow = 20, ncol = 1), DCE_d2 = as.matrix(round(coverage_dced1[2,
    ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(coverage_dced1[3, ], 2), nrow = 20, ncol = 1), DCE_d4 = as.matrix(round(coverage_dced1[4, ], 2), nrow = 20,
    ncol = 1), DCE_d5 = as.matrix(round(coverage_dced1[5, ], 2), nrow = 20, ncol = 1))

DCE_coverage_scenario1_nocov <- data.frame(DCE_ITT = as.matrix(round(coverage_dceitt1_nocov, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(coverage_dcend1_nocov,
    2), nrow = 20, ncol = 1), DCE_D = as.matrix(round(coverage_dceD1_nocov, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(coverage_dced1_nocov[1, ], 2), nrow = 20,
    ncol = 1), DCE_d2 = as.matrix(round(coverage_dced1_nocov[2, ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(coverage_dced1_nocov[3, ], 2), nrow = 20, ncol = 1),
    DCE_d4 = as.matrix(round(coverage_dced1_nocov[4, ], 2), nrow = 20, ncol = 1), DCE_d5 = as.matrix(round(coverage_dced1_nocov[5, ], 2), nrow = 20, ncol = 1))

DCE_coverage_scenario2 <- data.frame(DCE_ITT = as.matrix(round(coverage_dceitt2, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(coverage_dcend2, 2), nrow = 20,
    ncol = 1), DCE_D = as.matrix(round(coverage_dceD2, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(coverage_dced2[1, ], 2), nrow = 20, ncol = 1), DCE_d2 = as.matrix(round(coverage_dced2[2,
    ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(coverage_dced2[3, ], 2), nrow = 20, ncol = 1), DCE_d4 = as.matrix(round(coverage_dced2[4, ], 2), nrow = 20,
    ncol = 1), DCE_d5 = as.matrix(round(coverage_dced2[5, ], 2), nrow = 20, ncol = 1))

DCE_coverage_scenario2_nocov <- data.frame(DCE_ITT = as.matrix(round(coverage_dceitt2_nocov, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(coverage_dcend2_nocov,
    2), nrow = 20, ncol = 1), DCE_D = as.matrix(round(coverage_dceD2_nocov, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(coverage_dced2_nocov[1, ], 2), nrow = 20,
    ncol = 1), DCE_d2 = as.matrix(round(coverage_dced2_nocov[2, ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(coverage_dced2_nocov[3, ], 2), nrow = 20, ncol = 1),
    DCE_d4 = as.matrix(round(coverage_dced2_nocov[4, ], 2), nrow = 20, ncol = 1), DCE_d5 = as.matrix(round(coverage_dced2_nocov[5, ], 2), nrow = 20, ncol = 1))

# DCE - width
DCE_avg_width_scenario1 <- data.frame(DCE_ITT = as.matrix(round(avg_width_dceitt1, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(avg_width_dcend1, 2), nrow = 20,
    ncol = 1), DCE_D = as.matrix(round(avg_width_dceD1, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(avg_width_dced1[1, ], 2), nrow = 20, ncol = 1), DCE_d2 = as.matrix(round(avg_width_dced1[2,
    ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(avg_width_dced1[3, ], 2), nrow = 20, ncol = 1), DCE_d4 = as.matrix(round(avg_width_dced1[4, ], 2), nrow = 20,
    ncol = 1), DCE_d5 = as.matrix(round(avg_width_dced1[5, ], 2), nrow = 20, ncol = 1))

DCE_avg_width_scenario1_nocov <- data.frame(DCE_ITT = as.matrix(round(avg_width_dceitt1_nocov, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(avg_width_dcend1_nocov,
    2), nrow = 20, ncol = 1), DCE_D = as.matrix(round(avg_width_dceD1_nocov, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(avg_width_dced1_nocov[1, ], 2), nrow = 20,
    ncol = 1), DCE_d2 = as.matrix(round(avg_width_dced1_nocov[2, ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(avg_width_dced1_nocov[3, ], 2), nrow = 20,
    ncol = 1), DCE_d4 = as.matrix(round(avg_width_dced1_nocov[4, ], 2), nrow = 20, ncol = 1), DCE_d5 = as.matrix(round(avg_width_dced1_nocov[5, ], 2), nrow = 20,
    ncol = 1))

DCE_avg_width_scenario2 <- data.frame(DCE_ITT = as.matrix(round(avg_width_dceitt2, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(avg_width_dcend2, 2), nrow = 20,
    ncol = 1), DCE_D = as.matrix(round(avg_width_dceD2, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(avg_width_dced2[1, ], 2), nrow = 20, ncol = 1), DCE_d2 = as.matrix(round(avg_width_dced2[2,
    ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(avg_width_dced2[3, ], 2), nrow = 20, ncol = 1), DCE_d4 = as.matrix(round(avg_width_dced2[4, ], 2), nrow = 20,
    ncol = 1), DCE_d5 = as.matrix(round(avg_width_dced2[5, ], 2), nrow = 20, ncol = 1))

DCE_avg_width_scenario2_nocov <- data.frame(DCE_ITT = as.matrix(round(avg_width_dceitt2_nocov, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(avg_width_dcend2_nocov,
    2), nrow = 20, ncol = 1), DCE_D = as.matrix(round(avg_width_dceD2_nocov, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(avg_width_dced2_nocov[1, ], 2), nrow = 20,
    ncol = 1), DCE_d2 = as.matrix(round(avg_width_dced2_nocov[2, ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(avg_width_dced2_nocov[3, ], 2), nrow = 20,
    ncol = 1), DCE_d4 = as.matrix(round(avg_width_dced2_nocov[4, ], 2), nrow = 20, ncol = 1), DCE_d5 = as.matrix(round(avg_width_dced2_nocov[5, ], 2), nrow = 20,
    ncol = 1))

# DCE - bias
DCE_bias_scenario1 <- data.frame(DCE_ITT = as.matrix(round(dceitt_mean_bias1, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(dce.nd_mean_bias1, 2), nrow = 20,
    ncol = 1), DCE_D = as.matrix(round(dce.D_mean_bias1, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(dce.d_mean_bias1[1, ], 2), nrow = 20, ncol = 1), DCE_d2 = as.matrix(round(dce.d_mean_bias1[2,
    ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(dce.d_mean_bias1[3, ], 2), nrow = 20, ncol = 1), DCE_d4 = as.matrix(round(dce.d_mean_bias1[4, ], 2), nrow = 20,
    ncol = 1), DCE_d5 = as.matrix(round(dce.d_mean_bias1[5, ], 2), nrow = 20, ncol = 1))

DCE_bias_scenario1_nocov <- data.frame(DCE_ITT = as.matrix(round(dceitt_mean_bias1_nocov, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(dce.nd_mean_bias1_nocov,
    2), nrow = 20, ncol = 1), DCE_D = as.matrix(round(dce.D_mean_bias1_nocov, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(dce.d_mean_bias1_nocov[1, ], 2),
    nrow = 20, ncol = 1), DCE_d2 = as.matrix(round(dce.d_mean_bias1_nocov[2, ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(dce.d_mean_bias1_nocov[3, ], 2),
    nrow = 20, ncol = 1), DCE_d4 = as.matrix(round(dce.d_mean_bias1_nocov[4, ], 2), nrow = 20, ncol = 1), DCE_d5 = as.matrix(round(dce.d_mean_bias1_nocov[5, ], 2),
    nrow = 20, ncol = 1))

DCE_bias_scenario2 <- data.frame(DCE_ITT = as.matrix(round(dceitt_mean_bias2, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(dce.nd_mean_bias2, 2), nrow = 20,
    ncol = 1), DCE_D = as.matrix(round(dce.D_mean_bias2, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(dce.d_mean_bias2[1, ], 2), nrow = 20, ncol = 1), DCE_d2 = as.matrix(round(dce.d_mean_bias2[2,
    ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(dce.d_mean_bias2[3, ], 2), nrow = 20, ncol = 1), DCE_d4 = as.matrix(round(dce.d_mean_bias2[4, ], 2), nrow = 20,
    ncol = 1), DCE_d5 = as.matrix(round(dce.d_mean_bias2[5, ], 2), nrow = 20, ncol = 1))

DCE_bias_scenario2_nocov <- data.frame(DCE_ITT = as.matrix(round(dceitt_mean_bias2_nocov, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(dce.nd_mean_bias2_nocov,
    2), nrow = 20, ncol = 1), DCE_D = as.matrix(round(dce.D_mean_bias2_nocov, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(dce.d_mean_bias2_nocov[1, ], 2),
    nrow = 20, ncol = 1), DCE_d2 = as.matrix(round(dce.d_mean_bias2_nocov[2, ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(dce.d_mean_bias2_nocov[3, ], 2),
    nrow = 20, ncol = 1), DCE_d4 = as.matrix(round(dce.d_mean_bias2_nocov[4, ], 2), nrow = 20, ncol = 1), DCE_d5 = as.matrix(round(dce.d_mean_bias2_nocov[5, ], 2),
    nrow = 20, ncol = 1))

# RMSTE - coverage
RMSTE_coverage_scenario1 <- data.frame(RMSTE_ITT = as.matrix(round(coverage_rmsteitt1, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(coverage_rmstend1, 2),
    nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(coverage_rmsteD1, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(coverage_rmsted1[1, ], 2), nrow = 20,
    ncol = 1), RMSTE_d2 = as.matrix(round(coverage_rmsted1[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(coverage_rmsted1[3, ], 2), nrow = 20, ncol = 1),
    RMSTE_d4 = as.matrix(round(coverage_rmsted1[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(coverage_rmsted1[5, ], 2), nrow = 20, ncol = 1))

RMSTE_coverage_scenario1_nocov <- data.frame(RMSTE_ITT = as.matrix(round(coverage_rmsteitt1_nocov, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(coverage_rmstend1_nocov,
    2), nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(coverage_rmsteD1_nocov, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(coverage_rmsted1_nocov[1, ],
    2), nrow = 20, ncol = 1), RMSTE_d2 = as.matrix(round(coverage_rmsted1_nocov[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(coverage_rmsted1_nocov[3,
    ], 2), nrow = 20, ncol = 1), RMSTE_d4 = as.matrix(round(coverage_rmsted1_nocov[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(coverage_rmsted1_nocov[5,
    ], 2), nrow = 20, ncol = 1))

RMSTE_coverage_scenario2 <- data.frame(RMSTE_ITT = as.matrix(round(coverage_rmsteitt2, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(coverage_rmstend2, 2),
    nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(coverage_rmsteD2, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(coverage_rmsted2[1, ], 2), nrow = 20,
    ncol = 1), RMSTE_d2 = as.matrix(round(coverage_rmsted2[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(coverage_rmsted2[3, ], 2), nrow = 20, ncol = 1),
    RMSTE_d4 = as.matrix(round(coverage_rmsted2[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(coverage_rmsted2[5, ], 2), nrow = 20, ncol = 1))

RMSTE_coverage_scenario2_nocov <- data.frame(RMSTE_ITT = as.matrix(round(coverage_rmsteitt2_nocov, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(coverage_rmstend2_nocov,
    2), nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(coverage_rmsteD2_nocov, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(coverage_rmsted2_nocov[1, ],
    2), nrow = 20, ncol = 1), RMSTE_d2 = as.matrix(round(coverage_rmsted2_nocov[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(coverage_rmsted2_nocov[3,
    ], 2), nrow = 20, ncol = 1), RMSTE_d4 = as.matrix(round(coverage_rmsted2_nocov[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(coverage_rmsted2_nocov[5,
    ], 2), nrow = 20, ncol = 1))

# RMSTE - width
RMSTE_avg_width_scenario1 <- data.frame(RMSTE_ITT = as.matrix(round(avg_width_rmsteitt1, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(avg_width_rmstend1,
    2), nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(avg_width_rmsteD1, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(avg_width_rmsted1[1, ], 2), nrow = 20,
    ncol = 1), RMSTE_d2 = as.matrix(round(avg_width_rmsted1[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(avg_width_rmsted1[3, ], 2), nrow = 20, ncol = 1),
    RMSTE_d4 = as.matrix(round(avg_width_rmsted1[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(avg_width_rmsted1[5, ], 2), nrow = 20, ncol = 1))

RMSTE_avg_width_scenario1_nocov <- data.frame(RMSTE_ITT = as.matrix(round(avg_width_rmsteitt1_nocov, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(avg_width_rmstend1_nocov,
    2), nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(avg_width_rmsteD1_nocov, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(avg_width_rmsted1_nocov[1,
    ], 2), nrow = 20, ncol = 1), RMSTE_d2 = as.matrix(round(avg_width_rmsted1_nocov[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(avg_width_rmsted1_nocov[3,
    ], 2), nrow = 20, ncol = 1), RMSTE_d4 = as.matrix(round(avg_width_rmsted1_nocov[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(avg_width_rmsted1_nocov[5,
    ], 2), nrow = 20, ncol = 1))

RMSTE_avg_width_scenario2 <- data.frame(RMSTE_ITT = as.matrix(round(avg_width_rmsteitt2, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(avg_width_rmstend2,
    2), nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(avg_width_rmsteD2, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(avg_width_rmsted2[1, ], 2), nrow = 20,
    ncol = 1), RMSTE_d2 = as.matrix(round(avg_width_rmsted2[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(avg_width_rmsted2[3, ], 2), nrow = 20, ncol = 1),
    RMSTE_d4 = as.matrix(round(avg_width_rmsted2[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(avg_width_rmsted2[5, ], 2), nrow = 20, ncol = 1))

RMSTE_avg_width_scenario2_nocov <- data.frame(RMSTE_ITT = as.matrix(round(avg_width_rmsteitt2_nocov, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(avg_width_rmstend2_nocov,
    2), nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(avg_width_rmsteD2_nocov, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(avg_width_rmsted2_nocov[1,
    ], 2), nrow = 20, ncol = 1), RMSTE_d2 = as.matrix(round(avg_width_rmsted2_nocov[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(avg_width_rmsted2_nocov[3,
    ], 2), nrow = 20, ncol = 1), RMSTE_d4 = as.matrix(round(avg_width_rmsted2_nocov[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(avg_width_rmsted2_nocov[5,
    ], 2), nrow = 20, ncol = 1))

# RMSTE - bias
RMSTE_bias_scenario1 <- data.frame(RMSTE_ITT = as.matrix(round(rmsteitt_mean_bias1, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(rmste.nd_mean_bias1, 2),
    nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(rmste.D_mean_bias1, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(rmste.d_mean_bias1[1, ], 2), nrow = 20,
    ncol = 1), RMSTE_d2 = as.matrix(round(rmste.d_mean_bias1[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(rmste.d_mean_bias1[3, ], 2), nrow = 20, ncol = 1),
    RMSTE_d4 = as.matrix(round(rmste.d_mean_bias1[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(rmste.d_mean_bias1[5, ], 2), nrow = 20, ncol = 1))

RMSTE_bias_scenario1_nocov <- data.frame(RMSTE_ITT = as.matrix(round(rmsteitt_mean_bias1_nocov, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(rmste.nd_mean_bias1_nocov,
    2), nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(rmste.D_mean_bias1_nocov, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(rmste.d_mean_bias1_nocov[1,
    ], 2), nrow = 20, ncol = 1), RMSTE_d2 = as.matrix(round(rmste.d_mean_bias1_nocov[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(rmste.d_mean_bias1_nocov[3,
    ], 2), nrow = 20, ncol = 1), RMSTE_d4 = as.matrix(round(rmste.d_mean_bias1_nocov[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(rmste.d_mean_bias1_nocov[5,
    ], 2), nrow = 20, ncol = 1))

RMSTE_bias_scenario2 <- data.frame(RMSTE_ITT = as.matrix(round(rmsteitt_mean_bias2, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(rmste.nd_mean_bias2, 2),
    nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(rmste.D_mean_bias2, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(rmste.d_mean_bias2[1, ], 2), nrow = 20,
    ncol = 1), RMSTE_d2 = as.matrix(round(rmste.d_mean_bias2[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(rmste.d_mean_bias2[3, ], 2), nrow = 20, ncol = 1),
    RMSTE_d4 = as.matrix(round(rmste.d_mean_bias2[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(rmste.d_mean_bias2[5, ], 2), nrow = 20, ncol = 1))

RMSTE_bias_scenario2_nocov <- data.frame(RMSTE_ITT = as.matrix(round(rmsteitt_mean_bias2_nocov, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(rmste.nd_mean_bias2_nocov,
    2), nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(rmste.D_mean_bias2_nocov, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(rmste.d_mean_bias2_nocov[1,
    ], 2), nrow = 20, ncol = 1), RMSTE_d2 = as.matrix(round(rmste.d_mean_bias2_nocov[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(rmste.d_mean_bias2_nocov[3,
    ], 2), nrow = 20, ncol = 1), RMSTE_d4 = as.matrix(round(rmste.d_mean_bias2_nocov[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(rmste.d_mean_bias2_nocov[5,
    ], 2), nrow = 20, ncol = 1))

Table4 <- DCE_coverage_scenario1
Table5 <- DCE_coverage_scenario2
Table6 <- DCE_coverage_scenario1_nocov
Table7 <- DCE_coverage_scenario2_nocov

saveRDS(Table4, file = "tables/FULL/Table4.rds")
saveRDS(Table5, file = "tables/FULL/Table5.rds")
saveRDS(Table6, file = "tables/FULL/Table6.rds")
saveRDS(Table7, file = "tables/FULL/Table7.rds")

Table8 <- DCE_avg_width_scenario1
Table9 <- DCE_avg_width_scenario2
Table10 <- DCE_avg_width_scenario1_nocov
Table11 <- DCE_avg_width_scenario2_nocov

saveRDS(Table8, file = "tables/FULL/Table8.rds")
saveRDS(Table9, file = "tables/FULL/Table9.rds")
saveRDS(Table10, file = "tables/FULL/Table10.rds")
saveRDS(Table11, file = "tables/FULL/Table11.rds")

Table12 <- DCE_bias_scenario1
Table13 <- DCE_bias_scenario2
Table14 <- DCE_bias_scenario1_nocov
Table15 <- DCE_bias_scenario2_nocov

saveRDS(Table12, file = "tables/FULL/Table12.rds")
saveRDS(Table13, file = "tables/FULL/Table13.rds")
saveRDS(Table14, file = "tables/FULL/Table14.rds")
saveRDS(Table15, file = "tables/FULL/Table15.rds")

}else{
  load("final_simdata_intermediate.RData")
  nsim <- 2

  jump1_nocov <- NULL
  jump2_nocov <- NULL
  
  for (i in 1:nsim) {
    jump1_nocov <- rbind(jump1_nocov, is.numeric(sim_SCENARIO1_nocov[[i]]$ID_save) == TRUE)
    jump2_nocov <- rbind(jump2_nocov, is.numeric(sim_SCENARIO2_nocov[[i]]$ID_save) == TRUE)
  }
  
  s1_nocov <- which(jump1_nocov == 1)
  s2_nocov <- which(jump2_nocov == 1)
  
  # SCENARIO 1 w/o covariates
  ### PROPORTION OF ND patients
  ID_save1_nocov <- matrix(NA, nsim, 1)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    ID_save1_nocov[j] <- sim_SCENARIO1_nocov[[i]]$ID_save_nocov
  }
  pc_ID1_nocov <- mean(ID_save1_nocov)
  pc_ID1_nocov
  
  ### DISTRIBUTIONAL CAUSAL EFFECTS: CI width, coverage, and bias
  
  # D for each d
  coverage_dced1_nocov <- matrix(NA, 5, 20)
  for (d in 1:5) {
    dceHPD_nocov_include1 <- matrix(NA, nsim, 20)
    j <- NULL
    for (i in s1_nocov) {
      j <- sum(j, 1)
      dceHPD_nocov_include1[j, ] <- sim_SCENARIO1_nocov[[i]]$dced_HPD_nocov_include[d, ]
    }
    coverage_dced1_nocov[d, ] <- colMeans(dceHPD_nocov_include1)  # coverage
  }
  
  avg_width_dced1_nocov <- matrix(NA, 5, 20)
  for (d in 1:5) {
    dceHPD_width1_nocov <- matrix(NA, nsim, 20)
    j <- NULL
    for (i in s1_nocov) {
      j <- sum(j, 1)
      dceHPD_width1_nocov[j, ] <- sim_SCENARIO1_nocov[[i]]$width_hdi_dced_nocov[d, ]
    }
    avg_width_dced1_nocov[d, ] <- colMeans(dceHPD_width1_nocov)  # width
  }
  
  dce.d_mean_bias1_nocov <- matrix(NA, 5, 20)
  for (d in 1:5) {
    dce.d_bias1_nocov <- matrix(NA, nrow = nsim, ncol = 20)
    j <- NULL
    for (i in s1_nocov) {
      j <- sum(j, 1)
      dce.d_bias1_nocov[j, ] <- as.vector(sim_SCENARIO1_nocov[[i]]$dced_nocov_bias[d, ])
    }
    dce.d_mean_bias1_nocov[d, ] <- colMeans(dce.d_bias1_nocov)  # bias
  }
  
  # average D
  dceDHPD_nocov_include1 <- matrix(NA, nsim, 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    dceDHPD_nocov_include1[j, ] <- sim_SCENARIO1_nocov[[i]]$dceD_HPD_nocov_include
  }
  coverage_dceD1_nocov <- colMeans(dceDHPD_nocov_include1)  # coverage
  
  dceDHPD_width1_nocov <- matrix(NA, nsim, 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    dceDHPD_width1_nocov[j, ] <- sim_SCENARIO1_nocov[[i]]$width_hdi_dceD_nocov
  }
  avg_width_dceD1_nocov <- colMeans(dceDHPD_width1_nocov)  # width
  
  dce.D_bias1_nocov <- matrix(NA, nrow = nsim, ncol = 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    dce.D_bias1_nocov[j, ] <- as.vector(sim_SCENARIO1_nocov[[i]]$dceD_nocov_bias)
  }
  dce.D_mean_bias1_nocov <- colMeans(dce.D_bias1_nocov)  # bias
  
  # ND
  dcendHPD_nocov_include1 <- matrix(NA, nsim, 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    dcendHPD_nocov_include1[j, ] <- sim_SCENARIO1_nocov[[i]]$dcend_HPD_nocov_include
  }
  coverage_dcend1_nocov <- colMeans(dcendHPD_nocov_include1)  # coverage
  
  dcendHPD_width1_nocov <- matrix(NA, nsim, 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    dcendHPD_width1_nocov[j, ] <- sim_SCENARIO1_nocov[[i]]$width_hdi_dcend_nocov
  }
  avg_width_dcend1_nocov <- colMeans(dcendHPD_width1_nocov)  # width
  
  dce.nd_bias1_nocov <- matrix(NA, nrow = nsim, ncol = 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    dce.nd_bias1_nocov[j, ] <- as.vector(sim_SCENARIO1_nocov[[i]]$dcend_nocov_bias)
  }
  dce.nd_mean_bias1_nocov <- colMeans(dce.nd_bias1_nocov)  # bias
  
  # ITT
  dceittHPD_nocov_include1 <- matrix(NA, nsim, 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    dceittHPD_nocov_include1[j, ] <- sim_SCENARIO1_nocov[[i]]$dceitt_HPD_nocov_include
  }
  coverage_dceitt1_nocov <- colMeans(dceittHPD_nocov_include1)  # coverage
  
  dceittHPD_width1_nocov <- matrix(NA, nsim, 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    dceittHPD_width1_nocov[j, ] <- sim_SCENARIO1_nocov[[i]]$width_hdi_dceitt_nocov
  }
  avg_width_dceitt1_nocov <- colMeans(dceittHPD_width1_nocov)  # width
  
  dceitt_bias1_nocov <- matrix(NA, nrow = nsim, ncol = 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    dceitt_bias1_nocov[j, ] <- as.vector(sim_SCENARIO1_nocov[[i]]$dceitt_nocov_bias)
  }
  dceitt_mean_bias1_nocov <- colMeans(dceitt_bias1_nocov)  # bias
  
  ### RESTRICTED MEAN SURVIVAL TIME: CI width, coverage, and bias
  
  # D for each d
  coverage_rmsted1_nocov <- matrix(NA, 5, 20)
  for (d in 1:5) {
    rmsteHPD_nocov_include1 <- matrix(NA, nsim, 20)
    j <- NULL
    for (i in s1_nocov) {
      j <- sum(j, 1)
      rmsteHPD_nocov_include1[j, ] <- sim_SCENARIO1_nocov[[i]]$rmsted_HPD_nocov_include[d, ]
    }
    coverage_rmsted1_nocov[d, ] <- colMeans(rmsteHPD_nocov_include1)  # coverage
  }
  
  avg_width_rmsted1_nocov <- matrix(NA, 5, 20)
  for (d in 1:5) {
    rmsteHPD_width1_nocov <- matrix(NA, nsim, 20)
    j <- NULL
    for (i in s1_nocov) {
      j <- sum(j, 1)
      rmsteHPD_width1_nocov[j, ] <- sim_SCENARIO1_nocov[[i]]$width_hdi_rmsted_nocov[d, ]
    }
    avg_width_rmsted1_nocov[d, ] <- colMeans(rmsteHPD_width1_nocov)  # width
  }
  
  rmste.d_mean_bias1_nocov <- matrix(NA, 5, 20)
  for (d in 1:5) {
    rmste.d_bias1_nocov <- matrix(NA, nrow = nsim, ncol = 20)
    j <- NULL
    for (i in s1_nocov) {
      j <- sum(j, 1)
      rmste.d_bias1_nocov[j, ] <- as.vector(sim_SCENARIO1_nocov[[i]]$rmsted_nocov_bias[d, ])
    }
    rmste.d_mean_bias1_nocov[d, ] <- colMeans(rmste.d_bias1_nocov)  # bias
  }
  
  # average D
  rmsteDHPD_nocov_include1 <- matrix(NA, nsim, 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    rmsteDHPD_nocov_include1[j, ] <- sim_SCENARIO1_nocov[[i]]$rmsteD_HPD_nocov_include
  }
  coverage_rmsteD1_nocov <- colMeans(rmsteDHPD_nocov_include1)  # coverage
  
  rmsteDHPD_width1_nocov <- matrix(NA, nsim, 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    rmsteDHPD_width1_nocov[j, ] <- sim_SCENARIO1_nocov[[i]]$width_hdi_rmsteD_nocov
  }
  avg_width_rmsteD1_nocov <- colMeans(rmsteDHPD_width1_nocov)  # width
  
  rmste.D_bias1_nocov <- matrix(NA, nrow = nsim, ncol = 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    rmste.D_bias1_nocov[j, ] <- as.vector(sim_SCENARIO1_nocov[[i]]$rmsteD_nocov_bias)
  }
  rmste.D_mean_bias1_nocov <- colMeans(rmste.D_bias1_nocov)  # bias
  
  # ND
  rmstendHPD_nocov_include1 <- matrix(NA, nsim, 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    rmstendHPD_nocov_include1[j, ] <- sim_SCENARIO1_nocov[[i]]$rmstend_HPD_nocov_include
  }
  coverage_rmstend1_nocov <- colMeans(rmstendHPD_nocov_include1)  # coverage
  
  rmstendHPD_width1_nocov <- matrix(NA, nsim, 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    rmstendHPD_width1_nocov[j, ] <- sim_SCENARIO1_nocov[[i]]$width_hdi_rmstend_nocov
  }
  avg_width_rmstend1_nocov <- colMeans(rmstendHPD_width1_nocov)  # width
  
  rmste.nd_bias1_nocov <- matrix(NA, nrow = nsim, ncol = 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    rmste.nd_bias1_nocov[j, ] <- as.vector(sim_SCENARIO1_nocov[[i]]$rmstend_nocov_bias)
  }
  rmste.nd_mean_bias1_nocov <- colMeans(rmste.nd_bias1_nocov)  # bias
  
  # ITT
  rmsteittHPD_nocov_include1 <- matrix(NA, nsim, 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    rmsteittHPD_nocov_include1[j, ] <- sim_SCENARIO1_nocov[[i]]$rmsteitt_HPD_nocov_include
  }
  coverage_rmsteitt1_nocov <- colMeans(rmsteittHPD_nocov_include1)  # coverage
  
  rmsteittHPD_width1_nocov <- matrix(NA, nsim, 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    rmsteittHPD_width1_nocov[j, ] <- sim_SCENARIO1_nocov[[i]]$width_hdi_rmsteitt_nocov
  }
  avg_width_rmsteitt1_nocov <- colMeans(rmsteittHPD_width1_nocov)  # width
  
  rmsteitt_bias1_nocov <- matrix(NA, nrow = nsim, ncol = 20)
  j <- NULL
  for (i in s1_nocov) {
    j <- sum(j, 1)
    rmsteitt_bias1_nocov[j, ] <- as.vector(sim_SCENARIO1_nocov[[i]]$rmsteitt_nocov_bias)
  }
  rmsteitt_mean_bias1_nocov <- colMeans(rmsteitt_bias1_nocov)  # bias
  
  # SCENARIO 2 w/o covariates
  
  n2_nocov <- length(s2_nocov)
  
  ### PROPORTION OF ND patients
  ID_save2_nocov <- matrix(NA, n2_nocov, 1)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    ID_save2_nocov[j] <- sim_SCENARIO2_nocov[[i]]$ID_save_nocov
  }
  pc_ID2_nocov <- mean(ID_save2_nocov)
  pc_ID2_nocov
  
  ### DISTRIBUTIONAL CAUSAL EFFECTS: CI width, coverage, and bias
  
  # D for each d
  coverage_dced2_nocov <- matrix(NA, 5, 20)
  for (d in 1:5) {
    dceHPD_nocov_include2 <- matrix(NA, n2_nocov, 20)
    j <- NULL
    for (i in s2_nocov) {
      j <- sum(j, 1)
      dceHPD_nocov_include2[j, ] <- sim_SCENARIO2_nocov[[i]]$dced_HPD_nocov_include[d, ]
    }
    coverage_dced2_nocov[d, ] <- colMeans(dceHPD_nocov_include2)  # coverage
  }
  
  avg_width_dced2_nocov <- matrix(NA, 5, 20)
  for (d in 1:5) {
    dceHPD_width2_nocov <- matrix(NA, n2_nocov, 20)
    j <- NULL
    for (i in s2_nocov) {
      j <- sum(j, 1)
      dceHPD_width2_nocov[j, ] <- sim_SCENARIO2_nocov[[i]]$width_hdi_dced_nocov[d, ]
    }
    avg_width_dced2_nocov[d, ] <- colMeans(dceHPD_width2_nocov)  # width
  }
  
  dce.d_mean_bias2_nocov <- matrix(NA, 5, 20)
  for (d in 1:5) {
    dce.d_bias2_nocov <- matrix(NA, nrow = n2_nocov, ncol = 20)
    j <- NULL
    for (i in s2_nocov) {
      j <- sum(j, 1)
      dce.d_bias2_nocov[j, ] <- as.vector(sim_SCENARIO2_nocov[[i]]$dced_nocov_bias[d, ])
    }
    dce.d_mean_bias2_nocov[d, ] <- colMeans(dce.d_bias2_nocov)  # bias
  }
  
  # average D
  dceDHPD_nocov_include2 <- matrix(NA, n2_nocov, 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    dceDHPD_nocov_include2[j, ] <- sim_SCENARIO2_nocov[[i]]$dceD_HPD_nocov_include
  }
  coverage_dceD2_nocov <- colMeans(dceDHPD_nocov_include2)  # coverage
  
  dceDHPD_width2_nocov <- matrix(NA, n2_nocov, 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    dceDHPD_width2_nocov[j, ] <- sim_SCENARIO2_nocov[[i]]$width_hdi_dceD_nocov
  }
  avg_width_dceD2_nocov <- colMeans(dceDHPD_width2_nocov)  # width
  
  dce.D_bias2_nocov <- matrix(NA, nrow = n2_nocov, ncol = 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    dce.D_bias2_nocov[j, ] <- as.vector(sim_SCENARIO2_nocov[[i]]$dceD_nocov_bias)
  }
  dce.D_mean_bias2_nocov <- colMeans(dce.D_bias2_nocov)  # bias 
  
  # ND
  dcendHPD_nocov_include2 <- matrix(NA, n2_nocov, 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    dcendHPD_nocov_include2[j, ] <- sim_SCENARIO2_nocov[[i]]$dcend_HPD_nocov_include
  }
  coverage_dcend2_nocov <- colMeans(dcendHPD_nocov_include2)  # coverage
  
  dcendHPD_width2_nocov <- matrix(NA, n2_nocov, 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    dcendHPD_width2_nocov[j, ] <- sim_SCENARIO2_nocov[[i]]$width_hdi_dcend_nocov
  }
  avg_width_dcend2_nocov <- colMeans(dcendHPD_width2_nocov)  # width
  
  dce.nd_bias2_nocov <- matrix(NA, nrow = n2_nocov, ncol = 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    dce.nd_bias2_nocov[j, ] <- as.vector(sim_SCENARIO2_nocov[[i]]$dcend_nocov_bias)
  }
  dce.nd_mean_bias2_nocov <- colMeans(dce.nd_bias2_nocov)  # bias
  
  # ITT
  dceittHPD_nocov_include2 <- matrix(NA, n2_nocov, 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    dceittHPD_nocov_include2[j, ] <- sim_SCENARIO2_nocov[[i]]$dceitt_HPD_nocov_include
  }
  coverage_dceitt2_nocov <- colMeans(dceittHPD_nocov_include2)  # coverage
  
  dceittHPD_width2_nocov <- matrix(NA, n2_nocov, 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    dceittHPD_width2_nocov[j, ] <- sim_SCENARIO2_nocov[[i]]$width_hdi_dceitt_nocov
  }
  avg_width_dceitt2_nocov <- colMeans(dceittHPD_width2_nocov)  # width 
  
  dceitt_bias2_nocov <- matrix(NA, nrow = n2_nocov, ncol = 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    dceitt_bias2_nocov[j, ] <- as.vector(sim_SCENARIO2_nocov[[i]]$dceitt_nocov_bias)
  }
  dceitt_mean_bias2_nocov <- colMeans(dceitt_bias2_nocov)  # bias
  
  ### RESTRICTED MEAN SURVIVAL TIME: CI width, coverage, and bias
  
  # D for each d
  coverage_rmsted2_nocov <- matrix(NA, 5, 20)
  for (d in 1:5) {
    rmsteHPD_nocov_include2 <- matrix(NA, n2_nocov, 20)
    j <- NULL
    for (i in s2_nocov) {
      j <- sum(j, 1)
      rmsteHPD_nocov_include2[j, ] <- sim_SCENARIO2_nocov[[i]]$rmsted_HPD_nocov_include[d, ]
    }
    coverage_rmsted2_nocov[d, ] <- colMeans(rmsteHPD_nocov_include2)  # coverage
  }
  
  avg_width_rmsted2_nocov <- matrix(NA, 5, 20)
  for (d in 1:5) {
    rmsteHPD_width2_nocov <- matrix(NA, n2_nocov, 20)
    j <- NULL
    for (i in s2_nocov) {
      j <- sum(j, 1)
      rmsteHPD_width2_nocov[j, ] <- sim_SCENARIO2_nocov[[i]]$width_hdi_rmsted_nocov[d, ]
    }
    avg_width_rmsted2_nocov[d, ] <- colMeans(rmsteHPD_width2_nocov)  # width
  }
  
  rmste.d_mean_bias2_nocov <- matrix(NA, 5, 20)
  for (d in 1:5) {
    rmste.d_bias2_nocov <- matrix(NA, nrow = n2_nocov, ncol = 20)
    j <- NULL
    for (i in s2_nocov) {
      j <- sum(j, 1)
      rmste.d_bias2_nocov[j, ] <- as.vector(sim_SCENARIO2_nocov[[i]]$rmsted_nocov_bias[d, ])
    }
    rmste.d_mean_bias2_nocov[d, ] <- colMeans(rmste.d_bias2_nocov)  # bias
  }
  
  # average D
  rmsteDHPD_nocov_include2 <- matrix(NA, n2_nocov, 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    rmsteDHPD_nocov_include2[j, ] <- sim_SCENARIO2_nocov[[i]]$rmsteD_HPD_nocov_include
  }
  coverage_rmsteD2_nocov <- colMeans(rmsteDHPD_nocov_include2)  # coverage
  
  rmsteDHPD_width2_nocov <- matrix(NA, n2_nocov, 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    rmsteDHPD_width2_nocov[j, ] <- sim_SCENARIO2_nocov[[i]]$width_hdi_rmsteD_nocov
  }
  avg_width_rmsteD2_nocov <- colMeans(rmsteDHPD_width2_nocov)  # width
  
  rmste.D_bias2_nocov <- matrix(NA, nrow = n2_nocov, ncol = 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    rmste.D_bias2_nocov[j, ] <- as.vector(sim_SCENARIO2_nocov[[i]]$rmsteD_nocov_bias)
  }
  rmste.D_mean_bias2_nocov <- colMeans(rmste.D_bias2_nocov)  # bias
  
  # ND
  rmstendHPD_nocov_include2 <- matrix(NA, n2_nocov, 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    rmstendHPD_nocov_include2[j, ] <- sim_SCENARIO2_nocov[[i]]$rmstend_HPD_nocov_include
  }
  coverage_rmstend2_nocov <- colMeans(rmstendHPD_nocov_include2)  # coverage
  
  rmstendHPD_width2_nocov <- matrix(NA, n2_nocov, 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    rmstendHPD_width2_nocov[j, ] <- sim_SCENARIO2_nocov[[i]]$width_hdi_rmstend_nocov
  }
  avg_width_rmstend2_nocov <- colMeans(rmstendHPD_width2_nocov)  # width
  
  rmste.nd_bias2_nocov <- matrix(NA, nrow = n2_nocov, ncol = 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    rmste.nd_bias2_nocov[j, ] <- as.vector(sim_SCENARIO2_nocov[[i]]$rmstend_nocov_bias)
  }
  rmste.nd_mean_bias2_nocov <- colMeans(rmste.nd_bias2_nocov)  # bias
  
  # ITT
  rmsteittHPD_nocov_include2 <- matrix(NA, n2_nocov, 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    rmsteittHPD_nocov_include2[j, ] <- sim_SCENARIO2_nocov[[i]]$rmsteitt_HPD_nocov_include
  }
  coverage_rmsteitt2_nocov <- colMeans(rmsteittHPD_nocov_include2)  # coverage
  
  rmsteittHPD_width2_nocov <- matrix(NA, n2_nocov, 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    rmsteittHPD_width2_nocov[j, ] <- sim_SCENARIO2_nocov[[i]]$width_hdi_rmsteitt_nocov
  }
  avg_width_rmsteitt2_nocov <- colMeans(rmsteittHPD_width2_nocov)  # width
  
  rmsteitt_bias2_nocov <- matrix(NA, nrow = n2_nocov, ncol = 20)
  j <- NULL
  for (i in s2_nocov) {
    j <- sum(j, 1)
    rmsteitt_bias2_nocov[j, ] <- as.vector(sim_SCENARIO2_nocov[[i]]$rmsteitt_nocov_bias)
  }
  rmsteitt_mean_bias2_nocov <- colMeans(rmsteitt_bias2_nocov)  # bias
  
  # TABLES DCE - coverage
  DCE_coverage_scenario1_nocov <- data.frame(DCE_ITT = as.matrix(round(coverage_dceitt1_nocov, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(coverage_dcend1_nocov,
                                                                                                                                                  2), nrow = 20, ncol = 1), DCE_D = as.matrix(round(coverage_dceD1_nocov, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(coverage_dced1_nocov[1, ], 2), nrow = 20,
                                                                                                                                                                                                                                                                       ncol = 1), DCE_d2 = as.matrix(round(coverage_dced1_nocov[2, ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(coverage_dced1_nocov[3, ], 2), nrow = 20, ncol = 1),
                                             DCE_d4 = as.matrix(round(coverage_dced1_nocov[4, ], 2), nrow = 20, ncol = 1), DCE_d5 = as.matrix(round(coverage_dced1_nocov[5, ], 2), nrow = 20, ncol = 1))
  
  DCE_coverage_scenario2_nocov <- data.frame(DCE_ITT = as.matrix(round(coverage_dceitt2_nocov, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(coverage_dcend2_nocov,
                                                                                                                                                  2), nrow = 20, ncol = 1), DCE_D = as.matrix(round(coverage_dceD2_nocov, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(coverage_dced2_nocov[1, ], 2), nrow = 20,
                                                                                                                                                                                                                                                                       ncol = 1), DCE_d2 = as.matrix(round(coverage_dced2_nocov[2, ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(coverage_dced2_nocov[3, ], 2), nrow = 20, ncol = 1),
                                             DCE_d4 = as.matrix(round(coverage_dced2_nocov[4, ], 2), nrow = 20, ncol = 1), DCE_d5 = as.matrix(round(coverage_dced2_nocov[5, ], 2), nrow = 20, ncol = 1))
  
  # DCE - width
  DCE_avg_width_scenario1_nocov <- data.frame(DCE_ITT = as.matrix(round(avg_width_dceitt1_nocov, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(avg_width_dcend1_nocov,
                                                                                                                                                    2), nrow = 20, ncol = 1), DCE_D = as.matrix(round(avg_width_dceD1_nocov, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(avg_width_dced1_nocov[1, ], 2), nrow = 20,
                                                                                                                                                                                                                                                                          ncol = 1), DCE_d2 = as.matrix(round(avg_width_dced1_nocov[2, ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(avg_width_dced1_nocov[3, ], 2), nrow = 20,
                                                                                                                                                                                                                                                                                                                                                                                       ncol = 1), DCE_d4 = as.matrix(round(avg_width_dced1_nocov[4, ], 2), nrow = 20, ncol = 1), DCE_d5 = as.matrix(round(avg_width_dced1_nocov[5, ], 2), nrow = 20,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    ncol = 1))
  
  DCE_avg_width_scenario2_nocov <- data.frame(DCE_ITT = as.matrix(round(avg_width_dceitt2_nocov, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(avg_width_dcend2_nocov,
                                                                                                                                                    2), nrow = 20, ncol = 1), DCE_D = as.matrix(round(avg_width_dceD2_nocov, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(avg_width_dced2_nocov[1, ], 2), nrow = 20,
                                                                                                                                                                                                                                                                          ncol = 1), DCE_d2 = as.matrix(round(avg_width_dced2_nocov[2, ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(avg_width_dced2_nocov[3, ], 2), nrow = 20,
                                                                                                                                                                                                                                                                                                                                                                                       ncol = 1), DCE_d4 = as.matrix(round(avg_width_dced2_nocov[4, ], 2), nrow = 20, ncol = 1), DCE_d5 = as.matrix(round(avg_width_dced2_nocov[5, ], 2), nrow = 20,
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    ncol = 1))
  
  # DCE - bias
  DCE_bias_scenario1_nocov <- data.frame(DCE_ITT = as.matrix(round(dceitt_mean_bias1_nocov, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(dce.nd_mean_bias1_nocov,
                                                                                                                                               2), nrow = 20, ncol = 1), DCE_D = as.matrix(round(dce.D_mean_bias1_nocov, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(dce.d_mean_bias1_nocov[1, ], 2),
                                                                                                                                                                                                                                                                      nrow = 20, ncol = 1), DCE_d2 = as.matrix(round(dce.d_mean_bias1_nocov[2, ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(dce.d_mean_bias1_nocov[3, ], 2),
                                                                                                                                                                                                                                                                                                                                                                                               nrow = 20, ncol = 1), DCE_d4 = as.matrix(round(dce.d_mean_bias1_nocov[4, ], 2), nrow = 20, ncol = 1), DCE_d5 = as.matrix(round(dce.d_mean_bias1_nocov[5, ], 2),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        nrow = 20, ncol = 1))
  
  DCE_bias_scenario2_nocov <- data.frame(DCE_ITT = as.matrix(round(dceitt_mean_bias2_nocov, 2), nrow = 20, ncol = 1), DCE_ND = as.matrix(round(dce.nd_mean_bias2_nocov,
                                                                                                                                               2), nrow = 20, ncol = 1), DCE_D = as.matrix(round(dce.D_mean_bias2_nocov, 2), nrow = 20, ncol = 1), DCE_d1 = as.matrix(round(dce.d_mean_bias2_nocov[1, ], 2),
                                                                                                                                                                                                                                                                      nrow = 20, ncol = 1), DCE_d2 = as.matrix(round(dce.d_mean_bias2_nocov[2, ], 2), nrow = 20, ncol = 1), DCE_d3 = as.matrix(round(dce.d_mean_bias2_nocov[3, ], 2),
                                                                                                                                                                                                                                                                                                                                                                                               nrow = 20, ncol = 1), DCE_d4 = as.matrix(round(dce.d_mean_bias2_nocov[4, ], 2), nrow = 20, ncol = 1), DCE_d5 = as.matrix(round(dce.d_mean_bias2_nocov[5, ], 2),
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        nrow = 20, ncol = 1))
  
  # RMSTE - coverage
  RMSTE_coverage_scenario1_nocov <- data.frame(RMSTE_ITT = as.matrix(round(coverage_rmsteitt1_nocov, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(coverage_rmstend1_nocov,
                                                                                                                                                          2), nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(coverage_rmsteD1_nocov, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(coverage_rmsted1_nocov[1, ],
                                                                                                                                                                                                                                                                                           2), nrow = 20, ncol = 1), RMSTE_d2 = as.matrix(round(coverage_rmsted1_nocov[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(coverage_rmsted1_nocov[3,
                                                                                                                                                                                                                                                                                           ], 2), nrow = 20, ncol = 1), RMSTE_d4 = as.matrix(round(coverage_rmsted1_nocov[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(coverage_rmsted1_nocov[5,
                                                                                                                                                                                                                                                                                           ], 2), nrow = 20, ncol = 1))
  
  RMSTE_coverage_scenario2_nocov <- data.frame(RMSTE_ITT = as.matrix(round(coverage_rmsteitt2_nocov, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(coverage_rmstend2_nocov,
                                                                                                                                                          2), nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(coverage_rmsteD2_nocov, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(coverage_rmsted2_nocov[1, ],
                                                                                                                                                                                                                                                                                           2), nrow = 20, ncol = 1), RMSTE_d2 = as.matrix(round(coverage_rmsted2_nocov[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(coverage_rmsted2_nocov[3,
                                                                                                                                                                                                                                                                                           ], 2), nrow = 20, ncol = 1), RMSTE_d4 = as.matrix(round(coverage_rmsted2_nocov[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(coverage_rmsted2_nocov[5,
                                                                                                                                                                                                                                                                                           ], 2), nrow = 20, ncol = 1))
  
  # RMSTE - width
  RMSTE_avg_width_scenario1_nocov <- data.frame(RMSTE_ITT = as.matrix(round(avg_width_rmsteitt1_nocov, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(avg_width_rmstend1_nocov,
                                                                                                                                                            2), nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(avg_width_rmsteD1_nocov, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(avg_width_rmsted1_nocov[1,
                                                                                                                                                            ], 2), nrow = 20, ncol = 1), RMSTE_d2 = as.matrix(round(avg_width_rmsted1_nocov[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(avg_width_rmsted1_nocov[3,
                                                                                                                                                            ], 2), nrow = 20, ncol = 1), RMSTE_d4 = as.matrix(round(avg_width_rmsted1_nocov[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(avg_width_rmsted1_nocov[5,
                                                                                                                                                            ], 2), nrow = 20, ncol = 1))
  
  RMSTE_avg_width_scenario2_nocov <- data.frame(RMSTE_ITT = as.matrix(round(avg_width_rmsteitt2_nocov, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(avg_width_rmstend2_nocov,
                                                                                                                                                            2), nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(avg_width_rmsteD2_nocov, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(avg_width_rmsted2_nocov[1,
                                                                                                                                                            ], 2), nrow = 20, ncol = 1), RMSTE_d2 = as.matrix(round(avg_width_rmsted2_nocov[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(avg_width_rmsted2_nocov[3,
                                                                                                                                                            ], 2), nrow = 20, ncol = 1), RMSTE_d4 = as.matrix(round(avg_width_rmsted2_nocov[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(avg_width_rmsted2_nocov[5,
                                                                                                                                                            ], 2), nrow = 20, ncol = 1))
  
  # RMSTE - bias
  RMSTE_bias_scenario1_nocov <- data.frame(RMSTE_ITT = as.matrix(round(rmsteitt_mean_bias1_nocov, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(rmste.nd_mean_bias1_nocov,
                                                                                                                                                       2), nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(rmste.D_mean_bias1_nocov, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(rmste.d_mean_bias1_nocov[1,
                                                                                                                                                       ], 2), nrow = 20, ncol = 1), RMSTE_d2 = as.matrix(round(rmste.d_mean_bias1_nocov[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(rmste.d_mean_bias1_nocov[3,
                                                                                                                                                       ], 2), nrow = 20, ncol = 1), RMSTE_d4 = as.matrix(round(rmste.d_mean_bias1_nocov[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(rmste.d_mean_bias1_nocov[5,
                                                                                                                                                       ], 2), nrow = 20, ncol = 1))
  
  RMSTE_bias_scenario2_nocov <- data.frame(RMSTE_ITT = as.matrix(round(rmsteitt_mean_bias2_nocov, 2), nrow = 20, ncol = 1), RMSTE_ND = as.matrix(round(rmste.nd_mean_bias2_nocov,
                                                                                                                                                       2), nrow = 20, ncol = 1), RMSTE_D = as.matrix(round(rmste.D_mean_bias2_nocov, 2), nrow = 20, ncol = 1), RMSTE_d1 = as.matrix(round(rmste.d_mean_bias2_nocov[1,
                                                                                                                                                       ], 2), nrow = 20, ncol = 1), RMSTE_d2 = as.matrix(round(rmste.d_mean_bias2_nocov[2, ], 2), nrow = 20, ncol = 1), RMSTE_d3 = as.matrix(round(rmste.d_mean_bias2_nocov[3,
                                                                                                                                                       ], 2), nrow = 20, ncol = 1), RMSTE_d4 = as.matrix(round(rmste.d_mean_bias2_nocov[4, ], 2), nrow = 20, ncol = 1), RMSTE_d5 = as.matrix(round(rmste.d_mean_bias2_nocov[5,
                                                                                                                                                       ], 2), nrow = 20, ncol = 1))
  Table6 <- DCE_coverage_scenario1_nocov
  Table7 <- DCE_coverage_scenario2_nocov
  
  saveRDS(Table6, file = "tables/INTERMEDIATE/Table6.rds")
  saveRDS(Table7, file = "tables/INTERMEDIATE/Table7.rds")

  Table10 <- DCE_avg_width_scenario1_nocov
  Table11 <- DCE_avg_width_scenario2_nocov

  saveRDS(Table10, file = "tables/INTERMEDIATE/Table10.rds")
  saveRDS(Table11, file = "tables/INTERMEDIATE/Table11.rds")

  Table14 <- DCE_bias_scenario1_nocov
  Table15 <- DCE_bias_scenario2_nocov
  
  saveRDS(Table14, file = "tables/INTERMEDIATE/Table14.rds")
  saveRDS(Table15, file = "tables/INTERMEDIATE/Table15.rds")
  
}

