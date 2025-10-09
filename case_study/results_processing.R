################################################################################ 
# Evaluating causal effects on time-to-event outcomes in an RCT in Oncology 
# with treatment discontinuation
#
# Authors: V. Ballerini, B., Bornkamp, F. Mealli, C. Wang, Y. Zhang, A. Mattei
#
# Replication Package for results in Section 6 
# In particular: Posterior Predictive checks, posterior values of the parameters,
# proportion of ND patients.
#
# Code authors: Veronica Ballerini
# Last modified: September 27, 2025 
################################################################################ 

# Load results (if you started a new session)
# and initial dataset
load("final_casestudy.RData")
observed.data <- read.csv("synthetic_data.csv")

# Load functions and libraries
source("ce_apply_expexp.R")
source("apply_cov.R")
source("Functions.R")
source("RepData.R")

# set the number of iterations, burnin and
# thinning
niter <- 50000
burn <- 30000
thin <- 10

THETA_expexp <- rbind(chain_expexp$Theta[seq(burn,
    niter, by = thin), ])

# Parameters' mean
t_expexp_mean <- as.matrix(apply(THETA_expexp, 2,
    mean), ncol = 1)
t_expexp_hdi <- as.matrix(t(apply(THETA_expexp,
    2, hdi)))
t_expexp <- cbind(t_expexp_mean, t_expexp_hdi)
colnames(t_expexp) <- c("posterior_mean", "HDI_lower",
    "HDI_upper")
saveRDS(t_expexp, file = "in_text_results/t_expexp.rds")

# Parameters' trace plots
for (i in 1:dim(THETA_expexp)[2]) {
    jpeg(paste0("figures/trace/ExpExp/", colnames(THETA_expexp)[i],
        ".jpeg"), width = 3500, height = 2500, units = "px",
        res = 300)
    plot(THETA_expexp[, i], type = "l", xlab = "iter",
        ylab = paste(colnames(THETA_expexp)[i]))
    dev.off()
}

# Store principal strata membership posterior
# in lists
DID_expexp <- array(NA, dim = c(dim(chain_expexp$DID[[1]])[1],
    dim(chain_expexp$DID[[1]])[2], dim(THETA_expexp)[1]))
j <- NULL
for (i in seq(burn, niter, by = thin)) {
    j <- sum(j, 1)
    DID_expexp[, , j] <- chain_expexp$DID[[i]]
}

# Sore parameters' posterior in lists
THETA_list_expexp <- list()
for (i in 1:dim(THETA_expexp)[1]) {
    THETA_list_expexp[[i]] <- list(eta.0 = THETA_expexp[i,
        "eta.0"], eta.1 = THETA_expexp[i, "eta.1"],
        eta.2 = THETA_expexp[i, "eta.2"], eta.3 = THETA_expexp[i,
            "eta.3"], beta.d = THETA_expexp[i, "beta.d"],
        eta.d1 = THETA_expexp[i, "eta.d1"], eta.d2 = THETA_expexp[i,
            "eta.d2"], eta.d3 = THETA_expexp[i,
            "eta.d3"], beta.y1nd = THETA_expexp[i,
            "beta.y1nd"], beta.y1d = THETA_expexp[i,
            "beta.y1d"], beta.y0nd = THETA_expexp[i,
            "beta.y0nd"], beta.y0d = THETA_expexp[i,
            "beta.y0d"], eta.y1 = THETA_expexp[i,
            "eta.y1"], eta.y2 = THETA_expexp[i,
            "eta.y2"], eta.y3 = THETA_expexp[i,
            "eta.y3"], delta = THETA_expexp[i, "delta"],
        D = DID_expexp[, , i][, 1], ID = DID_expexp[,
            , i][, 2], data = observed.data)
}

## KM
km_post <- lapply(THETA_list_expexp, kaplanMeier_expexp)

# ppp
rep_data_expexp <- lapply(THETA_list_expexp, rep.data)
km_rep <- lapply(rep_data_expexp, kaplanMeier_expexp)

Diff.G.nd <- NULL
Diff.G.d <- NULL
GD.d1 <- NULL

for (i in 1:length(THETA_list_expexp)) {
    Diff.G.nd <- rbind(Diff.G.nd, km_rep[[i]]$Diff.G.nd -
        km_post[[i]]$Diff.G.nd > 0)
    Diff.G.d <- rbind(Diff.G.d, km_rep[[i]]$Diff.G.d -
        km_post[[i]]$Diff.G.d > 0)
    GD.d1 <- rbind(GD.d1, km_rep[[i]]$GD.d1 - km_post[[i]]$GD.d1 >
        0)
}

pppv.Diff.G.nd <- colMeans(Diff.G.nd)
pppv.Diff.G.d <- colMeans(Diff.G.d)
pppv.Diff.GD.d1 <- colMeans(GD.d1)

PPP <- cbind(mean(pppv.Diff.G.nd), mean(pppv.Diff.G.d),
    mean(pppv.Diff.GD.d1))
colnames(PPP) <- c("KMdm_NDpatients", "KMdm_Dpatients",
    "KMdm_discontinuation")
saveRDS(PPP, file = "in_text_results/PPP.rds")

### PROPORTION OF ND
ND_mat <- DID_expexp[, 2, ]
pc_ND <- mean(apply(ND_mat, 2, mean))
hdi_ND <- hdi(apply(ND_mat, 2, mean))

p_ND <- cbind(pc_ND, t(hdi_ND))
colnames(p_ND) <- c("prop_NDpatients", "HDI_lower",
    "HDI_upper")

saveRDS(p_ND, file = "in_text_results/p_ND.rds")

### DISTRIBUTIONAL CAUSAL EFFECTS & RESTRICTED MEAN SURVIVAL
rmste <- lapply(THETA_list_expexp, rmste_apply_expexp)

dce.itt <- matrix(NA, dim(THETA_expexp)[1], 20)
for (i in 1:dim(THETA_expexp)[1]) {
    dce.itt[i, ] <- rmste[[i]]$DCE
}

dce.nd <- matrix(NA, dim(THETA_expexp)[1], 20)
for (i in 1:dim(THETA_expexp)[1]) {
    dce.nd[i, ] <- rmste[[i]]$DCE_ND
}

dce.D <- matrix(NA, dim(THETA_expexp)[1], 20)
for (i in 1:dim(THETA_expexp)[1]) {
    dce.D[i, ] <- rmste[[i]]$DCE_D
}

dce.d_1 <- matrix(NA, dim(THETA_expexp)[1], 20)
for (i in 1:dim(THETA_expexp)[1]) {
    dce.d_1[i, ] <- rmste[[i]]$DCE_Y_D[1, ]
}

dce.d_2 <- matrix(NA, dim(THETA_expexp)[1], 20)
for (i in 1:dim(THETA_expexp)[1]) {
    dce.d_2[i, ] <- rmste[[i]]$DCE_Y_D[2, ]
}

dce.d_3 <- matrix(NA, dim(THETA_expexp)[1], 20)
for (i in 1:dim(THETA_expexp)[1]) {
    dce.d_3[i, ] <- rmste[[i]]$DCE_Y_D[3, ]
}

dce.d_4 <- matrix(NA, dim(THETA_expexp)[1], 20)
for (i in 1:dim(THETA_expexp)[1]) {
    dce.d_4[i, ] <- rmste[[i]]$DCE_Y_D[4, ]
}

rmste.itt <- matrix(NA, dim(THETA_expexp)[1], 20)
for (i in 1:dim(THETA_expexp)[1]) {
    rmste.itt[i, ] <- rmste[[i]]$RMSTE
}

rmste.nd <- matrix(NA, dim(THETA_expexp)[1], 20)
for (i in 1:dim(THETA_expexp)[1]) {
    rmste.nd[i, ] <- rmste[[i]]$RMSTE_ND
}

rmste.d_1 <- matrix(NA, dim(THETA_expexp)[1], 20)
for (i in 1:dim(THETA_expexp)[1]) {
    rmste.d_1[i, ] <- rmste[[i]]$RMSTE_d[1, ]
}

rmste.d_2 <- matrix(NA, dim(THETA_expexp)[1], 20)
for (i in 1:dim(THETA_expexp)[1]) {
    rmste.d_2[i, ] <- rmste[[i]]$RMSTE_d[2, ]
}

rmste.d_3 <- matrix(NA, dim(THETA_expexp)[1], 20)
for (i in 1:dim(THETA_expexp)[1]) {
    rmste.d_3[i, ] <- rmste[[i]]$RMSTE_d[3, ]
}

rmste.d_4 <- matrix(NA, dim(THETA_expexp)[1], 20)
for (i in 1:dim(THETA_expexp)[1]) {
    rmste.d_4[i, ] <- rmste[[i]]$RMSTE_d[4, ]
}

rmste.d_5 <- matrix(NA, dim(THETA_expexp)[1], 20)
for (i in 1:dim(THETA_expexp)[1]) {
    rmste.d_5[i, ] <- rmste[[i]]$RMSTE_d[5, ]
}

rmste.D <- matrix(NA, dim(THETA_expexp)[1], 20)
for (i in 1:dim(THETA_expexp)[1]) {
    rmste.D[i, ] <- rmste[[i]]$RMSTE_D
}

cov_distr <- lapply(THETA_list_expexp, apply_cov)

#### Save results
save(THETA_expexp, cov_distr, rmste, file = "results.RData")
