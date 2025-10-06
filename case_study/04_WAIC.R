################################################################################ 
# Evaluating causal effects on time-to-event outcomes in an RCT in Oncology 
# with treatment discontinuation
#
# Authors: V. Ballerini, B., Bornkamp, F. Mealli, C. Wang, Y. Zhang, A. Mattei
#
# Replication Package for results in Section 6 
# In particular: WAIC for model selection
#
# Code authors: Veronica Ballerini
# Last modified: October 6, 2025 
################################################################################ 

# Weibull-Weibull model: Weibull distributions for the potential time to 
# discontinuation and the potential PFS
source("MCMC_weibweib.R")
source("CompleteLogPost_weibweib.R")
source("DataAugmentation_weibweib.R")
source("MCMC_initialization_weibweib.R")

set.seed(474756)

chain_weibweib <- mcmc.tdiscontinue_weibweib(niter = niter, burn = burn, 
                                             thin = thin, 
                                             par.start = par.start_weibweib, 
                                             proposal = proposal_weibweib, 
                                             theta.prior = theta.prior_weibweib, 
                                             observed.data = observed.data)

# Exp-Weibull model: Exponential distribution for the potential time to 
# discontinuation and Weibull distributions for the potential PFS
source("MCMC_expweib.R")
source("CompleteLogPost_expweib.R")
source("DataAugmentation_expweib.R")
source("MCMC_initialization_expweib.R")

set.seed(474756)
chain_expweib <- mcmc.tdiscontinue_expweib(niter = niter, burn = burn, 
                                           thin = thin, 
                                           par.start = par.start_expweib, 
                                           proposal = proposal_expweib, 
                                           theta.prior = theta.prior_expweib, 
                                           observed.data = observed.data)
#### Save results
save(chain_weibweib, chain_expweib, "other_models.RData")

# Results
load("other_models.RData")

source("ce_apply_weibweib.R")
source("ce_apply_expweib.R")

THETA_weibweib <- rbind(chain_weibweib$Theta[seq(burn, niter, by = thin), ])
THETA_expweib <- rbind(chain_expweib$Theta[seq(burn, niter, by = thin), ])

# Parameters' mean
t_weibweib_mean<-as.matrix(apply(THETA_weibweib, 2, mean),ncol=1)
t_weibweib_hdi<-as.matrix(t(apply(THETA_weibweib, 2, hdi)))
t_weibweib<-cbind(t_weibweib_mean,t_weibweib_hdi)
colnames(t_weibweib)<-c("posterior_mean","HDI_lower","HDI_upper")
saveRDS(t_weibweib, file = "tables/extra/t_weibweib.rds")

t_expweib_mean<-as.matrix(apply(THETA_expweib, 2, mean),ncol=1)
t_expweib_hdi<-as.matrix(t(apply(THETA_expweib, 2, hdi)))
t_expweib<-cbind(t_expweib_mean,t_expweib_hdi)
colnames(t_expweib)<-c("posterior_mean","HDI_lower","HDI_upper")
saveRDS(t_expweib, file = "tables/extra/t_expweib.rds")

saveRDS(t_expexp, file = "tables/extra/t_expexp.rds")

# Parameters' trace plots
for (i in 1:dim(THETA_weibweib)[2]) {
  jpeg(paste0("figures/trace/WeibWeib/", colnames(THETA_weibweib)[i], ".jpeg"), width = 3500, height = 2500, units = "px", res = 300)
  plot(THETA_weibweib[, i], type = "l", xlab = "iter", 
       ylab = paste(colnames(THETA_weibweib)[i]))
  dev.off()
}
for (i in 1:dim(THETA_expweib)[2]) {
  jpeg(paste0("figures/trace/ExpWeib/", colnames(THETA_weibweib)[i], ".jpeg"), width = 3500, height = 2500, units = "px", res = 300)
  plot(THETA_expweib[, i], type = "l", xlab = "iter", 
       ylab = paste(colnames(THETA_expweib)[i]))
  dev.off()
}

# Store principal strata membership posterior in lists
DID_weibweib <- array(NA, dim = c(dim(chain_weibweib$DID[[1]])[1], 
                                  dim(chain_weibweib$DID[[1]])[2], 
                                  dim(THETA_weibweib)[1]))
DID_expweib <- array(NA, dim = c(dim(chain_expweib$DID[[1]])[1], 
                                 dim(chain_expweib$DID[[1]])[2], 
                                 dim(THETA_expweib)[1]))

j <- NULL
for (i in seq(burn, niter, by = thin)) {
  j <- sum(j, 1)
  DID_weibweib[, , j] <- chain_weibweib$DID[[i]]
  DID_expweib[, , j] <- chain_expweib$DID[[i]]
}

# Sore parameters' posterior in lists
THETA_list_weibweib <- list()
for (i in 1:dim(THETA_weibweib)[1]) {
  THETA_list_weibweib[[i]] <- list(eta.0 = THETA_weibweib[i, "eta.0"], 
                                   eta.1 = THETA_weibweib[i, "eta.1"], 
                                   eta.2 = THETA_weibweib[i, "eta.2"], 
                                   eta.3 = THETA_weibweib[i, "eta.3"], 
                                   alpha.d = THETA_weibweib[i, "alpha.d"], 
                                   beta.d = THETA_weibweib[i, "beta.d"], 
                                   eta.d1 = THETA_weibweib[i, "eta.d1"], 
                                   eta.d2 = THETA_weibweib[i, "eta.d2"], 
                                   eta.d3 = THETA_weibweib[i, "eta.d3"], 
                                   alpha.y1nd = THETA_weibweib[i, "alpha.y1nd"],
                                   beta.y1nd = THETA_weibweib[i, "beta.y1nd"], 
                                   alpha.y1d = THETA_weibweib[i, "alpha.y1d"], 
                                   beta.y1d = THETA_weibweib[i, "beta.y1d"], 
                                   alpha.y0nd = THETA_weibweib[i, "alpha.y0nd"], 
                                   beta.y0nd = THETA_weibweib[i, "beta.y0nd"], 
                                   alpha.y0d = THETA_weibweib[i, "alpha.y0d"], 
                                   beta.y0d = THETA_weibweib[i, "beta.y0d"], 
                                   eta.y1 = THETA_weibweib[i, "eta.y1"], 
                                   eta.y2 = THETA_weibweib[i, "eta.y2"], 
                                   eta.y3 = THETA_weibweib[i, "eta.y3"], 
                                   delta = THETA_weibweib[i, "delta"], 
                                   D = DID_weibweib[, , i][, 1], 
                                   ID = DID_weibweib[, , i][, 2], 
                                   data = observed.data)
}

THETA_list_expweib <- list()
for (i in 1:dim(THETA_expweib)[1]) {
  THETA_list_expweib[[i]] <- list(eta.0 = THETA_expweib[i, "eta.0"], 
                                  eta.1 = THETA_expweib[i, "eta.1"], 
                                  eta.2 = THETA_expweib[i, "eta.2"], 
                                  eta.3 = THETA_expweib[i, "eta.3"], 
                                  beta.d = THETA_expweib[i, "beta.d"], 
                                  eta.d1 = THETA_expweib[i, "eta.d1"], 
                                  eta.d2 = THETA_expweib[i, "eta.d2"], 
                                  eta.d3 = THETA_expweib[i, "eta.d3"], 
                                  alpha.y1nd = THETA_expweib[i, "alpha.y1nd"], 
                                  beta.y1nd = THETA_expweib[i, "beta.y1nd"], 
                                  alpha.y1d = THETA_expweib[i, "alpha.y1d"], 
                                  beta.y1d = THETA_expweib[i, "beta.y1d"], 
                                  alpha.y0nd = THETA_expweib[i, "alpha.y0d"], 
                                  beta.y0nd = THETA_expweib[i, "beta.y0nd"], 
                                  alpha.y0d = THETA_expweib[i, "alpha.y0d"], 
                                  beta.y0d = THETA_expweib[i, "beta.y0d"], 
                                  eta.y1 = THETA_expweib[i, "eta.y1"], 
                                  eta.y2 = THETA_expweib[i, "eta.y2"], 
                                  eta.y3 = THETA_expweib[i, "eta.y3"], 
                                  delta = THETA_expweib[i, "delta"], 
                                  D = DID_expweib[, , i][, 1], 
                                  ID = DID_expweib[, , i][, 2], 
                                  data = observed.data)
}

# WAIC computation
ll_weibweib <- lapply(THETA_list_weibweib, waicfun_weibweib)
ll_mat_weibweib <- matrix(NA, nrow = dim(observed.data)[1], 
                          ncol = length(THETA_list_weibweib))
for (s in 1:length(THETA_list_weibweib)) {
  ll_mat_weibweib[, s] <- ll_weibweib[[s]]
}

ll_expweib <- lapply(THETA_list_expweib, waicfun_expweib)
ll_mat_expweib <- matrix(NA, nrow = dim(observed.data)[1], 
                         ncol = length(THETA_list_expweib))
for (s in 1:length(THETA_list_expweib)) {
  ll_mat_expweib[, s] <- ll_expweib[[s]]
}

ll_expexp <- lapply(THETA_list_expexp, waicfun_expexp)
ll_mat_expexp <- matrix(NA, nrow = dim(observed.data)[1], 
                        ncol = length(THETA_list_expexp))
for (s in 1:length(THETA_list_expexp)) {
  ll_mat_expexp[, s] <- ll_expexp[[s]]
}

# Log pointwise predictive density
log_mean_weibweib <- log(rowMeans(exp(ll_mat_weibweib)))
lppd_weibweib <- sum(log_mean_weibweib)
log_mean_expweib <- log(rowMeans(exp(ll_mat_expweib)))
lppd_expweib <- sum(log_mean_expweib)
log_mean_expexp <- log(rowMeans(exp(ll_mat_expexp)))
lppd_expexp <- sum(log_mean_expexp)

# Effective number of parameters
p_waic_weibweib <- sum(apply(ll_mat_weibweib, 1, var))
p_waic_expweib <- sum(apply(ll_mat_expweib, 1, var))
p_waic_expexp <- sum(apply(ll_mat_expexp, 1, var))

# WAIC
waic_weibweib <- -2 * (lppd_weibweib - p_waic_weibweib)
waic_expweib <- -2 * (lppd_expweib - p_waic_expweib)
waic_expexp <- -2 * (lppd_expexp - p_waic_expexp)

WAIC<-cbind(waic_weibweib,waic_expweib,waic_expexp)

saveRDS(WAIC, file = "in_text_results/WAIC.rds")
