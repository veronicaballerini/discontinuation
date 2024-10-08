# Evaluating causal effects on time-to-event outcomes in an RCT in Oncology with treatment discontinuation
Authors: Veronica Ballerini, Björn Bornkamp, Alessandra Mattei, Fabrizia Mealli, Craig Wang, Yufen Zhang.

arXiv: https://arxiv.org/abs/2310.06653v1

This code reproduces the results for Scenario 1 and Scenario 2 in Section 5.4. 

Scenario 1 depicts a situation in which the principal causal effects are positive in all the latent strata, reflecting the efficacy of the treatment. 
Scenario 2 represents a more challenging case of a positive overall effect, i.e., $ITT > 0$, but the treatment assignment has no effect for D patients; $ACE_{\text{D}} = 0$. 
In other words, we mimic a situation where the treatment does not show efficacy due to discontinuation.
The summary statistics of the data simulated under such scenarios are very close between them and similar to the summaries of the real data. 

Data simulation and model estimation are included in the folders.

In both scenarios, we make some assumptions on discontinuation for synthetic data simulation.
- The higher the value of $X_1$, the lower the probability of being an ND patient.
- Patients with $X_2$, $X_3$ equal to $1$ are more likely to be ND patients, i.e., patients who experience progression-free survival (PFS) without discontinuing.
- For D patients, higher-risk patients are more likely to discontinue sooner.

For convenience, we standardise the continuous covariate $X_1$, and we use the standardised version in the data-generating process and in the estimation. 
With a little abuse of notation, we continue denoting $\mathbf{X}$ the vector of covariates including the standardised $X_1$.
