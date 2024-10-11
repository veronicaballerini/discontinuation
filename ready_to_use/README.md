# Your own analysis

This folder contains scripts that will allow you to do your own analysis. 

- Read the "README.txt" file to learn more about the required data structure.
- "example_data.txt" is an example of how your dataset should look like.
- "01_MCMC_Analysis.R" is the first script to run. You can change the number of iterations, burnin, and thinning there.
- "02_results.R" must be run to obtain a table of results and as a necessary step before running scripts "03_graphs.R" and "04_covariates.R".
- "03_graphs.R" is to obtain $ACE_D(d)$, $DCE_{ND}$, and $DCE_D(d)$ graphs (see Figures 6-8 in the paper). It requires "02_results.R".
- "04_covariates.R" is to obtain covariates' distribution by stratum (see Figures 12-14 in the paper). It requires "02_results.R".
