# SCENARIO 1

We generate data such that there is a positive treatment effect for all latent strata and, thus, a positive effect of the treatment for the ND patients.

- Generate data using file 00_datageneration.R; the output files are: simdata1_4816.txt, realata1_4816.txt, thetatrue_4816.txt. It prints KM figures.
- Initialize the analysis running file 01_MCMC_Analysis.R. The input files are: MCMC.R, CompleteLogPost.R, DataAugmentation.R, Functions.R. It saves results in a .RData file reporting date and time
- Check results running file 02_results.R (you will need to change the name of the .RData file). It requires the additional input files ce_apply.R and Estimands.R.
- Visualize results running file 03_graphs.R immediately after 02_results.R. It prints ACE_d.jpeg, DCE_ND.jpeg and DCE_D.jpeg.
