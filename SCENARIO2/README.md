# SCENARIO 2

We simulate the data such that the principal causal effects for D patients are zero. 
However, the treatment has a positive effect on ND patients; the overall ITT effect is positive.

- Generate data using file 00_datageneration_SC2.R; the output files are: simdata2_4800.txt, realata2_4800.txt, thetatrue2_4800.txt. It prints KM figures.
- Initialize the analysis running file 01_MCMC_Analysis_SC2.R. The input files are contained in folder "SCENARIO 1": SCENARIO1/MCMC.R,  SCENARIO1/CompleteLogPost.R,  SCENARIO1/DataAugmentation.R,  SCENARIO1/Functions.R. It saves results in a .RData file reporting date and time
- Check results running file 02_results.R (you will need to change the name of the .RData file). It requires the additional input files  SCENARIO1/ce_apply.R and  SCENARIO1/Estimands.R.
- Visualize results running file 03_graphs.R immediately after 02_results.R. It prints ACE_d.jpeg, DCE_ND.jpeg and DCE_D.jpeg.
