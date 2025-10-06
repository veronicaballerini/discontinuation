This folder contains scripts to reproduce the results of the simulation study in Section 5. 
You can choose to reproduce either the full simulation study or an "intermediate" one, by setting the option "full" (line 38 of the main script 01_MCMC_Analysis.R) equal to TRUE or FALSE, respectively. 
If full == FALSE, you will use 2 datasets (instead of 150 of the full study) and implement the model for both scenarios without covariates.
Also, change option at line 25 of 03\_results.R.

The content of the folder follows.

- "01\_MCMC\_Analysis.R" is the master script that runs the simulation study. 
  Inputs: 
  - "Functions.R" -> contains all useful functions
  - "estimands\_apply.R" -> function to compute the posterior of the effects; models with covariates
  - "estimands\_apply\_nocov.R" -> function to compute the posterior of the effects; models w/o covariates
  - "Estimands.R" -> function to compute the "true" values of the effects
  - "applyfun\_1.R" -> function that i) generates data under SCENARIO I with covariates, ii) estimate the joint    
     posterior distribution using the model correctly specified, iii) compute the coverage, width, and bias of the 
     model
  - "applyfun\_1nocov.R" -> function that i) generates data under SCENARIO I with covariates, ii) estimate the joint
     posterior distribution using the model correctly specified with respect to the discontinuation but w/o covariates, 
     iii) compute the coverage, width, and bias of the model
  - "applyfun\_2.R" -> function that i) generates data under SCENARIO II with covariates, ii) estimate the joint 
     posterior distribution using the misspecified model with covariates, iii) compute the coverage, width, and bias of 
     the model
  - "applyfun\_2nocov.R" -> function that i) generates data under SCENARIO II with covariates, ii) estimate the joint 
     posterior distribution using the misspecified model w/o covariates, iii) compute the coverage, width, and bias of 
     the model
  - "MCMC.R"-> contains the MCMC described by Algorithm 1 in Appendix B. 
       Inputs: "CompleteLogPost.R" and "DataAugmentation.R"
  - "MCMC\_nocov.R" -> contains the MCMC described by Algorithm 1 in Appendix B, but w/o covariates. 
       Inputs: "CompleteLogPost\_nocov.R" and "DataAugmentation\_nocov.R"
  - "MCMC\_initialization.R" -> sets the values for the prior, the proposal and the starting values of the parameters
  - "seeds1.txt", "seeds1\_nocov.txt", "seeds2.txt", "seeds2\_nocov.txt": files that fix the seeds for reproducibility. 
     If full == FALSE, only the first 2 seeds for each scenario are used (only without-covariates-specification). 
  Outputs:
  - Files .txt named "realdata\*SCENARIO\*\_\*SEED\*" and "simdata*SCENARIO*\_\*SEED\*" (in "sim" folder)
  - Files .RData named "thetatrue\*SCENARIO\*\_\*SEED\*" (in "sim" folder)
  - Files .txt named "jump*SCENARIO\*\_\*SEED*" that check the acceptance rates are optimal and the chain jumps
  - Final results file: "final\_simdata.RData" if full == TRUE; "final\_simdata\_intermediate.RData" if full == FALSE.    
    These files contain information on the execution time of the simulation study, and results for Scenario I and II, 
    with and without covariates if full == TRUE (objects sim\_SCENARIO1, sim\_SCENARIO1\_nocov, 
     sim\_SCENARIO2, sim\_SCENARIO2\_nocov) and only without covariates if full == FALSE (objects sim\_SCENARIO1\_nocov, 
     and sim\_SCENARIO2\_nocov). Each of these objects is a list of number of elements equal to the used 
     datasets, and each element of such a list is, in turn, a list of 25 element that are used to compute coverage, 
     width, and bias. 
- "02\_KM.R" reproduces Figures 3-6 in Section 5 (to run after the simulation study ends)
  Inputs: 
  - realdata + observed data .txt files printed during the (full) simulation study. We used seed 4798 as an example for 
     Scenario I, and seed 4800 for Scenario II. 
  Outputs (in "plots" folder):
   - Figure3a, Figure3b, Figure4, Figure5a, Figure5b, Figure6 (.jpeg files).
- "03\_results.R" makes the graphs in Figures 7-12, and tables in Appendix C. Results for full/intermediate are
   different; set the option at line 25.
   Inputs: 
      - "tables.R" -> compute coverage, width, and bias. 
        Inputs: 
        - "final\_simdata.RData" if full == TRUE, or "final\_simdata\_intermediate.RData" otherwise.
        Outputs (in "tables" folder, subfolder "FULL" if full == TRUE, "INTERMEDIATE" otherwise): 
        - TableX.rds", X = 4, …, 15
   Outputs (in "plots" folder -> subfolder "FULL" if full == TRUE, "INTERMEDIATE" otherwise):
   - Figure7, Figure8, Figure9, FIgure10, Figure11, Figure12 (.jpeg files).
