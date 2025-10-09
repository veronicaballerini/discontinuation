# Case study, Sections 2 and 6
This folder contains scripts that reproduce the results of the case study in Sections 2 and 6.
Run the following scripts in order (from 01 to 04, detailed below) to reproduce all results in Sections 2 and 6. 
If you run the scripts from 01 to 03, you will reproduce all results for the chosen model, namely the Exponential-Exponential, except for the WAIC value. To obtain the WAIC values, you need to run the "04_WAIC.R" file, which implements the two competing models. It takes > 2 hours to run.

- "01_MCMC_Analysis.R" is the master script that runs the MCMC with Exponential-Exponential specification. It takes about 1 hour to run.
Inputs: 
  - "synthetic_data.csv" -> dataset containing synthetic data based on a real case study, simulated using the "synthpop" package. It contains the following variables:
    	- Z: treatment assignment indicator
    	- RD.obs: indicator taking value 1 if discontinuation is observed during the follow-up
  		- D: time-to-discontinuation (if RD.obs == 1) or time-to-censoring (if RD.obs == 0)
    	- RY: indicator taking value 1 if the primary event is observed during the follow-up
    	- Y: PFS/time to the primary event (if RY == 1) or time-to-censoring (if RY == 0)
    	- ID.obs: it takes value 1 if the patient is observed to be a ND patient, namely RD.obs = 0 & RY = 1
    	- x1, x2, x3: covariates
  - "Functions.R" -> contains all useful functions
  - "descriptive.R" -> compute the descriptive statistics and produce Tables 1 and 2, and Figure 2 of Section 2.
  - "MCMC_expexp.R"-> contains the MCMC described by Algorithm 1 in Appendix B, but with exponential specifications. 
	Inputs: 
          - "CompleteLogPost_expexp.R" 
          - "DataAugmentation_expexp.R"
  "MCMC_initialization_expexp.R" -> sets the values for the prior, the proposal, and the starting values of the parameters
Outputs:
   - Final results file: "final_casestudy.RData" (in the working directory)
   - Figure 2 (.jpeg file, in the "figures" folder) 
   - Table1, Table2_left, and Table2_right (.rds files, in the "tables" folder) 

- "02_results.R" reproduces Figures 13-15 in Section 6. It takes a few minutes to run.
Inputs: 
   - "results_processing.R": processes the results.
	Inputs:
	  - "final_casestudy.RData"
          - "ce_apply_expexp.R" -> function to compute the posterior of the effects, the WAIC, and the Kaplan-Meier discrepancies
             for the Exponential-Exponential specifications
          - "RepData.R" -> function that draws replicated data from the joint posterior (for KM discrepancies computation)
          - "apply_cov.R" -> compute the posterior of the covariates' average value in each group - ND, late D, and early D, namely
             it takes the covariates' average value within strata of patients for each posterior draw of the latent memberships.
        Outputs:
          - "results.RData" -> contains the posterior and the distribution of the covariates.
          - "p_ND.rds" -> posterior mean and 95% HDI of the proportion of ND patients. In-text results in Section 6 (.rds file in "in_text_results" folder)
          - "PPP.rds" -> mean Bayesian PPP for KMdm. In-text results in Section 6.1 (.rds file in "in_text_results" folder)
          - "t_expexp.rds" -> posterior mean and 95% HDI of model's parameters. Only the delta parameter as in-text result, Section 6 (.rds file in "in_text_results" folder)
          - traceplots in folder "figures/trace/ExpExp" - not in the manuscript
Outputs:
   - Figures 13-15 (.jpeg files, in folder "figures").   

- "03_covariates -> to produce Figures 16-18. It takes a few seconds to run.
Input: 
   - "results.RData"
Outputs: 
   - Figures 16-18 (.jpeg files, in folder "figures").

- "04_WAIC.R" -> to obtain the WAIC values for model selection. It implements the Weibull-Weibull and Exponential-Weibull models. It takes > 2 hours to run.
Inputs: 
  - "synthetic_data.csv" -> dataset containing synthetic data based on the real case study, simulated using the "synthpop" package.
  - "Functions.R" -> contains all useful functions
  - "MCMC_weibeib.R"-> contains the MCMC described by Algorithm 1 in Appendix B 
	Inputs: 
          - "CompleteLogPost_weibweib.R" 
          - "DataAugmentation_weibweib.R"
  - "MCMC_expweib.R"-> contains the MCMC described by Algorithm 1 in Appendix B, but with exponential-weibull specification
	Inputs: 
          - "CompleteLogPost_expweib.R" 
          - "DataAugmentation_expweib.R"
  - "MCMC_initialization_weibweib.R" -> sets the values for the prior, the proposal, and the starting values of the parameters for the weib-weib model
  - "MCMC_initialization_expweib.R" -> sets the values for the prior, the proposal, and the starting values of the parameters for the exp-weib model
  - "ce_apply_weibweib.R" -> function to compute the posterior of the effects and the WAIC for the Weibull-Weibull specifications
  - "ce_apply_expweib.R" -> function to compute the posterior of the effects and the WAIC for the exponential-Weibull specifications
Outputs:
  - Results file: "other_models.RData" (in the working directory)
  - WAIC.rds ("in_text_results" folder)
  - "t_weibweib.rds", "t_expweib.rds", "t_expexp.rds" -> posterior mean and 95% HDI of model's parameters for each specification (tables/extra folder) - not in the manuscript
  - traceplots in folders "figures/trace/WeibWeib" and "figures/trace/ExpWeib" - not in the manuscript
