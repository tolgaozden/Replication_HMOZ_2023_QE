# Replication Package for "Behavioral Learning Equilibria in New Keynesian Models." by C. Hommes, K. Mavromatis, T. Ozden, and M. Zhu. 

Working paper link (2022): 

https://www.bankofcanada.ca/2022/09/staff-working-paper-2022-42/

### Contact information: 

Latest update on 04.04.2023. Contact tozden@bankofcanada.ca for any questions and inquiries. 

### Software requirements:

The codes in this package have been tested on a Linux server with MATLAB 2021a using global optimization, parallel computing and symbolic math toolboxes. Parts of the package utilize the routines in Dynare software version 4.6.1 (For instructions on how to install Dynare, please visit https://www.dynare.org/).  

### Databases: 

There are two key databases used for estimations and simulations in the study, which can be found in main_databases subfolder. These variables are also provided in a spreadsheet in main_databases/database.csv. Below we provide details of all variables, mnemonics and data sources. 

1. "slobodyan_dataset.mat" contains publicly available quarterly macroeconomic time series for U.S. The database is taken from the study of Slobodyan and Wouters (2012). See the replication package of the article for all data sources and transformations: 

https://www.aeaweb.org/articles?id=10.1257/mac.4.2.65

In this study we use q/q CPI inflation as our measure of inflation. This is retrieved from FRED database (Consumer Price Index for All Urban Consumers (FRED mnemonic: CPIAUCSL):

https://fred.stlouisfed.org/series/CPIAUCSL

The variable mnemonics contained in the database are as follows: 
    
    [cpi_quarterly]: U.S. CPI inflation (quarter/quarter change). 

    [dc]: U.S. consumption growth rate(quarter/quarter change). 

    [dinve]: U.S. business investment growth rate (quarter/quarter change). 

    [dw]: U.S. wage growth rate (quarter/quarter change). 

    [dy]: U.S. output growth rate (quarter/quarter change). 

    [hours]: U.S. average hours worked. 

    [labobs]: U.S. average hours worked (normalized). 

    [pinfobs]: U.S. GDP deflator inflation (quarter/quarter change). 

    [robs]: U.S. FED funds rate (quarterly rate). 

2. "expectations_dataset.mat" contains the same variables as above, with some additional time series on inflation expectations. These are taken from the Survey of Professional Forecasters of Philadelphia Fed, see https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/survey-of-professional-forecasters for 
further details. The additional variable mnemonics related to inflation expectations are as follows: 

    [cpi_annualized]: Annualized rate of cpi_quarterly (i.e. 4 * cpi_quarterly). 
    
    [inf_exp_1_step]: 1-quarter ahead CPI inflation expectations from Survey of Professional Forecasters. Database: Median Forecast Data for Levels, CPI Table, Mnemonic "CPI3". 

    [inf_exp_2_step]: 2-quarter ahead CPI inflation expectations from Survey of Professional Forecasters. Database: Median Forecast Data for Levels, CPI Table, Mnemonic "CPI4". 

    [inf_exp_3_step]: 3-quarter ahead CPI inflation expectations from Survey of Professional Forecasters. Database: Median Forecast Data for Levels, CPI Table, Mnemonic "CPI5". 

    [inf_exp_4_step]: 4-quarter ahead CPI inflation expectations from Survey of Professional Forecasters. Database: Median Forecast Data for Levels, CPI Table, Mnemonic "CPI6". 



--------------------------
###### Data citation information 
--------------------------

[1] Slobodyan, S., & Wouters, R. (2012). Learning in a medium-scale DSGE model with expectations based on small
forecasting models. American Economic Journal: Macroeconomics, 4(2), 65-101.

Online Appendix publicly available at 

https://www.aeaweb.org/articles?id=10.1257/mac.4.2.65.

[2] Croushore, D. D. (1993). Introducing: the survey of professional forecasters. Business Review-Federal 
Reserve Bank of Philadelphia, 6, 3.

Survey results publicly available at 

https://www.philadelphiafed.org/surveys-and-data/real-time-data-research/survey-of-professional-forecasters



### Replication scripts: 

The package is accompanied by a set of wrapper scripts to directly reproduce the figures and tables in the paper. 

**replicate_figure_1.m**: uses routines in subfolder "3equation" to replicate figure 1. 

**replicate_figure_2.m**: uses routines in subfolder "3euquation" to replicate figure 2. 

**replicate_table_1_3.m**: uses routines in subfoldfer "MCMC_results" to replicate Tables 1 and 3. Note that this script does not re-run posterior mode estimation and MCMC simulations. It uses the provided databases in the following subfolders to compute the posterior moments and create the corresponding tables: 

    posterior_mode/estimations_baseline/results/[model_name]_estimation_results.mat (posterior mode of baseline models in table 1)

    posterior_mode/estimation_with_exp/results/[model_name]_estimation_results.mat (posterior mode of models with inflation expectations in table 3)

    posterior_mode/estimations_baseline/dynare_initial_beliefs/SW_Estimation_REE_results.mat (posterior mode with Dynare implementation of baseline RE model)

    posterior_mode/estimations_with_exp/dynare_initial_beliefs/SW_Estimation_REE_expData.mat (posterior mode with Dynare implementation of RE model with inflation exp.)

    MCMC_results/mcmc_baseline/mcmc_[model_name].mat (posterior distribution of baseline models in table 1)

    MCMC_results/mcmc_expectations/mcmc_[model_name].mat (posterior distribution of models with inflation expectations in table 3)

The routines used in the posterior mode estimation and MCMC simulations are explained in further detail below. 

**replicate_table_2.m**: uses routines in subfolder "posterior_mode" to replicate Table 2. Note that this script does not re-run the pseudo-out-of-sample forecasting exercise. It uses the provided databases in the following subfolder to compute the forecast errors, underlying RMSE's and the comparison table: 

posterior_mode/estimations_baseline/forecast_summary/forecast_output_[model_name].mat 

The routines used in generating forecasts are explained in further detail below. 

**replicate_figure_3_table_4.m**: uses routines in posterior_mode/estimations_with_exp and MCMC_samplers/MCMC_sampler_with_expectation to replicate figure 3 and the correlations in table 4. This script makes use of the same databases that are highlighted for tables 1 and 3. Given the posterior distribution databases of each model, the Kalman filter is re-run 1000 times with posterior parameter draws to obtain HPD bands of inflation expectations. The correlations are generated using the point estimated at the posterior mode. 

**replicate_figure_4.m**: uses routines in posterior_mode/estimations_with_exp to replicate figure 4. This script uses the posterior mode estimations of BLE and SAC models, same databases for tables 1 and 3, to generate the figure. 

**replicate_figure_5_table_5**: uses routines in optimal_policy subfolder to replicate table 5 and figure 5. This script uses the posterior mode estimations, same databases used for tables 1 and 3, to simulate each model over a grid of smoothing parameter and compute the optimal parameter values given the objective function. 


### Main sections:

The package consists of 6 main sections. 

(1) **3equation**: contains routines for computing some moments of the simple 3-equation New Keynesian model.

**varratio.m**: computes and plots the variance ratios reported in figure 2 of the paper. 

**blew.m**: computes the beta functions reported in figure 1 of the paper. 

**beta2star_rho.m**: computes and plots and persistence as a function of shock persistence reported in figure 2 of the paper. 

(2) **posterior_mode**: contains routines for posterior mode search of all models. This contains two subfolders: 

2.1) estimations_baseline: contains routines for mode search of the baseline model specifications. 
    
2.2) estimations_with_exp: contains routines for mode search of altenrative model specifications that use survey data on inflation expectations as an input. 

Each subfolder contains the following sections: 

-results: contains the estimation and Kalman filter output files at the posterior mode in .mat format for each model under consideration. The estimation result files also contain information about which options should be used in the estimation routines to get the corresponding output. 
 
-main_functions: contains the main routines to carry out the posterior mode search of the models under consideration. Further information about the contents of these routines are provided in "Background functions" section. 
 
-inputs: contains .mat files that are used for initializing the estimation routines. 
 
-helpers: contains some helper functions for optimization, reporting and specifying priors of estimated parameters. 
 
-forecast_summary: contains .mat files with summary statistics of pseudo out-of-sample forecasting exercise. These databases also contain information about which options should be used in the estimation routines to get the corresponding output. 
 
-dynare_initial_beliefs: contains Dynare routines to estimate the baseline rational expectation model. The results of the rational expectation model are used for initializing the beliefs of some of the adaptive learning models as well. 
 
-Matfiles: contains some output files of the csminwel optimization routine. 
 
 ##### Background files: 
 
dynare_initial_beliefs/slobodyan_dataset.mat: main database used in estimations with historical data on U.S. This database is taken from the replication package of Slobodyan and Wouters (2012, AEJ:Macro). See the section on data information for further details. 
 
dynare_initial_beliefs/SW_Estimation_REE.mod: Dynare file for estimating the baseline rational expectations model. The results of this file are also used for initializing beliefs for some of the adaptive learning models. 
 
dynare_initial_beliefs/beliefs_initialization_[model_name].m: the scripts use the rational expectation model as an input and run an OLS regression from the RE-based simulation results to generate initial beliefs consistent with the underlying forecasting rule for each adaptive learning model. 
 
**main.m**: this is the main wrapper for carrying out all posterior mode estimations. BLE and all adaptive learning models, as well as all underlying options related to belief initialization, timing of expectations, the use of projection facilities, which optimizer to use etc. are specified in this file. When estimation procedure is finished, the results are saved into results/estimation_results.mat under the default options. 

For the results reported in the paper, the underlying posterior modes are saved under the names results/[model_name]estimation_results.mat for each model. These databases also contain the set of options that must be specified in main.m in order to obtain the corresponding results. 
 
**forecast_main.m**: this is the main wrapper to carry out pseudo-out-of-sample forecasting exercises for all models. For each selected model, an estimation is carried out at each quarter to obtain the corresponding posterior mode. The forecasts are then generated under the obtained set of parameter values at each quarter. The estimation results at each quarter are saved into a database in subfolder auxiliary_files. The underlying posterior mode databases used in the paper are not provided due to space limitations. 

The summary statistics and key inputs needed for pseudo-out-of-sample forecasting of each model are stored in subfolder forecast_summary under the name forecast_output_[model_name].mat. 

**forecast_tables.m**: produces forecast errors, RMSE's and the corresponding tables given the forecast output files generated by forecast_main.m. 

**update_matrices.m**: Function that uses the system matrices and beliefs to calculate the recursive matrices of the model.  

**update_beliefs.m**: Function that applies the belief updating step for a specified information set. 

**sysmat_SW_model.m**: Script that specifies that symbolic system matrices associated with Smets Wouters model. 

**symbolic_to_matrix_MSV.m**: Function that returns the symbolic system matrices associated with Smets Wouters model for a given set of model equations

**SW_sysmat_MSV_filter.m**: Function s`that returns the system matrices associated with Smets Wouters model for a given set of parameter values.

**SW_prior.m**: Function that specifies prior distributions of all estimated parameters for Smets Wouters model. 

**SW_fixedPoint.m**: Function that uses the model object and finds the Behavioral Learning Equilibrium fixed-point associated with the model. 

**REE_solve_uhlig.m**: Function that uses the model object and finds the Rational Expectations Equilibrium fixed-point associated with the model.

**ree_auxiliary_retrieve_moments.m**: an auxiliary function to retrieve moments of Rational Expectations model for each forecast vintage. 

**projection_facility.m**: Computes the projection facility for a model. Our version follows Slobodyan-Wouters approach here (which in turn is 
based on previous learning literature), which amounts to ignoring the data from latest period if it leads to explosive dynamics. If the 
largest eigenvalue exceeds 1, leave all variables at their previous value (i.e. ignore the last data) and flag the projection_facility=1. 
Otherwise keep the updated values, and leave the flag unchanged (the flag may have been changed to 1 in previous steps if something else went wrong.)

**point_forecast.m**: Returns up to xxx-step ahead forecasts given the model. 

**plot_alphas_betas.m**: For estimated BLE and SAC-learning models, this script plots the estimated alphas and betas over history.

**param_set.m**: Auxiliary file to do some re-ordering on the parameters and also set the fixed parameters. 

**MSV_regression.m**: Runs a multicollinearity test for a specified matrix of data. 

**msv_learning.m**: Generic function that applies the constant gain learning step for recursive least squares models for specified PLM. 

**likelihood.m**: Calculates the likelihood of the model given the Kalman filter output and the prior distributions. 

**kalmanSW.m**: Function that takes in data and model object, finds the recursive system associated with the equilibrium/learning model,
runs it through the Kalman filter and returns the likelihood value. 

**KF_output_MCMC.m**: Runs the Kalman filter for specified model for posterior distribution of estimated parameters. 

**KF_output.m**: Same as KalmanSW.m in script form instead of function. 

**cgl_learning_recursive.m**: This function is used for the updating step of the SAC-learning algorithm with AR(1) rule. 

**initialize_BLE_fixed_model.m**: Function to specify initial beliefs associated with the estimated BLE model. 

**AR1_BLE_initial_beliefs.m**: same as initialize_BLE_fixed_model.m

 "helpers" subfolder contains a set of generic functions that are used for calculating moments/hessian matrices or prior distributions. 
 
 
(2) MCMC_samplers: contains MCMC simulation codes for the models. Each model has a MetropolisHastings_[model_name].m script that starts the MCMC simulation. These scripts use the estimated posterior mode from (1) as candidate density to initialize the MCMC sampling. The covariance matrix of the proposal density is typically initialized with a diagonal matrix. This is followed by a short MCMC chain (5000 draws) to obtain a more reasonable covariance matrix, which is used for initializing a long MCMC chain that is used for computing the posterior moments reported in the paper.

The simulations for rational expectation models are carried out using Dynare routines, which are contained in subfolder dynare_initial_beliefs. 


(3) MCMC_results: contains the MCMC simulation results of all models under consideration. These results are provided in .mat databases since the sampling step can take several days to weeks depending on the number of simulations. 

A wrapper script run_me.m collects all MCMC simulations and computes posterior moments. This script also extracts 1000 posterior draws from each model. These draws are used as input to generate Kalman filter results with the posterior distribution. These distributions are used as an input to compute model-implied inflation expectations and the corresponding uncertainty band around it in Section 3.4 of the paper. 

(4) optimal_policy: 
 
### Background functions: 


[1.6] generate_figure_1.m: Executes relevant functions and scripts to replicate Figure 1 in the study. 

Format: MATLAB file

[1.7] generate_figure_2.m: Executes relevant functions and scripts to replicate Figure 1 in the study. 

Format: MATLAB file

[1.8] generate_figure_3.m: Executes relevant functions and scripts to replicate Figure 1 in the study. 

Format: MATLAB file

[1.9] generate_figure_4.m: Executes relevant functions and scripts to replicate Figure 1 in the study. 

Format: MATLAB file

[1.10] generate_figure_5.m: Executes relevant functions and scripts to replicate Figure 1 in the study. 

Format: MATLAB file

[1.11] table1_metadata.xlsx: Contains posterior moments that are used for generating Table 1 in the study. 

Format: XLSX file

[1.12] table3_metadata.xlsx: Contains posterior moments that are used for generating Table 3 in the study. 

Format: XLSX file 

[1.13] table_5_metadata.csv: Contains simulation results that are used for generating Table 5 in the study. 

[2] "codes" folder includes all required MATLAB programs to replicate the results in the study. 

[2.1] beliefs_initialization_AR1.m: This file takes as input the estimated Rational Expectations model from SW_Estimation_REE.mod and uses the 
recursive system to estimate an AR(1) OLS model to initialize beliefs for SAC-learning.

Format: MATLAB file

[2.2] beliefs_initialization_AR2.m: This file takes as input the estimated Rational Expectations model from SW_Estimation_REE.mod and uses the 
recursive system to estimate an AR(2) OLS model to initialize beliefs for AR(2) constant gain learning. 

Format: MATLAB file

[2.3] beliefs_initialization_MSV.m: This file takes as input the estimated Rational Expectations model from SW_Estimation_REE.mod and uses the 
recursive system to estimate an MSV-consistent model to initialize beliefs for MSV-learning model. 

Format: MATLAB file

[2.4] beliefs_initialization_MSV_t_1.m: This file is a replicate of beliefs_initialization_MSV.m but estimates a model consistent with t-1
timing of expectations.

Format: MATLAB file


[2.5] beliefs_initialization_ree.m: This file takes as input the estimated Rational Expectations model from SW_Estimation_REE.mod and uses the 
recursive system to estimate a REE-consistent model to initialize beliefs for a REE model without cross equation restrictions.  

Format: MATLAB file

[2.6] beliefs_initization_VAR.m: This file takes as input the estimated Rational Expectations model from SW_Estimation_REE.mod and uses the 
recursive system to estimate a VAR(1) OLS model to initialize beliefs for VAR(1) constant gain learning.

Format: MATLAB file


[2.7] beta2star_rho.m: This file calculates and returns the equilibrium beta value in the small New Keynesian model as a function of the exogenous
shock persistence. This file is used for generating figure 2(a) in the study. 

Format: MATLAB file

[2.8] blew.m: This file calculates the theoretical values of beta_1 as a function of beta_2, and beta_2 as a function of beta_1. This file is used 
for generating figure 1 in the study. 

Format: MATLAB file

[2.9] cgl_learning_recursive.m: This function is used for the updating step of the SAC-learning algorithm with AR(1) rule. 

Format: MATLAB file

[2.10] expectations_statistics_mcmc.m: This file is used for calculating and plotting model-implied inflation expectations for a given posterior 
distribution of estimated parameter values. 

Format: MATLAB file

[2.11] expectatios_statistics_mode.m: This file is used for calculating and plotting model-implied inflation expectations for a given point estimate
of the posterior model. 

Format: MATLAB file

[2.12] forecast_comparisons.m: Given the forecast_output_xxx.mat files, it reproduces the out-of-sample forecast comparisons tables and 
plots used in the table. 

Format: MATLAB file

[2.13] forecast_evaluation.m: Returns the output variable that includes the out-of-sample statistics RMSE,MSE,MAE and FE,uncentered log determinant
under consideration. 

Format: MATLAB file

[2.14] forecast_main.m: Carries out all necessary steps, i.e. estimation of the model at each time step//forecasting//saving the results given the 
model and the forecast time period. The auxiliary files for the initial beliefs at each step need to be pre-computed and saved. 
after the estimation at each step, calls functions: point_forecast.m to obtain the period-specific forecasts, then saves it to 
forecast_output_xxx.mat file to be used in the forecast_evaluation.m file.

Format: MATLAB file

[2.15] forecast_retrieve_results.m: Re-evaluates forecasts and forecast errors for given model objects, without recursively re-estimating the model 
at every period. 

Format: MATLAB file

[2.16] forecast_tables.m: Generates the point forecast comparison tables and saves it to a .csv/.xlsx file.

Format: MATLAB file

[2.17] generate_figures.m: Generic function to plot the figure as a .pdf to a specified location.

Format: MATLAB file

[2.18] initial_belief_special_cases.m: Generates initial beliefs for several cases not covered in initialize_beliefs.m

Format: MATLAB file

[2.19] initialize_beliefs.m: Function to specify initial beliefs for learning models. 

Format: MATLAB file

[2.20] initialize_BLE_fixed_model.m: Function to specify initial beliefs associated with the estimated BLE model. 

Format: MATLAB file

[2.21] kalmanSW.m: Function that takes in data and model object, finds the recursive system associated with the equilibrium/learning model,
runs it through the Kalman filter and returns the likelihood value. 

Format: MATLAB file

[2.22] KF_output.m: Same as KalmanSW.m in script form instead of function. 

Format: MATLAB file

[2.23] KF_output_MCMC.m: Runs the Kalman filter for specified model for posterior distribution of estimated parameters. 

Format: MATLAB file

[2.24] likelihood.m: Calculates the likelihood of the model given the Kalman filter output and the prior distributions. 

Format: MATLAB file

[2.25] main.m: Wrapper script to launch estimation of a specified model. 

Format: MATLAB file

MCMC_figures_tables.m: Plots and returns the MCMC tables of an estimated model. 

[2.26] MetropolisHastings.m: Wrapper scripts to launch the Metropolis Hastings algorithm that runs an MCMC simulation and returns the posterior 
distribution of the model. 

Format: MATLAB file

[2.27] msv_learning.m: Generic function that applies the constant gain learning step for recursive least squares models for specified PLM. 

Format: MATLAB file

[2.28] MSV_regression.m: Runs a multicollinearity test for a specified matrix of data. 

Format: MATLAB file

[2.29] param_set.m: Auxiliary file to do some re-ordering on the parameters and also set the fixed parameters. 

Format: MATLAB file

[2.30] plot_alphas_betas.m: For estimated BLE and SAC-learning models, this script plots the estimated alphas and betas over history.

Format: MATLAB file

[2.31] point_forecast.m: Returns up to xxx-step ahead forecasts given the model. 

Format: MATLAB file

[2.32] projection_facility.m: Computes the projection facility for a model. Our version follows Slobodyan-Wouters approach here (which in turn is 
based on previous learning literature), which amounts to ignoring the data from latest period if it leads to explosive dynamics. If the 
largest eigenvalue exceeds 1, leave all variables at their previous value (i.e. ignore the last data) and flag the projection_facility=1. 
Otherwise keep the updated values, and leave the flag unchanged (the flag may have been changed to 1 in previous steps if something else went wrong.)

Format: MATLAB file

[2.33] REE_solve_uhlig.m: Function that uses the model object and finds the Rational Expectations Equilibrium fixed-point associated with the model.

Format: MATLAB file

[2.34] SW_fixedPoint.m: Function that uses the model object and finds the Behavioral Learning Equilibrium fixed-point associated with the model. 

Format: MATLAB file

[2.35] SW_prior.m: Function that specifies prior distributions of all estimated parameters for Smets Wouters model. 

Format: MATLAB file

[2.36] SW_sysmat_MSV_filter.m: Function that returns the system matrices associated with Smets Wouters model for a given set of parameter values.

Format: MATLAB file

[2.37] symbolic_to_matrix_MSV.m: Function that returns the symbolic system matrices associated with Smets Wouters model for a given set of model equations.

Format: MATLAB file

[2.38] sysmat_SW_model.m: Script that specifies that symbolic system matrices associated with Smets Wouters model. 

Format: MATLAB file

[2.39] update_beliefs.m: Function that applies the belief updating step for a specified information set. 

Format: MATLAB file

[2.40] update_matrices.m: Function that uses the system matrices and beliefs to calculate the recursive matrices of the model. 

Format: MATLAB file

[2.41] varratio.m: Function that calculates that variance ratios in the small New Keynesian model. This file is used for generating Figure 2(a) in the study. 

Format: MATLAB file

[2.42] SW_Estimation_REE.mod: Dynare file that specifies the baseline Rational Expectations version of the model. 

Format: MATLAB file


Format: MATLAB file




