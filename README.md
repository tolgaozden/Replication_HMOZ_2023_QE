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

(2) posterior_mode: contains routines for posterior mode search of all models. This contains two subfolders: 

2.1) estimations_baseline: contains routines for mode search of the baseline model specifications. 
    
2.2) estimations_with_exp: contains routines for mode search of altenrative model specifications that use survey data on inflation expectations as an input. 

Each subfolder contains the following sections: 

- results: contains the estimation and Kalman filter output files at the posterior mode in .mat format for each model under consideration. The estimation result files also contain information about which options should be used in the estimation routines to get the corresponding output. 
 
- main_functions: contains the main routines to carry out the posterior mode search of the models under consideration. Further information about the contents of these routines are provided in "Background functions" section. 
 
- inputs: contains .mat files that are used for initializing the estimation routines. 
 
- helpers: contains some helper functions for optimization, reporting and specifying priors of estimated parameters. 
 
- forecast_summary: contains .mat files with summary statistics of pseudo out-of-sample forecasting exercise. These databases also contain information about which options should be used in the estimation routines to get the corresponding output. 
 
- dynare_initial_beliefs: contains Dynare routines to estimate the baseline rational expectation model. The results of the rational expectation model are used for initializing the beliefs of some of the adaptive learning models as well. 
 
- Matfiles: contains some output files of the csminwel optimization routine. 
 
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
 
 The subfolder estimations_with_exp/main_functions additionally contains the following functions and scripts related to inflation expectaitons: 
 
 **sysmat_SW_model_[model_name].m** and **symbolic_to_matrix_[model_name].m**: functions to retrieve system matrices for each model under consideration. These have the same functionality as sysmat_SW_model.m and symbolit_to_matrix_MSV.m, but but there is a different function for each model since the extended state space including expectations are needed when inflation expectations are included as observables. 
 
 **KF_output_MCMC_[mode_name].m**: script to retrieve Kalman filter output using parameter draws simulated from the posterior distribution of each model. 
 
 **expectations_statistics_mode_[model_name].m**: script to compute inflation expectations/inflation expectation errors and related statistics such as correlations for each model, at the estimation posterior mode. 
 
 **expectations_statistics_MCMC_[model_name].m**: same as above, using parameter draws simulated from the posterior distribution to obtain confidence bands around model-implied inflation expectations. 
 
(2) MCMC_samplers: contains MCMC simulation codes for the models. Each model has a MetropolisHastings_[model_name].m script that starts the MCMC simulation. These scripts use the estimated posterior mode from (1) as candidate density to initialize the MCMC sampling. The covariance matrix of the proposal density is typically initialized with a diagonal matrix. This is followed by a short MCMC chain (5000 draws) to obtain a more reasonable covariance matrix, which is used for initializing a long MCMC chain that is used for computing the posterior moments reported in the paper.

The simulations for rational expectation models are carried out using Dynare routines, which are contained in subfolder dynare_initial_beliefs. 

 ##### Background files: 

The subfolders in MCMC_sampler_expectation/codes and MCMC_samplers_no_expectation/codes contain the same routines as in the posterior_mode folder.  

The subfolders in MCMC_sampler_expectation/dynare_initial_beliefs and MCMC_samplers_no_expectation/dynare_initial_beliefs contain Dynare routines to launch MCMC simulations for the underlying Rational Expectations models. 

**MetropolisHastings_[model_name].m**: launches the MCMC simulation for each model under consideration. The estimated posterior modes are used as the candidate density to initialize the simulations. The covariance matrix of the proposal density is obtained by running a short MCMC sample with 5000 draws, using a diagonal Hessian matrix. The results of this short MCMC sample are used to compute a covariance matrix, which is used in the main MCMC chain. 

**MCMC_sampler_expectation/mcmc_for_ree_dynare/KF_output_REE.m**: retrieves 1000 draws from the posterior distribution of the RE model, and runs the Kalman filter 1000 times with the posterior draws to obtain model-implied inflation expectations and the corresponding HPD bands around it. 

(3) MCMC_results: contains the MCMC simulation results of all models under consideration. These results are provided in .mat databases since the sampling step can take several days to weeks depending on the number of simulations. 

A wrapper script run_me_mcmc_results.m collects all MCMC simulations and computes posterior moments. This script also extracts 1000 posterior draws from each model. These draws are used as input to generate Kalman filter results with the posterior distribution. These distributions are used as an input to compute model-implied inflation expectations and the corresponding uncertainty band around it in Section 3.4 of the paper. 

#### Background functions: 


**MHM_estimator.m**: computes the Modified Harmonic Mean estimator given the full posterior distribution of a model. Also retrieves a specified number of parameters draws from the posterior distribution. The thinning factor for these parameter draws are set to minimize the autocorrelation in the resulting draws. 

**compute_MHM.m**: collects the MCMC distributions of all models and returns MHM estimator along with 1000 posterior draws to be used in Kalman filter re-runs. 

**mcmc_summary_stat.m**: computes posterior moments of a model, given the MCMC distribution. 

**mcmc_summary_stat_all_models.m**: collects summary statistics and moments of MCMC draws from all models and writes the results into a spreadsheet. 

(4) optimal_policy: 
 
#### Background functions: 




