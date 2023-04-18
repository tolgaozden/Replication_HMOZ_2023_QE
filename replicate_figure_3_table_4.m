clear;
cd('posterior_mode/estimations_with_exp');
addpath(genpath(cd));

warning('off','all');

%for each model: run Kalman filter 1000 times with posterior distribution
%first, then generate HPD bands for inflation expectations.


%inputs: 

%estimation results at the posterior mode. The results reported in the
%paper are based on the files in: 

%posterior_mode/estimations_with_exp/results/[model_name]_estimation_results.mat

% a database containing 1000 posterior parameter draws for each model. 
% the results reported in the paper are based on the database in 

%posterior_mode/estimations_with_exp/inputs/mcmc_summary_all_models.mat 

% this database in an output of run_me_mcmc_results.m, which uses the full
% posterior distribution of each model to extract 1000 equally spaced mcmc
% draws to minimize the autocorrelation across the parameter draws. 


%run Kalman filter with MCMC draws. 
%store results into BLE_posteriorDist.mat
KF_output_MCMC_ble

%use BLE_posteriorDist.mat (for posterior dist.) and ble_kf_output.mat (for posterior mode)
%to generate inflation expectation statistics
expectations_statistics_mcmc_ble;
expectations_statistics_mode_ble;

%run Kalman filter with MCMC draws. 
%store results into SAC_posteriorDist.mat
KF_output_MCMC_sac;

%use SAC_posteriorDist.mat (for posterior dist.) and sac_kf_output.mat (for posterior mode)
%to generate inflation expectation statistics
expectations_statistics_mcmc_sac;
expectations_statistics_mode_sac;


%run Kalman filter with MCMC draws. 
%store results into MSV_posteriorDist.mat
KF_output_MCMC_msv;

%use MSV_posteriorDist.mat (for posterior dist.) and msv_kf_output.mat (for posterior mode)
%to generate inflation expectation statistics
expectations_statistics_mcmc_msv;
expectations_statistics_mode_msv;

%run Kalman filter with MCMC draws. 
%store results into VAR_posteriorDist.mat
KF_output_MCMC_var1;

%use VAR_posteriorDist.mat (for posterior dist.) and var_kf_output.mat (for posterior mode)
%to generate inflation expectation statistics
expectations_statistics_mcmc_var1;
expectations_statistics_mode_var1;


%run Kalman filter with MCMC draws. 
%store results into AR2_posteriorDist.mat
KF_output_MCMC_ar2;

%use AR2_posteriorDist.mat (for posterior dist.) and ar2_kf_output.mat (for posterior mode)
%to generate inflation expectation statistics
expectations_statistics_mcmc_ar2;
expectations_statistics_mode_ar2;

%generate RE-implied inflation expectations at posterior mode with dynare
expectations_statistics_ree;
cd('../../MCMC_samplers/MCMC_sampler_expectation/mcmc_for_ree_dynare/');

%run dynare with MCMC draws for posterior distribution and generate plots
%for RE 
KF_output_REE;
