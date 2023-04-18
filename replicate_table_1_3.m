clear;
clc;
addpath(genpath(cd));

warning('off','all');

%inputs: 

%estimation results at the posterior mode. 
%the results reported in the paper are based on .mat files stored in: 

%posterior_mode/estimations_baseline/results/[model_name]_estimation_results.mat (for table 1)
%posterior_mode/estimations_with_exp/results/[model_name]_estimation_results.mat (for table 3)

%MCMC results
%the results reported in the paper are based on .mat files stored in: 

%MCMC_results/mcmc_baseline/mcmc_[model_name].mat (for table 1)
%MCMC_results/mcmc_expectations/mcmc_[model_name].mat (for table 3)




%% the results will be written into the following spreadsheets: 
% mcmc_results.xlsx (table 1)
% mcmc_results_exp.xlsx (table 3)

% input files: 

% MCMC_results/mcmc_baseline/mcmc_[model_name].mat
% MCMC_results/mcmc_expectations/mcmc_[model_name].mat

run_me_mcmc_results;
% mcmc_summary_stat_all_models;



