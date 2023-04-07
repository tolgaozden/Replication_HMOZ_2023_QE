%% this script returns posterior moments of all models and writes 
% them into a spreadsheet. 

%% input files: 

%MCMC simulation results stored in subfolder "mcmc_baseline" for models 
% reported in Table 1.

%MCMC simulation results stored in subfolder "mcmc_expectations" for models 
% reported in Table 3. 

clear;
clc;
close all;

addpath(genpath(cd));

%compute covariance matrix and MHM estimator for all models 
compute_MHM;

%compute posterior moments of parameter and write into an excel file 
mcmc_summary_stat_all_models;

