clear;
cd('posterior_mode/estimations_with_exp');
addpath(genpath(cd));

%for each model: run Kalman filter 1000 times with posterior distribution
%first, then generate HPD bands for inflation expectations.

KF_output_MCMC_ble;
expectations_statistics_mcmc_ble;
expectations_statistics_mode_ble;

KF_output_MCMC_sac;
expectations_statistics_mcmc_sac;
expectations_statistics_mode_sac;

KF_output_MCMC_msv;
expectations_statistics_mcmc_msv;
expectations_statistics_mode_msv;


KF_output_MCMC_var1;
expectations_statistics_mcmc_var1;
expectations_statistics_mode_var1;

KF_output_MCMC_ar2;
expectations_statistics_mcmc_ar2;
expectations_statistics_mode_ar2;


expectations_statistics_ree;


cd('../..');