clear;
clc;

%% models with inf expectations

num_for_simul=1000;

file_name='mcmc_sac_exp.mat';
input_draws=500000;
[mcmc_sac_exp]=MHM_estimator(file_name,input_draws,num_for_simul);

file_name='mcmc_ble_exp.mat';
input_draws=500000;
[mcmc_ble_exp]=MHM_estimator(file_name,input_draws,num_for_simul);

file_name='mcmc_msv_exp.mat';
input_draws=500000;
[mcmc_msv_exp]=MHM_estimator(file_name,input_draws,num_for_simul);

file_name='mcmc_var1_exp.mat';
input_draws=500000;
[mcmc_var1_exp]=MHM_estimator(file_name,input_draws,num_for_simul);

file_name='mcmc_ar2_exp.mat';
input_draws=500000;
[mcmc_ar2_exp]=MHM_estimator(file_name,input_draws,num_for_simul);


TT=table(mcmc_sac_exp.marginal(:,2),mcmc_ble_exp.marginal(:,2),mcmc_msv_exp.marginal(:,2),mcmc_var1_exp.marginal(:,2),mcmc_ar2_exp.marginal(:,2));
TT.Properties.VariableNames={'sac','ble','msv','var1','ar2'};
disp(TT);





%% models without inf expectations 

file_name='mcmc_sac.mat';
input_draws=500000;
[mcmc_sac]=MHM_estimator(file_name,input_draws,num_for_simul);

file_name='mcmc_ble.mat';
input_draws=500000;
[mcmc_ble]=MHM_estimator(file_name,input_draws,num_for_simul);

file_name='mcmc_msv.mat';
input_draws=500000;
[mcmc_msv]=MHM_estimator(file_name,input_draws,num_for_simul);

file_name='mcmc_var1.mat';
input_draws=500000;
[mcmc_var1]=MHM_estimator(file_name,input_draws,num_for_simul);

file_name='mcmc_ar2.mat';
input_draws=500000;
[mcmc_ar2]=MHM_estimator(file_name,input_draws,num_for_simul);


TT2=table(mcmc_sac.marginal(:,2),mcmc_ble.marginal(:,2),mcmc_msv.marginal(:,2),mcmc_var1.marginal(:,2),mcmc_ar2.marginal(:,2));
TT2.Properties.VariableNames={'sac','ble','msv','var1','ar2'};
disp(TT2);


save mcmc_summary_all_models.mat mcmc_sac mcmc_ble mcmc_msv mcmc_var1 mcmc_ar2 ... 
    mcmc_sac_exp mcmc_ble_exp mcmc_msv_exp mcmc_var1_exp mcmc_ar2_exp;