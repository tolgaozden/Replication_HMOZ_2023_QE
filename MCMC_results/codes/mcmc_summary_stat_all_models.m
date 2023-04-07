clear;
clc;

%% write no-expectation models to excel

load mcmc_summary_all_models.mat;

file_name='mcmc_sac.mat';
input_draws=500000;
[mcmc_summary_sac]=mcmc_summary_stat(file_name,input_draws);

file_name='mcmc_ble.mat';
input_draws=500000;
[mcmc_summary_ble]=mcmc_summary_stat(file_name,input_draws);

file_name='mcmc_msv.mat';
input_draws=500000;
[mcmc_summary_msv]=mcmc_summary_stat(file_name,input_draws);

file_name='mcmc_var1.mat';
input_draws=500000;
[mcmc_summary_var1]=mcmc_summary_stat(file_name,input_draws);

file_name='mcmc_ar2.mat';
input_draws=500000;
[mcmc_summary_ar2]=mcmc_summary_stat(file_name,input_draws);


writetable(mcmc_summary_sac,'mcmc_results.xlsx','sheet','sac');
writematrix(mcmc_sac.marginal(end,2),'mcmc_results.xlsx','sheet','sac','range','A40');

writetable(mcmc_summary_ble,'mcmc_results.xlsx','sheet','ble');
writematrix(mcmc_ble.marginal(end,2),'mcmc_results.xlsx','sheet','ble','range','A40');

writetable(mcmc_summary_msv,'mcmc_results.xlsx','sheet','msv');
writematrix(mcmc_msv.marginal(end,2),'mcmc_results.xlsx','sheet','msv','range','A40');

writetable(mcmc_summary_var1,'mcmc_results.xlsx','sheet','var1');
writematrix(mcmc_var1.marginal(end,2),'mcmc_results.xlsx','sheet','var1','range','A40');

writetable(mcmc_summary_ar2,'mcmc_results.xlsx','sheet','ar2');
writematrix(mcmc_ar2.marginal(end,2),'mcmc_results.xlsx','sheet','ar2','range','A40');


%% write expectation models to excel 

file_name='mcmc_sac_exp.mat';
input_draws=500000;
[mcmc_summary_exp_sac]=mcmc_summary_stat(file_name,input_draws);

file_name='mcmc_ble_exp.mat';
input_draws=500000;
[mcmc_summary_exp_ble]=mcmc_summary_stat(file_name,input_draws);

file_name='mcmc_msv_exp.mat';
input_draws=500000;
[mcmc_summary_exp_msv]=mcmc_summary_stat(file_name,input_draws);

file_name='mcmc_var1_exp.mat';
input_draws=500000;
[mcmc_summary_exp_var1]=mcmc_summary_stat(file_name,input_draws);

file_name='mcmc_ar2_exp.mat';
input_draws=500000;
[mcmc_summary_exp_ar2]=mcmc_summary_stat(file_name,input_draws);


writetable(mcmc_summary_exp_sac,'mcmc_results_exp.xlsx','sheet','sac');
writematrix(mcmc_sac_exp.marginal(end,2),'mcmc_results_exp.xlsx','sheet','sac','range','A40');

writetable(mcmc_summary_exp_ble,'mcmc_results_exp.xlsx','sheet','ble');
writematrix(mcmc_ble_exp.marginal(end,2),'mcmc_results_exp.xlsx','sheet','ble','range','A40');

writetable(mcmc_summary_exp_msv,'mcmc_results_exp.xlsx','sheet','msv');
writematrix(mcmc_msv_exp.marginal(end,2),'mcmc_results_exp.xlsx','sheet','msv','range','A40');

writetable(mcmc_summary_exp_var1,'mcmc_results_exp.xlsx','sheet','var1');
writematrix(mcmc_var1_exp.marginal(end,2),'mcmc_results_exp.xlsx','sheet','var1','range','A40');

writetable(mcmc_summary_exp_ar2,'mcmc_results_exp.xlsx','sheet','ar2');
writematrix(mcmc_ar2_exp.marginal(end,2),'mcmc_results_exp.xlsx','sheet','ar2','range','A40');
%% rational expectations models -- output files are from Dynare and have a different structure

ree_noExp = load('mcmc_ree.mat');

writematrix(table2array(struct2table(ree_noExp.oo_.posterior_mean.parameters))',...
    'mcmc_results.xlsx','sheet','ree','range','B2');
writematrix(table2array(struct2table(ree_noExp.oo_.posterior_mean.shocks_std))',...
    'mcmc_results.xlsx','sheet','ree','range','B28');

writematrix(table2array(struct2table(ree_noExp.oo_.posterior_hpdinf.parameters))',...
    'mcmc_results.xlsx','sheet','ree','range','C2');
writematrix(table2array(struct2table(ree_noExp.oo_.posterior_hpdinf.shocks_std))',...
    'mcmc_results.xlsx','sheet','ree','range','C28');

writematrix(table2array(struct2table(ree_noExp.oo_.posterior_hpdsup.parameters))',...
    'mcmc_results.xlsx','sheet','ree','range','D2');
writematrix(table2array(struct2table(ree_noExp.oo_.posterior_hpdsup.shocks_std))',...
    'mcmc_results.xlsx','sheet','ree','range','D28');

mhm=ree_noExp.oo_.MarginalDensity.ModifiedHarmonicMean;
writematrix(mhm,'mcmc_results.xlsx','sheet','ree','range','A40');

%% with exp

ree_withExp = load('mcmc_ree_exp.mat');

writematrix(table2array(struct2table(ree_withExp.oo_.posterior_mean.parameters))',...
    'mcmc_results_exp.xlsx','sheet','ree','range','B2');
writematrix(table2array(struct2table(ree_withExp.oo_.posterior_mean.shocks_std))',...
    'mcmc_results_exp.xlsx','sheet','ree','range','B28');

writematrix(table2array(struct2table(ree_withExp.oo_.posterior_hpdinf.parameters))',...
    'mcmc_results_exp.xlsx','sheet','ree','range','C2');
writematrix(table2array(struct2table(ree_withExp.oo_.posterior_hpdinf.shocks_std))',...
    'mcmc_results_exp.xlsx','sheet','ree','range','C28');

writematrix(table2array(struct2table(ree_withExp.oo_.posterior_hpdsup.parameters))',...
    'mcmc_results_exp.xlsx','sheet','ree','range','D2');
writematrix(table2array(struct2table(ree_withExp.oo_.posterior_hpdsup.shocks_std))',...
    'mcmc_results_exp.xlsx','sheet','ree','range','D28');

mhm=ree_withExp.oo_.MarginalDensity.ModifiedHarmonicMean;
writematrix(mhm,'mcmc_results_exp.xlsx','sheet','ree','range','A40');
