clear;clc;close all;
 forecast.num_periods=96;


for jj=1:forecast.num_periods
file_name=['auxiliary_forecast_model_ree' num2str(jj) '.mat'];
location=[cd,'\auxiliary_files\'];
string=[location file_name];

load(string);
save estimation_results;
KF_output;
second_moments_init=...
diag(gam1(model.FL_indices,model.FL_indices))./...
diag(gam0(model.FL_indices,model.FL_indices));
file_name=['ble_init_auxiliary_moments' num2str(jj) '.mat'];
location=[cd,'\auxiliary_files\'];
string=[location file_name];
clearvars -except string second_moments_init;
save(string);
end