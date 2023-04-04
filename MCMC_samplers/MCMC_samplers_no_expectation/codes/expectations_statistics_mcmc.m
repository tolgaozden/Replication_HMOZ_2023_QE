
clear;
clc;
% close all;
% load kf_output_estimation_results_ble_t_1_full.mat;
% load kf_output;
load BLE_posteriorDist.mat;

% autocorr(S(2:end,10) - beta_tt(10,10) * S(1:end-1,10))

%filtered data for forward-looking variables 


%period covered: 1965Q1-2008Q4 (176 obs)
S_FL = S_mcmc(:,model.FL_indices,:);

%retrieve betas for expectations 




%Perceived inflation for period t = beta * S_{t-1} 
S_exp       = nan(size(S_FL));
exp_error   = nan(size(S_FL));

for mcmc=1:size(S_mcmc,3)
disp(mcmc)
for tt=3:size(S_mcmc,1)
S_exp(tt,:,mcmc) = beta_BLE(:,mcmc) .* beta_BLE(:,mcmc) .* S_FL(tt-2,:,mcmc)' ;

%forecast errors 
exp_error(tt,:,mcmc) = S_FL(tt,:,mcmc) - S_exp(tt,:,mcmc);
end


end

% %omit the first year (remaining sample is 1966Q1-2008Q4)
% S_exp = S_exp(5:end,:);
% exp_error=exp_error(5:end,:);


%look at the period for which SPF is available: 1981Q4 omwards 

% S_exp = S_exp(64:end,:);
% exp_error= exp_error(64:end,:);


%focus on inflation expectations (inflation index=6 in forward looking
%variables)
 
infl_exp        = squeeze(S_exp(:,6,:));
infl_exp_error  =squeeze(exp_error(:,6,:));


% figure;
% for jj=1:7
%     subplot(2,4,jj)
%     autocorr(exp_error(:,jj));
%     title(names_endo(model.FL_indices(jj)));
% end

% pi_bar = x(15) ; %intercept of inflation in the measurement equation

% clearvars -except S S_exp infl_exp infl_exp_error pi_bar; 

load SPF_results.mat;



%% derive mcmc bands

infl_exp_error_sorted = sort(infl_exp_error,2);

%derive the HPD band
mc_lb = 0.001;
mc_ub = 1;
mc_med = 0.5;

mc_lb = round(mc_lb * mcmc);
mc_ub = round(mc_ub * mcmc);
mc_med = round(mc_med * mcmc);

infl_exp_lb = infl_exp_error_sorted(:,mc_lb);
infl_exp_med = infl_exp_error_sorted(:,mc_med);
infl_exp_ub = infl_exp_error_sorted(:,mc_ub);

infl_data = [infl_exp_med infl_exp_ub infl_exp_lb];


%% compute RMSE statistics for inflation expectations 


% diff_exp =mean(sqrt((infl_exp_error(3:end,:) - exp_errors_empirical1(3:end,:)/4).^2),1);
% figure('Name','RMSE distribution of Inflation Expectations');
% histfit(diff_exp)

%% figures 
figure('Name','inflation expectations, model-implied vs. empirical');
plotx1(infl_data);
hold on;
plot(exp_errors_empirical1/4);
legend('model-implied, HPD interval','model-implied, median','empirical');
print -dpdf figures/exp_errors_inflation_mcmc


figure('Name','autocorrelation of exp errors');
acf_model = [autocorr(infl_exp_med) autocorr(infl_exp_ub) autocorr(infl_exp_lb)];
acf_emp = autocorr(exp_errors_empirical1/4);
plotx1(acf_model(2:end,:));
hold on;
plot(acf_emp(2:end));
legend('model-implied, HPD interval','model-implied, median','empirical');
print -dpdf figures/sacf_exp_errors_inflation_mcmc;


% figure('Name','inflation expectations 1-step ahead');
% plot((infl_exp+ pi_bar)*4);
% hold on;
% plot(inf_exp_1_step);%this is annualized, divide by 4
% legend('model-implied','empirical');


% 
% % clearvars -except infl_exp infl_exp_error S_exp exp_error beta_BLE;
% 
% %load the results from SPF, these are computed in the folder
% %"survey_expectations" in simple_regressions.m
% 
% % load SPF_results.mat;
% 
% diff_exp = pinfobs_exp_median - mean(pinfobs_exp_median) - 4*infl_exp;
% acf_diff_exp = autocorr(diff_exp);
% 
% diff_exp_error = exp_errors_empirical - 4*infl_exp_error;
% acf_diff_exp_error = autocorr(diff_exp_error);
% 
% figure;
% subplot(2,1,1)
% plot(pinfobs_exp_median - mean(pinfobs_exp_median),'color','black');
% hold on;
% plot(4*(infl_exp),'color','red');
% legend('SPF - mean inflation expectations','BLE model-implied expectations');
% subplot(2,1,2);
% plot((acf_diff_exp(2:end)),'color','black');
% legend('ACF of the Difference');
% 
% figure;
% autocorr(diff_exp);
% 
% 
% figure;
% subplot(2,1,1);
% plot(exp_errors_empirical,'color','black');
% hold on;
% plot(4*infl_exp_error,'color','red');
% legend('SPF - expectation errors','BLE model-implied expectation errors');
% subplot(2,1,2);
% plot((acf_diff_exp_error(2:end)),'color','black');
% legend('ACF of the Difference');
% 
% figure;
% autocorr(diff_exp_error);
% 
% %calculate RMSE 
% 
% rmse_ble = sqrt(mean(diff_exp.^2))
% 
% 
