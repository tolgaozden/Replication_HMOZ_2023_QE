clear;
clc;

load results/ar1_kf_output.mat;

%period covered: 1965Q1-2008Q4 (176 obs)
S_FL = S(:,model.FL_indices);

%retrieve betas for expectations 
switch model.learning 
    case 0 % retrieve BLE
beta_BLE = diag(beta_tt(model.FL_indices,model.FL_indices));
% beta_BLE = beta_init;

    case 1
end


S_exp       = nan(size(S_FL));
exp_error   = nan(size(S_FL));
for tt=3:length(S_FL)
    
    switch model.learning 
        case 0
S_exp(tt,:) = beta_BLE .* beta_BLE .* S_FL(tt-2,:)' ;
% S_exp_2(tt,:) = diag((diag(beta_BLE).^3)) .* S_FL(tt-2,:)' ;
% S_exp_3(tt,:) = diag((diag(beta_BLE).^4)) .* S_FL(tt-2,:)' ;
% S_exp_4(tt,:) = diag((diag(beta_BLE).^5)) .* S_FL(tt-2,:)' ;
% S_exp_LR(tt,:) = diag((diag(beta_BLE).^40)) .* S_FL(tt-2,:)';

        case 1
            switch model.learning_algo 
                case 'ar(1),sac'
S_exp(tt,:) = alpha_all(:,tt-1) +( beta_all(:,tt-1) .* beta_all(:,tt-1) .* (S_FL(tt-2,:)' - alpha_all(:,tt-1)));
                case 'ar(1),rls'
S_exp(tt,:) = alpha_all(:,tt-1) +( beta_all(:,tt-1) .* beta_all(:,tt-1) .* S_FL(tt-2,:)' + beta_all(:,tt-1) .* alpha_all(:,tt-1));
                case 'var(1)'
                    beta_tmp = zeros(model.numEndo);
                    beta_tmp(model.FL_indices,model.BL_indices) = beta_all(:,:,tt-1); 
                    beta_tmp = beta_tmp^2;
                    beta_tmp = beta_tmp(model.FL_indices,model.BL_indices);
                    
S_exp(tt,:) = alpha_all(:,tt-1) +  beta_tmp *S(tt,model.BL_indices)' + beta_tmp* alpha_all(:,tt-1) ;
            case 'msv'
                
                
%                 a_tmp = alpha_all(6,tt-1);
%                 b_tmp = beta_all(6,:,tt-1);
%                 d_tmp = dd_all(6,:,tt);
%                 S_BL = S(tt-1,model.BL_indices)';
%                 S_shk = S(tt,model.shock_indices)';
%                 
%                 pi_fore = a_tmp  + b_tmp * S_BL + (d_tmp .* diag(sysmat.RHO)') *  S_shk;
%                 
%                 S_exp(tt,6) = pi_fore;

a_tmp = zeros(17,1);
a_tmp(model.FL_indices) = alpha_all(:,tt-2);
bb_tmp = zeros(17,17);
bb_tmp(model.FL_indices,model.BL_indices) = beta_all(:,:,tt-2);
bb_tmp_sq = bb_tmp^2; 
d_tmp = zeros(17,7);
d_tmp(model.FL_indices,:) = dd_all(:,:,tt-2);
S_BL = S(tt-1,model.BL_indices)';
S_shk = S(tt,model.shock_indices)';  
                
fore_tmp = a_tmp + bb_tmp *  S(tt-1,1:17)' + (d_tmp * sysmat.RHO) * S_shk;
% fore_tmp = (a_tmp + bb_tmp * a_tmp)  + bb_tmp_sq * S(tt-1,1:17)' + (d_tmp * sysmat.RHO + bb_tmp * d_tmp) * S_shk;
S_exp(tt,:) = fore_tmp(model.FL_indices);
                
   

            end
    end
%forecast errors 
exp_error(tt,:) = S_FL(tt,:) - S_exp(tt,:);
end




% %omit the first year (remaining sample is 1966Q1-2008Q4)
% S_exp = S_exp(5:end,:);
% exp_error=exp_error(5:end,:);


%look at the period for which SPF is available: 1981Q4 omwards 

% S_exp = S_exp(64:end,:);
% exp_error= exp_error(64:end,:);


%focus on inflation expectations (inflation index=6 in forward looking
%variables)


 
infl_exp         = S_exp(:,6);
% infl_exp_2        = S_exp_2(:,6);
% infl_exp_3        = S_exp_3(:,6);
% infl_exp_4        = S_exp_4(:,6);

infl_exp_error  =exp_error(:,6);


%% comparison with survey of professional forecasters

load SPF_results.mat;

% adjust the sample size if full sample was used 

ll_spf = length(exp_errors_empirical1);

if length(infl_exp) > ll_spf 
    
    infl_exp = infl_exp(end-ll_spf:end-1);
    infl_exp_error = infl_exp_error(end-ll_spf:end-1);
    
end

rmse =mean(sqrt((infl_exp_error(3:end,:) - exp_errors_empirical1(3:end,:)/4).^2),1);
rmse_acf =mean(sqrt((autocorr(infl_exp_error) - autocorr(exp_errors_empirical1/4)).^2),1);



% figure('Name','sample autocorrelation of expectation erorrs');
% for jj=1:7
%     subplot(2,4,jj)
%     autocorr(exp_error(:,jj));
%     title(names_endo(model.FL_indices(jj)));
% end
% print -dpdf figures/sacf_expectation_errors;

pi_bar = x(15) ; %intercept of inflation in the measurement equation



figure('Name','inflation expectations 1-step ahead');
plot((infl_exp+ pi_bar)*4,'b');
hold on;
plot(inf_exp_1_step);%this is annualized, divide by 4
legend('model-implied','empirical','r','lineStyle','--');
% ylim([-1.5 1]);
ylabel('Inflation Expectation')
print -dpdf figures/ar1_exp_1step_inflation

corr_inf= corr((infl_exp+pi_bar)*4,inf_exp_1_step);
disp('correlation between model-implied and empirical inflation expectations (1-step ahead):');
disp(corr_inf);

figure('Name','exp errors');
plot(infl_exp_error,'b');
hold on;
plot(exp_errors_empirical1/4,'r','LineStyle','--');
legend('model-implied','empirical');
str_ = strcat(['rmse of inflation expectations: ', num2str(rmse)]);
title(str_);
% ylim([-1.5 1]);
ylabel('Inflation Expectation Errors')
print -dpdf figures/ar1_exp_errors_inflation


figure('Name','autocorrelation of exp errors');
acf_model = autocorr(infl_exp_error);
acf_emp = autocorr(exp_errors_empirical1/4);
plot(acf_model(2:end),'b');
hold on;
plot(acf_emp(2:end),'r','LineStyle','--');
legend('model-implied','empirical');
str_ = strcat(['rmse of exp error ACF: ', num2str(rmse_acf)]);
title(str_);
% ylim([-0.2 0.35]);
ylabel('ACF of Inflation Expectation Errors')
print -dpdf figures/ar1_sacf_infl_expectation_errors

% figure('Name','distributions of exp errors');
% histogram(infl_exp_error);
% hold on;
% histogram(exp_errors_empirical1/4);
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
