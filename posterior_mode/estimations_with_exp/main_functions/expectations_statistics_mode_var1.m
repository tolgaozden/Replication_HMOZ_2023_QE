clear;
clc;


load results/var1_kf_output.mat;

%period covered: 1965Q1-2008Q3 (175 obs)
S_FL = S(:,model.FL_indices);


% 1965Q1-2008Q3
counter=1:1:size(S_FL,1);
counter=counter';
xx1=datetime(1965,1,1);xx2=datetime(2008,9,1);date_tt=xx1:calmonths(3):xx2;
date_tt=date_tt';

%retrieve betas for expectations 
switch model.learning 
    case 0 % retrieve BLE
beta_BLE = diag(beta_tt(model.FL_indices,model.FL_indices));
% beta_BLE = beta_init;

    case 1
end

% beta_BLE(6)=0.1;

S_exp       = nan(size(S_FL));
exp_error   = nan(size(S_FL));
for tt=3:length(S_FL)
    
    switch model.learning 
        case 0
            
%expectations for period t, made at period t-1, with information available from t-2      
S_exp(tt,:) = beta_BLE .* beta_BLE .* S_FL(tt-2,:)' ;
% S_exp(tt,:) = S(tt-1,end-6:end);

        case 1
            switch model.learning_algo 
                case 'ar(1),sac'
S_exp(tt,:) = alpha_all(:,tt-2) +( beta_all(:,tt-2) .* beta_all(:,tt-2) .* (S_FL(tt-2,:)' - alpha_all(:,tt-2)));
                case 'ar(1),rls'
S_exp(tt,:) = alpha_all(:,tt-2) +( beta_all(:,tt-2) .* beta_all(:,tt-2) .* S_FL(tt-2,:)' + beta_all(:,tt-2) .* alpha_all(:,tt-2));
                case 'var(1)'
                a_tmp = alpha_all(:,tt-2);
                b_tmp = beta_all(:,:,tt-2);
                
  S_exp(tt,:) = a_tmp+ b_tmp * S(tt-2,model.BL_indices)' ;
  
  
                case 'ar(2),rls'
                    if tt>3
  
%            exp1 = -rkp + (alpha_rk + beta1_rk*alpha_rk) + (beta1_rk^2+beta2_rk) *  rkm + (beta1_rk*beta2_rk)*rkm2 == 0;         
a_tmp = squeeze(learning_parameters(tt-2,:,1))';
b1_tmp = squeeze(learning_parameters(tt-2,:,2))';
b2_tmp = squeeze(learning_parameters(tt-2,:,3))';

S_exp(tt,:)= (a_tmp + b1_tmp .* a_tmp) + (b1_tmp.^2 + b2_tmp) .* S_FL(tt-2,:)' + (b1_tmp .* b2_tmp) .* S_FL(tt-3,:)';
                    end
               
  
            case 'msv'
                
                a_tmp = alpha_all(:,tt-2);
                b_tmp = beta_all(:,:,tt-2);
                d_tmp = dd_all(:,1:end-1,tt-2); %exclude the iid exp. shock
                
  S_exp(tt,:) = a_tmp+ b_tmp * S(tt-2,model.BL_indices)' + d_tmp * S(tt-1,model.shock_indices)';              
                
% %                 a_tmp = alpha_all(6,tt-1);
% %                 b_tmp = beta_all(6,:,tt-1);
% %                 d_tmp = dd_all(6,:,tt);
% %                 S_BL = S(tt-1,model.BL_indices)';
% %                 S_shk = S(tt,model.shock_indices)';
% %                 
% %                 pi_fore = a_tmp  + b_tmp * S_BL + (d_tmp .* diag(sysmat.RHO)') *  S_shk;
% %                 
% %                 S_exp(tt,6) = pi_fore;
% 
% a_tmp = zeros(17,1);
% a_tmp(model.FL_indices) = alpha_all(:,tt-2);
% bb_tmp = zeros(17,17);
% bb_tmp(model.FL_indices,model.BL_indices) = beta_all(:,:,tt-2);
% bb_tmp_sq = bb_tmp^2; 
% d_tmp = zeros(17,7);
% d_tmp(model.FL_indices,:) = dd_all(:,:,tt-2);
% S_BL = S(tt-1,model.BL_indices)';
% S_shk = S(tt,model.shock_indices)';  
%                 
% fore_tmp = a_tmp + bb_tmp *  S(tt-1,1:17)' + (d_tmp * sysmat.RHO) * S_shk;
% % fore_tmp = (a_tmp + bb_tmp * a_tmp)  + bb_tmp_sq * S(tt-1,1:17)' + (d_tmp * sysmat.RHO + bb_tmp * d_tmp) * S_shk;
% S_exp(tt,:) = fore_tmp(model.FL_indices);
%                 
   

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
    
    infl_exp = infl_exp(end-ll_spf+1:end);
    infl_exp_error = infl_exp_error(end-ll_spf+1:end);
    
end

pi_bar = x(15) ; %intercept of inflation in the measurement equation

rmse_error =sqrt(mean(((infl_exp_error(1:end,:) - exp_errors_empirical1(1:end,:)/4).^2),1));
% rmse_acf =mean(sqrt((autocorr(infl_exp_error) - autocorr(exp_errors_empirical1/4)).^2),1);

rmse_inf_exp = sqrt(mean(((infl_exp(1:end,:) +pi_bar - inf_exp_1_step(1:end,:)/4 ).^2),1));

% figure('Name','sample autocorrelation of expectation erorrs');
% for jj=1:7
%     subplot(2,4,jj)
%     autocorr(exp_error(:,jj));
%     title(names_endo(model.FL_indices(jj)));
% end
% print -dpdf figures/sacf_expectation_errors;




% 
% figure('Name','inflation expectations 1-step ahead');
% plot((infl_exp(1:end)+ pi_bar)*4,'b');
% hold on;
% plot(inf_exp_1_step(1:end),'r','lineStyle','--');%this is annualized, divide by 4
% hold on;
% % yyaxis right;
% % plot(cpi_annualized(1:end),'--','Color','k');
% 
% legend('model-implied','empirical');
% % str_ = strcat(['rmse of inflation expectations: ', num2str(rmse_inf_exp)]);
% % title(str_);
% ylabel('Inflation Expectation')
% % print -dpdf figures/ble_exp_1step_inflation
% 
% generate_figures('ble_exp_1step_inflation','figures');

corr_inf= corr((infl_exp(1:end)+pi_bar)*4,inf_exp_1_step(1:end));
disp('correlation between model-implied and empirical inflation expectations (1-step ahead):');
disp(corr_inf);
% 
% figure('Name','exp errors');
% plot(infl_exp_error*4,'b');
% hold on;
% plot(exp_errors_empirical1,'r','LineStyle','--');
% legend('model-implied','empirical');
% str_ = strcat(['rmse of inflation expectation errors: ', num2str(rmse_error)]);
% title(str_);
% % ylim([-1.5 1]);
% ylabel('Inflation Expectation Errors')
% % print -dpdf figures/ble_exp_errors_inflation
% generate_figures('ble_exp_errors_inflation','figures');

corr_exp_error = corr(infl_exp_error,exp_errors_empirical1/4);
disp('correlation between model-implied and empirical inflation expectation errors (1-step ahead):');
disp(corr_exp_error);

 writematrix(corr_inf,...
     'inf_exp_correlations.xlsx','sheet','var1','range','B2');
 
 writematrix(corr_exp_error,...
     'inf_exp_correlations.xlsx','sheet','var1','range','B3');



% 
% figure('Name','autocorrelation of exp errors');
% acf_model = autocorr(infl_exp_error);
% acf_emp = autocorr(exp_errors_empirical1/4);
% plot(acf_model(2:end),'b');
% hold on;
% plot(acf_emp(2:end),'r','LineStyle','--');
% legend('model-implied','empirical');
% str_ = strcat(['rmse of exp error ACF: ', num2str(rmse_acf)]);
% title(str_);
% % ylim([-0.2 0.35]);
% ylabel('ACF of Inflation Expectation Errors')
% print -dpdf figures/ar1_sacf_infl_e[scale=0.75]xpectation_errors


% 
% figure('Name','Inflation beta for AR(1) learning model');
% plot(beta_all(6,end-104:end),'b');
% generate_figures('ar1_inflation_beta','figures');
