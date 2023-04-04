
clear;
clc;
% close all;
% load kf_output_estimation_results_ar2_t_1_full.mat;
% load kf_output;
load results/msv_kf_output;
load MSV_posteriorDist.mat;
% pi_bar = x(15);


% autocorr(S(2:end,10) - beta_tt(10,10) * S(1:end-1,10))

%filtered data for forward-looking variaar2s 

%period covered: 1965Q1-2008Q4 (176 obs)
S_FL = S_mcmc(:,model.FL_indices,:);

% 1965Q1-2008Q3
counter=1:1:size(S_FL,1);
counter=counter';
xx1=datetime(1965,1,1);xx2=datetime(2008,9,1);date_tt=xx1:calmonths(3):xx2;
date_tt=date_tt';

%retrieve betas for expectations 




%Perceived inflation for period t = beta * S_{t-1} 
S_exp       = nan(size(S_FL));
exp_error   = nan(size(S_FL));

for mcmc=1:size(S_mcmc,3)
disp(mcmc)



for tt=3:size(S_mcmc,1)
    switch model.learning 
        case 0
            
%expectations for period t, made at period t-1, with information availaar2 from t-2      
S_exp(tt,:,mcmc) = beta_BLE(:,mcmc) .* beta_BLE(:,mcmc) .* S_FL(tt-2,:,mcmc)' ;
% S_exp(tt,:) = S(tt-1,end-6:end);

        case 1
            switch model.learning_algo 
                case 'ar(1),sac'
S_exp(tt,:,mcmc) = alpha_all(:,tt-2,mcmc) +( beta_all(:,tt-2,mcmc) .* beta_all(:,tt-2,mcmc) .* (S_FL(tt-2,:,mcmc)' - alpha_all(:,tt-2,mcmc)));
                case 'ar(1),rls'
S_exp(tt,:,mcmc) = alpha_all(:,tt-2,mcmc) +( beta_all(:,tt-2,mcmc) .* beta_all(:,tt-2,mcmc) .* S_FL(tt-2,:,mcmc)' + beta_all(:,tt-2,mcmc) .* alpha_all(:,tt-2,mcmc));
                case 'var(1)'
                a_tmp = alpha_all(:,tt-2,mcmc);
                b_tmp = beta_all(:,:,tt-2,mcmc);
                
  S_exp(tt,:,mcmc) = a_tmp+ b_tmp * S_mcmc(tt-2,model.BL_indices,mcmc)' ;
  
  
                case 'ar(2),rls'
                    if tt>3 
  
%            exp1 = -rkp + (alpha_rk + beta1_rk*alpha_rk) + (beta1_rk^2+beta2_rk) *  rkm + (beta1_rk*beta2_rk)*rkm2 == 0;         
a_tmp = squeeze(learning_parameters(tt-2,:,1,mcmc))';
b1_tmp = squeeze(learning_parameters(tt-2,:,2,mcmc))';
b2_tmp = squeeze(learning_parameters(tt-2,:,3,mcmc))';

S_exp(tt,:,mcmc)= (a_tmp + b1_tmp .* a_tmp) + (b1_tmp.^2 + b2_tmp) .* S_FL(tt-2,:,mcmc)' + (b1_tmp .* b2_tmp) .* S_FL(tt-3,:,mcmc)';
                    end
               
  
            case 'msv'
                
                a_tmp = alpha_all(:,tt-2,mcmc);
                b_tmp = beta_all(:,:,tt-2,mcmc);
                d_tmp = dd_all(:,1:end-1,tt-2,mcmc); %exclude the iid exp. shock
                
  S_exp(tt,:,mcmc) = a_tmp+ b_tmp * S_mcmc(tt-2,model.BL_indices,mcmc)' + d_tmp * S_mcmc(tt-1,model.shock_indices,mcmc)';              



            end
    end
    
    
    
exp_error(tt,:,mcmc) = S_FL(tt,:,mcmc) - S_exp(tt,:,mcmc);
end



end

pi_bar_mcmc = param_mcmc(:,15);

% %omit the first year (remaining sample is 1966Q1-2008Q4)
% S_exp = S_exp(5:end,:);
% exp_error=exp_error(5:end,:);


%look at the period for which SPF is availaar2: 1981Q4 omwards 

% S_exp = S_exp(64:end,:);
% exp_error= exp_error(64:end,:);


%focus on inflation expectations (inflation index=6 in forward looking
%variaar2s)
 
infl_exp        = squeeze(S_exp(:,6,:));
infl_exp_error  =squeeze(exp_error(:,6,:));




load SPF_results.mat;



%% derive mcmc bands

%add back inflation trend 
for jj=1:mcmc
    infl_exp(:,jj)= infl_exp(:,jj) + pi_bar_mcmc(jj);
end

infl_exp_error_sorted = sort(infl_exp_error,2);
infl_exp_sorted =sort(infl_exp,2);


%derive the HPD band
mc_lb = 0.01;
mc_ub = 0.99;
mc_med = 0.5;

mc_lb = round(mc_lb * mcmc);
mc_ub = round(mc_ub * mcmc);
mc_med = round(mc_med * mcmc);

infl_exp_error_lb = infl_exp_error_sorted(:,mc_lb);
infl_exp_error_med = infl_exp_error_sorted(:,mc_med);
infl_exp_error_ub = infl_exp_error_sorted(:,mc_ub);

infl_exp_lb = infl_exp_sorted(:,mc_lb);
infl_exp_med = infl_exp_sorted(:,mc_med);
infl_exp_ub = infl_exp_sorted(:,mc_ub);

infl_data_exp_error = [infl_exp_error_med infl_exp_error_ub infl_exp_error_lb];
infl_data = 4*[infl_exp_med infl_exp_ub infl_exp_lb];
%% append empirical series with nans 

l_diff = size(infl_data,1) - size(exp_errors_empirical1,1);
exp_errors_empirical1 = [nan(l_diff,1);exp_errors_empirical1];
inf_exp_1_step = [nan(l_diff,1);inf_exp_1_step];
%% compute correlation distribution 

for jj=1:mcmc 
  corr_inf_exp_error(jj) =   corr(exp_errors_empirical1(l_diff+1:end),infl_exp_error(l_diff+1:end,jj));
  corr_inf_exp(jj) =   corr(inf_exp_1_step(l_diff+1:end),infl_exp(l_diff+1:end,jj));
end

corr_inf_exp_error = sort(corr_inf_exp_error);
corr_inf_exp_error = corr_inf_exp_error([mc_lb mc_med mc_ub]);

corr_inf_exp = sort(corr_inf_exp);
corr_inf_exp = corr_inf_exp([mc_lb mc_med mc_ub]);


% writematrix(corr_inf_exp(2),...
%     'inf_exp_correlations.xlsx','sheet','msv','range','B2');
% 
% writematrix(corr_inf_exp_error(2),...
%     'inf_exp_correlations.xlsx','sheet','msv','range','C2');



%% figures 


ff2=figure('Name','inflation expectations, model-implied vs. empirical');
% plotx1(infl_data_exp_error,'--','color','blue');
plotx1((infl_data));
hold on;
plot(date_tt,inf_exp_1_step,'color','red');
%legend('model-implied, HPD interval','model-implied, median','empirical');
% print -dpdf figures/ar2_inf_exp_mcmc
ylim([-1 16]);
ylabel('Annualized q/q %');
set(gca,'FontSize',15)
xlabel('Year');

generate_figures('msv_inf_exp_mcmc','figures');
