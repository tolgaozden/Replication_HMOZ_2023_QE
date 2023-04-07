clear;
clc;

load dynare_initial_beliefs/estimation_results_SW_Estimation_REE_expData.mat;



pi_bar = M_.params(7);

%model-implied 1-step expectations


%coefficients 

alpha_re = 0; %zero intercept
bb_y = 0.014120 ;
bb_r = 0.073144 ;
bb_kp =  -0.001427 ;
bb_c =  0.002143 ;
bb_inve = 0.001188  ;
bb_pinf =0.060088 ;
bb_w =  0.008444 ;

dd_a = 0.060554 ;
dd_b =   -0.000264  ;
dd_g =  -0.024590  ;
dd_i = 0.007394  ;
dd_r =  -0.003595 ;
dd_p=  0.738096  ;
dd_w=  1.117091   ;

sm_y = oo_.SmoothedVariables.y; 
sm_r = oo_.SmoothedVariables.r; 
sm_kp = oo_.SmoothedVariables.kp; 
sm_c = oo_.SmoothedVariables.c; 
sm_inve = oo_.SmoothedVariables.inve; 
sm_pinf = oo_.SmoothedVariables.pinf; 
sm_w = oo_.SmoothedVariables.w; 
sm_a= oo_.SmoothedVariables.eps_a;
sm_b= oo_.SmoothedVariables.eps_b;
sm_g= oo_.SmoothedVariables.eps_g;
sm_i= oo_.SmoothedVariables.eps_i;
sm_r= oo_.SmoothedVariables.eps_r;
sm_p= oo_.SmoothedVariables.eps_p;
sm_w= oo_.SmoothedVariables.eps_w;


%forecasts for t+1 made at t
pinf_exp_re = sm_y * bb_y + ...
          sm_r * bb_r + ... 
          sm_kp * bb_kp + ... 
          sm_c * bb_c + ... 
          sm_inve * bb_inve + ... 
          sm_pinf * bb_pinf + ... 
          sm_w * bb_w + ... 
          sm_a * dd_a + ... 
          sm_b * dd_b + ... 
          sm_g * dd_g + ... 
          sm_i * dd_i + ... 
          sm_r * dd_r + ... 
          sm_p * dd_p + ... 
          sm_w * dd_w ;

      
pinf_exp_re = 4*(pi_bar+pinf_exp_re);
      
% pinf_exp_re = (oo_.SmoothedVariables.pinf_exp-oo_.SmoothedShocks.eta_pi_exp+pi_bar)*4;
% pinf_exp_re = (oo_.SmoothedVariables.pinf_exp-oo_.SmoothedShocks.eta_pi_exp);
load SPF_results.mat;


ll_spf = length(exp_errors_empirical1);

pinf_exp_re = pinf_exp_re(end-ll_spf+1:end);


%rmse for inflation expectation errors 
pinf_model = 4*(pi_bar+oo_.SmoothedVariables.pinf(end-ll_spf+1:end));
%substract forecasts for period t made at t-1, from realized pinf at t
fore_error_model = pinf_exp_re(1:end-1) - pinf_model(2:end);


rmse_error =sqrt(mean(((fore_error_model(1:end,:)/4 - exp_errors_empirical1(1:end-1,:)/4).^2),1));
rmse_inf_exp=sqrt(mean(((pinf_exp_re(1:end,:)/4 - inf_fore_1(1:end,:)/4).^2),1));


% figure('Name','inflation expectations 1-step ahead');
% plot(pinf_exp_re,'b');
% hold on;
% plot(inf_fore_1);
% legend('model-implied','empirical');
% % str_ = strcat(['rmse of inflation expectations: ', num2str(rmse_inf_exp)]);
% % title(str_);
% generate_figures('ree_exp_1step_inflation','figures');
% 

corr_inf=corr(pinf_exp_re,inf_fore_1);
disp('correlation between model-implied and empirical inflation expectations (1-step ahead):');
disp(corr_inf);

% figure('Name','exp errors');
% plot(fore_error_model,'b');
% hold on;
% plot(exp_errors_empirical1(1:end-1),'r','LineStyle','--');
% str_ = strcat(['rmse of inflation expectation errors: ', num2str(rmse_error)]);
% title(str_);
% generate_figures('ree_exp_errors_inflation','figures');

corr_exp_error=corr(fore_error_model,exp_errors_empirical1(1:end-1));
disp('correlation between model-implied and empirical inflation expectation errors (1-step ahead):');
disp(corr_exp_error);



writematrix(corr_inf,...
    'inf_exp_correlations.xlsx','sheet','ree','range','B2');

writematrix(corr_exp_error,...
    'inf_exp_correlations.xlsx','sheet','ree','range','B3');


