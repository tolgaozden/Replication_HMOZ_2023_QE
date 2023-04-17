clear;
clc;

%addpath('/apps/matlab/matlab2020a/dynare-4.6.1/matlab');

addpath('SW_Estimation_REE_expData')
addpath('helpers')

mcmc_part1=load('SW_Estimation_REE_expData/metropolis/SW_Estimation_REE_expData_mh1_blck1.mat');
mcmc_part2=load('SW_Estimation_REE_expData/metropolis/SW_Estimation_REE_expData_mh2_blck1.mat');
info = load('SW_Estimation_REE_expData_mean.mat');
info2=load('SW_Estimation_REE_expData_results.mat');

size1= size(mcmc_part1.x2);
size2= size(mcmc_part2.x2);

param_dist = [mcmc_part1.x2;mcmc_part2.x2];


%discard the first half 

burn=250000;

param_dist = param_dist(burn+1:end,:);


keep=1000;
remaining=length(param_dist);
thin_factor = remaining/keep;


mcmc_for_simul = param_dist(1:thin_factor:end,:);

save ree_mcmc_for_simul.mat;

table(info.parameter_names,mcmc_for_simul(500,:)')
tic

pi_bar = info2.M_.params(7);

for jj=1:keep
% for jj=1:10
%     for jj=1:1
toc
disp(jj)
    %write new parameter vector into .mat file 
  load SW_Estimation_REE_expData_mode.mat;  
  xparam1 = mcmc_for_simul(jj,:)';
  
  save mcmc_param_values.mat;
  
  dynare SW_Estimation_REE_expData_forMCMC.mod noclearall nolog nograph nointeractive fast notmpterms;
 
%% manually compute expectations   
  
%   tmp_mat = oo_.dr.ghx;
%   tmp_mat_inv = tmp_mat(oo_.dr.inv_order_var,:); 
%   infl_coef(:,jj) = tmp_mat_inv(10,:);  
%   state_ind = oo_.dr.state_var; 
%   rhs_ind = M_.endo_names(state_ind);
%   ts_db=[];
%       for kk=1:length(rhs_ind)
%        ts_db=[ts_db, oo_.SmoothedVariables.(rhs_ind{kk})];
%       end
%     infl_exp(:,jj)  = ts_db * infl_coef(:,jj);


%% read from dynare variable 
infl_exp1(:,jj) = oo_.SmoothedVariables.pinf_exp_tt1;
infl_exp2(:,jj) = oo_.SmoothedVariables.pinf_exp_tt2;
infl_exp3(:,jj) = oo_.SmoothedVariables.pinf_exp_tt3;
infl_exp4(:,jj) = oo_.SmoothedVariables.pinf_exp_tt4;
    end

save ree_expectation_database.mat;
% load ree_expectation_database.mat;
infl_exp2 = 4*(infl_exp2 + pi_bar);
counter=1:1:size(infl_exp2,1);
counter=counter';
xx1=datetime(1965,1,1);xx2=datetime(2008,9,1);date_tt=xx1:calmonths(3):xx2;
date_tt=date_tt';

est_database = load('expectations_dataset.mat');

% save ree_expectation_database.mat;


%sort for each period derive bounds 

mc_lb = 0.1;
mc_ub = 0.9;
mc_med = 0.5;


[num_p,num_s]=size(infl_exp2);
mc_indices=round(num_s*[mc_med,mc_ub,mc_lb]);

for jj=1:num_p
   infl_exp_sorted(jj,:) = sort(infl_exp2(jj,:));
end

infl_exp_hpd = [nan(1,3);infl_exp_sorted(2:end,mc_indices)];

% infl_exp_hpd = [infl_exp infl_exp infl_exp];


load SPF_results.mat;

l_diff = [length(infl_exp_hpd) - length(inf_exp_2_step)];
% l_diff2= [length(infl_exp_hpd) - length(est_database.pinfobs)];


inf_exp_1_step=[nan(l_diff,1);inf_exp_1_step];
% pinfobs=[nan(l_diff2,1);est_database.pinfobs];


ff=figure('Name','inflation expectations 1-step ahead');
%ff.Position = [2067 179 709 632];
plotx1(infl_exp_hpd);
hold on;
pp=plot(date_tt,4*est_database.cpi_quarterly(end-175+1:end),'--','color','black','lineWidth',0.5)
hold on;
plot(date_tt,inf_exp_1_step,'color','red','lineWidth',1);
% hold on;

% hold on;
% yyaxis right;

ll=legend('Model-implied, HPD interval','Model-implied, Median','CPI Inflation','SPF Expectations');
%ll.Position=[0.457821830350871 0.771546422400264 0.445703508018660 0.154028441563225];
% legend('model-implied','empirical');
% yyaxis left ;ylim([-1 16]);
% yyaxis right ;
ylim([-1 16]);
ylabel('Annualized q/q %');
set(gca,'FontSize',15)
xlabel('Year');
% str_ = strcat(['rmse of inflation expectations: ', num2str(rmse_inf_exp)]);
% title(str_);
%generate_figures('ree_exp_1step_inflation','figures');

% corr(inf_fore_1(ll_dif+1:end),infl_exp(ll_dif+1:end))
