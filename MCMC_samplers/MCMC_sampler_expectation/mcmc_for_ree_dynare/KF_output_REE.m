clear;
clc;

addpath('SW_Estimation_REE_expData')

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

% for jj=1:keep
% for jj=1:10
    for jj=1:50
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
infl_exp(:,jj) = oo_.SmoothedVariables.pinf_exp_tt1;

    end


infl_exp = 4*(infl_exp + pi_bar);
counter=1:1:size(infl_exp,1);
counter=counter';
xx1=datetime(1965,1,1);xx2=datetime(2008,9,1);date_tt=xx1:calmonths(3):xx2;
date_tt=date_tt';

% save ree_expectation_database.mat;


%sort for each period derive bounds 

mc_lb = 0.01;
mc_ub = 0.99;
mc_med = 0.5;


[num_p,num_s]=size(infl_exp);
mc_indices=round(num_s*[mc_med,mc_ub,mc_lb]);

for jj=1:num_p
   infl_exp_sorted(jj,:) = sort(infl_exp(jj,:));
end

% infl_exp_hpd = [infl_exp_sorted(:,mc_indices)];

infl_exp_hpd = [infl_exp infl_exp infl_exp];


load SPF_results.mat;

l_diff = [length(infl_exp_hpd) - length(inf_exp_1_step)];

inf_exp_1_step=[nan(l_diff,1);inf_exp_1_step];



figure('Name','inflation expectations 1-step ahead');
plotx1(infl_exp_hpd);
hold on;
plot(date_tt,inf_exp_1_step,'color','red');
% legend('model-implied','empirical');
ylim([-1 16]);
% str_ = strcat(['rmse of inflation expectations: ', num2str(rmse_inf_exp)]);
% title(str_);
% generate_figures('ree_exp_1step_inflation','figures');

% corr(inf_fore_1(ll_dif+1:end),infl_exp(ll_dif+1:end))
