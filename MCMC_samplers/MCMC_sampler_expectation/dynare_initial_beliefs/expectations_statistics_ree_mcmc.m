clear;


dynare SW_Estimation_REE_expData.mod noclearall;

num_simul=1000;
pi_bar = M_.params(7);

% for jj=1:num_simul
for jj=1
    
    %set parameter values and load in the mod file --- to be done 
    
    tmp_mat = oo_.dr.ghx;
    
    %switch to declaration order. 10=inflation in tmp_mat_inv 
    tmp_mat_inv = tmp_mat(oo_.dr.inv_order_var,:); 
    
    infl_coef(:,jj) = tmp_mat_inv(10,:);
    
    %state variables in dynare order 
    state_ind = oo_.dr.state_var; 
    
    rhs_ind = M_.endo_names(state_ind);
    
    ts_db=[];
    for kk=1:length(rhs_ind)
        
        ts_db=[ts_db, oo_.SmoothedVariables.(rhs_ind{kk})];
        
    end
    
  infl_exp(:,jj)  = ts_db * infl_coef(:,jj);
    
end
%convert t+1 expectations to exp. for t

% infl_exp = infl_exp(1:end-1,:);
% 
% infl_exp = [nan(size(infl_exp,2),1);infl_exp];

%convert to annualized and add back infl trend

infl_exp = 4*(infl_exp + pi_bar);


counter=1:1:size(infl_exp,1);
counter=counter';
xx1=datetime(1965,1,1);xx2=datetime(2008,9,1);date_tt=xx1:calmonths(3):xx2;
date_tt=date_tt';

load SPF_results.mat;

ll_dif =size(infl_exp,1) - size(inf_fore_1,1);

inf_fore_1 = [nan(ll_dif,1);inf_fore_1];

figure('Name','inflation expectations 1-step ahead');
plot(date_tt,infl_exp,'b');
hold on;
plot(date_tt,inf_fore_1);
legend('model-implied','empirical');
% str_ = strcat(['rmse of inflation expectations: ', num2str(rmse_inf_exp)]);
% title(str_);
% generate_figures('ree_exp_1step_inflation','figures');

corr(inf_fore_1(ll_dif+1:end),infl_exp(ll_dif+1:end))