 clear;
 clc;
 close all;
 
 clear;clc;%close all;
names_endo=[{'mc' 'zcap' 'rk' 'k' 'q'  'Consumption' 'Investment' 'Output' 'lab' 'Inflation' 'w' 'r' 'kp'   'dy' 'dc' 'dinve' 'dw' 'eps_a'  'eps_b' 'eps_g' 'eps_i'  'eps_r'  'eps_p' 'eps_w' }];
 names_exo=[{ 'Productivity'  'Risk Premium' 'Government Spending' 'Investment'  'Monetary Policy'  'Price Mark-up' 'Wage Mark-up' }];

 load kf_output_msv_t_1_star;



%construct the dataset of RHS variables that go into MSV regressions
%this consists of backward looking variables and shocks 


%BL variables enter with a lag 
data_RHS = S(1:end-1,model.BL_indices)
%shocks enter contemp.
data_RHS = [data_RHS, S(2:end,model.shock_indices)];


%check the correlations 
collintest(data_RHS)

%suggests 1 2 3 and 10 exhibit multicollinearity. These are 
% (6) consumption, investment, output, government spending 

%try again excluding output

model.BL_indices2 = [ 6     7      10    11    12    13]

data_RHS2 = S(1:end-1,model.BL_indices2)
%shocks enter contemp.
data_RHS2 = [data_RHS2, S(2:end,model.shock_indices)];

collintest(data_RHS2)


