clear;
clc;


% inputs: 
%estimation results at the posterior mode
%should be stored in the subfolder optimal_policy/codes/results

% outputs: 

% the results with optimal policy parameters are printed out on the command
% window when simulations are finished running. 

% the figure with standard deviations as a function of smoothing is stored
% under optimal_policy/figures 

% the output files are stored in db_[model_name]_smoothing.mat in
% optimal_policy folder. 




cd('optimal_policy');
addpath(genpath(cd));

optimal_policy_smoothing_SAC;
optimal_policy_smoothing_BLE;
optimal_policy_smoothing_MSV;
optimal_policy_smoothing_REE;


reports_sac_smoothing;
reports_ble_smoothing;
reports_msv_smoothing;
reports_ree_smoothing;


smoothing_plots_all;

%% optimal parameters 

load T_BLE ;
load T_SAC;
load T_MSV ;
load T_REE;

disp('Optimal smoothing under BLE');
disp(T_BLE)

disp('Optimal smoothing under SAC');
disp(T_SAC)

disp('Optimal smoothing under MSV');
disp(T_MSV)

disp('Optimal smoothing under REE');
disp(T_REE)