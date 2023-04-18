clear;
clc;
tic

%% initialize simulation options

opt=struct();
opt.seed_ = 1;
rng(opt.seed_);
% opt.shock_type='all'; %current options: all, supply (only p-shock), demand (only b-shock)
opt.num_grid1 = 100;
opt.num_grid2= 1;
opt.numVar=24;
opt.sim_length= 5000;
opt.burn_in = round(0.8*opt.sim_length);
opt.param_grid1 =  1.4880 ;
opt.param_grid2= linspace(0.001,0.95,opt.num_grid1);
% opt.gain=0.02;  
% opt.parameters_path='/home/res/tolo/BLE_150222/estimations_new_data/results/ble_estimation_results.mat';
opt.parameters_path='results/msv_estimation_results.mat';
% opt.parameters_path='/home/res/tolo/BLE_150222/estimations_new_data/results/ble_estimation_results.mat';
% gain_grid = [0.001,0.005,0.01,0.02];
% gain_grid = [0.005,0.01,0.02];
gain_grid=0.000;
% gain_grid=[0.001,0.005];
% shock_types={'supply','demand','all'};
shock_types={'all'};


for gg=1:length(gain_grid)
    
    for ss=1:length(shock_types)
       
        opt.gain=gain_grid(gg);
        opt.shock_type=shock_types{ss};


%% model options 
opt.model=[];

opt.model.PLM='msv,t';%options: 'msv,t','msv,t-1',ar(1),t,rls','ar(1),t,sac','ar(1),t-1,rls','ar(1),t-1,sac'
opt.model.learning=0; %options: 1=adaptive learning invoked, 0=beliefs fixed at initial values 
opt.model.learning_algo='msv'; %options: 'msv','ar(1),rls','ar(1),sac';
opt.model.initial_beliefs='msv,ree-based';%options: 'msv,ree-based,fixed','ar(1),ree-based,fixed','ar(1),ree-based','ar(1),ble-based','ar(1),ble-based,fixed'
opt.model.PLM_timing='t';%this is repetitive, same info in model.PLM, fix it later
opt.model.projection_facility=0;%=1 impose PR facility, =0 do not
opt.model.projection_facility_PLM=0;%=1 imposes projection facility on the PLM on top of IALM
%

opt.model.l=7;%dataset # of variables
opt.model.burnIn=4;%presample data for KF initialization
opt.model.numVar=24;% #model size
opt.model.numShocks=7;% #number of shocks
opt.model.numEndo=17;% #number of endogenous variables
opt.model.numExo=7;% #number of exogenous variables
opt.model.numBackward=7;% #number of backward-looking variables
opt.model.numForward=7;% #number of forward-looking  variables
opt.model.BL_indices=[6 7 8 10 11 12 13];% backward-looking variables indices in PLM
opt.model.FL_indices=[3 5 6 7 9 10 11];% forward-looking  variables indices 
opt.model.shock_indices=18:1:24;
% model.shock_indices_lagged=[1];
opt.model.LB=[0.5 0.25 0.001 0.3 0.25 0.5 0.01 0.01 0.01 1 1 0.5 0.001 0.001 0.1 0.01 -10 0.1 0.01 0.01 0.01 0.019 0.01 0.01 ...
    0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.00001]';%lower bound on estimates
opt.model.UB=[15 3 0.999 0.95 10 0.99 0.99 0.99 1 3 3 0.999 0.5 0.5 2 2 10 0.8 1 0.9999 .9999 .9999 .9999 .9999 .9999 ...
    .9999 .9999 10 10 10 10 10 10 10 0.1]';%upper bound on estimates
opt.model.ridge_correction=0;%invokes ridge correction when updating beliefs if =1
opt.model.compute_init_hess=0;
opt.model.sim_length=opt.sim_length;
opt.model.num_simul=1;
opt.model.learning_mean=1;
opt.model.calculate_eig=0;
opt.model.burn_in=0000;
opt.model.hetero_gains=0;
opt.model.N_fixedPoint=200;
opt.model.optimizer=2; %options: 1=fminsearch, 2=csminwel 3=fmincon

%%



SS=zeros(opt.sim_length,opt.numVar);
SS_all=nan(opt.sim_length,opt.numVar,opt.num_grid1,opt.num_grid2);
var_all = nan(opt.numVar,opt.num_grid1,opt.num_grid2);
mean_all = nan(opt.numVar,opt.num_grid1,opt.num_grid2);


parfor nn=1:opt.num_grid1
% for nn=1:opt.num_grid1
    disp(nn);
    for mm=1:opt.num_grid2
 
       sim_index = nn;
        
        param1 = opt.param_grid1;
        param2 = opt.param_grid2(nn);

%       try
      [sim_results] = simulation_main_ree(param1,param2,opt); 
%       catch
%        sim_results=sim_results_prev;
%       end
      
      
results_tmp = sim_results; 


%       
SS_tmp = sim_results.S;


    SS_all(:,:,nn) = SS_tmp ; 
    var_all(:,nn) = var(SS_tmp(opt.burn_in:end,:));
    mean_all(:,nn) = mean(SS_tmp(opt.burn_in:end,:));
    sim_beta(:,nn) = diag(sim_results.beta(opt.model.FL_indices,opt.model.FL_indices));

  
sim_results_prev = sim_results;    
    
    end
end

inflation_var = var_all(10,:); 
output_gap_var =var_all(8,:);
interest_rate_var = var_all(12,:);



toc



%==================objective function 
% weight_pi_grid=[1 1 1 1 1 1 1 1 ];
% weight_y_grid=[0.1 0.048 0.25 0.2 0.15 0.1 0.05 0];
% weight_r_grid=[0.05 0.236 0 0.05 0.09 0.15 0.2 0.25];

% weight_pi_grid = [1  0  0    1   1      1    1     1  1 ];
% weight_y_grid =  [0  1  0    0.1  0.2    0.2  0.2   0.2 0.048];
% weight_pi_grid = [1 1 1 1 1 1 1 1 1 ];
% weight_y_grid = 0.048 * weight_pi_grid ; 
% weight_r_grid =  [0 0.001 0.05 0.1 0.2 0.236 0.5 0.75 1];

% num_policy = 20;

% weight_pi_grid = ones(num_policy,1);
% weight_y_grid = linspace(0,0.5,num_policy)';
% weight_r_grid = linspace(0,0,num_policy)';
% 
% weight_pi_grid = ones(num_policy,1);
% weight_y_grid = 0.048*ones(num_policy,1);
% weight_r_grid = linspace(0,0.5,num_policy)';


weight_pi_grid = [1    1         1 ]';
weight_y_grid = [0.048 0.048   0.1 ]';
weight_r_grid = [0     0.1      0.1 ]';
optimalPara_SAC = nan(1,length(weight_pi_grid));
welfare_SAC = nan(length(weight_pi_grid),1);

% output_gap_var = [nan diff(output_gap_var)];

for ll=1:length(weight_pi_grid)
    
   weight_pi=weight_pi_grid(ll);
   weight_r=weight_r_grid(ll);
   weight_y=weight_y_grid(ll);
   
   

objective_SAC=weight_pi * inflation_var + weight_y * output_gap_var + weight_r * interest_rate_var;


minSAC=min((objective_SAC));
[minSAC_ind]=find(objective_SAC==minSAC);

optimalPara_SAC(ll)=[opt.param_grid2(minSAC_ind(1))];

welfare_SAC(ll)=objective_SAC(minSAC_ind(1));



% table_output(ll,:)=[weight_pi,weight_y,weight_r,optimalPara_BLE,...
% inflation_var(minBLE_ind(1),minBLE_ind(2)),output_gap_var(minBLE_ind(1),minBLE_ind(2)),...
% interest_rate_var(minBLE_ind(1),minBLE_ind(2)),welfare_BLE(ll)];

end


        

disp([weight_pi_grid,weight_y_grid,weight_r_grid,optimalPara_SAC'])

T_REE = table(weight_pi_grid,weight_y_grid,weight_r_grid,optimalPara_SAC');
T_REE.Properties.VariableNames{4} = 'Optimal_Smoothing';
save T_REE.mat T_REE;


file_name=['codes/simulations/',opt.model.initial_beliefs,'_gain_',num2str(opt.gain),'_learning','_',num2str(opt.model.learning),'_',opt.shock_type,'_smoothing.mat']

save(file_name);

    end
end






% save opt_policy_SAC_1dim_rhor_gain001_final_rule2;
% save opt_policy_SAC_1dim_rhor_gain001_weighty_rule2;

% save opt_policy_AR1_1dim_rhor_3_final;
% save opt_policy_SAC_1dim_rhor_gain001_final_weighty;

% save opt_policy_SAC_1dim_rhor_gain0_01_final;
% save opt_policy_SAC_1dim_rhor_gain001_final_weighty2;

% save opt_policy_SAC_1dim_rhor_gain001_final_weighty2;
% save opt_policy_SAC_1dim_rhor_gain001_weighty2_rule2;

% figure;
% plot(inflation_var);

% figure;
% plot(optimalPara_SAC(1,:))


% figure('Name','variances');
% subplot(1,3,1);plot(rho_r_grid,inflation_var);title('inflation');
% subplot(1,3,2);plot(rho_r_grid,output_gap_var);title('output growth');
% subplot(1,3,3);plot(rho_r_grid,interest_rate_var);title('interest rate');

% figure('Name','averages');
% subplot(1,3,1);plot(rho_r_grid,inflation_mean);title('inflation');
% subplot(1,3,2);plot(rho_r_grid,output_gap_mean);title('output growth');
% subplot(1,3,3);plot(rho_r_grid,interest_rate_mean);title('interest rate');


% figure;
% subplot(1,3,1);plot(squeeze(SS_all(end,1,burn_in:end,10)));title('inflation');
% subplot(1,3,2);plot(squeeze(SS_all(end,1,burn_in:end,14)));title('output growth');
% subplot(1,3,3);plot(squeeze(SS_all(end,1,burn_in:end,12)));title('interest rate');