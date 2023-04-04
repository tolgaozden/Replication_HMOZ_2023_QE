clear;clc;close all;
clear;clc;%close all;
restoredefaultpath;
addpath('helpers');
addpath('Matfiles');
addpath('figures');
addpath('dynare_initial_beliefs');
addpath('inputs');
addpath('results');
tic
warning('off');
names=[{'phi','sigma_c' ,'lambda' ,'xi_w' ,'sigma_l', 'xi_p' , 'iota_w','iota_p',...
    'psi' ,'phi_p','r_pi', 'rho' ,'r_y', 'r_dy' ,...
           'pi_bar' ,'beta_const' ,'l_bar', 'gamma_bar' ,'alpha'...
           'rho_a', 'rho_b' ,'rho_g' ,'rho_i' ,'rho_r', 'rho_p', 'rho_w','rho_ga',...
           'eta_a', 'eta_b' ,'eta_g' ,'eta_i' ,'eta_r', 'eta_p', 'eta_w',...
           'gain'} ] ;  
%retrieve data
load('dynare_initial_beliefs/slobodyan_dataset.mat');
  %forecast-variables
  forecast=[];

 
  forecast.numVar=24;%same as model.numVar below
  forecast.numObs=7;
  forecast.first_forecast=167;%167;%151; %index of the first variable to forecast
  forecast.horizon=12; %how many steps ahead in each forecast--> h-step ahead gives [1,2...h]
  forecast.window_length=80;%sample size; rolling window length for the estimation period.
  forecast.num_periods=length(dy)-forecast.first_forecast+1;%how many periods do we want to forecast

  
  
  S_fore_all= nan(2*forecast.numVar,forecast.horizon,forecast.num_periods);%first element 2*forecast.numVar only if ar(2) with t-1
  obs_fore_all=nan(forecast.numObs,forecast.horizon,forecast.num_periods);

   %construct forecast dataset: this is the 1-step ahead dataset to be used in forecast evaluation
  %extend to h-step ahead later
  first_obs=forecast.first_forecast;
  last_obs=forecast.first_forecast+forecast.num_periods-1;
 forecast.dataset=[dy dc dinve dw cpi_quarterly robs labobs];
 forecast.dataset=forecast.dataset(first_obs:last_obs,:);
 
    save forecast_output.mat S_fore_all obs_fore_all forecast;
  %---------------------------------------------------------------------------
  %model variables that will not change throughout forecast
 model=[];
model.l=7;%dataset # of variables
model.burnIn=4;%presample data for KF initialization
model.numVar=24;% #model size
model.numShocks=7;% #number of shocks
model.numEndo=17;% #number of endogenous variables
model.numExo=7;% #number of exogenous variables
model.numBackward=7;% #number of backward-looking variables
model.numForward=7;% #number of forward-looking  variables
model.BL_indices=[6 7 8 10 11 12 13];% backward-looking variables indices in PLM
model.FL_indices=[3 5 6 7 9 10 11];% forward-looking  variables indices 
model.shock_indices=18:1:24;
model.shock_indices_lagged=[1 2 3 4 5 6 7];
model.LB=[0.5 0.1 0.001 0.5 0.25 0.5 0.01 0.01 0.01 1 1 0.5 0.001 0.001 0.1 0.01 -10 0.1 0.01 0.01 0.01 0.01 0.01 0.01 ...
    0 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.00001]';%lower bound on estimates
model.UB=[15 3 0.99 0.95 10 0.99 0.99 0.99 1 3 3 0.999 0.5 0.5 2 2 10 0.8 1 0.999 .999 .999 .999 .999 .999 ...
    .999 .999 10 10 10 10 10 10 10 0.1]';%upper bound on estimates
model.PLM='ar(2),t-1,rls';%options: 'msv,t','msv,t-1',ar(1),t,rls','ar(1),t,sac','ar(1),t-1,rls','ar(1),t-1,sac','ar(2),t-1,rls'
model.learning=1; %options: 1=adaptive learning invoked, 0=beliefs fixed at initial values 
model.learning_algo='ar(2),rls'; %options: 'msv','ar(1),rls','ar(1),sac','msv,restricted';
model.initial_beliefs='ar(2),ree-based,fixed';%options: 'msv,ree-based,fixed','ar(1),ree-based,fixed','ar(1),ree-based','ar(1),ble-based','ar(1),ble-based,fixed'
model.PLM_timing='t-1';%this is repetitive, same info in model.PLM, fix it later
model.projection_facility=0;%=1 impose PR facility, =0 do not
model.projection_facility_PLM=0;%=1 imposes projection facility on the PLM on top of IALM
model.N_fixedPoint=100;
model.optimizer=4; %options: 1=fminsearch, 2=csminwel 3=fmincon 4=patternsearch
model.ridge_correction=0;%invokes ridge correction when updating beliefs if =1
model.compute_init_hess=0;%=0 no hessian, =1 only diagonal, =2 full hessian
model.hetero_gains=0;%!!!!!keep at =0 for now, =1 not working
model.decreasing_gain=0;
model.rolling_horizon=1;
%---------------------------------------------------------------------------

 %---------------------------------------------------------------------------
 for jj=1:forecast.num_periods
% for jj=1:22
% for jj=1:1

% for jj=1:5


%      file_name=['auxiliary_forecast_model_ble' num2str(jj) '.mat']
% location=[cd,'\auxiliary_files\']
% string=[location file_name]
% load(string);
% param_init=x;
% param_init(25)=0.01;
% save param_init.mat param_init;
% init_H=H;
% if jj>7
% hh1=init_H(1:5,1:5);hh2=init_H(6:23,6:23);hh3=init_H(24:end,24:end);
% init_H(1:5,1:5)=hh1;init_H(6,6)=1;init_H(7:24,7:24)=hh2;init_H(25,25)=1;init_H(26:34,26:34)=hh3;init_H(35,35)=1;
% save init_H.mat init_H;
% end
% 
% forecast.horizon=12;
     %large price stickiness, small indexation &persistence to make sure
     %eqm exists
%       load inputs/param_init;
%       param_init(6)=0.75;
%      param_init(25)=0.2;
%      param_init(8)=param_init(8)-0.05;
%      save param_init.mat param_init;
     %=============================
% %      %retrieve initial moment files for fixed-point iteration in BLE
%==============================================
% file_name=['ble_init_auxiliary_moments' num2str(jj) '.mat'];
% location=[cd,'\auxiliary_files\'];
% string=[location file_name];
% load(string);
%      model.beta_init_FP=0*eye(model.numVar);
% model.beta_init_FP(model.FL_indices,model.FL_indices)=...
%    0.99*diag(second_moments_init);


%%
%initialize belief files for learning models 

if model.learning == 1 
    
    switch model.initial_beliefs 
        case 'ar(1),ble-based,fixed'
            
            file_name=['auxiliary_files/auxiliary_forecast_model_ble' num2str(jj) '.mat']
            db_tmp = load(file_name);
            alpha_init=db_tmp.alpha_tt(model.FL_indices);
            beta_init= diag(db_tmp.beta_tt(model.FL_indices,model.FL_indices));
            rr_init = (db_tmp.rr_tt(:,:,model.FL_indices));
            save AR1_BLE_initial_beliefs.mat beta_init rr_init;
            
      case 'msv,ree-based,fixed'
            
%             file_name=['auxiliary_files/auxiliary_forecast_model_ree' num2str(jj) '.mat']
%             db_tmp = load(file_name);
%             [initb]=beliefs_initialization_ree(db_tmp);
%             
%             rr_init=initb.rr_tt ;
%             alpha_init=initb.alpha_tt;
%             beta_init=initb.beta_tt;
%             dd_init=initb.dd_tt;
%             
%         ind_ = [1 2 3 4 5 6 7 8 9 10 11 13 14 15];    
%     rr_init=rr_init(ind_,ind_);
%     beta_init=beta_init(model.FL_indices,model.BL_indices);
%      dd_init=dd_init(model.FL_indices,:);
%      dd_init(4,:) = zeros(7,1);
     
%        save('dynare_initial_beliefs/MSV_initial_beliefs.mat');
     
                 case 'var(1),ree-based,fixed'
%                  
%             file_name=['auxiliary_files/auxiliary_forecast_model_ree' num2str(jj) '.mat']
%             db_tmp = load(file_name);
%             [initb]=beliefs_initialization_ree(db_tmp);
%             
%             rr_init=initb.rr_tt ;
%             alpha_init=initb.alpha_tt;
%             beta_init=initb.beta_tt;
%             dd_init=initb.dd_tt;
%             
%         ind_ = [1 2 3 4 5 6 7 8 ];    
%     rr_init=rr_init(ind_,ind_);
%     beta_init=beta_init(model.FL_indices,model.BL_indices);
%      dd_init=0*dd_init(model.FL_indices,:);
%             
            
%  save 'dynare_initial_beliefs/var1_initial_beliefs.mat' rr_init beta_init dd_init;
    
    end
    
end

model.forecast_index= jj;

%%
%---------------------------------------------------------------------------
last_obs=forecast.first_forecast-1 + (jj-1);%estimate up to first variable you want to forecast -1
first_obs=last_obs-forecast.window_length+1;%keeping sample size fixed==>rolling window

%load('full_dataset.mat');first_obs=106;last_obs=length(dy);
model.dataset=[dy dc dinve dw cpi_quarterly robs labobs];
model.dataset=model.dataset(first_obs:last_obs,:);
model.N=length(model.dataset);%dataset length
model.forecast_index=jj;
%set the objective function

load param_init; 
param_init=reshape(param_init,[1,length(param_init)]);

objective=@(x) likelihood(x,model);

if jj==1
options=optimset('Display','iter','MaxIter',500,'UseParallel',true);
else
 options=optimset('Display','iter','MaxIter',100,'UseParallel',true);
end
options_ps =optimoptions('particleswarm','SwarmSize',500,...
    'InitialSwarmMatrix',repmat(param_init,[500 1]),'UseParallel',true);
% options=optimset('Display','iter','MaxIter',20000,'UseParallel',true);

%optimization initial values

%init_H=nhess_diagonal(@likelihood,param_init');
%init_H=inv(init_H);init_H=diag(diag(init_H)); 
%save init_H;

%---------------------------------------------------------------------------

%---------------------------------------------------------------------------
%optimization input values for model

switch model.compute_init_hess
    case 2
        init_H=nhess(objective,param_init');
        init_H=inv(init_H);
        save init_H.mat init_H;
    case 1
        init_H=nhess_diagonal(objective,param_init');
        init_H=inv(init_H);init_H=diag(diag(init_H)); 
        save init_H.mat init_H;
    case 0
        load init_H.mat;
end

save options_initial.mat model param_init init_H;



switch model.optimizer
    case 1
[x fh]=fminsearch(objective,param_init',options);
    case 2
[fh,x,gh,H,itct,fcount,retcodeh] = csminwel(objective,param_init,init_H,[] ,10^(-6),9999);
    case 3
[x,fh,EXITFLAG,OUTPUT,LAMBDA,GRAD,H] =fmincon(objective,param_init,[],[],[],[],[],[],[],options);
    case 4
[x,fh,EXITFLAG,OUTPUT]=patternsearch(objective,param_init,[],[],[],[],[],[],[],options);     
    case 5
     x=   particleswarm(objective,length(param_init),model.LB,model.UB,options_ps);
end

x=reshape(x,[length(x),1]);
param_init=x;save param_init.mat param_init;
    if model.optimizer == 2
st_dev=sqrt(diag(H));
init_H=H;save init_H.mat init_H;
% laplace1=laplace_approximator(fh,x(model.est_indices),H(model.est_indices,model.est_indices))
laplace1=laplace_approximator(fh,x,H);
results=table(names',x,st_dev)
disp(laplace1);
    else 
results=table(names',x)   
    end
save estimation_results.mat;
%========================================================
%========================================================
kalman_retrieve_output;
%h-step ahead forecasts for current period
[S_fore,obs_fore]=point_forecast(S,gamma1,gamma2,gamma3,sysmat,model,forecast);
%load previous forecasts and add this on top
load forecast_output;
S_fore_all(:,:,jj)=S_fore;
obs_fore_all(:,:,jj)=obs_fore;
% likl_all(jj)=laplace1;
mode_all(jj) = fh;
estim_data_length(jj)= size(model.dataset,1);
save forecast_output.mat S_fore_all obs_fore_all forecast mode_all estim_data_length;

file_name=['auxiliary_forecast_model_ar2' num2str(jj) '.mat']
location=[cd,'/auxiliary_files/']
string=[location file_name]

save(string);
 end
 
 toc