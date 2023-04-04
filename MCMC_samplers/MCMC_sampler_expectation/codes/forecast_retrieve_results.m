clear;clc;close all;
tic
names=[{'phi','sigma_c' ,'lambda' ,'xi_w' ,'sigma_l', 'xi_p' , 'iota_w','iota_p',...
    'psi' ,'phi_p','r_pi', 'rho' ,'r_y', 'r_dy' ,...
           'pi_bar' ,'beta_const' ,'l_bar', 'gamma_bar' ,'alpha'...
           'rho_a', 'rho_b' ,'rho_g' ,'rho_i' ,'rho_r', 'rho_p', 'rho_w','rho_ga',...
           'eta_a', 'eta_b' ,'eta_g' ,'eta_i' ,'eta_r', 'eta_p', 'eta_w',...
           'gain'} ] ;  
%retrieve data
load('slobodyan_dataset.mat');
  %forecast-variables
  forecast=[];

 
  forecast.numVar=24;%same as model.numVar below
  forecast.numObs=7;
  forecast.first_forecast=151; %index of the first variable to forecast
  forecast.horizon=12; %how many steps ahead in each forecast--> h-step ahead gives [1,2...h]
  forecast.window_length=80;%sample size; rolling window length for the estimation period.
  forecast.num_periods=length(dy)-forecast.first_forecast+1;%how many periods do we want to forecast



   %construct forecast dataset: this is the 1-step ahead dataset to be used in forecast evaluation
  %extend to h-step ahead later
  first_obs=forecast.first_forecast;
  last_obs=forecast.first_forecast+forecast.num_periods-1;
 forecast.dataset=[dy dc dinve dw pinfobs robs labobs];
 forecast.dataset=forecast.dataset(first_obs:last_obs,:);

for jj=1:forecast.num_periods

file_name=['auxiliary_forecast_model_msv' num2str(jj) '.mat']
location=[cd,'\auxiliary_files\']
string=[location file_name]
load(string);

forecast.horizon=12;%loading the above overwrites the previous setup...
  S_fore_all= nan(forecast.numVar,forecast.horizon,forecast.num_periods);
  obs_fore_all=nan(forecast.numObs,forecast.horizon,forecast.num_periods);
  
save estimation_results.mat;
kalman_retrieve_output;
%h-step ahead forecasts for current period
[S_fore,obs_fore]=point_forecast(S,gamma1,gamma2,gamma3,sysmat,model,forecast);
%load previous forecasts and add this on top
load forecast_output;
S_fore_all(:,:,jj)=S_fore;
obs_fore_all(:,:,jj)=obs_fore;

%  laplace1=laplace_approximator(fh,x([1:5,7:24,26:34]),(H([1:5,7:24,26:34],[1:5,7:24,26:34])))
laplace1=laplace_approximator(fh,x,H);

likl_all(jj)=laplace1
estim_data_length(jj)= size(model.dataset,1);
save forecast_output.mat S_fore_all obs_fore_all forecast likl_all estim_data_length;

end