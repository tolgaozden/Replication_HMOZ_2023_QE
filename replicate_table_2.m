clear;
clc;

warning('off','all');


% inputs: 
%.mat files containing forecasts for each model. The results reported in
%the paper are based on files stored in
%posterior_mode/estimations_baseline/forecast_summary/forecast_output_[model_name].mat


%outputs: 

%the model comparison of RMSE's are exported into a spreadsheet called
%forecasts.xlsx 



cd('posterior_mode/estimations_baseline');

run('forecast_tables.m');
