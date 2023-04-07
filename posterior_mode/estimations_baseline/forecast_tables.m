

clear;clc;%close all;
addpath(genpath(cd));
load forecast_output_ble.mat;%retrieve forecast related variables-->should be the same in all files
% clearvars -except forecast;%
%dataset order: [dy dc dinve dw pinfobs robs labobs];
ree_result=forecast_evaluation('forecast_output_ree.mat');
msv_result=forecast_evaluation('forecast_output_msv.mat');
ble_result=forecast_evaluation('forecast_output_ble.mat');
sac_result=forecast_evaluation('forecast_output_sac.mat');
ar2_result=forecast_evaluation('forecast_output_ar2.mat');
var1_result=forecast_evaluation('forecast_output_var1.mat');



 output_file='forecast_template.csv';
 output_sheet='forecast_template';

names=[{'REE','MSV' ,'BLE' ,'SAC','VAR(1)','AR(2)'} ] ;  
horizons=[1 2 4 8 12]';
variables=[{' horizon' 'dy' 'dc' 'dinve' 'dw' 'cpi infl.' 'robs' 'labobs' 'overall'}];

%percentage gains over ree
msv_perc_gains=100*(ree_result.RMSE(:,horizons)'-msv_result.RMSE(:,horizons)')./ree_result.RMSE(:,horizons)';
sac_perc_gains=100*(ree_result.RMSE(:,horizons)'-sac_result.RMSE(:,horizons)')./ree_result.RMSE(:,horizons)';
ble_perc_gains=100*(ree_result.RMSE(:,horizons)'-ble_result.RMSE(:,horizons)')./ree_result.RMSE(:,horizons)';
var1_perc_gains=100*(ree_result.RMSE(:,horizons)'-var1_result.RMSE(:,horizons)')./ree_result.RMSE(:,horizons)';
ar2_perc_gains=100*(ree_result.RMSE(:,horizons)'-ar2_result.RMSE(:,horizons)')./ree_result.RMSE(:,horizons)';

msv_overall_perc_gains=100*(ree_result.uncentered_log_det(horizons)-msv_result.uncentered_log_det(horizons))/14;
sac_overall_perc_gains=100*(ree_result.uncentered_log_det(horizons)-sac_result.uncentered_log_det(horizons))/14;
ble_overall_perc_gains=100*(ree_result.uncentered_log_det(horizons)-ble_result.uncentered_log_det(horizons))/14;
var1_overall_perc_gains=100*(ree_result.uncentered_log_det(horizons)-var1_result.uncentered_log_det(horizons))/14;
ar2_overall_perc_gains=100*(ree_result.uncentered_log_det(horizons)-ar2_result.uncentered_log_det(horizons))/14;
% 
% %REE table
% disp('REE')
% disp('RMSE')
% disp(variables);
% disp('============================================================')
% disp([horizons,ree_result.RMSE(:,horizons)',ree_result.uncentered_log_det(:,horizons)']);
% %===================================================================
% %MSV table
% disp('MSV')
% disp('RMSE')
% disp(variables);
% disp('============================================================')
% disp([horizons,msv_result.RMSE(:,horizons)',msv_result.uncentered_log_det(:,horizons)']);
% %===================================================================
% %BLE table
% disp('BLE')
% disp('RMSE')
% disp(variables);
% disp('============================================================')
% disp([horizons,ble_result.RMSE(:,horizons)',ble_result.uncentered_log_det(:,horizons)']);
% %===================================================================
% %SAC table
% disp('SAC')
% disp('RMSE')
% disp(variables);
% disp('============================================================')
% disp([horizons,sac_result.RMSE(:,horizons)',sac_result.uncentered_log_det(:,horizons)']);
% %=====================================================================
% %VAR(1) table
% disp('VAR(1)')
% disp('RMSE')
% disp(variables);
% disp('============================================================')
% disp([horizons,var1_result.RMSE(:,horizons)',var1_result.uncentered_log_det(:,horizons)']);
% %=====================================================================
% %VAR(1) table
% disp('AR(2)')
% disp('RMSE')
% disp(variables);
% disp('============================================================')
% disp([horizons,ar2_result.RMSE(:,horizons)',ar2_result.uncentered_log_det(:,horizons)']);


load forecast_output_ble.mat;
forecast_ble_point=obs_fore_all;
load forecast_output_sac.mat;
forecast_sac_point=obs_fore_all;
load forecast_output_ree.mat;
forecast_ree_point=obs_fore_all;
load forecast_output_msv.mat;
forecast_msv_point=obs_fore_all;
load forecast_output_var1.mat;
forecast_var1_point=obs_fore_all;
load forecast_output_ar2.mat;
forecast_ar2_point=obs_fore_all;


%=======================================
disp('MSV percentage gains')
disp(variables)
[horizons msv_perc_gains msv_overall_perc_gains']
%=====================================
disp('SAC percentage gains')
disp(variables)
[horizons sac_perc_gains sac_overall_perc_gains']
%=========================================
disp('BLE percentage gains')
disp(variables)
[horizons ble_perc_gains ble_overall_perc_gains']
%=========================================
disp('VAR percentage gains')
disp(variables)
[horizons var1_perc_gains var1_overall_perc_gains']
%=========================================
disp('AR(2) percentage gains')
disp(variables)
[horizons ar2_perc_gains ar2_overall_perc_gains']

%============================================
% forecast_comb_perc_gains=100*(ree_result.RMSE(:,horizons)'-forecast_comb3_result.RMSE(:,horizons)')./ree_result.RMSE(:,horizons)';
% forecast_comb_overall_perc_gains=100*(ree_result.uncentered_log_det(horizons)-forecast_comb3_result.uncentered_log_det(horizons))/14;

% disp('Forecast combination percentage gains')
% disp(variables)
% [horizons forecast_comb_perc_gains forecast_comb_overall_perc_gains']
%===================================================================
%SAC table
%====================================================================
%====================================================================
%====================================================================
%excel outputs
%====================================================================
output=round([ree_result.RMSE(:,horizons)',ree_result.uncentered_log_det(:,horizons)'],2);
 writematrix(output,'forecasts.xlsx','sheet',output_sheet,'range','B3');
 

 output=round([ble_result.RMSE(:,horizons)',ble_result.uncentered_log_det(:,horizons)'],2);
 writematrix(output,'forecasts.xlsx','sheet',output_sheet,'range','B11');
 
  output=round([msv_result.RMSE(:,horizons)',msv_result.uncentered_log_det(:,horizons)'],2);
 writematrix(output,'forecasts.xlsx','sheet',output_sheet,'range','B19');
 
  output=round([sac_result.RMSE(:,horizons)',sac_result.uncentered_log_det(:,horizons)'],2);
writematrix(output,'forecasts.xlsx','sheet',output_sheet,'range','B27');
 
   output=round([var1_result.RMSE(:,horizons)',var1_result.uncentered_log_det(:,horizons)'],2);
 writematrix(output,'forecasts.xlsx','sheet',output_sheet,'range','B35');
 
    output=round([ar2_result.RMSE(:,horizons)',ar2_result.uncentered_log_det(:,horizons)'],2);
 writematrix(output,'forecasts.xlsx','sheet',output_sheet,'range','B43');
 
  output_sheet='forecast_template2';
 
 output=round([ msv_perc_gains msv_overall_perc_gains'],2);
 writematrix(output,'forecasts.xlsx','sheet',output_sheet,'range','B3');
  
   output=round([ ble_perc_gains ble_overall_perc_gains'],2);
 writematrix(output,'forecasts.xlsx','sheet',output_sheet,'range','B11');
  
   output=round([ sac_perc_gains sac_overall_perc_gains'],2);
writematrix(output,'forecasts.xlsx','sheet',output_sheet,'range','B19');
  
     output=round([ var1_perc_gains var1_overall_perc_gains'],2);
  writematrix(output,'forecasts.xlsx','sheet',output_sheet,'range','B27');
  
  output=round([ ar2_perc_gains ar2_overall_perc_gains'],2);
 writematrix(output,'forecasts.xlsx','sheet',output_sheet,'range','B35');
  
