

clear;clc;%close all;
load forecast_output_ble.mat;%retrieve forecast related variables-->should be the same in all files
clearvars -except forecast;%
%dataset order: [dy dc dinve dw pinfobs robs labobs];
ree_result=forecast_evaluation('forecast_output_ree.mat');
msv_result=forecast_evaluation('forecast_output_msv.mat');
ble_result=forecast_evaluation('forecast_output_ble.mat');
sac_result=forecast_evaluation('forecast_output_sac.mat');


  output_file='forecast_template.csv';
 output_sheet='forecast_template';

names=[{'REE','MSV' ,'BLE' ,'SAC' } ] ;  
horizons=[1 2 4 8 12]';
variables=[{' horizon' 'dy' 'dc' 'dinve' 'dw' 'pinfobs' 'robs' 'labobs' 'overall'}];

%percentage gains over ree
msv_perc_gains=100*(ree_result.RMSE(:,horizons)'-msv_result.RMSE(:,horizons)')./ree_result.RMSE(:,horizons)';
sac_perc_gains=100*(ree_result.RMSE(:,horizons)'-sac_result.RMSE(:,horizons)')./ree_result.RMSE(:,horizons)';
ble_perc_gains=100*(ree_result.RMSE(:,horizons)'-ble_result.RMSE(:,horizons)')./ree_result.RMSE(:,horizons)';
msv_overall_perc_gains=100*(ree_result.uncentered_log_det(horizons)-msv_result.uncentered_log_det(horizons))/14;
sac_overall_perc_gains=100*(ree_result.uncentered_log_det(horizons)-sac_result.uncentered_log_det(horizons))/14;
ble_overall_perc_gains=100*(ree_result.uncentered_log_det(horizons)-ble_result.uncentered_log_det(horizons))/14;


%REE table
disp('REE')
disp('RMSE')
disp(variables);
disp('============================================================')
disp([horizons,ree_result.RMSE(:,horizons)',ree_result.uncentered_log_det(:,horizons)']);
%===================================================================
%MSV table
disp('MSV')
disp('RMSE')
disp(variables);
disp('============================================================')
disp([horizons,msv_result.RMSE(:,horizons)',msv_result.uncentered_log_det(:,horizons)']);
%===================================================================
%BLE table
disp('BLE')
disp('RMSE')
disp(variables);
disp('============================================================')
disp([horizons,ble_result.RMSE(:,horizons)',ble_result.uncentered_log_det(:,horizons)']);
%===================================================================
%SAC table
disp('SAC')
disp('RMSE')
disp(variables);
disp('============================================================')
disp([horizons,sac_result.RMSE(:,horizons)',sac_result.uncentered_log_det(:,horizons)']);
%=====================================================================




% [sac_result.RMSE;ree_result.RMSE;ble_result.RMSE;msv_result.RMSE]
% %=======================================
% disp('MSE-SAC,MSE-REE,MSE-BLE,MSE-MSV');
% [sac_result.MSE;ree_result.MSE;ble_result.MSE;msv_result.MSE]
% %=======================================
% disp('MAE-SAC,MAE-REE,MAE-BLE,MAE-MSV');
% [sac_result.MAE;ree_result.MAE;ble_result.MAE;msv_result.MAE]
% %=======================================
% disp('LOGDET-SAC,LOGDET-REE,LOGDET-BLE');
% [sac_result.uncentered_log_det;ree_result.uncentered_log_det;...
%     ble_result.uncentered_log_det;msv_result.uncentered_log_det]
% %=======================================
figure('Name','rolling window likelihood','units','normalized','outerposition',[0 0 1 0.5]);
plot(sac_result.likl,'lineWidth',3);
hold on;
plot(ree_result.likl,'lineWidth',3);
hold on;
plot(ble_result.likl,'lineWidth',3);
hold on;
plot(msv_result.likl,'lineWidth',3);
legend('sac','ree','ble','msv');
%====================================
  fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print(fig,'rolling_window_likelihoods','-dpdf'); 
% 
% %===================================
% %forecast combinations
load forecast_output_ble.mat;
forecast_ble_point=obs_fore_all;
load forecast_output_sac.mat;
forecast_sac_point=obs_fore_all;
load forecast_output_ree.mat;
forecast_ree_point=obs_fore_all;
load forecast_output_msv.mat;
forecast_msv_point=obs_fore_all;
% 
% obs_fore_all=(forecast_ble_point+forecast_sac_point)/2;
% likl_all=0;
% save forecast_comb1.mat obs_fore_all forecast likl_all;
% 
% obs_fore_all=(forecast_ble_point+forecast_ree_point+forecast_sac_point)/3;
% save forecast_comb2.mat obs_fore_all forecast likl_all;
% 
 obs_fore_all=(forecast_ble_point+forecast_sac_point+forecast_ree_point+forecast_msv_point)/4;
%obs_fore_all=(forecast_ble_point+forecast_sac_point+forecast_msv_point)/3;
save forecast_comb3.mat obs_fore_all forecast likl_all;
% 
% forecast_comb1_result=forecast_evaluation('forecast_comb1.mat');
% forecast_comb2_result=forecast_evaluation('forecast_comb2.mat');
forecast_comb3_result=forecast_evaluation('forecast_comb3.mat');

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
%============================================
forecast_comb_perc_gains=100*(ree_result.RMSE(:,horizons)'-forecast_comb3_result.RMSE(:,horizons)')./ree_result.RMSE(:,horizons)';
forecast_comb_overall_perc_gains=100*(ree_result.uncentered_log_det(horizons)-forecast_comb3_result.uncentered_log_det(horizons))/14;

disp('Forecast combination percentage gains')
disp(variables)
[horizons forecast_comb_perc_gains forecast_comb_overall_perc_gains']
%===================================================================
%SAC table
%====================================================================
%====================================================================
%====================================================================
%excel outputs
%====================================================================
output=round([ree_result.RMSE(:,horizons)',ree_result.uncentered_log_det(:,horizons)'],2);
 xlswrite(output_file,output,output_sheet,'b3');

 output=round([ble_result.RMSE(:,horizons)',ble_result.uncentered_log_det(:,horizons)'],2);
 xlswrite(output_file,output,output_sheet,'b11');
 
  output=round([msv_result.RMSE(:,horizons)',msv_result.uncentered_log_det(:,horizons)'],2);
 xlswrite(output_file,output,output_sheet,'b19');
 
  output=round([sac_result.RMSE(:,horizons)',sac_result.uncentered_log_det(:,horizons)'],2);
 xlswrite(output_file,output,output_sheet,'b27');
 
 output=round([ msv_perc_gains msv_overall_perc_gains'],2);
  xlswrite(output_file,output,output_sheet,'b35');
  
   output=round([ ble_perc_gains ble_overall_perc_gains'],2);
  xlswrite(output_file,output,output_sheet,'b43');
  
   output=round([ sac_perc_gains sac_overall_perc_gains'],2);
  xlswrite(output_file,output,output_sheet,'b51');
  
     output=round([ forecast_comb_perc_gains forecast_comb_overall_perc_gains'],2);
  xlswrite(output_file,output,output_sheet,'b59');