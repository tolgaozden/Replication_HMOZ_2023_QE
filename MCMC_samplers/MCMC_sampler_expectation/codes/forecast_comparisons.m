
clear;clc;%close all;
load forecast_output_ble.mat;%retrieve forecast related variables-->should be the same in all files
clearvars -except forecast;%
%dataset order: [dy dc dinve dw pinfobs robs labobs];
ree_result=forecast_evaluation('forecast_output_ree.mat');
msv_result=forecast_evaluation('forecast_output_msv.mat');
ble_result=forecast_evaluation('forecast_output_ble.mat');
sac_result=forecast_evaluation('forecast_output_sac.mat');
ar2_result=forecast_evaluation('forecast_output_ar2_t_1_final.mat');


startDate=datenum('01-12-1984');
endDate = datenum('01-09-2008');
Date=linspace(startDate,endDate,length(forecast.dataset));




%first forecast model: 1965:I-1984:IV (71-150), first year
%presample. Forecast from 151 onwards. 
%... last model:       1988:IV-2008:III (166:245). Forecast 246 onw. 



names=[{'REE','MSV' ,'BLE' ,'SAC','AR(2)' } ] ;  
horizons=[1 2 4 8 12]';
variables=[{' horizon' 'dy' 'dc' 'dinve' 'dw' 'pinfobs' 'robs' 'labobs' 'overall'}];

%percentage gains over ree
msv_perc_gains=100*(ree_result.RMSE(:,horizons)'-msv_result.RMSE(:,horizons)')./ree_result.RMSE(:,horizons)';
sac_perc_gains=100*(ree_result.RMSE(:,horizons)'-sac_result.RMSE(:,horizons)')./ree_result.RMSE(:,horizons)';
ble_perc_gains=100*(ree_result.RMSE(:,horizons)'-ble_result.RMSE(:,horizons)')./ree_result.RMSE(:,horizons)';
ar2_perc_gains=100*(ree_result.RMSE(:,horizons)'-ar2_result.RMSE(:,horizons)')./ree_result.RMSE(:,horizons)';

msv_overall_perc_gains=100*(ree_result.uncentered_log_det(horizons)-msv_result.uncentered_log_det(horizons))/14;
sac_overall_perc_gains=100*(ree_result.uncentered_log_det(horizons)-sac_result.uncentered_log_det(horizons))/14;
ble_overall_perc_gains=100*(ree_result.uncentered_log_det(horizons)-ble_result.uncentered_log_det(horizons))/14;
ar2_overall_perc_gains=100*(ree_result.uncentered_log_det(horizons)-ar2_result.uncentered_log_det(horizons))/14;

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
%AR(2) table
disp('AR(2)')
disp('RMSE')
disp(variables);
disp('============================================================')
disp([horizons,ar2_result.RMSE(:,horizons)',ar2_result.uncentered_log_det(:,horizons)']);



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
figure('Name','rolling window likelihood','units','normalized','outerposition',[0 0 0.5 0.5]);

plot(Date,-sac_result.likl,'lineWidth',3,'LineStyle','-','color','red');
  % xlim([startDate endDate]);datetick('x','yyyy','keeplimits');
hold on;
plot(Date,-ree_result.likl,'lineWidth',3,'LineStyle','-.','color','black');

   % xlim([startDate endDate]);datetick('x','yyyy','keeplimits');
hold on;
plot(Date,-ble_result.likl,'lineWidth',3,'LineStyle','--','color','blue');
   % xlim([startDate endDate]);datetick('x','yyyy','keeplimits');
hold on;
plot(Date,-msv_result.likl,'lineWidth',3,'LineStyle',':','color','green');
   % xlim([startDate endDate]);datetick('x','yyyy','keeplimits');
hold on;
plot(Date,-ar2_result.likl,'lineWidth',3,'LineStyle','-.','color','yellow');
    xlim([startDate endDate]);datetick('x','yy','keeplimits');
    fig_=gca;
 fig_.XTickLabelRotation=90;
legend('SAC','REE','BLE','MSV','AR(2)');

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
load forecast_output_ar2_t_1_final.mat;
forecast_ar2_point=obs_fore_all;
% 
% obs_fore_all=(forecast_ble_point+forecast_sac_point)/2;
% likl_all=0;
% save forecast_comb1.mat obs_fore_all forecast likl_all;
% 
% obs_fore_all=(forecast_ble_point+forecast_ree_point+forecast_sac_point)/3;
% save forecast_comb2.mat obs_fore_all forecast likl_all;
% 
 obs_fore_all=(forecast_ble_point+forecast_sac_point+forecast_ree_point+forecast_msv_point+forecast_ar2_point)/5;
%obs_fore_all=(forecast_ble_point+forecast_sac_point+forecast_msv_point)/3;
save forecast_comb3.mat obs_fore_all forecast likl_all;

 obs_fore_all=(forecast_ble_point+forecast_sac_point+forecast_msv_point+forecast_ar2_point)/4;
%obs_fore_all=(forecast_ble_point+forecast_sac_point+forecast_msv_point)/3;
save forecast_comb2.mat obs_fore_all forecast likl_all;


% 
% forecast_comb1_result=forecast_evaluation('forecast_comb1.mat');
% forecast_comb2_result=forecast_evaluation('forecast_comb2.mat');
forecast_comb3_result=forecast_evaluation('forecast_comb3.mat');
forecast_comb2_result=forecast_evaluation('forecast_comb2.mat');

%forecast comb. table
disp('forecast comb-all, excluding ree model')
disp('RMSE')
disp(variables);
disp('============================================================')
disp([horizons,forecast_comb2_result.RMSE(:,horizons)',forecast_comb2_result.uncentered_log_det(:,horizons)']);

%===================================================================
%SAC table
disp('forecast comb-all')
disp('RMSE')
disp(variables);
disp('============================================================')
disp([horizons,forecast_comb3_result.RMSE(:,horizons)',forecast_comb3_result.uncentered_log_det(:,horizons)']);


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
disp('AR(2) percentage gains')
disp(variables)
[horizons ar2_perc_gains ar2_overall_perc_gains']
%============================================
forecast_comb3_perc_gains=100*(ree_result.RMSE(:,horizons)'-forecast_comb3_result.RMSE(:,horizons)')./ree_result.RMSE(:,horizons)';
forecast_comb3_overall_perc_gains=100*(ree_result.uncentered_log_det(horizons)-forecast_comb3_result.uncentered_log_det(horizons))/14;

forecast_comb2_perc_gains=100*(ree_result.RMSE(:,horizons)'-forecast_comb3_result.RMSE(:,horizons)')./ree_result.RMSE(:,horizons)';
forecast_comb2_overall_perc_gains=100*(ree_result.uncentered_log_det(horizons)-forecast_comb2_result.uncentered_log_det(horizons))/14;


disp('Forecast combination percentage gains (all included)')
disp(variables)
[horizons forecast_comb3_perc_gains forecast_comb3_overall_perc_gains']

disp('Forecast combination percentage gains (ree excluded)')
disp(variables)
[horizons forecast_comb2_perc_gains forecast_comb2_overall_perc_gains']

%bayes factors
%msv
msv_bayes=log10(exp(-msv_result.likl+ree_result.likl));
sac_bayes=log10(exp(-sac_result.likl+ree_result.likl));
ble_bayes=log10(exp(-ble_result.likl+ree_result.likl));
ar2_bayes=log10(exp(-ar2_result.likl+ree_result.likl));



%all pairwise combinations
%relative to ree
bayes_all(:,1)=log10(exp(-msv_result.likl+ree_result.likl));
bayes_all(:,2)=log10(exp(-sac_result.likl+ree_result.likl));
bayes_all(:,3)=log10(exp(-ble_result.likl+ree_result.likl));
bayes_all(:,4)=log10(exp(-ar2_result.likl+ree_result.likl));
%relative to msv
bayes_all(:,5)=log10(exp(-sac_result.likl+msv_result.likl));
bayes_all(:,6)=log10(exp(-ble_result.likl+msv_result.likl));
bayes_all(:,7)=log10(exp(-ar2_result.likl+msv_result.likl));
%relative to sac
bayes_all(:,8)=log10(exp(-ble_result.likl+sac_result.likl));
bayes_all(:,9)=log10(exp(-ar2_result.likl+sac_result.likl));
%relative to ble
bayes_all(:,10)=log10(exp(-ar2_result.likl+ble_result.likl));

comp_names = {'MSV vs. REE','SAC vs. REE','BLE vs. REE','AR(2) vs. REE','SAC vs. MSV','BLE vs. MSV','AR(2) vs. MSV','BLE vs. SAC','AR(2) vs. SAC','AR(2) vs. BLE'}



% msv_bayes=(exp(-msv_result.likl+ree_result.likl));
% sac_bayes=(exp(-sac_result.likl+ree_result.likl));
%  ble_bayes=(exp(-ble_result.likl+ree_result.likl));



figure('Name','Bayes Factors','units','normalized','outerposition',[0 0 0.5 0.5]);
subplot(2,2,1);
plot(Date,msv_bayes,'lineWidth',3,'color','black');
hold on;
plot(Date,zeros(length(ble_bayes),1),'--','color','black');

hold on
plot(Date,2*ones(length(ble_bayes),1),'--','color','black');
hold on
plot(Date,-2*ones(length(ble_bayes),1),'--','color','black');


xlim([0 length(msv_bayes)]);
title('MSV');
    xlim([startDate endDate]);datetick('x','yyyy','keeplimits');
  fig_=gca;
 fig_.XTickLabelRotation=90;
 ylim([-20 20]);
 
subplot(2,2,2);
plot(Date,sac_bayes,'lineWidth',3,'color','black');
hold on;
plot(Date,0*ones(length(ble_bayes),1),'--','color','black');

hold on
plot(Date,2*ones(length(ble_bayes),1),'--','color','black');
hold on
plot(Date,-2*ones(length(ble_bayes),1),'--','color','black');

xlim([0 length(msv_bayes)]);
title('SAC');
    xlim([startDate endDate]);datetick('x','yyyy','keeplimits');
 fig_=gca;
 fig_.XTickLabelRotation=90;
 ylim([-20 20]);
subplot(2,2,3);
plot(Date,ble_bayes,'lineWidth',3,'color','black');
hold on;
plot(Date,0*ones(length(ble_bayes),1),'--','color','black');

hold on
plot(Date,2*ones(length(ble_bayes),1),'--','color','black');
hold on
plot(Date,-2*ones(length(ble_bayes),1),'--','color','black');

xlim([0 length(msv_bayes)]);
title('BLE');
    xlim([startDate endDate]);datetick('x','yyyy','keeplimits');
   fig_=gca;
 fig_.XTickLabelRotation=90;
 ylim([-20 20]);
subplot(2,2,4);
plot(Date,ar2_bayes,'lineWidth',3,'color','black');
%legend('msv','sac','ble','ar(2)');
hold on;
plot(Date,0*ones(length(ble_bayes),1),'--','color','black');

hold on
plot(Date,2*ones(length(ble_bayes),1),'--','color','black');
hold on
plot(Date,-2*ones(length(ble_bayes),1),'--','color','black');

title('AR(2)');
xlim([0 length(msv_bayes)]);
    xlim([startDate endDate]);datetick('x','yyyy','keeplimits');
 fig_=gca;
 fig_.XTickLabelRotation=90;
ylim([-20 20]);

  fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print(fig,'bayes_factor_comparisons','-dpdf'); 



figure('Name','Bayes Factors','units','normalized','outerposition',[0 0 0.5 1]);
for jj=1:10
    subplot(5,2,jj);
plot(Date,bayes_all(:,jj),'lineWidth',3,'color','black');
hold on;
plot(Date,zeros(length(ble_bayes),1),'--','color','black');
hold on
plot(Date,2*ones(length(ble_bayes),1),'--','color','black');
hold on
plot(Date,-2*ones(length(ble_bayes),1),'--','color','black');

xlim([0 length(ble_bayes)]);
title(comp_names(jj));
    xlim([startDate endDate]);datetick('x','yy','keeplimits');
  fig_=gca;
 fig_.XTickLabelRotation=90;
 ylim([-20 20]);
 
end
 

  fig = gcf;
fig.PaperPositionMode = 'auto'
fig_pos = fig.PaperPosition;
fig.PaperSize = [fig_pos(3) fig_pos(4)];

print(fig,'bayes_factor_comparisons_all','-dpdf'); 

