function [output]=forecast_evaluation(output_name)
%obs_fore_all structure: 
% #numObs x #Horizon x #Forecast
load(output_name);
output=[];
output.FE=zeros(forecast.num_periods,forecast.numObs,forecast.horizon);
output.FEsqrt=zeros(forecast.num_periods,forecast.numObs,forecast.horizon);
%loop for horizon here%%to be written

for jj=1:forecast.horizon
actual=forecast.dataset(jj:end,:);
% fore=reshape(obs_fore_all(:,horizon,:),[forecast.num_periods,forecast.numObs]);
fore=squeeze(obs_fore_all(:,jj,1:end+1-jj))';

% % RMSE
output.RMSE(:,jj)=sqrt(mean((actual-fore).^2));
% % MSE 
output.MSE(:,jj)=(mean((actual-fore).^2));
% % MAE
output.MAE(:,jj)=mean(abs(actual-fore));
% % squared forecast error matrix
output.FE(1:end+1-jj,:,jj)=(actual-fore);
output.FEsqrt(1:end+1-jj,:,jj)=(actual-fore).^2;

%  output.uncentered_log_det(:,jj)= log(det(cov(output.FE(1:end+1-jj,:,jj))));

output.uncentered_log_det(:,jj)= log(det(   (output.FE(1:end+1-jj,:,jj))' * (output.FE(1:end+1-jj,:,jj)) ));
end
  %overall_alternative(j)=100*(varCovarREE(j)-varCovarBLE(j))/14;
output.likl=likl_all;
output.actual_data=actual;
output.fore=fore;
end