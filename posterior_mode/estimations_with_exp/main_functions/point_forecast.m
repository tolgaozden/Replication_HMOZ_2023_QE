function [S_fore,obs_fore]=point_forecast(S,gamma1,gamma2,gamma3,sysmat,model,forecast)
S_end=S(end,:)';%retrieve state variables in the last period as a vector

S_fore=nan(length(S_end),forecast.horizon);
obs_fore=nan(size(model.dataset,2),forecast.horizon);

%first forecast
S_fore(:,1)=gamma1*S_end+gamma2;

if strcmp(model.PLM,'ar(2),t-1,rls')==1
F2=[sysmat.F,zeros(model.l,model.numVar)];
end


% % if strcmp(model.PLM,'ar(2),t-1,rls')==1
% %         for ii=2:forecast.horizon
% %         S_fore(:,ii)=gamma1*S_fore(:,ii-1)+gamma2;
% %         obs_fore(:,ii)=sysmat.E+F2*S_fore(:,ii);
% %         end
% % else
obs_fore(:,1)=sysmat.E + sysmat.F*S_fore(:,1);

        for ii=2:forecast.horizon
        S_fore(:,ii)=gamma1*S_fore(:,ii-1)+gamma2;
        obs_fore(:,ii)=sysmat.E+sysmat.F*S_fore(:,ii);
        end

% % end
end