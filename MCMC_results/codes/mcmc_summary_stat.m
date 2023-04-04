function [TT] = mcmc_summary_stat(file_name,input_draws)

load(file_name);

Ndraws=input_draws;
Nburn=round(Ndraws *0.5);

 disp('POSTERIOR MEAN')
  mh_conf_sig=0.95;
 posteriorMean=mean(posteriorDraws(Nburn+1:Ndraws,:))
posteriorDist=posteriorDraws(Nburn+1:Ndraws,:);
posteriorDist=sort(posteriorDist);
hpdDraws=round((1-mh_conf_sig)*length(posteriorDist));
kk=zeros(hpdDraws,13);


for paraInd=1:length(posteriorMean)
    jj=size(posteriorDist,1)-hpdDraws-2;

for ii=1:hpdDraws
    kk(ii,paraInd)=posteriorDist(jj,paraInd)-posteriorDist(ii,paraInd);
    jj=jj+1;
end
[kmin,idx]=min(kk(:,paraInd));
hpd_interval(paraInd,:)=[posteriorDist(idx,paraInd) posteriorDist(idx,paraInd)+kmin];
post_deciles(paraInd,:)=posteriorDist([round(0.05*length(posteriorDist(:,paraInd)))...
    round(0.2*length(posteriorDist(:,paraInd)))...
    round(0.3*length(posteriorDist(:,paraInd)))...
    round(0.4*length(posteriorDist(:,paraInd)))...
    round(0.5*length(posteriorDist(:,paraInd)))...
    round(0.6*length(posteriorDist(:,paraInd)))...
    round(0.7*length(posteriorDist(:,paraInd)))...
    round(0.8*length(posteriorDist(:,paraInd)))...
    round(0.95*length(posteriorDist(:,paraInd)))],paraInd);
end


TT= table(names',posteriorMean',hpd_interval);
TT.Properties.VariableNames={'Variable Names','Posterior Mean','90% HPD Interval'};
disp(TT)
