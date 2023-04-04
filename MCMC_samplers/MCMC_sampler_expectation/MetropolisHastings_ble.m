clear;clc;close all;
clear;clc;%close all;
restoredefaultpath;
addpath('helpers');
addpath('Matfiles');
addpath('figures');
addpath('dynare_initial_beliefs');
addpath('inputs');
addpath('results');
% seed=round(rand*100000)
seed=111;
rng(seed);
tic
load results/ble_estimation_results.mat;

delete MCMC_in_progress.log;
diary MCMC_in_progress.log;

%  model.decreasing_gain=0;
mode=x;
% Sigma=diag(diag(H));

%estimate initial hessian at mode
% init_H=nhess_diagonal(objective,mode);
% init_H = nearestSPD(init_H);
% init_H=inv(init_H);
% Sigma=init_H;

Sigma = 0.001*eye(36);


%============================
Ndraws_initialize=5000;
Ndraws=500000;
numVar=length(mode);
c=0.01; %
c_eps=0.2;
recursiveAverages=nan(Ndraws,numVar);
acc_lb=0.2;
acc_ub=0.4;
Nburn=round(Ndraws/2);
% Nburn=100000;
posteriorDraws=nan(Ndraws,numVar);
currentDraw=mvnrnd(mode,c*Sigma);
% posteriorDraws(1,:)=currentDraw;
posteriorDraws(1,:) = mode;
accept=0;
objective=-likelihood(posteriorDraws(1,:),model);
counter=0;
logposterior=nan*ones(Ndraws,1);


%first run a small chain to initialize the Covariance matrix 
posteriorDraws_initialize = nan(Ndraws_initialize,numVar);
logposterior_initialize = nan(Ndraws_initialize,1);
posteriorDraws_initialize(1,:) = mode;
for i=1:Ndraws_initialize
    
    currentDraw=mvnrnd(posteriorDraws_initialize(i,:),c*Sigma);
    objectiveNew=-likelihood(currentDraw,model);
    alpha=min(1,exp(objectiveNew-objective));
    u=rand(1);

    if u<=alpha
        posteriorDraws_initialize(i+1,:)=currentDraw;
        accept=accept+1;
        objective=objectiveNew;
        logposterior_initialize(i+1)=objectiveNew;
    else
        posteriorDraws_initialize(i+1,:)=posteriorDraws_initialize(i,:);
        logposterior_initialize(i+1)=objective;
    end

acceptanceRate=accept/i; 

counter=counter+1;
if counter==50
    toc
     disp(['Acceptance Rate: ', num2str(acceptanceRate)]);
    disp(['Remaining Draws: ', num2str(Ndraws_initialize-i)]);
     disp(['Scale Coefficient: ', num2str(c)]);
    counter=0;
end

    
end


Sigma = cov(posteriorDraws_initialize);
Sigma=nearestSPD(Sigma);
c=0.1;
accept=0;

for i=1:Ndraws
    
    currentDraw=mvnrnd(posteriorDraws(i,:),c*Sigma);
    objectiveNew=-likelihood(currentDraw,model);
    alpha=min(1,exp(objectiveNew-objective));
    u=rand(1);

    if u<=alpha
        posteriorDraws(i+1,:)=currentDraw;
        accept=accept+1;
        objective=objectiveNew;
        logposterior(i+1)=objectiveNew;
    else
        posteriorDraws(i+1,:)=posteriorDraws(i,:);
        logposterior(i+1)=objective;
    end

acceptanceRate=accept/i;

% if acceptanceRate<acc_lb
%     c=(1-c_eps)*c;
% elseif acceptanceRate>acc_ub
%         c=(1+c_eps)*c;
%  end

counter=counter+1;
if counter==500
    
    save mcmc_ble_220322.mat;
    
    toc
     disp(['Acceptance Rate: ', num2str(acceptanceRate)]);
    disp(['Remaining Draws: ', num2str(Ndraws-i)]);
     disp(['Scale Coefficient: ', num2str(c)]);
    counter=0;
end


end

save mcmc_ble_220322_final.mat;

% 
% 
% figure('Name','beta-level-sac','units','normalized','outerposition',[0 0 1 0.35]);
% for i=1:numVar;
% subplot(7,5,i);hist(posteriorDraws(Nburn:end,i));hold on;plot([mode(i) mode(i)],[ylim],'r'); 
% end
% 
% 
% 
% 
% figure;
% title('Recursive Averages');
% 
%      disp('calculating recursive averages for cusum plots...');
% %      disp(i)
% for j=1:Ndraws
%     disp(j)
%     recursiveAverages(j,:)=mean(posteriorDraws(1:j,:));
%     
% end
% 
%  for i=1:numVar;
%    subplot(7,5,i);
%    plot(recursiveAverages(:,i),'lineWidth',2);
%   
%  end
% 
%  disp('POSTERIOR MEAN')
%   mh_conf_sig=0.95;
%  posteriorMean=mean(posteriorDraws(Nburn+1:end,:))
% posteriorDist=posteriorDraws(Nburn+1:end,:);
% posteriorDist=sort(posteriorDist);
% hpdDraws=round((1-mh_conf_sig)*length(posteriorDist));
% kk=zeros(hpdDraws,13);
% 
% 
% for paraInd=1:length(mode)
%     jj=size(posteriorDist,1)-hpdDraws-2;
% 
% for ii=1:hpdDraws
%     kk(ii,paraInd)=posteriorDist(jj,paraInd)-posteriorDist(ii,paraInd);
%     jj=jj+1;
% end
% [kmin,idx]=min(kk(:,paraInd));
% hpd_interval(paraInd,:)=[posteriorDist(idx,paraInd) posteriorDist(idx,paraInd)+kmin];
% post_deciles(paraInd,:)=posteriorDist([round(0.05*length(posteriorDist(:,paraInd)))...
%     round(0.2*length(posteriorDist(:,paraInd)))...
%     round(0.3*length(posteriorDist(:,paraInd)))...
%     round(0.4*length(postescp -r [marx folder location] tolo@edith-login1:[edith location]riorDist(:,paraInd)))...
%     round(0.5*length(posteriorDist(:,paraInd)))...
%     round(0.6*length(posteriorDist(:,paraInd)))...
%     round(0.7*length(posteriorDist(:,paraInd)))...
%     round(0.8*length(posteriorDist(:,paraInd)))...
%     round(0.95*length(posteriorDist(:,paraInd)))],paraInd);
% end
%  
%  %=========================================================================\
% %modified harmonic mean estimator
% %posteriorDraws-->likl values
