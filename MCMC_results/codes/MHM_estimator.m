function [MCMC]= MHM_estimator(file_name,input_draws,num_for_simul)

load(file_name);

%overwrite #draws and #burn_in

Ndraws=input_draws;
Nburn=round(Ndraws *0.5);




posterior_covar=zeros(numVar,numVar);
posterior_length=Ndraws-Nburn+1;


posteriorMean=mean(posteriorDraws(Nburn+1:Ndraws,:));



for ii=Nburn+1:Ndraws
    disp(ii)
    
XX=(posteriorDraws(ii,:)-posteriorMean)';    
    
 posterior_covar= posterior_covar + XX * XX';
    
end
    
posterior_covar=posterior_covar/posterior_length;   
detCovar=det(posterior_covar);
invCovar=inv(posterior_covar);
marginal=[];
deviation=zeros(posterior_length,1);

for pp=0.1:0.1:0.9
    index=0;
     critval = qchisq(pp,numVar);
     tmp=0;
     
     
    for ii=Nburn+1:Ndraws
     index=index+1;
     XX=(posteriorDraws(ii,:)-posteriorMean);   
     deviation(index)=XX*invCovar*XX';
        if deviation(index) <=critval
         lftheta= -log(pp)-(numVar*log(2*pi)+log(detCovar)+deviation(index))/2;
         tmp=tmp+exp(lftheta-logposterior(ii)-fh);
        end
    end

    
    
    marginal=cat(1,marginal,[pp,-fh-log(tmp/posterior_length)]);
    
end
format long;
disp(marginal(:,2));

MCMC.marginal = marginal;
if strcmp(file_name,'mcmc_sac.mat')==1
    MCMC.marginal = MCMC.marginal + 11.4026 %adjust for gain
end
    
MCMC.Ndraws = Ndraws;
MCMC.posterior_covar = posterior_covar;


%prepare thinned draws for simulations 

%draws after burn-in 

mcmc_exclude_burnin = posteriorDraws(Nburn+1:Ndraws,:); 

thin_factor = round(size(mcmc_exclude_burnin,1)/num_for_simul)


mcmc_for_simul = mcmc_exclude_burnin(1:thin_factor:end,:); 

mcmc_for_simul = mcmc_for_simul(1:min(num_for_simul,size(mcmc_for_simul,1)),:);

MCMC.mcmc_for_simul = mcmc_for_simul;

end
