clear;clc;%close all;
restoredefaultpath;
addpath(genpath(cd));
names_endo=[{'mc' 'zcap' 'rk' 'k' 'q'  'Consumption' 'Investment' 'Output' 'lab' 'Inflation' 'w' 'r' 'kp'   'dy' 'dc' 'dinve' 'dw' 'eps_a'  'eps_b' 'eps_g' 'eps_i'  'eps_r'  'eps_p' 'eps_w' }];
 names_exo=[{ 'Productivity'  'Risk Premium' 'Government Spending' 'Investment'  'Monetary Policy'  'Price Mark-up' 'Wage Mark-up' }];
tic


load results/ar2_estimation_results.mat;
load('inputs/mcmc_summary_all_models.mat');
param_mcmc = mcmc_ar2_exp.mcmc_for_simul;
num_simul = size(param_mcmc,1);

for mcmc = 1:num_simul
% for mcmc=563:num_simul
disp(mcmc)
    
param=param_mcmc(mcmc,:);
[parameters,gain,sigma]=param_set(param,model);

eig_warn = zeros(model.N,1);

switch model.initial_beliefs 
    case {'ar(1),ble-based','msv,ree-based'}
    [A, B, C, D, Et ,RHO ,Ft, G ,E ,F]=SW_sysmat_nonreduced(parameters);
    
sysmat=[];sysmat.A=A;sysmat.B=B;sysmat.C=C;sysmat.D=D;sysmat.Et=Et;sysmat.RHO=RHO;sysmat.Ft=Ft;sysmat.G=G;sysmat.E=E;sysmat.F=F;sysmat.sigma=sigma;
[alpha_tt,beta_tt,dd_tt,ee_tt rr_tt,indicator]=initialize_beliefs_new(parameters,model,sysmat);
%     case 'msv,ree-based,fixed'
%     [A, B, C, D, Et ,RHO ,Ft, G ,E ,F]=SW_sysmat_nonreduced(parameters);
%     
% sysmat=[];sysmat.A=A;sysmat.B=B;sysmat.C=C;sysmat.D=D;sysmat.Et=Et;sysmat.RHO=RHO;sysmat.Ft=Ft;sysmat.G=G;sysmat.E=E;sysmat.F=F;sysmat.sigma=sigma;
% [alpha_tt,beta_tt,dd_tt,ee_tt rr_tt,indicator]=initialize_beliefs_new(parameters,model,sysmat);
    otherwise


[alpha_tt,beta_tt,dd_tt,ee_tt rr_tt,indicator]=initialize_beliefs(parameters,model);
end

bmat = [];
bmat.alpha = alpha_tt;
bmat.beta = beta_tt;
bmat.dd = dd_tt;
bmat.ee = ee_tt;
bmat.rr = rr_tt;


[A, B , C, D , E , F]=SW_sysmat_reducedform(parameters,model,bmat);
sysmat=[];
sysmat.A=A;sysmat.B=B;sysmat.C=C;sysmat.D=D;sysmat.E=E;sysmat.F=F;
sysmat.sigma=sigma;
[gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_tt,dd_tt,ee_tt);

% initial_belief_special_cases;


if strcmp(model.PLM,'ar(2),t-1,rls')==1
S(1,:)=zeros(2*model.numVar,1);P(:,:,1)=10*eye(2*model.numVar);%use companion form if AR(2) rule is used
S(2,:)=zeros(2*model.numVar,1);P(:,:,2)=10*eye(2*model.numVar);%use companion form if AR(2) rule is used
start_index=3;
    [fdim1 fdim2] = size(F);
F=[F,zeros(fdim1,fdim2)];

elseif strcmp(model.PLM,'ar(2),t,rls')==1
S(1,:)=zeros(model.numVar,1);P(:,:,1)=10*eye(model.numVar);%use companion form if AR(2) rule is used
S(2,:)=zeros(model.numVar,1);P(:,:,2)=10*eye(model.numVar);%use companion form if AR(2) rule is used
start_index=3;


else
S(1,:)=zeros(model.numVar,1);P(:,:,1)=10*eye(model.numVar);
likl=nan(model.N,1);
pr_flag=zeros(model.N,1);
start_index=2;
end

likl=nan(model.N,1);

pr_flag=zeros(model.N,1);



try

for i=start_index:model.N
    
    model.tt= i;
    
    if model.decreasing_gain==1
        gain=1/i;
    end
    
    
    
%     if model.tt<model.exp_start
%     data_i=model.dataset(i,1:end-1)';
%     model.l=length(data_i);
%     else
%     data_i=model.dataset(i,:)';
%     model.l=length(data_i);
%     end


data_i = model.dataset(i,:)';
if isnan(data_i(end))==1
    data_i = data_i(1:end-1);
end

model.l = length(data_i);
   
 [E, F] = update_measurement( parameters,model);
 sysmat.E = E;
 sysmat.F = F;
 
 if strcmp(model.PLM,'ar(2),t-1,rls')==1
    [fdim1 fdim2] = size(F);
F=[F,zeros(fdim1,fdim2)];
sysmat.F = F;

end
   
    Sprime=gamma1*S(i-1,:)'+gamma2;
    Pprime=gamma1*P(:,:,i-1)*gamma1'+gamma3*sigma*gamma3';
     kGain=Pprime*F'*(F*Pprime*F')^(-1);
    S(i,:)=Sprime+kGain*(data_i-E-F*Sprime);
    P(:,:,i)=Pprime-kGain*(F*Pprime);
    v=data_i-E-F*Sprime;
       Fe=(F*Pprime*F');
likl(i)=-0.5*model.l*log(2*pi)-0.5*log(det(Fe))-0.5*v'*((Fe)\v);
      

 switch model.learning
     case 1
 [alpha_tt beta_tt dd_tt ee_tt rr_tt largest_eig_plm(i) pr_flag(i)...
    alpha_old beta_old dd_old ee_old rr_old]=update_beliefs(parameters,model,sysmat,...
    alpha_tt,beta_tt,dd_tt,ee_tt,rr_tt,S,i,gain);
% 

try
bmat = [];bmat.alpha = alpha_tt;bmat.beta = beta_tt;bmat.dd = dd_tt;bmat.ee = ee_tt;bmat.rr = rr_tt;
[A, B , C, D , E , F]=SW_sysmat_reducedform(parameters,model,bmat);
sysmat=[];sysmat.A=A;sysmat.B=B;sysmat.C=C;sysmat.D=D;sysmat.E=E;sysmat.F=F;sysmat.sigma=sigma;
[gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_tt,dd_tt,ee_tt);

                    switch model.projection_facility
                    case 1
                 [pr_flag(i),largest_eig(i),alpha_tt,beta_tt,dd_tt,ee_tt,rr_tt]=...
                projection_facility(gamma1,alpha_old,beta_old,dd_old,ee_old,rr_old,alpha_tt,beta_tt,dd_tt,ee_tt,rr_tt,pr_flag(i),model);
            
            
bmat = [];bmat.alpha = alpha_tt;bmat.beta = beta_tt;bmat.dd = dd_tt;bmat.ee = ee_tt;bmat.rr = rr_tt;
[A, B , C, D , E , F]=SW_sysmat_reducedform(parameters,model,bmat);
sysmat=[];sysmat.A=A;sysmat.B=B;sysmat.C=C;sysmat.D=D;sysmat.E=E;sysmat.F=F;sysmat.sigma=sigma;

if strcmp(model.PLM,'ar(2),t-1,rls')==1
    [fdim1 fdim2] = size(F);
F=[F,zeros(fdim1,fdim2)];
sysmat.F = F; 
end

                [gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_tt,dd_tt,ee_tt);
                
                
                
                        otherwise 
                           
                     largest_eig(i)=eigs(gamma1,1);       
                            
                     end
catch
    largest_eig(i) = 999;
    eig_warn(i) = 1;
    likl(i)=-999;
    
end
    
bmat = [];bmat.alpha = alpha_tt;bmat.beta = beta_tt;bmat.dd = dd_tt;bmat.ee = ee_tt;bmat.rr = rr_tt;
[A, B , C, D , E , F]=SW_sysmat_reducedform(parameters,model,bmat);
sysmat=[];sysmat.A=A;sysmat.B=B;sysmat.C=C;sysmat.D=D;sysmat.E=E;sysmat.F=F;sysmat.sigma=sigma;

if strcmp(model.PLM,'ar(2),t-1,rls')==1
    [fdim1 fdim2] = size(F);
F=[F,zeros(fdim1,fdim2)];
sysmat.F=F; 
end

[gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_tt,dd_tt,ee_tt);
 
     case 0
         

 [E, F] = update_measurement( parameters,model);
 sysmat.E = E; 
 sysmat.F = F;
 
 if strcmp(model.PLM,'ar(2),t-1,rls')==1
    [fdim1 fdim2] = size(F);
F=[F,zeros(fdim1,fdim2)];
sysmat.F = F;
end

 end
 
 
 
 
 
 if model.learning==1
     
     switch model.learning_algo
         case 'ar(1),rls'  
     
%  learning_parameters(i,:,:)=[alpha_tt(model.FL_indices)...
%    beta_tt(model.FL_indices,model.BL_indices) dd_tt(model.FL_indices,:) ee_tt(model.FL_indices,:)];
  learning_parameters(i,:,:,mcmc)=[alpha_tt(model.FL_indices) diag(beta_tt(model.FL_indices,model.FL_indices)) ];
   beta_all(:,i,mcmc) = diag(beta_tt(model.FL_indices,model.FL_indices));
  alpha_all(:,i,mcmc) = learning_parameters(i,:,1,mcmc);
   %  learning_covar(i,:,:,:)=rr_tt(:,:,model.FL_indices);
   
         case 'var(1)'  
     
  learning_parameters(i,:,:,mcmc)=[alpha_tt(model.FL_indices) (beta_tt(model.FL_indices,model.BL_indices)) ];
   beta_all(:,:,i,mcmc) = (beta_tt(model.FL_indices,model.BL_indices));
  alpha_all(:,i,mcmc) = learning_parameters(i,:,1,mcmc);
   %  learning_covar(i,:,:,:)=rr_tt(:,:,model.FL_indices);

         case 'ar(1),sac'
             
                  
%  learning_parameters(i,:,:)=[alpha_tt(model.FL_indices)...
%    beta_tt(model.FL_indices,model.BL_indices) dd_tt(model.FL_indices,:) ee_tt(model.FL_indices,:)];
  learning_parameters(i,:,:,mcmc)=[alpha_tt(model.FL_indices) diag(beta_tt(model.FL_indices,model.FL_indices)) ];
   beta_all(:,i,mcmc) = diag(beta_tt(model.FL_indices,model.FL_indices));
   alpha_all(:,i,mcmc) = learning_parameters(i,:,1,mcmc);
%  learning_covar(i,:,:,:)=rr_tt(:,:,model.FL_indices);


         case 'ar(2),rls'
%    learning_parameters(i,:,:)=[alpha_tt(model.FL_indices) ... 
%      diag(beta_tt(model.FL_indices,model.FL_indices)) ... 
%        diag(beta_tt(model.FL_indices+model.numEndo,model.FL_indices+model.numEndo)) ];
   
   learning_parameters(i,:,:,mcmc) = [alpha_tt(model.FL_indices), beta_tt(1:7), beta_tt(8:14)];
%    beta_all(:,i) = diag(beta_tt(model.FL_indices,model.FL_indices));   
         case 'msv'
             learning_parameters(i,:,:,mcmc)  =[alpha_tt(model.FL_indices) ...
                 (beta_tt(model.FL_indices,model.BL_indices)) ... 
                 dd_tt(model.FL_indices,:)];
                 
             alpha_all(:,i,mcmc) = alpha_tt(model.FL_indices);
             beta_all(:,:,i,mcmc) = beta_tt(model.FL_indices,model.BL_indices);
             dd_all(:,:,i,mcmc) = dd_tt(model.FL_indices,:);
         case 'msv,restricted'
             learning_parameters(i,:,:,mcmc)  =[alpha_tt(model.FL_indices) ...
                 (beta_tt(model.FL_indices,model.BL_indices)) ... 
                 dd_tt(model.FL_indices,:)];
                 
             alpha_all(:,i,mcmc) = alpha_tt(model.FL_indices);
             beta_all(:,:,i,mcmc) = beta_tt(model.FL_indices,model.BL_indices);
%              dd_all(:,:,i) = dd_tt(model.FL_indices,:);

     end
 
end
 
 
end

catch
%     
%                learning_parameters(:,:,:,mcmc)  =learning_parameters(:,:,:,mcmc-1);
%                alpha_all(:,i,mcmc)=alpha_all(:,i,mcmc-1);
%                beta_all(:,:,i,mcmc)=beta_all(:,:,i,mcmc-1);
%                dd_all(:,:,i,mcmc)=dd_all(:,:,i,mcmc-1);
%                S=S_mcmc(:,:,mcmc-1);
end
 


if model.learning==0
    largest_eig_gamma1= abs(eigs(gamma1,1));
%     largest_eig_gamma2 = abs(eigs(gamma3,1));
end


%=====================
if indicator==0
 likl=-sum(likl(1+model.burnIn:model.N));
elseif indicator==1
  likl=-sum(likl(1+model.burnIn:model.N)) - model.penalty;
end




if model.learning == 0 && strcmp(model.initial_beliefs,'ar(1),ble-based')==1
beta_BLE(:,mcmc) = diag(beta_tt(model.FL_indices,model.FL_indices));

end

S_mcmc(:,:,mcmc) = S;
toc

end


% save BLE_posteriorDist.mat beta_BLE S_mcmc param_mcmc model sysmat ;
% save MSV_posteriorDist.mat learning_parameters alpha_all beta_all dd_all S_mcmc param_mcmc  model sysmat;
% save SAC_posteriorDist.mat learning_parameters alpha_all beta_all  S_mcmc param_mcmc  model sysmat;
save AR2_posteriorDist.mat learning_parameters S_mcmc param_mcmc model sysmat;

