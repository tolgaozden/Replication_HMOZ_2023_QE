clear;clc;%close all;

addpath('helpers');
addpath('Matfiles');
addpath('figures');
addpath('dynare_initial_beliefs');
addpath('inputs');
addpath('results');

%% 
names_endo=[{'mc' 'zcap' 'rk' 'k' 'q'  'Consumption' 'Investment' 'Output' 'lab' 'Inflation' 'w' 'r' 'kp'   'dy' 'dc' 'dinve' 'dw' 'eps_a'  'eps_b' 'eps_g' 'eps_i'  'eps_r'  'eps_p' 'eps_w' }];
 names_exo=[{ 'Productivity'  'Risk Premium' 'Government Spending' 'Investment'  'Monetary Policy'  'Price Mark-up' 'Wage Mark-up' }];

string_ ='results/sac_estimation_results.mat';
load(string_);
% 
% model.PLM='ar(2),t-1,rls';%options: 'msv,t','msv,t-1',ar(1),t,rls','ar(1),t,sac','ar(1),t-1,rls','ar(1),t-1,sac','ar(2),t-1,rls'
% model.learning=0; %options: 1=adaptive learning invoked, 0=beliefs fixed at initial values 
% model.learning_algo='ar(2)'; %options: 'msv','ar(1),rls','ar(1),sac','msv,restricted';
% model.initial_beliefs='ar(2),ree-based,fixed';%options: 'msv,ree-based,fixed','ar(1),ree-based,fixed','ar(1),ree-based','ar(1),ble-based','ar(1),ble-based,fixed'
% model.PLM_timing='t-1';%this is repetitive, same info in model.PLM, fix it later
% 

param=x;
[parameters,gain,sigma]=param_set(param,model);

eig_warn = zeros(model.N,1);


[A, B, C, D, Et ,RHO ,Ft, G ,E ,F]=SW_sysmat_MSV_filter(parameters);

sysmat=[];
sysmat.A=A;sysmat.B=B;sysmat.C=C;sysmat.D=D;sysmat.Et=Et;
sysmat.RHO=RHO;sysmat.Ft=Ft;sysmat.G=G;sysmat.E=E;sysmat.F=F;
sysmat.sigma=sigma;



[alpha_tt,beta_tt,dd_tt,ee_tt rr_tt,indicator]=initialize_beliefs(parameters,model,sysmat);
[gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_tt,dd_tt,ee_tt);

initial_belief_special_cases;


if strcmp(model.PLM,'ar(2),t-1,rls')==1
S(1,:)=zeros(2*model.numVar,1);P(:,:,1)=10*eye(2*model.numVar);%use companion form if AR(2) rule is used
S(2,:)=zeros(2*model.numVar,1);P(:,:,2)=10*eye(2*model.numVar);%use companion form if AR(2) rule is used
start_index=3;
F=[F,zeros(model.l,model.numVar)];

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





for i=start_index:model.N
    
    if model.decreasing_gain==1
        gain=1/i;
    end
    
    data_i=model.dataset(i,:)';
   
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
[gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_tt,dd_tt,ee_tt);
                    switch model.projection_facility
                    case 1
                 [pr_flag(i),largest_eig(i),alpha_tt,beta_tt,dd_tt,ee_tt,rr_tt]=...
                projection_facility(gamma1,alpha_old,beta_old,dd_old,ee_old,rr_old,alpha_tt,beta_tt,dd_tt,ee_tt,rr_tt,pr_flag(i),model);
                [gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_tt,dd_tt,ee_tt);
                        otherwise 
                           
                     largest_eig(i)=eigs(gamma1,1);       
                            
                     end
catch
    largest_eig(i) = 999;
    eig_warn(i) = 1;
    likl(i)=-999;
    
end
    

[gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_tt,dd_tt,ee_tt);
 end
 
 
 if model.learning==1
     
     switch model.learning_algo
         case 'ar(1),rls'  
     
%  learning_parameters(i,:,:)=[alpha_tt(model.FL_indices)...
%    beta_tt(model.FL_indices,model.BL_indices) dd_tt(model.FL_indices,:) ee_tt(model.FL_indices,:)];
  learning_parameters(i,:,:)=[alpha_tt(model.FL_indices) diag(beta_tt(model.FL_indices,model.FL_indices)) ];
   beta_all(:,i) = diag(beta_tt(model.FL_indices,model.FL_indices));
  alpha_all(:,i) = learning_parameters(i,:,1);
   %  learning_covar(i,:,:,:)=rr_tt(:,:,model.FL_indices);
   
         case 'var(1)'  
     
  learning_parameters(i,:,:)=[alpha_tt(model.FL_indices) (beta_tt(model.FL_indices,model.BL_indices)) ];
   beta_all(:,:,i) = (beta_tt(model.FL_indices,model.BL_indices));
  alpha_all(:,i) = learning_parameters(i,:,1);
   %  learning_covar(i,:,:,:)=rr_tt(:,:,model.FL_indices);

         case 'ar(1),sac'
             
                  
%  learning_parameters(i,:,:)=[alpha_tt(model.FL_indices)...
%    beta_tt(model.FL_indices,model.BL_indices) dd_tt(model.FL_indices,:) ee_tt(model.FL_indices,:)];
  learning_parameters(i,:,:)=[alpha_tt(model.FL_indices) diag(beta_tt(model.FL_indices,model.FL_indices)) ];
   beta_all(:,i) = diag(beta_tt(model.FL_indices,model.FL_indices));
   alpha_all(:,i) = learning_parameters(i,:,1);
%  learning_covar(i,:,:,:)=rr_tt(:,:,model.FL_indices);


         case 'ar(2),rls'
   learning_parameters(i,:,:)=[alpha_tt(model.FL_indices) ... 
     diag(beta_tt(model.FL_indices,model.FL_indices)) ... 
       diag(beta_tt(model.FL_indices+model.numEndo,model.FL_indices+model.numEndo)) ];
%    beta_all(:,i) = diag(beta_tt(model.FL_indices,model.FL_indices));   
         case 'msv'
             learning_parameters(i,:,:)  =[alpha_tt(model.FL_indices) ...
                 (beta_tt(model.FL_indices,model.BL_indices)) ... 
                 dd_tt(model.FL_indices,:)];
                 
             alpha_all(:,i) = alpha_tt(model.FL_indices);
             beta_all(:,:,i) = beta_tt(model.FL_indices,model.BL_indices);
             dd_all(:,:,i) = dd_tt(model.FL_indices,:);
         case 'msv,restricted'
             learning_parameters(i,:,:)  =[alpha_tt(model.FL_indices) ...
                 (beta_tt(model.FL_indices,model.BL_indices)) ... 
                 dd_tt(model.FL_indices,:)];
                 
             alpha_all(:,i) = alpha_tt(model.FL_indices);
             beta_all(:,:,i) = beta_tt(model.FL_indices,model.BL_indices);
%              dd_all(:,:,i) = dd_tt(model.FL_indices,:);

     end

% beta1_tt=beta_tt(1:model.numEndo,1:model.numEndo);
% beta2_tt=beta_tt(model.numEndo+1:end,model.numEndo+1:end);
%  learning_parameters(i,:,:)=[alpha_tt(model.FL_indices) diag(beta1_tt(model.FL_indices,model.FL_indices)) diag(beta2_tt(model.FL_indices,model.FL_indices))];
% % % 
% 
gamma1_tv(i,:,:)=gamma1;
gamma2_tv(i,:,:)=gamma2;
gamma3_tv(i,:,:)=gamma3;

 
 end
 

end
%=====================
if indicator==0
 likl=-sum(likl(1+model.burnIn:model.N));
elseif indicator==1
  likl=-sum(likl(1+model.burnIn:model.N)) - model.penalty;
end
%=====================
% save kf_output_estimation_results_ble_t_1_full.mat;
% save kf_output_ar1_rls_baseline_070122.mat;
% save kf_output_ble_baseline_070122.mat;
% save kf_output_ble_fullSample_070122.mat;
save results/sac_kf_output;

%figures------
% startDate=datenum('01-01-1961');
% endDate = datenum('01-12-2008');
% Date=linspace(startDate,endDate,length(model.dataset));


% beta1_average=mean(learning_parameters(end-40:end-1,:,2));
% beta2_average=mean(learning_parameters(end-40:end-1,:,3));


% % save kf_output_workspace_sac.mat;
% if save_output==1
% string_uu=['kf_output_' string_];
% save(string_uu);
% end

% % % ar(2)
%  figure;
%  index=0;
%  for ii=1:7
% %      for jj=1:2
% jj=1;
%          index=index+1;
%          subplot(7,1,index);
%          plot(learning_parameters(3:end,ii,jj),'b');
% %      end;
%  end
%  
%   figure;
%  index=0;
%  for ii=1:7
% %      for jj=1:2
% jj=2;
%          index=index+1;
%          subplot(7,1,index);
%          plot(learning_parameters(3:end,ii,jj),'b');
% %      end;
%  end


% % ar(1)
%  figure;
%  index=0;
%  for ii=1:7;
%      for jj=1:2;
%          index=index+1;
%          subplot(7,2,index);
%          plot(learning_parameters(3:end,ii,jj));
%      end;
%  end