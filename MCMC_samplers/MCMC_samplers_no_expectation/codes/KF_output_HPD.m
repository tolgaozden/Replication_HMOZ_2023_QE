clear;clc;%close all;
names_endo=[{'mc' 'zcap' 'rk' 'k' 'q'  'Consumption' 'Investment' 'Output' 'lab' 'Inflation' 'w' 'r' 'kp'   'dy' 'dc' 'dinve' 'dw' 'eps_a'  'eps_b' 'eps_g' 'eps_i'  'eps_r'  'eps_p' 'eps_w' }];
 names_exo=[{ 'Productivity'  'Risk Premium' 'Government Spending' 'Investment'  'Monetary Policy'  'Price Mark-up' 'Wage Mark-up' }];
string_=('MH_results_AR2');
load(string_);

%repeat for HPD 5% and HPD 95% 
for jj=1:3

model.projection_facility=1;
model.projection_facility_PLM=1;
model.hetero_gains=0;
model.decreasing_gain=0;
save_output=1;
% model.initial_beliefs='ar(1),ble-based,fixed';
% model.learning=1;
% model.learning_algo='ar(2),rls';
% model.initial_beliefs= 'ar(2),ree-based,fixed';
% model.PLM='ar(2),t-1,rls';
% model.projection_facility=0;
% model.projection_facility_PLM=0;


% load param_init;param=param_init;
%   load estimation_results_ble_t_1;
%   load AR1_initial_beliefs.mat;
% load estimation_results_sac_t_1_ble_based_star;
% load estimation_results_msv_t_star;
% model.beta_init_FP(model.FL_indices,model.FL_indices)=0.9*diag(beta_init);
%   load param_init;
%   param=param_init;
%  load estimation_results_msv_t;
%  model.initial_beliefs='msv,ree-based'
% load estimation_results_ree_t;
% load estimation_results_rls_t_1_star.mat;
% load estimation_results_msv_t_diffuse.mat;
% load estimation_results_sac_t_1_ble_based_star;
% load estimation_results_rls_t_1_ble_based_star;
% load param_init
%  load estimation_results;
%  load estimation_results_ble_t_1_high_xi_p;
% load estimation_results_sac_t_1_ble_based_star;
% load estimation_results_rls_t_1.mat;
 %model.projection_facility_PLM=1;
%   model.beta_init_FP=0.01*eye(model.numVar);
%   model.beta_init_FP=diag(rand(model.numVar,1));
if jj==3
    param=posteriorMean;
else
 param=hpd_interval(:,jj);
end
%  param(25)=0;
%   param(6)=0.9;
%  param(25)=0;
% param_init(6)=0.8;
% param_init(25)=0.03;
% param=param_init;
% load options_initial;
% x(6)=0.69;
% param=param_init;
%  param(7)=0; param(8)=0; param(3)=0.8;
%           model.N_fixedPoint=200
[parameters,gain,sigma]=param_set(param,model);



[A, B, C, D, Et ,RHO ,Ft, G ,E ,F]=SW_sysmat_MSV_filter(parameters);

sysmat=[];
sysmat.A=A;sysmat.B=B;sysmat.C=C;sysmat.D=D;sysmat.Et=Et;
sysmat.RHO=RHO;sysmat.Ft=Ft;sysmat.G=G;sysmat.E=E;sysmat.F=F;
sysmat.sigma=sigma;

[alpha_tt,beta_tt,dd_tt,ee_tt rr_tt,indicator]=initialize_beliefs(parameters,model,sysmat);
[gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_tt,dd_tt,ee_tt);

switch model.initial_beliefs
    case 'msv,ree-based'
 [rr_tt,gam0,gam1,gam_0x,gam_0xeps,gam_0eps,gam_1x,gam_1xeps,gam_1eps]...
     =second_moments2(gamma1,gamma3,model,sysmat,...
     beta_tt(model.FL_indices,model.BL_indices),...
     dd_tt(model.FL_indices,:),ee_tt(model.FL_indices,:));
    case 'ar(1),ree-based'
        %first retrieve the auxiliary REE solution--initialize
        beta_aux=0*eye(model.numEndo);
alpha_aux=zeros(model.numEndo,1);
dd_aux=0*ones(model.numEndo,model.numShocks);
ee_aux=0*ones(model.numEndo,model.numShocks);
%-------------------------------
       [beta_init]=REE_solve_uhlig(parameters,model,sysmat);
beta_aux(model.FL_indices,model.BL_indices)=beta_init;
ee_aux=(sysmat.A-sysmat.C*beta_aux)\sysmat.Et;
vec_dd=(kron(eye(model.numExo),sysmat.A-sysmat.C*beta_aux)-...
    kron(sysmat.RHO',sysmat.C))\vec(sysmat.D+sysmat.C*ee_aux);
dd_aux=reshape(vec_dd,[model.numEndo,model.numExo]); 
      %retrieve the reduced-form matrices
      [gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_aux,dd_aux,ee_aux);
    %obtain second moments consistent with REE
   [rr_aux,gam0,gam1,gam_0x,gam_0xeps,gam_0eps,gam_1x,gam_1xeps,gam_1eps]...
     =second_moments2(gamma1,gamma3,model,sysmat,...
     beta_aux(model.FL_indices,model.BL_indices),...
     dd_aux(model.FL_indices,:),ee_aux(model.FL_indices,:));
        %take the relevant moments for AR(1) rule
        rr_tt(2,2,model.FL_indices)=diag(gam0(model.FL_indices,model.FL_indices));
        beta_tt(model.FL_indices,model.FL_indices)=...
            diag(diag(gam1(model.FL_indices,model.FL_indices))./diag(gam0(model.FL_indices,model.FL_indices)));
        
        %re-compute reduced-form matrices, now consistent with AR(1)
        [gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_tt,dd_tt,ee_tt);
        
end

%  [S(1,:),P(1,:,:)]=kalman_init(gamma1,gamma2,gamma3,sigma);
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
%     gain=1/i;
%     if model.initial_beliefs=='msv,diffuse'
%      S(model.burnIn+1,:)=zeros(model.numVar,1);   
%     end
    
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
[gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_tt,dd_tt,ee_tt);
                    switch model.projection_facility
                    case 1
                 [pr_flag(i),largest_eig(i),alpha_tt,beta_tt,dd_tt,ee_tt,rr_tt]=...
                projection_facility(gamma1,alpha_old,beta_old,dd_old,ee_old,rr_old,alpha_tt,beta_tt,dd_tt,ee_tt,rr_tt,pr_flag(i),model);
                [gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_tt,dd_tt,ee_tt);
                        otherwise 
                     largest_eig(i)=eigs(gamma1,1);       
                            
                     end


 end
 
 if model.learning==1
  learning_parameters(i,:,:)=[alpha_tt(model.FL_indices)...
    beta_tt(model.FL_indices,model.BL_indices) dd_tt(model.FL_indices,:) ee_tt(model.FL_indices,:)];
% %  learning_parameters(i,:,:)=[alpha_tt(model.FL_indices) diag(beta_tt(model.FL_indices,model.FL_indices)) ];
%  learning_covar(i,:,:,:)=rr_tt(:,:,model.FL_indices);

% beta1_tt=beta_tt(1:model.numEndo,1:model.numEndo);
% beta2_tt=beta_tt(model.numEndo+1:end,model.numEndo+1:end);
%  learning_parameters(i,:,:)=[alpha_tt(model.FL_indices) diag(beta1_tt(model.FL_indices,model.FL_indices)) diag(beta2_tt(model.FL_indices,model.FL_indices))];
% % % 
% 
% gamma1_tv(i,:,:)=gamma1;
% gamma2_tv(i,:,:)=gamma2;
% gamma3_tv(i,:,:)=gamma3;
 end
 
end
%=====================
% if indicator==0
 likl=-sum(likl(model.burnIn:model.N));
% else
%   likl=-Inf;
% end
%=====================
% save kf_output_workspace_sac_t_1.mat;



%figures------
startDate=datenum('01-01-1961');
endDate = datenum('01-12-2008');
Date=linspace(startDate,endDate,length(model.dataset));


% beta1_average=mean(learning_parameters(end-40:end-1,:,2));
% beta2_average=mean(learning_parameters(end-40:end-1,:,3));


% save kf_output_workspace_sac.mat;
if save_output==1
string_uu=['kf_output_hpd' num2str(jj) string_];
save(string_uu);
end

% % ar(2)
%  figure;
%  index=0;
%  for ii=1:7;
%      for jj=1:3;
%          index=index+1;
%          subplot(7,3,index);
%          plot(learning_parameters(3:end,ii,jj));
%      end;
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

end