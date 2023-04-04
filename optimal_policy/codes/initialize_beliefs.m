function [alpha_tt beta_tt dd_tt ee_tt rr_tt,indicator]=initialize_beliefs(parameters,model,sysmat)
switch model.PLM
    case 'ar(2),t-1,rls'
        
beta1_tt=0*eye(model.numEndo);
beta2_tt=0*eye(model.numEndo);
alpha_tt=zeros(model.numEndo,1);
dd_tt=0*ones(model.numEndo,model.numShocks);
ee_tt=0*ones(model.numEndo,model.numShocks);
indicator=0;       

    case 'ar(2),t,rls'
        
beta1_tt=0*eye(model.numEndo);
beta2_tt=0*eye(model.numEndo);
alpha_tt=zeros(model.numEndo,1);
dd_tt=0*ones(model.numEndo,model.numShocks);
ee_tt=0*ones(model.numEndo,model.numShocks);
indicator=0;  


    otherwise
beta_tt=0*eye(model.numEndo);
alpha_tt=zeros(model.numEndo,1);
dd_tt=0*ones(model.numEndo,model.numShocks);
ee_tt=0*ones(model.numEndo,model.numShocks);
indicator=0;
end
%=========================================================================
%=========================================================================
switch model.initial_beliefs
    case 'msv,ree-based'
%           load('MSV_initial_beliefs.mat');
rr_tt=nan;
[beta_init]=REE_solve_uhlig(parameters,model,sysmat);
beta_tt(model.FL_indices,model.BL_indices)=beta_init;
ee_tt=(sysmat.A-sysmat.C*beta_tt)\sysmat.Et;
vec_dd=(kron(eye(model.numExo),sysmat.A-sysmat.C*beta_tt)-...
    kron(sysmat.RHO',sysmat.C))\vec(sysmat.D+sysmat.C*ee_tt);
dd_tt=reshape(vec_dd,[model.numEndo,model.numExo]);



% load('MSV_initial_beliefs.mat');
%   rr_tt=1*rr_init;
%=========================================================================
%=========================================================================  
    case 'msv,ree-based,fixed'
  load('MSV_initial_beliefs.mat');
  rr_tt=1*rr_init;
beta_tt(model.FL_indices,model.BL_indices)=0.5*beta_init;
dd_tt(model.FL_indices,:)=dd_init;
% ee_tt(model.FL_indices,:)=ee_init;
%=========================================================================
%=========================================================================  
    case 'ar(1),ree-based,fixed'
     load('AR1_initial_beliefs.mat');
    rr_tt=repmat(eye(2),[1 1 model.numVar]);
for jj=1:length(beta_init)
   rr_tt(:,:,model.FL_indices(jj))=rr_init(:,:,jj);
    beta_tt(model.FL_indices(jj),model.FL_indices(jj))=beta_init(jj);

end   
%=========================================================================
%=========================================================================  
    case 'ar(2),ree-based,fixed'
     load('AR2_initial_beliefs.mat');
    rr_tt=repmat(eye(3),[1 1 model.numVar]);
for jj=1:length(beta1_init)
   rr_tt(:,:,model.FL_indices(jj))=rr_init(:,:,jj);
    beta1_tt(model.FL_indices(jj),model.FL_indices(jj))=beta1_init(jj);
    beta2_tt(model.FL_indices(jj),model.FL_indices(jj))=beta2_init(jj);
 
end
beta_tt=([beta1_tt,zeros(model.numEndo,model.numEndo);zeros(model.numEndo,model.numEndo),beta2_tt]);

    case 'ar(2),ble-based,fixed'
        load AR1_BLE_initial_beliefs.mat;
 rr_tt=repmat(eye(3),[1 1 model.numVar]);
beta1_tt(model.FL_indices,model.FL_indices)=diag(beta_init);
rr_tt(2,2,model.FL_indices)=rr_init(2,2,:);
rr_tt(3,3,model.FL_indices)=rr_init(2,2,:);

beta_tt=([beta1_tt,zeros(model.numEndo,model.numEndo);zeros(model.numEndo,model.numEndo),beta2_tt]);
%=========================================================================


    case 'ar(1),ree-based'
       rr_tt=repmat(eye(2),[1 1 model.numVar]);
        %do nothing here, this is taken care of inside main loop
%=========================================================================  
    case 'ar(1),ble-based'
   
        try
[beta_init,indicator,vec0,vec1,model]=SW_fixedPoint(model,sysmat);
        catch 
  indicator=1;
  beta_init=2*eye(length(beta_tt));
  vec0=nan;
  vec1=nan;
        end

% beta_init=diag(beta_init);
 for jj=1:length(model.FL_indices)
   rr_tt=repmat(eye(2),[1 1 model.numVar]);
   rr_tt(2,2,model.FL_indices)=diag(vec0(model.FL_indices,model.FL_indices));
    beta_tt(model.FL_indices(jj),model.FL_indices(jj))=beta_init(model.FL_indices(jj),model.FL_indices(jj));
 end
 %========================================================================
    case 'ar(1),ble-based,fixed'
        load AR1_BLE_initial_beliefs.mat;
rr_tt=repmat(eye(2),[1 1 model.numVar]);
beta_tt(model.FL_indices,model.FL_indices)=0.5*diag(beta_init);
rr_tt(:,:,model.FL_indices)=0.01*rr_init;
               
%=========================================================================
    case 'msv,ree-based,fixed,forecast'
        jj=model.forecast_index;
file_name=['auxiliary_forecast_model_ree' num2str(jj) '.mat'];
location=[cd,'\auxiliary_files\'];
string=[location file_name];
load(string,'alpha_tt','beta_tt','dd_tt','ee_tt','rr_tt');
 %========================================================================
    case 'ar(1),ble-based,fixed,forecast'
        jj=model.forecast_index;
file_name=['auxiliary_forecast_model_ble' num2str(jj) '.mat'];
location=[cd,'\auxiliary_files\'];
string=[location file_name];
load(string,'beta_tt','rr_tt');
               
 %========================================================================
    case 'ar(2),ble-based,fixed,forecast'
        jj=model.forecast_index;
file_name=['auxiliary_forecast_model_ble' num2str(jj) '.mat'];
location=[cd,'\auxiliary_files\'];
string=[location file_name];
load(string,'beta_tt','rr_tt');

beta_init=beta_tt;
rr_init=rr_tt;
rr_tt=repmat(eye(3),[1 1 model.numVar]);
beta1_tt(model.FL_indices,model.FL_indices)=(beta_init(model.FL_indices,model.FL_indices));
rr_tt(2,2,model.FL_indices)=rr_init(2,2,model.FL_indices);
rr_tt(3,3,model.FL_indices)=rr_init(2,2,model.FL_indices);

beta_tt=([beta1_tt,zeros(model.numEndo,model.numEndo);zeros(model.numEndo,model.numEndo),beta2_tt]);
%=========================================================================
    case 'ar(1),diffuse'
          rr_tt=0.01*repmat(eye(2),[1 1 model.numVar]);
               %do nothing, leave everything at zeros
%=========================================================================
    case 'msv,diffuse'
        load msv_initial_beliefs.mat;
        rr_tt=0.1*eye(length(rr_init));
               %do nothing, leave everything at zeros
%=========================================================================  
    otherwise
        error('error: specified model.initial_beliefs does not exist');
end



% % %  load initial_rr_tt;
%   beta_tt(model.FL_indices,model.BL_indices)=beta_init;
%   dd_tt(model.FL_indices,:)=dd_init;
%   ee_tt(model.FL_indices,:)=ee_init;
% beta_prev=beta_tt(model.FL_indices,model.BL_indices);
% save beta_prev.mat beta_prev;

end