function  [alpha_tt beta_tt dd_tt ee_tt rr_tt largest_eig pr_flag...
    alpha_old beta_old dd_old ee_old rr_old]=update_beliefs(parameters,model,sysmat,...
    alpha_tt,beta_tt,dd_tt,ee_tt,rr_tt,S,i,gain)




pr_flag=0;
beta_old=beta_tt;
alpha_old=alpha_tt;
dd_old=dd_tt;
ee_old=ee_tt;
rr_old=rr_tt;
  

    switch model.learning_algo
        
        case 'msv'
            
            
%   thetaOld=[alpha_tt(model.FL_indices)...
%    beta_tt(model.FL_indices,model.BL_indices) dd_tt(model.FL_indices,:) ee_tt(model.FL_indices,model.shock_indices_lagged(1))];

  thetaOld=[alpha_tt(model.FL_indices)...
   beta_tt(model.FL_indices,model.BL_indices) dd_tt(model.FL_indices,[1 2 3 5 6 7])];

regressor=[1;S(i-1,model.BL_indices)';S(i,model.shock_indices([1 2 3 5 6 7]))'];
 regressand=S(i,model.FL_indices)'; 
 
[theta rr_tt largest_eig pr_flag] =msv_learning(regressand,regressor,thetaOld,rr_old,gain,...
    model.FL_indices,model.BL_indices,model.numEndo,model.ridge_correction);

alpha_tt(model.FL_indices,1)=theta(1,:)';
beta_tt(model.FL_indices,model.BL_indices)=theta(2:model.numBackward+1,:)';
dd_tt(model.FL_indices,[1 2 3 5 6 7])=theta(model.numBackward+2:model.numBackward+model.numShocks,:)';
% ee_tt(model.FL_indices,model.shock_indices_lagged(1))=theta(model.numBackward+model.numShocks+2:end,:)';
%check eigenvalue of PLM
switch model.projection_facility_PLM
    case 1
    
            if pr_flag==0
                largest_eig=abs(eigs(beta_tt,1));
                if largest_eig>1
                     pr_flag=1;
                end
            end

end
%=========================================================================
%=========================================================================
        case 'ar(1),sac'
            
           
            switch model.hetero_gains
                case 0
                        for jj=model.FL_indices
    [alpha_tt(jj),beta_tt(jj,jj),rr_tt(2,2,jj)]=...
    cgl_learning_recursive(S(i,jj),S(i-1,jj),alpha_tt(jj),beta_tt(jj,jj),rr_tt(2,2,jj),gain);
largest_eig=0;
pr_flag=0;
                           end
                        
                case 1
                    index_gain=0;
                           for jj=model.FL_indices
                               index_gain=index_gain+1;
    [alpha_tt(jj),beta_tt(jj,jj),rr_tt(2,2,jj)]=...
    cgl_learning_recursive(S(i,jj),S(i-1,jj),alpha_tt(jj),beta_tt(jj,jj),rr_tt(2,2,jj),gain(index_gain),model);
largest_eig=0;
pr_flag=0;
                           end
            end
                    
                        
  %check eigenvalue of PLM                      
switch model.projection_facility_PLM
    case 1
    
            if pr_flag==0
                largest_eig=abs(eigs(beta_tt,1));
                if largest_eig>1
                     pr_flag=1;
                end
            end

end
%=========================================================================
%=========================================================================

        case 'ar(1),rls'
            
%             rr_old=rr_tt;
                           for jj=model.FL_indices
                                    thetaOld=[alpha_tt(jj) beta_tt(jj,jj)]; 
%                                     thetaOld=[alpha_tt(jj)-beta_tt(jj,jj)*alpha_tt(jj) beta_tt(jj,jj)]; 
%                                     regressor=[1-S(i-1,jj);S(i-1,jj)'];
                                    regressor=[1;S(i-1,jj)'];
                                    regressand=S(i,jj)'; 
                    [theta rr_tt(:,:,jj) largest_eig pr_flag] =msv_learning(regressand,regressor,thetaOld,rr_old(:,:,jj),gain,...
    model.FL_indices,model.BL_indices,model.numEndo,model.ridge_correction);
alpha_tt(jj)=theta(1);
%alpha_tt(jj)=theta(1)-theta(1)*theta(2);
% alpha_tt(jj)=theta(1)/(1-theta(2));
beta_tt(jj,jj)=theta(2);
                           end
                               
  %check eigenvalue of PLM                      
switch model.projection_facility_PLM
    case 1
    
            if pr_flag==0
                largest_eig=abs(eigs(beta_tt,1));
                if largest_eig>1
                     pr_flag=1;
                end
            end

end

%===================================================
%===================================================

        case 'ar(2),rls'
            
%             rr_old=rr_tt;
root1=[];
root2=[];
beta1_tt=beta_tt(1:model.numEndo,1:model.numEndo);
beta2_tt=beta_tt(model.numEndo+1:end,model.numEndo+1:end);

                           for jj=model.FL_indices
                                     thetaOld=[alpha_tt(jj) beta1_tt(jj,jj) beta2_tt(jj,jj)]; 
%                                     thetaOld=[alpha_tt(jj)-beta_tt(jj,jj)*alpha_tt(jj) beta_tt(jj,jj)]; 
%                                     regressor=[1-S(i-1,jj);S(i-1,jj)'];
                                    regressor=[1;S(i-1,jj)'; S(i-2,jj)'];
                                    regressand=S(i,jj)'; 
                    [theta rr_tt(:,:,jj) largest_eig pr_flag] =msv_learning(regressand,regressor,thetaOld,rr_old(:,:,jj),gain,...
    model.FL_indices,model.BL_indices,model.numEndo,model.ridge_correction);
alpha_tt(jj)=theta(1);
%alpha_tt(jj)=theta(1)-theta(1)*theta(2);
% alpha_tt(jj)=theta(1)/(1-theta(2));
beta1_tt(jj,jj)=theta(2);
beta2_tt(jj,jj)=theta(3);

root1(jj)=(beta1_tt(jj,jj)+sqrt(beta1_tt(jj,jj)^2+4*beta2_tt(jj,jj)))/2;
root2(jj)=(beta1_tt(jj,jj)-sqrt(beta1_tt(jj,jj)^2+4*beta2_tt(jj,jj)))/2;

beta_tt=([beta1_tt,zeros(model.numEndo,model.numEndo);zeros(model.numEndo,model.numEndo),beta2_tt]);%convenience form
                           end
                               
  %check eigenvalue of PLM                      
switch model.projection_facility_PLM
    case 1
    
            if pr_flag==0
root1_max=max(root1);
root2_max=max(root2);
largest_eig=max(root1_max,root2_max);
                
                
                
                if largest_eig>1
                     pr_flag=1;
                end
            end

end

%========================================================================


        otherwise
        error('error: specified model.learning_algo does not exist');
    end
    
% beta_prev=beta_tt(model.FL_indices,model.BL_indices);
% save beta_prev.mat beta_prev;
end