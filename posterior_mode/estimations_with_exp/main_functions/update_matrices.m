function [gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_tt,dd_tt,ee_tt)
%beta_tt is stacked with beta1_tt and beta2_tt if AR(2)
%normal sizes: gamma1 N x N, gamma2 N x 1, gamma3 N x M 
%sizes with AR(2) gamma1 2N x 2N, gamma2 2N x 1, gamma3 2N x M
%N=numVar,M=numExo





switch model.PLM
    
    case 'msv,t'
                gamma1 = sysmat.A \ sysmat.B; 
        gamma2 = sysmat.A \ sysmat.C; 
        gamma3 = sysmat.A \ sysmat.D;  
        
        
        
%         aux_matrix= [(sysmat.A-sysmat.C*beta_tt), -(sysmat.C * dd_tt + sysmat.D);...
%                     zeros(model.numExo,model.numEndo),eye(model.numExo)];
%         gamma1 = aux_matrix\[sysmat.B,zeros(model.numEndo,model.numExo);zeros(model.numExo,model.numEndo), sysmat.RHO];
%         gamma2=  aux_matrix\[sysmat.C*alpha_tt;zeros(model.numExo,1)];
%         gamma3=aux_matrix\[zeros(model.numEndo,model.numExo);sysmat.Ft];

%=========================================================================
%=========================================================================
    case 'var(1),t-1'
        
        gamma1 = sysmat.A \ sysmat.B; 
        gamma2 = sysmat.A \ sysmat.C; 
        gamma3 = sysmat.A \ sysmat.D;  
        
%         
%     aux_matrix=[sysmat.A,-(sysmat.C*(beta_tt*dd_tt+dd_tt*sysmat.RHO+ee_tt)+sysmat.D);...
%         zeros(model.numExo,model.numEndo),eye(model.numExo)];
%     gamma1=aux_matrix\[sysmat.B+sysmat.C*beta_tt^2,sysmat.C*beta_tt*ee_tt+sysmat.Et;...
%         zeros(model.numExo,model.numEndo),sysmat.RHO];
%     gamma3=aux_matrix\[zeros(model.numEndo,model.numExo);sysmat.Ft];
%     gamma2=aux_matrix\[sysmat.C*(alpha_tt+beta_tt*alpha_tt);zeros(model.numExo,1)];  



    case 'msv,t-1'
        
             gamma1 = sysmat.A \ sysmat.B; 
        gamma2 = sysmat.A \ sysmat.C; 
        gamma3 = sysmat.A \ sysmat.D; 
        
%         
%     aux_matrix=[sysmat.A,-(sysmat.C*(beta_tt*dd_tt+dd_tt*sysmat.RHO+ee_tt)+sysmat.D);...
%         zeros(model.numExo,model.numEndo),eye(model.numExo)];
%     gamma1=aux_matrix\[sysmat.B+sysmat.C*beta_tt^2,sysmat.C*beta_tt*ee_tt+sysmat.Et;...
%         zeros(model.numExo,model.numEndo),sysmat.RHO];
%     gamma3=aux_matrix\[zeros(model.numEndo,model.numExo);sysmat.Ft];
%     gamma2=aux_matrix\[sysmat.C*(alpha_tt+beta_tt*alpha_tt);zeros(model.numExo,1)];  
   
% %=========================================================================
%=========================================================================    
%     case 'ar(1),t,sac'
%         aux_matrix=[sysmat.A-sysmat.C*beta_tt,-sysmat.D;zeros(model.numExo,model.numEndo),eye(model.numExo)];
%         gamma1=aux_matrix\[sysmat.B,sysmat.Et;zeros(model.numExo,model.numEndo),sysmat.RHO];
%         gamma3=aux_matrix\[zeros(model.numEndo,model.numExo);sysmat.Ft];
%         gamma2=aux_matrix\[sysmat.C*(alpha_tt-beta_tt*alpha_tt);zeros(model.numExo,1)];
%=========================================================================
%=========================================================================       
        
    case 'ar(1),t-1,sac'
        
        gamma1 = sysmat.A \ sysmat.B; 
        gamma2 = sysmat.A \ sysmat.C; 
        gamma3 = sysmat.A \ sysmat.D; 
        
%         aux_matrix=[sysmat.A,-sysmat.D;zeros(model.numExo,model.numEndo),eye(model.numExo)];
%         gamma1=aux_matrix\[sysmat.B+sysmat.C*beta_tt^2,sysmat.Et;zeros(model.numExo,model.numEndo),sysmat.RHO];
%         gamma3=aux_matrix\[zeros(model.numEndo,model.numExo);sysmat.Ft];
%         gamma2=aux_matrix\[sysmat.C*(alpha_tt-beta_tt^2*alpha_tt);zeros(model.numExo,1)];
% gamma2=aux_matrix\[sysmat.C*(alpha_tt);zeros(model.numExo,1)];
%=========================================================================
%=========================================================================        
%     case 'ar(1),t,rls'
%       aux_matrix=[sysmat.A-sysmat.C*beta_tt,-sysmat.D;zeros(model.numExo,model.numEndo),eye(model.numExo)];
%       gamma1=aux_matrix\[sysmat.B,sysmat.Et;zeros(model.numExo,model.numEndo),sysmat.RHO];
%      gamma3=aux_matrix\[zeros(model.numEndo,model.numExo);sysmat.Ft];
%      gamma2=aux_matrix\[sysmat.C*(alpha_tt);zeros(model.numExo,1)];        
%=========================================================================
%=========================================================================        
    case 'ar(1),t-1,rls'
        
        gamma1 = sysmat.A \ sysmat.B; 
        gamma2 = sysmat.A \ sysmat.C; 
        gamma3 = sysmat.A \ sysmat.D; 
        
       
    
       
 %=========================================================================
%=========================================================================
    case 'ar(2),t-1,rls'
        
       gamma1 = sysmat.A \ sysmat.B; 
        gamma2 = sysmat.A \ sysmat.C; 
        gamma3 = sysmat.A \ sysmat.D;  
        
% beta1_tt=beta_tt(1:model.numEndo,1:model.numEndo);
% beta2_tt=beta_tt(model.numEndo+1:end,model.numEndo+1:end);      
%         
%          aux_matrix=[sysmat.A,-sysmat.D;zeros(model.numExo,model.numEndo),eye(model.numExo)];  
% gamma1_tilde=aux_matrix\[sysmat.B+sysmat.C*(beta1_tt^2+beta2_tt),sysmat.Et;zeros(model.numExo,model.numEndo),sysmat.RHO];
% gamma2_tilde=aux_matrix\[sysmat.C*(alpha_tt+beta1_tt*alpha_tt);zeros(model.numExo,1)];  %intercept term    
% gamma3_tilde=aux_matrix\[zeros(model.numEndo,model.numExo);sysmat.Ft];
% gamma4_tilde=aux_matrix\[sysmat.C*beta1_tt*beta2_tt, zeros(model.numEndo,model.numExo);...
%     zeros(model.numExo,model.numEndo),zeros(model.numExo,model.numExo)];
% 
% gamma1=[gamma1_tilde,gamma4_tilde;eye(model.numVar,model.numVar),zeros(model.numVar,model.numVar)];
% gamma2=[gamma2_tilde;zeros(model.numVar,1)];
% gamma3=[gamma3_tilde;zeros(model.numVar,model.numExo)];


%======================================
%=========================================================================
    case 'ar(2),t,rls'
        
% beta1_tt=beta_tt(1:model.numEndo,1:model.numEndo);
% beta2_tt=beta_tt(model.numEndo+1:end,model.numEndo+1:end);      
%         
%   aux_matrix=[sysmat.A-sysmat.C*beta1_tt,-sysmat.D;zeros(model.numExo,model.numEndo),eye(model.numExo)];
% 
%       gamma1=aux_matrix\[sysmat.B+sysmat.C*beta2_tt,sysmat.Et;zeros(model.numExo,model.numEndo),sysmat.RHO];
%      gamma3=aux_matrix\[zeros(model.numEndo,model.numExo);sysmat.Ft];
%      gamma2=aux_matrix\[sysmat.C*(alpha_tt);zeros(model.numExo,1)]; 

%======================================

   

%====================================
    otherwise
        error('error: specified model.PLM does not exist');
end



end



