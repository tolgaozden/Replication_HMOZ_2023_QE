   function [BLE,indicator,vec0,vec1,model]=SW_fixedPoint(model,sysmat)

indicator=0;

eps_conv=10e-8;
%matrices with period t-1 dating
aux=[sysmat.A,-sysmat.D;zeros(model.numExo,model.numEndo),eye(model.numExo)];
gamma1=aux\[sysmat.B,sysmat.Et;zeros(model.numExo,model.numEndo),sysmat.RHO];
gamma2=aux\[sysmat.C,zeros(model.numEndo,model.numExo);zeros(model.numExo,model.numEndo),zeros(model.numExo,model.numExo)];
gamma3=aux\[zeros(model.numEndo,model.numExo);sysmat.Ft];
gamma3_tilde=gamma3;


sigma_vec=reshape(sysmat.sigma,[length(sysmat.sigma)^2,1]);
beta=zeros(model.numVar,model.numVar,model.N_fixedPoint);
%beta(:,:,1)=model.beta_init_FP;
beta(:,:,1)=0.5*eye(model.numVar);
dist=1;
% i=0;
% while dist>eps_conv
msgstr=[];
msgid=[];
lastwarn('Success');

for i=1:model.N_fixedPoint
    
%  disp(i)
    
    
    
% dist_total(i)=dist;
%     i=i+1;
    betaAux = beta(:,:,i);
    
    %period t-1 dating: 
    
    switch model.PLM_timing
        case 't-1'
    
    M=gamma1+gamma2*betaAux^2;
         case 't'    
    M=(eye(model.numVar)-gamma2*betaAux)\gamma1;
    gamma3=(eye(model.numVar)-gamma2*betaAux)\gamma3_tilde;
    otherwise
        error('error: specified model.PLM_timing does not exist');
            
            
    end
    
    %period t dating: 
%     M=(eye(numVar)-gamma2*betaAux)\gamma1;
%     gamma3=(eye(numVar)-gamma2*betaAux)\gamma3;
    

%% calculate moments - bruce force 
% vec0=(eye(model.numVar^2)-kron(M,M))\kron(gamma3,gamma3)*sigma_vec;
%     vec1=(kron(eye(model.numVar),gamma1)+kron(eye(model.numVar),gamma2*betaAux^2))*vec0;
%      for j=1:model.numVar
%      beta(j,j,i+1)=vec1( (j-1)*model.numVar+j)/vec0( (j-1)*model.numVar+j);
% %             if beta(j,j,i+1)>.9999
% %          beta(j,j,i+1)=rand*scale;
% %             end
%      end
%         
%% calculate moments -- lyapunov equation 
lyap_a = M;
% lyap_b = gamma3 * sqrt(varCovar);
lyap_b = gamma3 * sysmat.sigma * gamma3';
% lyap_b= gamma3;
% try
% gamm0_chol= dlyapchol(lyap_a,lyap_b);
gamm0= dlyap(lyap_a,lyap_b);
gamm1 = M * gamm0;
beta(:,:,i+1)=diag(diag(gamm1)./diag(gamm0));
% catch
% disp('No stable solution at current parameter values, discarding the draw.');
% lyap_a=0.99*eye(numVar);%-rand*eye(numVar);
% gamm0= dlyap(lyap_a,lyap_b);
% end
% gamm0 = gamm0_chol' * gamm0_chol;






%%
     
     
%     if max(eigs(beta(:,:,i+1)))>.99999
%         draw=rand(24,1)*scale;
%         beta(:,:,i+1)=diag(draw);
%     end
   
    if strcmp(msgid,'MATLAB:nearlySingularMatrix')==1
indicator=1;
 break;
%     elseif sum((abs(diag(beta(:,:,i))>1)))>0
%         indicator=1;
%         break;
%  
    end
   
end




% if indicator==0
%  BLE=reshape(beta(:,:,end),[model.numVar model.numVar]);
% vec0=reshape(vec0,[sqrt(length(vec0)),sqrt(length(vec0))]);
BLE = squeeze(beta(:,:,end));
vec0=gamm0;
vec1=gamm1;
 
%    model.beta_init_FP=beta(:,:,end);
%    model.vec0_init=vec0;
   
   
% else
%    BLE=model.beta_init_FP;
%    vec0=model.vec0_init;
% end;
  
% varMat=reshape(vec0,[numVar,numVar]);

 
 
%  indicator=0;
%  if sum(diag(BLE)>1)>0
%   BLE(BLE>1)=.9999;
%   indicator=1;
%  elseif sum(diag(BLE)<-1)>0 
%   BLE(BLE<-1)=-.9999;
%   indicator=1;
%  end

    
% indicator=0;
% for j=1:numVar
%     if beta(j,j,end)>1
%         indicator=1;
%     end
% end
% if indicator==1
%     BLE=0.5*eye(numVar,numVar);
% else
% BLE=beta(:,:,end);BLE=reshape(BLE,[numVar,numVar]);
% end
% BLE(3,3)=0.9;BLE(5,5)=0.85 ;BLE(6,6)=0.9;
% BLE(7,7)=0.9;BLE(9,9)=0.9;BLE(10,10)=0.4;BLE(11,11)=0.9;

end