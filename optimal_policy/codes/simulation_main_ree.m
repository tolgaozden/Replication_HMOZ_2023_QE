function [sim_results] = simulation_main(param1, param2,opt)
rng(opt.seed_);
tmp = load(opt.parameters_path);


sim_results.names_endo=[{'mc' 'zcap' 'rk' 'k' 'q'  'Consumption' 'Investment' 'Output' 'lab' 'Inflation' 'w' 'r' 'kp'   'dy' 'dc' 'dinve' 'dw' 'eps_a'  'eps_b' 'eps_g' 'eps_i'  'eps_r'  'eps_p' 'eps_w' }];
sim_results.names_exo=[{ 'Productivity'  'Risk Premium' 'Government Spending' 'Investment'  'Monetary Policy'  'Price Mark-up' 'Wage Mark-up' }];

sim_results.param=tmp.x;%for reporting only
param = tmp.x; %this vector feeds into simulation - identical to sim_results.param
% 
% if strcmp(opt.model.learning_algo,'msv')==1  %&& opt.model.learning==0
  param(1)= 5.4590   ;
   param(2)=  1.2980    ;
   param(3)= 0.7761    ;
   param(4)=  0.6917    ; 
   param(5)= 1.1086     ;
   param(6)= 0.5825    ;
   param(7)= 0.2798    ;
  param(8)=  0.1810    ;
   param(9)= 0.5341    ;
   param(10)=  1.6690   ;
   param(11)=   1.4851   ; 
    param(12)= 0.8585    ;
    param(13)= 0.0943    ; 
    param(14)= 0.1458    ;
    param(15)= 0.6602   ; 
    param(16)= 0.1498   ;  
    param(17)=  1.4755   ;  
    param(18)= 0.4012    ; 
    param(19)=  0.1727   ;  
    param(20)= 0.9108    ;
    param(21)= 0.2673     ;
     param(22)= 0.9900    ;
    param(23)= 0.8027   ;
    param(24)= 0.0623    ;
    param(25)= 0.5810   ;
    param(26)= 0.8844   ;
    param(27)= 0.5007     ;
  
    param(28)= (0.4422)    ;
    param(29)= (2.4899 )    ;
    param(30)= (0.5599 )    ;
    param(31)= (0.3782 )    ;
    param(32)= (0.2186  ) ;
    param(33)= (0.2075  ) ;
     param(34)= (0.1010 )  ;
% end

sim_results.names_param=[{'phi','sigma_c' ,'lambda' ,'xi_w' ,'sigma_l', 'xi_p' , 'iota_w','iota_p',...
    'psi' ,'phi_p','r_pi', 'rho' ,'r_y', 'r_dy' ,...
           'pi_bar' ,'beta_const' ,'l_bar', 'gamma_bar' ,'alpha'...
           'rho_a', 'rho_b' ,'rho_g' ,'rho_i' ,'rho_r', 'rho_p', 'rho_w','rho_ga',...
           'eta_a', 'eta_b' ,'eta_g' ,'eta_i' ,'eta_r', 'eta_p', 'eta_w',...
           'gain'} ] ;  
       

%---------------------------------------------------------------------------
%====================================================



model = opt.model;

%if optimize over smoothing
param(11) = param1;
param(12) = param2;
%if optimize of phi_pi
% param(11)=param2;
% param(12)=param1;

param(13)=0.2;%r_y
param(14)=0.2;%r_dy

%re-order parameters
[parameters,gain,sigma]=param_set(param);

% switch opt.shock_type
%     case 'supply'
% % tmp=sigma(6,6);
% % sigma=0 * sigma;
% % sigma(6,6)=tmp;
%     case 'demand'
% % tmp=sigma(2,2);
% % sigma=0 * sigma;
% % sigma(2,2)=tmp;
%     case 'all'
% end

%over-write the gain from estimation path
gain = opt.gain ;

%====================================================
%set the system matrices: A X_t = B X_{t-1} + C X_{t+1}^e + D eps_t + Et eps_{t-1};
%eps_t = RHO eps_{t-1} + Ft eta_t; G eta_{t-1}; 
%E and F are for measurement equations
[A, B, C, D, Et ,RHO ,Ft, G ,E ,F]=SW_sysmat_MSV_filter(parameters);

sysmat=[];
sysmat.A=A;sysmat.B=B;sysmat.C=C;sysmat.D=D;sysmat.Et=Et;
sysmat.RHO=RHO;sysmat.Ft=Ft;sysmat.G=G;sysmat.E=E;sysmat.F=F;
sysmat.sigma=sigma;
%====================================================
[alpha_init,beta_init,dd_init,ee_init rr_init]=initialize_beliefs(parameters,model,sysmat);
% beta_init = 0 * beta_init;
alpha_tt=alpha_init;beta_tt=beta_init;dd_tt=dd_init;ee_tt=ee_init;rr_tt=rr_init;
[gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_tt,dd_tt,ee_tt);
%====================================================
S=zeros(model.sim_length,model.numVar);

% eta=mvnrnd(zeros(model.numShocks,1),sigma,model.sim_length);

tmp1=mvnrnd(zeros(model.numShocks,1),sigma,model.sim_length);

eta = 0 * tmp1; 

switch opt.shock_type
    case 'supply'
eta(:,6) = 1*tmp1(:,6);
    case 'demand'
  eta(:,[2]) = 1*tmp1(:,[2]);
    case 'all'
  eta = tmp1;
end

pr_flag=zeros(model.sim_length,1);
largest_eig = zeros(model.sim_length,1);

for i=2:model.sim_length
%     disp(i)
%     gain=1/i;
% eta(1000:1009,1)=3*ones(10,1);
if strcmp(model.learning_algo,'msv')==1
      learning_parameters(i,:,:)=[alpha_tt(model.FL_indices)...
    beta_tt(model.FL_indices,model.BL_indices) dd_tt(model.FL_indices,:) ee_tt(model.FL_indices,1)];
else
 learning_parameters(i,:,:)=[alpha_tt(model.FL_indices) diag(beta_tt(model.FL_indices,model.FL_indices)) ];
end
 %-------------------
    S(i,:)=gamma2+ gamma1*S(i-1,:)'+gamma3*eta(i,:)';
    
     switch model.learning
     case 1
 [alpha_tt beta_tt dd_tt ee_tt rr_tt largest_eig(i) pr_flag(i)...
    alpha_old beta_old dd_old ee_old rr_old]=update_beliefs(parameters,model,sysmat,...
    alpha_tt,beta_tt,dd_tt,ee_tt,rr_tt,S,i,gain);
% 
[gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_tt,dd_tt,ee_tt);
%=======================================
                    switch model.projection_facility
                    case 1
                 [pr_flag(i),largest_eig(i),alpha_tt,beta_tt,dd_tt,ee_tt,rr_tt]=...
                projection_facility(gamma1,alpha_old,beta_old,dd_old,ee_old,rr_old,alpha_tt,beta_tt,dd_tt,ee_tt,rr_tt,pr_flag(i));
                [gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_tt,dd_tt,ee_tt);
                        otherwise 
                            %======================
                            switch model.calculate_eig 
                                case 1
                     largest_eig(i)=eigs(gamma1,1);       
                                case 0
                                    largest_eig(i)=0;
                            end
                            %======================
                            
                    end
                     %==================================
end


end


sim_results.S = S ; 
sim_results.learning_parameters = learning_parameters; 
sim_results.beta = beta_tt;
sim_results.pr_flag = pr_flag;
sim_results.largest_eig = largest_eig;



