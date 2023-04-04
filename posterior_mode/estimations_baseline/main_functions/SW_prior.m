function [prior_]=SW_prior(param,model)
% param=param_init;
% names={'phi','sigma_c' ,'lambda' ,'xi_w' ,'sigma_l', 'xi_p' , 'iota_w','iota_p',...
%     'psi' ,'phi_p','r_pi', 'rho' ,'r_y', 'r_dy' ,...
%            'pi_bar' ,'beta_const' ,'l_bar', 'gamma_bar' ,'alpha'...
%            'rho_a', 'rho_b' ,'rho_g' ,'rho_i' ,'rho_r', 'rho_p', 'rho_w','mu_p','mu_w','rho_ga',...
%            'eta_a', 'eta_b' ,'eta_g' ,'eta_i' ,'eta_r','eta_p', 'eta_w',...
%            'gain'}  ;    

prior(1)=normpdf(param(1),4,1.5);%phi
prior(2)=normpdf(param(2),1.5,0.375);%sigma_c
% prior2=normpdf(param(2),1.5,0.15);%sigma_c

mu=0.7;sigma2=0.1^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior(3)= betapdf(param(3),a,b);%lambda
% mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
% prior3= betapdf(param(3),a,b);%lambda

mu=0.5;sigma2=0.1^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior(4)= betapdf(param(4),a,b);%xi_w
% mu=0.75;sigma2=0.1^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
% prior4= betapdf(param(4),a,b);%xi_w

prior(5)=normpdf(param(5),2,0.75);%sigma_l

%higher prior for BLE model, everything else uses standard prior.
if strcmp(model.initial_beliefs, 'ar(1),ble-based')==1 && model.learning ==0
 mu=0.75;sigma2=0.05^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
 prior(6)= betapdf(param(6),a,b);%xi_p
else
mu=0.5;sigma2=0.1^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior(6)= betapdf(param(6),a,b);%xi_p
end;

% prior for BLE model
%  mu=0.75;sigma2=0.05^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
%  prior(6)= betapdf(param(6),a,b);%xi_p


%prior for REE model & and others except for BLE
% mu=0.5;sigma2=0.1^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
% prior6= betapdf(param(6),a,b);%xi_p

mu=0.5;sigma2=0.15^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior(7)= betapdf(param(7),a,b);%iota_w

mu=0.5;sigma2=0.15^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior(8)= betapdf(param(8),a,b);%iota_p



mu=0.5;sigma2=0.15^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior(9)= betapdf(param(9),a,b);%psi

prior(10)=normpdf(param(10),1.25,0.125);%phi_p

prior(11)=normpdf(param(11),1.5,0.25);%r_pi

mu=0.75;sigma2=0.1^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior(12)= betapdf(param(12),a,b);%rho

prior(13)=normpdf(param(13),0.125,0.05);%r_y
prior(14)=normpdf(param(14),0.125,0.05);%r_dy

mu=0.625;sigma2=0.1^2;b  = sigma2/mu;a  = mu/b;
prior(15)=gampdf(param(15),a,b);

mu=0.25;sigma2=0.1^2;b  = sigma2/mu;a  = mu/b;
prior(16)=gampdf(param(16),a,b);

prior(17)=normpdf(param(17),0,2);
prior(18)=normpdf(param(18),0.4,0.1);
prior(19)=normpdf(param(19),0.3,0.05);

mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior(20)= betapdf(param(20),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior(21)= betapdf(param(21),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior(22)= betapdf(param(22),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior(23)= betapdf(param(23),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior(24)= betapdf(param(24),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior(25)= betapdf(param(25),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
prior(26)= betapdf(param(26),a,b);
mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
% prior27= betapdf(param(27),a,b);
% mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
% prior28= betapdf(param(28),a,b);

% mu=0.5;sigma2=0.2^2;a = (1-mu)*mu*mu/sigma2-mu;b = a*(1/mu-1);
% prior(27)= betapdf(param(27),a,b);
prior(27)=1;

[s,nu] = inverse_gamma_specification(0.1,4, 0,1, false, 'name');
% [s,nu] = inverse_gamma_specification(0.5,16, 0,1, false, 'name');
prior(28)= exp(lpdfig1(param(28),s,nu));
prior(29)= exp(lpdfig1(param(29),s,nu));
prior(30)= exp(lpdfig1(param(30),s,nu));
prior(31)= exp(lpdfig1(param(31),s,nu));
prior(32)= exp(lpdfig1(param(32),s,nu));

prior(33)= exp(lpdfig1(param(33),s,nu));
prior(34)= exp(lpdfig1(param(34),s,nu));

 mu=0.035;sigma2=0.015^2;b  = sigma2/mu;a  = mu/b;
 
 
switch model.hetero_gains
    case 0
        %same gain for everything
prior(35)=gampdf(param(35),a,b);%gain
% prior(35)=1;
prior(36)=1;
prior(37)=1;
prior(38)=1;
prior(39)=1;
prior(40)=1;
prior(41)=1;


    case 1
        %7 different gains
 prior(35)=gampdf(param(35),a,b);%gain       
 prior(36)=gampdf(param(36),a,b);%gain       
 prior(37)=gampdf(param(37),a,b);%gain       
 prior(38)=gampdf(param(38),a,b);%gain       
 prior(39)=gampdf(param(39),a,b);%gain       
 prior(40)=gampdf(param(40),a,b);%gain       
 prior(41)=gampdf(param(41),a,b);%gain       
end

%discard the gain prior if the model doesn't have learning.
if model.learning == 0 
    prior(35) = 1; 
    prior(36) = 1; 
    prior(37) = 1; 
    prior(38) = 1; 
    prior(39) = 1; 
    prior(40) = 1; 
    prior(41) = 1; 
end
 
 
 
 prior_=prior(1)*prior(2)*prior(3)*prior(4)*prior(5)*prior(6)*prior(7)*prior(8)*...
    prior(9)*prior(10)*prior(11)*prior(12)*...
    prior(13)*prior(14)*prior(15)*prior(16)*prior(17)*...
    prior(18)*prior(19)*prior(20)*prior(21)*prior(22)*prior(23)*prior(24)*...
    prior(25)*prior(26)*prior(27)*prior(28)*prior(29)*prior(30)*prior(31)*...
    prior(32)*prior(33)*prior(34)*prior(35)*prior(36)*prior(37)*prior(38)*...
    prior(39)*prior(40)*prior(41);
 
 
% prior_excl=1;
% for jj=1:length(model.excl_indices);
%     prior_excl=prior_excl*prior(jj);
% end
% model.prior_excl=-log(prior_excl);
        
        

prior_=-log(prior_);
end