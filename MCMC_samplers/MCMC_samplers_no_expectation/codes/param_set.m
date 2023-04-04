function [parameters,gain,sigma]=param_set(param,model)

parameters(1)   = [0.025];     %delta
parameters(2)   = [0.18];    % G 
parameters(3)   =[1.5];       % phi_w 
parameters(4)   = [10];        %curv_p
parameters(5)   = [10];         %curv_w
 
 
%nominal and real frictions
parameters(6)  = param(1) ;    %phi
parameters(7)  = param(2) ;      %sigma_c
parameters(8)  = param(3);   %lambda 
parameters(9)  = param(4); %xi_w
parameters(10) = param(5);        %sigma_l 
parameters(11) = param(6);      %xi_p 
parameters(12) = param(7) ;     %iota_w
parameters(13) = param(8);         %iota_p
parameters(14) = param(9);   %psi
parameters(15) = param(10); %phi_p

%policy related parameters

parameters(16)    =  param(11);   %r_pi
parameters(17)    =  param(12); %rho
parameters(18)    =  param(13);%0.0746;    %r_y
parameters(19)    =  param(14);     %r_dy

%SS related parameters
parameters(20)    = param(15);    %pi_bar
parameters(21)    = param(16);       %beta_const
parameters(22)    = param(17);     %l_bar
parameters(23)    = param(18);         %gamma_bar
parameters(24)    = param(19);    %alpha

% %shock persistence
parameters(25) =param(20);      %rho_a
parameters(26) =param(21);  %rho_b
parameters(27) =param(22);   %rho_g
parameters(28) =param(23);   %rho_i
parameters(29) =param(24);   %rho_r
parameters(30) =param(25);  %rho_p
parameters(31) =param(26);%rho_w 

% parameters(32) =param(27);    %mu_p 
% parameters(33) =param(28);    %mu_w
parameters(32) =param(27) ;  %rho_ga

%shock standard deviations
parameters(33)=param(28);  %sigma_a
parameters(34)=param(29);  %sigma_b
parameters(35)=param(30); %sigma_g
parameters(36)=param(31) ;   %sigma_i
parameters(37)=param(32);   %sigma_r
parameters(38)=param(33);  %sigma_p
parameters(39)=param(34);   %sigma_w

switch model.hetero_gains
        case 0
        gain=param(35);
        case 1
        gain=param(35:41);
end
sigma=diag(parameters(33:39))^2;
end