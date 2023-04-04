function [ AA,BB, CC ,DD, E, F] = SW_sysmat_reducedform( parameters,model,bmat )
 delta   = parameters(1);
 G     = parameters(2);
 phi_w   =parameters(3); 
 curvp  = parameters(4);
 curvw  = parameters(5);
 
phi     =parameters(6) ;    
sigma_c =  parameters(7) ;
lambda   =  parameters(8)   ;
xi_w     =  parameters(9)    ;
sigma_l  = parameters(10);
xi_p     = parameters(11) ;
iota_w   =    parameters(12);
iota_p   =  parameters(13);
psi      =   parameters(14) ;
phi_p    =  parameters(15)  ;

r_pi     =     parameters(16);
rho      =  parameters(17) ;
r_y      = parameters(18) ;
r_dy=     parameters(19);




pi_bar   = parameters(20);
beta1_const=   parameters(21);
l_bar    =parameters(22) ;
gamma_bar =parameters(23);
alpha    =   parameters(24);

rho_a= parameters(25);
rho_b =  parameters(26);
rho_g=    parameters(27);
rho_i=  parameters(28); 
rho_r =parameters(29) ;
rho_p=  parameters(30) ;
rho_w= parameters(31);
rho_ga=     parameters(32)  ;
mu_p=0;mu_w=0;





switch model.PLM
    case 'ar(1),t-1,rls'

        
alpha_tt = bmat.alpha(model.FL_indices);
beta_tt = diag(bmat.beta(model.FL_indices,model.FL_indices));



alpha_rk=alpha_tt(1);
alpha_q=alpha_tt(2); 
alpha_c=alpha_tt(3);
alpha_i=alpha_tt(4);
alpha_l=alpha_tt(5);
 alpha_pi=alpha_tt(6);
alpha_w=alpha_tt(7);
beta_rk=beta_tt(1);
beta_q=beta_tt(2);
beta_c=beta_tt(3); 
beta_i=beta_tt(4); 
beta_l=beta_tt(5);
beta_pi=beta_tt(6);
beta_w=beta_tt(7);         
        
sysmat_SW_model_AR1;

    case 'ar(1),t-1,sac'

        
alpha_tt = bmat.alpha(model.FL_indices);
beta_tt = diag(bmat.beta(model.FL_indices,model.FL_indices));



alpha_rk=alpha_tt(1);
alpha_q=alpha_tt(2); 
alpha_c=alpha_tt(3);
alpha_i=alpha_tt(4);
alpha_l=alpha_tt(5);
 alpha_pi=alpha_tt(6);
alpha_w=alpha_tt(7);
beta_rk=beta_tt(1);
beta_q=beta_tt(2);
beta_c=beta_tt(3); 
beta_i=beta_tt(4); 
beta_l=beta_tt(5);
beta_pi=beta_tt(6);
beta_w=beta_tt(7);         
        
sysmat_SW_model_SAC;


    case 'ar(2),t-1,rls'
        
       

        
alpha_tt = bmat.alpha(model.FL_indices);
beta_tt = bmat.beta;



alpha_rk=alpha_tt(1);
alpha_q=alpha_tt(2); 
alpha_c=alpha_tt(3);
alpha_i=alpha_tt(4);
alpha_l=alpha_tt(5);
 alpha_pi=alpha_tt(6);
alpha_w=alpha_tt(7);


beta1_rk=beta_tt(1);
beta1_q=beta_tt(2);
beta1_c=beta_tt(3); 
beta1_i=beta_tt(4); 
beta1_l=beta_tt(5);
beta1_pi=beta_tt(6);
beta1_w=beta_tt(7);         

beta2_rk=beta_tt(8);
beta2_q=beta_tt(9);
beta2_c=beta_tt(10); 
beta2_i=beta_tt(11); 
beta2_l=beta_tt(12);
beta2_pi=beta_tt(13);
beta2_w=beta_tt(14);     
        
sysmat_SW_model_AR2;        



    case {'msv,t-1','var(1),t-1'}
alpha_tt = bmat.alpha(model.FL_indices);
beta_tt = (bmat.beta(model.FL_indices,model.BL_indices));
dd_tt = bmat.dd(model.FL_indices,1:end-1);%exclude iid exp shock

alpha_rk=alpha_tt(1);
 alpha_q=alpha_tt(2);
alpha_c=alpha_tt(3);
alpha_i=alpha_tt(4);
alpha_l =alpha_tt(5);
alpha_pi=alpha_tt(6);
alpha_w =alpha_tt(7);

 beta_11=beta_tt(1,1);
 beta_12=beta_tt(1,2);
 beta_13=beta_tt(1,3);
 beta_14=beta_tt(1,4);
 beta_15 =beta_tt(1,5);
 beta_16=beta_tt(1,6);
 beta_17 =beta_tt(1,7);
     beta_21=beta_tt(2,1);
     beta_22=beta_tt(2,2);
     beta_23 =beta_tt(2,3);
     beta_24 =beta_tt(2,4);
     beta_25 =beta_tt(2,5);
     beta_26 =beta_tt(2,6);
     beta_27 =beta_tt(2,7);
     beta_31=beta_tt(3,1);
     beta_32=beta_tt(3,2);
     beta_33=beta_tt(3,3);
     beta_34=beta_tt(3,4);
     beta_35=beta_tt(3,5);
     beta_36 =beta_tt(3,6);
     beta_37 =beta_tt(3,7);
     beta_41=beta_tt(4,1);
     beta_42=beta_tt(4,2);
     beta_43=beta_tt(4,3);
     beta_44=beta_tt(4,4);
     beta_45 =beta_tt(4,5);
     beta_46 =beta_tt(4,6);
     beta_47 =beta_tt(4,7);
     beta_51=beta_tt(5,1);
     beta_52 =beta_tt(5,2);
     beta_53=beta_tt(5,3);
     beta_54=beta_tt(5,4);
     beta_55=beta_tt(5,5);
     beta_56=beta_tt(5,6);
     beta_57 =beta_tt(5,7);
     beta_61=beta_tt(6,1);
     beta_62 =beta_tt(6,2);
     beta_63 =beta_tt(6,3);
     beta_64 =beta_tt(6,4);
     beta_65 =beta_tt(6,5);
     beta_66 =beta_tt(6,6);
     beta_67 =beta_tt(6,7);
     beta_71 =beta_tt(7,1);
     beta_72=beta_tt(7,2);
     beta_73 =beta_tt(7,3);
     beta_74 =beta_tt(7,4);
     beta_75=beta_tt(7,5);
     beta_76=beta_tt(7,6);
     beta_77=beta_tt(7,7);

 %d_mat, 7x7 
 
 dd_11=dd_tt(1,1);
 dd_12=dd_tt(1,2);
 dd_13=dd_tt(1,3);
 dd_14 =dd_tt(1,4);
 dd_15=dd_tt(1,5);
 dd_16=dd_tt(1,6);
 dd_17 =dd_tt(1,7);
     dd_21=dd_tt(2,1);
     dd_22=dd_tt(2,2);
     dd_23=dd_tt(2,3);
     dd_24 =dd_tt(2,4);
     dd_25 =dd_tt(2,5);
     dd_26 =dd_tt(2,6);
     dd_27 =dd_tt(2,7);
     dd_31 = dd_tt(3,1);
     dd_32 = dd_tt(3,2);
     dd_33= dd_tt(3,3);
     dd_34= dd_tt(3,4);
     dd_35= dd_tt(3,5);
     dd_36= dd_tt(3,6);
     dd_37 = dd_tt(3,7);
     dd_41  = dd_tt(4,1);
     dd_42 = dd_tt(4,2);
     dd_43 = dd_tt(4,3);
     dd_44 = dd_tt(4,4);
     dd_45 = dd_tt(4,5);
     dd_46 = dd_tt(4,6);
     dd_47  = dd_tt(4,7);
     dd_51 = dd_tt(5,1);
     dd_52= dd_tt(5,2);
     dd_53= dd_tt(5,3);
     dd_54= dd_tt(5,4);
     dd_55 = dd_tt(5,5);
     dd_56 = dd_tt(5,6);
     dd_57 = dd_tt(5,7);
     dd_61= dd_tt(6,1);
     dd_62= dd_tt(6,2);
     dd_63= dd_tt(6,3);
     dd_64= dd_tt(6,4);
     dd_65= dd_tt(6,5);
     dd_66 = dd_tt(6,6);
     dd_67 = dd_tt(6,7);
     dd_71 = dd_tt(7,1);
     dd_72 = dd_tt(7,2);
     dd_73 = dd_tt(7,3);
     dd_74 = dd_tt(7,4);
     dd_75 = dd_tt(7,5);
     dd_76 = dd_tt(7,6);
     dd_77 = dd_tt(7,7); 



sysmat_SW_model_MSV;


end