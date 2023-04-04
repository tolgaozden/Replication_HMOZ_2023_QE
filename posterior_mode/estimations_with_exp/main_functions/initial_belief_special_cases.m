%this script calculates the second moments of rational expectations
%equilibrium. This is needed for initializing beliefs based on REE

switch model.initial_beliefs
    case 'msv,ree-based'
        
        
       
        
 [gam0,gam1,gam_0x,gam_0xeps,gam_0eps,gam_1x,gam_1xeps,gam_1eps]...
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