function [initb] = beliefs_initialization_ree(db)

param=db.x;
model=db.model;

[parameters,gain,sigma]=param_set(param,model);
[A, B, C, D, Et ,RHO ,Ft, G ,E ,F]=SW_sysmat_MSV_filter(parameters);
sysmat=[];
sysmat.A=A;sysmat.B=B;sysmat.C=C;sysmat.D=D;sysmat.Et=Et;
sysmat.RHO=RHO;sysmat.Ft=Ft;sysmat.G=G;sysmat.E=E;sysmat.F=F;
sysmat.sigma=sigma;

[alpha_tt,beta_tt,dd_tt,ee_tt rr_tt,indicator]=initialize_beliefs(parameters,model,sysmat);
[gamma1 gamma2 gamma3]=update_matrices(model,sysmat,alpha_tt,beta_tt,dd_tt,ee_tt);


 [rr_tt,gam0,gam1,gam_0x,gam_0xeps,gam_0eps,gam_1x,gam_1xeps,gam_1eps]...
     =second_moments2(gamma1,gamma3,model,sysmat,...
     beta_tt(model.FL_indices,model.BL_indices),...
     dd_tt(model.FL_indices,:),ee_tt(model.FL_indices,:));
 


initb.rr_tt = rr_tt;
initb.alpha_tt=alpha_tt;
initb.beta_tt=beta_tt;
initb.dd_tt=dd_tt;

end