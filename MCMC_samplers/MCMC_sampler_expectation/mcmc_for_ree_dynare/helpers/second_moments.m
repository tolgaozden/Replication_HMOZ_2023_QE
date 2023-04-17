function [gam_0eps,gam_1eps,gam_0xeps,gam_0x,gam_1x,init_R]=...
    second_moments(beta_tt,dd_tt,ee_tt,sysmat,model)
RHO=sysmat.RHO;
Et=sysmat.Et(model.BL_indices,:);
Ft=sysmat.Ft;
sigma=sysmat.sigma;
vec_sigma=vec(sigma);

vec_0eps=(eye(length(RHO)^2)-kron(RHO,RHO))\(kron(Ft,Ft)*vec_sigma);
ll=sqrt(length(vec_0eps));
gam_0eps=reshape(vec_0eps,[ll ll]);

vec_1eps=(kron(eye(length(RHO)),RHO))*vec_0eps;
ll=sqrt(length(vec_1eps));
gam_1eps=reshape(vec_1eps,[ll ll]);

vec_0xeps=(eye(length(RHO)^2)-kron(RHO,beta_tt))\(kron(eye(length(dd_tt)),dd_tt)+kron(RHO,ee_tt))*vec_0eps;
ll=sqrt(length(vec_0xeps));
gam_0xeps=reshape(vec_0xeps,[ll ll]);

vec_0x=(eye(length(beta_tt)^2)-kron(beta_tt,beta_tt))...
    \(kron((dd_tt*RHO+ee_tt),beta_tt)*vec_0xeps+...
    kron(eye(length(dd_tt)),dd_tt)*vec(gam_0xeps')+...
    kron(beta_tt,ee_tt)*vec(gam_0xeps')+kron(dd_tt,ee_tt)*vec(gam_1eps')+...
    kron(ee_tt,ee_tt)*vec(gam_0eps'));
ll=sqrt(length(vec_0x));
gam_0x=reshape(vec_0x,[ll ll]);

gam_1x=beta_tt*gam_0x + (dd_tt*RHO+ee_tt)*gam_0xeps;

init_R=[1,zeros(1,size(beta_tt,2)),zeros(1,size(dd_tt,2)),zeros(1,size(ee_tt,2));
    zeros(length(beta_tt),1),gam_0x, gam_0xeps*RHO',gam_0xeps;
    zeros(length(dd_tt),1), RHO*gam_0xeps',gam_0eps,gam_1eps;
    zeros(length(ee_tt),1),gam_0xeps',gam_1eps',gam_0eps];

end