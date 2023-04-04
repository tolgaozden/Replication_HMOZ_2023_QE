function [gam0,gam1,gam_0x,gam_0xeps,gam_0eps,gam_1x,gam_1xeps,gam_1eps]=...
    second_moments2(gamma1,gamma3,model,sysmat,beta_tt,dd_tt,ee_tt)


%gam0: covariance matrix of endogenous+exogenous variables
%gam1: first-order autocovariance matrix of endogenous+exogenous variables
%

sigma=sysmat.sigma;
RHO=sysmat.RHO;
vec0=(eye(length(gamma1)^2)-kron(gamma1,gamma1))\...
kron(gamma3,gamma3)*vec(sigma);

ll=sqrt(length(vec0));
gam0=reshape(vec0,[ll,ll]);
gam1=gamma1*gam0;

gam_0x=gam0(1:model.numEndo,1:model.numEndo);
gam_0x=gam_0x(model.BL_indices,model.BL_indices);

gam_0xeps=gam0(1:model.numEndo,model.numEndo+1:end);
gam_0xeps=gam_0xeps(model.BL_indices,:);

gam_0eps=gam0(model.numEndo+1:end,model.numEndo+1:end);

gam_1x=gam1(1:model.numEndo,1:model.numEndo);
gam_1x=gam_1x(model.FL_indices,model.BL_indices);

gam_1xeps=gam1(1:model.numEndo,model.numEndo+1:end);
gam_1xeps=gam_1xeps(model.BL_indices,:);

gam_1eps=gam1(model.numEndo+1:end,model.numEndo+1:end);


%if you have eps_{t-1} only

rr=[1,zeros(1,size(beta_tt,2)),zeros(1,size(dd_tt,2)),0;
    zeros(length(beta_tt),1),gam_0x, gam_0xeps*RHO',gam_0xeps(:,1);
    zeros(length(dd_tt),1), RHO*gam_0xeps',gam_0eps,gam_1eps(:,1);
    0,gam_0xeps(:,1)',gam_1eps(:,1)',gam_0eps(1,1)];


%if you have all eps_{t-1}
% rr=[1,zeros(1,size(beta_tt,2)),zeros(1,size(dd_tt,2)),zeros(1,size(ee_tt,2));
%     zeros(length(beta_tt),1),gam_0x, gam_0xeps*RHO',gam_0xeps;
%     zeros(length(dd_tt),1), RHO*gam_0xeps',gam_0eps,gam_1eps;
%     zeros(length(ee_tt),1),gam_0xeps',gam_1eps',gam_0eps];

end
