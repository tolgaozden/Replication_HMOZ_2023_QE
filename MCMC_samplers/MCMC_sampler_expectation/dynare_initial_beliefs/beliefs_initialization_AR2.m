clear;clc;close all;
% dynare SW_Estimation_REE;
load initial_ree_estimation;
dataset=oo_.endo_simul;
clearvars -except dataset;
forward_indices=[3 5 6 7 9 10 11];

second_lag=dataset(forward_indices,1:end-2)';
first_lag=dataset(forward_indices,2:end-1)';
contemp=dataset(forward_indices,3:end)';
const=ones(length(first_lag),1);
beta1_init=zeros(length(forward_indices),1);
beta2_init=zeros(length(forward_indices),1);
for jj=1:length(forward_indices)
    regr=([const first_lag(:,jj) second_lag(:,jj)]'*[const first_lag(:,jj) second_lag(:,jj)])^(-1)*...
        ([const first_lag(:,jj) second_lag(:,jj)]'*contemp(:,jj));
    beta1_init(jj)=regr(2);
    beta2_init(jj)=regr(3);
    rr_init(:,:,jj)=([const first_lag(:,jj) second_lag(:,jj)]'*[const first_lag(:,jj) second_lag(:,jj)])/length(const);
end

save AR2_initial_beliefs.mat beta1_init beta2_init rr_init;
