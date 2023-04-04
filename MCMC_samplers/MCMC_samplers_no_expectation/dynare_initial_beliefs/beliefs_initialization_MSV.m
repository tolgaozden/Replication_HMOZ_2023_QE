clear;clc;close all;
%  dynare SW_Estimation_REE;
 %dynare T_RE;
%  save('initial_ree_estimation.mat');
  load('initial_ree_estimation.mat');
%dataset=oo_.endo_simul;
% clearvars -except dataset;

% forward_indices=[3 5 6 7 9 10 11];
% shock_indices=[14 15 16 17 18 19 20];
% backward_indices=[6 7 8 10 12 13];
% lagged=dataset(backward_indices,1:end-1)';
% contemp=dataset(forward_indices,2:end)';
% shocks=dataset(shock_indices,2:end)';
burn_in=1;
shock_indices_lagged=[1 2 3 4 5 6 7];
lagged=[ c inve y pinf w r kp];
lagged=lagged(1+burn_in:end-1,:);
contemp=[rk q c inve lab pinf w];
%contemp=[rk pk c inve lab pinf w];
contemp=contemp(2+burn_in:end,:);
shocks=[eps_a eps_b eps_g eps_i eps_r eps_p eps_w];
%shocks=[a b g qs ms spinf sw];
shocks_lagged=shocks(1+burn_in:end-1,shock_indices_lagged);
shocks=shocks(2+burn_in:end,:);

const=ones(length(lagged),1);
clearvars -except lagged contemp shocks shocks_lagged const shock_indices_lagged;
regressor=[const lagged shocks(:,[1 2 3 5 6 7])];
regressand=contemp;

regression=(regressor'*regressor)^(-1) * (regressor'*regressand);
rr_init=(regressor'*regressor) * length(contemp)^(-1);

alpha_init=regression(1,:)';
beta_init=regression(2:8,:)';

dd_init=zeros(7,7);
dd_init(:,[1 2 3 5 6 7])=regression(9:14,:)';

% ee_init=zeros(size(dd_init));
% ee_init=zeros(size(dd_init,2),1)';
% ee_init(1,shock_indices_lagged)=regression(16,:)';
save MSV_initial_beliefs alpha_init beta_init dd_init  rr_init;%ee_init


