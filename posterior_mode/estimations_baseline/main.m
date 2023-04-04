clear;clc;%close all;
restoredefaultpath;
% addpath('helpers');
% addpath('Matfiles');
% addpath('figures');
% addpath('dynare_initial_beliefs');
% addpath('inputs');
% addpath('results');
addpath(genpath(cd));

warning('off');
names=[{'phi','sigma_c' ,'lambda' ,'xi_w' ,'sigma_l', 'xi_p' , 'iota_w','iota_p',...
    'psi' ,'phi_p','r_pi', 'rho' ,'r_y', 'r_dy' ,...
           'pi_bar' ,'beta_const' ,'l_bar', 'gamma_bar' ,'alpha'...
           'rho_a', 'rho_b' ,'rho_g' ,'rho_i' ,'rho_r', 'rho_p', 'rho_w','rho_ga',...
           'eta_a', 'eta_b' ,'eta_g' ,'eta_i' ,'eta_r', 'eta_p', 'eta_w',...
           'gain'} ] ;  

%---------------------------------------------------------------------------
%optimization initial values
load param_init; 
param_init=reshape(param_init,[1,length(param_init)]);
%init_H=nhess_diagonal(@likelihood,param_init');
%init_H=inv(init_H);init_H=diag(diag(init_H)); 
%save init_H;

%---------------------------------------------------------------------------
%---------------------------------------------------------------------------
%optimization input values for model
% model=[];
% load('slobodyan_dataset.mat');%this dataset ends in 2008Q4
% % first_obs=71;%1965Q1
% first_obs = 141;
% last_obs=length(dy)-1;%2008Q4
% %load('full_dataset.mat');first_obs=106;last_obs=length(dy);
% % model.dataset=[dy dc dinve dw pinfobs robs labobs];
% model.dataset=[dy dc dinve dw cpi_quarterly robs labobs];
% model.dataset=model.dataset(first_obs:last_obs,:);


model=[];
load('dynare_initial_beliefs/slobodyan_dataset.mat');
first_obs=71;%1965Q1
last_obs=length(dy);%2008Q4
%load('full_dataset.mat');first_obs=106;last_obs=length(dy);
model.dataset=[dy dc dinve dw cpi_quarterly robs labobs];
model.dataset=model.dataset(first_obs:last_obs,:);


%this is data between 1965Q1 2008Q4 
%expectations data will start from 1982Q3. Shorten the sample to
%accommodate that
% model.dataset = model.dataset(71:end,:);

% load inflation_expectations_cleaned;
% 
% cpi_quarterly = cpi_annualized(1:106)/4 ; % quarterterly rate, 82Q3-2008Q4
% %switch pinfobs (gdp deflator) with cpi inflation 
% model.dataset(:,5) = cpi_quarterly;


model.N=length(model.dataset);%dataset length
model.l=7;%dataset # of variables
model.burnIn=4;%presample data for KF initialization
model.numVar=24;% #model size
model.numShocks=7;% #number of shocks
model.numEndo=17;% #number of endogenous variables
model.numExo=7;% #number of exogenous variables
model.numBackward=7;% #number of backward-looking variables
model.numForward=7;% #number of forward-looking  variables
model.BL_indices=[6 7 8 10 11 12 13];% backward-looking variables indices in PLM
model.FL_indices=[3 5 6 7 9 10 11];% forward-looking  variables indices 
model.shock_indices=18:1:24;
model.shock_indices_lagged=[1 2 3 4 5 6 7];
model.LB=[0.5 0.1 0.001 0.5 0.25 0.5 0.01 0.01 0.01 1 1 0.5 0.001 0.001 0.1 0.01 -10 0.1 0.01 0.01 0.01 0.01 0.01 0.01 ...
    0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.01 0.00001]';%lower bound on estimates
model.UB=[15 3 0.99 0.99 10 0.99 0.99 0.99 1 3 3 0.999 0.5 0.5 2 2 10 0.8 1 0.999 .999 .999 .999 .999 .999 ...
    .999 .999 10 10 10 10 10 10 10 0.2]';%upper bound on estimates
model.PLM='ar(1),t-1,sac';%options: 'msv,t','msv,t-1',ar(1),t,rls','ar(1),t,sac','ar(1),t-1,rls','ar(1),t-1,sac','ar(2),t-1,rls'
model.learning=1; %options: 1=adaptive learning invoked, 0=beliefs fixed at initial values 
model.learning_algo='ar(1),sac'; %options: 'msv','ar(1),rls','ar(1),sac','msv,restricted';
model.initial_beliefs='ar(1),ble-based,fixed';%options: 'msv,ree-based,fixed','ar(1),ree-based,fixed','ar(1),ree-based','ar(1),ble-based','ar(1),ble-based,fixed'
model.PLM_timing='t-1';%this is repetitive, same info in model.PLM, fix it later
model.projection_facility=0;%=1 impose PR facility, =0 do not
model.projection_facility_PLM=0;%=1 imposes projection facility on the PLM on top of IALM
model.N_fixedPoint=100;
model.optimizer=1; %options: 1=fminsearch, 2=csminwel 3=fmincon,4=patternsearch,5=particleswarm
model.ridge_correction=0;%invokes ridge correction when updating beliefs if =1
model.compute_init_hess=0;%=0 no hessian, =1 only diagonal, =2 full hessian
model.hetero_gains=0;%!!!!!keep at =0 for now, =1 not working
model.decreasing_gain=0;
model.penalty = 0;
% model.est_indices=1:1:34;
% model.excl_indices=35;
%---------------------------------------------------------------------------
% load AR1_initial_beliefs.mat;
% model.beta_init_FP=0*eye(model.numVar);
% % model.beta_init_FP(model.FL_indices,model.FL_indices)=...
% %     diag([0.98,0.52,0.98,0.98,0.98,0.92,0.98]);
% model.beta_init_FP(model.FL_indices,model.FL_indices)=...
%     diag([0.98,0.52,0.98,0.98,0.98,0.85,0.98]);
%---------------------------------------------------------------------------
objective=@(x) likelihood(x,model);
% options=optimset('Display','iter','PlotFcns',@optimplotfval,'MaxIter',20000,'UseParallel',true);
options=optimset('Display','iter','MaxIter',1000,'UseParallel',true);
options_ps =optimoptions('particleswarm','SwarmSize',500,...
    'InitialSwarmMatrix',repmat(param_init,[500 1]),'UseParallel',true);

switch model.compute_init_hess
    case 2
        init_H=nhess(objective,param_init');
        init_H=inv(init_H);
        save init_H.mat init_H;
    case 1
        init_H=nhess_diagonal(objective,param_init');
        init_H=inv(init_H);init_H=diag(diag(init_H)); 
        save init_H.mat init_H;
    case 0
        load init_H.mat;
end


save options_initial.mat model param_init init_H;


switch model.optimizer
    case 1
[x fh]=fminsearch(objective,param_init',options);
    case 2
[fh,x,gh,H,itct,fcount,retcodeh] = csminwel(objective,param_init,init_H,[] ,10^(-6),9999);
    case 3
[x,fh,EXITFLAG,OUTPUT,LAMBDA,GRAD,H] =fmincon(objective,param_init,[],[],[],[],[],[],[],options);
    case 4
[x,fh,EXITFLAG,OUTPUT]=patternsearch(objective,param_init,[],[],[],[],[],[],[],options);     
    case 5
     x=   particleswarm(objective,length(param_init),model.LB,model.UB,options_ps);
end

% [~,model]=likelihood(x,model);
% fh=fh-model.prior_excl;

x=reshape(x,[length(x),1]);
param_init=x;save param_init.mat param_init;
    if model.optimizer == 2
st_dev=sqrt(diag(H));
init_H=H;save init_H.mat init_H;
% laplace1=laplace_approximator(fh,x(model.est_indices),H(model.est_indices,model.est_indices))
laplace1=laplace_approximator(fh,x,H);
results=table(names',x,st_dev)
disp(laplace1);
    else 
results=table(names',x)   
    end


save results/estimation_results.mat;



%========================================================
%========================================================

% Sigma=nhess(@likelihood,x);
% Sigma=inv(Sigma);
% laplace2=laplace_approximator(fh,x,Sigma)
% %smoothed_variables(x);
% smoother;

