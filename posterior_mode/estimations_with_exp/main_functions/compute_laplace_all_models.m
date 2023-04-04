clear;
clc;

db_ble=load('ble_estimation_results.mat');
db_sac=load('sac_estimation_results.mat');
db_ar2=load('ar2_estimation_results.mat');
db_var=load('var1_estimation_results.mat');
db_msv=load('msv_estimation_results.mat');


models={'db_ble','db_sac','db_ar2','db_var','db_msv'};
% models={'db_ble'};

for jj=1:length(models)
    
    disp(jj)
    
    fh_tmp=eval([models{jj},'.fh']);
    x_tmp=eval([models{jj},'.x']);
    objective=eval([models{jj},'.objective']);
    model=eval([models{jj},'.model']);
    
    init_H=nhess_diagonal(objective,x_tmp);
    % init_H=nearestSPD(init_H);
    
    if jj==1
    temp_H=inv(init_H(1:end-1,1:end-1));
    temp_H(end+1,end+1)=1;
    init_H=temp_H;
    else
      init_H=inv(init_H);
    end
        

    assignin('base',['H_',models{jj}],init_H);
    
    if jj==1
    laplace1=laplace_approximator(fh_tmp,x_tmp(1:end-1),init_H(1:end-1,1:end-1))
    else
    laplace1=laplace_approximator(fh_tmp,x_tmp(1:end),init_H(1:end,1:end))
    end
    
    assignin('base',['laplace_',models{jj}],laplace1);
end


save results/laplace_approximations.mat;

