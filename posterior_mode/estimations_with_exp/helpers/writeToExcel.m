clear;clc;close all;

%AR(1) learning
filename='ar1_estimation_results';
filepath = strcat('results/',filename);
load(filepath);



output_ar1=round(x,3);
output_ar1 = [output_ar1;-fh];

%SAC learning
filename='sac_estimation_results';
filepath = strcat('results/',filename);
load(filepath);


output_sac=round(x,3);
output_sac= [output_sac;-fh];



%VAR(1) learning
filename='var1_estimation_results';
filepath = strcat('results/',filename);
load(filepath);



output_var1=round(x,3);
output_var1 = [output_var1;-fh];


%MSV learning
filename='msv_estimation_results';
filepath = strcat('results/',filename);
load(filepath);



output_msv=round(x,3);
output_msv= [output_msv;-fh];

%BLE
filename='ble_estimation_results';
filepath = strcat('results/',filename);
load(filepath);

output_ble=round(x,3);
output_ble= [output_ble;-fh];

% %  

names_final=[names,{'likelihood at mode'}];

TT = table(names_final',output_ar1,output_sac, output_ble,output_var1,output_msv,'VariableNames',{'Parameter','AR(1)','SAC','BLE','VAR(1)','MSV'});
output_file='SW_Estimation_Results.csv';
writetable(TT,output_file);




%REE -- write to separate file 

filename='dynare_initial_beliefs/SW_Estimation_REE_expData_mode.mat';
filepath=filename;
load(filepath);

names_ree = parameter_names;
names_ree=[names_ree;'likelihood at mode'];

output_ree= xparam1;
output_ree=[output_ree;-fval];

TT_RE = table(names_ree,output_ree,'VariableNames',{'Parameter','REE'});
output_file='SW_Estimation_Results_REE.csv';
writetable(TT_RE,output_file);