
clear;
clc;


opt.model.initial_beliefs='msv,ree-based,fixed';
opt.model.learning=1;



gain_grid = [0.001];
% shock_types={'supply','demand','all'};
shock_types={'all'};

index=0;
for gg=1:length(gain_grid);
    for ss=1:length(shock_types)
        index=index+1;
        

        opt.gain = gain_grid(gg);
        opt.shock_type=shock_types{ss};

file_name=['codes/simulations/',opt.model.initial_beliefs,'_gain_',num2str(opt.gain),'_learning','_',num2str(opt.model.learning),'_',opt.shock_type,'_smoothing.mat']

db_tmp{index}=load(file_name);

label{index}=[opt.shock_type,'_',num2str(opt.gain)]

    end
end

label = label';

%parse relevant variables into vector 

for jj=1:index
    
    simul_var(jj,:,:)=db_tmp{jj}.var_all';
    
 
end




%% plots 
var_list_long={'mc' 'zcap' 'rk' 'k' 'q'  'Consumption' 'Investment' 'Output' 'lab' 'Inflation' 'w' 'r' 'kp'   'dy' 'dc' 'dinve' 'dw' 'eps_a'  'eps_b' 'eps_g' 'eps_i'  'eps_r'  'eps_p' 'eps_w' };
% var_list={'q'  'Consumption' 'Investment' 'Output' 'lab' 'Inflation' 'w' 'r' };
var_list={'Output','Inflation','r'};


db_small=[];
for jj=1:length(var_list)
   
ind_= strcmp(var_list_long,var_list(jj));
ind_tmp = find(ind_==1);
db_small.(var_list{jj})=simul_var(:,:,ind_tmp)';
end

for jj=1:100
pi_tt = db_tmp{1}.SS_all(4002:end,10,jj);
r_tt1 = db_tmp{1}.SS_all(4001:end-1,12,jj);

corr_rpi(jj)=corr(pi_tt,r_tt1);

end

save db_msv_smoothing.mat;
% 
% % figure('Name','corr rpi');plot(corr_rpi);
% %%
% 
% % for jj=1:size(db_tmp{1}.sim_beta,1)
% 
% % ind_all=[1;2;3];
% ind_all=1;
% var_list_name={'\Delta y_t','\pi_t','r_t','E[\pi_t r_{t-1}]'};
% for jj=1:1
%     figure('Name',shock_types{jj});
%     for kk=1: length(fieldnames(db_small))
% %     subplot(3,5,kk);
% subplot(2,2,kk)
%     
%     plot(db_tmp{1}.opt.param_grid2(2:end-5),100*diff(log((db_small.(var_list{kk})(1:end-5,ind_all(jj,:))))),'lineWidth',2,'color','black');
%     if kk==length(fieldnames(db_small))
% %     legend(arrayfun(@num2str,gain_grid,'UniformOutput',0),'location','northwest');
%     end
%     title(var_list_name{kk});
%     end   
%         generate_figures(['ble_',shock_types{jj}],'figures');
% end
% subplot(2,2,4);
%  plot(db_tmp{1}.opt.param_grid2(1:end-5),((corr_rpi(1:end-5))),'lineWidth',2,'color','black');
%   title(var_list_name{4})
%   
%   
%   
%   % for jj=6
% %     
% %     name=['beta_',var_list_long{db_tmp{1}.opt.model.FL_indices(jj)}];
% %     
% %     for kk=1:index
% %     
% %     db_small.(name)(:,kk) = db_tmp{kk}.sim_beta(jj,:)';
% %     end
% % 
% %     var_list=[var_list name];
% % 
% % end
