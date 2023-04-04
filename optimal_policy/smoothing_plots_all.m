clear;

ble_db = load('db_ble_smoothing.mat');
sac_db = load('db_sac_smoothing.mat');
msv_db = load('db_msv_smoothing.mat');
ree_db = load('db_ree_smoothing.mat');

para_grid = ble_db.db_tmp{1}.opt.param_grid2;
trim=3;

db_names={'ble_db','sac_db','msv_db','ree_db'};
plot_names={'BLE','SAC','MSV','REE'};

fig=figure;
fig.Position=[3000 600 800 600];
var_list={'Output','Inflation','r'};
var_list_plot={' y_t','\pi_t','r_t','E[\pi_t r_{t-1}]'};


    for ii=1:length(var_list)   
          subplot(2,2,ii)
        for jj=1:length(db_names)
            
            dat_tmp=100*diff(log(sqrt(eval([db_names{jj},'.db_small.',var_list{ii}]))));
            if jj==1
             plot(para_grid(2:end-trim),dat_tmp(1:end-trim),'LineWidth',2,'color','red');
            elseif jj==2 
               plot(para_grid(2:end-trim),dat_tmp(1:end-trim),'--','LineWidth',2,'color','red');   
            elseif jj==3
                 plot(para_grid(2:end-trim),dat_tmp(1:end-trim),'--','LineWidth',2,'color','black');   
            elseif jj==4
              plot(para_grid(2:end-trim),dat_tmp(1:end-trim),'LineWidth',2,'color','black');   
            end
    hold on;
    plot(0:0.01:1,0*ones(101,1),'-','LineWidth',0.5,'color','black');
        end
%         if ii==length(var_list)
%         legend(plot_names);
%         end
        
        title(var_list_plot(ii));
        ylabel('% change in st. dev');
        xlabel('Smoothing \rho');
    end
    
    
subplot(2,2,4);
for jj=1:length(db_names)
    dat_tmp=eval([db_names{jj},'.','corr_rpi']);
    if jj==1
            plot(para_grid(1:end-trim),dat_tmp(1:end-trim),'lineWidth',2,'color','red');
    elseif jj==2 
        plot(para_grid(1:end-trim),dat_tmp(1:end-trim),'--','lineWidth',2,'color','red'); 
    elseif jj==3
          plot(para_grid(1:end-trim),dat_tmp(1:end-trim),'--','lineWidth',2,'color','black'); 
    elseif jj==4 
       plot(para_grid(1:end-trim),dat_tmp(1:end-trim),'lineWidth',2,'color','black');
    end
    hold on;
end
ylabel('Correlation');
 xlabel('Smoothing \rho');
 hl=legend(plot_names,'Location','SouthWest','Orientation','horizontal','Position', [ 0.336388888888889       0.00351851978511722         0.392083338244756         0.031666667620341]);

title(var_list_plot{4});
 generate_figures('smoothing_figures_all','figures');