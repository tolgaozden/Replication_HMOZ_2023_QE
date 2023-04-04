clear;

ble_database = load('posterior_mode/estimations_with_exp/results/ble_kf_output.mat');
sac_database = load('posterior_mode/estimations_with_exp/results/sac_kf_output.mat');

beta_ble = ble_database.beta_tt(10,10);
beta_sac = sac_database.beta_all(6,:);

alpha_ble = 0 * ble_database.beta_tt(10,10);
alpha_sac = sac_database.alpha_all(6,:);

counter=1:1:size(beta_sac,1);
counter=counter';
xx1=datetime(1965,1,1);xx2=datetime(2008,9,1);date_tt=xx1:calmonths(3):xx2;
date_tt=date_tt';

ff=figure('Name','sac-betas');

plot(date_tt(2:end),beta_sac(2:end),'color','blue','lineWidth',2);
hold on;
plot(date_tt(2:end),beta_ble*ones(length(date_tt(2:end)),1),'color','red','lineWidth',2);
xlabel('Year');
ylabel('Perceived Persistence');
% legend('SAC-learning','BLE');
set(gca,'FontSize',15);
% generate_figures('estimates_betas','figures');


ff=figure('Name','sac-alphas');

plot(date_tt(2:end),alpha_sac(2:end),'color','blue','lineWidth',2);
hold on;
plot(date_tt(2:end),alpha_ble*ones(length(date_tt(2:end)),1),'color','red','lineWidth',2);
xlabel('Year');
ylabel('Perceived Mean');
legend('SAC-learning','BLE');
set(gca,'FontSize',15);
% generate_figures('estimates_alphas','figures');
