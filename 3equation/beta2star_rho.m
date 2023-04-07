% close all 
clear all
clf

%%%%%%%%%% For contemporanous Taylor rule

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Our own calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0, sigma2/sigma1=0.5; 
% 
% % (\beta^*_1, \beta^*_2)=(0, 0)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)

% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.1,
% % sigma2/sigma1=0.5; 
% beta2=[0.1, 0.1118];
% T2=[0.1085, 0.1088];
% beta1=[0.1085, 0.1088];
% T1=[0.1118,0.1118];
% % (\beta^*_1, \beta^*_2)=(0.1088, 0.1118)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.2,
% % sigma2/sigma1=0.5; 
% beta2=[0.2,0.2626];
% T2=[0.2388,0.2428];
% beta1=[0.2388,0.2428];
% T1=[0.2626,0.2625];
% % (\beta^*_1, \beta^*_2)=(0.2428, 0.2625)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% % 

% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.25,
% % sigma2/sigma1=0.5; 
% beta2=[0.3,0.3679,0.3681];
% T2=[0.3185,0.3256,0.3256];
% beta1=[0.3185,0.3256];
% T1=[0.3681,0.3679];
% % (\beta^*_1, \beta^*_2)=(0.3256, 0.3679)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% % 


% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.3,
% % sigma2/sigma1=0.5; 
% beta2=[0.35,0.5179,0.5199];
% T2=[0.4025,0.4284,0.4288];
% beta1=[0.4025,0.4288];
% T1=[0.5199,0.5179];
% % (\beta^*_1, \beta^*_2)=(0.4284, 0.5179)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% 


% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.32,
% % sigma2/sigma1=0.5; 
% beta2=[0.55,0.601,0.6025];
% T2=[0.4684,0.4809,0.4813];
% beta1=[0.4684,0.4809,0.4813];
% T1=[0.6025,0.601,0.601];
% % (\beta^*_1, \beta^*_2)=(0.4809, 0.601)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% 
% 

% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.35,
% % sigma2/sigma1=0.5; 
% beta2=[0.6,0.7516,0.7638];
% T2=[0.5305,0.5851,0.5913];
% beta1=[0.5305,0.5851,0.5913];
% T1=[0.7638,0.753,0.7516];
% % (\beta^*_1, \beta^*_2)=(0.5857, 0.7529)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.38,
% % sigma2/sigma1=0.5; 
% beta2=[0.8,0.8607,0.8705];
% T2=[0.6571,0.7009,0.7098];
% beta1=[0.6571,0.702,0.7098];
% T1=[0.8705,0.8623,0.8607];
% % (\beta^*_1, \beta^*_2)=(0.7023, 0.8622)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)


% % % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.4,
% % % sigma2/sigma1=0.5; 
% beta2=[0.82,0.8962,0.9074];
% T2=[0.6968,0.759,0.7714];
% beta1=[0.6968,0.759,0.7714];
% T1=[0.9074,0.8983,0.8962];
% % (\beta^*_1, \beta^*_2)=(0.7609, 0.898)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% % % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.45,
% % % sigma2/sigma1=0.5; 
% beta2=[0.92,0.9396,0.9425];
% T2=[0.8288,0.8508,0.8544];
% beta1=[0.8288,0.8508,0.8544];
% T1=[0.9425,0.9401,0.9396];
% % (\beta^*_1, \beta^*_2)=(0.8513, 0.94)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% 

% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.5,
% % sigma2/sigma1=0.5; 
% beta2=[0.94,0.959,0.9609];
% T2=[0.8796,0.8997,0.9019];
% beta1=[0.8796,0.8997,0.9019];
% T1=[0.9609,0.9592,0.959];
% % (\beta^*_1, \beta^*_2)=(0.8999, 0.9592)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.55,
% % sigma2/sigma1=0.5; 
% beta2=[0.96,0.9711];
% T2=[0.9192,0.9301];
% beta1=[0.9192,0.9301];
% T1=[0.9711,0.9704];
% % (\beta^*_1, \beta^*_2)=(0.9294, 0.9704)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.6,
% % sigma2/sigma1=0.5; 
% beta2=[0.97,0.9783];
% T2=[0.9422,0.9493];
% beta1=[0.9422,0.9493];
% T1=[0.9783,0.9779];
% % (\beta^*_1, \beta^*_2)=(0.949, 0.9779)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% % 

% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.65,
% % sigma2/sigma1=0.5; 
% beta2=[0.98,0.9834];
% T2=[0.9604,0.9629];
% beta1=[0.9604,0.9629];
% T1=[0.9834,0.9833];
% % (\beta^*_1, \beta^*_2)=(0.9628, 0.9833)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% 

% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.7,
% % sigma2/sigma1=0.5; 
% beta2=[0.9874,0.989];
% T2=[0.973,0.974];
% beta1=[0.973,0.974];
% T1=[0.9874,0.9874];
% % (\beta^*_1, \beta^*_2)=(0.973, 0.9874)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.75,
% % sigma2/sigma1=0.5; 
% beta2=[0.989,0.9906];
% T2=[0.98,0.9808];
% beta1=[0.98,0.9808];
% T1=[0.9906,0.9906];
% % (\beta^*_1, \beta^*_2)=(0.9808, 0.9906)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% 
% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.8,
% % sigma2/sigma1=0.5; 
% beta2=[0.988,0.9932];
% T2=[0.9849,0.9868];
% beta1=[0.9849,0.9868];
% T1=[0.9932,0.9932];
% % (\beta^*_1, \beta^*_2)=(0.9868, 0.9932)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.85,
% % sigma2/sigma1=0.5; 
% beta2=[0.995,0.9953];
% T2=[0.9914,0.9915];
% beta1=[0.9914,0.9915];
% T1=[0.9953,0.9953];
% % (\beta^*_1, \beta^*_2)=(0.9915, 0.9953)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.9,
% % sigma2/sigma1=0.5; 
% beta2=[0.995,0.9971];
% T2=[0.9948,0.9952];
% beta1=[0.9948,0.9952];
% T1=[0.9971,0.9971];
% % (\beta^*_1, \beta^*_2)=(0.9952, 0.9971)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.04, lambda=0.99, rho=0.95,
% % sigma2/sigma1=0.5; 
% beta2=[0.998,0.9987];
% T2=[0.998,0.998];
% beta1=[0.998];
% T1=[0.9987];
% % (\beta^*_1, \beta^*_2)=(0.998, 0.9987)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)



% % % 
% % % 
% figure(1)
% phipi=[0.8:0.1:1.2, 1.4:0.2:2];0.9929, 0.9962
beta1=0:0.05:1;
rho=[0,0.1,0.2:0.05:0.3,0.32,0.35,0.38,0.4,0.45,0.5:0.05:0.95,0.9999];
beta1star=[0, 0.1088,  0.2428,  0.3256, 0.4284,0.4809, 0.5857,0.7023, 0.7609,0.8513, 0.8999,0.9294, 0.949, 0.9628, 0.973, 0.9808,  0.9868,0.9915,  0.9952,0.998,1 ];
beta2star=[0, 0.1118,  0.2625,   0.3679, 0.5179,0.601, 0.7529, 0.8622,0.898, 0.94,   0.9592, 0.9704, 0.9779,0.9833,0.9874,0.9906,0.9932,0.9953, 0.9971,0.9987,1];
% plot(rho, beta1star,'b.-',beta1,beta1,'r--')
% xlabel('\bf\rho') 
% ylabel('\bf\beta_1^*','Rotation',0)
% hold on ,
% % 
% figure(2)
% phipi=[0.8:0.1:1.2, 1.4:0.2:2];
% beta1=0:0.05:1;
% rho=[0:0.1:0.9, 0.999];
rhoq=0:0.01:1;
beta1starq=interp1(rho, beta1star, rhoq,'pchip');
beta2starq=interp1(rho, beta2star, rhoq,'pchip');


plot(rhoq, beta1starq,'b-',rhoq, beta2starq,'g-',beta1,beta1,'r--')
xlabel('\rho') 
ylabel('\beta_i^*','Rotation',0)
% legend('\beta_1^* for output gap','\beta_2^* for inflation')
% hold on 

% 


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Our own calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % phiy=0.5, phipi=1.5, varphi=1, gamma=0.075, lambda=0.99, rho=0, sigmapi^2/sigmay^2=1; 
% % 
% % % (\beta^*_1, \beta^*_2)=(0, 0)*****
% % plot(beta1,T1,'r-',T2,beta2,'r-')
% % % axis([0,1,0,1])
% % xlabel('\bf\beta_1') 
% % ylabel('\bf\beta_2','Rotation',0)
% 
% % %%phiy=0.5, phipi=1.5, varphi=1, gamma=0.075, lambda=0.99, rho=0.1, sigmapi^2/sigmay^2=1; 
% % beta2=[0.1, 0.12];
% % T2=[0.1093, 0.1103];
% % beta1=[0.1, 0.12];
% % T1=[0.1114,0.1112];
% % % (\beta^*_1, \beta^*_2)=(0.1099, 0.1113)*****
% % plot(beta1,T1,'r-',T2,beta2,'r-')
% % % axis([0,1,0,1])
% % xlabel('\bf\beta_1') 
% % ylabel('\bf\beta_2','Rotation',0)
% % % hold on
% % % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % % axis([0,1,0,1].)
% 
% % %%phiy=0.5, phipi=1.5, varphi=1, gamma=0.075, lambda=0.99, rho=0.4, sigmapi^2/sigmay^2=1; 
% % beta2=[0.7,0.82, 0.9];
% % T2=[0.6899,0.7626,0.8299];
% % beta1=[0.7, 0.775,0.84];
% % T1=[0.8464,0.8203,0.7896];
% % % (\beta^*_1, \beta^*_2)=(0.7656, 0.8236)*****
% % plot(beta1,T1,'r-',T2,beta2,'r-')
% % % axis([0,1,0,1])
% % xlabel('\bf\beta_1') 
% % ylabel('\bf\beta_2','Rotation',0)
% % % hold on
% % % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % % axis([0,1,0,1].)
% % 
% % %%phiy=0.5, phipi=1.5, varphi=1, gamma=0.075, lambda=0.99, rho=0.5, sigmapi^2/sigmay^2=1; 
% % beta2=[0.9,0.925 0.94];
% % T2=[0.8822,0.8985,0.9092];
% % beta1=[0.85, 0.9,0.902,0.91];
% % T1=[0.9372,0.9284,0.928,0.9263];
% % % (\beta^*_1, \beta^*_2)=(0.9008, 0.9282)*****
% % plot(beta1,T1,'r-',T2,beta2,'r-')
% % % axis([0,1,0,1])
% % xlabel('\bf\beta_1') 
% % ylabel('\bf\beta_2','Rotation',0)
% % % hold on
% % % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % % axis([0,1,0,1].)
% 
% % %%phiy=0.5, phipi=1.5, varphi=1, gamma=0.075, lambda=0.99, rho=0.6, sigmapi^2/sigmay^2=1; 
% % beta2=[0.95,0.962,0.99];
% % T2=[0.9421,0.9484,0.965];
% % beta1=[0.92, 0.948,0.96];
% % T1=[0.965,0.962,0.9605];
% % % (\beta^*_1, \beta^*_2)=(0.9484, 0.962)*****
% % plot(beta1,T1,'r-',T2,beta2,'r-')
% % % axis([0,1,0,1])
% % xlabel('\bf\beta_1') 
% % ylabel('\bf\beta_2','Rotation',0)
% % % hold on
% % % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % % axis([0,1,0,1].)
% % 
% 
% % %%phiy=0.5, phipi=1.5, varphi=1, gamma=0.075, lambda=0.99, rho=0.8, sigmapi^2/sigmay^2=1; 
% % beta2=[0.98,0.995];
% % T2=[0.9833,0.9872];
% % beta1=[0.96,0.99];
% % T1=[0.9895,0.9884];
% % % (\beta^*_1, \beta^*_2)=(0.9855, 0.9886)*****
% % plot(beta1,T1,'r-',T2,beta2,'r-')
% % % axis([0,1,0,1])
% % xlabel('\bf\beta_1') 
% % ylabel('\bf\beta_2','Rotation',0)
% % % hold on
% % % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % % axis([0,1,0,1].)
% 
% 
% % %%phiy=0.5, phipi=1.5, varphi=1, gamma=0.075, lambda=0.99, rho=0.999, sigmapi^2/sigmay^2=1; 
% % beta2=[0.99,0.9999];
% % T2=[0.9999,1];
% % beta1=[0.99,0.9999];
% % T1=[1,1];
% % % (\beta^*_1, \beta^*_2)=(1, 1)*****
% % plot(beta1,T1,'r-',T2,beta2,'r-')
% % % axis([0,1,0,1])
% % xlabel('\bf\beta_1') 
% % ylabel('\bf\beta_2','Rotation',0)
% % % hold on
% % % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % % axis([0,1,0,1].)
% % 
% 
% 
% % 
% % 
% % 
% % figure(1)
% % phipi=[0.8:0.1:1.2, 1.4:0.2:2];0.9929, 0.9962
% beta1=0:0.05:1;
% rho=[0,0.1,0.4,0.5,0.6,0.8,0.999];
% beta1star=[0,0.1099, 0.7656,0.9008, 0.9484,0.9855, 1];
% beta2star=[0,0.1113, 0.8236,0.9282, 0.962,0.9886,1];
% % plot(rho, beta1star,'b.-',beta1,beta1,'r--')
% % xlabel('\bf\rho') 
% % ylabel('\bf\beta_1^*','Rotation',0)
% % hold on 
% % % 
% % figure(2)
% % phipi=[0.8:0.1:1.2, 1.4:0.2:2];
% % beta1=0:0.05:1;
% % rho=[0:0.1:0.9, 0.999];
% 
% plot(rho, beta1star,'b.--',rho, beta2star,'g.-',beta1,beta1,'r--')
% xlabel('\rho') 
% ylabel('\beta_i^*','Rotation',0)
% legend('\beta_1^* for output gap','\beta_2^* for inflation')
% % hold on 




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Clarida et al's calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% phiy=0.5, phipi=1.5, varphi=1, gamma=0.3, lambda=0.99, rho=0, sigmapi^2/sigmay^2=0.1; 
% beta2=[0.95, 1];
% T2=[0.9932, 1];
% beta1=[0.95, 1];
% T1=[0.9998, 1];
% % (\beta^*_1, \beta^*_2)=(0, 0)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% hold on
% plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% %%phiy=0.5, phipi=1.5, varphi=1, gamma=0.3, lambda=0.99, rho=0.1, sigmapi^2/sigmay^2=0.1; 
% beta2=[0.1, 0.12];
% T2=[0.1056, 0.1056];
% beta1=[0.1, 0.12];
% T1=[0.1122, 0.1124];
% % (\beta^*_1, \beta^*_2)=(0.1056, 0.1123)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% %%phiy=0.5, phipi=1.5, varphi=1, gamma=0.3, lambda=0.99, rho=0.2, sigmapi^2/sigmay^2=0.1; 
% beta2=[0.2, 0.28];
% T2=[0.2244, 0.2236];
% beta1=[0.2, 0.24];
% T1=[0.2606, 0.2618];
% % (\beta^*_1, \beta^*_2)=(0.2238, 0.2613)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% 
% %%phiy=0.5, phipi=1.5, varphi=1, gamma=0.3, lambda=0.99, rho=0.3, sigmapi^2/sigmay^2=0.1; 
% beta2=[0.3, 0.5];
% T2=[0.3572, 0.3536];
% beta1=[0.3, 0.4];
% T1=[0.4898, 0.4955];
% % (\beta^*_1, \beta^*_2)=(0.3537, 0.4929)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% %%phiy=0.5, phipi=1.5, varphi=1, gamma=0.3, lambda=0.99, rho=0.4, sigmapi^2/sigmay^2=0.1; 
% beta2=[0.5, 0.8];
% T2=[0.4976, 0.4927];
% beta1=[0.4, 0.5];
% T1=[0.7931, 0.7937];
% % (\beta^*_1, \beta^*_2)=(0.4928, 0.7937)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% %%phiy=0.5, phipi=1.5, varphi=1, gamma=0.3, lambda=0.99, rho=0.5, sigmapi^2/sigmay^2=0.1; 
% beta2=[0.9, 0.92];
% T2=[0.6377, 0.6433];
% beta1=[0.6, 0.64];
% T1=[0.9076, 0.9066];
% % (\beta^*_1, \beta^*_2)=(0.6396, 0.9066)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
%

% %%phiy=0.5, phipi=1.5, varphi=1, gamma=0.3, lambda=0.99, rho=0.6, sigmapi^2/sigmay^2=0.1; 
% beta2=[0.94, 0.97];
% T2=[0.7572, 0.7729];
% beta1=[0.74, 0.77];
% T1=[0.9465, 0.9456];
% % (\beta^*_1, \beta^*_2)=(0.7603, 0.9459)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% 
% %%phiy=0.5, phipi=1.5, varphi=1, gamma=0.3, lambda=0.99, rho=0.7, sigmapi^2/sigmay^2=0.1; 
% beta2=[0.95, 0.98];
% T2=[0.8407, 0.8535];
% beta1=[0.8, 0.86];
% T1=[0.9676, 0.9661];
% % (\beta^*_1, \beta^*_2)=(0.8477, 0.9664)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% %%phiy=0.5, phipi=1.5, varphi=1, gamma=0.3, lambda=0.99, rho=0.8, sigmapi^2/sigmay^2=0.1; 
% beta2=[0.96, 0.99];
% T2=[0.9038, 0.9135];
% beta1=[0.9, 0.92];
% T1=[0.9803, 0.9799];
% % (\beta^*_1, \beta^*_2)=(0.9103, 0.9801)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% 
% %%phiy=0.5, phipi=1.5, varphi=1, gamma=0.3, lambda=0.99, rho=0.9, sigmapi^2/sigmay^2=0.1; 
% beta2=[0.98, 0.995];
% T2=[0.955, 0.9583];
% beta1=[0.94, 0.98];
% T1=[0.9909, 0.9903];
% % (\beta^*_1, \beta^*_2)=(0.9573, 0.9906)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% %%phiy=0.5, phipi=1.5, varphi=1, gamma=0.3, lambda=0.99, rho=0.999, sigmapi^2/sigmay^2=0.1; 
% beta2=[0.99, 0.9999];
% T2=[0.9995, 0.9995];
% beta1=[0.99, 0.9999];
% T1=[0.9999, 0.9999];
% % (\beta^*_1, \beta^*_2)=(0.9995, 0.9999)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% % 
% % 
% % 
% % figure(1)
% % phipi=[0.8:0.1:1.2, 1.4:0.2:2];0.9929, 0.9962
% beta1=0:0.05:1;
% rho=[0:0.1:0.9, 0.999];
% beta1star=[0,0.1056,0.2238, 0.3537,0.4928, 0.6396, 0.7603,0.8477, 0.9103,0.9573,0.9995];
% plot(rho, beta1star,'b.-',beta1,beta1,'r--')
% xlabel('\bf\rho') 
% ylabel('\bf\beta_1^*','Rotation',0)
% hold on 
% % 
% % figure(2)
% % phipi=[0.8:0.1:1.2, 1.4:0.2:2];
% beta1=0:0.05:1;
% rho=[0:0.1:0.9, 0.999];
% beta2star=[0, 0.1123, 0.2613,0.4929,0.7937,0.9066, 0.9459,0.9664,0.9801,0.9906, 0.9999];
% plot(rho, beta2star,'g.-',beta1,beta1,'r--')
% xlabel('\bf\rho') 
% ylabel('\bf\beta_i^*','Rotation',0)
% % hold on 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Woodford's calibration %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% phiy=0.5, phipi=1.5, varphi=6.3694, gamma=0.024, lambda=0.99, rho=0, sigmapi^2/sigmay^2=1; 
% beta2=[0, 1];
% T2=[0, 1];
% beta1=[0, 1];
% T1=[0, 1];
% % (\beta^*_1, \beta^*_2)=(0, 0)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% 
% %% phiy=0.5, phipi=1.5, varphi=6.3694, gamma=0.024, lambda=0.99, rho=0.1, sigmapi^2/sigmay^2=1; 
% beta2=[0.1, 0.12];
% T2=[0.1055, 0.107];
% beta1=[0.1, 0.12];
% T1=[0.1125, 0.1124];
% % (\beta^*_1, \beta^*_2)=(0.1064, 0.1125)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% %% phiy=0.5, phipi=1.5, varphi=6.3694, gamma=0.024, lambda=0.99, rho=0.2, sigmapi^2/sigmay^2=1; 
% beta2=[0.25, 0.27];
% T2=[0.2297, 0.2331];
% beta1=[0.2, 0.23, 0.24];
% T1=[0.2644, 0.2641, 0.264];
% % (\beta^*_1, \beta^*_2)=(0.2321, 0.2641)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% % 
% 
% %% phiy=0.5, phipi=1.5, varphi=6.3694, gamma=0.024, lambda=0.99, rho=0.3,
%  %% sigmapi^2/sigmay^2=1; 
% beta2=[0.4, 0.55];
% T2=[0.3723, 0.4213];
% beta1=[0.3, 0.4, 0.42];
% T1=[0.5273, 0.5242, 0.5234];
% % (\beta^*_1, \beta^*_2)=(0.4127, 0.5237)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% % %% phiy=0.5, phipi=1.5, varphi=6.3694, gamma=0.024, lambda=0.99, rho=0.4,
% %  %% sigmapi^2/sigmay^2=1; 
% beta2=[0.8, 0.92];
% T2=[0.6579, 0.7752];
% beta1=[0.5, 0.75, 0.78];
% T1=[0.9181, 0.9105, 0.9093];
% % (\beta^*_1, \beta^*_2)=(0.7653, 0.9099)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)

% % %% phiy=0.5, phipi=1.5, varphi=6.3694, gamma=0.024, lambda=0.99, rho=0.5,
%  %% sigmapi^2/sigmay^2=1; 
% beta2=[0.9, 0.98];
% T2=[0.8162, 0.9169];
% beta1=[0.8, 0.93];
% T1=[0.9693, 0.9676];
% % (\beta^*_1, \beta^*_2)=(0.9018, 0.968)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% %

% %%% phiy=0.5, phipi=1.5, varphi=6.3694, gamma=0.024, lambda=0.99, rho=0.6,
%  %% sigmapi^2/sigmay^2=1; 
% beta2=[0.96, 0.99];
% T2=[0.9217, 0.9572];
% beta1=[0.92, 0.95];
% T1=[0.984, 0.9837];
% % (\beta^*_1, \beta^*_2)=(0.9518, 0.9837) previous (0.9497, 0.9837)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% % %

% %%% phiy=0.5, phipi=1.5, varphi=6.3694, gamma=0.024, lambda=0.99, rho=0.7,
%  %% sigmapi^2/sigmay^2=1; 
% beta2=[0.98, 0.99,0.992];
% T2=[0.9653,0.9737,0.9755];
% beta1=[0.97, 0.98];
% T1=[0.9911,0.9911];
% % (\beta^*_1, \beta^*_2)=(0.9747,0.9911)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% %

% %%% phiy=0.5, phipi=1.5, varphi=6.3694, gamma=0.024, lambda=0.99, rho=0.8,
%  %% sigmapi^2/sigmay^2=1; 
% beta2=[0.98, 0.999];
% T2=[0.9814, 0.9907];
% beta1=[0.98, 0.99];
% T1=[0.9954, 0.9953];
% % (\beta^*_1, \beta^*_2)=(0.9889, 0.9953)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% % %

% %%% phiy=0.5, phipi=1.5, varphi=6.3694, gamma=0.024, lambda=0.99, rho=0.9,
%  %% sigmapi^2/sigmay^2=1; 
% beta2=[0.99,0.999];
% T2=[0.9948,0.9966];
% beta1=[0.99, 0.996,0.997];
% T1=[0.9981,0.9981,0.9981];
% % (\beta^*_1, \beta^*_2)=(0.9964,0.9981)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)


% %%% phiy=0.5, phipi=1.5, varphi=6.3694, gamma=0.024, lambda=0.99, rho=0.999,
%  %% sigmapi^2/sigmay^2=1; 
% beta2=[0.99, 0.999];
% T2=[1, 1];
% beta1=[0.99, 0.99];
% T1=[1, 1];
% % (\beta^*_1, \beta^*_2)=(1, 1)*****
% plot(beta1,T1,'r-',T2,beta2,'r-')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1].)
% % %


% 
% 
% % figure(1)
% % phipi=[0.8:0.1:1.2, 1.4:0.2:2];0.9929, 0.9962
% beta1=0:0.05:1;
% rho=[0:0.1:0.6,0.7,0.8, 0.9,0.999];
% beta1star=[0,0.1064,0.2321, 0.4127, 0.7653, 0.9018, 0.9518,0.9747,0.9889,0.9964,1 ];
% plot(rho, beta1star,'b.-',beta1,beta1,'r--')
% xlabel('\rho') 
% ylabel('\beta_1^*','Rotation',0)
% hold on 
% % 
% % figure(2)
% % phipi=[0.8:0.1:1.2, 1.4:0.2:2];
% beta1=0:0.05:1;
% rho=[0:0.1:0.6,0.7, 0.8,0.9, 0.999];
% beta2star=[0, 0.1125,0.2641, 0.5237,0.9099,0.968,0.9837,0.9911,0.9953,0.9981,1];
% plot(rho, beta2star,'g.-',beta1,beta1,'r--')
% xlabel('\rho') 
% ylabel('\beta_i^*','Rotation',0)
% % hold on 
% 
