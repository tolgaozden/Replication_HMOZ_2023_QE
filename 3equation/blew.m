close all 
clear all
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Our own calibration
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % phipi=1.5,phiy=0.5; lambda=0.99, rho=0.5, stvepi=0.5;stvey=1;
% gamma=0.04;varphi=1;
 % 
beta2=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9592,0.9999];
beta1star=[0.7073,0.7081,0.7107,0.7153,0.7221,0.7318,0.7456,0.7654,0.7955,0.8467,0.9,0.9576];
beta1=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95,0.9999];
beta2star=[0.9776,0.9775,0.9772,0.9767,0.9758,0.9746,0.9729,0.9702,0.9661,0.9592,0.9538,0.9466];
% % % (\beta^*_1, \beta^*_2)=(0.9518, 0.9837)
% % 
beta2q=0:0.01:1;
beta1starq=interp1(beta2, beta1star, beta2q,'pchip');

beta1q=0:0.01:1;
beta2starq=interp1(beta1, beta2star, beta1q,'pchip');

plot(beta1q,beta2starq,'b-', beta1starq,beta2q,'b-',0.9, 0.9592,'b.', 0.5, 0.5,'r.')
% axis([0,1,0,1])
xlabel('\beta_y') 
ylabel('\beta_\pi','Rotation',0)




%%%%%%%%%%%%%%%%%%%%%%%%% Previous calibration %%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% % % phipi=1.5,phiy=0.5; lambda=0.99, rho=0.9, stvepi=1;stvey=1;
% % gamma=0.024;varphi=1/0.157;
%  % 
% beta2=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95,0.999];
% beta1star=[0.6596,0.6615,0.6674,0.6773,0.6916,0.7109,0.736,0.7685,0.8112,0.8704,0.9117,0.97];
% beta1=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9,0.95, 0.999];
% beta2star=[0.9864,0.9864,0.9863,0.9862,0.986,0.9858,0.9855,0.9851,0.9847,0.9841,0.9837,0.9834];
% % % % (\beta^*_1, \beta^*_2)=(0.9518, 0.9837)
% % % 
% plot(beta1,beta2star,'b-', beta1star,beta2,'b--',0.9518, 0.9837,'r.')
% % axis([0,1,0,1])
% xlabel('\bf\beta_1') 
% ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'r-', 0.94635, 0.99128, 'r.')
% % % axis([0,1,0,1].)

% hold on 

% % phipi=1.2, lambda=0.99, rho=0.8, n2=0.8; 
% beta1=[0.9, 0.95];
% T1=[0.97795, 0.9795];
% beta2=[0.95, 1];
% T2=[0.9267, 0.9214];
% % % (\beta^*_1, \beta^*_2)=(0.9005, 0.958)
% 
% plot(beta1,T1,'r-',T2,beta2,'r-', 0.9237, 0.9787, 'b.')


% hold on 
% 
% % phipi=1.2, lambda=0.99, rho=0.8, n2=0.6; 
% beta1=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
% T1=[0.9278, 0.9283, 0.9298, 0.9321, 0.9351, 0.9388, 0.9431, 0.9476, 0.9527, 0.958, 0.9636];
% beta2=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
% T2=[0.9096, 0.9095, 0.9094, 0.9091, 0.9088, 0.9082, 0.9075, 0.9064, 0.9049, 0.9026, 0.899];
% % % (\beta^*_1, \beta^*_2)=(0.9005, 0.958)
% % 
% plot(beta1,T1,'c-',T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % axis([0,1,0,1])
% % xlabel('\bf\beta_1') 
% % ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'c-', 0.9005, 0.958, 'c.')
% % % axis([0,1,0,1].)
% 
% hold on 

% % phipi=1.2, lambda=0.99, rho=0.8, n2=0.5; 
% beta1=[0.85, 0.9, 0.92];
% T1=[0.9411, 0.9443, 0.9446];
% beta2=[0.9, 0.94, 0.95];
% T2=[0.8869, 0.88607, 0.8858];
% % % (\beta^*_1, \beta^*_2)=(0.9005, 0.958)
% 
% plot(beta1,T1,'r-',T2,beta2,'r-', 0.88596, 0.9434, 'b.')

% hold on
% 
% % phipi=1.2, lambda=0.99, rho=0.8, n2=0.3; 
% beta1=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
% T1=[0.8664, 0.867, 0.8686, 0.8714, 0.8751, 0.8799, 0.8852, 0.8913, 0.898, 0.9051, 0.9125];
% beta2=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
% T2=[0.8561, 0.8561, 0.85603, 0.8559, 0.8557, 0.8554, 0.855, 0.8546, 0.854, 0.8532, 0.8522];
% % % (\beta^*_1, \beta^*_2)=(0.8531, 0.9019)
% % 
% plot(beta1,T1,'b-', T2,beta2,'b-', 0.8531, 0.9019, 'b.')
% % % axis([0,1,0,1])
% % xlabel('\bf\beta_1') 
% % ylabel('\bf\beta_2','Rotation',0)
% % hold on
% % plot(T2,beta2,'b-', 0.8531, 0.9019, 'b.')
% % % axis([0,1,0,1].)
% hold on

% % phipi=1.2, lambda=0.99, rho=0.8, n2=0.1; 
% beta1=[0.8, 0.82, 0.85];
% T1=[0.8382, 0.839, 0.8402];
% beta2=[0.8,0.85];
% T2=[0.81805, 0.81795];
% % % (\beta^*_1, \beta^*_2)=(0.9005, 0.958)
% 
% plot(beta1,T1,'r-',T2,beta2,'r-', 0.818, 0.8389, 'b.')

% hold on

% % phipi=1.2, lambda=0.99, rho=0.8, n2=0.01; 
% beta1=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
% T1=[0.80221, 0.80225, 0.80232, 0.80248, 0.80269, 0.80294, 0.80325, 0.80362, 0.80405, 0.80454, 0.80509];
% beta2=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1];
% T2=[0.80187, 0.80187, 0.80186, 0.80185, 0.80184, 0.80183, 0.80182, 0.80181, 0.8018, 0.80179, 0.80176];
% % % (\beta^*_1, \beta^*_2)=(0.8018, 0.80405)
% 
% 
% plot(beta1,T1,'m-', T2,beta2,'m-', , 0.80405, 'm.')

% hold off
% hold off
% hold off
% 
% 
% % 
% figure(2)
% beta1=[0.8018,0.818,0.8531, 0.88596, 0.9005,  0.9237,0.94635]
% beta2=[0.80405,0.8389, 0.9019, 0.9434, 0.958,  0.9787, 0.99128];
% plot(beta1, beta2, 'b-', 0.8:0.01:1, 0.8:0.01:1, 'r-.')
% xlabel('\beta_1^*') 
% ylabel('\beta_2^*','Rotation',0)
