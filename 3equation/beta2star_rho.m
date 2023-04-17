% close all 
clear all
clf

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code is used to calculate Figure 2a.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



beta1=0:0.05:1;
rho=[0,0.1,0.2:0.05:0.3,0.32,0.35,0.38,0.4,0.45,0.5:0.05:0.95,0.9999];
%The following values for beta1star and beta2star are obtained with the same method as Figure 1, i.e. by running the code blew.m with different rho %
beta1star=[0, 0.1088,  0.2428,  0.3256, 0.4284,0.4809, 0.5857,0.7023, 0.7609,0.8513, 0.8999,0.9294, 0.949, 0.9628, 0.973, 0.9808,  0.9868,0.9915,  0.9952,0.998,1 ];
beta2star=[0, 0.1118,  0.2625,   0.3679, 0.5179,0.601, 0.7529, 0.8622,0.898, 0.94,   0.9592, 0.9704, 0.9779,0.9833,0.9874,0.9906,0.9932,0.9953, 0.9971,0.9987,1];

rhoq=0:0.01:1;
beta1starq=interp1(rho, beta1star, rhoq,'pchip');
beta2starq=interp1(rho, beta2star, rhoq,'pchip');


plot(rhoq, beta1starq,'b-',rhoq, beta2starq,'g-',beta1,beta1,'r--')
xlabel('\rho') 
ylabel('\beta_i^*','Rotation',0)




