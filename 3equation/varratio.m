clear all
% close all
hold off


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code is used to calculate the ratio of varances of output
% gap and inflation at BLE and REE in Figure 2(b).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stvepi=0.5;
stvey=1; 
gamma=0.04;
varphi=1;
phipi=1.5;
phiy=0.5;
lambda=0.99;

rhov=[0.35,0.38,0.4,0.45,0.5:0.05:0.95,0.9999];
%The followong beta^* are obtained from beta2star_rho.m%
beta1starv=[ 0.5857,0.7023, 0.7609,0.8513, 0.8999,0.9294, 0.949, 0.9628, 0.973, 0.9808,  0.9868,0.9915,  0.9952,0.998,1 ];
beta2starv=[0.7529, 0.8622,0.898, 0.94,   0.9592, 0.9704, 0.9779,0.9833,0.9874,0.9906,0.9932,0.9953, 0.9971,0.9987,1];

B=[1, varphi*(1-lambda*phipi); gamma, gamma*varphi+lambda*(1+varphi*phiy)]/(1+gamma*varphi*phipi+varphi*phiy);
C=[1, -varphi*phipi; gamma, 1+varphi*phiy]/(1+gamma*varphi*phipi+varphi*phiy);
Sigmae=[stvey^2,0;0,stvepi^2];
 
for i=1:15;
rho=rhov(i);
beta1=beta1starv(i);
beta2=beta2starv(i);
beta=[beta1, 0; 0, beta2]; 

Sigma1=1/(1-rho^2)*inv(eye(2)-rho*B)*C*Sigmae*(inv(eye(2)-rho*B)*C)';  % variance at REE

%%The following expresssions are based on the theoretical results, which are used to calculate the variances of the system at BLE%%
lambda1=(beta1^2+(gamma*varphi+lambda+lambda*varphi*phiy)*beta2^2+sqrt((beta1^2+(gamma*varphi+lambda+lambda*varphi*phiy)*beta2^2)^2-4*lambda*beta1^2*beta2^2*(1+gamma*varphi*phipi+varphi*phiy)))/(2*(1+gamma*varphi*phipi+varphi*phiy));
lambda2=(beta1^2+(gamma*varphi+lambda+lambda*varphi*phiy)*beta2^2-sqrt((beta1^2+(gamma*varphi+lambda+lambda*varphi*phiy)*beta2^2)^2-4*lambda*beta1^2*beta2^2*(1+gamma*varphi*phipi+varphi*phiy)))/(2*(1+gamma*varphi*phipi+varphi*phiy));

delta=varphi*(1-lambda*phipi)*beta2^2*(lambda2-lambda1);
mhat=delta*(rho-lambda1)*(rho-lambda2);
khat=1/((mhat*(1+gamma*varphi*phipi+varphi*phiy)*delta)^2);

t=(beta1^2+(gamma*varphi+lambda*(1+varphi*phiy))*beta2^2)/(1+gamma*varphi*phipi+varphi*phiy);
d=lambda*beta1^2*beta2^2/(1+gamma*varphi*phipi+varphi*phiy);
d1=rho+t;
d2=rho*d;
d3=(rho*t+d);

b1=1;
b2=lambda*beta2^2;
c1=varphi*phipi;
c2=varphi*beta2^2;
dd2=stvey^2*(b1^2+b2^2-2*b1*b2*d1+(b1^2+b2^2)*d3-d2*((b1^2+b2^2)*d1-2*b1*b2*d3+(b1^2+b2^2)*d2))+stvepi^2*(c1^2+c2^2-2*c1*c2*d1+(c1^2+c2^2)*d3-d2*((c1^2+c2^2)*d1-2*c1*c2*d3+(c1^2+c2^2)*d2));

b1=-gamma;
b2=0;
c1=1+varphi*phiy;
c2=beta1^2;
dd22=stvey^2*(b1^2+b2^2-2*b1*b2*d1+(b1^2+b2^2)*d3-d2*((b1^2+b2^2)*d1-2*b1*b2*d3+(b1^2+b2^2)*d2))+stvepi^2*(c1^2+c2^2-2*c1*c2*d1+(c1^2+c2^2)*d3-d2*((c1^2+c2^2)*d1-2*c1*c2*d3+(c1^2+c2^2)*d2));

vary=khat*(delta^4*(rho-lambda2)^2*(rho-lambda1)^2/((1-rho^2)*(1-rho*lambda1)*(1-lambda1^2)*(1-rho*lambda2)*(1-lambda2^2)*(1-lambda1*lambda2)))*dd2;
varpi=khat*(delta^4*(rho-lambda2)^2*(rho-lambda1)^2/((1-rho^2)*(1-rho*lambda1)*(1-lambda1^2)*(1-rho*lambda2)*(1-lambda2^2)*(1-lambda1*lambda2)))*dd22;

%%The following two expressions are used to calculate the ratio of variances at BLE and at REE%%
varratioy_BLEvsREE(i)=vary/Sigma1(1,1);
varratiopi_BLEvsREE(i)=varpi/Sigma1(2,2);
end

rhoq=0.35:0.01:0.9999;
Ratio_vary=interp1(rhov, varratioy_BLEvsREE, rhoq,'pchip');
Ratio_varpi=interp1(rhov, varratiopi_BLEvsREE, rhoq,'pchip');


plot(rhoq, Ratio_vary,'b-',rhoq, Ratio_varpi,'r-',rhoq,ones(length(rhoq)),'g--')
xlabel('\rho') 
ylabel('{\sigma^2_{i,BLE}}/{\sigma^2_{i,REE}}')
axis([0.35 1 0.5 4.5])







