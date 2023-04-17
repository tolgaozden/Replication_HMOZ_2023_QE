close all 
clear all
clf


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%This code is used to calculate Figure 1.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

stvepi=0.5;
stvey=1; 
gamma=0.04;
varphi=1;
phipi=1.5;
phiy=0.5;
lambda=0.99;
rho=0.5;
beta2v=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.9592,0.9999];
beta1v=[0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.95,0.9999];  

for i=1:12
beta2=beta2v(i);
syms beta1

%%The following expresssions are based on the theoretical results, which are shown in Appendix A.5%%
t=(beta1^2+(gamma*varphi+lambda*(1+varphi*phiy))*beta2^2)/(1+gamma*varphi*phipi+varphi*phiy);
d=lambda*beta1^2*beta2^2/(1+gamma*varphi*phipi+varphi*phiy);
d1=rho+t;
d2=rho*d;
d3=(rho*t+d);

b1=1;
b2=lambda*beta2^2;
c1=varphi*phipi;
c2=varphi*beta2^2;
dd1=stvey^2*((b1*d1-b2)*(b1-b2*d1)+(b2*d3-b1*d2)*(b1*d3-b2*d2))+stvepi^2*((c1*d1-c2)*(c1-c2*d1)+(c2*d3-c1*d2)*(c1*d3-c2*d2));
dd2=stvey^2*(b1^2+b2^2-2*b1*b2*d1+(b1^2+b2^2)*d3-d2*((b1^2+b2^2)*d1-2*b1*b2*d3+(b1^2+b2^2)*d2))+stvepi^2*(c1^2+c2^2-2*c1*c2*d1+(c1^2+c2^2)*d3-d2*((c1^2+c2^2)*d1-2*c1*c2*d3+(c1^2+c2^2)*d2));
corrytheo=dd1/dd2;

beta1star(i)=single(vpasolve(corrytheo==beta1,beta1,[0 1]));
end

for i=1:12
beta1=beta1v(i);
syms beta2

%%The following expresssions are based on the theoretical results, which are shown in Appendix A.5%%
t=(beta1^2+(gamma*varphi+lambda*(1+varphi*phiy))*beta2^2)/(1+gamma*varphi*phipi+varphi*phiy);
d=lambda*beta1^2*beta2^2/(1+gamma*varphi*phipi+varphi*phiy);
d1=rho+t;
d2=rho*d;
d3=(rho*t+d);
b1=-gamma;
b2=0;
c1=1+varphi*phiy;
c2=beta1^2;
dd11=stvey^2*((b1*d1-b2)*(b1-b2*d1)+(b2*d3-b1*d2)*(b1*d3-b2*d2))+stvepi^2*((c1*d1-c2)*(c1-c2*d1)+(c2*d3-c1*d2)*(c1*d3-c2*d2));
dd22=stvey^2*(b1^2+b2^2-2*b1*b2*d1+(b1^2+b2^2)*d3-d2*((b1^2+b2^2)*d1-2*b1*b2*d3+(b1^2+b2^2)*d2))+stvepi^2*(c1^2+c2^2-2*c1*c2*d1+(c1^2+c2^2)*d3-d2*((c1^2+c2^2)*d1-2*c1*c2*d3+(c1^2+c2^2)*d2));
corrpitheo=dd11/dd22;

beta2star(i)=single(vpasolve(corrpitheo==beta2,beta2,[0 1]));
end

 
beta2q=0:0.01:1;
beta1starq=interp1(beta2v, beta1star, beta2q,'pchip');

beta1q=0:0.01:1;
beta2starq=interp1(beta1v, beta2star, beta1q,'pchip');

plot(beta1q,beta2starq,'b-', beta1starq,beta2q,'b-',0.9, 0.9592,'b.', 0.5, 0.5,'r.')
% axis([0,1,0,1])
xlabel('\beta_y') 
ylabel('\beta_\pi','Rotation',0)




