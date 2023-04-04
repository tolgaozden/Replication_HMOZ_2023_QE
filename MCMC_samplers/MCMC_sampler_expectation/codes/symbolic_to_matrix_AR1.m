clear;clc;close all;

syms  G curvw curvp rho_ga delta phi_w l_bar pi_bar beta_const alpha gamma_bar mu_w mu_p...  
psi  phi lambda phi_p iota_w xi_w iota_p  xi_p sigma_c sigma_l...
r_pi r_dy r_y rho beta ...
rho_a  rho_b rho_p  rho_w rho_i rho_r rho_ga rho_g mu_w mu_p  ...
cpie beta_bar cr crk cw cikbar cik clk cky ciy ccy crkky cwhlc cwly...
r_bar  gamma ...
alpha_rk alpha_q alpha_c alpha_i alpha_l alpha_pi alpha_w ... 
beta_rk beta_q beta_c beta_i beta_l beta_pi beta_w  ;


syms  mc zcap rk k1 q       c inve y lab pinf w r kp eps_a  eps_b eps_g eps_i  eps_r  eps_p eps_w  labobs robs pinfobs dy dc dinve dw ...  
     mcm zcapm rkm k1m qm    cm invem ym labm pinfm wm rm kpm eps_am  eps_bm eps_gm eps_im  eps_rm  eps_pm eps_wm dym dcm dinvem dwm... 
     mcp zcapp rkp k1p qp   cp invep yp labp pinfp wp rp kpp dyp dcp dinvep dwp ... 
     eta_a eta_b eta_g eta_i eta_r eta_p eta_w eta_pi_exp...
     eta_am eta_bm eta_gm eta_im eta_rm eta_pm eta_wm eta_pi_expm... 
     dy_obs dc_obs dinve_obs dw_obs pinfobs robs labobs pinfobs_exp...
     rkp_m qp_m cp_m invep_m labp_m  pinfp_m wp_m;


beta=0.9995;
cpie=1+pi_bar/100;
gamma=1+gamma_bar/100;
beta_bar=beta*gamma^(-sigma_c);
cr=cpie/(beta*gamma^(-sigma_c));
crk=(beta^(-1))*(gamma^sigma_c) - (1-delta);
cw = (alpha^alpha*(1-alpha)^(1-alpha)/(phi_p*crk^alpha))^(1/(1-alpha));
cikbar=(1-(1-delta)/gamma);
cik=(1-(1-delta)/gamma)*gamma;
clk=((1-alpha)/alpha)*(crk/cw);
cky=phi_p*(clk)^(alpha-1);
ciy=cik*cky;
ccy=1-G-cik*cky;
crkky=crk*cky;
cwhlc=(1/phi_w)*(1-alpha)/alpha*crk*cky/ccy;
cwly=1-crk*cky;
r_bar=(cr-1)*100;



fs1=    -mc +  alpha*rk+(1-alpha)*(w) - 1*eps_a - 0*(1-alpha)*eps_a==0 ;
fs2=	-zcap +  (1/(psi/(1-psi)))* rk ==0;
fs3=	-rk +  w+lab-k1 ==0;
fs4=	-k1 +  kpm+zcap ==0;
fs5=	-inve + (1/(1+beta_bar*gamma))* (  invem + beta_bar*gamma*invep+(1/(gamma^2*phi))*q ) +eps_i ==0;
fs6=    -q + -r+pinfp-1*eps_b +0*(1/((1-lambda/gamma)/(sigma_c*(1+lambda/gamma))))*eps_b + (crk/(crk+(1-delta)))*rkp +  ((1-delta)/(crk+(1-delta)))*qp ==0;
fs7=	-c + (lambda/gamma)/(1+lambda/gamma)*cm + (1/(1+lambda/gamma))*cp +((sigma_c-1)*cwhlc/(sigma_c*(1+lambda/gamma)))*(lab-labp) - (1-lambda/gamma)/(sigma_c*(1+lambda/gamma))*(r-pinfp + 1*eps_b) +0*eps_b ==0;
fs8=	-y + ccy*c+ciy*inve+eps_g  +  1*crkky*zcap ==0 ;
fs9=	-y + phi_p*( alpha*k1+(1-alpha)*lab +eps_a ) ==0;
fs10=	-pinf +  (1/(1+beta_bar*gamma*iota_p)) * ( beta_bar*gamma*pinfp +iota_p*pinfm... 
               +((1-xi_p)*(1-beta_bar*gamma*xi_p)/xi_p)/((phi_p-1)*curvp+1)*(mc)  )  + eps_p ==0; 
fs11=	-w +  (1/(1+beta_bar*gamma))*wm...
               +(beta_bar*gamma/(1+beta_bar*gamma))*wp...
               +(iota_w/(1+beta_bar*gamma))*pinfm...
               -(1+beta_bar*gamma*iota_w)/(1+beta_bar*gamma)*pinf...
               +(beta_bar*gamma)/(1+beta_bar*gamma)*pinfp...
               +(1-xi_w)*(1-beta_bar*gamma*xi_w)/((1+beta_bar*gamma)*xi_w)*(1/((phi_w-1)*curvw+1))*...
               (sigma_l*lab + (1/(1-lambda/gamma))*c - ((lambda/gamma)/(1-lambda/gamma))*cm -w) ...
               + 1*eps_w ==0;
fs12=	-r +  r_pi*(1-rho)*pinf...
               +r_y*(1-rho)*(y-phi_p*eps_a)     ...
               +r_dy*(y-phi_p*eps_a-ym+phi_p*eps_am)...
               +rho*rm...
               +eps_r ==0 ;
	         
fs13=	-kp +  (1-cikbar)*kpm+cikbar*inve + cikbar*gamma^2*phi*eps_i ==0;

fs14=   -dy + y - ym ==0;
fs15=   -dc + c - cm ==0;
fs16=   -dw + w - wm ==0;
fs17=   -dinve+inve-invem==0;


err1=   -eps_a + rho_a*eps_am  + eta_a==0;
err2=   -eps_b + rho_b*eps_bm + eta_b==0;
err3=   -eps_g + rho_g*(eps_gm) + eta_g+ rho_ga*eta_a==0;
err4=   -eps_i + rho_i*eps_im + eta_i==0;
err5=   -eps_r + rho_r*eps_rm + eta_r==0;
err6=   -eps_p + rho_p*eps_pm + eta_p ==0;%- mu_p*eta_pm==0;
err7=   -eps_w + rho_w*eps_wm + eta_w ==0;%- mu_w*eta_wm==0 ;          


meas1= -dy_obs + dy + gamma_bar==0;
meas2= -dc_obs + dc + gamma_bar==0;
meas3= -dinve_obs + dinve + gamma_bar==0;
meas4= -dw_obs  + dw + gamma_bar==0;
meas5= -pinfobs + pinf + pi_bar==0;
meas6= -robs + r + r_bar  ==0;
meas7= -labobs +lab + l_bar ==0;
meas8= -pinfobs_exp + pinfp + pi_bar ==0;

exp1 = -rkp + (alpha_rk + beta_rk*alpha_rk) + beta_rk^2 * rkm == 0;
exp2 = -qp  + (alpha_q + beta_q*alpha_q) + beta_q^2 * qm == 0;
exp3 = -cp   + (alpha_c + beta_c*alpha_c) + beta_c^2 * cm == 0;
exp4 = -invep   + (alpha_i + beta_i*alpha_i) + beta_i^2 * invem == 0;
exp5 = -labp  + (alpha_l + beta_l*alpha_l) + beta_l^2 * labm == 0;
exp6 = -pinfp  + (alpha_pi + beta_pi*alpha_pi) + beta_pi^2 * pinfm  + eta_pi_exp == 0;
exp7 = -wp + (alpha_w + beta_w*alpha_w) + beta_w^2 * wm == 0;





equations_endo = [fs1,fs2,fs3,fs4,fs5,fs6,fs7,fs8,fs9,fs10,fs11,fs12,fs13,fs14,fs15,fs16,fs17,...
    err1,err2,err3,err4,err5,err6,err7,...
    exp1 exp2 exp3 exp4 exp5 exp6 exp7];


meas_equations = [meas1 meas2 meas3 meas4 meas5 meas6 meas7 meas8];

%17 endogenous, 7 exogenous, 7 expectations = 31 variables
Contemp_endo = [mc,zcap,rk,k1,q,c,inve,y,lab,pinf,w,r,kp,dy,dc,dinve,dw, eps_a,eps_b,eps_g,eps_i,eps_r,eps_p,eps_w,rkp,qp,cp,invep,labp, pinfp,wp]; ...
    
observables = [dy_obs dc_obs dinve_obs dw_obs pinfobs robs labobs pinfobs_exp];

%17 endogenous variables and 7 endogenous, 7 expectations
Lagged_endo=[mcm zcapm rkm k1m   qm cm invem ym labm pinfm wm rm kpm dym,dcm,dwm,dinvem,eps_am  eps_bm eps_gm eps_im  eps_rm  eps_pm eps_wm rkp_m qp_m cp_m invep_m labp_m  pinfp_m wp_m;];...

%7 iid shocks
Shocks = [eta_a ,eta_b ,eta_g ,eta_i ,eta_r ,eta_p ,eta_w,eta_pi_exp];
     
AA = equationsToMatrix(equations_endo, Contemp_endo);   
AA=-AA;    
BB = equationsToMatrix(equations_endo, Lagged_endo);
DD = equationsToMatrix(equations_endo, Shocks);

CC_tmp = [ (alpha_rk + beta_rk*alpha_rk) ;...
                (alpha_q + beta_q*alpha_q); ...
                (alpha_c + beta_c*alpha_c); ...
                (alpha_i + beta_i*alpha_i); ...
                (alpha_l + beta_l*alpha_l); ...
                (alpha_pi + beta_pi*alpha_pi); ...
                (alpha_w + beta_w*alpha_w)];

CC = [zeros(24,1); CC_tmp];


E =[gamma_bar;gamma_bar;gamma_bar;gamma_bar;pi_bar;r_bar;l_bar;pi_bar];
F = equationsToMatrix(meas_equations,[Contemp_endo]);

%  F=[ 0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0;
%      0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0;
%      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0;
%      0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0;
%      0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0;
%      0,0,0,0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0;
%      0,0,0,0,0,0,0,0,1,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0];
%  
 
fid = fopen('sysmat_SW_model_AR1.m','wt');
fprintf(fid,'%s\n','AA=...');
fprintf(fid,'%s',char(AA), ';');
fprintf(fid,'\n');
fprintf(fid,'%s\n','BB=...');
fprintf(fid,'%s',char(BB), ';');
fprintf(fid,'\n');
fprintf(fid,'%s\n','CC=...');
fprintf(fid,'%s',char(CC), ';');
fprintf(fid,'\n');
fprintf(fid,'%s\n','DD=...')
fprintf(fid,'%s',char(DD), ';')
fprintf(fid,'\n')
fprintf(fid,'%s\n','E=...')
fprintf(fid,'%s',char(E), ';')
fprintf(fid,'\n')
fprintf(fid,'%s\n','F=...')
fprintf(fid,'%s',char(F), ';');
fclose(fid);
 
 
 
% clearvars -except AA BB CC DD EE RHO FF GG E F;

