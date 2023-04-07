var   
 mc zcap rk k q  c inve y lab pinf w r kp eps_a  eps_b eps_g eps_i  eps_r  eps_p eps_w  
labobs robs cpi_quarterly dy dc dinve dw ;    
 
varexo eta_a eta_b eta_g eta_i  eta_r  eta_p eta_w  ;  
 

parameters curvw curvp rho_ga delta phi_w l_bar pi_bar beta_const alpha gamma_bar mu_w mu_p  
psi  phi lambda phi_p iota_w xi_w iota_p  xi_p sigma_c sigma_l
r_pi r_dy r_y rho  G 

rho_a  rho_b rho_p  rho_w rho_i rho_r rho_g ; 

//cpie beta_bar cr crk cw cikbar cik clk cky ciy ccy crkky cwhlc cwly r_bar  gamma   beta 

// fixed parameters
delta=.025;
phi_w=1.5;
G=0.18;
curvp=10;
curvw=10;

// estimated parameters initialisation
alpha=.24;

sigma_c=1.5;

phi_p=1.5;
rho_ga=0.51;

phi= 6.0144;
lambda=    0.6361;    
xi_w=   0.8087;
sigma_l=    1.9423;
xi_p=   0.6;
iota_w=    0.3243;
iota_p=    0.47;
psi=    0.2696;
r_pi=     1.488;
rho=      0.8762;
r_y=      0.0593;
r_dy=     0.2347;

rho_a=    0.9977;
rho_b=    0.5799;
rho_g=    0.9957;
rho_i=   0.7165;

rho_r=0;
rho_p=0;
rho_w=0;
mu_p = 0;
mu_w  = 0;

// derived from steady state






model(linear); 




	

// sticky price - wage economy

# cpie=pi_bar/100+1;
# gamma=gamma_bar/100+1;
# beta = 1/(1 + beta_const/100);
# cr=cpie/(beta*gamma^(-sigma_c));
# r_bar=(cr-1)*100;
# beta_bar=beta*gamma^(-sigma_c);
# crk=(beta^(-1))*(gamma^sigma_c) - (1-delta);
# cw = (alpha^alpha*(1-alpha)^(1-alpha)/(phi_p*crk^alpha))^(1/(1-alpha));
# cikbar=(1-(1-delta)/gamma);
# cik=(1-(1-delta)/gamma)*gamma;
# clk=((1-alpha)/alpha)*(crk/cw);
# cky=phi_p*(clk)^(alpha-1);
# ciy=cik*cky;
# ccy=1-G-cik*cky;
# crkky=crk*cky;
# cwhlc=(1/phi_w)*(1-alpha)/alpha*crk*cky/ccy;
# cwly=1-crk*cky;



	      mc =  alpha*rk+(1-alpha)*(w) - 1*eps_a - 0*(1-alpha)*eps_a ;
	      zcap =  (1/(psi/(1-psi)))* rk ;
	      rk =  w+lab-k ;
	      k =  kp(-1)+zcap ;
	      inve = (1/(1+beta_bar*gamma))* (  inve(-1) + beta_bar*gamma*inve(1)+(1/(gamma^2*phi))*q ) +eps_i ;
          q = -r+pinf(1)-1*eps_b +0*(1/((1-lambda/gamma)/(sigma_c*(1+lambda/gamma))))*eps_b + (crk/(crk+(1-delta)))*rk(1) +  ((1-delta)/(crk+(1-delta)))*q(1) ;
	      c = (lambda/gamma)/(1+lambda/gamma)*c(-1) + (1/(1+lambda/gamma))*c(+1) +((sigma_c-1)*cwhlc/(sigma_c*(1+lambda/gamma)))*(lab-lab(+1)) - (1-lambda/gamma)/(sigma_c*(1+lambda/gamma))*(r-pinf(+1) + 1*eps_b) +0*eps_b ;
	      y = ccy*c+ciy*inve+eps_g  +  1*crkky*zcap ;
	      y = phi_p*( alpha*k+(1-alpha)*lab +eps_a );
	      pinf =  (1/(1+beta_bar*gamma*iota_p)) * ( beta_bar*gamma*pinf(1) +iota_p*pinf(-1) 
               +((1-xi_p)*(1-beta_bar*gamma*xi_p)/xi_p)/((phi_p-1)*curvp+1)*(mc)  )  + eps_p; 
	      w =  (1/(1+beta_bar*gamma))*w(-1)
               +(beta_bar*gamma/(1+beta_bar*gamma))*w(1)
               +(iota_w/(1+beta_bar*gamma))*pinf(-1)
               -(1+beta_bar*gamma*iota_w)/(1+beta_bar*gamma)*pinf
               +(beta_bar*gamma)/(1+beta_bar*gamma)*pinf(1)
               +(1-xi_w)*(1-beta_bar*gamma*xi_w)/((1+beta_bar*gamma)*xi_w)*(1/((phi_w-1)*curvw+1))*
               (sigma_l*lab + (1/(1-lambda/gamma))*c - ((lambda/gamma)/(1-lambda/gamma))*c(-1) -w) 
               + 1*eps_w ;
	      r =  r_pi*(1-rho)*pinf
               +r_y*(1-rho)*(y-phi_p*eps_a)     
               +r_dy*(y-phi_p*eps_a-y(-1)+phi_p*eps_a(-1))
               +rho*r(-1)
               +eps_r  ;

	      eps_a = rho_a*eps_a(-1)  + eta_a;
	      eps_b = rho_b*eps_b(-1) + eta_b;
	      eps_g = rho_g*(eps_g(-1)) + eta_g ;%+ rho_ga*eta_a;
	      eps_i = rho_i*eps_i(-1) + eta_i;
	      eps_r = rho_r*eps_r(-1) + eta_r;
	      eps_p = rho_p*eps_p(-1) + eta_p ;%- mu_p*eta_p(-1);
	      eps_w = rho_w*eps_w(-1) + eta_w ;%- mu_w*eta_w(-1) ;
	         
	      kp =  (1-cikbar)*kp(-1)+cikbar*inve + cikbar*gamma^2*phi*eps_i ;

// measurement equations

dy=y-y(-1)+gamma_bar;
dc=c-c(-1)+gamma_bar;
dinve=inve-inve(-1)+gamma_bar;
dw=w-w(-1)+gamma_bar;
cpi_quarterly = 1*(pinf) + pi_bar;
robs =    1*(r) + r_bar;
labobs = lab + l_bar;

end; 




shocks;
var eta_a;
stderr 0.4618;
var eta_b;
stderr 1.8513;
var eta_g;
stderr 0.6090;
var eta_i;
stderr 0.6017;
var eta_r;
stderr 0.2397;
var eta_p;
stderr 0.1455;
var eta_w;
stderr 0.2089;
end;



estimated_params;


phi,5.2,2,15,NORMAL_PDF,4,1.5;
sigma_c,1.36,0.25,3,NORMAL_PDF,1.50,0.375;
lambda,0.75,0.001,0.99,BETA_PDF,0.7,0.1;
xi_w,0.59,0.3,0.95,BETA_PDF,0.5,0.1;
sigma_l,0.5,0.25,10,NORMAL_PDF,2,0.75;
xi_p,0.7813,0.5,0.99,BETA_PDF,0.5,0.10;
iota_w,0.36,0.01,0.99,BETA_PDF,0.5,0.15;
iota_p,0.58,0.01,0.99,BETA_PDF,0.5,0.15;
psi,0.45,0.01,1,BETA_PDF,0.5,0.15;
phi_p,1.61,1.0,3,NORMAL_PDF,1.25,0.125;
r_pi,1.88,1.0,3,NORMAL_PDF,1.5,0.25;
rho,0.84,0.5,0.975,BETA_PDF,0.75,0.10;
r_y,0.11,0.001,0.5,NORMAL_PDF,0.125,0.05;
r_dy,0.16,0.001,0.5,NORMAL_PDF,0.125,0.05;
pi_bar,0.67,0.1,2.0,GAMMA_PDF,0.625,0.1;//20;
beta_const,0.16,0.01,2.0,GAMMA_PDF,0.25,0.1;//0.20;
l_bar,1.18,-10.0,10.0,NORMAL_PDF,0.0,2.0;
gamma_bar,0.41,0.1,0.8,NORMAL_PDF,0.4,0.10;

alpha,0.19,0.01,1.0,NORMAL_PDF,0.3,0.05;

rho_a,0.92 ,.01,.9999,BETA_PDF,0.5,0.20;
rho_b,0.21,.01,.9999,BETA_PDF,0.5,0.20;
rho_g,0.99,.01,.9999,BETA_PDF,0.5,0.20;
rho_i,0.77,.01,.9999,BETA_PDF,0.5,0.20;
rho_r,0.06,.01,.9999,BETA_PDF,0.5,0.20;
rho_p,0.09,.01,.9999,BETA_PDF,0.5,0.20;
rho_w,0.94,.001,.9999,BETA_PDF,0.5,0.20;
rho_ga,0.53,0.01,2.0,NORMAL_PDF,0.5,0.25;
//mu_p,.7652,0.01,.9999,BETA_PDF,0.5,0.2;
//mu_w,.8936,0.01,.9999,BETA_PDF,0.5,0.2;


stderr eta_a,0.4618,0.01,3,INV_GAMMA_PDF,0.1,2;
stderr eta_b,2.47,0.025,5,INV_GAMMA_PDF,0.1,2;
stderr eta_g,0.5,0.01,3,INV_GAMMA_PDF,0.1,2;
stderr eta_i,0.4,0.01,3,INV_GAMMA_PDF,0.1,2;
stderr eta_r,0.21,0.01,3,INV_GAMMA_PDF,0.1,2;
stderr eta_p,0.16,0.01,3,INV_GAMMA_PDF,0.1,2;
stderr eta_w,0.09,0.01,3,INV_GAMMA_PDF,0.1,2;

end;

varobs dy dc dinve labobs cpi_quarterly dw robs;


estimation(optim=('MaxIter',500),
mode_file=SW_Estimation_REE_mode,
datafile=slobodyan_dataset,
nograph,
mode_compute=0,
first_obs=71,
//nobs=141,
presample=4,
kalman_algo=1,
cova_compute=0,
lik_init=2,
prefilter=0,
mh_replic=00000,
mh_nblocks=1,
mh_jscale=0.3,
mh_drop=0.5,
irf=100);

//estimation(optim=('MaxIter',500),mode_file=SW_Estimation_REE_mode,datafile=slobodyan_dataset,nograph,mode_compute=1,first_obs=71,presample=4,kalman_algo=1,lik_init=2,prefilter=0,mh_replic=00000,mh_nblocks=1,mh_jscale=0.3,mh_drop=0.5,irf=100);//long training sample
//estimation(optim=('MaxIter',500),datafile=raf_dataset,nograph,mode_compute=1,first_obs=75,presample=4,kalman_algo=1,lik_init=2,prefilter=0,mh_replic=0000,mh_nblocks=1,mh_jscale=0.3,mh_drop=0.2,irf=100);//long training sample
//estimation(optim=('MaxIter',500),datafile=slobodyan_dataset,nograph,mode_compute=1,first_obs=141,presample=4,kalman_algo=2,lik_init=2,prefilter=0,mh_replic=500000,mh_nblocks=1,mh_jscale=0.3,mh_drop=0.5,irf=100);//long training sample

 stoch_simul(periods=50000,nograph);

    