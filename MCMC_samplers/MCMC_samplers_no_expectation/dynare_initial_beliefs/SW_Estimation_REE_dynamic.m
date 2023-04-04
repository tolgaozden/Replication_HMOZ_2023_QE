function [residual, g1, g2, g3] = SW_Estimation_REE_dynamic(y, x, params, steady_state, it_)
%
% Status : Computes dynamic model for Dynare
%
% Inputs :
%   y         [#dynamic variables by 1] double    vector of endogenous variables in the order stored
%                                                 in M_.lead_lag_incidence; see the Manual
%   x         [nperiods by M_.exo_nbr] double     matrix of exogenous variables (in declaration order)
%                                                 for all simulation periods
%   steady_state  [M_.endo_nbr by 1] double       vector of steady state values
%   params    [M_.param_nbr by 1] double          vector of parameter values in declaration order
%   it_       scalar double                       time period for exogenous variables for which to evaluate the model
%
% Outputs:
%   residual  [M_.endo_nbr by 1] double    vector of residuals of the dynamic model equations in order of 
%                                          declaration of the equations.
%                                          Dynare may prepend auxiliary equations, see M_.aux_vars
%   g1        [M_.endo_nbr by #dynamic variables] double    Jacobian matrix of the dynamic model equations;
%                                                           rows: equations in order of declaration
%                                                           columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g2        [M_.endo_nbr by (#dynamic variables)^2] double   Hessian matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%   g3        [M_.endo_nbr by (#dynamic variables)^3] double   Third order derivative matrix of the dynamic model equations;
%                                                              rows: equations in order of declaration
%                                                              columns: variables in order stored in M_.lead_lag_incidence followed by the ones in M_.exo_names
%
%
% Warning : this file is generated automatically by Dynare
%           from model file (.mod)

%
% Model equations
%

residual = zeros(29, 1);
cpie__ = 1+params(7)/100;
gamma__ = 1+params(10)/100;
beta__ = 1/(1+params(8)/100);
cr__ = cpie__/(beta__*gamma__^(-params(21)));
r_bar__ = 100*(cr__-1);
beta_bar__ = beta__*gamma__^(-params(21));
crk__ = beta__^(-1)*gamma__^params(21)-(1-params(4));
cw__ = (params(9)^params(9)*(1-params(9))^(1-params(9))/(params(16)*crk__^params(9)))^(1/(1-params(9)));
cikbar__ = 1-(1-params(4))/gamma__;
cik__ = gamma__*(1-(1-params(4))/gamma__);
clk__ = (1-params(9))/params(9)*crk__/cw__;
cky__ = params(16)*clk__^(params(9)-1);
ciy__ = cik__*cky__;
ccy__ = 1-params(27)-cik__*cky__;
crkky__ = crk__*cky__;
cwhlc__ = cky__*crk__*(1-params(9))*1/params(5)/params(9)/ccy__;
T88 = 1/(params(13)/(1-params(13)));
T103 = 1/(1+gamma__*beta_bar__);
T108 = gamma__^2;
T126 = params(15)/gamma__;
T130 = (1-T126)/(params(21)*(1+T126));
T153 = (params(21)-1)*cwhlc__/(params(21)*(1+T126));
T184 = 1/(1+gamma__*beta_bar__*params(19));
T199 = (1-params(20))*(1-gamma__*beta_bar__*params(20))/params(20)/(1+(params(16)-1)*params(2));
T208 = gamma__*beta_bar__/(1+gamma__*beta_bar__);
T235 = (1-params(18))*(1-gamma__*beta_bar__*params(18))/((1+gamma__*beta_bar__)*params(18))*1/(1+(params(5)-1)*params(1));
lhs =y(17);
rhs =params(9)*y(19)+(1-params(9))*y(27)-y(30);
residual(1)= lhs-rhs;
lhs =y(18);
rhs =y(19)*T88;
residual(2)= lhs-rhs;
lhs =y(19);
rhs =y(27)+y(25)-y(20);
residual(3)= lhs-rhs;
lhs =y(20);
rhs =y(18)+y(7);
residual(4)= lhs-rhs;
lhs =y(23);
rhs =T103*(y(2)+gamma__*beta_bar__*y(49)+1/(T108*params(14))*y(21))+y(33);
residual(5)= lhs-rhs;
lhs =y(21);
rhs =(-y(28))+y(51)-y(31)+crk__/(1-params(4)+crk__)*y(46)+(1-params(4))/(1-params(4)+crk__)*y(47);
residual(6)= lhs-rhs;
lhs =y(22);
rhs =T126/(1+T126)*y(1)+1/(1+T126)*y(48)+T153*(y(25)-y(50))-T130*(y(31)+y(28)-y(51));
residual(7)= lhs-rhs;
lhs =y(24);
rhs =ccy__*y(22)+y(23)*ciy__+y(32)+y(18)*crkky__;
residual(8)= lhs-rhs;
lhs =y(24);
rhs =params(16)*(y(30)+params(9)*y(20)+(1-params(9))*y(25));
residual(9)= lhs-rhs;
lhs =y(26);
rhs =T184*(gamma__*beta_bar__*y(51)+params(19)*y(4)+y(17)*T199)+y(35);
residual(10)= lhs-rhs;
lhs =y(27);
rhs =T103*y(5)+T208*y(52)+y(4)*params(17)/(1+gamma__*beta_bar__)-y(26)*(1+gamma__*beta_bar__*params(17))/(1+gamma__*beta_bar__)+y(51)*T208+T235*(y(25)*params(22)+y(22)*1/(1-T126)-y(1)*T126/(1-T126)-y(27))+y(36);
residual(11)= lhs-rhs;
lhs =y(28);
rhs =y(26)*params(23)*(1-params(26))+(1-params(26))*params(25)*(y(24)-params(16)*y(30))+params(24)*(y(24)-params(16)*y(30)-y(3)+params(16)*y(8))+params(26)*y(6)+y(34);
residual(12)= lhs-rhs;
lhs =y(30);
rhs =y(8)*params(28)+x(it_, 1);
residual(13)= lhs-rhs;
lhs =y(31);
rhs =params(29)*y(9)+x(it_, 2);
residual(14)= lhs-rhs;
lhs =y(32);
rhs =params(34)*y(10)+x(it_, 3)+x(it_, 1)*params(3);
residual(15)= lhs-rhs;
lhs =y(33);
rhs =params(32)*y(11)+x(it_, 4);
residual(16)= lhs-rhs;
lhs =y(34);
rhs =params(33)*y(12)+x(it_, 5);
residual(17)= lhs-rhs;
lhs =y(35);
rhs =params(30)*y(13)+x(it_, 6)-params(12)*y(15);
residual(18)= lhs-rhs;
lhs =y(36);
rhs =params(31)*y(14)+x(it_, 7)-params(11)*y(16);
residual(19)= lhs-rhs;
lhs =y(29);
rhs =y(7)*(1-cikbar__)+y(23)*cikbar__+y(33)*params(14)*T108*cikbar__;
residual(20)= lhs-rhs;
lhs =y(40);
rhs =params(10)+y(24)-y(3);
residual(21)= lhs-rhs;
lhs =y(41);
rhs =params(10)+y(22)-y(1);
residual(22)= lhs-rhs;
lhs =y(42);
rhs =params(10)+y(23)-y(2);
residual(23)= lhs-rhs;
lhs =y(43);
rhs =params(10)+y(27)-y(5);
residual(24)= lhs-rhs;
lhs =y(39);
rhs =params(7)+y(26);
residual(25)= lhs-rhs;
lhs =y(38);
rhs =y(28)+r_bar__;
residual(26)= lhs-rhs;
lhs =y(37);
rhs =y(25)+params(6);
residual(27)= lhs-rhs;
lhs =y(44);
rhs =x(it_, 6);
residual(28)= lhs-rhs;
lhs =y(45);
rhs =x(it_, 7);
residual(29)= lhs-rhs;
if nargout >= 2,
  g1 = zeros(29, 59);

  %
  % Jacobian matrix
  %

  g1(1,17)=1;
  g1(1,19)=(-params(9));
  g1(1,27)=(-(1-params(9)));
  g1(1,30)=1;
  g1(2,18)=1;
  g1(2,19)=(-T88);
  g1(3,19)=1;
  g1(3,20)=1;
  g1(3,25)=(-1);
  g1(3,27)=(-1);
  g1(4,18)=(-1);
  g1(4,20)=1;
  g1(4,7)=(-1);
  g1(5,21)=(-(T103*1/(T108*params(14))));
  g1(5,2)=(-T103);
  g1(5,23)=1;
  g1(5,49)=(-(gamma__*beta_bar__*T103));
  g1(5,33)=(-1);
  g1(6,46)=(-(crk__/(1-params(4)+crk__)));
  g1(6,21)=1;
  g1(6,47)=(-((1-params(4))/(1-params(4)+crk__)));
  g1(6,51)=(-1);
  g1(6,28)=1;
  g1(6,31)=1;
  g1(7,1)=(-(T126/(1+T126)));
  g1(7,22)=1;
  g1(7,48)=(-(1/(1+T126)));
  g1(7,25)=(-T153);
  g1(7,50)=T153;
  g1(7,51)=(-T130);
  g1(7,28)=T130;
  g1(7,31)=T130;
  g1(8,18)=(-crkky__);
  g1(8,22)=(-ccy__);
  g1(8,23)=(-ciy__);
  g1(8,24)=1;
  g1(8,32)=(-1);
  g1(9,20)=(-(params(9)*params(16)));
  g1(9,24)=1;
  g1(9,25)=(-((1-params(9))*params(16)));
  g1(9,30)=(-params(16));
  g1(10,17)=(-(T184*T199));
  g1(10,4)=(-(params(19)*T184));
  g1(10,26)=1;
  g1(10,51)=(-(gamma__*beta_bar__*T184));
  g1(10,35)=(-1);
  g1(11,1)=(-(T235*(-(T126/(1-T126)))));
  g1(11,22)=(-(T235*1/(1-T126)));
  g1(11,25)=(-(T235*params(22)));
  g1(11,4)=(-(params(17)/(1+gamma__*beta_bar__)));
  g1(11,26)=(1+gamma__*beta_bar__*params(17))/(1+gamma__*beta_bar__);
  g1(11,51)=(-T208);
  g1(11,5)=(-T103);
  g1(11,27)=1-(-T235);
  g1(11,52)=(-T208);
  g1(11,36)=(-1);
  g1(12,3)=params(24);
  g1(12,24)=(-((1-params(26))*params(25)+params(24)));
  g1(12,26)=(-(params(23)*(1-params(26))));
  g1(12,6)=(-params(26));
  g1(12,28)=1;
  g1(12,8)=(-(params(16)*params(24)));
  g1(12,30)=(-((1-params(26))*params(25)*(-params(16))+params(24)*(-params(16))));
  g1(12,34)=(-1);
  g1(13,8)=(-params(28));
  g1(13,30)=1;
  g1(13,53)=(-1);
  g1(14,9)=(-params(29));
  g1(14,31)=1;
  g1(14,54)=(-1);
  g1(15,10)=(-params(34));
  g1(15,32)=1;
  g1(15,53)=(-params(3));
  g1(15,55)=(-1);
  g1(16,11)=(-params(32));
  g1(16,33)=1;
  g1(16,56)=(-1);
  g1(17,12)=(-params(33));
  g1(17,34)=1;
  g1(17,57)=(-1);
  g1(18,13)=(-params(30));
  g1(18,35)=1;
  g1(18,58)=(-1);
  g1(18,15)=params(12);
  g1(19,14)=(-params(31));
  g1(19,36)=1;
  g1(19,59)=(-1);
  g1(19,16)=params(11);
  g1(20,23)=(-cikbar__);
  g1(20,7)=(-(1-cikbar__));
  g1(20,29)=1;
  g1(20,33)=(-(params(14)*T108*cikbar__));
  g1(21,3)=1;
  g1(21,24)=(-1);
  g1(21,40)=1;
  g1(22,1)=1;
  g1(22,22)=(-1);
  g1(22,41)=1;
  g1(23,2)=1;
  g1(23,23)=(-1);
  g1(23,42)=1;
  g1(24,5)=1;
  g1(24,27)=(-1);
  g1(24,43)=1;
  g1(25,26)=(-1);
  g1(25,39)=1;
  g1(26,28)=(-1);
  g1(26,38)=1;
  g1(27,25)=(-1);
  g1(27,37)=1;
  g1(28,58)=(-1);
  g1(28,44)=1;
  g1(29,59)=(-1);
  g1(29,45)=1;

if nargout >= 3,
  %
  % Hessian matrix
  %

  g2 = sparse([],[],[],29,3481);
if nargout >= 4,
  %
  % Third order derivatives
  %

  g3 = sparse([],[],[],29,205379);
end
end
end
end
