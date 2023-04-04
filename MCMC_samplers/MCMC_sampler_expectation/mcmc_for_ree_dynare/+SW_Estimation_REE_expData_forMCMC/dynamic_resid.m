function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
% function residual = dynamic_resid(T, y, x, params, steady_state, it_, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double   vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double   vector of endogenous variables in the order stored
%                                                     in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double   matrix of exogenous variables (in declaration order)
%                                                     for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double   vector of steady state values
%   params        [M_.param_nbr by 1]        double   vector of parameter values in declaration order
%   it_           scalar                     double   time period for exogenous variables for which
%                                                     to evaluate the model
%   T_flag        boolean                    boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = SW_Estimation_REE_expData_forMCMC.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(34, 1);
lhs = y(16);
rhs = params(9)*y(18)+(1-params(9))*y(26)-y(29);
residual(1) = lhs - rhs;
lhs = y(17);
rhs = y(18)*1/(params(13)/(1-params(13)));
residual(2) = lhs - rhs;
lhs = y(18);
rhs = y(26)+y(24)-y(19);
residual(3) = lhs - rhs;
lhs = y(19);
rhs = y(17)+y(7);
residual(4) = lhs - rhs;
lhs = y(22);
rhs = 1/(1+T(3)*T(11))*(y(2)+T(3)*T(11)*y(53)+1/(T(3)^2*params(14))*y(20))+y(32);
residual(5) = lhs - rhs;
lhs = y(20);
rhs = y(36)-y(27)-y(30)+T(5)/(1-params(4)+T(5))*y(50)+(1-params(4))/(1-params(4)+T(5))*y(51);
residual(6) = lhs - rhs;
lhs = y(21);
rhs = params(15)/T(3)/(1+params(15)/T(3))*y(1)+1/(1+params(15)/T(3))*y(52)+(params(21)-1)*T(12)/(params(21)*(1+params(15)/T(3)))*(y(24)-y(54))-(1-params(15)/T(3))/(params(21)*(1+params(15)/T(3)))*(y(30)+y(27)-y(36));
residual(7) = lhs - rhs;
lhs = y(23);
rhs = T(10)*y(21)+y(22)*T(13)+y(31)+y(17)*T(14);
residual(8) = lhs - rhs;
lhs = y(23);
rhs = params(16)*(y(29)+params(9)*y(19)+(1-params(9))*y(24));
residual(9) = lhs - rhs;
lhs = y(25);
rhs = 1/(1+T(3)*T(11)*params(19))*(T(3)*T(11)*y(36)+params(19)*y(4)+y(16)*(1-params(20))*(1-T(3)*T(11)*params(20))/params(20)/(1+(params(16)-1)*params(2)))+y(34);
residual(10) = lhs - rhs;
lhs = y(26);
rhs = 1/(1+T(3)*T(11))*y(5)+T(3)*T(11)/(1+T(3)*T(11))*y(56)+y(4)*params(17)/(1+T(3)*T(11))-y(25)*(1+T(3)*T(11)*params(17))/(1+T(3)*T(11))+y(36)*T(3)*T(11)/(1+T(3)*T(11))+(1-params(18))*(1-T(3)*T(11)*params(18))/((1+T(3)*T(11))*params(18))*1/(1+(params(5)-1)*params(1))*(y(24)*params(22)+y(21)*1/(1-params(15)/T(3))-y(1)*params(15)/T(3)/(1-params(15)/T(3))-y(26))+y(35);
residual(11) = lhs - rhs;
lhs = y(27);
rhs = y(25)*params(23)*(1-params(26))+(1-params(26))*params(25)*(y(23)-params(16)*y(29))+params(24)*(y(23)-params(16)*y(29)-y(3)+params(16)*y(8))+params(26)*y(6)+y(33);
residual(12) = lhs - rhs;
lhs = y(29);
rhs = y(8)*params(28)+x(it_, 1);
residual(13) = lhs - rhs;
lhs = y(30);
rhs = params(29)*y(9)+x(it_, 2);
residual(14) = lhs - rhs;
lhs = y(31);
rhs = params(34)*y(10)+x(it_, 3);
residual(15) = lhs - rhs;
lhs = y(32);
rhs = params(32)*y(11)+x(it_, 4);
residual(16) = lhs - rhs;
lhs = y(33);
rhs = params(33)*y(12)+x(it_, 5);
residual(17) = lhs - rhs;
lhs = y(34);
rhs = params(30)*y(13)+x(it_, 6);
residual(18) = lhs - rhs;
lhs = y(35);
rhs = params(31)*y(14)+x(it_, 7);
residual(19) = lhs - rhs;
lhs = y(28);
rhs = y(7)*(1-T(15))+y(22)*T(15)+y(32)*params(14)*T(3)^2*T(15);
residual(20) = lhs - rhs;
lhs = y(36);
rhs = x(it_, 8)+y(59);
residual(21) = lhs - rhs;
lhs = y(45);
rhs = y(59);
residual(22) = lhs - rhs;
lhs = y(46);
rhs = y(15);
residual(23) = lhs - rhs;
lhs = y(40);
rhs = params(10)+y(23)-y(3);
residual(24) = lhs - rhs;
lhs = y(41);
rhs = params(10)+y(21)-y(1);
residual(25) = lhs - rhs;
lhs = y(42);
rhs = params(10)+y(22)-y(2);
residual(26) = lhs - rhs;
lhs = y(43);
rhs = params(10)+y(26)-y(5);
residual(27) = lhs - rhs;
lhs = y(39);
rhs = params(7)+y(25);
residual(28) = lhs - rhs;
lhs = y(38);
rhs = y(27)+T(16);
residual(29) = lhs - rhs;
lhs = y(37);
rhs = y(24)+params(6);
residual(30) = lhs - rhs;
lhs = y(44);
rhs = params(7)+y(36);
residual(31) = lhs - rhs;
lhs = y(47);
rhs = y(55);
residual(32) = lhs - rhs;
lhs = y(48);
rhs = y(57);
residual(33) = lhs - rhs;
lhs = y(49);
rhs = y(58);
residual(34) = lhs - rhs;

end
