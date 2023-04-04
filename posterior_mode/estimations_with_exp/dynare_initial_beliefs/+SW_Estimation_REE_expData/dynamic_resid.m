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
    T = SW_Estimation_REE_expData.dynamic_resid_tt(T, y, x, params, steady_state, it_);
end
residual = zeros(29, 1);
lhs = y(15);
rhs = params(9)*y(17)+(1-params(9))*y(25)-y(28);
residual(1) = lhs - rhs;
lhs = y(16);
rhs = y(17)*T(17);
residual(2) = lhs - rhs;
lhs = y(17);
rhs = y(25)+y(23)-y(18);
residual(3) = lhs - rhs;
lhs = y(18);
rhs = y(16)+y(7);
residual(4) = lhs - rhs;
lhs = y(21);
rhs = T(18)*(y(2)+T(3)*T(11)*y(47)+1/(T(19)*params(14))*y(19))+y(31);
residual(5) = lhs - rhs;
lhs = y(19);
rhs = y(35)-y(26)-y(29)+T(5)/(1-params(4)+T(5))*y(44)+(1-params(4))/(1-params(4)+T(5))*y(45);
residual(6) = lhs - rhs;
lhs = y(20);
rhs = T(20)/(1+T(20))*y(1)+1/(1+T(20))*y(46)+T(22)*(y(23)-y(48))-T(21)*(y(29)+y(26)-y(35));
residual(7) = lhs - rhs;
lhs = y(22);
rhs = T(10)*y(20)+y(21)*T(13)+y(30)+y(16)*T(14);
residual(8) = lhs - rhs;
lhs = y(22);
rhs = params(16)*(y(28)+params(9)*y(18)+(1-params(9))*y(23));
residual(9) = lhs - rhs;
lhs = y(24);
rhs = T(23)*(T(3)*T(11)*y(35)+params(19)*y(4)+y(15)*T(24))+y(33);
residual(10) = lhs - rhs;
lhs = y(25);
rhs = T(18)*y(5)+T(25)*y(50)+y(4)*params(17)/(1+T(3)*T(11))-y(24)*(1+T(3)*T(11)*params(17))/(1+T(3)*T(11))+y(35)*T(25)+T(26)*(y(23)*params(22)+y(20)*1/(1-T(20))-y(1)*T(20)/(1-T(20))-y(25))+y(34);
residual(11) = lhs - rhs;
lhs = y(26);
rhs = y(24)*params(23)*(1-params(26))+(1-params(26))*params(25)*(y(22)-params(16)*y(28))+params(24)*(y(22)-params(16)*y(28)-y(3)+params(16)*y(8))+params(26)*y(6)+y(32);
residual(12) = lhs - rhs;
lhs = y(28);
rhs = y(8)*params(28)+x(it_, 1);
residual(13) = lhs - rhs;
lhs = y(29);
rhs = params(29)*y(9)+x(it_, 2);
residual(14) = lhs - rhs;
lhs = y(30);
rhs = params(34)*y(10)+x(it_, 3)+x(it_, 1)*params(3);
residual(15) = lhs - rhs;
lhs = y(31);
rhs = params(32)*y(11)+x(it_, 4);
residual(16) = lhs - rhs;
lhs = y(32);
rhs = params(33)*y(12)+x(it_, 5);
residual(17) = lhs - rhs;
lhs = y(33);
rhs = params(30)*y(13)+x(it_, 6);
residual(18) = lhs - rhs;
lhs = y(34);
rhs = params(31)*y(14)+x(it_, 7);
residual(19) = lhs - rhs;
lhs = y(27);
rhs = y(7)*(1-T(15))+y(21)*T(15)+y(31)*params(14)*T(19)*T(15);
residual(20) = lhs - rhs;
lhs = y(35);
rhs = y(49)+x(it_, 8);
residual(21) = lhs - rhs;
lhs = y(39);
rhs = params(10)+y(22)-y(3);
residual(22) = lhs - rhs;
lhs = y(40);
rhs = params(10)+y(20)-y(1);
residual(23) = lhs - rhs;
lhs = y(41);
rhs = params(10)+y(21)-y(2);
residual(24) = lhs - rhs;
lhs = y(42);
rhs = params(10)+y(25)-y(5);
residual(25) = lhs - rhs;
lhs = y(38);
rhs = params(7)+y(24);
residual(26) = lhs - rhs;
lhs = y(37);
rhs = y(26)+T(16);
residual(27) = lhs - rhs;
lhs = y(36);
rhs = y(23)+params(6);
residual(28) = lhs - rhs;
lhs = y(43);
rhs = params(7)+y(35);
residual(29) = lhs - rhs;

end
