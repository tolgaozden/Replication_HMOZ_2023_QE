function residual = static_resid(T, y, x, params, T_flag)
% function residual = static_resid(T, y, x, params, T_flag)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%                                              to evaluate the model
%   T_flag    boolean                 boolean  flag saying whether or not to calculate temporary terms
%
% Output:
%   residual
%

if T_flag
    T = SW_Estimation_REE_expData_forMCMC.static_resid_tt(T, y, x, params);
end
residual = zeros(34, 1);
lhs = y(1);
rhs = params(9)*y(3)+(1-params(9))*y(11)-y(14);
residual(1) = lhs - rhs;
lhs = y(2);
rhs = y(3)*1/(params(13)/(1-params(13)));
residual(2) = lhs - rhs;
lhs = y(3);
rhs = y(11)+y(9)-y(4);
residual(3) = lhs - rhs;
lhs = y(4);
rhs = y(2)+y(13);
residual(4) = lhs - rhs;
lhs = y(7);
rhs = 1/(1+T(3)*T(11))*(y(7)+y(7)*T(3)*T(11)+1/(T(3)^2*params(14))*y(5))+y(17);
residual(5) = lhs - rhs;
lhs = y(5);
rhs = y(21)-y(12)-y(15)+y(3)*T(5)/(1-params(4)+T(5))+y(5)*(1-params(4))/(1-params(4)+T(5));
residual(6) = lhs - rhs;
lhs = y(6);
rhs = y(6)*params(15)/T(3)/(1+params(15)/T(3))+y(6)*1/(1+params(15)/T(3))-(1-params(15)/T(3))/(params(21)*(1+params(15)/T(3)))*(y(15)+y(12)-y(21));
residual(7) = lhs - rhs;
lhs = y(8);
rhs = T(10)*y(6)+y(7)*T(12)+y(16)+y(2)*T(13);
residual(8) = lhs - rhs;
lhs = y(8);
rhs = params(16)*(y(14)+params(9)*y(4)+(1-params(9))*y(9));
residual(9) = lhs - rhs;
lhs = y(10);
rhs = 1/(1+T(3)*T(11)*params(19))*(T(3)*T(11)*y(21)+y(10)*params(19)+y(1)*(1-params(20))*(1-T(3)*T(11)*params(20))/params(20)/(1+(params(16)-1)*params(2)))+y(19);
residual(10) = lhs - rhs;
lhs = y(11);
rhs = y(11)*1/(1+T(3)*T(11))+y(11)*T(3)*T(11)/(1+T(3)*T(11))+y(10)*params(17)/(1+T(3)*T(11))-y(10)*(1+T(3)*T(11)*params(17))/(1+T(3)*T(11))+y(21)*T(3)*T(11)/(1+T(3)*T(11))+(1-params(18))*(1-T(3)*T(11)*params(18))/((1+T(3)*T(11))*params(18))*1/(1+(params(5)-1)*params(1))*(y(9)*params(22)+y(6)*1/(1-params(15)/T(3))-y(6)*params(15)/T(3)/(1-params(15)/T(3))-y(11))+y(20);
residual(11) = lhs - rhs;
lhs = y(12);
rhs = y(10)*params(23)*(1-params(26))+(1-params(26))*params(25)*(y(8)-params(16)*y(14))+params(24)*(params(16)*y(14)+y(8)-params(16)*y(14)-y(8))+y(12)*params(26)+y(18);
residual(12) = lhs - rhs;
lhs = y(14);
rhs = y(14)*params(28)+x(1);
residual(13) = lhs - rhs;
lhs = y(15);
rhs = y(15)*params(29)+x(2);
residual(14) = lhs - rhs;
lhs = y(16);
rhs = y(16)*params(34)+x(3);
residual(15) = lhs - rhs;
lhs = y(17);
rhs = y(17)*params(32)+x(4);
residual(16) = lhs - rhs;
lhs = y(18);
rhs = y(18)*params(33)+x(5);
residual(17) = lhs - rhs;
lhs = y(19);
rhs = y(19)*params(30)+x(6);
residual(18) = lhs - rhs;
lhs = y(20);
rhs = y(20)*params(31)+x(7);
residual(19) = lhs - rhs;
lhs = y(13);
rhs = y(13)*(1-T(14))+y(7)*T(14)+y(17)*params(14)*T(3)^2*T(14);
residual(20) = lhs - rhs;
lhs = y(21);
rhs = y(10)+x(8);
residual(21) = lhs - rhs;
lhs = y(30);
rhs = y(10);
residual(22) = lhs - rhs;
lhs = y(31);
rhs = y(30);
residual(23) = lhs - rhs;
lhs = y(25);
rhs = params(10);
residual(24) = lhs - rhs;
lhs = y(26);
rhs = params(10);
residual(25) = lhs - rhs;
lhs = y(27);
rhs = params(10);
residual(26) = lhs - rhs;
lhs = y(28);
rhs = params(10);
residual(27) = lhs - rhs;
lhs = y(24);
rhs = params(7)+y(10);
residual(28) = lhs - rhs;
lhs = y(23);
rhs = y(12)+T(15);
residual(29) = lhs - rhs;
lhs = y(22);
rhs = y(9)+params(6);
residual(30) = lhs - rhs;
lhs = y(29);
rhs = params(7)+y(21);
residual(31) = lhs - rhs;
lhs = y(32);
rhs = y(10);
residual(32) = lhs - rhs;
lhs = y(33);
rhs = y(10);
residual(33) = lhs - rhs;
lhs = y(34);
rhs = y(10);
residual(34) = lhs - rhs;
if ~isreal(residual)
  residual = real(residual)+imag(residual).^2;
end
end
