function T = static_resid_tt(T, y, x, params)
% function T = static_resid_tt(T, y, x, params)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T         [#temp variables by 1]  double   vector of temporary terms to be filled by function
%   y         [M_.endo_nbr by 1]      double   vector of endogenous variables in declaration order
%   x         [M_.exo_nbr by 1]       double   vector of exogenous variables in declaration order
%   params    [M_.param_nbr by 1]     double   vector of parameter values in declaration order
%
% Output:
%   T         [#temp variables by 1]  double   vector of temporary terms
%

assert(length(T) >= 15);

T(1) = 1+params(7)/100;
T(2) = 0.9995;
T(3) = 1+params(10)/100;
T(4) = T(1)/(T(2)*T(3)^(-params(21)));
T(5) = T(2)^(-1)*T(3)^params(21)-(1-params(4));
T(6) = (params(9)^params(9)*(1-params(9))^(1-params(9))/(params(16)*T(5)^params(9)))^(1/(1-params(9)));
T(7) = (1-params(9))/params(9)*T(5)/T(6);
T(8) = T(3)*(1-(1-params(4))/T(3));
T(9) = params(16)*T(7)^(params(9)-1);
T(10) = 1-params(27)-T(8)*T(9);
T(11) = T(2)*T(3)^(-params(21));
T(12) = T(8)*T(9);
T(13) = T(5)*T(9);
T(14) = 1-(1-params(4))/T(3);
T(15) = 100*(T(4)-1);

end
