function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
% function T = dynamic_resid_tt(T, y, x, params, steady_state, it_)
%
% File created by Dynare Preprocessor from .mod file
%
% Inputs:
%   T             [#temp variables by 1]     double  vector of temporary terms to be filled by function
%   y             [#dynamic variables by 1]  double  vector of endogenous variables in the order stored
%                                                    in M_.lead_lag_incidence; see the Manual
%   x             [nperiods by M_.exo_nbr]   double  matrix of exogenous variables (in declaration order)
%                                                    for all simulation periods
%   steady_state  [M_.endo_nbr by 1]         double  vector of steady state values
%   params        [M_.param_nbr by 1]        double  vector of parameter values in declaration order
%   it_           scalar                     double  time period for exogenous variables for which
%                                                    to evaluate the model
%
% Output:
%   T           [#temp variables by 1]       double  vector of temporary terms
%

assert(length(T) >= 26);

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
T(12) = T(9)*T(5)*(1-params(9))*1/params(5)/params(9)/T(10);
T(13) = T(8)*T(9);
T(14) = T(5)*T(9);
T(15) = 1-(1-params(4))/T(3);
T(16) = 100*(T(4)-1);
T(17) = 1/(params(13)/(1-params(13)));
T(18) = 1/(1+T(3)*T(11));
T(19) = T(3)^2;
T(20) = params(15)/T(3);
T(21) = (1-T(20))/(params(21)*(1+T(20)));
T(22) = (params(21)-1)*T(12)/(params(21)*(1+T(20)));
T(23) = 1/(1+T(3)*T(11)*params(19));
T(24) = (1-params(20))*(1-T(3)*T(11)*params(20))/params(20)/(1+(params(16)-1)*params(2));
T(25) = T(3)*T(11)/(1+T(3)*T(11));
T(26) = (1-params(18))*(1-T(3)*T(11)*params(18))/((1+T(3)*T(11))*params(18))*1/(1+(params(5)-1)*params(1));

end
