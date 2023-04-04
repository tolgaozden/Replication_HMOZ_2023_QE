function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
% function g1 = dynamic_g1(T, y, x, params, steady_state, it_, T_flag)
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
%   g1
%

if T_flag
    T = SW_Estimation_REE_expData_forMCMC.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(34, 67);
g1(1,16)=1;
g1(1,18)=(-params(9));
g1(1,26)=(-(1-params(9)));
g1(1,29)=1;
g1(2,17)=1;
g1(2,18)=(-(1/(params(13)/(1-params(13)))));
g1(3,18)=1;
g1(3,19)=1;
g1(3,24)=(-1);
g1(3,26)=(-1);
g1(4,17)=(-1);
g1(4,19)=1;
g1(4,7)=(-1);
g1(5,20)=(-(1/(1+T(3)*T(11))*1/(T(3)^2*params(14))));
g1(5,2)=(-(1/(1+T(3)*T(11))));
g1(5,22)=1;
g1(5,53)=(-(T(3)*T(11)*1/(1+T(3)*T(11))));
g1(5,32)=(-1);
g1(6,50)=(-(T(5)/(1-params(4)+T(5))));
g1(6,20)=1;
g1(6,51)=(-((1-params(4))/(1-params(4)+T(5))));
g1(6,27)=1;
g1(6,30)=1;
g1(6,36)=(-1);
g1(7,1)=(-(params(15)/T(3)/(1+params(15)/T(3))));
g1(7,21)=1;
g1(7,52)=(-(1/(1+params(15)/T(3))));
g1(7,24)=(-((params(21)-1)*T(12)/(params(21)*(1+params(15)/T(3)))));
g1(7,54)=(params(21)-1)*T(12)/(params(21)*(1+params(15)/T(3)));
g1(7,27)=(1-params(15)/T(3))/(params(21)*(1+params(15)/T(3)));
g1(7,30)=(1-params(15)/T(3))/(params(21)*(1+params(15)/T(3)));
g1(7,36)=(-((1-params(15)/T(3))/(params(21)*(1+params(15)/T(3)))));
g1(8,17)=(-T(14));
g1(8,21)=(-T(10));
g1(8,22)=(-T(13));
g1(8,23)=1;
g1(8,31)=(-1);
g1(9,19)=(-(params(9)*params(16)));
g1(9,23)=1;
g1(9,24)=(-((1-params(9))*params(16)));
g1(9,29)=(-params(16));
g1(10,16)=(-(1/(1+T(3)*T(11)*params(19))*(1-params(20))*(1-T(3)*T(11)*params(20))/params(20)/(1+(params(16)-1)*params(2))));
g1(10,4)=(-(params(19)*1/(1+T(3)*T(11)*params(19))));
g1(10,25)=1;
g1(10,34)=(-1);
g1(10,36)=(-(T(3)*T(11)*1/(1+T(3)*T(11)*params(19))));
g1(11,1)=(-((1-params(18))*(1-T(3)*T(11)*params(18))/((1+T(3)*T(11))*params(18))*1/(1+(params(5)-1)*params(1))*(-(params(15)/T(3)/(1-params(15)/T(3))))));
g1(11,21)=(-((1-params(18))*(1-T(3)*T(11)*params(18))/((1+T(3)*T(11))*params(18))*1/(1+(params(5)-1)*params(1))*1/(1-params(15)/T(3))));
g1(11,24)=(-((1-params(18))*(1-T(3)*T(11)*params(18))/((1+T(3)*T(11))*params(18))*1/(1+(params(5)-1)*params(1))*params(22)));
g1(11,4)=(-(params(17)/(1+T(3)*T(11))));
g1(11,25)=(1+T(3)*T(11)*params(17))/(1+T(3)*T(11));
g1(11,5)=(-(1/(1+T(3)*T(11))));
g1(11,26)=1+(1-params(18))*(1-T(3)*T(11)*params(18))/((1+T(3)*T(11))*params(18))*1/(1+(params(5)-1)*params(1));
g1(11,56)=(-(T(3)*T(11)/(1+T(3)*T(11))));
g1(11,35)=(-1);
g1(11,36)=(-(T(3)*T(11)/(1+T(3)*T(11))));
g1(12,3)=params(24);
g1(12,23)=(-((1-params(26))*params(25)+params(24)));
g1(12,25)=(-(params(23)*(1-params(26))));
g1(12,6)=(-params(26));
g1(12,27)=1;
g1(12,8)=(-(params(16)*params(24)));
g1(12,29)=(-((1-params(26))*params(25)*(-params(16))+params(24)*(-params(16))));
g1(12,33)=(-1);
g1(13,8)=(-params(28));
g1(13,29)=1;
g1(13,60)=(-1);
g1(14,9)=(-params(29));
g1(14,30)=1;
g1(14,61)=(-1);
g1(15,10)=(-params(34));
g1(15,31)=1;
g1(15,62)=(-1);
g1(16,11)=(-params(32));
g1(16,32)=1;
g1(16,63)=(-1);
g1(17,12)=(-params(33));
g1(17,33)=1;
g1(17,64)=(-1);
g1(18,13)=(-params(30));
g1(18,34)=1;
g1(18,65)=(-1);
g1(19,14)=(-params(31));
g1(19,35)=1;
g1(19,66)=(-1);
g1(20,22)=(-T(15));
g1(20,7)=(-(1-T(15)));
g1(20,28)=1;
g1(20,32)=(-(params(14)*T(3)^2*T(15)));
g1(21,36)=1;
g1(21,67)=(-1);
g1(21,59)=(-1);
g1(22,45)=1;
g1(22,59)=(-1);
g1(23,15)=(-1);
g1(23,46)=1;
g1(24,3)=1;
g1(24,23)=(-1);
g1(24,40)=1;
g1(25,1)=1;
g1(25,21)=(-1);
g1(25,41)=1;
g1(26,2)=1;
g1(26,22)=(-1);
g1(26,42)=1;
g1(27,5)=1;
g1(27,26)=(-1);
g1(27,43)=1;
g1(28,25)=(-1);
g1(28,39)=1;
g1(29,27)=(-1);
g1(29,38)=1;
g1(30,24)=(-1);
g1(30,37)=1;
g1(31,36)=(-1);
g1(31,44)=1;
g1(32,55)=(-1);
g1(32,47)=1;
g1(33,57)=(-1);
g1(33,48)=1;
g1(34,58)=(-1);
g1(34,49)=1;

end
