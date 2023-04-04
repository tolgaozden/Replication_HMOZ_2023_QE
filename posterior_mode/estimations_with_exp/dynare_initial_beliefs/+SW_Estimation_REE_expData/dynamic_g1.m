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
    T = SW_Estimation_REE_expData.dynamic_g1_tt(T, y, x, params, steady_state, it_);
end
g1 = zeros(29, 58);
g1(1,15)=1;
g1(1,17)=(-params(9));
g1(1,25)=(-(1-params(9)));
g1(1,28)=1;
g1(2,16)=1;
g1(2,17)=(-T(17));
g1(3,17)=1;
g1(3,18)=1;
g1(3,23)=(-1);
g1(3,25)=(-1);
g1(4,16)=(-1);
g1(4,18)=1;
g1(4,7)=(-1);
g1(5,19)=(-(T(18)*1/(T(19)*params(14))));
g1(5,2)=(-T(18));
g1(5,21)=1;
g1(5,47)=(-(T(3)*T(11)*T(18)));
g1(5,31)=(-1);
g1(6,44)=(-(T(5)/(1-params(4)+T(5))));
g1(6,19)=1;
g1(6,45)=(-((1-params(4))/(1-params(4)+T(5))));
g1(6,26)=1;
g1(6,29)=1;
g1(6,35)=(-1);
g1(7,1)=(-(T(20)/(1+T(20))));
g1(7,20)=1;
g1(7,46)=(-(1/(1+T(20))));
g1(7,23)=(-T(22));
g1(7,48)=T(22);
g1(7,26)=T(21);
g1(7,29)=T(21);
g1(7,35)=(-T(21));
g1(8,16)=(-T(14));
g1(8,20)=(-T(10));
g1(8,21)=(-T(13));
g1(8,22)=1;
g1(8,30)=(-1);
g1(9,18)=(-(params(9)*params(16)));
g1(9,22)=1;
g1(9,23)=(-((1-params(9))*params(16)));
g1(9,28)=(-params(16));
g1(10,15)=(-(T(23)*T(24)));
g1(10,4)=(-(params(19)*T(23)));
g1(10,24)=1;
g1(10,33)=(-1);
g1(10,35)=(-(T(3)*T(11)*T(23)));
g1(11,1)=(-(T(26)*(-(T(20)/(1-T(20))))));
g1(11,20)=(-(T(26)*1/(1-T(20))));
g1(11,23)=(-(T(26)*params(22)));
g1(11,4)=(-(params(17)/(1+T(3)*T(11))));
g1(11,24)=(1+T(3)*T(11)*params(17))/(1+T(3)*T(11));
g1(11,5)=(-T(18));
g1(11,25)=1+T(26);
g1(11,50)=(-T(25));
g1(11,34)=(-1);
g1(11,35)=(-T(25));
g1(12,3)=params(24);
g1(12,22)=(-((1-params(26))*params(25)+params(24)));
g1(12,24)=(-(params(23)*(1-params(26))));
g1(12,6)=(-params(26));
g1(12,26)=1;
g1(12,8)=(-(params(16)*params(24)));
g1(12,28)=(-((1-params(26))*params(25)*(-params(16))+params(24)*(-params(16))));
g1(12,32)=(-1);
g1(13,8)=(-params(28));
g1(13,28)=1;
g1(13,51)=(-1);
g1(14,9)=(-params(29));
g1(14,29)=1;
g1(14,52)=(-1);
g1(15,10)=(-params(34));
g1(15,30)=1;
g1(15,51)=(-params(3));
g1(15,53)=(-1);
g1(16,11)=(-params(32));
g1(16,31)=1;
g1(16,54)=(-1);
g1(17,12)=(-params(33));
g1(17,32)=1;
g1(17,55)=(-1);
g1(18,13)=(-params(30));
g1(18,33)=1;
g1(18,56)=(-1);
g1(19,14)=(-params(31));
g1(19,34)=1;
g1(19,57)=(-1);
g1(20,21)=(-T(15));
g1(20,7)=(-(1-T(15)));
g1(20,27)=1;
g1(20,31)=(-(params(14)*T(19)*T(15)));
g1(21,49)=(-1);
g1(21,35)=1;
g1(21,58)=(-1);
g1(22,3)=1;
g1(22,22)=(-1);
g1(22,39)=1;
g1(23,1)=1;
g1(23,20)=(-1);
g1(23,40)=1;
g1(24,2)=1;
g1(24,21)=(-1);
g1(24,41)=1;
g1(25,5)=1;
g1(25,25)=(-1);
g1(25,42)=1;
g1(26,24)=(-1);
g1(26,38)=1;
g1(27,26)=(-1);
g1(27,37)=1;
g1(28,23)=(-1);
g1(28,36)=1;
g1(29,35)=(-1);
g1(29,43)=1;

end
