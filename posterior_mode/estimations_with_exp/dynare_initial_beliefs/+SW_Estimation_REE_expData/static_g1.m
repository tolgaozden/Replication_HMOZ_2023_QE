function g1 = static_g1(T, y, x, params, T_flag)
% function g1 = static_g1(T, y, x, params, T_flag)
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
%   g1
%

if T_flag
    T = SW_Estimation_REE_expData.static_g1_tt(T, y, x, params);
end
g1 = zeros(29, 29);
g1(1,1)=1;
g1(1,3)=(-params(9));
g1(1,11)=(-(1-params(9)));
g1(1,14)=1;
g1(2,2)=1;
g1(2,3)=(-T(16));
g1(3,3)=1;
g1(3,4)=1;
g1(3,9)=(-1);
g1(3,11)=(-1);
g1(4,2)=(-1);
g1(4,4)=1;
g1(4,13)=(-1);
g1(5,5)=(-(T(17)*1/(T(18)*params(14))));
g1(5,17)=(-1);
g1(6,3)=(-(T(5)/(1-params(4)+T(5))));
g1(6,5)=1-(1-params(4))/(1-params(4)+T(5));
g1(6,12)=1;
g1(6,15)=1;
g1(6,21)=(-1);
g1(7,6)=1-(T(19)/(1+T(19))+1/(1+T(19)));
g1(7,12)=T(20);
g1(7,15)=T(20);
g1(7,21)=(-T(20));
g1(8,2)=(-T(13));
g1(8,6)=(-T(10));
g1(8,7)=(-T(12));
g1(8,8)=1;
g1(8,16)=(-1);
g1(9,4)=(-(params(9)*params(16)));
g1(9,8)=1;
g1(9,9)=(-((1-params(9))*params(16)));
g1(9,14)=(-params(16));
g1(10,1)=(-(T(21)*T(22)));
g1(10,10)=1-params(19)*T(21);
g1(10,19)=(-1);
g1(10,21)=(-(T(3)*T(11)*T(21)));
g1(11,6)=(-(T(24)*(1/(1-T(19))-T(19)/(1-T(19)))));
g1(11,9)=(-(T(24)*params(22)));
g1(11,10)=(-(params(17)/(1+T(3)*T(11))-(1+T(3)*T(11)*params(17))/(1+T(3)*T(11))));
g1(11,11)=1-(T(17)+T(23)-T(24));
g1(11,20)=(-1);
g1(11,21)=(-T(23));
g1(12,8)=(-((1-params(26))*params(25)));
g1(12,10)=(-(params(23)*(1-params(26))));
g1(12,12)=1-params(26);
g1(12,14)=(-((1-params(26))*params(25)*(-params(16))));
g1(12,18)=(-1);
g1(13,14)=1-params(28);
g1(14,15)=1-params(29);
g1(15,16)=1-params(34);
g1(16,17)=1-params(32);
g1(17,18)=1-params(33);
g1(18,19)=1-params(30);
g1(19,20)=1-params(31);
g1(20,7)=(-T(14));
g1(20,13)=1-(1-T(14));
g1(20,17)=(-(params(14)*T(18)*T(14)));
g1(21,10)=(-1);
g1(21,21)=1;
g1(22,25)=1;
g1(23,26)=1;
g1(24,27)=1;
g1(25,28)=1;
g1(26,10)=(-1);
g1(26,24)=1;
g1(27,12)=(-1);
g1(27,23)=1;
g1(28,9)=(-1);
g1(28,22)=1;
g1(29,21)=(-1);
g1(29,29)=1;
if ~isreal(g1)
    g1 = real(g1)+2*imag(g1);
end
end
