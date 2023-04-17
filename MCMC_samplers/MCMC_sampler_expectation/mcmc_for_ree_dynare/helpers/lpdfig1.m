function [ldens] = lpdfig1(x,s,nu)
% Evaluates the logged INVERSE-GAMMA-1 PDF at x.
%
% X ~ IG1(s,nu) if X = sqrt(Y) where Y ~ IG2(s,nu) and Y = inv(Z) with Z ~ G(nu/2,2/s) (Gamma distribution)
%
% See L. Bauwens, M. Lubrano and J-F. Richard [1999, appendix A] for more details.
%
%
% INPUTS
%    x     [double]  m*n matrix of locations,
%    s     [double]  m*n matrix or scalar, First INVERSE-GAMMA-1 distribution parameters,
%    nu    [double]  m*n matrix or scalar, Second INVERSE-GAMMA-1 distribution parameters.
%
% OUTPUTS
%    ldens [double]  m*n matrix of logged INVERSE-GAMMA-1 densities evaluated at x.



ldens = -Inf( size(x) ) ;
idx = find( x>0 ) ;

if length(s)==1
    ldens(idx) = log(2) - gammaln(.5*nu) - .5*nu*(log(2)-log(s)) - (nu+1)*log(x(idx)) - .5*s./(x(idx).*x(idx)) ;
else
    ldens(idx) = log(2) - gammaln(.5*nu(idx)) - .5*nu(idx).*(log(2)-log(s(idx))) - (nu(idx)+1).*log(x(idx)) - .5*s(idx)./(x(idx).*x(idx)) ;
end


%ldens=exp(ldens);
end
