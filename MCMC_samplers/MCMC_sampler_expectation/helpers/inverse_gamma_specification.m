function [s,nu] = inverse_gamma_specification(mu, sigma2, lb, type, use_fzero_flag, name) % --*-- Unitary tests --*--

% Computes the inverse Gamma hyperparameters from the prior mean and standard deviation.
%
% INPUTS
% - mu               [double]   scalar, prior mean.
% - sigma2           [double]   positive scalar, prior variance.
%lb???

% - type             [integer]  scalar equal to 1 or 2, type of the inverse gamma distribution
% - use_fzero_flag   [logical]  scalar, Use (matlab/octave's implementation of) fzero to solve for nu if true, use
%                               dynare's implementation of the secant method otherwise.
% - name             [string]   name of the parameter or random variable.
%
% OUTPUS
% - s                [double]    scalar, first hyperparameter.
% - nu               [double]    scalar, second hyperparameter.
%
% REMARK
% The call to the matlab's implementation of the secant method is here for testing purpose and should not be used. This routine fails
% more often in finding an interval for nu containing a signe change because it expands the interval on both sides and eventually
% violates  the condition nu>2.

% Copyright (C) 2003-2017 Dynare Team
%
% This file is part of Dynare.
%
% Dynare is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% Dynare is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Dynare.  If not, see <http://www.gnu.org/licenses/>.

if nargin<4
    error('At least four input arguments are required!')
end

if ~isnumeric(mu) || ~isscalar(mu) || ~isreal(mu)
    error('First input argument must be a real scalar!')
end

if ~isnumeric(sigma2) || ~isscalar(sigma2) || ~isreal(sigma2) || sigma2<=0
    error('Second input argument must be a real positive scalar!')
end

if ~isnumeric(lb) || ~isscalar(lb) || ~isreal(lb)
    error('Third input argument must be a real scalar!')
end

if ~isnumeric(type) || ~isscalar(type) || ~ismember(type, [1, 2])
    error('Fourth input argument must be equal to 1 or 2!')
end

if nargin==4 || isempty(use_fzero_flag)
    use_fzero_flag = false;
else
    if ~isscalar(use_fzero_flag) || ~islogical(use_fzero_flag)
        error('Fourth input argument must be a scalar logical!')
    end
end

if nargin>5 && (~ischar(name) || size(name, 1)>1)
    error('Sixth input argument must be a string!')
else
    name = '';
end

if ~isempty(name)
    name = sprintf(' for %s', name);
end

if mu<=lb
    error('The prior mean%s (%f) must be above the lower bound (%f)of the Inverse Gamma (type %d) prior distribution!', mu, lb, name, type);
end

check_solution_flag = true;
s = [];
nu = [];

sigma = sqrt(sigma2);
mu2 = mu*mu;

if type == 2       % Inverse Gamma 2
    nu   = 2*(2+mu2/sigma2);
    s    = 2*mu*(1+mu2/sigma2);
elseif type == 1   % Inverse Gamma 1
    if sigma2 < Inf
        nu = sqrt(2*(2+mu2/sigma2));
        if use_fzero_flag
            nu = fzero(@(nu)ig1fun(nu,mu2,sigma2),nu);
        else
            nu2 = 2*nu;
            nu1 = 2;
            err  = ig1fun(nu,mu2,sigma2);
            err2 = ig1fun(nu2,mu2,sigma2);
            if err2 > 0         % Too short interval.
                while nu2 < 1e12 % Shift the interval containing the root.
                    nu1  = nu2;
                    nu2  = nu2*2;
                    err2 = ig1fun(nu2,mu2,sigma2);
                    if err2<0
                        break
                    end
                end
                if err2>0
                    error('inverse_gamma_specification:: Failed in finding an interval containing a sign change! You should check that the prior variance is not too small compared to the prior mean...');
                end
            end
            % Solve for nu using the secant method.
            while abs(nu2/nu1-1) > 1e-14
                if err > 0
                    nu1 = nu;
                    if nu < nu2
                        nu = nu2;
                    else
                        nu = 2*nu;
                        nu2 = nu;
                    end
                else
                    nu2 = nu;
                end
                nu =  (nu1+nu2)/2;
                err = ig1fun(nu,mu2,sigma2);
            end
        end
        s = (sigma2+mu2)*(nu-2);
        if check_solution_flag
            if abs(log(mu)-log(sqrt(s/2))-gammaln((nu-1)/2)+gammaln(nu/2))>1e-7
                error('inverse_gamma_specification:: Failed in solving for the hyperparameters!');
            end
            if abs(sigma-sqrt(s/(nu-2)-mu*mu))>1e-7
                error('inverse_gamma_specification:: Failed in solving for the hyperparameters!');
            end
        end
    else
        nu  = 2;
        s   = 2*mu2/pi;
    end
else
    error('inverse_gamma_specification: unkown type')
end

%@test:1
%$ try
%$    [s, nu] = inverse_gamma_specification(.5, .05, 0, 1);
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    t(2) = abs(0.5-sqrt(.5*s)*gamma(.5*(nu-1))/gamma(.5*nu))<1e-12;
%$    t(3) = abs(0.05-s/(nu-2)+.5^2)<1e-12;
%$ end
%$ T = all(t);
%@eof:1

%@test:2
%$ try
%$    [s, nu] = inverse_gamma_specification(.5, .05, 0, 2);
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    t(2) = abs(0.5-s/(nu-2))<1e-12;
%$    t(3) = abs(0.05-2*.5^2/(nu-4))<1e-12;
%$ end
%$ T = all(t);
%@eof:2

%@test:3
%$ try
%$    [s, nu] = inverse_gamma_specification(.5, Inf, 0, 1);
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    t(2) = abs(0.5-sqrt(.5*s)*gamma(.5*(nu-1))/gamma(.5*nu))<1e-12;
%$    t(3) = isequal(nu, 2); %abs(0.05-2*.5^2/(nu-4))<1e-12;
%$ end
%$ T = all(t);
%@eof:3

%@test:4
%$ try
%$    [s, nu] = inverse_gamma_specification(.5, Inf, 0, 2);
%$    t(1) = true;
%$ catch
%$    t(1) = false;
%$ end
%$
%$ if t(1)
%$    t(2) = abs(0.5-s/(nu-2))<1e-12;
%$    t(3) = isequal(nu, 4);
%$ end
%$ T = all(t);
%@eof:4