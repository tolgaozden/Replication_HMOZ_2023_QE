function [alpha beta r] =cgl_learning_recursive(x,xOld,alphaOld,betaOld,rOld,gain)
%constant gain sample autocorrelation learning with AR(1) rule, recursive
%form.


% alphaOld=0;

r = rOld+gain*((x-alphaOld)^2-rOld);

alpha= alphaOld+gain*(x-alphaOld);

% alpha=0;

beta= betaOld+ gain*(pinv(r))*((x-alphaOld)*(xOld-alphaOld)-betaOld*(x-alphaOld)^2);



end