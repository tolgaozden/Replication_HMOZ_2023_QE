function [hssn]=nhess(fcn,para,varargin);

if nargout>2
	error('maximum number of output arguments is 2.');
end

% printing format
fmt='%5.0f ';		

npara = length(para);

% Compute Hessian, element by element, fine tune with dxscale
ndx = 6;
h0  = exp(-(6:2:(6+2*(ndx-1))));	% step size

hssn     = zeros(npara,npara);
hessdiag = zeros(ndx,1);

dxscale  = ones(npara,1);			% specify different scales across parameters
dxscale  = sparse(1:npara,1:npara,dxscale,npara,npara);

fx    = feval(fcn,para,varargin{:});	% evaluate function 

disp(' ');
disp(' ');
disp('Computing Hessian..');
disp('--------------------------------------------------');
disp(sprintf(['  diagonal elements     :   ' fmt ],npara));


% Compute diagonal elements first
for seli=1:npara

	h = dxscale(:,seli)*h0;
	
	for i=1:ndx

		% forward point
		paradx = para + h(:,i);

		% backward point
		parady = para - h(:,i);
		
		% evaluate function at forward and backward points
		fdx   = feval(fcn,paradx,varargin{:});
		fdy   = feval(fcn,parady,varargin{:});

		% Hessian
		hessdiag(i) = -(2*fx-fdx-fdy)/(h(seli,i))^2; 

    end

    hssn(seli,seli) = 0.5*(hessdiag(3)+hessdiag(4));
	disp(sprintf(['                            ' fmt ],seli));
end


% Now compute off-diagonal elements
% Make sure that correlations are between -1 and 1
% errorij contains the index of elements that are invalid


hssn = real(hssn);

% varargout = {errorij};
