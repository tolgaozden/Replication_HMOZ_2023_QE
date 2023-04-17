function[laplace]=laplace_approximator(likl,mode,sigma)

n=length(mode);

laplace= likl-(n*log(2*pi)+log(det(sigma)))/2;

end



