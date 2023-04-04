function [vecX]= vec(X)

vecX= reshape(X,[size(X,1)*size(X,2),1]);

end