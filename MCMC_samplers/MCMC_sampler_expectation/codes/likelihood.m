function [likelihood]= likelihood(param,model)

[pr]=SW_prior(param,model);
try
[likl]=kalmanSW(param,model);
likelihood=likl+pr;
catch
    likelihood = Inf;
end

param=reshape(param,[length(param),1]);
if sum(double(param>model.UB))>0 || sum(double(param<model.LB))>0
likelihood=Inf;
end


if isnan(likelihood)==1
    likelihood=Inf;
end

end