string_='estimation_results_ree_t_1.mat';
index_exclude=35;
indices_all=1:1:35;
prior_=SW_prior_output(string_);

for jj=1:length(index_exclude);
    
    prior_exclude(jj)=prior_(index_exclude(jj));
    indices_all(index_exclude(jj))=[];
end

prior_exclude=sum(prior_exclude);

fh_final=fh+prior_exclude;

param_final=x(indices_all);

%H=nhess_diagonal(objective,x);
H=nearestSPD(H);

H_final=H(indices_all,indices_all);
laplace_final=laplace_approximator(fh_final,param_final,H_final);


FLAG_adjusted_likelihood=1;
save string_;

