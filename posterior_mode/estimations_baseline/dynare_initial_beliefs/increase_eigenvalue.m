function[newm]=increase_eigenvalue(m,eig_crit)

    [VV,DD]=eig(m);
    DD=diag(DD);
    ii=find(DD<eig_crit);
    DD(ii)=eig_crit;

newm=real(VV*diag(DD)*pinv(VV));


end