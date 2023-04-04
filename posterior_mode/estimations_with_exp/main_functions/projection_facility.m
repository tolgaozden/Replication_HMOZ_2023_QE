
function [pr_flag,largest_eig,alpha_tt,beta_tt,dd_tt,ee_tt,rr_tt]=...
    projection_facility(gamma1,alpha_old,beta_old,dd_old,ee_old,rr_old,alpha_tt,beta_tt,dd_tt,ee_tt,rr_tt,pr_flag_old,model)

%projection facilitytry largest_eig(i)=max(abs(eigs(gamma1,1)),abs(eigs(beta_tt,1)));
  


try largest_eig=abs(eigs(gamma1,1));
catch
    largest_eig=1.01;
end;


if pr_flag_old==1
    
   pr_flag=1;
    alpha_tt=alpha_old;
    beta_tt=beta_old;
    dd_tt=dd_old;
    ee_tt=ee_old;
    rr_tt=rr_old;

else  

if abs(largest_eig)>1% if eig(gamma1)>1 holds, projection facility always binds (IALM supercedes PLM)
    pr_flag=1;
    alpha_tt=alpha_old;
    beta_tt=beta_old;
    dd_tt=dd_old;
    ee_tt=ee_old;
    rr_tt=rr_old;
else
    pr_flag=pr_flag_old;
end
    
    
end


end




