function [theta r largest_eig projectionFlag] =msv_learning(x,regressor,thetaOld,rOld,gain,...
    forward_indices,backward_indices,numEndo,ridge_correction)
msgstr=[];
msgid=[];
numBackward=length(backward_indices);
lastwarn('Success');
lambda=10e-8;
projectionFlag=0;
largest_eig=0;
eig_crit=10e-5;

 yy=regressor;
 
 r_low = 0;

%r=increase_eigenvalue(r,lambda,0);
%r=triu(r);
%r=(r+r')/2;
% r=rOld;

%  try
     r=rOld+gain*(yy*yy'-rOld);
     
     if abs(min(eigs(r)))<eig_crit
         r_low = 1;
     end
%      if ridge_correction==1
%      r=increase_eigenvalue(r,eig_crit);
%      end

if r_low == 0
theta=thetaOld'+gain*(pinv(r)*yy)*(x-thetaOld*yy)';
elseif r_low ==1 
theta=thetaOld'+gain*(pinv(r+eig_crit*eye(length(r)))*yy)*(x-thetaOld*yy)';
end
%  catch
% %     projectionFlag=1;
%      r=rOld;
%      theta=thetaOld';
%      projectionFlag=1;
%  end
 [msgstr, msgid] = lastwarn; 
 if strcmp(msgid,'MATLAB:nearlySingularMatrix')==1
 r=rOld;
 theta=thetaOld';
 projectionFlag=1;
 end




%   auxiliary_r=100*eye(length(r));  
%   auxiliary_r=triu(auxiliary_r);
%   auxiliary_r=(auxiliary_r+auxiliary_r')/2;
% % auxiliary_r=nearestSPD(r+eye(length(r))*lambda);
% r=rOld;
%  theta=thetaOld'+gain*(pinv(auxiliary_r)*yy)*(x-thetaOld*yy)';
% end
    
% beta_tt=zeros(numEndo,numEndo);
% beta_tt(forward_indices,backward_indices)=theta(2:numBackward+1,:)';


% gamma1_tilde=(AA-CC*beta_tt)\BB;
% 
% try ev=abs(eigs(gamma1_tilde,1)); catch; ev=1.01;end;
% largest_eig=ev;




% [msgstr, msgid] = lastwarn;
% 
% if strcmp(msgid,'MATLAB:nearlySingularMatrix')==1
% r=rOld;
% theta=thetaOld';
% projectionFlag=1;
% end


% if largest_eig1>1
%     r=rOld;
% theta=thetaOld';
% projectionFlag=1;
% 
% elseif largest_eig2>1
%     
%         r=rOld;
% theta=thetaOld';
% projectionFlag=1;
% 
% end

% if eig_active>1;
%        r=rOld;
% theta=thetaOld';
% projectionFlag=1;
% end

  
% if largest_eig>1
%     r=rOld;
%     theta=thetaOld';
%     projectionFlag=1;
% end



end



