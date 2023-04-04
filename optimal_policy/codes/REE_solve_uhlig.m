function[bb]=REE_solve_uhlig(parameters,model,sysmat)





m_states=13;
Sigma=sysmat.sigma;
AA=sysmat.A(1:m_states,1:m_states);
BB=sysmat.B(1:m_states,1:m_states);
CC=sysmat.C(1:m_states,1:m_states);
DD=sysmat.D(1:m_states,:);


PHI=CC;
LAMBDA=AA;
THETA=-BB;
MANUAL_ROOTS=0;
Xi_manual=1;
DISPLAY_IMMEDIATELY=1;
TOL=1e-6;
warnings='';

Xi_mat=[LAMBDA,THETA;eye(m_states),zeros(m_states,m_states)];
Delta_mat=[PHI,zeros(m_states,m_states);zeros(m_states,m_states),eye(m_states)];




 [Xi_eigvec,Xi_eigval] = eig(Xi_mat,Delta_mat);
     if rank(Xi_eigvec)<m_states,
        message = ['SOLVE.M: Sorry! Xi is not diagonalizable! Cannot solve for PP.         '
                   '         Try to run your program again with DO_QZ = 1.                 '];
        if DISPLAY_IMMEDIATELY, disp(message); end;
        warnings = [warnings;message];
     else
       [Xi_sortabs,Xi_sortindex] = sort(abs(diag(Xi_eigval)));
       Xi_sortvec = Xi_eigvec(1:2*m_states,Xi_sortindex);
       Xi_sortval = diag(Xi_eigval(Xi_sortindex,Xi_sortindex));
       Xi_select = 1 : m_states;
       if imag(Xi_sortval(m_states))~=0,
         if (abs( Xi_sortval(m_states) - conj(Xi_sortval(m_states+1)) ) < TOL),
         % NOTE: THIS LAST LINE MIGHT CREATE PROBLEMS, IF THIS EIGENVALUE OCCURS MORE THAN ONCE!!
         % IF YOU HAVE THAT PROBLEM, PLEASE TRY MANUAL ROOT SELECTION.  
           drop_index = 1;
           while (abs(imag(Xi_sortval(drop_index)))>TOL) & (drop_index < m_states),
             drop_index = drop_index + 1;
           end;
           if drop_index >= m_states,
             message = ['SOLVE.M: You are in trouble. You have complex eigenvalues, and I cannot'
                        '   find a real eigenvalue to drop to only have conjugate-complex pairs.'
                        '   Put differently: your PP matrix will contain complex numbers. Sorry!'
                        '   Try increasing the dimension of your state space. You may then get  '
                        '   sunspots, too.                                                      '];
             if DISPLAY_IMMEDIATELY, disp(message); end;
             warnings = [warnings;message];
           else
             message = ['SOLVE.M: I will drop the lowest real eigenvalue to get real PP.        '
                        '         I hope that is ok. You may have sunspots.                     ']; 
             if DISPLAY_IMMEDIATELY, disp(message); end;
             warnings = [warnings;message];
             Xi_select = [ 1: (drop_index-1), (drop_index+1):(m_states+1)];
           end; % if drop_index >= m_states,
         end; % if (abs( Xi_sortval(m_states) - ...
       end; % if imag(Xi_sortval(m_states))~=0,
       if MANUAL_ROOTS,
         message = ['SOLVE.M: You have chosen to select roots manually.  I am crossing my   '
                    '         fingers that you are doing it correctly.  In particular,      '
                    '         you should have defined Xi_manual.  Type help solve           '
                    '         and inspect SOLVE.M to get further information on how to do it'];
         if DISPLAY_IMMEDIATELY, disp(message); end;
         warnings = [warnings;message];
         if exist('Xi_manual'),
            Xi_select = Xi_manual;
         else
            message = ['SOLVE.M: You have not defined Xi_manual.  Either define it or turn off '
                       '         the manual roots selection procedure with                     '
                       '         MANUAL_ROOTS = 0                                              '
                       '         Right now, I better let your calculations crash - sorry!      '
                       '         If you get results, they are based on previous calculations.  '];
            disp(message);
            warnings = [warnings;message];
         end; % if exist('Xi_manual'),
       else
         if max(Xi_select) < 2*m_states,
           if Xi_sortabs(max(Xi_select)+1) < 1 - TOL,
             message = ['SOLVE.M: You may be in trouble. There are stable roots NOT used for PP.'
                        '         I have used the smallest roots: I hope that is ok.            '  
                        '         If not, try manually selecting your favourite roots.          '
                        '         For manual root selection, take a look at the file solve.m    '
                        '         Watch out for sunspot solutions.                              '
                        '         Better yet: move the time index of some endogenous variables  '
                        '         back by one and turn them into (predetermined) state variables'];
             if DISPLAY_IMMEDIATELY, disp(message); end;
             warnings = [warnings;message];
           end; % if Xi_sortabs(max(Xi_select)+1) < 1 - TOL,
         end; % if max(Xi_select) < 2*m_states,
       end; % if MANUAL_ROOTS,
       if max(abs(Xi_sortval(Xi_select)))  > 1 + TOL,
         message = ['SOLVE.M: You may be in trouble.  There are unstable roots used for PP. '
                    '         Keep your fingers crossed or change your model.               '];
         if DISPLAY_IMMEDIATELY, disp(message); end;
         warnings = [warnings;message];
       end; % if max(abs(Xi_sortval(Xi_select))) ... 
       if abs( max(abs(Xi_sortval(Xi_select))) - 1  ) < TOL,
         message = ['SOLVE.M: Your matrix PP contains a unit root. You probably do not have '
                    '         a unique steady state, do you?  Should not be a problem, but  '
                    '         you do not have convergence back to steady state after a shock'
                    '         and you should better not trust long simulations.             '];
         if DISPLAY_IMMEDIATELY, disp(message); end;
         warnings = [warnings;message];
       end; % if abs( max(abs(Xi_sortval(Xi_select))) - 1 ... 
       Lambda_mat = diag(Xi_sortval(Xi_select));
    
       
       Omega_mat  = [Xi_sortvec((m_states+1):(2*m_states),Xi_select)];
       
    

       if rank(Omega_mat)<m_states,
         message = 'SOLVE.M: Sorry! Omega is not invertible. Cannot solve for PP.          ';
         if DISPLAY_IMMEDIATELY, disp(message); end;
         warnings = [warnings;message];
       else
         PP = Omega_mat*Lambda_mat/Omega_mat;
         PP_imag = imag(PP);
         PP = real(PP);
         if sum(sum(abs(PP_imag))) / sum(sum(abs(PP))) > .000001,
           message = ['SOLVE.M: PP is complex.  I proceed with the real part only.            '  
                      '         Hope that is ok, but you are probably really in trouble!!     '
                      '         You should better check everything carefully and be           '
                      '         distrustful of all results which follow now.                  '];
           if DISPLAY_IMMEDIATELY, disp(message); end;
          % warnings = [warnings;message];
         end; % if sum(sum(abs(PP_imag)))
      end; % if rank(Omega_mat)<m_states,
      % End of calculating the PP matrix.  Now comes the rest.
 
    end; % if rank(Xi_eigvec)<m_states,



    
       
%        try

      bb=PP;
      bb=bb(model.FL_indices,model.BL_indices);
% %        catch
%       load beta_prev.mat;
%       bb=beta_prev;
%        end



end

