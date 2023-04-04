%this script is used to generate initial beliefs based on BLE model for
%AR(1) learning models 

%run after KF_output.m with BLE 


beta_init = diag(beta_tt(model.FL_indices,model.FL_indices))
rr_init = rr_tt(:,:,model.FL_indices);
save AR1_BLE_initial_beliefs.mat beta_init rr_init;