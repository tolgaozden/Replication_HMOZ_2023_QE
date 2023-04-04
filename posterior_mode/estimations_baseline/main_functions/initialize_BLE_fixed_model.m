clear;
load kf_output;

rr_init=(rr_tt(:,:,model.FL_indices))
beta_init=diag(beta_tt(model.FL_indices,model.FL_indices))

clearvars -except rr_init beta_init;


save inputs/AR1_BLE_initial_beliefs.mat;
