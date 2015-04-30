function [W_m, W_p] ...
        = jump_W_PWC(parameter_vector, dimN, dimM, v_vec, roots_0, root_vec, weights_0, weights)
    
    % - 29.6.06 - V is NxM - %

for loopN = 1:dimN
    
    for loopM=1:dimM
        
        [W_m(loopN,loopM), W_p(loopN,loopM)] = ...
            ...
            ip_wv_gen_PWC(parameter_vector, roots_0(loopN,:), root_vec(loopN,:), v_vec(loopM,:),...
            ...
            weights_0(loopN,:), weights(loopN,:));
        
    end
    
end