function C = mat_C_dr_PWC(parameter_vector, dim, root_vec, mu_0, mu_1,...
    wt_vec_dr, dum_beta)

kappa = parameter_vector(6);            

dum_H = parameter_vector(10);

%--------------------------------------%

v_0 = vec_v_new_PWC(parameter_vector, mu_0, dim, root_vec, wt_vec_dr);

v_1 = vec_v_new_PWC(parameter_vector, mu_1, dim, root_vec, wt_vec_dr);

%--------------------------------------%

C = zeros(dim + 2);

C(1:dim,1:dim) = eye(dim);

%--------------------------------------%

for loop = 1:dim
 
    C(dim+1, loop) = wt_vec_dr(loop,1)*root_vec(loop)*sinh(root_vec(loop)*dum_H)/kappa;
    C(dim+2, loop) = -dum_beta * (root_vec(loop)^2) * C(dim+1, loop);
    
end

C(1:dim, dim+1) = v_0;
C(1:dim, dim+2) = v_1;

C(dim+1 , dim+1:dim+2) = ones(1,2);

C(dim+2, dim+1) = -dum_beta*(mu_0^2);
C(dim+2, dim+2) = -dum_beta*(mu_1^2);
