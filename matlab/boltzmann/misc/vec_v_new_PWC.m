function v = vec_v_new_PWC(parameter_vector, u, dim, root_vec_dr, wt_vec_dr)

dum_beta = parameter_vector(12);
dum_H = parameter_vector(10);

%---------------------------------------------------------------%

A = matrix_A_PWC(parameter_vector, dim, root_vec_dr, wt_vec_dr, dum_H);

K = diag(root_vec_dr(1:dim));

S = diag(sinh(root_vec_dr(1:dim)*dum_H));

M = diag(wt_vec_dr);

o = ones(dim,1);

v = -dum_beta.*inv(A)*((K^2)+(u^2).*eye(dim))*K*M*S*o;


