function [c_4, c_2, c_0] = quartic_eqn_new_PWC(parameter_vector, dim, root_vec_dr, wt_vec_dr)

%----------------------------------------------------------------%

dum_H = parameter_vector(10);

%------------------------------------------------------------%

one_vec = ones(dim,1);

W = diag(wt_vec_dr);

R = diag(root_vec_dr);

A = matrix_A_PWC(parameter_vector, dim, root_vec_dr, wt_vec_dr, dum_H);

C = diag(cosh(root_vec_dr.*dum_H));

S = diag(sinh(root_vec_dr.*dum_H));

c_4 = 1;

c_2 = one_vec'*R*W*S*inv(A)*W*C*one_vec;

c_0 = one_vec'*(R^3)*W*S*inv(A)*W*C*one_vec;






