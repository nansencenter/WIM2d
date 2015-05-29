function F = dr_pi_ii_PWC(w, parameter_vector)

H = parameter_vector(10);

kappa = parameter_vector(6);

F = kappa/tan(w*H);