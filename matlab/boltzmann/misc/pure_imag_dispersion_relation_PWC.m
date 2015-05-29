function dr = pure_imag_dispersion_relation_PWC(w, parameter_vector)

kappa = parameter_vector(6);
dummy_alpha = parameter_vector(11);
dummy_beta = parameter_vector(12);
dummy_cap_H = parameter_vector(10);

dr = (1 - kappa*dummy_alpha + dummy_beta* (w^4) ) * w  + kappa* cot(w* dummy_cap_H) ;