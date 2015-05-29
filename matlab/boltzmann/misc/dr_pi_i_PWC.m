function f = dr_pi_i_PWC(w, parameter_vector)

al = parameter_vector(11);
be = parameter_vector(12);

f = (1 - parameter_vector(6)*al + be*(w^4))*w;