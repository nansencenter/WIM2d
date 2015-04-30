function dr = pimag_disprel_FS_PWC(w, parameter_vector)

dummy_h = parameter_vector(7);

dr = w  + parameter_vector(6)* cot(w* dummy_h) ;