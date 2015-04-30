function [wm, wp] = ...
    ip_wv_gen_PWC(parameter_vector, roots_0_i, roots_i, roots_j, weights_0_i, weights_i)


[ cosh_cosh, waste_i, waste_ii] = hyp_inner_prods_gen_PWC(parameter_vector, roots_i(1,1), roots_j(1,1), parameter_vector(10), 1, 0, 0);

wp = cosh_cosh*weights_i(1,1);


[ cosh_cosh, waste_i, waste_ii] = hyp_inner_prods_gen_PWC(parameter_vector, roots_0_i(1,1), roots_j(1,1), parameter_vector(10), 1, 0, 0);

wm = cosh_cosh*weights_0_i(1,1);



 
