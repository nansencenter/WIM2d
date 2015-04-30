function m = normalising_weight_PWC(parameter_vector, root)

%-------------------------------------%

[ cosh_cosh, waste1, waste2] = hyp_inner_prods_gen_PWC(root, root, parameter_vector(10), 1, 0, 0);

m = sqrt(1 / cosh_cosh);

