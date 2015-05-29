function m = normalising_weight_0_PWC(bed, root)

[cosh_cosh, waste_1, waste_2] = hyp_inner_prods_gen_PWC(root, root, bed, 1, 0, 0);

m = sqrt(1/cosh_cosh);