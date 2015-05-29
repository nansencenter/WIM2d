function a_ij = a_ij_gen_PWC(parameter_vector, root_i, root_j, m_i, m_j, water_depth)

%-------------------------------------------------------------------------------
   
[ cosh_cosh, no_use, no_use1] = hyp_inner_prods_gen_PWC(root_i, root_j, water_depth, 1, 0, 0);

a_ij = m_i*m_j*cosh_cosh;


