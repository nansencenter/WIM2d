function A = matrix_A_PWC(parameter_vector, dim, root_vec, weight_vec, water_depth)

A = zeros(dim);

for loop = 1:dim
    
    for loop2 = 1:dim
        
        A(loop, loop2) = a_ij_gen_PWC(parameter_vector,...
            root_vec(loop), root_vec(loop2),...
            weight_vec(loop), weight_vec(loop2), water_depth);
        
    end
    
end
