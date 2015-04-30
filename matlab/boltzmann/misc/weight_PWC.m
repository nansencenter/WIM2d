function w = weight_PWC(parameter_vector, root)

H = parameter_vector(10);

w =  sech(root*H);