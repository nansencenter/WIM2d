function alpha = alpha_fn_PWC(parameter_vector)

%---------------------------------------

liquid_density = parameter_vector(1);
ice_density = parameter_vector(2);

%-------------------------------------

% required function

alpha = (ice_density * parameter_vector(9)) / liquid_density ;