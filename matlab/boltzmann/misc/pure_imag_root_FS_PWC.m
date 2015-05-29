function w = pure_imag_root_FS_PWC(parameter_vector, no, Tol_vec)

tol_x = Tol_vec(2);

dum_h = parameter_vector(7);

up = ( pi * no / dum_h) - (tol_x / dum_h);

low = ( pi * (no - 1) / dum_h) + (tol_x / dum_h);

while pimag_disprel_FS_PWC( up, parameter_vector) > 0
    
    tol_x = tol_x*10;
    
    up = ( pi * no / dum_h) - (tol_x / dum_h);
    
end

tol_x = Tol_vec(2);

while pimag_disprel_FS_PWC(low, parameter_vector) < 0
    
    tol_x = tol_x*10;
    
    low = ( pi * (no - 1) / dum_h) + (tol_x / dum_h);
    
end
 
tol_x = Tol_vec(2);

interval = [low, up];

w = fzero('pimag_disprel_FS_PWC', interval, ...
    optimset('tolX', tol_x),parameter_vector);

return