function r = GetRootsMMA_FS_PWC(parameter_vector, Dim, Tol_vec)

if Dim==1
    
    r = CalcRealRoot_PWC([parameter_vector(7), parameter_vector(6)], ...
            'FS_DispRel_PWC', 'UppLimReal_FS_PWC', Tol_vec); 
    else

    r = 1i*pure_imag_root_FS_PWC(parameter_vector, Dim-1, Tol_vec);
    
end
    
return