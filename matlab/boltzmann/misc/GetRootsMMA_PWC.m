function r = GetRootsMMA_PWC(parameter_vector, Dim, Tol_vec)

if Dim==1
    
    r = CalcRealRoot_PWC([parameter_vector(11), parameter_vector(12), ...
        parameter_vector(10), parameter_vector(6)], ...
        'ice_DispRel_PWC', 'UppLimReal_PWC', Tol_vec);
    
    else

    r = pure_imag_root_PWC(parameter_vector, Dim-1, Tol_vec);
    
end

return
    


            
    
    
