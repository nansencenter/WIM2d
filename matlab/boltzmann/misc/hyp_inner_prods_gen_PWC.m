function [cosh_cosh, cosh_sinh, sinh_sinh] = hyp_inner_prods_gen_PWC(root_i, root_j, water_depth,...
								 yn1, yn2, yn3)

% - 05.09.05 - %

%-------------------------------------------------------------------------------

% - i \neq j - %

if root_i ~= root_j
 
    sinh_i = sinh(root_i*water_depth);
    sinh_j = sinh(root_j*water_depth);
    
    cosh_i = cosh(root_i*water_depth);
    cosh_j = cosh(root_j*water_depth);
    
    quot = root_i^2 - root_j^2;
    quot = 1 / quot;
    
    %-------------------------------------------%
    
    if yn1 == 0
        
        cosh_cosh = 0;
        
    else

        cosh_cosh = quot*(root_i*sinh_i*cosh_j - root_j*cosh_i*sinh_j);
        
    end
    
    if yn3 == 0
        
        sinh_sinh = 0;
        
    else

        sinh_sinh = quot*(root_i*sinh_j*cosh_i - root_j*cosh_j*sinh_i);
        
    end
    
    if yn2 == 0
        
        cosh_sinh = 0;
        
    else
    
        cosh_sinh = quot*(root_j + root_i*sinh_i*sinh_j - root_j*cosh_i*cosh_j);
        
    end
    
end

% - i = j - %

if root_i == root_j
    
    if root_i == 0
        
        if yn1 == 0
        
        cosh_cosh = 0;
        
        else

        cosh_cosh = dum_H;
        
        end
        
        if yn3 == 0
        
        sinh_sinh = 0;
        
        else

        sinh_sinh = 0;
        
        end
        
        if yn2 == 0
        
        cosh_sinh = 0;
        
        else
    
        cosh_sinh = 0;
        
        end
    
    else
    
    arg = 2*root_i*water_depth;
    
    quot = (4*root_i)^(-1);
    
    sinh2rH = sinh(arg);
    
    cosh2rH = cosh(arg);
    
    %---------------------------------------------%
    if yn1 == 0
        
        cosh_cosh = 0;
        
    else

        cosh_cosh = sinh2rH + arg;
    
        cosh_cosh = cosh_cosh * quot;
        
    end
    
    if yn3 == 0
        
        sinh_sinh = 0;
        
        else

        sinh_sinh = sinh2rH - arg;
    
        sinh_sinh = sinh_sinh * quot;
        
    end
    
    if yn2 == 0
        
        cosh_sinh = 0;
        
        else
    
        cosh_sinh = cosh2rH - 1;
    
        cosh_sinh = cosh_sinh * quot;
        
    end
    
end
    
end
    
    
    
    