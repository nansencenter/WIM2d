function w = pure_imag_root_PWC(parameter_vector, no, Tol_vec)

% - 28.09.05 - %

tol_x = Tol_vec(2);

w_hold = -2*tol_x;
mk=0;

kappa = parameter_vector(6);
dum_H = parameter_vector(10);

al = kappa*parameter_vector(11);
be = parameter_vector(12);

res = 100;

up = pi * no / dum_H;
mid = ((2*no-1) * pi ) / (2*dum_H);
low = ((no-1) * pi ) / (dum_H);

% ------------------------- %

% ------------------------- %      #1         % ------------------------- %

if and(dr_pi_i_PWC(mid, parameter_vector)>0, al <= 1) % no chance of multi roots
  
low0 = acot(-dr_pi_i_PWC(mid, parameter_vector)/kappa);
low0 = up + (low0/dum_H);

% ------------- %

    count = 1;
    
    w_hold  = low0;
    
    while count ~= 0
  
    interval = linspace(w_hold, up, res);
    
    ed = res;
    
    while dr_pi_ii_PWC(interval(ed), parameter_vector) >  0
        
        ed = ed - 1;
        
    end
    
    if dr_pi_ii_PWC(interval(ed), parameter_vector) + ...
            dr_pi_i_PWC(interval(ed), parameter_vector) < 0
        
        count = 0;
        
    else
        
        count = count + 1;
        
    end
    
    if and(count ~= 0, abs(interval(1)-interval(2)) < tol_x)
        
        mk = 1;
        
        break
        
    end
    
    w_hold = interval(ed);
    
    end
    
    vec = zeros(res,1);

    interval = linspace(low0, w_hold, res);
    
% ------------- %
  
for loop = 2:res
    
    if (dr_pi_i_PWC(interval(loop-1), parameter_vector)+dr_pi_ii_PWC(interval(loop-1), parameter_vector))*...
            (dr_pi_i_PWC(interval(loop), parameter_vector)+dr_pi_ii_PWC(interval(loop), parameter_vector)) < 0
        
        vec(loop)=1;
        
    end
    
end
  
non_zero = find(vec); 

step = 1;

while step <= length(non_zero)
    
    w(step) = i*fzero('pure_imag_dispersion_relation_PWC', [interval(non_zero(step)-1), interval(non_zero(step))],...
        optimset('tolX', tol_x) , parameter_vector);
    
    step = step + 1;
    
end

if mk == 1
    
    w(step) = i*w_hold;
    
end
    
% ------------------------- %      #2         % ------------------------- %

elseif and(dr_pi_i_PWC(mid, parameter_vector)>0, al>1) % - chance of multi roots
    
    min = (al - 1) / (5*be);
    min = min^(1/4);
 
    low0 = acot(-dr_pi_i_PWC(min, parameter_vector)/kappa);
    low0 = low + (low0/dum_H);

      vec = zeros(res,1);

      l_interval = linspace(low0, mid, res);

      for loop = 2:res
    
            if (dr_pi_i_PWC(l_interval(loop-1), parameter_vector)+dr_pi_ii_PWC(l_interval(loop-1), parameter_vector))*...
                (dr_pi_i_PWC(l_interval(loop), parameter_vector)+dr_pi_ii_PWC(l_interval(loop), parameter_vector)) < 0
        
                vec(loop)=1;
        
            end
    
        end

        l_non_zero = find(vec); 
        
        for loop=1:length(l_non_zero)

            w(loop) = i*fzero('pure_imag_dispersion_relation_PWC', [l_interval(l_non_zero(loop)-1), l_interval(l_non_zero(loop))],...
                    optimset('tolX', tol_x), parameter_vector);
            
        end
            
            % ------------------------ %
      
        low0 = acot(-dr_pi_i_PWC(mid, parameter_vector)/kappa);
        low0 = up + (low0/dum_H);

% ------------- %

    count = 1;
    
    w_hold  = mid;
    
    while count ~= 0
  
    r_interval = linspace(w_hold, up, res);
    
    ed = res;
    
    while dr_pi_ii_PWC(r_interval(ed), parameter_vector) >  0
        
        ed = ed - 1;
        
    end
    
    if dr_pi_ii_PWC(r_interval(ed), parameter_vector) + ...
            dr_pi_i_PWC(r_interval(ed), parameter_vector) < 0
        
        count = 0;
        
    else
        
        count = count + 1;
        
    end
    
    if and(count ~= 0, abs(r_interval(1)-r_interval(2)) < tol_x)
        
        mk = 1;
        
        break
        
    end
    
    w_hold = r_interval(ed);
    
    end
    
    vec = zeros(res,1);

    r_interval = linspace(mid, w_hold, res);
    
        % ----------------------- %

        for loop = 2:res
    
            if (dr_pi_i_PWC(r_interval(loop-1), parameter_vector)+dr_pi_ii_PWC(r_interval(loop-1), parameter_vector))*...
                (dr_pi_i_PWC(r_interval(loop), parameter_vector)+dr_pi_ii_PWC(r_interval(loop), parameter_vector)) < 0
        
                vec(loop)=1;
        
            end
    
        end

    r_non_zero = find(vec); 

step = 1;

while step <= length(r_non_zero)
    
    w(length(l_non_zero) + step) = i*fzero('pure_imag_dispersion_relation_PWC', [r_interval(r_non_zero(step)-1), r_interval(r_non_zero(step))],...
           optimset('tolX', tol_x) , parameter_vector);
    
       step = step + 1;
    
end

if mk == 1
    
    w(length(l_non_zero) + step) = i*w_hold;
    
end
       
% ------------------------- %      #3         % ------------------------- %

    elseif and(dr_pi_i_PWC(mid, parameter_vector)<0, al>1)
        
        min = (al - 1) / (5*be);
        min = min^(1/4);
 
        low0 = acot(-dr_pi_i_PWC(min, parameter_vector)/kappa);
        low0 = low + (low0/dum_H);

      vec = zeros(res,1);

      l_interval = linspace(low0, mid, res);

      for loop = 2:res
    
            if (dr_pi_i_PWC(l_interval(loop-1), parameter_vector)+dr_pi_ii_PWC(l_interval(loop-1), parameter_vector))*...
                (dr_pi_i_PWC(l_interval(loop), parameter_vector)+dr_pi_ii_PWC(l_interval(loop), parameter_vector)) < 0
        
                vec(loop)=1;
        
            end
    
      end

        l_non_zero = find(vec); 
        
        for loop=1:length(l_non_zero)

            w(loop) = i*fzero('pure_imag_dispersion_relation_PWC', [l_interval(l_non_zero(loop)-1), l_interval(l_non_zero(loop))],...
                    optimset('tolX', tol_x), parameter_vector);
            
        end 
        
        % - ------------- - %
        
        max_deriv = 1 - al + 5*be*(up^4);
        
        if max_deriv > 0
        
% ------------- %

    w_hold  = mid;
    
    count = 1;
    
    while count ~= 0
  
    r_interval = linspace(w_hold, up, res);
    
    ed = res;
    
    while dr_pi_ii_PWC(r_interval(ed), parameter_vector) >  0
        
        ed = ed - 1;
        
    end
    
    if dr_pi_ii_PWC(r_interval(ed), parameter_vector) + ...
            dr_pi_i_PWC(r_interval(ed), parameter_vector) < 0
        
        count = 0;
        
    else
        
        count = count + 1;
        
    end
    
    if and(count ~= 0, abs(r_interval(1)-r_interval(2)) < tol_x)
        
        mk = 1;
        
        break
        
    end
    
    w_hold = r_interval(ed);
    
    end
    
    vec = zeros(res,1);

    r_interval = linspace(mid, w_hold, res);
    
% ------------- %

        for loop = 2:res
    
            if (dr_pi_i_PWC(r_interval(loop-1), parameter_vector)+dr_pi_ii_PWC(r_interval(loop-1), parameter_vector))*...
                (dr_pi_i_PWC(r_interval(loop), parameter_vector)+dr_pi_ii_PWC(r_interval(loop), parameter_vector)) < 0
        
                vec(loop)=1;
        
            end
    
        end

    r_non_zero = find(vec); 

step = 1;

while step <= length(r_non_zero)
    
    w(length(l_non_zero) + step) = i*fzero('pure_imag_dispersion_relation_PWC', [r_interval(r_non_zero(step)-1), r_interval(r_non_zero(step))],...
           optimset('tolX', tol_x), parameter_vector);
    
       step = step + 1;
    
end

if mk == 1
    
    w(length(l_non_zero) + step) = i*w_hold;
    
end
    
    end

end
    
    