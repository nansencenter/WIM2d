function [mu_0, mu_1] = mu_new_PWC(parameter_vector, dim, root_vec_dr, wt_vec_dr)

% - 18.05.05  - %


%----------------------------------------------------------------%

[c_4, c_2, c_0] = quartic_eqn_new_PWC(parameter_vector, dim, root_vec_dr,  wt_vec_dr);

dum_alpha = parameter_vector(6)*parameter_vector(11);
dum_beta = parameter_vector(12);

% - 22.12.04 - %

v_sq = ((1-dum_alpha) / (c_4*dum_beta))...
     ...
     + ( c_0 / c_4 )...
     ...
     - ( ( c_2 / (2*c_4)) ^ 2 );
 
% - either v_sq is +ve or -ve and real - %
 
if v_sq > 0
    
    % - want the complex root in the upp left quad - %
    
    dummy = 1;
    
    v = sqrt(v_sq);
    
else
    
    % - want the upper half p.i. root or +ve real root - %
    
    dummy = 0;
    
    v = sqrt(-v_sq)*i;
    
end
    
u = c_2 / (2*c_4);

lam_0 = -u + i*v;

% - either lam_0 is on the real line or it is in the upper half plane - %

mod_lam0 = abs(lam_0);
arg_lam0 = angle(lam_0); % \geq 0 \leq pi

lam_1 = -u - i*v;

% - either lam_1 is on the real line or it is in the lower half plane - %

mod_lam1 = abs(lam_1);
arg_lam1 = angle(lam_1); % \geq pi \leq 2*pi

mu_0 = sqrt(mod_lam0)*exp(i*(arg_lam0/2));

mu_1 = sqrt(mod_lam1)*exp(i*(arg_lam1/2) + i*dummy*pi);