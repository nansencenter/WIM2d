function [W_m, W_p] ...
 = jump_W_PWC_wtd(parameter_vector, dimN, rts_v, rts_0, rts, wts_v, wts_0, wts)

% - 03.06.09 - rewritten for PWC case with weights for interfacial fns - %

water_depth=parameter_vector(10); clear parameter_vector

W_m = complex(zeros(dimN,dimN));
W_p = complex(zeros(dimN,dimN));

for loop_i = 1:dimN
 for loop_j=1:dimN
  
  W_m(loop_i,loop_j) = fn_cosh_cosh([], rts_0(loop_i), rts_v(loop_j), water_depth);
  W_p(loop_i,loop_j) = fn_cosh_cosh([], rts(loop_i), rts_v(loop_j), water_depth);
  
  W_m(loop_i,loop_j) = wts_v(loop_j)*wts_0(loop_i)*W_m(loop_i,loop_j);
  W_p(loop_i,loop_j) = wts_v(loop_j)*wts(loop_i)*W_p(loop_i,loop_j);
  
 end
end

%% - Subfns
%% - taken from old code:

function cosh_cosh = fn_cosh_cosh(parameter_vector, root_i, root_j, water_depth)

clear parameter_vector

% - i \neq j - %

if root_i ~= root_j
 
 sinh_i = sinh(root_i*water_depth);
 sinh_j = sinh(root_j*water_depth);
 
 cosh_i = cosh(root_i*water_depth);
 cosh_j = cosh(root_j*water_depth);
 
 quot = root_i^2 - root_j^2;
 quot = 1 / quot;
 
 %-------------------------------------------%
 
 cosh_cosh = quot*(root_i*sinh_i*cosh_j - root_j*cosh_i*sinh_j);
 
 % - i = j - %
elseif root_i == root_j
 
 if root_i == 0
  
  cosh_cosh = dum_H;
  
 else
  
  arg = 2*root_i*water_depth;
  
  quot = (4*root_i)^(-1);
  
  sinh2rH = sinh(arg);
  
  cosh2rH = cosh(arg);
  
  cosh_cosh = sinh2rH + arg;
  
  cosh_cosh = cosh_cosh * quot;
 end
end

return