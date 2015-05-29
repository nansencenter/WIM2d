% function w = fn_ImagRoot_ice(no, al, be, dum_H, kappa, Tol_vec)
%
% L Bennetts Sept 2013 / Adelaide
%
% Last updated:  Mar 2015

function w = fn_ImagRoot_ice(no, al, be, dum_H, kappa, Tol_vec)

res=100;

be_hat = be/(dum_H^4);
al_hat = 1-kappa*al;
kap_hat = kappa*dum_H;

% interval
i0=(no-1)*pi; %i1=no*pi;

% calculate minimum of polynomial part of disp rel:
mn = -al_hat / 5 / be_hat; mn = mn^(1/4);

if and(mn<=i0,fn_poly_dr((no-1)*pi,al_hat,be_hat)>=0)  
 %% No chance of multiple roots:
 
 low = (no-0.5)*pi;
 upp = (no+0.0)*pi;
 
 w = fzero('fn_ImDispRel_ice', [low,upp], optimset('tolX', Tol_vec(2)),...
   [al_hat,be_hat,kap_hat])/dum_H;
 
else
 %% Chance of multiple roots:
 
 mid  = (no-0.5)*pi;
 
 %low0 = acot(-fn_poly_dr(mn,al_hat,be_hat)/kap_hat);
 low0 = (no-1)*pi; % + low0;
 
 % calculate any roots in interval to mid-pt:
 vec = zeros(res,1);
 
 l_interval = linspace(low0, mid, res);
 
 for loop = 2:res
  
%   if (fn_poly_dr(l_interval(loop-1),al_hat,be_hat)+fn_cot_dr(l_interval(loop-1), kap_hat))*...
%     (fn_poly_dr(l_interval(loop),al_hat,be_hat)+fn_cot_dr(l_interval(loop), kap_hat)) < 0

  if fn_ImDispRel_ice(l_interval(loop-1),[al_hat,be_hat,kap_hat])*...
   fn_ImDispRel_ice(l_interval(loop),[al_hat,be_hat,kap_hat]) < 0
   
   vec(loop)=1;
   
  end
  
 end
 
 l_non_zero = find(vec); clear vec
 
 w=[];
 
 for loop=1:length(l_non_zero)
  
  w(loop) = fzero('fn_ImDispRel_ice', ...
   [l_interval(l_non_zero(loop)-1), l_interval(l_non_zero(loop))],...
   optimset('tolX', Tol_vec(2)), [al_hat,be_hat,kap_hat])/dum_H;
  
 end
 
 % 
 
 upp  = no*pi;
 
 %max_deriv = al_hat + 5*be_hat*(upp^4);
 
 % If Poly is negative at RH end => no roots in 2nd half of interval
 if fn_poly_dr(upp,al_hat,be_hat) > 0
  
  % V1:
  
  N = length(w);
  
  r_interval = linspace(mid, upp, res);
  
  for loop = 2:res
   if fn_ImDispRel_ice(r_interval(loop-1),[al_hat,be_hat,kap_hat])*...
    fn_ImDispRel_ice(r_interval(loop),[al_hat,be_hat,kap_hat]) < 0
    vec(loop)=1;
   end
  end
  
  r_non_zero = find(vec); clear vec
 
  for loop=1:length(r_non_zero)
   w(N+loop) = fzero('fn_ImDispRel_ice', ...
    [r_interval(r_non_zero(loop)-1), r_interval(r_non_zero(loop))],...
    optimset('tolX', Tol_vec(2)), [al_hat,be_hat,kap_hat])/dum_H;
  end
  
  % V2
  
%   w_hold  = mid;
%   
%   count = 1;
%   
%   while count ~= 0
%    
%    r_interval = linspace(w_hold, upp, res);
%    
%    ed = res;
%    
%    while fn_cot_dr(r_interval(ed), kap_hat) >  0
%     
%     ed = ed - 1;
%     
%    end
%    
%    if fn_cot_dr(r_interval(ed),kap_hat) + ...
%      fn_poly_dr(r_interval(ed),al_hat,be_hat) < 0
%     
%     count = 0;
%     
%    else
%     
%     count = count + 1;
%     
%    end
%    
%    if and(count ~= 0, abs(r_interval(1)-r_interval(2)) < Tol_vec(1))
%     
%     mk = 1;
%     
%     break
%     
%    end
%    
%    w_hold = r_interval(ed);
%    
%   end
%   
%   vec = zeros(res,1);
%   
%   r_interval = linspace(mid, w_hold, res);
%   
%   % ------------- %
%   
%   for loop = 2:res
%    
%    if (fn_poly_dr(r_interval(loop-1),al_hat,be_hat)+fn_cot_dr(r_interval(loop-1),kap_hat))*...
%      (fn_poly_dr(r_interval(loop),al_hat,be_hat)+fn_cot_dr(r_interval(loop),kap_hat)) < 0
%     
%     vec(loop)=1;
%     
%    end
%    
%   end
%   
%   r_non_zero = find(vec);
%   
%   step = 1;
%   
%   while step <= length(r_non_zero)
%    
%    w(length(l_non_zero) + step) = fzero('fn_ImDispRel_ice', ...
%     [r_interval(r_non_zero(step)-1), r_interval(r_non_zero(step))],...
%     optimset('tolX', Tol_vec(2)), [al_hat,be_hat,kap_hat])/dum_H;
%    
%    step = step + 1;
%    
%   end
%   
%   if mk == 1
%    
%    w(length(l_non_zero) + step) = 1i*w_hold;
%    
%   end
  
 end
 
 % Arbitrarily choose the first root
 
 w = w(1);
 
end

return

% %%%%%% SUBFUNCTIONS %%%%% %%

function out = fn_poly_dr(z,al_hat,be_hat)

out = (al_hat + be_hat*(z.^4)).*z;

return

%

function out = fn_cot_dr(z, kap)

out = kap*cot(z);

return