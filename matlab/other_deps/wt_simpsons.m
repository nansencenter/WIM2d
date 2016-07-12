function wt = wt_simpsons(nw,dom)
%% weights for Simpson's rule

if nw==1
   disp(['Number of integration points: ',num2str(nw)]);
   error('Too few points for integration with Simpson''s rule');
elseif mod(nw,2)==0
   disp(['Number of integration points: ',num2str(nw)]);
   error('Need an odd number of points for integration with Simpson''s rule');
else
   wt_simp            = 2+zeros(nw,1);
   wt_simp([1 end])   = 1;
   wt_simp(2:2:end-1) = 4;
   %%
   wt = dom/3*wt_simp;
end
