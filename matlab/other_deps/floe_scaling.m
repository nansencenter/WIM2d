function Dave = floe_scaling(f,xi,Dmin,Dmax)
% This function computes the average floe size within a grid cell as a
% function of the maximum floe size using a bounded fractal renormalization
% groud method.

% We suggest to use Dmin >= 20. Below that value, there is no scattering
% by the floes for periods larger than 6 s (it's probably viscous however).

Dave     = Dmax;
Dthresh  = 200;

for j=1:length(Dmax(:));
   dmax  = Dmax(j);
   M     = floor(log2(dmax/Dmin));
   if (isfinite(M) && M > 0)&(dmax<Dthresh) %>Dmin & <Dthresh
       m       = 0:M;
       N       = (1 - f).*(xi.^2.*f).^m;
       ND      = N.*dmax./(xi.^m);
       Dave(j) = sum(ND)./sum(N);
   elseif (dmax<Dmin)
       Dave(j) = Dmin;
   end
end

Dave = max(Dave,Dmin);
