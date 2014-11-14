function Dave = floe_scaling(f,xi,Dmin,Dmax)
% This function computes the average floe size within a grid cell as a
% function of the maximum floe size using a bounded fractal renormalization
% groud method.

% We suggest to use Dmin >= 20. Below that value, there is no scattering
% by the floes for periods larger than 6 s (it's probably viscous however).

M  = floor(log2(Dmax/Dmin));

if isfinite(M) && M > 0
    m    = 0:M;
    N    = (1 - f).*(xi.^2.*f).^m;
    ND   = N.*Dmax./(xi.^m);
    Dave = sum(ND)./sum(N);
else
    Dave = Dmin;
end

Dave = max(Dave,Dmin);
